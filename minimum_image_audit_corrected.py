#!/usr/bin/env python3
"""
minimum_image_audit_corrected.py

Audit whether any non-excluded peptide atom pair is evaluated through a
non-zero periodic image inside the production real-space cutoff.

This replaces the peptide-extent/half-box heuristic.  For every selected frame
the script:

1. reconstructs the peptide as a whole molecule using its bond graph;
2. removes the actual OpenMM NonbondedForce exception pairs;
3. compares unwrapped direct pair distances with periodic pair distances;
4. identifies pairs whose nearest representation uses a non-zero image; and
5. tests whether those image-mediated distances are below the 1.0 nm cutoff.

The canonical window is 20--200 ns and the sampling stride is fixed at 100
saved frames.  This script is deliberately a SCREEN only.  A no-hit result is
reported as screening evidence and never as a definitive all-frame pass.

This is a direct-space audit only.  A clear result does not remove PME
reciprocal-space finite-size effects or the neutralizing-background concern in
net-charged systems.

Examples
--------
    # Fixed stride-100 screen of all systems
    python minimum_image_audit_corrected.py

    # Fixed stride-100 screen of selected systems
    python minimum_image_audit_corrected.py \
        --systems GGE_reline GGE_water

    # Also save per-sampled-frame data for correlation with SASA/Rg
    python minimum_image_audit_corrected.py --write-frame-data

Outputs
-------
    extension/analysis/audit/minimum_image_audit.csv
    extension/analysis/audit/minimum_image_frames/*.csv.gz  (optional)

Simulation files are opened read-only.  Only the audit outputs are written.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from collections import deque
from itertools import product
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd


SCRIPT_VERSION = "1.1.0"
DEFAULT_PROJECT_DIR = Path("~/des-peptide-study").expanduser()

SYSTEMS = [
    "GGE_reline", "GGE_glyceline", "GGE_water",
    "CME_reline", "CME_glyceline", "CME_water",
    "YIY_reline", "YIY_glyceline", "YIY_water",
]

DEFAULT_CUTOFF_NM = 1.0
DEFAULT_START_NS = 20.0
DEFAULT_END_NS = 200.0
DEFAULT_FRAME_INTERVAL_PS = 1.0
SCREEN_STRIDE = 100

SHIFTS_27 = np.asarray(list(product((-1, 0, 1), repeat=3)), dtype=np.int8)
ZERO_SHIFT_INDEX = int(np.flatnonzero(np.all(SHIFTS_27 == 0, axis=1))[0])
SHIFTS_26 = np.delete(SHIFTS_27, ZERO_SHIFT_INDEX, axis=0)


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def find_prmtop(project_dir: Path, system: str) -> Path:
    candidates = (
        project_dir / "systems" / system / f"{system}.prmtop",
        project_dir / "systems" / f"{system}.prmtop",
    )
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(f"no prmtop found for {system}")


def parse_run_name(name: str) -> tuple[str, int | None]:
    match = re.search(r"_(compact|mid|open|extended)_r(\d+)$", name)
    if match:
        return match.group(1), int(match.group(2))
    return name, None


def find_runs(project_dir: Path, system: str) -> list[dict]:
    root = project_dir / "extension" / "trajectories_extended"
    runs = []
    for directory in sorted(root.glob(f"{system}_*")):
        dcd = directory / f"{directory.name}.dcd"
        if not dcd.exists():
            continue
        start, replicate = parse_run_name(directory.name)
        metadata = directory / f"{directory.name}_metadata.json"
        runs.append({
            "run_id": directory.name,
            "start": start,
            "replicate": replicate,
            "dcd": dcd,
            "metadata": metadata if metadata.exists() else None,
        })
    return runs


def read_metadata(path: Path | None) -> dict:
    if path is None:
        return {}
    with path.open() as handle:
        return json.load(handle)


def dcd_frame_count(path: Path) -> int:
    with md.open(str(path)) as handle:
        return int(len(handle))


def atom_symbol(atom) -> str:
    if atom.element is not None:
        return str(atom.element.symbol).upper()
    letters = "".join(character for character in atom.name if character.isalpha())
    return letters[:1].upper()


def is_heavy(atom) -> bool:
    return atom_symbol(atom) != "H"


def atom_label(atom) -> str:
    residue = atom.residue
    return f"{residue.name}{residue.resSeq}:{atom.name}[{atom.index}]"


def peptide_component_indices(topology) -> np.ndarray:
    """Return protein atoms plus covalently connected ACE/NME cap atoms."""
    seeds = set(int(index) for index in topology.select("protein"))
    if not seeds:
        raise ValueError("MDTraj protein selection is empty")

    adjacency = [set() for _ in range(topology.n_atoms)]
    for bond in topology.bonds:
        i, j = bond.atom1.index, bond.atom2.index
        adjacency[i].add(j)
        adjacency[j].add(i)

    selected = set(seeds)
    queue = deque(seeds)
    while queue:
        atom_index = queue.popleft()
        for neighbour in adjacency[atom_index]:
            if neighbour not in selected:
                selected.add(neighbour)
                queue.append(neighbour)
    return np.asarray(sorted(selected), dtype=np.int32)


def unwrapping_tree(peptide_topology) -> tuple[list[int], list[tuple[int, int]]]:
    adjacency = [set() for _ in range(peptide_topology.n_atoms)]
    for bond in peptide_topology.bonds:
        i, j = bond.atom1.index, bond.atom2.index
        adjacency[i].add(j)
        adjacency[j].add(i)

    roots = []
    edges = []
    visited = set()
    for root in range(peptide_topology.n_atoms):
        if root in visited:
            continue
        roots.append(root)
        visited.add(root)
        queue = deque([root])
        while queue:
            parent = queue.popleft()
            for child in sorted(adjacency[parent]):
                if child in visited:
                    continue
                visited.add(child)
                edges.append((parent, child))
                queue.append(child)

    if len(roots) != 1:
        raise ValueError(
            f"peptide selection has {len(roots)} bonded components; "
            "cannot define one whole peptide safely"
        )
    return roots, edges


def make_peptide_whole(
    xyz_nm: np.ndarray,
    box_vectors_nm: np.ndarray,
    roots: list[int],
    edges: list[tuple[int, int]],
) -> tuple[np.ndarray, np.ndarray]:
    """Make a bonded peptide whole using each frame's complete box matrix.

    The 27 neighboring lattice translations are tested explicitly for every
    bond.  This works for orthorhombic and triclinic unit cells.
    """
    xyz = np.asarray(xyz_nm, dtype=np.float64)
    boxes = np.asarray(box_vectors_nm, dtype=np.float64)
    if xyz.ndim != 3 or boxes.shape != (len(xyz), 3, 3):
        raise ValueError("coordinate or box array has an unexpected shape")
    if not np.isfinite(xyz).all() or not np.isfinite(boxes).all():
        raise ValueError("trajectory contains NaN/Inf coordinates or box vectors")
    if np.any(np.abs(np.linalg.det(boxes)) < 1e-10):
        raise ValueError("trajectory contains a singular periodic box")

    lattice_vectors = np.einsum("si,fij->fsj", SHIFTS_27, boxes)
    whole = np.empty_like(xyz)
    split_frame = np.zeros(len(xyz), dtype=bool)

    for root in roots:
        whole[:, root, :] = xyz[:, root, :]

    frame_indices = np.arange(len(xyz))
    for parent, child in edges:
        raw_delta = xyz[:, child, :] - xyz[:, parent, :]
        candidates = raw_delta[:, None, :] - lattice_vectors
        squared = np.einsum("fsi,fsi->fs", candidates, candidates)
        best_shift = np.argmin(squared, axis=1)
        split_frame |= best_shift != ZERO_SHIFT_INDEX
        whole[:, child, :] = (
            whole[:, parent, :] + candidates[frame_indices, best_shift]
        )

    # Translation is immaterial to pair distances; centering keeps numbers small.
    whole -= whole.mean(axis=1, keepdims=True)
    return whole, split_frame


def openmm_allowed_pairs(
    prmtop_path: Path,
    peptide_indices: np.ndarray,
    peptide_topology,
    cutoff_nm: float,
) -> tuple[np.ndarray, np.ndarray, bool]:
    """Remove every pair listed as an OpenMM NonbondedForce exception."""
    from openmm import NonbondedForce, unit
    from openmm.app import AmberPrmtopFile, HBonds, PME

    amber = AmberPrmtopFile(str(prmtop_path))
    system = amber.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=cutoff_nm * unit.nanometer,
        constraints=HBonds,
    )
    nonbonded = next(
        force for force in system.getForces()
        if isinstance(force, NonbondedForce)
    )

    global_to_local = {
        int(global_index): local_index
        for local_index, global_index in enumerate(peptide_indices)
    }
    exceptions = set()
    for exception_index in range(nonbonded.getNumExceptions()):
        atom_i, atom_j, *_ = nonbonded.getExceptionParameters(exception_index)
        global_i, global_j = int(atom_i), int(atom_j)
        if global_i in global_to_local and global_j in global_to_local:
            local_i = global_to_local[global_i]
            local_j = global_to_local[global_j]
            exceptions.add(tuple(sorted((local_i, local_j))))

    pairs = np.asarray([
        (i, j)
        for i in range(peptide_topology.n_atoms)
        for j in range(i + 1, peptide_topology.n_atoms)
        if (i, j) not in exceptions
    ], dtype=np.int32)
    if not len(pairs):
        raise ValueError("no non-exception peptide atom pairs were found")

    heavy_atoms = np.asarray([
        is_heavy(atom) for atom in peptide_topology.atoms
    ])
    heavy_pair_mask = heavy_atoms[pairs[:, 0]] & heavy_atoms[pairs[:, 1]]
    exceptions_use_pbc = bool(
        nonbonded.getExceptionsUsePeriodicBoundaryConditions()
    )
    return pairs, heavy_pair_mask, exceptions_use_pbc


def episode_statistics(flags: np.ndarray, sample_interval_ps: float) -> tuple[int, float]:
    flags = np.asarray(flags, dtype=bool)
    if not len(flags) or not flags.any():
        return 0, 0.0
    padded = np.concatenate(([False], flags, [False])).astype(np.int8)
    changes = np.diff(padded)
    starts = np.flatnonzero(changes == 1)
    ends = np.flatnonzero(changes == -1)
    lengths = ends - starts
    return int(len(lengths)), float(lengths.max() * sample_interval_ps / 1000.0)


def finite_summary(values_nm: np.ndarray, prefix: str) -> dict:
    finite = np.asarray(values_nm, dtype=float)
    finite = finite[np.isfinite(finite)] * 10.0
    if not len(finite):
        return {
            f"{prefix}_min_A": np.nan,
            f"{prefix}_p05_A": np.nan,
            f"{prefix}_median_A": np.nan,
        }
    return {
        f"{prefix}_min_A": float(np.min(finite)),
        f"{prefix}_p05_A": float(np.percentile(finite, 5)),
        f"{prefix}_median_A": float(np.median(finite)),
    }


def recurrent_pair_labels(
    counts: np.ndarray,
    pairs: np.ndarray,
    topology,
    limit: int = 5,
) -> str:
    nonzero = np.flatnonzero(counts)
    if not len(nonzero):
        return ""
    ranked = nonzero[np.argsort(counts[nonzero])[::-1]][:limit]
    labels = []
    for pair_index in ranked:
        atom_i, atom_j = pairs[pair_index]
        labels.append(
            f"{atom_label(topology.atom(int(atom_i)))}--"
            f"{atom_label(topology.atom(int(atom_j)))}:"
            f"{int(counts[pair_index])}"
        )
    return " | ".join(labels)


def audit_run(
    system: str,
    run: dict,
    prmtop: Path,
    peptide_indices: np.ndarray,
    peptide_topology,
    pairs: np.ndarray,
    heavy_pair_mask: np.ndarray,
    roots: list[int],
    unwrap_edges: list[tuple[int, int]],
    exceptions_use_pbc: bool,
    start_ns: float,
    end_ns: float,
    frame_interval_ps: float,
    cutoff_nm: float,
    chunk_size: int,
    write_frame_data: bool,
    frame_output_dir: Path,
) -> dict:
    stride = SCREEN_STRIDE
    dcd = run["dcd"]
    metadata = read_metadata(run["metadata"])
    total_frames = dcd_frame_count(dcd)

    start_frame = int(round(start_ns * 1000.0 / frame_interval_ps))
    end_frame_exclusive = int(round(end_ns * 1000.0 / frame_interval_ps))
    if start_frame < 0 or end_frame_exclusive <= start_frame:
        raise ValueError("invalid analysis time window")
    if total_frames < end_frame_exclusive:
        raise ValueError(
            f"requested {end_ns:g} ns requires {end_frame_exclusive} frames, "
            f"but {run['run_id']} has {total_frames}"
        )

    metadata_frames = metadata.get("n_frames", metadata.get("expected_frames"))
    if metadata_frames is not None and int(metadata_frames) != total_frames:
        raise ValueError(
            f"DCD has {total_frames} frames but metadata reports {metadata_frames}"
        )

    all_shortest = []
    heavy_shortest = []
    all_contact_counts = []
    heavy_contact_counts = []
    raw_split_flags = []
    box_lengths = []
    box_volumes = []
    raw_frame_indices = []
    times_ns = []
    nonexception_extent = []
    pair_contact_totals = np.zeros(len(pairs), dtype=np.int64)

    global_minimum_nm = np.inf
    global_pair_index = -1
    global_frame = -1
    sampled_read = 0

    iterator = md.iterload(
        str(dcd),
        top=str(prmtop),
        atom_indices=peptide_indices,
        skip=start_frame,
        stride=stride,
        chunk=chunk_size,
    )

    for trajectory in iterator:
        n_loaded = len(trajectory)
        frame_indices = start_frame + (
            sampled_read + np.arange(n_loaded, dtype=np.int64)
        ) * stride
        sampled_read += n_loaded

        keep = frame_indices < end_frame_exclusive
        if not np.any(keep):
            break
        reached_end = not np.all(keep)
        if reached_end:
            trajectory = trajectory[keep]
            frame_indices = frame_indices[keep]

        boxes = trajectory.unitcell_vectors
        if boxes is None:
            raise ValueError("DCD does not contain periodic box vectors")
        boxes = np.asarray(boxes, dtype=np.float64)

        # The nearest nonzero lattice-vector length must exceed 2*cutoff for a
        # unique real-space image under the minimum-image convention.
        lattice_vectors = np.einsum("si,fij->fsj", SHIFTS_26, boxes)
        shortest_lattice = np.linalg.norm(lattice_vectors, axis=2).min(axis=1)
        if np.any(shortest_lattice <= 2.0 * cutoff_nm):
            raise ValueError("a box vector is not greater than twice the cutoff")

        whole, split_flags = make_peptide_whole(
            trajectory.xyz, boxes, roots, unwrap_edges
        )
        trajectory.xyz = whole.astype(np.float32)

        direct = md.compute_distances(
            trajectory, pairs, periodic=False, opt=True
        ).astype(np.float64, copy=False)
        periodic = md.compute_distances(
            trajectory, pairs, periodic=True, opt=True
        ).astype(np.float64, copy=False)

        uses_nonzero_image = periodic < (direct - 1e-6)
        image_distances = np.where(uses_nonzero_image, periodic, np.inf)
        contacts = uses_nonzero_image & (periodic < cutoff_nm)

        frame_pair = np.argmin(image_distances, axis=1)
        frame_minimum = image_distances[np.arange(len(trajectory)), frame_pair]
        no_image_shortcut = ~np.isfinite(frame_minimum)
        frame_pair[no_image_shortcut] = -1

        heavy_indices = np.flatnonzero(heavy_pair_mask)
        heavy_images = image_distances[:, heavy_indices]
        heavy_minimum = heavy_images.min(axis=1)
        heavy_contacts = contacts[:, heavy_indices]

        all_shortest.append(frame_minimum)
        heavy_shortest.append(heavy_minimum)
        all_contact_counts.append(contacts.sum(axis=1))
        heavy_contact_counts.append(heavy_contacts.sum(axis=1))
        raw_split_flags.append(split_flags)
        pair_contact_totals += contacts.sum(axis=0)
        nonexception_extent.append(direct.max(axis=1))

        lengths = np.linalg.norm(boxes, axis=2)
        box_lengths.append(lengths)
        box_volumes.append(np.abs(np.linalg.det(boxes)))
        raw_frame_indices.append(frame_indices)
        nominal_times = frame_indices * frame_interval_ps / 1000.0
        times_ns.append(nominal_times)

        local_frame = int(np.argmin(frame_minimum))
        local_value = float(frame_minimum[local_frame])
        if np.isfinite(local_value) and local_value < global_minimum_nm:
            global_minimum_nm = local_value
            global_pair_index = int(frame_pair[local_frame])
            global_frame = int(frame_indices[local_frame])

        if reached_end:
            break

    if not all_shortest:
        raise ValueError("no trajectory frames remained in the requested window")

    all_shortest = np.concatenate(all_shortest)
    heavy_shortest = np.concatenate(heavy_shortest)
    all_contact_counts = np.concatenate(all_contact_counts)
    heavy_contact_counts = np.concatenate(heavy_contact_counts)
    raw_split_flags = np.concatenate(raw_split_flags)
    nonexception_extent = np.concatenate(nonexception_extent)
    box_lengths = np.concatenate(box_lengths)
    box_volumes = np.concatenate(box_volumes)
    raw_frame_indices = np.concatenate(raw_frame_indices)
    times_ns = np.concatenate(times_ns)

    all_contact_frames = all_contact_counts > 0
    heavy_contact_frames = heavy_contact_counts > 0
    sample_interval_ps = frame_interval_ps * stride
    all_episodes, longest_all_episode_ns = episode_statistics(
        all_contact_frames, sample_interval_ps
    )
    heavy_episodes, longest_heavy_episode_ns = episode_statistics(
        heavy_contact_frames, sample_interval_ps
    )

    if all_contact_frames.any():
        status = "CONTACTS_DETECTED_IN_STRIDE_100_SCREEN"
    else:
        status = "NO_HITS_IN_STRIDE_100_SCREEN"

    if global_pair_index >= 0:
        atom_i, atom_j = pairs[global_pair_index]
        global_pair_label = (
            f"{atom_label(peptide_topology.atom(int(atom_i)))}--"
            f"{atom_label(peptide_topology.atom(int(atom_j)))}"
        )
        global_time_ns = global_frame * frame_interval_ps / 1000.0
        global_distance_A = global_minimum_nm * 10.0
    else:
        global_pair_label = ""
        global_time_ns = np.nan
        global_distance_A = np.nan

    eligible_frames = end_frame_exclusive - start_frame
    frames_checked = len(all_shortest)
    record = {
        "script_version": SCRIPT_VERSION,
        "command": " ".join(sys.argv),
        "system": system,
        "run_id": run["run_id"],
        "start": run["start"],
        "replicate": run["replicate"],
        "status": status,
        "dcd": str(dcd),
        "dcd_size_bytes": dcd.stat().st_size,
        "dcd_mtime_utc": pd.Timestamp(
            dcd.stat().st_mtime, unit="s", tz="UTC"
        ).isoformat(),
        "prmtop": str(prmtop),
        "prmtop_sha256": sha256_file(prmtop),
        "metadata": str(run["metadata"]) if run["metadata"] else "not_found",
        "metadata_verified": metadata.get("verified", np.nan),
        "metadata_frames": metadata_frames,
        "metadata_seed": metadata.get("seed", np.nan),
        "total_dcd_frames": total_frames,
        "eligible_saved_frames": eligible_frames,
        "frames_checked": frames_checked,
        "saved_frame_coverage_pct": 100.0 * frames_checked / eligible_frames,
        "analysis_start_ns": start_ns,
        "analysis_end_ns": end_ns,
        "frame_interval_ps": frame_interval_ps,
        "stride": stride,
        "sample_interval_ps": sample_interval_ps,
        "screening_only": True,
        "cutoff_A": cutoff_nm * 10.0,
        "n_nonexception_pairs": len(pairs),
        "n_heavy_nonexception_pairs": int(heavy_pair_mask.sum()),
        "exceptions_use_periodic_boundary_conditions": exceptions_use_pbc,
        **finite_summary(all_shortest, "periodic_shortcut_all"),
        **finite_summary(heavy_shortest, "periodic_shortcut_heavy"),
        "frames_with_all_atom_contact": int(all_contact_frames.sum()),
        "pct_frames_with_all_atom_contact": float(
            100.0 * all_contact_frames.mean()
        ),
        "all_atom_contact_episodes": all_episodes,
        "longest_sampled_all_atom_episode_ns": longest_all_episode_ns,
        "frames_with_heavy_atom_contact": int(heavy_contact_frames.sum()),
        "pct_frames_with_heavy_atom_contact": float(
            100.0 * heavy_contact_frames.mean()
        ),
        "heavy_atom_contact_episodes": heavy_episodes,
        "longest_sampled_heavy_atom_episode_ns": longest_heavy_episode_ns,
        "mean_contact_pairs_per_frame": float(all_contact_counts.mean()),
        "max_contact_pairs_in_frame": int(all_contact_counts.max()),
        "raw_coordinate_split_frames": int(raw_split_flags.sum()),
        "pct_raw_coordinate_split_frames": float(100.0 * raw_split_flags.mean()),
        "global_minimum_shortcut_A": global_distance_A,
        "global_minimum_frame": global_frame,
        "global_minimum_time_ns": global_time_ns,
        "global_minimum_pair": global_pair_label,
        "recurrent_contact_pairs": recurrent_pair_labels(
            pair_contact_totals, pairs, peptide_topology
        ),
        "nonexception_extent_mean_A": float(
            nonexception_extent.mean() * 10.0
        ),
        "nonexception_extent_max_A": float(
            nonexception_extent.max() * 10.0
        ),
        "box_a_min_A": float(box_lengths[:, 0].min() * 10.0),
        "box_a_median_A": float(np.median(box_lengths[:, 0]) * 10.0),
        "box_a_max_A": float(box_lengths[:, 0].max() * 10.0),
        "box_b_min_A": float(box_lengths[:, 1].min() * 10.0),
        "box_b_median_A": float(np.median(box_lengths[:, 1]) * 10.0),
        "box_b_max_A": float(box_lengths[:, 1].max() * 10.0),
        "box_c_min_A": float(box_lengths[:, 2].min() * 10.0),
        "box_c_median_A": float(np.median(box_lengths[:, 2]) * 10.0),
        "box_c_max_A": float(box_lengths[:, 2].max() * 10.0),
        "volume_min_nm3": float(box_volumes.min()),
        "volume_median_nm3": float(np.median(box_volumes)),
        "volume_max_nm3": float(box_volumes.max()),
    }

    if write_frame_data:
        frame_output_dir.mkdir(parents=True, exist_ok=True)
        frame_path = frame_output_dir / f"{run['run_id']}_image_frames.csv.gz"
        pd.DataFrame({
            "raw_frame": raw_frame_indices,
            "nominal_time_ns": times_ns,
            "periodic_shortcut_all_A": np.where(
                np.isfinite(all_shortest), all_shortest * 10.0, np.nan
            ),
            "periodic_shortcut_heavy_A": np.where(
                np.isfinite(heavy_shortest), heavy_shortest * 10.0, np.nan
            ),
            "all_atom_contact_pairs": all_contact_counts,
            "heavy_atom_contact_pairs": heavy_contact_counts,
            "raw_coordinates_split": raw_split_flags,
            "box_a_A": box_lengths[:, 0] * 10.0,
            "box_b_A": box_lengths[:, 1] * 10.0,
            "box_c_A": box_lengths[:, 2] * 10.0,
            "volume_nm3": box_volumes,
        }).to_csv(frame_path, index=False, compression="gzip")
        record["frame_data"] = str(frame_path)
    else:
        record["frame_data"] = "not_written"

    return record


def prepare_system(
    project_dir: Path,
    system: str,
    cutoff_nm: float,
) -> tuple[Path, list[dict], np.ndarray, object, np.ndarray, np.ndarray, list, list, bool]:
    prmtop = find_prmtop(project_dir, system)
    runs = find_runs(project_dir, system)
    if not runs:
        raise FileNotFoundError(f"no DCD runs found for {system}")

    first_frame = md.load_frame(
        str(runs[0]["dcd"]), 0, top=str(prmtop)
    )
    peptide_indices = peptide_component_indices(first_frame.topology)
    peptide_frame = first_frame.atom_slice(peptide_indices)
    peptide_topology = peptide_frame.topology
    roots, unwrap_edges = unwrapping_tree(peptide_topology)
    pairs, heavy_pair_mask, exceptions_use_pbc = openmm_allowed_pairs(
        prmtop,
        peptide_indices,
        peptide_topology,
        cutoff_nm,
    )
    return (
        prmtop,
        runs,
        peptide_indices,
        peptide_topology,
        pairs,
        heavy_pair_mask,
        roots,
        unwrap_edges,
        exceptions_use_pbc,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", type=Path, default=DEFAULT_PROJECT_DIR)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--systems", nargs="*", default=SYSTEMS)
    parser.add_argument("--chunk", type=int, default=100)
    parser.add_argument("--start-ns", type=float, default=DEFAULT_START_NS)
    parser.add_argument("--end-ns", type=float, default=DEFAULT_END_NS)
    parser.add_argument(
        "--frame-interval-ps",
        type=float,
        default=DEFAULT_FRAME_INTERVAL_PS,
        help="time represented by one saved DCD frame",
    )
    parser.add_argument("--cutoff-nm", type=float, default=DEFAULT_CUTOFF_NM)
    parser.add_argument("--write-frame-data", action="store_true")
    parser.add_argument(
        "--allow-partial",
        action="store_true",
        help="exit zero even if a requested system or run fails",
    )
    args = parser.parse_args()

    if args.chunk < 1:
        parser.error("--chunk must be a positive integer")
    if args.frame_interval_ps <= 0 or args.cutoff_nm <= 0:
        parser.error("frame interval and cutoff must be positive")
    if args.end_ns <= args.start_ns:
        parser.error("--end-ns must be greater than --start-ns")

    project_dir = args.project_dir.expanduser().resolve()
    output_dir = (
        args.out_dir.expanduser().resolve()
        if args.out_dir is not None
        else project_dir / "extension" / "analysis" / "audit"
    )
    frame_output_dir = output_dir / "minimum_image_frames"

    rows = []
    errors = []

    print(
        f"Direct periodic-image audit: {args.start_ns:g}--{args.end_ns:g} ns, "
        f"cutoff={args.cutoff_nm * 10.0:g} A, stride={SCREEN_STRIDE}"
    )
    print("Mode: STRIDE-100 SCREEN (not a full-frame pass)\n")

    for system in args.systems:
        try:
            (
                prmtop,
                runs,
                peptide_indices,
                peptide_topology,
                pairs,
                heavy_pair_mask,
                roots,
                unwrap_edges,
                exceptions_use_pbc,
            ) = prepare_system(project_dir, system, args.cutoff_nm)
        except Exception as error:
            message = f"{system}: {error}"
            errors.append(message)
            rows.append({
                "system": system,
                "status": "ERROR",
                "error": repr(error),
            })
            print(f"  ERROR {message}", file=sys.stderr)
            continue

        for run in runs:
            try:
                record = audit_run(
                    system=system,
                    run=run,
                    prmtop=prmtop,
                    peptide_indices=peptide_indices,
                    peptide_topology=peptide_topology,
                    pairs=pairs,
                    heavy_pair_mask=heavy_pair_mask,
                    roots=roots,
                    unwrap_edges=unwrap_edges,
                    exceptions_use_pbc=exceptions_use_pbc,
                    start_ns=args.start_ns,
                    end_ns=args.end_ns,
                    frame_interval_ps=args.frame_interval_ps,
                    cutoff_nm=args.cutoff_nm,
                    chunk_size=args.chunk,
                    write_frame_data=args.write_frame_data,
                    frame_output_dir=frame_output_dir,
                )
                rows.append(record)
                minimum = record["global_minimum_shortcut_A"]
                minimum_text = f"{minimum:.3f} A" if np.isfinite(minimum) else "none"
                print(
                    f"  {system:<14} {run['start']:<9} "
                    f"{record['status']:<29} "
                    f"min shortcut={minimum_text:<12} "
                    f"contact frames={record['frames_with_all_atom_contact']}/"
                    f"{record['frames_checked']} "
                    f"({record['pct_frames_with_all_atom_contact']:.3f}%)"
                )
            except Exception as error:
                message = f"{system}/{run['run_id']}: {error}"
                errors.append(message)
                rows.append({
                    "system": system,
                    "run_id": run["run_id"],
                    "start": run["start"],
                    "status": "ERROR",
                    "error": repr(error),
                })
                print(f"  ERROR {message}", file=sys.stderr)

    if not rows:
        print("No runs were found.", file=sys.stderr)
        return 2

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "minimum_image_audit.csv"
    table = pd.DataFrame(rows)
    table.to_csv(output_path, index=False)

    successful = table[table["status"] != "ERROR"]
    print("\n" + "=" * 78)
    if not successful.empty:
        contacts = successful[successful["frames_with_all_atom_contact"] > 0]
        print(f"Runs audited successfully: {len(successful)}")
        print(f"Runs with sampled direct contacts: {len(contacts)}")
        print(
            "A no-hit result is stride-100 screening evidence only; this "
            "script does not perform an all-frame audit."
        )
        print(
            "This audit does not assess reciprocal-space PME or net-charge "
            "finite-size effects."
        )
    print(f"Saved: {output_path}")

    if errors and not args.allow_partial:
        print("\nAudit completed with errors:", file=sys.stderr)
        for message in errors:
            print(f"  - {message}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
