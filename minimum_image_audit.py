#!/usr/bin/env python3
"""
minimum_image_audit.py — the definitive periodic-image test.

The earlier extent-vs-half-box heuristic is NOT diagnostic when the solute is
large relative to the box (which is the case here): it fires on genuine peptide
size, not on wrapping. This script instead computes the quantity that actually
matters:

    the minimum distance between any peptide atom and any atom of the peptide's
    own PERIODIC IMAGES (i.e. excluding the central copy).

Interpretation:
  d_img >= 2 * r_cut   -> solute cannot see its image at all (safe)
  d_img >= r_cut       -> no direct nonbonded solute-image interaction (acceptable)
  d_img <  r_cut       -> the peptide DIRECTLY interacts with its own image
                          through the real-space nonbonded terms (artefact)
  d_img <  ~4 A        -> images are in van der Waals contact (severe)

r_cut = 10 A (PME real-space cutoff used in production).

It also reports the fraction of frames below each threshold, per run, so the
artefact can be quantified rather than asserted.

Usage
-----
    python minimum_image_audit.py                     # all systems, stride 100
    python minimum_image_audit.py --stride 50
    python minimum_image_audit.py --systems GGE_reline

Writes: extension/analysis/audit/minimum_image_audit.csv
Read-only otherwise.
"""

import argparse
from pathlib import Path
from itertools import product

import numpy as np
import pandas as pd
import mdtraj as md

PROJECT_DIR = Path("~/des-peptide-study").expanduser()
EXT_DIR = PROJECT_DIR / "extension" / "trajectories_extended"
OUT_DIR = PROJECT_DIR / "extension" / "analysis" / "audit"

SYSTEMS = ["GGE_reline", "GGE_glyceline", "GGE_water",
           "CME_reline", "CME_glyceline", "CME_water",
           "YIY_reline", "YIY_glyceline", "YIY_water"]

R_CUT = 10.0   # A, PME real-space cutoff used in production


def find_prmtop(system):
    for p in (PROJECT_DIR / "systems" / system / f"{system}.prmtop",
              PROJECT_DIR / "systems" / f"{system}.prmtop"):
        if p.exists():
            return p
    raise FileNotFoundError(f"no prmtop for {system}")


def find_runs(system):
    runs = []
    for d in sorted(EXT_DIR.glob(f"{system}_*")):
        dcd = d / f"{d.name}.dcd"
        if dcd.exists():
            runs.append((d.name.split("_")[-2], dcd))
    return runs


# the 26 neighbouring image translations (exclude 0,0,0 = the central copy)
SHIFTS = np.array([s for s in product((-1, 0, 1), repeat=3) if s != (0, 0, 0)],
                  dtype=float)


def min_image_distance(xyz, box):
    """Minimum distance from the peptide to any of its 26 periodic images.

    xyz : (N, 3) peptide coords, Angstrom, molecule intact (unwrapped)
    box : (3,)   box lengths, Angstrom  (orthorhombic)

    Returns the minimum over all image translations of the minimum
    atom-atom distance between the central copy and that image.
    """
    # centre the molecule so translations are about its own centroid
    best = np.inf
    for s in SHIFTS:
        shifted = xyz + s * box            # one periodic image
        d2 = ((xyz[:, None, :] - shifted[None, :, :]) ** 2).sum(-1)
        m = np.sqrt(d2.min())
        if m < best:
            best = m
    return float(best)


def audit(system, stride=100):
    prmtop = find_prmtop(system)
    rows = []
    for label, dcd in find_runs(system):
        # box can fluctuate under the barostat -> read it per frame
        t = md.load(str(dcd), top=str(prmtop), stride=stride)
        sel = t.topology.select("protein")
        pep = t.atom_slice(sel)
        xyz = pep.xyz * 10.0                       # (F, N, 3) Angstrom
        boxes = t.unitcell_lengths * 10.0          # (F, 3) Angstrom

        d = np.array([min_image_distance(xyz[i], boxes[i])
                      for i in range(len(xyz))])

        # peptide extent per frame, for context
        ext = np.array([np.sqrt((((xyz[i][:, None, :] - xyz[i][None, :, :]) ** 2)
                                 .sum(-1)).max()) for i in range(len(xyz))])

        rows.append({
            "system": system, "start": label,
            "frames": len(d), "stride_ps": stride,
            "box_mean_A": round(float(boxes[:, 0].mean()), 2),
            "extent_mean_A": round(float(ext.mean()), 2),
            "extent_max_A": round(float(ext.max()), 2),
            "d_image_min_A": round(float(d.min()), 2),
            "d_image_mean_A": round(float(d.mean()), 2),
            "d_image_p5_A": round(float(np.percentile(d, 5)), 2),
            "pct_frames_below_cutoff": round(100.0 * (d < R_CUT).mean(), 2),
            "pct_frames_below_2cutoff": round(100.0 * (d < 2 * R_CUT).mean(), 2),
            "pct_frames_vdw_contact": round(100.0 * (d < 4.0).mean(), 2),
        })

        r = rows[-1]
        status = ("SAFE" if r["d_image_min_A"] >= 2 * R_CUT else
                  "ACCEPTABLE" if r["d_image_min_A"] >= R_CUT else
                  "** DIRECT SOLUTE-IMAGE INTERACTION")
        print(f"  {system:<14} {label:<9} "
              f"d_img min {r['d_image_min_A']:6.2f} A  mean {r['d_image_mean_A']:6.2f}  "
              f"| <{R_CUT:.0f}A in {r['pct_frames_below_cutoff']:5.1f}% frames  "
              f"| vdW contact {r['pct_frames_vdw_contact']:5.1f}%   {status}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--systems", nargs="*", default=SYSTEMS)
    ap.add_argument("--stride", type=int, default=100)
    args = ap.parse_args()

    print(f"Minimum peptide-to-own-image distance (PME real-space cutoff = {R_CUT} A)\n")
    all_rows = []
    for s in args.systems:
        try:
            all_rows.extend(audit(s, args.stride))
        except Exception as e:
            print(f"  !! {s}: {e}")

    if not all_rows:
        return
    df = pd.DataFrame(all_rows)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_DIR / "minimum_image_audit.csv", index=False)

    print("\n" + "=" * 70)
    worst = df["d_image_min_A"].min()
    n_bad = (df["d_image_min_A"] < R_CUT).sum()
    print(f"Runs with any frame below the {R_CUT:.0f} A cutoff: {n_bad} / {len(df)}")
    print(f"Closest solute-image approach anywhere: {worst:.2f} A")
    frac = df["pct_frames_below_cutoff"].mean()
    print(f"Mean % of frames with direct solute-image interaction: {frac:.1f}%")
    print(f"\nsaved: {OUT_DIR / 'minimum_image_audit.csv'}")


if __name__ == "__main__":
    main()
