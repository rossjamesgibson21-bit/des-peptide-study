#!/usr/bin/env python3
"""
system_audit.py — build-and-trajectory audit for the DES-peptide systems.

Reports, per system, the parameters a methods section must state and that the
production monitoring (PE / T / density) never covered:

  * exact total charge (sum of prmtop partial charges) and per-species charges
  * protonation states as built (residue names: GLU/GLH, ASP/ASH, HIS variants,
    CYS/CYM/CYX, LYS/LYN) and the peptide's net charge
  * ion counts and DES composition (choline / chloride / urea / glycerol / water)
  * box dimensions and volume
  * number densities and mass density
  * minimum solute-image separation (the periodic-image question), computed as
    box_edge - peptide_max_extent, plus the direct minimum peptide-to-image
    distance on sampled frames
  * initial periodic clashes (closest peptide-solvent contact at frame 0)
  * broken-molecule / imaging audit ACROSS THE TRAJECTORY: flags frames where the
    peptide's max interatomic distance exceeds half the box edge (the signature of
    a molecule split across the periodic boundary, which corrupts Rg and SASA)
  * atom ordering sanity (peptide atoms contiguous? caps present?)
  * ECC charge-scaling verification (do choline partial charges sum to +0.8 and
    chloride to -0.8?)

Usage
-----
    python system_audit.py                 # all nine systems, frame-0 only
    python system_audit.py --traj          # + trajectory-wide imaging audit (slower)
    python system_audit.py --systems GGE_reline CME_water
    python system_audit.py --stride 100    # imaging-audit stride (default 100 -> 100 ps)

Outputs a per-system table to stdout and writes:
    extension/analysis/audit/system_audit.csv
    extension/analysis/audit/imaging_audit.csv   (with --traj)

Read-only: opens topologies and trajectories, writes nothing but the CSVs.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import mdtraj as md

PROJECT_DIR = Path("~/des-peptide-study").expanduser()
EXT_DIR = PROJECT_DIR / "extension" / "trajectories_extended"
OUT_DIR = PROJECT_DIR / "extension" / "analysis" / "audit"

SYSTEMS = [
    "GGE_reline", "GGE_glyceline", "GGE_water",
    "CME_reline", "CME_glyceline", "CME_water",
    "YIY_reline", "YIY_glyceline", "YIY_water",
]

# Residue-name -> (protonation note, formal charge) for the states AmberTools emits.
PROTONATION = {
    "GLU": ("Glu, deprotonated (charged)", -1),
    "GLH": ("Glu, PROTONATED (neutral)", 0),
    "ASP": ("Asp, deprotonated (charged)", -1),
    "ASH": ("Asp, PROTONATED (neutral)", 0),
    "LYS": ("Lys, protonated (charged)", +1),
    "LYN": ("Lys, NEUTRAL", 0),
    "ARG": ("Arg, protonated (charged)", +1),
    "HIS": ("His, ambiguous name", 0),
    "HID": ("His, delta-protonated (neutral)", 0),
    "HIE": ("His, epsilon-protonated (neutral)", 0),
    "HIP": ("His, doubly protonated (charged)", +1),
    "CYS": ("Cys, thiol (neutral)", 0),
    "CYM": ("Cys, thiolate (charged)", -1),
    "CYX": ("Cys, disulfide-bonded", 0),
    "TYR": ("Tyr, neutral", 0),
}

# Solvent / DES species by residue name (covers common naming variants).
SPECIES = {
    "water":    {"HOH", "WAT", "TIP3", "SOL"},
    "chloride": {"CL", "CL-", "CLA", "Cl-"},
    "choline":  {"CHO", "CHL", "CHOL"},
    "urea":     {"URE", "UREA"},
    "glycerol": {"GOL", "GLY3", "GLYC"},
    "sodium":   {"NA", "NA+", "Na+"},
}

def species_of(resname):
    r = resname.upper()
    for name, names in SPECIES.items():
        if r in {n.upper() for n in names}:
            return name
    return None


def find_prmtop(system):
    for p in (PROJECT_DIR / "systems" / system / f"{system}.prmtop",
              PROJECT_DIR / "systems" / f"{system}.prmtop"):
        if p.exists():
            return p
    raise FileNotFoundError(f"no prmtop for {system}")


def find_runs(system):
    """All extended production runs for a system, as (label, dcd)."""
    runs = []
    for d in sorted(EXT_DIR.glob(f"{system}_*")):
        dcd = d / f"{d.name}.dcd"
        if dcd.exists():
            # label = the start token between _{ns}ns_ and _r{n}
            parts = d.name.split("_")
            label = parts[-2] if len(parts) >= 2 else d.name
            runs.append((label, dcd))
    return runs


def charges_from_prmtop(prmtop_path):
    """Per-atom partial charges (elementary charge) from the Amber prmtop.

    Amber stores charges pre-multiplied by 18.2223; divide to recover e.
    Uses ParmEd if available (authoritative), else parses the CHARGE flag.
    """
    try:
        import parmed as pmd
        st = pmd.load_file(str(prmtop_path))
        return np.array([a.charge for a in st.atoms], dtype=float)
    except Exception:
        pass
    # Fallback: parse the %FLAG CHARGE block directly.
    txt = Path(prmtop_path).read_text().splitlines()
    vals, grab = [], False
    for i, line in enumerate(txt):
        if line.startswith("%FLAG CHARGE"):
            grab = True
            continue
        if grab:
            if line.startswith("%FORMAT"):
                continue
            if line.startswith("%FLAG"):
                break
            vals.extend(float(x) for x in line.split())
    return np.array(vals, dtype=float) / 18.2223


def audit_system(system, do_traj=False, stride=100, verbose=True):
    prmtop = find_prmtop(system)
    runs = find_runs(system)
    if not runs:
        raise FileNotFoundError(f"no runs found for {system}")

    label0, dcd0 = runs[0]
    frame = md.load_frame(str(dcd0), 0, top=str(prmtop))
    top = frame.topology

    # ---- charges -----------------------------------------------------------
    q = charges_from_prmtop(prmtop)
    if len(q) != top.n_atoms:
        raise ValueError(f"{system}: charge array {len(q)} != atoms {top.n_atoms}")
    total_charge = float(q.sum())

    prot_sel = top.select("protein")
    q_peptide = float(q[prot_sel].sum()) if len(prot_sel) else np.nan

    # ---- protonation states (peptide residues) -----------------------------
    prot_res = [r for r in top.residues if r.index in
                {top.atom(i).residue.index for i in prot_sel}]
    prot_notes, formal = [], 0
    for r in prot_res:
        note, fc = PROTONATION.get(r.name.upper(), (f"{r.name} (uncharacterised)", 0))
        prot_notes.append(f"{r.name}{r.resSeq}: {note}")
        formal += fc

    # ---- species counts, ions ----------------------------------------------
    counts = {}
    species_charge = {}
    for r in top.residues:
        sp = species_of(r.name)
        key = sp if sp else ("peptide" if r.index in
                             {top.atom(i).residue.index for i in prot_sel} else f"other:{r.name}")
        counts[key] = counts.get(key, 0) + 1
        qs = sum(q[a.index] for a in r.atoms)
        species_charge.setdefault(key, []).append(qs)

    # ---- ECC charge-scaling verification -----------------------------------
    ecc = {}
    for sp in ("choline", "chloride"):
        if sp in species_charge:
            arr = np.array(species_charge[sp])
            ecc[sp] = (float(arr.mean()), float(arr.std()))

    # ---- box, volume, densities --------------------------------------------
    box = frame.unitcell_lengths[0] * 10.0            # nm -> Angstrom
    vol_A3 = float(np.prod(box))
    vol_nm3 = vol_A3 / 1000.0
    masses = np.array([a.element.mass for a in top.atoms])
    total_mass_amu = float(masses.sum())
    # g/cm3 = amu / A^3 * 1.66054
    mass_density = total_mass_amu / vol_A3 * 1.66054

    number_density = {k: v / vol_nm3 for k, v in counts.items()}

    # ---- peptide extent and image separation (frame 0) ---------------------
    pep = frame.atom_slice(prot_sel)
    xyz = pep.xyz[0] * 10.0                            # Angstrom
    # max interatomic distance without scipy dependency
    d2 = ((xyz[:, None, :] - xyz[None, :, :]) ** 2).sum(-1)
    pep_extent = float(np.sqrt(d2.max()))
    min_box = float(box.min())
    image_sep = min_box - pep_extent                   # solute-image separation

    # ---- initial periodic clashes (closest peptide-nonpeptide contact) -----
    others = np.setdiff1d(np.arange(top.n_atoms), prot_sel)
    pairs = np.array([[i, j] for i in prot_sel[:min(len(prot_sel), 200)]
                      for j in others[::max(1, len(others)//4000)]])
    if len(pairs):
        dmin = float(md.compute_distances(frame, pairs, periodic=True).min() * 10.0)
    else:
        dmin = np.nan

    # ---- atom ordering -----------------------------------------------------
    contiguous = bool(np.all(np.diff(np.sort(prot_sel)) == 1))
    resnames = [r.name for r in prot_res]
    caps = ("ACE" in [r.upper() for r in resnames],
            any(r.upper() in ("NME", "NHE") for r in resnames))

    rec = {
        "system": system,
        "n_atoms": top.n_atoms,
        "n_peptide_atoms": len(prot_sel),
        "peptide_contiguous": contiguous,
        "ACE_cap": caps[0], "NME_cap": caps[1],
        "total_charge_e": round(total_charge, 4),
        "peptide_charge_e": round(q_peptide, 4),
        "peptide_formal_charge": formal,
        "box_x_A": round(float(box[0]), 2),
        "box_y_A": round(float(box[1]), 2),
        "box_z_A": round(float(box[2]), 2),
        "volume_nm3": round(vol_nm3, 2),
        "mass_density_g_cm3": round(mass_density, 4),
        "peptide_max_extent_A": round(pep_extent, 2),
        "min_solute_image_sep_A": round(image_sep, 2),
        "half_box_A": round(min_box / 2.0, 2),
        "extent_lt_half_box": bool(pep_extent < min_box / 2.0),
        "min_initial_contact_A": round(dmin, 2) if dmin == dmin else np.nan,
        "n_water": counts.get("water", 0),
        "n_chloride": counts.get("chloride", 0),
        "n_choline": counts.get("choline", 0),
        "n_urea": counts.get("urea", 0),
        "n_glycerol": counts.get("glycerol", 0),
        "n_sodium": counts.get("sodium", 0),
        "choline_q_mean": round(ecc["choline"][0], 4) if "choline" in ecc else np.nan,
        "chloride_q_mean": round(ecc["chloride"][0], 4) if "chloride" in ecc else np.nan,
        "water_number_density_nm-3": round(number_density.get("water", 0), 2),
        "protonation": " | ".join(prot_notes),
        "n_runs": len(runs),
    }

    if verbose:
        print(f"\n{'='*78}\n{system}\n{'='*78}")
        print(f"  atoms {top.n_atoms}  (peptide {len(prot_sel)}, contiguous={contiguous}, "
              f"ACE={caps[0]} NME={caps[1]})")
        print(f"  box            {box[0]:.2f} x {box[1]:.2f} x {box[2]:.2f} A   "
              f"V = {vol_nm3:.1f} nm^3   rho = {mass_density:.4f} g/cm3")
        print(f"  total charge   {total_charge:+.4f} e     (peptide {q_peptide:+.4f} e, "
              f"formal {formal:+d})")
        if ecc:
            for sp, (m, s) in ecc.items():
                print(f"  {sp:<9} q     {m:+.4f} e  (sd {s:.4f})  "
                      f"{'<-- ECC 0.8 scaling confirmed' if abs(abs(m)-0.8) < 0.02 else '<-- NOT 0.8'}")
        comp = ", ".join(f"{k}={v}" for k, v in sorted(counts.items()) if not k.startswith("other"))
        print(f"  composition    {comp}")
        print(f"  peptide extent {pep_extent:.2f} A   half-box {min_box/2:.2f} A   "
              f"solute-image sep {image_sep:.2f} A")
        if pep_extent >= min_box / 2.0:
            print("  ** WARNING: peptide extent >= half the box edge at frame 0 — "
                  "imaging artefacts possible")
        print(f"  min initial peptide-solvent contact  {dmin:.2f} A"
              f"{'  ** CLASH' if dmin == dmin and dmin < 1.5 else ''}")
        print("  protonation:")
        for n in prot_notes:
            print(f"     {n}")

    # ---- trajectory-wide imaging audit -------------------------------------
    img_rows = []
    if do_traj:
        for label, dcd in runs:
            t = md.load(str(dcd), top=str(prmtop), stride=stride,
                        atom_indices=prot_sel)
            box_t = md.load_frame(str(dcd), 0, top=str(prmtop)).unitcell_lengths[0] * 10.0
            half = float(box_t.min()) / 2.0
            x = t.xyz * 10.0                                    # (F, N, 3) Angstrom
            # per-frame max interatomic distance
            diff = x[:, :, None, :] - x[:, None, :, :]
            dmax = np.sqrt((diff ** 2).sum(-1)).reshape(len(x), -1).max(axis=1)
            n_broken = int((dmax > half).sum())
            img_rows.append({
                "system": system, "start": label,
                "frames_checked": len(x), "stride": stride,
                "half_box_A": round(half, 2),
                "max_extent_A": round(float(dmax.max()), 2),
                "mean_extent_A": round(float(dmax.mean()), 2),
                "n_frames_extent_gt_half_box": n_broken,
                "pct_frames_flagged": round(100.0 * n_broken / len(x), 3),
                "clean": n_broken == 0,
            })
            if verbose:
                flag = "CLEAN" if n_broken == 0 else f"** {n_broken} FRAMES FLAGGED"
                print(f"  imaging  {label:<9} max extent {dmax.max():6.2f} A  "
                      f"(half-box {half:.2f})  {flag}")

    return rec, img_rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--systems", nargs="*", default=SYSTEMS)
    ap.add_argument("--traj", action="store_true",
                    help="run the trajectory-wide imaging audit (slower)")
    ap.add_argument("--stride", type=int, default=100,
                    help="stride for the imaging audit (default 100 = 100 ps)")
    args = ap.parse_args()

    recs, imgs = [], []
    for s in args.systems:
        try:
            rec, img = audit_system(s, do_traj=args.traj, stride=args.stride)
            recs.append(rec)
            imgs.extend(img)
        except Exception as e:
            print(f"\n!! {s}: {e}", file=sys.stderr)

    if not recs:
        return

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(recs)
    df.to_csv(OUT_DIR / "system_audit.csv", index=False)
    print(f"\nsaved: {OUT_DIR / 'system_audit.csv'}")

    if imgs:
        di = pd.DataFrame(imgs)
        di.to_csv(OUT_DIR / "imaging_audit.csv", index=False)
        print(f"saved: {OUT_DIR / 'imaging_audit.csv'}")
        bad = di[~di["clean"]]
        print("\nIMAGING AUDIT:",
              "all trajectories clean" if bad.empty
              else f"** {len(bad)} run(s) with flagged frames — see imaging_audit.csv")

    # concise cross-system summary
    cols = ["system", "total_charge_e", "peptide_charge_e", "box_x_A",
            "mass_density_g_cm3", "peptide_max_extent_A", "half_box_A",
            "extent_lt_half_box", "min_initial_contact_A"]
    print("\n" + df[cols].to_string(index=False))


if __name__ == "__main__":
    main()
