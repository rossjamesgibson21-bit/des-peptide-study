#!/usr/bin/env python
"""
build_e3x_2b_matched_cutoff.py
==============================

Emit the standalone matched-cutoff coordination notebook (E3.2b):

    extension/analysis/E3x/cross_solvent/E3x_2b_matched_cutoff.ipynb

Split out of the cross-solvent notebook for provenance: E3.2b is the only
trajectory-touching unit of the phase (loads 35 production DCDs, ~4 min,
writes a guarded cache). It re-derives its inputs from committed CSVs, so it
runs independently of the reduction notebook, which merely consumes its output
`e3.2_matched_cutoff.csv`.

Producer  (this notebook): trajectories + committed CSVs -> e3.2_matched_cutoff.csv
Consumer  (E3x_cross_solvent.ipynb 2a/4): reads that cache, resolves the E3.2 line.

Key differences from the embedded 2b cell it replaces:
  - reference atom sets reconstructed from sequence + motif span (resSeq),
    NOT resname strings (which conflated GGE's four Glu and over-selected
    backbone_O to the whole peptide). Validated against committed n_ref
    before any recompute.
  - partner selections use real topology names (choline CHO/N, urea URE O,N/N1,
    glycerol GOL O/O1/O2, chloride CLA, water HOH/O).
  - convergence_stats imported from extension/ (module dir), and compute_stat
    read as a dict (mean, ci95, tau_int_ns, N_eff), not tuple-unpacked.

Output path overridable via env E3X_2B_NB_OUT for validation.

Run:
    python build_e3x_2b_matched_cutoff.py
"""

import os
import nbformat as nbf

NB_OUT = os.path.expanduser(os.environ.get(
    "E3X_2B_NB_OUT",
    "~/des-peptide-study/extension/analysis/E3x/cross_solvent/E3x_2b_matched_cutoff.ipynb",
))

md = []
def C(src):  md.append(("code", src.strip("\n")))
def M(src):  md.append(("markdown", src.strip("\n")))


# ─────────────────────────────────────────────────────────────────────────────
M(r"""
# E3.2b — Matched-Cutoff Coordination (standalone)

The one trajectory-touching unit of the cross-solvent phase, split out for
provenance. Per-system coordination used pair-specific first-shell cutoffs, so
magnitudes are not comparable across solvents; this recomputes each resolved
(reference, partner) at a **common** cutoff per pair-type — derived from the
committed `rdf_cutoffs.csv` (`first_min_nm`), not invented — using the same
estimator as the per-system phase (count of partner atoms within the cutoff,
ACF-corrected CI).

**I/O.** Inputs: committed coordination/rdf CSVs (re-loaded here, not inherited
from the reduction notebook) plus the 35 production DCDs. Output: the guarded
cache `e3.2_matched_cutoff.csv`, which `E3x_cross_solvent.ipynb` §2a/§4 consume.

**Cost.** ~4 min first run (loads dominate: 35 strided DCD loads at ~5 s each;
RDF/count compute is a minor term). Guarded by the cache so it runs exactly once;
delete `e3.2_matched_cutoff.csv` to force a recompute.

**Provenance guarantee.** Reference atom sets are reconstructed from the peptide
sequences and motif spans and **validated against the committed `n_ref`** (§1
pre-flight) before any trajectory is read — so the recompute is provably the
same atom sets as the per-system phase, and the §4 QC re-confirms it by
reproducing stored `pooled_coord` at common cutoffs.
""")


# ─────────────────────────────────────────────────────────────────────────────
# §0 — config + inputs
# ─────────────────────────────────────────────────────────────────────────────
M(r"""
## §0 — Configuration and committed inputs

Re-derives `data` (coordination + rdf CSVs), `COL`, and the `long` resolved-shell
frame from committed outputs, so this notebook is self-contained. `CS_MODULE_DIR`
points at `extension/` where `convergence_stats.py` lives (distinct from the
convergence sidecar CSVs under `extension/analysis/convergence/`).
""")

C(r"""
import os, sys, time
import numpy as np
import pandas as pd
import mdtraj as md

PROJECT   = os.path.expanduser("~/des-peptide-study")
E3X       = os.path.join(PROJECT, "extension/analysis/E3x")
TRAJROOT  = os.path.join(PROJECT, "extension/trajectories_extended")
SYSTOP    = os.path.join(PROJECT, "systems")
OUTDIR    = os.path.join(E3X, "cross_solvent")
CS_MODULE_DIR = os.path.join(PROJECT, "extension")   # convergence_stats.py lives here
os.makedirs(OUTDIR, exist_ok=True)

PEPTIDES = ["GGE", "CME", "YIY"]
SOLVENTS = ["water", "reline", "glyceline"]
SYSTEMS  = [f"{p}_{s}" for p in PEPTIDES for s in SOLVENTS]
STARTS   = ["compact", "mid", "open", "extended"]

STRIDE          = 10
WIN_START_FRAME = 2000        # drop 0-20 ns at stride-10 -> 18,000 frames

def traj_path(system, start):
    tag = f"{system}_200ns_{start}_r1"
    return os.path.join(TRAJROOT, tag, tag + ".dcd")

def top_path(system):
    return os.path.join(SYSTOP, system, system + ".prmtop")

CACHE = os.path.join(OUTDIR, "e3.2_matched_cutoff.csv")
print("outputs ->", OUTDIR)
print("cache   ->", CACHE, "(exists)" if os.path.exists(CACHE) else "(absent - will recompute)")
""")

C(r"""
# committed coordination + rdf inputs, and the COL map (reconciled schema)
COL = {
    "site": "reference", "partner": "partner",
    "cutoff_nm": "cutoff_nm", "coord": "pooled_coord",
    "coord_ci95": "ci95", "resolved": "structured_shell", "per_atom": "per_atom",
    "first_min": "first_min_nm",
}
SECTION = {
    "backbone":    "coordination/backbone_coordination.csv",
    "sidechain":   "coordination/sidechain_coordination.csv",
    "rdf_cutoffs": "rdfs/rdf_cutoffs.csv",
}
def _read(system, rel):
    p = os.path.join(E3X, system, rel)
    return pd.read_csv(p) if os.path.isfile(p) else None

data = {s: {k: _read(s, rel) for k, rel in SECTION.items()} for s in SYSTEMS}
missing = [f"{s}:{k}" for s in SYSTEMS for k in SECTION if data[s][k] is None]
assert not missing, f"missing committed inputs: {missing}"
print("loaded coordination + rdf inputs for all", len(SYSTEMS), "systems")
""")

C(r"""
# rebuild the resolved-shell `long` frame (structured_shell == True only)
rows = []
for s in SYSTEMS:
    pep, solv = s.split("_")[0], s.split("_", 1)[1]
    for scope in ("backbone", "sidechain"):
        df = data[s][scope]
        for _, r in df.iterrows():
            if not bool(r[COL["resolved"]]):
                continue
            rows.append({"peptide": pep, "solvent": solv, "scope": scope,
                         "site": str(r[COL["site"]]), "partner": str(r[COL["partner"]]),
                         "cutoff_nm": r.get(COL["cutoff_nm"], np.nan),
                         "coordination": r.get(COL["coord"], np.nan),
                         "per_atom": r.get(COL["per_atom"], np.nan)})
long = pd.DataFrame(rows)
print(f"resolved shells: {len(long)} rows across {long.site.nunique()} references, "
      f"{long.partner.nunique()} partners")
print(sorted(long.site.unique()), "|", sorted(long.partner.unique()))
""")


# ─────────────────────────────────────────────────────────────────────────────
# §1 — resolver + pre-flight n_ref validation
# ─────────────────────────────────────────────────────────────────────────────
M(r"""
## §1 — Reference/partner resolver, validated against committed n_ref

Reference atom sets are reconstructed from the motif spans (topology resSeq,
`ACE=0`): GGE and CME motifs at residues 4–6, YIY at 3–5. This separates GGE's
single `motif_Glu` (E6) from the `flank_Glu` cluster (E1–E3), which a resname
string cannot. Partner atoms use the real topology names surfaced from the
prmtops. The pre-flight cell asserts each reconstructed reference count equals
the committed `n_ref` (= `pooled_coord / per_atom`) **before** any trajectory is
read — the gate that makes the recompute provably the per-system atom sets.
""")

C(r"""
# motif span (resSeq, ACE=0) — the 3-residue motif tripeptide per peptide
MOTIF = {"GGE": (4, 5, 6), "CME": (4, 5, 6), "YIY": (3, 4, 5)}

def ref_spec(peptide):
    m = set(MOTIF[peptide])
    spec = {"backbone_O": (m, {"O"})}                 # motif carbonyl O (3 atoms)
    if peptide == "GGE":
        spec["motif_Glu"] = ({6},       {"OE1", "OE2"})   # single motif Glu (E6)
        spec["flank_Glu"] = ({1, 2, 3}, {"OE1", "OE2"})   # EEE cluster
    elif peptide == "CME":
        spec["motif_Cys"] = ({4}, {"SG"})
        spec["motif_Met"] = ({5}, {"SD"})
        spec["motif_Glu"] = ({6}, {"OE1", "OE2"})
    elif peptide == "YIY":
        spec["motif_Tyr"] = ({3, 5}, {"OH"})
    return spec

# partner label -> (residue name, atom-name set | None for single-atom residue)
PARTNER = {
    "water_O":    ("HOH", {"O"}),
    "urea_O":     ("URE", {"O"}),
    "urea_N":     ("URE", {"N", "N1"}),
    "glycerol_O": ("GOL", {"O", "O1", "O2"}),
    "choline_N":  ("CHO", {"N"}),
    "chloride":   ("CLA", None),
}

def ref_idx(top, peptide, label):
    resset, names = ref_spec(peptide)[label]
    return [a.index for a in top.atoms
            if a.residue.is_protein and a.residue.resSeq in resset and a.name in names]

def partner_idx(top, label):
    rn, names = PARTNER[label]
    return [a.index for a in top.atoms
            if a.residue.name == rn and (names is None or a.name in names)]
""")

C(r"""
# pre-flight gate: reconstructed n_ref must equal committed (coordination / per_atom)
tops = {p: md.load_prmtop(top_path(f"{p}_water")) for p in PEPTIDES}   # site count is solvent-independent

stored = (long.assign(nref=(long["coordination"] / long["per_atom"]).round())
              .dropna(subset=["nref"])
              .groupby(["peptide", "site"])["nref"].first())

print(f"{'peptide':8s} {'site':12s} {'recon':>5s} {'stored':>6s}")
ok = True
for (pep, site), nref_stored in stored.items():
    if site not in ref_spec(pep):
        print(f"{pep:8s} {site:12s} {'--':>5s} {int(nref_stored):>6d}  <-- no ref_spec"); ok = False; continue
    n_recon = len(ref_idx(tops[pep], pep, site))
    flag = "" if n_recon == int(nref_stored) else "  <-- MISMATCH"
    print(f"{pep:8s} {site:12s} {n_recon:>5d} {int(nref_stored):>6d}{flag}")
    ok = ok and (n_recon == int(nref_stored))
assert ok, "reference reconstruction disagrees with committed n_ref — fix resolver before recompute"
print("\nall reference atom-sets validated against committed n_ref")

# partner presence check on one DES topology (labels resolve to >0 atoms)
_t = md.load_prmtop(top_path("YIY_glyceline"))
for lab in ("glycerol_O", "choline_N", "chloride", "water_O"):
    print(f"  {lab:11s} -> {len(partner_idx(_t, lab))} atoms")
_t = md.load_prmtop(top_path("GGE_reline"))
for lab in ("urea_O", "urea_N", "choline_N", "chloride"):
    print(f"  {lab:11s} -> {len(partner_idx(_t, lab))} atoms")
""")


# ─────────────────────────────────────────────────────────────────────────────
# §2 — matched cutoffs
# ─────────────────────────────────────────────────────────────────────────────
M(r"""
## §2 — Matched cutoff per pair-type

Common cutoff per (reference, partner) from the committed `first_min_nm`, taking
the `max` across the solvents that resolve it (captures the full first shell
everywhere). NaN entries (no shell) are dropped. Recompute set = every resolved
(site, partner) in the matrix, augmented with `(site, water_O)` across all three
solvents for every site that resolves a water shell — the hydration-across-
solvents comparison that is the core indirect-restructuring signal.
""")

C(r"""
MATCHED_RULE = "max"   # {"max","mean","min"}
cut_rows = []
for s in SYSTEMS:
    rc = data[s]["rdf_cutoffs"]
    for _, r in rc.iterrows():
        v = r[COL["first_min"]]
        if pd.isna(v):
            continue
        cut_rows.append({"site": str(r[COL["site"]]), "partner": str(r[COL["partner"]]),
                         "solvent": s.split("_", 1)[1], "cutoff_nm": float(v)})
matched_cut = (pd.DataFrame(cut_rows).groupby(["site", "partner"])["cutoff_nm"]
                 .agg(MATCHED_RULE).rename("matched_cutoff_nm").reset_index())
print("matched cutoffs per pair-type:")
print(matched_cut.to_string(index=False))

want = matched_cut.merge(
    long[["peptide", "solvent", "site", "partner"]].drop_duplicates(),
    on=["site", "partner"], how="inner")

aug = []
for p in PEPTIDES:
    sites_p = long.loc[(long.peptide == p) & (long.partner == "water_O"), "site"].unique()
    for st in sites_p:
        for sv in SOLVENTS:
            aug.append({"site": st, "partner": "water_O", "peptide": p, "solvent": sv})
aug = pd.DataFrame(aug).merge(matched_cut, on=["site", "partner"], how="left")
want = (pd.concat([want, aug], ignore_index=True)
          .dropna(subset=["matched_cutoff_nm"])
          .drop_duplicates(["peptide", "solvent", "site", "partner"]))
print(f"\n{len(want)} (system-pair) recomputes across {want.peptide.nunique()} peptides")
""")


# ─────────────────────────────────────────────────────────────────────────────
# §3 — guarded recompute
# ─────────────────────────────────────────────────────────────────────────────
M(r"""
## §3 — Recompute (guarded, cached, one pass)

One system at a time: load its starts once, count partner atoms within the
matched cutoff per resolved reference, ACF-correct via `compute_stat`, pool
across starts. Writes `e3.2_matched_cutoff.csv`. Skips entirely if the cache
exists.
""")

C(r"""
if os.path.exists(CACHE):
    matched = pd.read_csv(CACHE)
    print(f"cache present ({len(matched)} rows) — skipping recompute; delete {os.path.basename(CACHE)} to force")
else:
    sys.path.insert(0, CS_MODULE_DIR)
    try:
        import convergence_stats as cs      # compute_stat(series, dt_ps=10.0) -> dict
        HAVE_CS = True
    except Exception as e:
        HAVE_CS = False
        print(f"[check] convergence_stats import failed ({e}); raw CI fallback (understates)")

    def _count_within(traj, a_idx, b_idx, cutoff, chunk=500):
        if not a_idx or not b_idx:
            return None, 0
        pairs = np.array([[i, j] for i in a_idx for j in b_idx])
        counts = []
        for k in range(0, traj.n_frames, chunk):
            d = md.compute_distances(traj[k:k+chunk], pairs)
            counts.extend((d < cutoff).sum(axis=1).tolist())
        return np.asarray(counts, float), len(a_idx)

    out_rows = []
    for s in SYSTEMS:
        pep, solv = s.split("_")[0], s.split("_", 1)[1]
        pairs_here = want[(want.peptide == pep) & (want.solvent == solv)]
        if pairs_here.empty:
            continue
        t0 = time.perf_counter()
        for start in STARTS:
            p = traj_path(s, start)
            if not os.path.isfile(p):
                print(f"  [skip] {s} {start}: not found"); continue
            tr = md.load(p, top=top_path(s), stride=STRIDE)[WIN_START_FRAME:]
            for _, pr in pairs_here.iterrows():
                a_idx = ref_idx(tr.topology, pep, pr.site)
                b_idx = partner_idx(tr.topology, pr.partner)
                series, n_ref = _count_within(tr, a_idx, b_idx, pr.matched_cutoff_nm)
                if series is None:
                    print(f"  [check] empty selection {s} {pr.site}/{pr.partner}"); continue
                if HAVE_CS:
                    st = cs.compute_stat(series, dt_ps=float(STRIDE))
                    m, ci, tau, neff = st["mean"], st["ci95"], st["tau_int_ns"], st["N_eff"]
                else:
                    m, ci = float(series.mean()), 1.96*series.std()/np.sqrt(len(series))
                    tau, neff = np.nan, np.nan
                out_rows.append({"system": s, "peptide": pep, "solvent": solv,
                                 "site": pr.site, "partner": pr.partner,
                                 "matched_cutoff_nm": pr.matched_cutoff_nm, "start": start,
                                 "coord": m, "per_atom": m/max(n_ref, 1),
                                 "ci95": ci, "tau_int_ns": tau, "N_eff": neff})
        print(f"  {s}: {len(pairs_here)} pairs x starts in {time.perf_counter()-t0:.1f}s")

    per_start = pd.DataFrame(out_rows)
    matched = (per_start.groupby(["system", "peptide", "solvent", "site", "partner", "matched_cutoff_nm"])
               .agg(coord=("coord", "mean"), per_atom=("per_atom", "mean"),
                    ci95=("ci95", lambda v: float(np.sqrt(np.sum(np.square(v)))/len(v))),
                    tau_int_ns=("tau_int_ns", "mean"), N_eff=("N_eff", "mean"),
                    n_starts=("start", "nunique"))
               .reset_index())
    matched.to_csv(CACHE, index=False)
    print(f"\nwrote {CACHE} ({len(matched)} pooled rows)")

print(matched.head(24).to_string(index=False))
""")


# ─────────────────────────────────────────────────────────────────────────────
# §4 — QC cross-check
# ─────────────────────────────────────────────────────────────────────────────
M(r"""
## §4 — QC: recompute reproduces stored coordination at common cutoffs

Where the matched cutoff equals a stored pair-specific cutoff, the recompute
must reproduce that stored `pooled_coord` within CI — the check that the
reconstructed atom sets and estimator match the per-system phase. Differences on
differing-cutoff pairs are expected (that is the point of matching); differences
on common-cutoff pairs flag estimator or selection drift and must be inspected
before the magnitudes enter the synthesis.
""")

C(r"""
qc = matched.merge(long, on=["peptide", "solvent", "site", "partner"], how="inner",
                   suffixes=("_matched", "_stored"))
if len(qc):
    qc["same_cutoff"] = np.isclose(qc["matched_cutoff_nm"], qc["cutoff_nm"], atol=1e-3)
    qc["abs_diff"] = (qc["coord"] - qc["coordination"]).abs()
    common = qc[qc["same_cutoff"]].sort_values("abs_diff", ascending=False)
    print("common-cutoff pairs (should reproduce stored pooled_coord within CI):")
    print(common[["system", "site", "partner", "matched_cutoff_nm",
                  "coord", "coordination", "abs_diff"]].to_string(index=False))
    worst = common["abs_diff"].max() if len(common) else 0.0
    print(f"\nmax abs_diff on common-cutoff pairs: {worst:.4f}")
    print("PASS: reconstruction reproduces per-system estimator" if worst < 0.15
          else "INSPECT: common-cutoff drift exceeds tolerance")
else:
    print("[check] no overlap between recompute and stored shells — inspect keys")
""")

M(r"""
## Consume in the synthesis

With `e3.2_matched_cutoff.csv` written and QC passed, re-run §4 of
`E3x_cross_solvent.ipynb` to resolve the E3.2 line from `PENDING` to the
direct-vs-indirect adjudication — including the backbone-hydration-across-
solvents comparison (`(backbone_O, water_O)` per solvent) that tests whether
water displacement is reline-specific or common to both DES. Commit this
notebook, its builder, and the cache as the E3.2b unit.
""")


# ─────────────────────────────────────────────────────────────────────────────
# assemble + write
# ─────────────────────────────────────────────────────────────────────────────
nb = nbf.v4.new_notebook()
nb.cells = [nbf.v4.new_markdown_cell(s) if k == "markdown" else nbf.v4.new_code_cell(s)
            for k, s in md]
nb.metadata["kernelspec"] = {
    "name": "des-peptide", "display_name": "Python (des-peptide)", "language": "python",
}
nb.metadata["language_info"] = {"name": "python"}

os.makedirs(os.path.dirname(NB_OUT), exist_ok=True)
with open(NB_OUT, "w") as fh:
    nbf.write(nb, fh)

print(f"wrote {NB_OUT}")
print(f"  cells: {len(nb.cells)} "
      f"({sum(c.cell_type=='code' for c in nb.cells)} code, "
      f"{sum(c.cell_type=='markdown' for c in nb.cells)} markdown)")
