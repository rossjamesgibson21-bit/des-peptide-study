#!/usr/bin/env python
"""
build_e3x_cross_solvent.py
==========================

Emit the single cross-solvent synthesis notebook for the E3.x phase:

    extension/analysis/E3x/cross_solvent/E3x_cross_solvent.ipynb

Mirrors the per-system `build_e3x_*.py` pattern (nbformat, one logical
concern per cell, run on the `des-peptide` kernel). The notebook is a
reduction over already-computed E3.x section outputs, with a single
trajectory-touching section (E3.2b, matched-cutoff coordination) that is
guarded, cached, and placed last so the cheap sections verify first.

Column names reconciled against the real committed CSVs (2024 schema):
summary -> sasa_pooled_nm2/sasa_ci95/rg_pooled_nm/rg_ci95;
coordination -> reference/partner/cutoff_nm/pooled_coord/ci95/N_eff_pooled/
per_atom/structured_shell/note; rdf_cutoffs -> first_min_nm/cutoff_source;
per-start hbonds -> start/n_episodes/n_bonds; hbond_life -> n_distinct_bonds.

Section order in the emitted notebook:
    0   Load + inventory (schema surfaced before any downstream consumption)
    1   E3.1  SASA significance      (reduction; tier-aware within-peptide row)
    2a  E3.2  qualitative matrix     (reduction; present/absent vs water null)
    3   E3.3  YIY backbone H-bonds   (reduction; per-start; direction robustness)
    4   Synthesis scaffold          (deterministic from the section outputs)
    2b  E3.2  matched-cutoff coord   (trajectory recompute; guarded + cached)

Output path is overridable via env E3X_NB_OUT for validation.

Run:
    python build_e3x_cross_solvent.py
"""

import os
import nbformat as nbf

NB_OUT = os.path.expanduser(os.environ.get(
    "E3X_NB_OUT",
    "~/des-peptide-study/extension/analysis/E3x/cross_solvent/E3x_cross_solvent.ipynb",
))

md = []   # (kind, source) tuples, assembled in notebook order
def C(src):  md.append(("code", src.strip("\n")))
def M(src):  md.append(("markdown", src.strip("\n")))


# ─────────────────────────────────────────────────────────────────────────────
# Title
# ─────────────────────────────────────────────────────────────────────────────
M(r"""
# E3.x — Cross-Solvent Synthesis

Reads the nine per-system E3.x summaries against one another, at matched cutoffs,
to frame the logged observations as **hypotheses supported by the data** — not as
conclusions. Every number is drawn from a committed E3.x section output (CSV or
convergence sidecar); nothing is parsed from the summary markdown, and the
Phase-3/Phase-4 strided SASA columns are not used here.

Three questions, one section each:

- **E3.1** — does the DES-induced SASA change survive corrected CIs across the matrix? (tier-aware, within-peptide row)
- **E3.2** — is the SASA change accompanied by direct DES contact, or consistent with indirect restructuring? (qualitative matrix + matched-cutoff magnitudes)
- **E3.3** — is the YIY reline backbone-H-bond pattern a genuine cross-solvent difference, and in which direction?

Run order: §0 → §1 → §2a → §3 → §4 are pure reductions and complete in a few
minutes. **§2b is the only trajectory-touching section** (guarded + cached,
~5–10 min on first run); it lives at the end so the cheap sections verify first.

Conventions carried from the phase: coordination on pair-specific RDF-derived
first-shell cutoffs with ACF-corrected 95% CIs (`N_eff = N/τ_int`); SASA/Rg from
the E2.2 pooled leave-one-out framework, reported per tier; backbone H-bonds at
1 ps resolution, backbone N→O scope, continuous episodes. No external citations —
`[REF]` placeholders for anything not in a committed output.
""")


# ─────────────────────────────────────────────────────────────────────────────
# §0 — config + paths
# ─────────────────────────────────────────────────────────────────────────────
M(r"""
## §0 — Configuration, load, inventory

Paths, the system×tier map, and the analysis-window parameters confirmed by the
calibration probe (stride-10, 18,000 frames retained after dropping 0–20 ns).
`RDF_CURVES_AVAILABLE` is `False` — the inventory found `rdfs/` holds derived
cutoffs only, no g(r) arrays — which is why the matched-cutoff magnitudes in §2b
require a trajectory recompute rather than a cheap re-integration.
""")

C(r"""
import os, sys, json, glob, time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PROJECT  = os.path.expanduser("~/des-peptide-study")
E3X      = os.path.join(PROJECT, "extension/analysis/E3x")
CONV     = os.path.join(PROJECT, "extension/analysis/convergence")
TRAJROOT = os.path.join(PROJECT, "extension/trajectories_extended")
SYSTOP   = os.path.join(PROJECT, "systems")                 # {SYSTEM}/{SYSTEM}.prmtop
OUTDIR   = os.path.join(E3X, "cross_solvent")
os.makedirs(OUTDIR, exist_ok=True)

PEPTIDES = ["GGE", "CME", "YIY"]
SOLVENTS = ["water", "reline", "glyceline"]
SYSTEMS  = [f"{p}_{s}" for p in PEPTIDES for s in SOLVENTS]
STARTS   = ["compact", "mid", "open", "extended"]

TIER = {
    "GGE_reline": "A", "GGE_glyceline": "A", "CME_water": "A",
    "YIY_water": "B", "YIY_reline": "B", "YIY_glyceline": "B",
    "CME_reline": "C",
    "GGE_water": "D", "CME_glyceline": "D",
}

# analysis window (confirmed by calibration): 10 ps/frame; drop 0-20 ns
STRIDE          = 10
WIN_START_FRAME = 2000        # 2000 strided frames == 20 ns -> 18,000 retained
RDF_CURVES_AVAILABLE = False  # inventory: rdfs/ holds cutoffs only

# trajectory path convention (confirmed by the size/glob probe)
def traj_path(system, start):
    tag = f"{system}_200ns_{start}_r1"
    return os.path.join(TRAJROOT, tag, tag + ".dcd")

def top_path(system):
    return os.path.join(SYSTOP, system, system + ".prmtop")

print("outputs ->", OUTDIR)
print("systems  :", ", ".join(SYSTEMS))
""")

M(r"""
### §0.1 — Load section outputs and build the manifest

One row per system; a `True`/`False` per expected section file. This is the
gate: any `False` is a gap to backfill before the synthesis leans on it.
""")

C(r"""
SECTION_FILES = {
    "summary":        "{S}_summary.csv",
    "backbone":       "coordination/backbone_coordination.csv",
    "sidechain":      "coordination/sidechain_coordination.csv",
    "rdf_cutoffs":    "rdfs/rdf_cutoffs.csv",
    "hbond_life":     "hbonds/hbond_lifetimes.csv",
    "hbond_perstart": "hbonds/per_start_hbond_summary.csv",
}
# tier-specific SASA decompositions (verify the location/name against disk)
TIER_FILES = {
    "D": "{S}_basin_sasa.csv",   # multi-basin: separated SASA states
    "C": "{S}_trap_sasa.csv",    # one-trap: pool vs consensus + inclusion shift
}

def _read(system, rel):
    p = os.path.join(E3X, system, rel)
    return pd.read_csv(p) if os.path.isfile(p) else None

data = {s: {} for s in SYSTEMS}
for s in SYSTEMS:
    for key, tmpl in SECTION_FILES.items():
        data[s][key] = _read(s, tmpl.format(S=s))
    tf = TIER_FILES.get(TIER[s])
    data[s]["tier_sasa"] = _read(s, tf.format(S=s)) if tf else None

manifest = pd.DataFrame([
    {"system": s, "tier": TIER[s],
     **{k: (data[s][k] is not None) for k in SECTION_FILES},
     "tier_sasa": (data[s]["tier_sasa"] is not None) if TIER[s] in TIER_FILES else "—"}
    for s in SYSTEMS
])
manifest.to_csv(os.path.join(OUTDIR, "manifest.csv"), index=False)
print(manifest.to_string(index=False))
assert manifest[list(SECTION_FILES)].all().all(), "missing section output(s) — backfill before continuing"
""")

M(r"""
### §0.2 — Surface the schema before consuming columns

The downstream cells reference columns through the `COL` map below. Print the
real columns first; if a name differs, fix it **once** in `COL` (or in `SEL`
for §2b) and the rest of the notebook follows. Nothing here is hard-wired to a
column name except through these maps.
""")

C(r"""
for key in SECTION_FILES:
    df = data["CME_water"][key]
    print(f"\n=== {key} (CME_water) ===")
    print("columns:", list(df.columns))
    print(df.head(3).to_string(index=False))

# central column map — reconciled to the committed CSV schema
COL = {
    # summary scalars
    "sasa_mean": "sasa_pooled_nm2", "sasa_ci95": "sasa_ci95",
    "rg_mean":   "rg_pooled_nm",    "rg_ci95":   "rg_ci95",
    # coordination CSVs
    "site":      "reference",       "partner":   "partner",
    "cutoff_nm": "cutoff_nm",       "coord":     "pooled_coord",
    "coord_ci95":"ci95",            "resolved":  "structured_shell",
    "per_atom":  "per_atom",        "neff":      "N_eff_pooled",
    # rdf cutoffs
    "first_min": "first_min_nm",
    # per-start hbond summary
    "start":     "start",           "episodes":  "n_episodes",  "bonds": "n_bonds",
    # hbond lifetimes (system level)
    "distinct_bonds": "n_distinct_bonds", "n_events": "n_events",
}

# cross-check anchors (committed exec summary) — QC tolerance only, NOT the source
REF_SASA = {  # (mean, ci95) nm^2
    "GGE_reline":(2.806,0.066), "GGE_glyceline":(2.612,0.064), "GGE_water":(2.590,0.066),
    "CME_water":(3.151,0.084),  "CME_reline":(3.399,0.079),    "CME_glyceline":(3.240,0.060),
    "YIY_water":(5.277,0.041),  "YIY_reline":(5.412,0.061),    "YIY_glyceline":(5.279,0.038),
}
""")


# ─────────────────────────────────────────────────────────────────────────────
# §1 — E3.1 SASA significance
# ─────────────────────────────────────────────────────────────────────────────
M(r"""
## §1 — E3.1: SASA significance (tier-aware, within-peptide row)

The comparison unit is DES-vs-water **within** each peptide; cross-peptide SASA
reflects size, not solvent. For each row we report ΔSASA (DES − water) with a
combined ACF-corrected 95% CI (SE-propagated from the per-system corrected CIs),
and flag whether that interval excludes zero — the corrected significance test.

Tier guards, applied not smoothed:
- **D** (GGE_water baseline; CME_glyceline): pooled mean is a diverse-start
  *ensemble*, not an equilibrium value — a Δ against it carries that caveat.
- **C** (CME_reline): entered both ways (full pool and consensus, with the
  +0.092 nm² inclusion shift).
- **A/B**: pooled mean is reportable directly.

`GGE_glyceline` is the load-bearing case (carried Phase-3 effect size 1.18); its
water baseline is tier-D and brackets the glyceline mean, so §1.2 tests it
against the pooled water mean *and* each separated water basin.
""")

C(r"""
def get_scalar(system, mean_key, ci_key):
    df = data[system]["summary"]
    row = df.iloc[0]
    return float(row[COL[mean_key]]), float(row[COL[ci_key]])

# QC: recomputed scalars should track the committed anchors within tolerance
sasa = {}
for s in SYSTEMS:
    m, ci = get_scalar(s, "sasa_mean", "sasa_ci95")
    sasa[s] = (m, ci)
    rm, rci = REF_SASA[s]
    if abs(m - rm) > 0.02:
        print(f"  [check] {s}: SASA {m:.3f} vs committed {rm:.3f} — provenance drift?")

sasa_tbl = pd.DataFrame([
    {"system": s, "peptide": s.split("_")[0], "solvent": s.split("_",1)[1],
     "tier": TIER[s], "SASA_nm2": sasa[s][0], "SASA_ci95": sasa[s][1],
     "se": sasa[s][1]/1.96}
    for s in SYSTEMS
])
print(sasa_tbl.to_string(index=False))
""")

C(r"""
# ΔSASA (DES - water) per peptide row, with SE-propagated corrected CI
def delta_row(pep, des):
    w = sasa_tbl.query("peptide==@pep and solvent=='water'").iloc[0]
    d = sasa_tbl.query("peptide==@pep and solvent==@des").iloc[0]
    delta = d.SASA_nm2 - w.SASA_nm2
    se    = np.hypot(d.se, w.se)
    ci    = 1.96 * se
    return {
        "peptide": pep, "contrast": f"{des} - water",
        "dSASA_nm2": round(delta, 3), "ci95": round(ci, 3),
        "lo": round(delta-ci, 3), "hi": round(delta+ci, 3),
        "excludes_zero": bool(abs(delta) > ci),
        "water_tier": TIER[f"{pep}_water"], "des_tier": TIER[f"{pep}_{des}"],
        "note": ("water baseline tier-D (ensemble, not equilibrium)"
                 if TIER[f"{pep}_water"]=="D" else
                 "des tier-D (ensemble, not equilibrium)"
                 if TIER[f"{pep}_{des}"]=="D" else ""),
    }

e31 = pd.DataFrame([delta_row(p, d) for p in PEPTIDES for d in ("reline","glyceline")])
e31.to_csv(os.path.join(OUTDIR, "e3.1_sasa_effects.csv"), index=False)
print(e31.to_string(index=False))
""")

C(r"""
# §1.2 — load-bearing GGE_glyceline against the tier-D water ensemble and its basins
gly_m, gly_ci = sasa["GGE_glyceline"]
basins = data["GGE_water"]["tier_sasa"]   # separated SASA states (high/low)
print(f"GGE_glyceline SASA = {gly_m:.3f} +/- {gly_ci:.3f} nm^2 (tier A)\n")

if basins is not None:
    print("GGE_water separated states (tier D — NOT an equilibrium mean):")
    print(basins.to_string(index=False))
    print("\nObservation: glyceline sits relative to each water basin as printed;")
    print("the pooled water Delta is reported in e3.1 but is against an ensemble baseline.")
else:
    print("[check] GGE_water basin file not found — set TIER_FILES/location in §0.1")
""")

C(r"""
# figure: SASA by system, grouped by peptide, tier-annotated, corrected CIs
fig, ax = plt.subplots(figsize=(8.5, 4.2))
col = {"water":"#2c7bb6","reline":"#d7191c","glyceline":"#1a9641"}
x = np.arange(len(PEPTIDES)); w = 0.26
for k, solv in enumerate(SOLVENTS):
    ys  = [sasa[f"{p}_{solv}"][0]  for p in PEPTIDES]
    es  = [sasa[f"{p}_{solv}"][1]  for p in PEPTIDES]
    trs = [TIER[f"{p}_{solv}"]     for p in PEPTIDES]
    bars = ax.bar(x + (k-1)*w, ys, w, yerr=es, capsize=3,
                  color=col[solv], label=solv, alpha=0.9)
    for b, tr in zip(bars, trs):
        ax.text(b.get_x()+b.get_width()/2, b.get_height(), tr,
                ha="center", va="bottom", fontsize=8, color="0.3")
ax.set_xticks(x); ax.set_xticklabels(PEPTIDES)
ax.set_ylabel("motif SASA (nm$^2$)")
ax.set_title("E3.1 — pooled SASA, ACF-corrected CIs (tier letters annotated)")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "e3.1_sasa_by_system.png"), dpi=200)
plt.show()
""")


# ─────────────────────────────────────────────────────────────────────────────
# §2a — E3.2 qualitative matrix
# ─────────────────────────────────────────────────────────────────────────────
M(r"""
## §2a — E3.2: qualitative present/absent matrix (cutoff-free)

The direct-vs-indirect adjudication needs no recomputation: against the water
null (no DES side-chain shells), a resolved DES shell is direct contact; SASA
moving without one is consistent with the indirect restructuring hypothesis.
Resolved shells are read via the `structured_shell` flag (reportable nulls have
`structured_shell == False` and a NaN cutoff, and are correctly excluded).
""")

C(r"""
def resolved_shells(system, which):
    df = data[system][which]
    if df is None or len(df) == 0:
        return []
    out = []
    for _, r in df.iterrows():
        ok = True
        if COL["resolved"] in df.columns:
            ok = bool(r[COL["resolved"]])
        elif COL["cutoff_nm"] in df.columns:
            ok = pd.notna(r[COL["cutoff_nm"]])
        if ok:
            out.append((str(r[COL["site"]]), str(r[COL["partner"]]),
                        r.get(COL["cutoff_nm"], np.nan),
                        r.get(COL["coord"], np.nan),
                        r.get(COL["per_atom"], np.nan)))
    return out

rows = []
for s in SYSTEMS:
    pep, solv = s.split("_")[0], s.split("_",1)[1]
    for which in ("backbone", "sidechain"):
        for site, partner, cut, coord, per_atom in resolved_shells(s, which):
            rows.append({"peptide": pep, "solvent": solv, "scope": which,
                         "site": site, "partner": partner, "cutoff_nm": cut,
                         "coordination": coord, "per_atom": per_atom})
long = pd.DataFrame(rows)

# present/absent matrix: partners coordinating each (peptide, scope, site) by solvent
def _agg(g):
    return ", ".join(sorted(set(g["partner"]))) if len(g) else "—"
matrix = (long.groupby(["peptide","scope","site","solvent"])
              .apply(_agg).unstack("solvent").reindex(columns=SOLVENTS).fillna("—"))
matrix.to_csv(os.path.join(OUTDIR, "e3.2_coordination_matrix.csv"))
print(matrix.to_string())
""")

M(r"""
Read this against the water column: entries there should be water/none. A DES
partner appearing in the reline/glyceline column where water shows `—` is the
direct-contact signal; its **magnitude** across solvents is deferred to §2b at a
matched cutoff. Qualitative observations already legible here (urea–backbone
present GGE/YIY absent CME; glycerol side-chain at CME Cys/Glu and YIY Tyr;
CME sulfur solvation cosolvent-specific) are recorded as such, framed as
hypotheses for the synthesis.
""")


# ─────────────────────────────────────────────────────────────────────────────
# §3 — E3.3 YIY H-bonds
# ─────────────────────────────────────────────────────────────────────────────
M(r"""
## §3 — E3.3: YIY backbone H-bonds (per-start; direction robustness)

The metric is already matched across the YIY row (1 ps, backbone N→O, continuous
episodes), so this is a per-start reduction. Episodes sum across starts; the
distinct-bond count ("turns") does **not** (a bond recurs across starts), so
turns is taken system-level from `hbond_lifetimes` (`n_distinct_bonds`) while
per-start `n_bonds` is used only to flag the zero-bond start. The test: does
"reline is the low-H-bond arm" survive removing that zero-bond start?
""")

C(r"""
YIY = [f"YIY_{s}" for s in SOLVENTS]

per_start = []
for s in YIY:
    df = data[s]["hbond_perstart"].copy()
    df["system"] = s; df["solvent"] = s.split("_",1)[1]
    per_start.append(df)
per_start = pd.concat(per_start, ignore_index=True)
keep = [COL["start"], COL["episodes"], COL["bonds"], "system", "solvent"]
per_start = per_start[[c for c in keep if c in per_start.columns]]
print(per_start.to_string(index=False))

# system-level distinct backbone bonds ("turns") from hbond_lifetimes
life = pd.concat([data[s]["hbond_life"].assign(solvent=s.split("_",1)[1]) for s in YIY],
                 ignore_index=True)
distinct = life.set_index("solvent")[COL["distinct_bonds"]].reindex(SOLVENTS)

def ep_totals(df):
    return df.groupby("solvent")[COL["episodes"]].sum().reindex(SOLVENTS)

full_ep = ep_totals(per_start)
noz_ep  = ep_totals(per_start[per_start[COL["bonds"]] > 0])

row = pd.DataFrame({
    "episodes_all":           full_ep,
    "episodes_excl_zerobond": noz_ep,
    "distinct_bonds":         distinct,
    "n_starts":       per_start.groupby("solvent")[COL["start"]].nunique().reindex(SOLVENTS),
    "zerobond_starts": (per_start[per_start[COL["bonds"]]==0]
                        .groupby("solvent")[COL["start"]].nunique()
                        .reindex(SOLVENTS).fillna(0).astype(int)),
})
print("\n", row.to_string())

order_full = list(full_ep.sort_values().index)
order_noz  = list(noz_ep.sort_values().index)
print(f"\nepisode ordering (all)      : {' < '.join(order_full)}")
print(f"episode ordering (excl. 0)  : {' < '.join(order_noz)}")
print(f"reline is low-arm (all)     : {order_full[0]=='reline'}")
print(f"reline is low-arm (excl. 0) : {order_noz[0]=='reline'}")

per_start.to_csv(os.path.join(OUTDIR, "e3.3_yiy_hbond_per_start.csv"), index=False)
row.to_csv(os.path.join(OUTDIR, "e3.3_yiy_row_totals.csv"))
""")

C(r"""
# figure: per-start episodes by solvent, distinct-bond count annotated per bar
fig, ax = plt.subplots(figsize=(8, 4))
col = {"water":"#2c7bb6","reline":"#d7191c","glyceline":"#1a9641"}
xs = np.arange(len(STARTS)); w = 0.26
for k, solv in enumerate(SOLVENTS):
    sub = per_start[per_start.solvent==solv].set_index(COL["start"]).reindex(STARTS)
    ep  = sub[COL["episodes"]].values
    nb  = sub[COL["bonds"]].values
    bars = ax.bar(xs + (k-1)*w, ep, w, color=col[solv], label=solv, alpha=0.9)
    for b, n in zip(bars, nb):
        if np.isfinite(b.get_height()):
            ax.text(b.get_x()+b.get_width()/2, b.get_height(),
                    f"{int(n) if np.isfinite(n) else 0}b",
                    ha="center", va="bottom", fontsize=7, color="0.3")
ax.set_xticks(xs); ax.set_xticklabels(STARTS)
ax.set_ylabel("backbone N->O episodes"); ax.set_xlabel("diverse start")
ax.set_title("E3.3 — YIY backbone H-bonds per start (distinct bonds annotated)")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "e3.3_yiy_hbond_per_start.png"), dpi=200)
plt.show()
""")


# ─────────────────────────────────────────────────────────────────────────────
# §4 — synthesis scaffold
# ─────────────────────────────────────────────────────────────────────────────
M(r"""
## §4 — Synthesis scaffold (deterministic from the section outputs)

Assembles the cross-solvent write-up from the §1–§3 tables (and §2b if cached),
as observations/hypotheses with `[REF]` placeholders. Values are inserted from
the computed frames so the document stays a function of the outputs, not authored
by hand. Prose stays hedged: hypotheses supported by the data, tier-aware, no
generalisation beyond the matrix. Edit the interpretive lines in place; re-running
regenerates the numeric slots only.
""")

C(r"""
matched_path = os.path.join(OUTDIR, "e3.2_matched_cutoff.csv")
matched_note = ("present (see §2b)" if os.path.exists(matched_path)
                else "PENDING — run §2b for matched-cutoff magnitudes")

lines = []
lines.append("# E3.x Cross-Solvent Synthesis\n")
lines.append("Observations and hypotheses read across the nine per-system E3.x "
             "summaries at matched cutoffs. Supported by the committed section "
             "outputs; not conclusions of the phase. External refs as [REF].\n")

lines.append("\n## E3.1 — SASA significance\n")
for _, r in e31.iterrows():
    verdict = "excludes zero" if r.excludes_zero else "includes zero"
    tail = f" ({r.note})" if r.note else ""
    lines.append(f"- {r.peptide} {r.contrast}: dSASA = {r.dSASA_nm2:+.3f} "
                 f"[{r.lo:+.3f}, {r.hi:+.3f}] nm^2, corrected 95% CI {verdict}{tail}.")
lines.append("\nLoad-bearing: GGE_glyceline is tested against a tier-D water "
             "baseline that brackets its mean; the corrected verdict is reported "
             "against the pooled ensemble and each separated basin (§1.2). [REF]\n")

lines.append("\n## E3.2 — direct contact vs indirect restructuring\n")
lines.append(f"- Matched-cutoff magnitudes: {matched_note}.")
lines.append("- Qualitative matrix (§2a): DES side-chain shells are absent in the "
             "water baseline; where a DES partner is resolved in reline/glyceline "
             "the contact is direct, otherwise the SASA shift is consistent with "
             "the indirect restructuring hypothesis. [REF]")

lines.append("\n## E3.3 — YIY backbone H-bond direction\n")
lines.append(f"- Episode ordering across the YIY row: {' < '.join(order_full)} "
             f"(all starts); {' < '.join(order_noz)} excluding zero-bond starts.")
lines.append(f"- Reline as the low-H-bond arm is {'robust to' if order_noz[0]=='reline' else 'sensitive to'} "
             "removing the zero-bond start; this inverts the Phase-3 region-only, "
             "strided framing. [REF]")

synth = "\n".join(lines)
with open(os.path.join(OUTDIR, "E3x_cross_solvent_synthesis.md"), "w") as fh:
    fh.write(synth + "\n")
print(synth)
""")


# ─────────────────────────────────────────────────────────────────────────────
# §2b — E3.2 matched-cutoff coordination (trajectory recompute; guarded, cached)
# ─────────────────────────────────────────────────────────────────────────────
M(r"""
## §2b — E3.2: matched-cutoff coordination (guarded, cached, one pass)

**The only trajectory-touching section.** Per-system coordination used
pair-specific first-shell cutoffs, so magnitudes are not comparable across
solvents. This recomputes each resolved (reference, partner) at a **common**
cutoff per pair-type — derived from the committed `rdf_cutoffs.csv`
(`first_min_nm`), not invented — using the same estimator as the per-system
phase (total count of partner atoms within the cutoff, ACF-corrected CI).

It additionally forces `(site, water_O)` across **all three solvents** for every
site that resolves a water shell in any solvent. That is the sharper E3.2 signal:
the indirect-restructuring hypothesis predicts hydration at a site changes across
solvents even where no DES shell forms, which is exactly a matched-cutoff
water-coordination comparison — not visible in the qualitative matrix.

Design, per the calibration (5.1 s/load, 899 waters, 0.8 s/water-RDF → full pass
≈ 5–10 min): **one system at a time, load its four starts once, compute every
matched-cutoff pair for that system in the same pass, write, free.** Loads
dominate and are fixed, so recomputing all pairs (plus the already-common-cutoff
pairs as QC) is the efficient choice. Guarded by the cache so it runs once.

**Verify before first run** (the four provenance items):
1. `SEL` atom-selection strings and reference/partner **labels** against a real
   reline/glyceline coordination CSV — the water system confirms `backbone_O`,
   `motif_Cys/Met/Glu`; `motif_Tyr` and the DES partner labels (`urea_O`,
   `urea_N`, `glycerol_O`, `chloride`) are inferred.
2. `MATCHED_RULE` — how the common cutoff is chosen per pair-type (default `max`,
   capturing the full first shell across solvents).
3. The `convergence_stats` import path and `compute_stat` signature.
4. That the per-system coordination used this same traj root / stride / window
   (the calibration confirmed format and window match).
""")

C(r"""
CACHE = os.path.join(OUTDIR, "e3.2_matched_cutoff.csv")

if os.path.exists(CACHE):
    matched = pd.read_csv(CACHE)
    print(f"loaded cached matched-cutoff coordination ({len(matched)} rows) — skipping recompute")
else:
    import mdtraj as _md
    sys.path.insert(0, CONV)
    try:
        import convergence_stats as cs   # E2.2 module: compute_stat(series) -> mean, ci95, ...
        HAVE_CS = True
    except Exception as e:
        HAVE_CS = False
        print(f"[check] convergence_stats import failed ({e}); falling back to raw CI (understates)")

    # ---- selections keyed by CSV reference/partner labels (verify vs prmtop) --
    SEL = {
        "backbone_O": "protein and name O",
        "motif_Tyr":  "resname TYR and name OH",
        "motif_Cys":  "resname CYS CYX and name SG",
        "motif_Met":  "resname MET and name SD",
        "motif_Glu":  "resname GLU and name OE1 OE2",
        "water_O":    "water and name O",
        "urea_O":     "resname URE and name O",
        "urea_N":     "resname URE and name N1 N2",
        "glycerol_O": "resname GOL and name O1 O2 O3",
        "chloride":   "resname CLA",
    }

    # ---- matched cutoff per (reference, partner) from committed first_min_nm --
    MATCHED_RULE = "max"   # {"max","mean","min"} — max captures the full first shell everywhere
    cut_rows = []
    for s in SYSTEMS:
        rc = data[s]["rdf_cutoffs"]
        if rc is None: continue
        for _, r in rc.iterrows():
            v = r[COL["first_min"]]
            if pd.isna(v): continue
            cut_rows.append({"site": str(r[COL["site"]]), "partner": str(r[COL["partner"]]),
                             "solvent": s.split("_",1)[1], "cutoff_nm": float(v)})
    cuts = pd.DataFrame(cut_rows)
    matched_cut = (cuts.groupby(["site","partner"])["cutoff_nm"].agg(MATCHED_RULE)
                        .rename("matched_cutoff_nm").reset_index())
    print("matched cutoffs per pair-type:")
    print(matched_cut.to_string(index=False))

    # ---- recompute set: every resolved (site,partner) in the matrix ...
    want = matched_cut.merge(
        long[["peptide","solvent","site","partner"]].drop_duplicates(),
        on=["site","partner"], how="inner")

    # ... plus (site, water_O) across ALL solvents, per peptide where water resolves
    aug = []
    for p in PEPTIDES:
        sites_p = long.loc[(long.peptide==p) & (long.partner=="water_O"), "site"].unique()
        for st in sites_p:
            for sv in SOLVENTS:
                aug.append({"site": st, "partner": "water_O", "peptide": p, "solvent": sv})
    aug = pd.DataFrame(aug).merge(matched_cut, on=["site","partner"], how="left")
    want = (pd.concat([want, aug], ignore_index=True)
              .dropna(subset=["matched_cutoff_nm"])
              .drop_duplicates(["peptide","solvent","site","partner"]))

    def _count_within(traj, sel_a, sel_b, cutoff, chunk=500):
        top = traj.topology
        a = top.select(sel_a); b = top.select(sel_b)
        if len(a)==0 or len(b)==0:
            return None, 0
        pairs = np.array([[i, j] for i in a for j in b])
        counts = []
        for k in range(0, traj.n_frames, chunk):
            d = _md.compute_distances(traj[k:k+chunk], pairs)
            counts.extend((d < cutoff).sum(axis=1).tolist())
        return np.asarray(counts, float), len(a)   # total per frame (matches pooled_coord); n_ref

    out_rows = []
    for s in SYSTEMS:
        pep, solv = s.split("_")[0], s.split("_",1)[1]
        pairs_here = want[(want.peptide==pep) & (want.solvent==solv)]
        if pairs_here.empty:
            continue
        t0 = time.perf_counter()
        for start in STARTS:
            p = traj_path(s, start)
            if not os.path.isfile(p):
                print(f"  [skip] {s} {start}: {p} not found"); continue
            tr = _md.load(p, top=top_path(s), stride=STRIDE)[WIN_START_FRAME:]
            for _, pr in pairs_here.iterrows():
                sa = SEL.get(pr.site); sb = SEL.get(pr.partner)
                if sa is None or sb is None:
                    print(f"  [check] no SEL for {pr.site}/{pr.partner}"); continue
                series, n_ref = _count_within(tr, sa, sb, pr.matched_cutoff_nm)
                if series is None:
                    continue
                if HAVE_CS:
                    m, ci = cs.compute_stat(series)[:2]
                else:
                    m, ci = float(series.mean()), 1.96*series.std()/np.sqrt(len(series))
                out_rows.append({"system": s, "peptide": pep, "solvent": solv,
                                 "site": pr.site, "partner": pr.partner,
                                 "matched_cutoff_nm": pr.matched_cutoff_nm,
                                 "start": start, "coord": m, "per_atom": m/max(n_ref,1),
                                 "ci95": ci})
        print(f"  {s}: {len(pairs_here)} pairs x starts in {time.perf_counter()-t0:.1f}s")

    per_start_m = pd.DataFrame(out_rows)
    # pool across starts: mean of per-start means; CI combined in quadrature / k
    matched = (per_start_m.groupby(["system","peptide","solvent","site","partner","matched_cutoff_nm"])
               .agg(coord=("coord","mean"), per_atom=("per_atom","mean"),
                    ci95=("ci95", lambda v: float(np.sqrt(np.sum(np.square(v)))/len(v))),
                    n_starts=("start","nunique"))
               .reset_index())
    matched.to_csv(CACHE, index=False)
    print(f"\nwrote {CACHE} ({len(matched)} pooled rows)")

print(matched.head(20).to_string(index=False))
""")

C(r"""
# QC: pairs already at a common cutoff should reproduce the stored pooled_coord
if "matched" in dir() and len(matched):
    qc = matched.merge(long, on=["peptide","solvent","site","partner"],
                       how="inner", suffixes=("_matched","_stored"))
    if len(qc):
        qc["abs_diff"] = (qc["coord"] - qc["coordination"]).abs()
        print("recompute vs stored (common-cutoff pairs are the meaningful QC):")
        print(qc[["system","site","partner","matched_cutoff_nm","coord","coordination","abs_diff"]]
              .sort_values("abs_diff", ascending=False).head(15).to_string(index=False))
        print("\nLarge diffs on differing-cutoff pairs are expected (that is the point);")
        print("large diffs on already-common-cutoff pairs flag estimator drift — inspect before synthesis.")
""")

M(r"""
### Fold §2b into the synthesis

With the cache written, re-run **§4** to insert the matched-cutoff magnitudes:
the pending line resolves to the direct-vs-indirect magnitude comparison per
observable — including the water-coordination-across-solvents test that is the
core indirect-restructuring signal. Every subsequent notebook run reads the
cache and skips the recompute.
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
