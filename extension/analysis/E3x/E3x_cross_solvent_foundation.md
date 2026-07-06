# Cross-Solvent Analysis — Thread Foundation

Seed document for the cross-solvent stage of the DES-peptide study. Paste or attach at
the start of the new thread. It defines the task, the evidence base, and the constraints;
it does not carry out the analysis.

## Task

Synthesise the nine E3.x per-system summaries into a cross-solvent interpretation. This is
the interpretation deliberately withheld from the per-system files: reading the systems
against one another, at matched cutoffs, to frame the logged observations as hypotheses
supported by the data. The per-system files assert nothing beyond their own tables; this
stage is where the comparative claims are made — and they are made as observations and
hypotheses, not established conclusions.

## Project context

MD study of how two choline-chloride deep eutectic solvents — reline (ChCl:urea) and
glyceline (ChCl:glycerol), both 20 wt% aqueous — modulate the local solvation and backbone
hydrogen bonding of three *Torreya grandis* peptide motifs (GGE, CME, YIY) relative to
water. Design is a 3 × 3 matrix (nine systems), four diverse-start runs each
(200 ns; canonical window 20–200 ns). Repo `rossjamesgibson21-bit/des-peptide-study`,
env `des-peptide`.

## Status entering this thread

The E3.x per-system pipeline is complete and committed (9/9). Evidence base:

- Nine `{SYSTEM}_summary.md` files at `extension/analysis/E3x/{SYSTEM}/`.
- `extension/analysis/E3x/E3x_executive_summary.md` — index and per-system synopsis.
- Per-system section outputs on disk: `coordination/{backbone,sidechain}_coordination.csv`,
  `rdfs/rdf_cutoffs.csv`, `hbonds/hbond_lifetimes.csv`, `hbonds/per_start_hbond_summary.csv`,
  and the E2.2 convergence sidecars at `extension/analysis/convergence/`.

The summaries carry tables, carried E2.2 verdicts, and computed caveats only — generated
through a single shared module, so each is a deterministic function of its section outputs.

## Inputs and I/O

Primary input: the nine summaries plus the executive summary — load these first; they are
the reduced evidence. For any matched-cutoff recomputation, the per-system coordination and
RDF CSVs are the secondary input. Output of this stage is a cross-solvent analysis
(tables/figures per observable, plus the written synthesis).

Compute note, flag before committing: the per-system coordination numbers use
*pair-specific* first-shell cutoffs, so they are not directly comparable in magnitude
across solvents. A valid cross-solvent magnitude comparison requires a common (matched)
cutoff per pair. If the stored RDF g(r) curves suffice, matched-cutoff coordination is a
cheap re-integration; if it requires re-deriving from trajectories, that is costly and
should be estimated first. Qualitative comparisons (shell present/absent, tier, verdict)
need no recomputation.

## Reporting conventions that constrain the numbers

- Coordination: pair-specific RDF-derived first-shell cutoffs; ACF-corrected 95% CIs
  (parametric normal approximation, N_eff = N/τ_int). Absolute magnitudes carry
  cutoff-dependence — hence the matched-cutoff requirement above.
- SASA and Rg: carried from the E2.2 pooled, leave-one-out framework. Reporting depends on
  the convergence tier.
- Backbone H-bonds: full 1 ps resolution, backbone N→O scope, continuous episodes.
  Stride-10 aliases ~87% of episodes; the corrected metric differs from the Phase-3
  region-only, strided metric.
- Statistics are autocorrelation-corrected throughout; corrected values are primary, full
  methodology to supplementary. Raw bootstrap CIs understate uncertainty by ~14–36× here.

## Convergence tiers (govern which comparisons are valid)

- **A — clean/ergodic** (GGE_reline, GGE_glyceline, CME_water): SASA a single pooled mean.
- **B — robust SASA / non-ergodic backbone** (YIY_water, YIY_reline, YIY_glyceline): SASA
  pooled and reportable; Rg may be trapped and carries a non-ergodic caveat.
- **C — one-trap** (CME_reline): SASA reported as full pool and consensus, with the
  trap-inclusion shift.
- **D — multi-basin** (GGE_water 2 states, CME_glyceline 4 states): SASA a diverse-start
  ensemble of separated states, not an equilibrium mean.

Consequences for cross-solvent comparison: do not treat a tier-D ensemble mean as an
equilibrium value; account for tier-C trap inclusion when comparing SASA; SASA
comparisons must respect the tier of each system rather than compare pooled means blindly.

## The three cross-solvent questions

- **E3.1 — SASA significance.** Does the DES-induced SASA change survive corrected CIs
  across the matrix? GGE_glyceline is the named significance test case (carried Phase-3
  effect size 1.18); the headline verdict against corrected CIs is settled here. Report
  effect sizes with corrected CIs, tier-aware.
- **E3.2 — indirect-restructuring hypothesis.** Is the SASA change accompanied by direct
  DES–side-chain contact, or is it consistent with an indirect (solvent-restructuring)
  mechanism? The per-system evidence is mixed: urea makes backbone contact (GGE, YIY, not
  CME); glycerol makes side-chain contact (CME Cys/Glu, YIY Tyr) and a weak YIY backbone
  contact. Adjudicate direct vs indirect against the water baseline (no DES side-chain
  shells there).
- **E3.3 — YIY backbone H-bond anomaly.** At matched 1 ps backbone-only scope, is the
  YIY reline pattern a genuine cross-solvent difference, and in which direction? The
  corrected metric appears to invert the Phase-3 framing (see observations).

## Observations to adjudicate (from the per-system files)

Recorded across the nine systems; the raw material for this stage. Each is an observation
supported by the tables, to be weighed at matched cutoffs and framed as a hypothesis — not
a conclusion of the per-system phase.

- Urea–backbone contact is present in GGE_reline and YIY_reline, absent in CME_reline
  (a CME-specific absence, not a GGE-only presence).
- Glycerol reach differs by peptide: a weak YIY backbone shell (~0.32/O), none for GGE;
  side-chain contacts at CME Cys/Glu and YIY Tyr.
- CME sulfur solvation is cosolvent-specific: no shell in water, urea coordinates both Cys
  and Met in reline, glycerol coordinates Cys only in glyceline.
- The YIY backbone H-bond pattern appears to invert the Phase-3 framing: at matched 1 ps
  backbone-only scope, reline is the low-H-bond arm (45,678 episodes, 3 turns, one start
  with none) relative to water (64,585, 4) and glyceline (62,681, 4).
- SASA–Rg decoupling is a YIY property with solvent-specific traps: SASA robust across the
  row while Rg is trapped, the dominating start differing by solvent (compact in reline,
  extended/mid in glyceline, none in water).

## Working principles

- Observations, not conclusions: interpretations are hypotheses supported by the data.
- No text that diverges from the data; comparative claims trace to the tables.
- Matched cutoffs for any cross-solvent magnitude comparison; qualitative comparisons
  (shell present/absent, tier, verdict) are cutoff-free.
- Corrected statistics primary; effect sizes with corrected CIs; do not overclaim or
  generalise beyond the matrix.
- No external citations; use [REF] placeholders for anything not in an uploaded document.
  Internal cross-references point to the committed E3.x outputs.
- Flag compute costs before committing; do not change working parameters without cause.

## Deferred beyond this stage

- Equilibrium SASA and basin populations for the multi-basin systems (GGE_water,
  CME_glyceline) — enhanced sampling, Paper 2.
- The backbone-H-bond persistence definition (continuous vs intermittent/ACF), if the
  cross-solvent reading requires it.
- E5 Methods rewrite incorporating the convergence framework (reconciliation note already
  committed at `extension/manuscript/E3x_methodology_reconciliation.md`).
