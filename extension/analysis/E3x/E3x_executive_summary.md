# E3.x — Extended MD Analysis: Executive Summary

Introduction and navigation guide to the nine per-system summaries. This document
orients a reader to the phase, the reporting conventions, and the observations logged
per system. It does not perform the cross-solvent synthesis: interpretation that weighs
the systems against one another is deferred to the cross-solvent stage, and the
statements below are observations supported by the tables, not conclusions established
by them.

## Scope

The E3.x phase characterises how two choline-chloride deep eutectic solvents — reline
(ChCl:urea) and glyceline (ChCl:glycerol), both at 20 wt% aqueous — modulate the local
solvation and backbone hydrogen bonding of three *Torreya grandis* peptide motifs
relative to water. The design is a 3 × 3 matrix (three peptides × three solvents),
giving nine systems, each analysed from four diverse-start production runs
(200 ns each; canonical window 20–200 ns).

The three analysis questions, one per section of every summary:

- **E3.1** — replicate-aware SASA and backbone coordination (convergence and reporting).
- **E3.2** — side-chain coordination (the indirect-restructuring question).
- **E3.3** — backbone hydrogen bonds (the YIY anomaly).

## How to read the per-system summaries

Each `{SYSTEM}_summary.md` reports tables, verdicts carried verbatim from the E2.2
convergence sidecars, and caveats computed from those tables. No interpretation is
authored into the per-system files; all nine are generated through a single shared
module, so a summary is a deterministic function of its persisted section outputs.

Quantities and their basis:

- **Coordination numbers** use pair-specific, RDF-derived first-shell cutoffs
  (not a uniform cutoff), with autocorrelation-corrected 95% confidence intervals
  (parametric normal approximation, N_eff = N/τ_int). A structured shell is reported
  only where the RDF resolves one; flat-RDF pairs are recorded as reportable nulls,
  not omitted.
- **SASA and Rg** are carried from the E2.2 pooled, leave-one-out framework. How SASA is
  reported depends on the convergence tier (below).
- **Backbone hydrogen bonds** are computed at full 1 ps resolution, backbone N→O scope
  only, as continuous episodes. Stride-10 sampling aliases roughly 87% of episodes on
  every system, so the 1 ps resolution is a requirement, not a preference.

## Convergence tiers

The tier sets how SASA is reported, reflecting what the four diverse starts support:

- **Tier A — clean / ergodic.** SASA carried as a single pooled mean.
  Systems: GGE_reline, GGE_glyceline, CME_water.
- **Tier B — robust SASA / non-ergodic backbone.** SASA pooled and reportable; Rg may be
  dominated by a trapped start and carries a non-ergodic caveat where so.
  Systems: YIY_water, YIY_reline, YIY_glyceline.
- **Tier C — one-trap.** One metastable start sits apart from the consensus; SASA reported
  both ways (full pool and consensus), with the inclusion shift stated.
  System: CME_reline.
- **Tier D — multi-basin.** SASA is a diverse-start ensemble of statistically separated
  states (adjacent-gap test against combined CIs), not an equilibrium mean; the true
  basin populations require enhanced sampling (Paper 2).
  Systems: GGE_water (2 states), CME_glyceline (4 states).

## System index

| System | Peptide | Solvent | Tier | Convergence character |
|---|---|---|---|---|
| GGE_reline | GGE | reline | A | clean / ergodic (3-start) |
| GGE_glyceline | GGE | glyceline | A | clean / ergodic; E3.1 significance test case |
| GGE_water | GGE | water | D | multi-basin (2 SASA states) |
| CME_water | CME | water | A | single broad basin (robust) |
| CME_reline | CME | reline | C | one-trap (metastable extended) |
| CME_glyceline | CME | glyceline | D | multi-basin (4 SASA states) |
| YIY_water | YIY | water | B | robust SASA / backbone converged (mildest YIY) |
| YIY_reline | YIY | reline | B | robust SASA / non-ergodic (Rg dominated by compact) |
| YIY_glyceline | YIY | glyceline | B | robust SASA / non-ergodic (Rg dominated by extended/mid) |

Each row links to its file: `GGE_reline/GGE_reline_summary.md`, and so on.

## Per-system synopsis

Anchoring figures are the pooled, ACF-corrected values reported in each file.

**GGE_reline** (A). SASA 2.806 ± 0.066 nm², Rg 0.764 ± 0.019 nm, both robust
(3-start system). Urea makes direct backbone contact: urea-O shell at 0.562 nm
(coordination 1.28), urea-N at 0.423 nm (0.92). Backbone H-bonds: 4 N→O turns,
36,061 episodes.

**GGE_glyceline** (A). The designated E3.1 replicate-significance test case (carried
Phase-3 effect size 1.18); its verdict against corrected CIs is the load-bearing E3.1
result for the campaign. Glycerol resolves no structured backbone shell (the RDF guard
rejects a spurious second-shell minimum); side-chain glycerol contact is at the Glu
carboxyl. Pooled SASA/Rg to be read from the file.

**GGE_water** (D, 2 states). SASA a two-state ensemble — high 2.844 ± 0.062 nm²
(extended, open) and low 2.336 ± 0.088 nm² (compact, mid), pooled 2.590 ± 0.066 nm² —
not an equilibrium mean. Rg 0.704 ± 0.018 nm, dominated. Backbone water shell at
0.333 nm (coordination 4.89). Backbone H-bonds: 8 turns, 146,682 episodes.

**CME_water** (A). SASA 3.151 ± 0.084 nm², Rg 0.696 ± 0.013 nm, both robust. The Glu
side chain is hydrated (5.13); the Cys and Met sulfurs resolve no structured shell in
water (reportable nulls). Backbone H-bonds: 11 turns, 307,222 episodes.

**CME_reline** (C, trap = extended). SASA full pool 3.399 ± 0.079 nm² vs consensus
3.307 ± 0.070 nm² (inclusion shift +0.092 nm²); Rg 0.734 ± 0.023 nm, robust. Urea makes
no structured backbone shell here — a contrast to GGE and YIY. Backbone H-bonds: 9 turns,
157,729 episodes.

**CME_glyceline** (D, 4 states). SASA four separated states (mid 2.883, extended 3.148,
compact 3.311, open 3.619 nm²), pooled 3.240 ± 0.060 nm²; the equilibrium mean and state
populations are unresolved (heaviest sampling caveat). Rg 0.691 ± 0.016 nm, dominated.
Glycerol contacts the Cys (0.38/O) and Glu (0.50/O) side chains but not Met. Backbone
H-bonds: 12 turns, 292,136 episodes.

**YIY_water** (B). SASA 5.277 ± 0.041 nm², Rg 0.637 ± 0.013 nm, both robust — the mildest
YIY backbone (converges at four starts where reline and glyceline are trapped). Backbone
water shell 0.338 nm (5.40); the two Tyr hydroxyls are strongly hydrated (6.69). Backbone
H-bonds: 4 turns, 64,585 episodes.

**YIY_reline** (B). SASA 5.412 ± 0.061 nm², robust; Rg 0.660 ± 0.010 nm, dominated by the
compact start. Urea contacts the YIY backbone (urea-O 0.562 nm, coordination 1.52 — close
to GGE_reline) and the Tyr hydroxyls (urea-N 0.55); a weak Tyr–chloride shell appears
(0.075). Backbone H-bonds: 3 turns, 45,678 episodes, with the mid start forming no
backbone turns (verified, not an artifact).

**YIY_glyceline** (B). SASA 5.279 ± 0.038 nm², robust; Rg 0.640 ± 0.013 nm, dominated by
the extended/mid pair. Glycerol contacts the Tyr hydroxyls strongly (0.538 nm,
coordination 2.34) and the backbone weakly (0.322 nm, 0.32) — a backbone contact GGE
glyceline lacked. Backbone H-bonds: 4 turns, 62,681 episodes, all starts populated.

## Cross-system observations (pending cross-solvent adjudication)

The following are observations recorded across the nine systems. They are the raw
material for the cross-solvent stage, where they will be weighed at matched cutoffs and
framed as hypotheses; they are not conclusions of this phase.

- **Urea–backbone contact is peptide-dependent.** Direct urea–backbone shells are present
  in GGE_reline and YIY_reline but absent in CME_reline — a CME-specific absence rather
  than a GGE-only presence.
- **Glycerol reach differs by peptide.** Glycerol resolves a weak YIY backbone shell
  (0.32/O) but none for GGE; its side-chain contacts are Cys and Glu in CME, and the Tyr
  hydroxyls in YIY.
- **CME sulfur solvation is cosolvent-specific.** The Cys and Met sulfurs resolve no shell
  in water; urea coordinates both in reline; glycerol coordinates Cys only in glyceline.
- **The YIY backbone H-bond pattern inverts the Phase-3 framing.** At matched 1 ps
  backbone-only scope, reline is the low-H-bond arm (45,678 episodes, 3 turns, one start
  with none) relative to water (64,585, 4) and glyceline (62,681, 4). The Phase-3
  "reline high, water/glyceline zero" reading was a region-only, strided artifact;
  the corrected metric reverses the direction.
- **SASA–Rg decoupling is a YIY property with solvent-specific traps.** Across the YIY row
  SASA is robust while Rg is trapped; the dominating start differs by solvent (compact in
  reline, extended/mid in glyceline, none in water).

## Deferred to later stages

- **Cross-solvent synthesis** — the interpretation withheld above, read across all nine
  summaries at matched cutoffs.
- **Equilibrium SASA for the multi-basin systems** (GGE_water, CME_glyceline) — basin
  populations require enhanced sampling (Paper 2).
- **Persistence definition for backbone H-bonds** — continuous-episode vs
  intermittent/ACF is fixed at the cross-solvent stage; the per-system medians are largely
  blind to slow exchange.

## Provenance

Summaries are regenerated from persisted section outputs through a single shared module
(`e3x_summary.py`) via `regenerate_summary_from_outputs.py`; no value is authored by hand.
Statistics are autocorrelation-corrected throughout (raw bootstrap CIs understate
uncertainty by roughly 14–36× at the effective sample sizes here). No external references
are cited; internal cross-references point to the project's own committed outputs.
