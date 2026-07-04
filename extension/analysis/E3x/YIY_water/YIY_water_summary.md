# E3.x system summary - YIY_water

Tier B | 4-start (compact, extended, mid, open) | robust SASA / backbone converged (mildest YIY)
Peptide YIY | solvent water

## E3.1 - Replicate-aware SASA & backbone coordination
- SASA (carried E2.2 pooled): 5.2769 +/- 0.0407 nm2  [robust]
- Rg   (carried E2.2 pooled): 0.6374 +/- 0.0125 nm  [robust]
- Convergence character: robust SASA / backbone converged (mildest YIY)
- Backbone-O coordination (pair-specific cutoffs, ACF-corrected CIs):
```
partner cutoff_nm pooled_coord   ci95
water_O     0.338       5.4004 0.0470
```
- Replicate significance for THIS system: achieved (both observables robust). The named E3.1 test case (GGE_glyceline, corrected ES 1.18) is a separate system; its significance verdict is deferred there.

## E3.2 - Side-chain coordination (indirect-restructuring test)
- Side-chain coordination (pair-specific cutoffs, ACF-corrected CIs):
```
reference partner cutoff_nm pooled_coord   ci95
motif_Tyr water_O     0.358       6.6924 0.0159
```
- CONTRIBUTES (baseline): water carries no glycerol, so it cannot test the glycerol-sidechain question directly. It supplies aqueous baseline that rules out a peptide-intrinsic explanation. The decisive glycerol comparison is deferred to the glyceline systems + water baseline.

## E3.3 - Backbone hydrogen bonds (anomaly test)
- 4 backbone N->O turns; 1 ps resolution, 64585 continuous episodes, median 2.0 ps (break_tol=0).
- Per-start (episodes/distinct-bonds): compact:34007ep/3b; extended:9761ep/2b; mid:11074ep/2b; open:9743ep/2b. A start-dependent backbone-H-bond pattern shows here (relevant to the YIY reline anomaly, which is start-specific).
- CONTRIBUTES (baseline): YIY, but this solvent is the anomaly-control arm. It supplies the water backbone-H-bond baseline (4 turns, 64585 episodes at 1 ps) against which YIY_reline is adjudicated. Note: at 1 ps backbone-only scope this is NOT zero, unlike the Phase-3 region-only metric - the anomaly is a cross-solvent difference, not presence/absence.
- Note: the continuous-episode median is largely blind to slow exchange; the persistence definition (continuous vs intermittent/ACF) is fixed at the cross-solvent stage.

## Carried caveats
- Chloride: no structured shell (backbone or side-chain) - reportable null.
- Side-chain cutoffs are novel (no E3.0 anchor); low reference-atom count -> convergence caveat.

## Question ledger
CLOSED by this system:
  - E3.1 convergence for this system (robust)
CONTRIBUTES toward (deferred to cross-solvent stage):
  - E3.1 headline significance (GGE_glyceline ES 1.18)
  - E3.2 indirect-mechanism adjudication (glyceline + water)
  - E3.3 anomaly adjudication (needs YIY_reline)
