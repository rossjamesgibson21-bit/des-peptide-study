# E3.x system summary - YIY_water

Tier B | 4-start (compact, extended, mid, open) | robust SASA / backbone converged (mildest YIY)
Peptide YIY | solvent water

## E3.1 - Replicate-aware SASA & backbone coordination
- SASA (carried): 5.2769 +/- 0.0407 nm2  [robust]
- Rg   (carried): 0.6374 +/- 0.0125 nm  [robust]
- Backbone-O coordination (pair-specific cutoffs, ACF-corrected CIs):
```
partner cutoff_nm pooled_coord   ci95
water_O     0.338       5.4004 0.0470
```

## E3.2 - Side-chain coordination
- Side-chain coordination (pair-specific cutoffs, ACF-corrected CIs):
```
reference partner cutoff_nm pooled_coord   ci95
motif_Tyr water_O     0.358       6.6924 0.0159
```
- Water only: aqueous side-chain baseline (no DES partners).

## E3.3 - Backbone hydrogen bonds
- 4 backbone N->O turns; 1 ps resolution, 64585 continuous episodes, median 2.0 ps (break_tol=0).
- Per-start (episodes/distinct-bonds): compact:34007ep/3b; extended:9761ep/2b; mid:11074ep/2b; open:9743ep/2b.
- Anomaly-control arm: water backbone-H-bond baseline for the YIY_reline comparison (1 ps backbone-only scope; not zero at this scope, unlike the Phase-3 region-only metric).

## Caveats
- Side-chain cutoffs are novel (no E3.0 anchor); low reference-atom count -> convergence caveat.

## Question ledger
CLOSED by this system:
  - E3.1 convergence for this system (both robust)
CONTRIBUTES toward (deferred to cross-solvent stage):
  - E3.1 headline significance (GGE_glyceline)
  - E3.2 indirect-mechanism adjudication (glyceline vs water)
  - E3.3 anomaly adjudication (needs YIY_reline)
