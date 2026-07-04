# E3.x system summary - CME_reline

Tier C | 4-start (compact, extended, mid, open) | one-trap (compact/mid/open consensus + metastable extended)
Peptide CME | solvent reline

## E3.1 - Replicate-aware SASA & backbone coordination
- SASA one-trap: one start (loo-flagged) sits apart from the consensus; reported both ways.
  - full pool (all starts): 3.399 +/- 0.0792 nm2
  - consensus (drop extended): 3.3073 +/- 0.0703 nm2
  - inclusion shift (full - consensus): +0.0917 nm2
  - sidecar 3.399 (delta 0.0)
- Rg (carried): 0.7335 +/- 0.0227 nm  [robust]
- Backbone-O coordination (pooled):
```
  partner cutoff_nm pooled_coord   ci95
  water_O     0.333       2.6493 0.1675
choline_N     0.407       0.1508 0.0187
 chloride  no shell       0.0000 0.0000
   urea_N  no shell       0.0000 0.0000
   urea_O  no shell       0.0000 0.0000
```

## E3.2 - Side-chain coordination
- Side-chain coordination (pair-specific cutoffs, ACF-corrected CIs):
```
reference   partner cutoff_nm pooled_coord   ci95
motif_Cys   water_O  no shell       0.0000 0.0000
motif_Cys choline_N     0.433       0.0739 0.0052
motif_Cys  chloride  no shell       0.0000 0.0000
motif_Cys    urea_N     0.447       0.4379 0.0193
motif_Cys    urea_O     0.377       0.0937 0.0053
motif_Met   water_O  no shell       0.0000 0.0000
motif_Met choline_N     0.447       0.1245 0.0074
motif_Met  chloride  no shell       0.0000 0.0000
motif_Met    urea_N     0.433       0.3793 0.0157
motif_Met    urea_O  no shell       0.0000 0.0000
motif_Glu   water_O     0.327       5.2250 0.0853
motif_Glu choline_N     0.527       0.4816 0.0283
motif_Glu  chloride  no shell       0.0000 0.0000
motif_Glu    urea_N     0.392       0.4192 0.0209
motif_Glu    urea_O  no shell       0.0000 0.0000
```
- reline DES contacts (choline backbone shoulder excluded): no direct backbone shell; side-chain shell at motif_Cys-choline_N, motif_Cys-urea_N, motif_Cys-urea_O, motif_Met-choline_N, motif_Met-urea_N, motif_Glu-choline_N, motif_Glu-urea_N.

## E3.3 - Backbone hydrogen bonds
- 9 backbone N->O turns; 1 ps resolution, 157729 continuous episodes, median 2.0 ps (break_tol=0).
- Per-start (episodes/distinct-bonds): compact:39786ep/4b; extended:14472ep/3b; mid:53091ep/6b; open:50380ep/6b.

## Caveats
- Chloride: no structured shell (backbone or side-chain) - reportable null.
- Choline-N -> backbone-O: weak feature, highest cutoff-sensitivity - shoulder, not an established shell.
- Side-chain cutoffs are novel (no E3.0 anchor); low reference-atom count -> convergence caveat.

## Question ledger
CLOSED by this system:
  - (none - validation/baseline system)
CONTRIBUTES toward (deferred to cross-solvent stage):
  - E3.1 SASA (one-trap; full-pool vs consensus reported)
  - E3.2 indirect-mechanism adjudication (glyceline vs water)
  - E3.3 anomaly adjudication (YIY systems)
