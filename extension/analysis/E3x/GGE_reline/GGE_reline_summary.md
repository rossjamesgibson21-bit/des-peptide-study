# E3.x system summary - GGE_reline

Tier A | 3-start (compact, extended, mid) | clean / ergodic
Peptide GGE | solvent reline

## E3.1 - Replicate-aware SASA & backbone coordination
- SASA (carried): 2.8056 +/- 0.0662 nm2  [robust]
- Rg   (carried): 0.7636 +/- 0.0189 nm  [robust]
- Backbone-O coordination (pair-specific cutoffs, ACF-corrected CIs):
```
  partner cutoff_nm pooled_coord   ci95
  water_O     0.333       4.6432 0.0976
choline_N     0.412       0.3969 0.0349
 chloride  no shell       0.0000 0.0000
   urea_N     0.423       0.9229 0.0454
   urea_O     0.562       1.2849 0.0565
```

## E3.2 - Side-chain coordination
- Side-chain coordination (pair-specific cutoffs, ACF-corrected CIs):
```
reference   partner cutoff_nm pooled_coord   ci95
motif_Glu   water_O     0.322       5.7587 0.0437
motif_Glu choline_N     0.548       0.6285 0.0351
motif_Glu  chloride  no shell       0.0000 0.0000
motif_Glu    urea_N     0.352       0.3087 0.0188
motif_Glu    urea_O  no shell       0.0000 0.0000
flank_Glu   water_O     0.327      17.9749 0.0530
flank_Glu choline_N     0.522       1.6996 0.0506
flank_Glu  chloride  no shell       0.0000 0.0000
flank_Glu    urea_N     0.392       1.3552 0.0437
flank_Glu    urea_O  no shell       0.0000 0.0000
```
- reline DES contacts (choline backbone shoulder excluded): backbone shell at urea_N, urea_O; side-chain shell at motif_Glu-choline_N, motif_Glu-urea_N, flank_Glu-choline_N, flank_Glu-urea_N.

## E3.3 - Backbone hydrogen bonds
- 4 backbone N->O turns; 1 ps resolution, 36061 continuous episodes, median 2.0 ps (break_tol=0).
- Per-start (episodes/distinct-bonds): compact:9124ep/1b; extended:7321ep/2b; mid:19616ep/2b.

## Caveats
- Chloride: no structured shell (backbone or side-chain) - reportable null.
- Choline-N -> backbone-O: weak feature, highest cutoff-sensitivity - shoulder, not an established shell.
- Side-chain cutoffs are novel (no E3.0 anchor); low reference-atom count -> convergence caveat.

## Question ledger
CLOSED by this system:
  - E3.1 convergence for this system (both robust)
CONTRIBUTES toward (deferred to cross-solvent stage):
  - E3.1 headline significance (GGE_glyceline)
  - E3.2 indirect-mechanism adjudication (glyceline vs water)
  - E3.3 anomaly adjudication (YIY systems)
