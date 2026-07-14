# E3.x system summary - GGE_glyceline

Tier A | 4-start (compact, extended, mid, open) | clean / ergodic
Peptide GGE | solvent glyceline

## E3.1 - Replicate-aware SASA & backbone coordination
- SASA (carried): 2.6122 +/- 0.0638 nm2  [robust]
- Rg   (carried): 0.7063 +/- 0.0176 nm  [robust]
- Backbone-O coordination (pair-specific cutoffs, ACF-corrected CIs):
```
   partner cutoff_nm pooled_coord   ci95
   water_O     0.333       3.9949 0.1450
 choline_N     0.412       0.2632 0.0235
  chloride  no shell       0.0000 0.0000
glycerol_O  no shell       0.0000 0.0000
```

## E3.2 - Side-chain coordination
- Side-chain coordination (pair-specific cutoffs, ACF-corrected CIs):
```
reference    partner cutoff_nm pooled_coord   ci95
motif_Glu    water_O     0.327       5.5090 0.0600
motif_Glu  choline_N     0.522       0.4435 0.0234
motif_Glu   chloride  no shell       0.0000 0.0000
motif_Glu glycerol_O     0.317       0.6309 0.0433
flank_Glu    water_O     0.327      16.7601 0.0917
flank_Glu  choline_N     0.527       1.3800 0.0452
flank_Glu   chloride  no shell       0.0000 0.0000
flank_Glu glycerol_O     0.317       1.8527 0.0718
```
- Glycerol side-chain first-shell contact at: motif_Glu (0.63/O), flank_Glu (1.85/O). No backbone-glycerol shell.

## E3.3 - Backbone hydrogen bonds
- 8 backbone N->O turns; 1 ps resolution, 112691 continuous episodes, median 2.0 ps (break_tol=0).
- Per-start (episodes/distinct-bonds): compact:19234ep/2b; extended:28315ep/5b; mid:23490ep/3b; open:41652ep/7b.

## Caveats
- Chloride: no structured shell (backbone or side-chain) - reportable null.
- Choline-N -> backbone-O: weak feature, highest cutoff-sensitivity - shoulder, not an established shell.
- Side-chain cutoffs are novel (no E3.0 anchor); low reference-atom count -> convergence caveat.

## Question ledger
CLOSED by this system:
  - E3.1 significance (named test case: GGE_glyceline)
  - E3.2 glycerol contact - side-chain and backbone (this solvent)
CONTRIBUTES toward (deferred to cross-solvent stage):
  - E3.3 anomaly adjudication (YIY systems)
