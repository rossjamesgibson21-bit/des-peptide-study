# E3.x system summary - CME_glyceline

Tier D | 4-start (compact, extended, mid, open) | multi-basin (4 separated SASA states; heaviest caveat)
Peptide CME | solvent glyceline

## E3.1 - Replicate-aware SASA & backbone coordination
- SASA: 4 statistically separated state(s) across the diverse starts (adjacent gap > combined 95% CI); diverse-start ensemble, not an equilibrium mean.
- Per-start SASA means (nm2):
```
  mid:2.8831(ci0.1014)
  extended:3.1479(ci0.1206)
  compact:3.3110(ci0.0667)
  open:3.6189(ci0.1020)
```
- Grouped states: mid:2.8831+/-0.1014 | extended:3.1479+/-0.1206 | compact:3.3110+/-0.0667 | open:3.6189+/-0.1020
- Pooled ensemble: 3.2402 +/- 0.0602 nm2 (sidecar 3.2402, delta 0.0)
- Rg (carried): 0.691 +/- 0.0155 nm  [dominated]
- Backbone-O coordination (pooled):
```
   partner cutoff_nm pooled_coord   ci95
   water_O     0.333       2.6139 0.1024
 choline_N     0.412       0.1519 0.0155
  chloride  no shell       0.0000 0.0000
glycerol_O  no shell       0.0000 0.0000
```

## E3.2 - Side-chain coordination
- Side-chain coordination (pair-specific cutoffs, ACF-corrected CIs):
```
reference    partner cutoff_nm pooled_coord   ci95
motif_Cys    water_O  no shell       0.0000 0.0000
motif_Cys  choline_N     0.442       0.0862 0.0060
motif_Cys   chloride  no shell       0.0000 0.0000
motif_Cys glycerol_O     0.412       0.3841 0.0186
motif_Met    water_O  no shell       0.0000 0.0000
motif_Met  choline_N     0.442       0.0804 0.0062
motif_Met   chloride  no shell       0.0000 0.0000
motif_Met glycerol_O  no shell       0.0000 0.0000
motif_Glu    water_O     0.327       4.2553 0.1104
motif_Glu  choline_N     0.498       0.2729 0.0198
motif_Glu   chloride  no shell       0.0000 0.0000
motif_Glu glycerol_O     0.317       0.4986 0.0414
```
- Glycerol side-chain first-shell contact at: motif_Cys (0.38/O), motif_Glu (0.50/O); no side-chain shell at motif_Met. No backbone-glycerol shell.

## E3.3 - Backbone hydrogen bonds
- 12 backbone N->O turns; 1 ps resolution, 292136 continuous episodes, median 2.0 ps (break_tol=0).
- Per-start (episodes/distinct-bonds): compact:43398ep/6b; extended:79424ep/7b; mid:122860ep/10b; open:46454ep/4b.

## Caveats
- Chloride: no structured shell (backbone or side-chain) - reportable null.
- Choline-N -> backbone-O: weak feature, highest cutoff-sensitivity - shoulder, not an established shell.
- Side-chain cutoffs are novel (no E3.0 anchor); low reference-atom count -> convergence caveat.
- Non-ergodic backbone: Rg dominated (trapped at 200 ns); SASA robust and pooled/reportable. SASA and Rg are decoupled for this system.

## Question ledger
CLOSED by this system:
  - E3.2 glycerol contact - side-chain and backbone (this solvent)
CONTRIBUTES toward (deferred to cross-solvent stage):
  - E3.1 SASA equilibrium (4 separated states; enhanced sampling -> Paper 2)
  - E3.3 anomaly adjudication (YIY systems)
