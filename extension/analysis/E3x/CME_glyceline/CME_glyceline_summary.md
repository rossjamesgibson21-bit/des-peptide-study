# E3.x system summary - CME_glyceline

Tier D | 4-start (compact, extended, mid, open) | multi-basin (4 separated SASA states; heaviest caveat)
Peptide CME | solvent glyceline

## E3.1 - Replicate-aware SASA & backbone coordination
- SASA is the DIVERGENT observable: 4 statistically separated state(s) across the diverse starts (adjacent gap > combined 95% CI). Reported as a diverse-start ENSEMBLE,
  NOT an equilibrium mean. Per-start means (the honest record):
```
  mid:2.8831(ci0.1014)
  extended:3.1479(ci0.1206)
  compact:3.3110(ci0.0667)
  open:3.6189(ci0.1020)
```
- Resolved states (grouped only where a gap test cannot separate): mid:2.8831+/-0.1014 | extended:3.1479+/-0.1206 | compact:3.3110+/-0.0667 | open:3.6189+/-0.1020
- Pooled ensemble (NOT equilibrium): 3.2402 +/- 0.0602 nm2
- E2.2 sidecar cross-check: 3.2402 (delta 0.0)
- Rg (carried E2.2 pooled): 0.691 +/- 0.0155 nm  [dominated]
- Convergence character: multi-basin (4 separated SASA states; heaviest caveat)
- Backbone-O coordination (pooled; local, state-independent):
```
   partner cutoff_nm pooled_coord   ci95
   water_O     0.333       2.6139 0.1024
 choline_N     0.412       0.1519 0.0155
  chloride  no shell       0.0000 0.0000
glycerol_O  no shell       0.0000 0.0000
```
- CONTRIBUTES: 4 separated SASA states, no two-basin structure; the equilibrium mean and state populations are unresolved (enhanced sampling -> Paper 2). Heaviest caveat.

## E3.2 - Side-chain coordination (indirect-restructuring test)
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
- CLOSES (in part): glycerol makes direct first-shell contact at motif_Cys (0.38/O), motif_Glu (0.50/O); no glycerol shell at motif_Met. Direct side-chain contact qualifies the purely-indirect reading; read against the water baseline (no side-chain shells there).
- Observation: motif vs flanking carboxylate coordinate near-uniformly per O (uniform acidic behaviour rather than local context) - a hypothesis the other systems test.

## E3.3 - Backbone hydrogen bonds (anomaly test)
- 12 backbone N->O turns; 1 ps resolution, 292136 continuous episodes, median 2.0 ps (break_tol=0).
- CONTRIBUTES: the H-bond anomaly is YIY reline; this system validated the H-bond machinery and established the 1 ps resolution requirement (stride=10 aliases ~87% of episodes). Anomaly adjudication deferred to the YIY systems.
- Note: the continuous-episode median is largely blind to slow exchange; the persistence definition (continuous vs intermittent/ACF) is fixed at the cross-solvent stage.

## Carried caveats
- Chloride: no structured shell (backbone or side-chain) - reportable null.
- Choline-N -> backbone-O registers a weak feature absent in the E3.0 preliminary and carrying the highest cutoff-sensitivity - treated as a shoulder, not an established shell.
- Side-chain cutoffs are novel (no E3.0 anchor); low reference-atom count -> convergence caveat.

## Question ledger
CLOSED by this system:
  - E3.2 glycerol-sidechain contact (this solvent)
CONTRIBUTES toward (deferred to cross-solvent stage):
  - E3.1 SASA equilibrium (4 separated states; enhanced sampling -> Paper 2)
  - E3.3 anomaly adjudication (YIY systems)
