# E3.x system summary - YIY_glyceline

Tier B | 4-start (compact, extended, mid, open) | robust SASA / non-ergodic backbone (Rg dominated by extended/mid)
Peptide YIY | solvent glyceline

## E3.1 - Replicate-aware SASA & backbone coordination
- SASA (carried): 5.2789 +/- 0.0378 nm2  [robust]
- Rg   (carried): 0.6397 +/- 0.0126 nm  [dominated]
- Backbone-O coordination (pair-specific cutoffs, ACF-corrected CIs):
```
   partner cutoff_nm pooled_coord   ci95
   water_O     0.338       4.2288 0.0530
 choline_N     0.423       0.2756 0.0297
  chloride  no shell       0.0000 0.0000
glycerol_O     0.322       0.3206 0.0183
```

## E3.2 - Side-chain coordination
- Side-chain coordination (pair-specific cutoffs, ACF-corrected CIs):
```
reference    partner cutoff_nm pooled_coord   ci95
motif_Tyr    water_O     0.358       5.4553 0.0262
motif_Tyr  choline_N  no shell       0.0000 0.0000
motif_Tyr   chloride     0.382       0.0569 0.0033
motif_Tyr glycerol_O     0.538       2.3387 0.0571
```
- Glycerol side-chain first-shell contact at: motif_Tyr (2.34/O). Backbone-glycerol shell: 0.32/O at 0.322.

## E3.3 - Backbone hydrogen bonds
- 4 backbone N->O turns; 1 ps resolution, 62681 continuous episodes, median 2.0 ps (break_tol=0).
- Per-start (episodes/distinct-bonds): compact:16932ep/3b; extended:16679ep/2b; mid:14495ep/2b; open:14575ep/3b.
- Anomaly-control arm: glyceline backbone-H-bond baseline for the YIY_reline comparison (1 ps backbone-only scope; not zero at this scope, unlike the Phase-3 region-only metric).

## Caveats
- Chloride: no backbone shell; structured shell at motif_Tyr-chloride (0.057) (low occupancy - verify).
- Choline-N -> backbone-O: weak feature, highest cutoff-sensitivity - shoulder, not an established shell.
- Side-chain cutoffs are novel (no E3.0 anchor); low reference-atom count -> convergence caveat.
- Non-ergodic backbone: Rg dominated (trapped at 200 ns); SASA robust and pooled/reportable. SASA and Rg are decoupled for this system.

## Question ledger
CLOSED by this system:
  - E3.2 glycerol contact - side-chain and backbone (this solvent)
CONTRIBUTES toward (deferred to cross-solvent stage):
  - E3.1 headline significance (GGE_glyceline)
  - E3.3 anomaly adjudication (needs YIY_reline)
