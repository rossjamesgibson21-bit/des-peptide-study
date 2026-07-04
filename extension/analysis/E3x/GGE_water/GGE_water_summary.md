# E3.x system summary - GGE_water

Tier D | 4-start (compact, extended, mid, open) | multi-basin (2+2 SASA split)
Peptide GGE | solvent water

## E3.1 - Replicate-aware SASA & backbone coordination
- SASA: 2 statistically separated state(s) across the diverse starts (adjacent gap > combined 95% CI); diverse-start ensemble, not an equilibrium mean.
- Per-start SASA means (nm2):
```
  compact:2.3824
  extended:2.7884
  mid:2.2891
  open:2.9003
```
- Grouped states: HIGH:2.8443+/-0.0617 (extended;open) | LOW:2.3357+/-0.0882 (compact;mid)
- Pooled ensemble: 2.59 +/- 0.0656 nm2 (sidecar 2.59, delta 0.0)
- Rg (carried): 0.704 +/- 0.0176 nm  [dominated]
- Backbone-O coordination (pooled):
```
partner cutoff_nm pooled_coord   ci95
water_O     0.333       4.8929 0.2550
```

## E3.2 - Side-chain coordination
- Side-chain coordination (pair-specific cutoffs, ACF-corrected CIs):
```
reference partner cutoff_nm pooled_coord   ci95
motif_Glu water_O     0.322       6.3564 0.0201
flank_Glu water_O     0.322      19.3362 0.0215
```
- Water only: aqueous side-chain baseline (no DES partners).

## E3.3 - Backbone hydrogen bonds
- 8 backbone N->O turns; 1 ps resolution, 146682 continuous episodes, median 2.0 ps (break_tol=0).
- Per-start (episodes/distinct-bonds): compact:51849ep/5b; extended:20630ep/3b; mid:56055ep/6b; open:18148ep/3b.

## Caveats
- Side-chain cutoffs are novel (no E3.0 anchor); low reference-atom count -> convergence caveat.
- Non-ergodic backbone: Rg dominated (trapped at 200 ns); SASA robust and pooled/reportable. SASA and Rg are decoupled for this system.

## Question ledger
CLOSED by this system:
  - (none - validation/baseline system)
CONTRIBUTES toward (deferred to cross-solvent stage):
  - E3.1 SASA equilibrium (2 separated states; enhanced sampling -> Paper 2)
  - E3.2 indirect-mechanism adjudication (glyceline vs water)
  - E3.3 anomaly adjudication (YIY systems)
