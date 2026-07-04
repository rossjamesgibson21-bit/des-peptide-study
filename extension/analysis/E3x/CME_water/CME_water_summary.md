# E3.x system summary - CME_water

Tier A | 4-start (compact, extended, mid, open) | single broad basin (robust)
Peptide CME | solvent water

## E3.1 - Replicate-aware SASA & backbone coordination
- SASA (carried): 3.1505 +/- 0.0841 nm2  [robust]
- Rg   (carried): 0.6959 +/- 0.013 nm  [robust]
- Backbone-O coordination (pair-specific cutoffs, ACF-corrected CIs):
```
partner cutoff_nm pooled_coord   ci95
water_O     0.333       2.4159 0.1875
```

## E3.2 - Side-chain coordination
- Side-chain coordination (pair-specific cutoffs, ACF-corrected CIs):
```
reference partner cutoff_nm pooled_coord   ci95
motif_Cys water_O  no shell       0.0000 0.0000
motif_Met water_O  no shell       0.0000 0.0000
motif_Glu water_O     0.327       5.1293 0.0985
```
- Water only: aqueous side-chain baseline (no DES partners).

## E3.3 - Backbone hydrogen bonds
- 11 backbone N->O turns; 1 ps resolution, 307222 continuous episodes, median 2.0 ps (break_tol=0).
- Per-start (episodes/distinct-bonds): compact:90772ep/8b; extended:52475ep/6b; mid:47095ep/4b; open:116880ep/10b.

## Caveats
- No structured shell against any partner at side-chain site(s): motif_Cys, motif_Met - reportable null.
- Side-chain cutoffs are novel (no E3.0 anchor); low reference-atom count -> convergence caveat.

## Question ledger
CLOSED by this system:
  - E3.1 convergence for this system (both robust)
CONTRIBUTES toward (deferred to cross-solvent stage):
  - E3.1 headline significance (GGE_glyceline)
  - E3.2 indirect-mechanism adjudication (glyceline vs water)
  - E3.3 anomaly adjudication (YIY systems)
