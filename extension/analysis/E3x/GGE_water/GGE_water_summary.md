# E3.x system summary - GGE_water

Tier D | 4-start (compact, extended, mid, open) | multi-basin (2+2 SASA split)
Peptide GGE | solvent water

## E3.1 - Replicate-aware SASA & backbone coordination
- SASA is the DIVERGENT observable (multi-basin) - reported as a diverse-start ENSEMBLE,
  NOT an equilibrium mean. Grouped-start basin means (corrected CIs):
  - HIGH basin (extended;open): 2.8443 +/- 0.0617 nm2
  - LOW basin  (compact;mid): 2.3357 +/- 0.0882 nm2
  - Pooled ensemble (NOT equilibrium): 2.59 +/- 0.0656 nm2
  - E2.2 sidecar cross-check: 2.59 (delta 0.0)
- Rg (carried E2.2 pooled): 0.704 +/- 0.0176 nm  [dominated]
- Convergence character: multi-basin (2+2 SASA split)
- Backbone-O coordination (pooled; local, basin-independent):
```
partner  cutoff_nm  pooled_coord  ci95
water_O      0.333        4.8929 0.255
```
- CONTRIBUTES: this water system supplies the aqueous SASA baseline, but as a multi-basin ensemble the equilibrium mean is unresolved; true basin populations need enhanced sampling.

## E3.2 - Side-chain coordination (indirect-restructuring test)
- Side-chain coordination (pair-specific cutoffs, ACF-corrected CIs):
```
reference partner  cutoff_nm  pooled_coord   ci95
motif_Glu water_O      0.322        6.3564 0.0201
flank_Glu water_O      0.322       19.3362 0.0215
```
- CONTRIBUTES (baseline): water carries no glycerol, so it cannot test the glycerol-sidechain question directly. It supplies the aqueous baseline that rules out a peptide-intrinsic explanation. The decisive glycerol comparison is deferred to the glyceline systems + water baseline.
- Observation: motif vs flanking carboxylate coordinate near-uniformly per O (uniform acidic behaviour rather than local context) - a hypothesis the other systems test.

## E3.3 - Backbone hydrogen bonds (anomaly test)
- 8 backbone i->i-3 turns; 1 ps resolution, 146682 continuous episodes, median 2.0 ps (break_tol=0).
- CONTRIBUTES: the H-bond anomaly is YIY reline; this system validated the H-bond machinery and established the 1 ps resolution requirement (stride=10 aliases ~87% of episodes). Anomaly adjudication deferred to the YIY systems.
- Note: the continuous-episode median is largely blind to slow exchange; the persistence definition (continuous vs intermittent/ACF) is fixed at the cross-solvent stage.

## Carried caveats
- Chloride: no structured shell (backbone or side-chain) - reportable null.
- Side-chain cutoffs are novel (no E3.0 anchor); low reference-atom count -> convergence caveat.

## Question ledger
CLOSED by this system:
  - (none - validation/contribution system)
CONTRIBUTES toward (deferred to cross-solvent stage):
  - E3.1 SASA equilibrium mean (multi-basin; enhanced sampling -> Paper 2)
  - E3.2 indirect-mechanism adjudication (glyceline + water)
  - E3.3 anomaly adjudication (YIY systems)
