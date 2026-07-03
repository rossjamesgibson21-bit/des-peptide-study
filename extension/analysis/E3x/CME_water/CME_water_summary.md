# E3.x system summary - CME_water

Tier A | 4-start (compact, extended, mid, open) | single broad basin (robust)
Peptide CME | solvent water

## E3.1 - Replicate-aware SASA & backbone coordination
- SASA (carried E2.2 pooled): 3.1505 +/- 0.0841 nm2  [robust]
- Rg   (carried E2.2 pooled): 0.6959 +/- 0.013 nm  [robust]
- Convergence character: single broad basin (robust)
- Backbone-O coordination (pair-specific cutoffs, ACF-corrected CIs):
```
partner cutoff_nm pooled_coord   ci95
water_O     0.333       2.4159 0.1875
```
- Replicate significance for THIS system: achieved (both observables robust). The named E3.1 test case (GGE_glyceline, corrected ES 1.18) is a separate system; its significance verdict is deferred there.

## E3.2 - Side-chain coordination (indirect-restructuring test)
- Side-chain coordination (pair-specific cutoffs, ACF-corrected CIs):
```
reference partner cutoff_nm pooled_coord   ci95
motif_Cys water_O  no shell       0.0000 0.0000
motif_Met water_O  no shell       0.0000 0.0000
motif_Glu water_O     0.327       5.1293 0.0985
```
- CONTRIBUTES (baseline): water carries no glycerol, so it cannot test the glycerol-sidechain question directly. It supplies the aqueous baseline that rules out a peptide-intrinsic explanation. The decisive glycerol comparison is deferred to the glyceline systems + water baseline.
- Observation: motif vs flanking carboxylate coordinate near-uniformly per O (uniform acidic behaviour rather than local context) - a hypothesis the other systems test.

## E3.3 - Backbone hydrogen bonds (anomaly test)
- 11 backbone i->i-3 turns; 1 ps resolution, 307222 continuous episodes, median 2.0 ps (break_tol=0).
- CONTRIBUTES: the H-bond anomaly is YIY reline; this system validated the H-bond machinery and established the 1 ps resolution requirement (stride=10 aliases ~87% of episodes). Anomaly adjudication deferred to the YIY systems.
- Note: the continuous-episode median is largely blind to slow exchange; the persistence definition (continuous vs intermittent/ACF) is fixed at the cross-solvent stage.

## Carried caveats
- Chloride: no structured shell (backbone or side-chain) - reportable null.
- No structured hydration/coordination shell (flat RDF) at side-chain site(s): motif_Cys, motif_Met - reportable null, not a failure (e.g. hydrophobic Cys/Met sulfur).
- Side-chain cutoffs are novel (no E3.0 anchor); low reference-atom count -> convergence caveat.

## Question ledger
CLOSED by this system:
  - E3.1 convergence for this system (robust)
CONTRIBUTES toward (deferred to cross-solvent stage):
  - E3.1 headline significance (GGE_glyceline ES 1.18)
  - E3.2 indirect-mechanism adjudication (glyceline + water)
  - E3.3 anomaly adjudication (YIY systems)
