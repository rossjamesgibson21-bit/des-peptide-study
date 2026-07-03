# E3.x system summary - GGE_reline

Tier A | 3-start (compact, extended, mid) | clean / ergodic
Peptide GGE | solvent reline

## E3.1 - Replicate-aware SASA & backbone coordination
- SASA (carried E2.2 pooled): 2.8056 +/- 0.0662 nm2  [robust]
- Rg   (carried E2.2 pooled): 0.7636 +/- 0.0189 nm  [robust]
- Convergence character: clean / ergodic
- Backbone-O coordination (pair-specific cutoffs, ACF-corrected CIs):
```
  partner  cutoff_nm  pooled_coord   ci95
  water_O      0.333        4.6432 0.0976
choline_N      0.412        0.3969 0.0349
 chloride        NaN        0.0000 0.0000
   urea_N      0.423        0.9229 0.0454
   urea_O      0.562        1.2849 0.0565
```
- Replicate significance for THIS system: achieved (both observables robust). The named E3.1 test case (GGE_glyceline, corrected ES 1.18) is a separate system; deferred there.

## E3.2 - Side-chain coordination (indirect-restructuring test)
- Side-chain coordination (pair-specific cutoffs, ACF-corrected CIs):
```
reference   partner  cutoff_nm  pooled_coord   ci95
motif_Glu   water_O      0.322        5.7587 0.0437
motif_Glu choline_N      0.548        0.6285 0.0351
motif_Glu  chloride        NaN        0.0000 0.0000
motif_Glu    urea_N      0.352        0.3087 0.0188
motif_Glu    urea_O        NaN        0.0000 0.0000
flank_Glu   water_O      0.327       17.9749 0.0530
flank_Glu choline_N      0.522        1.6996 0.0506
flank_Glu  chloride        NaN        0.0000 0.0000
flank_Glu    urea_N      0.392        1.3552 0.0437
flank_Glu    urea_O        NaN        0.0000 0.0000
```
- CONTRIBUTES (baseline): reline carries no glycerol, so it cannot test the glycerol-sidechain question directly. It supplies the urea direct-contact baseline (urea coordinates backbone AND carboxylate here). The decisive glycerol comparison is deferred to the glyceline systems + water baseline.
- Observation: motif vs flanking carboxylate coordinate near-uniformly per O (uniform acidic behaviour rather than local context) - a hypothesis the other systems test.

## E3.3 - Backbone hydrogen bonds (anomaly test)
- 4 backbone i->i-3 turns; 1 ps resolution, 36061 continuous episodes, median 2.0 ps (break_tol=0).
- CONTRIBUTES: the H-bond anomaly is YIY reline; this system validated the H-bond machinery and established the 1 ps resolution requirement (stride=10 aliases ~87% of episodes). Anomaly adjudication deferred to the YIY systems.
- Note: the continuous-episode median is largely blind to slow exchange; the persistence definition (continuous vs intermittent/ACF) is fixed at the cross-solvent stage.

## Carried caveats
- Chloride: no structured shell (backbone or side-chain) - reportable null.
- Choline-N -> backbone-O registers a weak feature absent in the E3.0 preliminary and carrying the highest cutoff-sensitivity - treated as a shoulder, not an established shell.
- Side-chain cutoffs are novel (no E3.0 anchor); low reference-atom count -> convergence caveat.

## Question ledger
CLOSED by this system:
  - E3.1 convergence for this system (robust)
CONTRIBUTES toward (deferred to cross-solvent stage):
  - E3.1 headline significance (GGE_glyceline ES 1.18)
  - E3.2 indirect-mechanism adjudication (glyceline + water)
  - E3.3 anomaly adjudication (YIY systems)
