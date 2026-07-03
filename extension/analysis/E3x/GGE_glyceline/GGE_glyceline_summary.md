# E3.x system summary - GGE_glyceline

Tier A | 4-start (compact, extended, mid, open) | clean / ergodic
Peptide GGE | solvent glyceline

## E3.1 - Replicate-aware SASA & backbone coordination
- SASA (carried E2.2 pooled): 2.6122 +/- 0.0638 nm2  [robust]
- Rg   (carried E2.2 pooled): 0.7063 +/- 0.0176 nm  [robust]
- Convergence character: clean / ergodic
- Backbone-O coordination (pair-specific cutoffs, ACF-corrected CIs):
```
   partner  cutoff_nm  pooled_coord   ci95
   water_O      0.333        3.9949 0.1450
 choline_N      0.412        0.2632 0.0235
  chloride        NaN        0.0000 0.0000
glycerol_O        NaN        0.0000 0.0000
```
- CLOSES: named E3.1 significance test case (Phase-3 corrected ES 1.18); the verdict against corrected CIs is the load-bearing E3.1 result for the campaign.

## E3.2 - Side-chain coordination (indirect-restructuring test)
- Side-chain coordination (pair-specific cutoffs, ACF-corrected CIs):
```
reference    partner  cutoff_nm  pooled_coord   ci95
motif_Glu    water_O      0.327        5.5090 0.0600
motif_Glu  choline_N      0.522        0.4435 0.0234
motif_Glu   chloride        NaN        0.0000 0.0000
motif_Glu glycerol_O      0.317        0.6309 0.0433
flank_Glu    water_O      0.327       16.7601 0.0917
flank_Glu  choline_N      0.527        1.3800 0.0452
flank_Glu   chloride        NaN        0.0000 0.0000
flank_Glu glycerol_O      0.317        1.8527 0.0718
```
- CLOSES (in part): glyceline system - glycerol-sidechain coordination is the direct test of the indirect hypothesis. Substantial glycerol-Tyr/Cys/Glu contact challenges purely-indirect; absence supports it. Read against the water baseline.
- Observation: motif vs flanking carboxylate coordinate near-uniformly per O (uniform acidic behaviour rather than local context) - a hypothesis the other systems test.

## E3.3 - Backbone hydrogen bonds (anomaly test)
- 8 backbone i->i-3 turns; 1 ps resolution, 112691 continuous episodes, median 2.0 ps (break_tol=0).
- CONTRIBUTES: the H-bond anomaly is YIY reline; this system validated the H-bond machinery and established the 1 ps resolution requirement (stride=10 aliases ~87% of episodes). Anomaly adjudication deferred to the YIY systems.
- Note: the continuous-episode median is largely blind to slow exchange; the persistence definition (continuous vs intermittent/ACF) is fixed at the cross-solvent stage.

## Carried caveats
- Chloride: no structured shell (backbone or side-chain) - reportable null.
- Choline-N -> backbone-O registers a weak feature absent in the E3.0 preliminary and carrying the highest cutoff-sensitivity - treated as a shoulder, not an established shell.
- Side-chain cutoffs are novel (no E3.0 anchor); low reference-atom count -> convergence caveat.

## Question ledger
CLOSED by this system:
  - E3.1 significance (named test case)
  - E3.2 glycerol-sidechain contact (this solvent)
CONTRIBUTES toward (deferred to cross-solvent stage):
  - E3.3 anomaly adjudication (YIY systems)
