# E3.x system summary - CME_reline

Tier C | 4-start (compact, extended, mid, open) | one-trap (compact/mid/open consensus + metastable extended)
Peptide CME | solvent reline

## E3.1 - Replicate-aware SASA & backbone coordination
- SASA is one-trap: a metastable start sits apart from the consensus. Reported both ways
  (the trap is a sampling outcome, not an equilibrium basin):
  - FULL pool (all starts): 3.399 +/- 0.0792 nm2
  - CONSENSUS (drop trap: extended): 3.3073 +/- 0.0703 nm2
  - inclusion shift (full - consensus): +0.0917 nm2
  - E2.2 sidecar cross-check: 3.399 (delta 0.0)
- Rg (carried E2.2 pooled): 0.7335 +/- 0.0227 nm  [robust]
- Convergence character: one-trap (compact/mid/open consensus + metastable extended)
- Backbone-O coordination (pooled; local, trap-independent):
```
  partner cutoff_nm pooled_coord   ci95
  water_O     0.333       2.6493 0.1675
choline_N     0.407       0.1508 0.0187
 chloride  no shell       0.0000 0.0000
   urea_N  no shell       0.0000 0.0000
   urea_O  no shell       0.0000 0.0000
```
- CONTRIBUTES: SASA carries a trap-inclusion sensitivity; the consensus value is the more defensible central estimate, with the full-pool shift reported alongside.

## E3.2 - Side-chain coordination (indirect-restructuring test)
- Side-chain coordination (pair-specific cutoffs, ACF-corrected CIs):
```
reference   partner cutoff_nm pooled_coord   ci95
motif_Cys   water_O  no shell       0.0000 0.0000
motif_Cys choline_N     0.433       0.0739 0.0052
motif_Cys  chloride  no shell       0.0000 0.0000
motif_Cys    urea_N     0.447       0.4379 0.0193
motif_Cys    urea_O     0.377       0.0937 0.0053
motif_Met   water_O  no shell       0.0000 0.0000
motif_Met choline_N     0.447       0.1245 0.0074
motif_Met  chloride  no shell       0.0000 0.0000
motif_Met    urea_N     0.433       0.3793 0.0157
motif_Met    urea_O  no shell       0.0000 0.0000
motif_Glu   water_O     0.327       5.2250 0.0853
motif_Glu choline_N     0.527       0.4816 0.0283
motif_Glu  chloride  no shell       0.0000 0.0000
motif_Glu    urea_N     0.392       0.4192 0.0209
motif_Glu    urea_O  no shell       0.0000 0.0000
```
- CONTRIBUTES (baseline): reline carries no glycerol, so it cannot test the glycerol-sidechain question directly. It supplies the reline direct-contact baseline: direct backbone contact (choline_N); side-chain contact at motif_Cys-choline_N, motif_Cys-urea_N, motif_Cys-urea_O, motif_Met-choline_N, motif_Met-urea_N, motif_Glu-choline_N, motif_Glu-urea_N. The decisive glycerol comparison is deferred to the glyceline systems + water baseline.
- Observation: motif vs flanking carboxylate coordinate near-uniformly per O (uniform acidic behaviour rather than local context) - a hypothesis the other systems test.

## E3.3 - Backbone hydrogen bonds (anomaly test)
- 9 backbone i->i-3 turns; 1 ps resolution, 157729 continuous episodes, median 2.0 ps (break_tol=0).
- CONTRIBUTES: the H-bond anomaly is YIY reline; this system validated the H-bond machinery and established the 1 ps resolution requirement (stride=10 aliases ~87% of episodes). Anomaly adjudication deferred to the YIY systems.
- Note: the continuous-episode median is largely blind to slow exchange; the persistence definition (continuous vs intermittent/ACF) is fixed at the cross-solvent stage.

## Carried caveats
- Chloride: no structured shell (backbone or side-chain) - reportable null.
- Choline-N -> backbone-O registers a weak feature absent in the E3.0 preliminary and carrying the highest cutoff-sensitivity - treated as a shoulder, not an established shell.
- Side-chain cutoffs are novel (no E3.0 anchor); low reference-atom count -> convergence caveat.

## Question ledger
CLOSED by this system:
  - (none - validation/contribution system)
CONTRIBUTES toward (deferred to cross-solvent stage):
  - E3.1 SASA (one-trap; consensus vs full-pool sensitivity reported)
  - E3.2 indirect-mechanism adjudication (glyceline + water)
  - E3.3 anomaly adjudication (YIY systems)
