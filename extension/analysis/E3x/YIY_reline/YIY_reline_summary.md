# E3.x system summary - YIY_reline

Tier B | 4-start (compact, extended, mid, open) | robust SASA / non-ergodic backbone (Rg dominated by compact)
Peptide YIY | solvent reline

## E3.1 - Replicate-aware SASA & backbone coordination
- SASA (carried E2.2 pooled): 5.4118 +/- 0.0607 nm2  [robust]
- Rg   (carried E2.2 pooled): 0.6603 +/- 0.0098 nm  [dominated]
- Convergence character: robust SASA / non-ergodic backbone (Rg dominated by compact)
- Backbone-O coordination (pair-specific cutoffs, ACF-corrected CIs):
```
  partner cutoff_nm pooled_coord   ci95
  water_O     0.338       4.2702 0.0495
choline_N     0.417       0.2229 0.0208
 chloride  no shell       0.0000 0.0000
   urea_N     0.407       0.8174 0.0345
   urea_O     0.562       1.5151 0.0693
```
- Replicate significance for THIS system: see convergence character. The named E3.1 test case (GGE_glyceline, corrected ES 1.18) is a separate system; its significance verdict is deferred there.

## E3.2 - Side-chain coordination (indirect-restructuring test)
- Side-chain coordination (pair-specific cutoffs, ACF-corrected CIs):
```
reference   partner cutoff_nm pooled_coord   ci95
motif_Tyr   water_O     0.358       5.3810 0.0234
motif_Tyr choline_N  no shell       0.0000 0.0000
motif_Tyr  chloride     0.377       0.0749 0.0040
motif_Tyr    urea_N     0.387       0.5546 0.0139
motif_Tyr    urea_O     0.312       0.0998 0.0053
```
- CONTRIBUTES (baseline): reline carries no glycerol, so it cannot test the glycerol-sidechain question directly. It supplies the reline direct-contact baseline: direct backbone contact (urea_N, urea_O); side-chain contact at motif_Tyr-chloride, motif_Tyr-urea_N, motif_Tyr-urea_O. The decisive glycerol comparison is deferred to the glyceline systems + water baseline.

## E3.3 - Backbone hydrogen bonds (anomaly test)
- 3 backbone N->O turns; 1 ps resolution, 45678 continuous episodes, median 2.0 ps (break_tol=0).
- Per-start (episodes/distinct-bonds): compact:22335ep/3b; extended:13773ep/2b; mid:0ep/0b; open:9570ep/2b. A start-dependent backbone-H-bond pattern shows here (relevant to the YIY reline anomaly, which is start-specific).
- CLOSES (in part): YIY_reline is the anomaly system. Backbone-H-bond counts here read against the YIY water and glyceline baselines (at matched 1 ps backbone-only scope, not the Phase-3 region-only metric); a start-dependent or solvent-specific pattern is the load-bearing E3.3 result.
- Note: the continuous-episode median is largely blind to slow exchange; the persistence definition (continuous vs intermittent/ACF) is fixed at the cross-solvent stage.

## Carried caveats
- Chloride: no backbone shell; structured shell at motif_Tyr-chloride (0.075) (low occupancy - verify).
- Choline-N -> backbone-O registers a weak feature absent in the E3.0 preliminary and carrying the highest cutoff-sensitivity - treated as a shoulder, not an established shell.
- Side-chain cutoffs are novel (no E3.0 anchor); low reference-atom count -> convergence caveat.
- Non-ergodic backbone: Rg dominated (trapped at 200 ns); SASA robust and pooled/reportable. SASA and Rg are decoupled for this system.

## Question ledger
CLOSED by this system:
  - E3.3 YIY reline H-bond anomaly (vs water/glyceline baselines)
CONTRIBUTES toward (deferred to cross-solvent stage):
  - E3.1 headline significance (GGE_glyceline ES 1.18)
  - E3.2 indirect-mechanism adjudication (glyceline + water)
