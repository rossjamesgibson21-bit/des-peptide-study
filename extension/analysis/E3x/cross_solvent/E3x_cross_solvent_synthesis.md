# E3.x Cross-Solvent Synthesis

Observations and hypotheses read across the nine per-system E3.x summaries at matched cutoffs. Supported by the committed section outputs; not conclusions of the phase. External refs as [REF].


## E3.1 — SASA significance

- Across the matrix, all three reline contrasts exclude zero (GGE +0.216, CME +0.248, YIY +0.135) and all three glyceline contrasts include it (GGE +0.022, CME +0.090, YIY +0.002): under corrected CIs the SASA-opening signal segregates with reline, not glyceline. [REF]
  - GGE reline - water: +0.216 [+0.122, +0.309] nm^2, corrected 95% CI excludes zero — water baseline tier-D (ensemble, not equilibrium).
  - GGE glyceline - water: +0.022 [-0.069, +0.114] nm^2, corrected 95% CI includes zero — water baseline tier-D (ensemble, not equilibrium).
  - CME reline - water: +0.248 [+0.133, +0.364] nm^2, corrected 95% CI excludes zero.
  - CME glyceline - water: +0.090 [-0.014, +0.193] nm^2, corrected 95% CI includes zero — des tier-D (ensemble, not equilibrium).
  - YIY reline - water: +0.135 [+0.062, +0.208] nm^2, corrected 95% CI excludes zero.
  - YIY glyceline - water: +0.002 [-0.054, +0.058] nm^2, corrected 95% CI includes zero.
- Load-bearing (GGE_glyceline): the carried Phase-3 effect size (1.18) does not survive — dSASA +0.022 [-0.069, +0.114] includes zero. The water baseline is tier-D bimodal (high 2.844±0.062, low 2.336±0.088); glyceline (2.612) falls between the basins, so the null reflects baseline non-equilibrium as much as any glyceline effect. Equilibrium water populations are deferred (Paper 2). [REF]
- GGE_reline (2.806) coincides with the high water basin (2.844), not the pooled mean, so its +0.216 against the pooled baseline is a direction-robust, magnitude-provisional estimate. [REF]

## E3.2 — direct contact versus indirect restructuring

- Backbone hydration is displaced by both deep eutectics, not reline alone, and by essentially equal amounts. At matched cutoff on the three-atom YIY motif backbone (waters per carbonyl O): water 1.80, reline 1.42, glyceline 1.41; the reline−glyceline difference +0.014 [±0.025] does not resolve a difference (statistically equal displacement). Backbone dehydration is a shared DES property. [REF: e3.2_matched_cutoff.csv]
- This decouples backbone dehydration from the SASA effect: reline opens YIY SASA (+0.135, excludes zero) while glyceline does not (+0.002, includes zero), yet both displace backbone water equally — the SASA opening is not explained by backbone water displacement. [REF]
- What distinguishes reline is direct urea occupancy of the vacated sites, not the vacancy. Urea contacts the motif backbone in GGE_reline (urea-O 0.428, urea-N 0.308 per O) and YIY_reline (urea-O 0.505); glycerol reaches the YIY backbone only weakly (0.107) and the GGE backbone not at all. Reline replaces backbone water with urea where glyceline leaves it unreplaced. [REF]
- The contact mode is not uniform, so no single mechanism accounts for the SASA opening. CME_reline opens SASA (+0.248, excludes zero) with no urea–backbone shell — urea instead occupies the CME side chains (Cys urea-N 0.438, urea-O 0.094; Met urea-N 0.379). Direct urea occupancy is the reline-specific feature, but its location (backbone in GGE/YIY, side chain in CME) differs by peptide. [REF]
- Read together these are mechanism-constraining observations, not an established mechanism: backbone dehydration is shared and does not track SASA opening; reline's distinguishing feature is direct urea occupancy; the occupancy site differs by peptide. Choline coordinates backbone and side chains pervasively in both DES (§2a) and is not a discriminating variable. Equilibrium and causal status are unresolved (tier caveats on GGE_water and the multi-basin systems stand); a causal test of urea occupancy against SASA is deferred (Paper 2). [REF]

## E3.3 — YIY backbone hydrogen bonding

- Across the YIY row, reline is the only solvent whose four diverse starts include one (mid) that forms no motif backbone N→O turns across the full 20–200 ns window. Water and glyceline populate every start (water 34,007 / 9,761 / 11,074 / 9,743 episodes; glyceline 16,932 / 16,679 / 14,495 / 14,575, the tightest spread of the three); reline reads 22,335 / 13,773 / 0 / 9,570. In aggregate reline totals 45,678 episodes over 3 distinct turns, against water 64,585 / 4 and glyceline 62,681 / 4 — lower on both counts and the only row member missing a turn from its inventory. [REF]

- The zero is a real conformational property of that trajectory rather than a detection artifact. An independent backbone hydrogen-bond pass over all 180,000 frames returned no donor–acceptor pair meeting the N→O geometry, which excludes the artifact classes that produce a spurious zero: a frame-loading failure truncates the count, stride aliasing spares long-lived bonds, and a scope or cutoff error misassigns rather than eliminates. Over the window the mid start occupies a backbone whose turn geometry is never satisfied. [REF]

- Because diverse starts are seeded independently per system, reline's mid start and the mid starts of water and glyceline are separately generated conformations sharing a label, not one structure held in three solvents. The supported statement is therefore that reline is the only solvent whose start ensemble includes a turn-free backbone basin, not that a shared conformation loses its turns in reline. [REF]

- Time-resolved backbone geometry characterises the state as neither a rigid single basin nor a directional drift, but a mobile turn-free ensemble. φ is held across all residues (per-residue circular dispersion ≤ 34°) while ψ at the motif positions TYR3 and ILE4 interconverts between two discrete non-turn states on a tens-of-ns timescale (ψ dispersion 71° and 87°; the third motif tyrosine TYR5 comparatively ordered at 37°). Backbone Rg spans 0.281 nm about a stable mean with no net drift (first- to second-half shift +0.013 nm), its compact and extended excursions coinciding in time with the ψ state-blocks. The zero turn count thus reflects occupancy of two specific turn-incompatible ψ basins, not backbone disorder. [REF: q3_yiy_reline_mid_backbone.png; q3_yiy_reline_mid_timeseries.npz]

- Scope and ceiling. The metric is matched across the YIY row (1 ps resolution, backbone N→O scope, continuous episodes), but the per-start episode ordering is not resolved at the start level and is fragile to the dominant start, so it is recorded as the categorical zero-turn observation above rather than a quantitative ranking. Relative to the Phase-3 region-scoped, strided reading, the 1 ps backbone-scope metric reveals a turn-free state the earlier metric could not resolve — this is not a correction of the earlier direction. All characterisation is on the 200 ns single-start timescale, and whether the turn-free ψ ensemble carries genuine equilibrium weight or is kinetically stranded is not separable from this data (non-ergodic backbone across the YIY row); that separation is deferred to enhanced sampling (Paper 2). [REF]

- The observation coheres with the one fully converged E3.1 contrast (YIY reline − water SASA +0.135 [+0.062, +0.208]): a reline-accessible, more open YIY backbone forming fewer intramolecular turns is consistent across two independent observables. Offered as a hypothesis, since the backbone signal rests on a single start. [REF]
