# DES Peptide Simulation Study — Project README

## Title

Solvent-Dependent Modulation of Bioactive Motif Accessibility in *Torreya grandis* Peptides: A Molecular Dynamics Study in Water and Choline Chloride-Based Deep Eutectic Solvents

## Author

Ross Gibson

## Last updated

April 2026 — Phase 3 complete

------------------------------------------------------------------------

## Background

This project is a computational follow-up to Durrani et al. (Zhejiang A&F University), who used ChCl:urea deep eutectic solvent to extract proteins from *Torreya grandis* seeds and identified 600 peptides via LC-MS/MS. Of 87 functional oligopeptides found in the DES extract (vs 41 in water), three distinct antioxidant motifs were highlighted: **GGE**, **CME**, and **YIY**.

This study uses molecular dynamics simulations to investigate *why* DES extraction yields more bioactive peptides — specifically, whether DES modulates the solvent accessibility and structural dynamics of these motifs compared to water.

------------------------------------------------------------------------

## Key decisions and rationale

### Peptide selection

**Problem:** The original study identified GGE, CME, and YIY as motifs of interest, but finding natural proteins containing these motifs with reliable structural predictions proved difficult.

**What was done:** - 87 peptide sequences from Durrani et al. were BLASTed across a phased taxonomic expansion (T. grandis → Taxaceae → Pinus → Plants) - \~30 candidate proteins were evaluated via AlphaFold pTM scores - Most T. grandis BLAST hits were from chloroplast/mitochondrial genome proteins or hypothetical proteins from the *Taxus chinensis* genome assembly (KAH9... accessions) — computationally predicted gene models with no experimental validation - The only natural, well-characterised protein containing GGE with a reliable AlphaFold structure was **A0A445BRA0** from *Arachis hypogaea* (peanut) — a 43-residue ribosomal protein (L2-GGE) - CME and YIY were only found in hypothetical *Taxus chinensis* proteins (KAH9325454.1 and KAH9327367.1) — 736 and 481 residues respectively, both too large for tractable MD on available hardware

**Decision:** Use **short capped peptide constructs** derived from the original *T. grandis* peptide sequences identified by LC-MS/MS, rather than full-length host proteins. This provides: - Compute feasibility (6–9 residues vs 400–700 residues) - Consistency across all three motifs (same system type) - Direct connection to the experimental source material - Scientific defensibility (studying the actual peptide sequences, not unrelated host proteins)

**Phase 1 pilot:** The L2-GGE full protein (43 residues) was simulated in water as a protocol validation exercise. This confirmed the simulation setup works and revealed indexing errors in the analysis pipeline that were subsequently corrected.

### Peptide constructs

All flanking residues are drawn from the original *T. grandis* peptide sequences in Durrani et al.:

| Motif | Source (Durrani) | Construct | Residues | Motif position (0-indexed in construct) |
|------------|------------|------------|------------|--------------------------|
| GGE | NTDEEEGGEIVFGGVDPK (Peptide 6) | Ace-EEEGGEIVF-NMe | 9 | 3–5 |
| CME | LYQCMEFVR (Peptide 5) | Ace-LYQCMEFVR-NMe | 9 | 3–5 |
| YIY | NPYIYK (Peptide 83) | Ace-NPYIYK-NMe | 6 | 2–4 |

Ace = acetyl cap (N-terminus), NMe = N-methylamide cap (C-terminus). Caps eliminate terminal charge artefacts so peptides behave like internal protein fragments.

### Force field choices

| Component | Force field | Charges | Justification |
|-----------------|-----------------|-----------------|----------------------|
| Peptides | Amber14SB (ff14SB) | Native | Standard for protein/peptide MD with Amber |
| Caps (ACE/NME) | Amber14SB | Native | Built-in support |
| Water | TIP3P | Native | Consistent with Amber14SB parameterisation |
| Choline cation | GAFF2 | AM1-BCC, scaled ×0.8 | GAFF2 is the standard small-molecule companion to Amber; charge scaling corrects for charge-transfer effects in DES |
| Chloride anion | Amber ions | Scaled to −0.8 | Consistent with choline scaling |
| Urea (reline) | GAFF2 | AM1-BCC, unscaled | Neutral molecule — scaling has no net effect |
| Glycerol (glyceline) | GAFF2 | AM1-BCC, unscaled | Neutral molecule |

### Charge scaling rationale (0.8)

Non-polarisable force fields with integer ionic charges (±1) for DES components systematically overstructure the solvent — producing viscosities an order of magnitude too high and self-diffusion coefficients that are far too low. This is because substantial charge transfer occurs between choline, chloride, and the hydrogen-bond donor (urea/glycerol) in the real system. The 0.8 scaling factor is the most widely validated correction for ChCl-based DES in the literature (Doherty & Bhatt; Perkins et al.), reproducing experimental densities, viscosities, and radial distribution functions.

### Solvent conditions

| Solvent | Composition | Validation target (density, 300 K) |
|----------------|----------------|-----------------------------------------|
| Water | Pure TIP3P | \~997 kg/m³ |
| Reline 20 wt% | ChCl:urea (1:2 mol) in TIP3P | \~1020–1040 kg/m³ |
| Glyceline 20 wt% | ChCl:glycerol (1:2 mol) in TIP3P | \~1030–1050 kg/m³ |

20 wt% aqueous DES was chosen over neat DES because: - More experimentally relevant (neat reline is too viscous for most bioprocessing) - Computationally faster (shorter equilibration, faster dynamics) - Consistent with typical DES dilutions used in protein extraction - Limitation: below eutectic composition, so characteristic DES H-bond cage is partially disrupted — this is acknowledged in the study design

### Simulation protocol

-   **Engine:** OpenMM 8.4+
-   **Integrator:** Langevin, 300 K, 1 ps⁻¹ friction, 2 fs timestep
-   **Barostat:** Monte Carlo, 1 bar
-   **Electrostatics:** PME, 1.0 nm cutoff
-   **Constraints:** HBonds (allows 2 fs step)
-   **Equilibration:** 1 ns water, 5 ns DES systems (NPT with position restraints → free NPT)
-   **Production:** 10 ns minimum per system
-   **Save interval:** 500 steps (1 ps)
-   **Hardware:** Apple M3 Max, CPU backend

### Analysis metrics

For each of the 9 production runs (3 peptides × 3 solvents):

1.  **Backbone RMSD** — equilibration assessment
2.  **Per-residue SASA** — protein residues only (Shrake-Rupley)
3.  **Motif SASA time series** — GGE/CME/YIY residues specifically
4.  **Solvent coordination number** — contacts within 0.35 nm of motif backbone O atoms
5.  **Backbone H-bond lifetimes** — ±3 residues around motif, using DCD save interval (not integrator step)
6.  **Bootstrap 95% CIs** — 2000 replicates on equilibrated portion (last ⅔ of trajectory)

Cross-solvent comparisons: ΔSASA, ΔHydration, ΔH-bond lifetime with propagated confidence intervals.

------------------------------------------------------------------------

## Project status

| Step | Description | Status |
|------------------|----------------------------------|---------------------|
| Phase 1 | L2-GGE pilot (water baseline) | Complete — corrected analysis notebook ready to run |
| 2.1 | Peptide construction (3 capped constructs) | Complete |
| 2.2 | Force field parameterisation (DES components) | Complete |
| 2.3 | Box construction (Packmol, 9 systems) | Complete |
| 2.4 | Equilibration + density validation | Complete |
| 2.5 | Production MD (9 runs) | Complete  |

======================================================================
PHASE 2 — COMPLETE VALIDATION SUMMARY
======================================================================
  System                  T (K)  ρ (kg/m³)   Target    Δρ%   Status
  ──────────────────────────────────────────────────────────────
  GGE_water               302.0       1008      997  +1.1%        ✓
  GGE_reline              300.3       1045     1030  +1.4%        ✓
  GGE_glyceline           300.6       1038     1040  -0.2%        ✓
  CME_water               301.7       1002      997  +0.5%        ✓
  CME_reline              300.6       1043     1030  +1.2%        ✓
  CME_glyceline           300.3       1032     1040  -0.7%        ✓
  YIY_water               302.2        996      997  -0.1%        ✓
  YIY_reline              300.7       1035     1030  +0.5%        ✓
  YIY_glyceline           300.3       1028     1040  -1.2%        ✓

  All systems validated from equilibration logs.
  
| 3.1 | Per-system MD analysis (9 systems) | Complete |
| 3.2 | Comparative analysis (deltas, effect sizes) | Complete |

======================================================================
PHASE 3 — COMPARATIVE MD ANALYSIS RESULTS
======================================================================

Analysis window: 1–10 ns (9,000 frames per system after 1 ns uniform
equilibration cutoff). All metrics with bootstrap 95% CIs (2,000
replicates).

### Motif SASA: DES vs water

  Peptide  Solvent      SASA_water  SASA_DES   ΔSASA     ΔSASA%   ES
  ─────────────────────────────────────────────────────────────────────
  GGE      Reline       0.450       0.920      +0.470    +104.3%  36.9
  GGE      Glyceline    0.450       0.894      +0.444     +98.6%  30.6
  CME      Reline       0.434       1.246      +0.812    +187.1%  50.4
  CME      Glyceline    0.434       1.492      +1.058    +243.6%  54.1
  YIY      Reline       0.593       1.304      +0.711    +119.8%  37.9
  YIY      Glyceline    0.593       1.963      +1.370    +230.8%  59.1

  SASA in nm². ES = effect size (|Δ| / pooled CI half-width; >2 =
  non-overlapping CIs). All comparisons ES > 30.

### Motif responsiveness ranking

  Rank  Motif  Mean ΔSASA%   Mean ES   Notes
  ─────────────────────────────────────────────────────────────────
  1     CME    +215.4%       52.2      Sulphur-containing residues
  2     YIY    +175.3%       48.5      Bulky aromatic residues
  3     GGE    +101.4%       33.7      Small, flexible residues

### Reline vs glyceline selectivity

  Motif  Glyceline − Reline ΔSASA   Direction
  ──────────────────────────────────────────────
  CME    +0.246 nm²                  Glyceline > Reline
  YIY    +0.659 nm²                  Glyceline >> Reline
  GGE    −0.026 nm²                  No significant preference

### Water coordination (contacts within 0.35 nm of motif backbone O)

  Peptide  Solvent      Water_w  Water_DES  ΔWater    ΔWater%
  ─────────────────────────────────────────────────────────────
  GGE      Reline       7.16     5.68       −1.48     −20.6%
  GGE      Glyceline    7.16     5.46       −1.70     −23.8%
  CME      Reline       5.89     4.72       −1.16     −19.7%
  CME      Glyceline    5.89     5.11       −0.77     −13.1%
  YIY      Reline       5.48     4.36       −1.11     −20.3%
  YIY      Glyceline    5.48     4.70       −0.77     −14.1%

### Key mechanistic observations

  - Urea is the only DES component with meaningful direct motif
    contact (0.32–0.38). Choline, chloride, glycerol are negligible.
    DES effect operates through solvent restructuring, not direct
    binding.
  - GGE in water has longest backbone H-bond lifetime (median 3 ps
    vs 2 ps elsewhere); events drop 26–48% in DES.
  - YIY: no backbone H-bonds in water or glyceline; 1,433 events
    in reline only (urea-stabilised transient conformations).
    
| 4.1 | In silico proteolysis | Not started |
| 4.2 | Bioactivity prediction | Not started |
| 4.3 | Molecular docking | Not started |
| 4.4 | ADMET prediction | Not started |
| 5 | Manuscript | Not started |



  Full results: analysis/Phase_3_Summary.md

------------------------------------------------------------------------

## Known issues and limitations

1.  **Phase 1 water baseline metrics are incorrect.** The original analysis used `g_idx=[0,1]` (MET, ASN) instead of the actual GGE motif at residues 8–10. Per-residue SASA included water molecules. The corrected analysis notebook fixes these issues but has not yet been rerun against the trajectory.

2.  **Sequence discrepancy in Phase 1 report.** The interim report stated the sequence as MNAVDHPHGGGEGRAPIGRKK (21 residues, with a typo at position 3). The actual simulated sequence was the full 43-mer MNPVDHPHGGGEGRAPIGRKKPATPWGYPALSKCFFFYYLNIQ.

3.  **H-bond lifetime timestep.** MDTraj trajectory objects report the DCD save interval as `traj.timestep`, not the integrator step. The corrected notebook uses the save interval explicitly (500 steps × 2 fs = 1 ps per frame).

4.  **No experimental DES density data at 20 wt%.** Literature values for neat reline (\~1.197 g/cm³ at 298 K) are well established, but 20 wt% aqueous mixture densities will need to be estimated and the validation tolerance set accordingly.

------------------------------------------------------------------------

## File inventory

### Reference documents

-   `Deep_Eutectic_Solvents_Manuscript_marked_Red.doc` — Durrani et al. original wet-lab manuscript
-   `peptides_t_grandis_homologs__1_.xlsx` — BLAST search results for 87 peptide sequences
-   `Torreya_Grandis_Motif_Homologs.pptx` — AlphaFold pTM score analysis of candidate proteins

### Workflow and analysis

-   `DES_Peptide_Project_Workflow.md` — Complete project workflow (this document's companion)
-   `corrected_water_baseline_analysis.ipynb` — Fixed analysis for Phase 1 pilot
-   `Step_2_1_Peptide_Construction.ipynb` — Build 3 capped peptide constructs

### Phase 3 analysis

-   `Step_3_1_Per_System_MD_Analysis.ipynb` — Per-system metrics for all 9 trajectories (~87 min runtime)
-   `Step_3_2_Comparative_Analysis.ipynb` — Cross-solvent deltas, effect sizes, motif ranking (<30 s runtime)
-   `analysis/Phase_3_Summary.md` — Complete Phase 3 findings summary
-   `analysis/per_system_metrics.csv` — All per-system metrics with bootstrap 95% CIs (36 rows × 27 columns)
-   `analysis/rmsd_equilibration_summary.csv` — Per-system RMSD equilibration cutoffs
-   `analysis/comparative_deltas.csv` — Cross-solvent deltas with propagated CIs and effect sizes
-   `analysis/cross_motif_summary.csv` — Motif responsiveness ranking
-   `analysis/figures/` — 11 diagnostic and publication figures

### Phase 1 pilot (L2-GGE)

-   `AF-A0A445BRA0-F1-model_v6.pdb` — AlphaFold structure (43 residues, *A. hypogaea*)
-   `GGE_MDS.ipynb` — Original simulation + analysis notebook (contains indexing errors)
-   `MD_Simulation_Interim_Report.pdf` — Water baseline results (metrics need recalculation)

------------------------------------------------------------------------

## Software environment

```         
Platform: Apple M3 Max (macOS, osx-arm64)
Python: 3.11 (Anaconda)
OpenMM: 8.4+
MDTraj: latest
AmberTools: required for tleap (conda install -c conda-forge ambertools)
Packmol: required for Step 2.3 (conda install -c conda-forge packmol)
```

------------------------------------------------------------------------

## Collaboration

This project involves collaboration with researchers at Zhejiang A&F University (Durrani et al. group). The computational study extends their experimental findings on DES-mediated protein extraction from *Torreya grandis*.
