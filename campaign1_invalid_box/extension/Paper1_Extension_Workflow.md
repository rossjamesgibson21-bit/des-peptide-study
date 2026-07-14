# Paper 1 — Computational Analysis Extension Workflow

## Molecular Dynamics Simulations Suggest Solvent-Dependent Modulation of Bioactive Motif Accessibility from Short Peptides of *Torreya grandis* in Aqueous Choline Chloride-Based DES Mixtures

**Author:** Ross Gibson
**Date:** June 2026
**Status:** Extension planning — awaiting execution

---

## Context

Phases 1–4 of the original workflow are complete. A first-draft manuscript has been shared with Rabia Durrani (Zhejiang A&F University), whose experimental dataset underpins this study. The extension work described here addresses three categories of limitation identified in the draft:

1. **Simulation sampling** — single 10 ns trajectories per system, insufficient for statistical robustness and vulnerable to kinetic trapping artefacts, particularly in glyceline.
2. **Coordination analysis scope** — restricted to backbone oxygen atoms; side-chain heteroatom contacts (Tyr –OH, Cys –SH, Met –SCH₃, Glu –COO⁻) were not assessed.
3. **Bioactivity prediction coverage** — limited to PeptideRanker, AntiCP 2.0, and BIOPEP-UWM; additional ML-based predictors (DeepMolecules) recommended by Rabia.

This document defines the extension phases, their inputs and outputs, dependencies on existing analyses, and connections to manuscript revision.

---

## Relationship to existing phases

| Existing phase | Description | Status | Extension dependency |
|---|---|---|---|
| Phase 1 | Pilot validation (L2-GGE water baseline) | Complete | None — superseded by Phase 2 |
| Phase 2 | Production simulations (9 systems, 10 ns each) | Complete | Phase E1 extends these trajectories; Phase E2 adds replicates |
| Phase 3 | MD analysis (SASA, coordination, H-bonds) | Complete | Phase E3 re-runs analysis on extended + replicate data |
| Phase 4 | Downstream pipeline (proteolysis, bioactivity, docking, ADMET) | Complete | Phase E4 adds DeepMolecules; docking/ADMET unchanged unless new fragments emerge |
| Draft manuscript | First draft shared with collaborator | Complete | Phase E5 revises manuscript with all new results |

---

## Extension Phase E1: Hardware Migration and Trajectory Extension

### E1.1 — Hardware benchmarking — COMPLETE

**Objective:** Establish per-system throughput on Apple M5 Max (128 GB) to calibrate runtime estimates for all subsequent phases.

**Environment:**
- OpenMM 8.5.1 (conda-forge, des-peptide environment)
- Python 3.11 (Anaconda)
- Available platforms: Reference, CPU, OpenCL (no native Metal backend in OpenMM 8.5.1; OpenCL uses Apple's OpenCL-to-Metal translation layer)
- Conda environment: `des-peptide`

**Benchmark protocol:** 5,000 steps (10 ps) per system after 500-step warmup and energy minimisation. Three representative systems tested on CPU and OpenCL platforms.

**Results (3 June 2026):**

| System | Atoms | CPU (ns/day) | OpenCL (ns/day) |
|---|---|---|---|
| GGE_water | 2,838 | 79.8 | 717.2 |
| CME_reline | 2,898 | 76.3 | 705.9 |
| YIY_glyceline | 2,876 | 77.8 | 711.5 |

**Key observations:**
- OpenCL throughput is consistent across solvent classes (~706–717 ns/day), as expected for similar atom counts on GPU
- OpenCL is ~9× faster than CPU on this hardware
- Compared to M3 Max CPU estimates (30–60 ns/day), the M5 Max OpenCL backend represents a ~12–20× throughput improvement
- All systems have similar atom counts (~2,850–2,900); DES systems are not slower than water on GPU, unlike the CPU-bound M3 Max where DES viscosity-related overhead was observed
- **Production platform: OpenCL** (all subsequent runs)

**Reference throughput for runtime estimation: 710 ns/day** (conservative average across all systems)

### E1.2 — Trajectory extension

**Objective:** Extend all 9 existing production trajectories to a target length sufficient to address kinetic trapping and autocorrelation concerns.

**Target length:** Minimum 100 ns per system. The decision on whether to extend further (to 500 ns) should be made after examining autocorrelation times from the 100 ns data. If glyceline τ_int values remain >5 ns at 100 ns total, further extension is warranted.

**Tasks:**
1. Resume each production trajectory from the final state of the Phase 2 run
2. Maintain identical simulation parameters (integrator, barostat, PME cutoff, save interval)
3. Save interval: 500 steps (1 ps), consistent with existing trajectories
4. Concatenate new trajectory frames with existing 10 ns data post-production
5. Validate continuity: check density, temperature, and RMSD across the junction

**Execution order:** Sequential, ascending estimated runtime (consistent with established practice):
1. GGE_water
2. CME_water
3. YIY_water
4. GGE_reline
5. CME_reline
6. YIY_reline
7. GGE_glyceline
8. CME_glyceline
9. YIY_glyceline

**Calibrated runtime estimate (from E1.1 benchmark, 710 ns/day on OpenCL):**
- Per system: 90 ns extension (10→100 ns) ≈ 3.0 hours
- Total for 9 systems: 810 ns ≈ 27 hours (~1.1 days)

**Inputs:** Final state files from Phase 2; validated M5 Max environment (OpenCL platform)
**Outputs:** 9 extended trajectories (≥100 ns each); concatenated .dcd files; per-system logs

**Commit strategy:** GitHub commit after each completed system, not in batch.

**Key decision point:** After all 9 extensions are complete, re-evaluate autocorrelation times. If glyceline systems still show τ_int > 5 ns, the trajectory target increases. This decision gates Phase E3.

---

## Extension Phase E2: Independent Replicates

**Objective:** Generate ≥3 independent velocity replicates per system to establish statistical robustness beyond block-bootstrap correction.

**Design:**
- Each replicate starts from the same equilibrated coordinates but with different initial velocities (drawn from the Maxwell-Boltzmann distribution at 300 K using a different random seed)
- Production length per replicate: matched to the extended trajectory length from Phase E1 (≥100 ns)
- Total new production runs: 9 systems × 2 additional replicates = 18 runs (the extended Phase E1 trajectory serves as replicate 1)

**Tasks:**
1. For each of the 9 systems, generate 2 new velocity assignments from the post-equilibration coordinates
2. Run production MD with identical parameters to Phase E1
3. Apply same validation checks (density, temperature, RMSD)

**Execution order:** System-by-system (all replicates of one system before moving to the next), ascending estimated runtime.

**Calibrated runtime estimate (from E1.1 benchmark, 710 ns/day on OpenCL):**
- Per replicate: 100 ns ≈ 3.4 hours
- Total for 18 replicates: 1,800 ns ≈ 2.5 days

**Inputs:** Post-equilibration coordinates from Phase 2; M5 Max (OpenCL platform)
**Outputs:** 18 new trajectories; per-system logs; validation summaries

**Commit strategy:** GitHub commit after each completed replicate.

**Note on starting conformations:** All replicates use extended starting conformations, consistent with the original Phase 2 simulations. Compact-start replicates have been deferred (see Resolved Design Decisions, point 2). After E2 completion, inter-replicate SASA convergence will be assessed. If glyceline replicates show poor convergence (suggesting kinetic trapping from the extended start), compact-start replicates will be run for the 3 glyceline systems only — adding 3–6 additional runs rather than redesigning the full replicate set.

---

## Extension Phase E3: Extended MD Analysis

### E3.0 — Preliminary RDF analysis (existing data) — COMPLETE

**Objective:** Compute atom-pair RDFs from the existing 10 ns trajectories to establish approximate coordination cutoffs and scope the side-chain analysis before extended data is available.

**Completed:** 4 June 2026
**Notebook:** `extension/analysis/preliminary_rdfs/E3_0_Preliminary_RDF_Analysis.ipynb`
**Outputs:** `extension/analysis/preliminary_rdfs/preliminary_rdfs_all_systems.png`, `preliminary_rdf_cutoffs.csv`

**Findings:**

| Atom pair | First minimum (nm) | Δ from 0.35 nm | RDF quality | Phase 3 impact |
|---|---|---|---|---|
| Water O → backbone O | 0.333–0.343 | −0.07 to −0.17 Å | Clean, well-defined | Negligible — slight overcounting |
| Urea N → backbone O | 0.358–0.412 | +0.08 to +0.62 Å | Moderate | Mild undercounting for CME/YIY |
| Urea O → backbone O | 0.538–0.543 | +1.88 to +1.93 Å | Clear first shell | **Substantial** — entire coordination shell missed at 0.35 nm |
| Choline N → backbone O | 0.387–0.403 | +0.37 to +0.53 Å | Diffuse, no clean trough (g(r_min) = 1.4–3.1) | No structured first shell; Phase 3 "negligible" finding confirmed |
| Chloride → backbone O | N/A | N/A | Flat, featureless | Confirmed negligible at any cutoff |
| Glycerol O → backbone O | N/A | N/A | No first-shell peak | **Key finding:** confirms no structured glycerol–backbone coordination; indirect mechanism hypothesis strengthened |

**Interpretation:**

1. **Urea O is the actionable finding.** The 0.35 nm cutoff captured urea N contacts (~0.36–0.41 nm first shell) but completely missed urea carbonyl O coordination (~0.54 nm first shell). Using pair-specific cutoffs will increase total urea coordination counts, strengthening the Phase 3 finding that urea is the dominant direct-contact DES component.

2. **Glycerol result is cutoff-independent.** No structured coordination at any distance — the "indirect mechanism" hypothesis (Section 4.2 of the draft) is supported by an independent analytical method. This is not a cutoff artefact; glycerol genuinely does not form first-shell contacts with backbone O atoms.

3. **Chloride and choline show no structured coordination.** Phase 3 conclusions for these components are confirmed.

4. **Qualitative conclusions from Phase 3 are robust.** The 0.35 nm cutoff was imperfect but did not produce misleading qualitative results. Pair-specific cutoffs will refine the quantitative coordination numbers without changing the mechanistic narrative.

**Preliminary cutoffs for E3.1/E3.2 (to be confirmed from 100 ns extended data):**

| Atom pair | Recommended cutoff (nm) | Notes |
|---|---|---|
| Water O → backbone O | 0.34 | Confirmed from RDF first minimum |
| Urea N → backbone O | 0.40 | Slight system dependence; refine from extended data |
| Urea O → backbone O | 0.54 | Consistent across all reline systems |
| Choline N → backbone O | — | No structured shell; report as unstructured |
| Chloride → backbone O | — | No structured shell; report as unstructured |
| Glycerol O → backbone O | — | No structured shell; report as unstructured |

**Caveat:** These RDFs were computed from 10 ns trajectories with small reference groups (3 backbone O atoms per system). Low-abundance species (choline, chloride) and systems with long autocorrelation times (glyceline) may show improved RDF convergence in the extended trajectories. The definitive cutoffs should be derived from the 100 ns data in Phase E3.1.

### E3.1 — Replicate-aware SASA and coordination analysis

**Objective:** Re-run the Phase 3 analysis pipeline on the extended and replicated trajectory data, with replicate-aware statistical treatment.

**Tasks:**
1. Apply equilibration cutoff to each trajectory independently (RMSD-based, as in Phase 3)
2. Compute motif SASA time series for all 27 trajectories (9 systems × 3 replicates)
3. Compute per-replicate means, then inter-replicate means and standard errors
4. Report both within-trajectory block-bootstrap CIs and inter-replicate standard errors
5. Re-compute autocorrelation times (τ_int) for extended trajectories — assess whether longer runs resolve the glyceline sampling concern
6. Re-compute all cross-solvent deltas (ΔSASA) with replicate-aware statistics

**Statistical framework:**
- Primary inference: inter-replicate mean ± SEM (n=3 per condition)
- Secondary: within-trajectory block-bootstrap CIs (as in Phase 3, but on longer data)
- Effect sizes: re-computed from replicate means
- The GGE glyceline comparison (corrected ES = 1.18 in Phase 3) is the key test case — does it achieve significance with replicates?

**Inputs:** 27 trajectories from Phases E1 + E2; existing Phase 3 analysis notebooks as templates
**Outputs:** Updated per_system_metrics.csv (expanded to 27 rows); updated comparative_deltas.csv; updated figures

### E3.2 — Side-chain coordination analysis

**Objective:** Extend the solvent coordination analysis from backbone oxygen atoms to side-chain heteroatoms, addressing the limitation identified in manuscript Sections 3.4 and 4.2.

**Target atoms by motif:**
- **GGE:** Glu side-chain carboxyl oxygens (OE1, OE2)
- **CME:** Cys thiol sulphur (SG), Met thioether sulphur (SD), Glu carboxyl oxygens (OE1, OE2)
- **YIY:** Tyr phenolic oxygen (OH) — both tyrosine residues

**Coordination partners (within cutoff):**
- Water oxygen (OH2)
- Choline nitrogen (N)
- Chloride (Cl⁻)
- Urea nitrogen and oxygen (reline systems)
- Glycerol oxygen atoms (glyceline systems)

**Coordination cutoff protocol (RDF-derived, pair-specific):**

The uniform 0.35 nm cutoff used in Phase 3 is calibrated for water-water first-shell distances and likely undercounts contacts from bulkier DES components (glycerol effective radius ~0.26 nm, choline ~0.35 nm). Phase E3.2 replaces this with pair-specific cutoffs derived from radial distribution functions:

1. **Preliminary RDFs (from existing 10 ns data):** Compute atom-pair RDFs for each DES component atom type vs backbone oxygen atoms across all existing trajectories. This provides approximate cutoff values and a sanity check before extended data is available. Note: RDFs may be noisy for low-abundance species (choline, chloride) at 10 ns.

2. **Definitive RDFs (from 100 ns extended data):** Re-compute all RDFs from the extended trajectories (E1.2), where statistics will be substantially better. Place the coordination cutoff at the first RDF minimum for each atom pair separately (e.g. urea-O to backbone-O; glycerol-O to backbone-O; water-O to backbone-O).

3. **Side-chain RDFs:** Apply the same protocol to side-chain heteroatom pairs — compute RDFs for each DES component vs each target side-chain atom (Tyr OH, Cys SG, Met SD, Glu OE1/OE2), and derive pair-specific cutoffs.

4. **Sensitivity analysis:** For each atom pair, recompute coordination numbers at the RDF-derived cutoff and at ±0.05 nm. Report whether qualitative conclusions (e.g. glycerol coordination negligible vs substantial) are robust to cutoff variation.

**Tasks:**
1. Compute preliminary RDFs from existing 10 ns data (backbone O pairs only — scoping step)
2. Identify atom indices for all target side-chain heteroatoms in each system topology
3. Compute definitive RDFs from extended trajectories for all backbone and side-chain atom pairs
4. Derive pair-specific cutoffs at first RDF minimum
5. Compute per-frame coordination numbers for each heteroatom–partner combination using pair-specific cutoffs
6. Run sensitivity analysis at ±0.05 nm around each cutoff
7. Report time-averaged coordination with bootstrap CIs and inter-replicate statistics
8. Compare backbone vs side-chain coordination profiles across solvents
9. Assess whether direct DES–side-chain contacts are present (particularly glycerol–Tyr OH and urea–Cys SH)

**Key question this answers:** Is the "indirect solvent restructuring" hypothesis (Section 4.2) supported or challenged when side-chain contacts are included? If glycerol makes substantial direct contact with Tyr or Cys side chains, the mechanism may be partially direct rather than purely indirect.

**Inputs:** 27 trajectories; topology files with atom name mappings
**Outputs:** Side-chain coordination tables; decomposed coordination bar charts; updated mechanistic interpretation

### E3.3 — Backbone hydrogen bond analysis (extended data)

**Objective:** Re-run H-bond lifetime analysis on extended trajectories to assess whether the YIY reline anomaly (1,433 events in reline, zero in water/glyceline) persists with longer sampling.

**Tasks:**
1. Compute backbone H-bond events and lifetime distributions for all 27 trajectories
2. Report per-replicate statistics; assess inter-replicate consistency of the YIY anomaly
3. If the anomaly is replicate-dependent (i.e. appears in some replicates but not others), this changes the mechanistic interpretation

**Inputs:** 27 trajectories
**Outputs:** Updated H-bond lifetime distributions; replicate consistency assessment

---

## Extension Phase E4: Bioactivity Prediction Expansion

**Objective:** Incorporate DeepMolecules (deepmolecules.org) predictions alongside the existing PeptideRanker/AntiCP/BIOPEP stack, per Rabia's recommendation.

**Tasks:**
1. Submit the 9 candidate peptide sequences to DeepMolecules
2. Record predicted bioactivity classes and confidence scores
3. Compare DeepMolecules predictions with existing PeptideRanker scores and BIOPEP motif matches
4. Assess whether the CME > YIY > GGE antioxidant ranking is supported by the additional predictor
5. Integrate results into the bioactivity summary table (phase4_integrated_summary.csv)

**Independence from simulation phases:** This phase has no dependency on Phases E1–E3 and can be executed in parallel.

**Inputs:** 9 candidate sequences from Phase 4; DeepMolecules web server
**Outputs:** Updated bioactivity_summary.csv; comparison table across predictors

---

## Extension Phase E5: Manuscript Revision

**Objective:** Revise the first-draft manuscript to incorporate all extension results and address structural issues identified in internal review.

**Structural issues to address (from internal review):**
1. **Figure numbering:** Figure 5 (backbone H-bond lifetime distributions) appears in the Discussion (Section 4.2) but is referenced from Section 3.5. Relocate to Results or renumber.
2. **Circularity caveat placement:** The caveat that antioxidant scoring weights were informed by the Durrani et al. dataset currently surfaces only in Section 3.7. State earlier, in Methods (Section 2.8).
3. **Docking limitations prominence:** Acknowledgement of uncontrolled zinc ion, protonation state, and co-crystallised ligand handling is buried in Section 2.9. Consider foregrounding in a limitations paragraph or moving to Discussion Section 4.5.
4. **Tonal alignment:** The email to Rabia described antioxidant scoring concordance as "perfectly aligned"; the manuscript is appropriately more circumspect. Ensure consistency throughout.

**Content updates required:**
1. Update all SASA, coordination, and H-bond results with extended trajectory + replicate data
2. Add side-chain coordination results and revised mechanistic discussion
3. Add DeepMolecules bioactivity predictions to Section 3.7
4. Revise statistical framework description (Section 2.6) to reflect replicate-aware analysis
5. Update limitations section (4.5) — some limitations may be resolved (trajectory length, single-trajectory statistics), new ones may emerge
6. Update simulation protocol (Section 2.5) with extended trajectory lengths and replicate design
7. Update hardware description if M5 Max is used for extension work
8. Compile reference list (currently [REF] placeholders)

**Inputs:** All outputs from Phases E1–E4; existing first-draft manuscript
**Outputs:** Revised manuscript; updated supplementary materials; updated figures

---

## Critical path and dependencies

```
Phase E1.1 (benchmark) ✓ COMPLETE
    │
    ▼
Phase E1.2 (extend 9 trajectories to 100 ns) — ~27 hours
    │
    ▼
CHECKPOINT: re-compute τ_int — is 100 ns sufficient?
    │                                          
    ├── Yes ──────────────────┐
    │                         │
    ▼                         ▼
Phase E2 (18 replicates     Phase E4 (DeepMolecules)  ← independent
 at 100 ns) — ~61 hours      │
    │                         │
    ▼                         │
Phase E3.1–E3.3              │
 (extended analysis)          │
    │                         │
    ├─────────────────────────┘
    ▼
Phase E5 (manuscript revision)

Phase E3.0 (preliminary RDFs) ✓ COMPLETE — confirms pair-specific cutoffs needed for urea O
Phase E4 (DeepMolecules)      ← independent, can run any time
```

The critical path runs E1 → checkpoint → E2 → E3 → E5. Phases E3.0 and E4 are off the critical path. If the checkpoint indicates 100 ns is insufficient for glyceline, E1.2 continues (extending to 200–500 ns) before E2 begins.

---

## Estimated timeline (calibrated from E1.1 benchmark: 710 ns/day on M5 Max OpenCL)

| Phase | Description | Estimated duration | Dependencies |
|---|---|---|---|
| E1.1 | Hardware benchmarking | **COMPLETE** (3 June 2026) | — |
| E1.2 | Trajectory extension (9 × 90 ns new) | ~1.1 days (~27 hours) | E1.1 |
| E3.0 | Preliminary RDF analysis (existing data) | **COMPLETE** (4 June 2026) | — |
| Checkpoint | Autocorrelation analysis of extended data | 1–2 days | E1.2 |
| E2 | Independent replicates (18 × 100 ns) | ~2.5 days (~61 hours) | Checkpoint passed |
| E3.1–E3.3 | Extended analysis (SASA, RDF, coordination, H-bonds) | 1–2 weeks | E1.2 + E2 complete |
| E4 | DeepMolecules predictions | 1–2 days | None |
| E5 | Manuscript revision | 2–4 weeks | E3 + E4 complete |
| **Total** | | **~5–8 weeks** | |

**Compute time summary:** Total MD production (E1.2 + E2) ≈ 2,610 ns ≈ 3.7 days wall-clock. The rate-limiting step is now analysis and manuscript revision, not simulation compute. If extension to 500 ns proves necessary, total compute increases to ~19 days — still feasible on the M5 Max without HPC access.

**Note:** All extensions (E1.2) complete before any replicates (E2) begin — this preserves the checkpoint logic. If the checkpoint reveals glyceline requires >100 ns, the replicate target length is adjusted before committing to 18 runs. Phase E4 can be completed at any time during E1.2 or E2.

---

## File structure (extension)

```
des-peptide-study/
├── ... (existing structure unchanged)
├── extension/
│   ├── Paper1_Extension_Workflow.md          ← this document
│   ├── benchmarks/
│   │   ├── benchmark_m5max.py               ← E1.1 benchmark script
│   │   └── m5_max_throughput.csv
│   ├── trajectories/
│   │   ├── GGE_water_ext/                   ← extended trajectory + 2 replicates
│   │   ├── GGE_water_rep2/
│   │   ├── GGE_water_rep3/
│   │   ├── ... (27 directories total)
│   ├── analysis/
│   │   ├── preliminary_rdfs/                ← E3.0: from existing 10 ns data
│   │   ├── extended_per_system_metrics.csv
│   │   ├── replicate_summary.csv
│   │   ├── sidechain_coordination/
│   │   ├── extended_hbond_analysis/
│   │   └── figures/
│   ├── bioactivity/
│   │   └── deepmolecules_predictions.csv
│   └── manuscript/
│       └── revision_notes.md
```

---

## Resolved design decisions

The following decisions were resolved during planning (June 2026):

1. **Target trajectory length:** 100 ns initial target. After all 9 extensions are complete (E1.2), re-compute τ_int values. If glyceline systems still show τ_int > 5 ns, extend further before committing to replicates. This is a checkpoint, not a pre-commitment to 500 ns.

2. **Compact-start replicates:** Deferred. All 27 trajectories use extended starting conformations with different velocity seeds. After E2 is complete, assess inter-replicate convergence. If glyceline replicates diverge (suggesting kinetic trapping), run compact-start replicates for glyceline systems only at that point. This avoids speculative compute while keeping the option open.

3. **Execution order:** Sequential — all 9 extensions first (E1.2), then autocorrelation checkpoint, then all 18 replicates (E2). This preserves the checkpoint logic: if 100 ns proves insufficient for glyceline, the replicate target length can be adjusted before committing to 18 runs at an inadequate length. Within each phase, systems are run sequentially in ascending estimated runtime (water → reline → glyceline).

4. **RDF-derived coordination cutoffs:** Yes — pair-specific cutoffs from RDF first minima, with sensitivity analysis at ±0.05 nm. **E3.0 preliminary results (4 June 2026) confirmed this is necessary:** urea carbonyl O has its first-shell minimum at ~0.54 nm, nearly 0.2 nm above the 0.35 nm default — the Phase 3 analysis missed this entire coordination shell. Water O first minima (0.333–0.343 nm) are close to the default. Glycerol O and chloride show no structured coordination at any distance, confirming Phase 3 findings independently. Definitive cutoffs will be derived from the 100 ns extended trajectories in Phase E3. The same protocol applies to both backbone oxygen and side-chain heteroatom coordination analyses.

5. **Docking re-run:** Conditional. The existing Phase 4.3 docking results stand unless the extended analysis produces new candidate fragments from revised proteolysis predictions. If new fragments emerge, re-dock those fragments only against the same three targets (Keap1, ACE, COX-2) using the established CB-Dock2 protocol.
