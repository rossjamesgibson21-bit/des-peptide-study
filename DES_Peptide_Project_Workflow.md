# DES Peptide Simulation Project — Complete Workflow

## Solvent-Dependent Modulation of Bioactive Motif Accessibility in *Torreya grandis* Peptides: A Molecular Dynamics Study in Water and Choline Chloride-Based Deep Eutectic Solvents

**Author:** Ross Gibson
**Created:** April 2026 · **Revised:** June 2026

> **Revision note (June 2026).** This overview has been updated to reflect the post-setup
> state of the project. Production was carried out at 100 ns rather than the original 10 ns,
> on an Apple M5 Max (OpenCL) rather than the M3 Max (CPU); the statistical treatment now
> uses autocorrelation-corrected confidence intervals; and a new **Phase E** documents the
> sampling-adequacy checkpoint and the extended-production design that follows from it. The
> setup phases (1–2) are retained as the build record. This document covers **Paper 1 (the MD
> study)**; the DFT multi-scale framework (Paper 2) is a separate, deferred track.

---

## Study overview

This study uses molecular dynamics (MD) simulations to investigate how two choline chloride-based deep eutectic solvents (reline and glyceline) modulate the solvent accessibility and structural dynamics of three bioactive peptide motifs (GGE, CME, YIY) identified in *Torreya grandis* seed proteins by Durrani et al. The aim is to provide a computational rationale for the enhanced bioactive peptide yield observed experimentally when DES is used for protein extraction.

### Peptide constructs

| Motif | Source peptide (Durrani et al.) | Construct sequence | Residues | Motif position (0-indexed) |
|-------|--------------------------------|-------------------|----------|---------------------------|
| GGE | NTDEEEGGEIVFGGVDPK (Peptide 6) | Ace-EEEGGEIVF-NMe | 9 | 3–5 |
| CME | LYQCMEFVR (Peptide 5) | Ace-LYQCMEFVR-NMe | 9 | 3–5 |
| YIY | NPYIYK (Peptide 83) | Ace-NPYIYK-NMe | 6 | 2–4 |

### Solvent conditions

| Solvent | Composition | Charge scaling | Notes |
|---------|------------|---------------|-------|
| Water | Pure TIP3P | N/A | Control |
| Reline (20 wt%) | ChCl:urea (1:2 molar) in TIP3P water | 0.8 on choline (+0.8) and chloride (−0.8) | Most effective DES in Durrani et al. |
| Glyceline (20 wt%) | ChCl:glycerol (1:2 molar) in TIP3P water | 0.8 on choline (+0.8) and chloride (−0.8) | Different H-bond character |

**Systems: 9** (3 peptides × 3 solvents).

---

## Phase 1: Pilot validation — COMPLETED

Water baseline MD of full-length L2-GGE protein; protocol validated (Amber14SB + TIP3P, Langevin NPT, 2 fs); corrected analysis notebook produced (GGE indexing fixed to residues 8–10).

---

## Phase 2: Production setup — COMPLETED

Steps 2.1–2.4 (peptide construction with Ace/NMe caps, GAFF2 + AM1-BCC parameterisation with 0.8 charge scaling on choline/chloride, Packmol box construction at 20 wt% DES, minimisation + NVT/NPT equilibration) were completed as originally specified. The nine equilibrated systems reside in `systems/{name}/` with Amber `.prmtop` topologies and `_prod.dcd` trajectories.

> **Coordinate-sourcing note.** The per-system `_equilibrated.pdb` files are mis-ordered
> relative to the `.prmtop`. Production runs therefore source starting coordinates from
> **frame 0 of `{system}_prod.dcd`** via `md.load_frame(..., top=prmtop)`, which enforces the
> topology atom order and yields physically sane start energies (~−31 to −37k kJ/mol). A
> start-PE guard (|PE| < 1e7) protects every run.

### Step 2.5 — Production MD (as run)

- Integrator: Langevin, 300 K, friction 1 ps⁻¹, 2 fs timestep
- Barostat: Monte Carlo, 1 bar, 300 K
- Nonbonded: PME, 1.0 nm cutoff · Constraints: HBonds
- Trajectory save: every 500 steps (1 ps) · OpenCL default precision
- **Run length: 100 ns** per trajectory (revised up from the original 10 ns)

Production and all extended sampling are governed by **Phase E** below.

---

## Phase 3: Analysis

Metrics per trajectory (MDTraj): backbone RMSD; per-residue and motif SASA (Shrake-Rupley); solvent coordination number (RDF-derived, pair-specific cutoffs with ±0.05 nm sensitivity — note the urea carbonyl-O first shell sits near 0.54 nm, missed by the default 0.35 nm; glycerol-O shows no structured coordination with backbone O); backbone H-bond lifetimes; secondary structure where relevant.

**Statistical treatment (revised).**
- Equilibration detected by an excursion-tolerant criterion (≥95% of remaining frames in-band, ±15% tolerance, capped at 50% of trajectory length) rather than a fixed last-⅔ rule.
- **Confidence intervals are autocorrelation-corrected.** Raw bootstrap CIs treat correlated frames as independent and are misleadingly tight; CIs are widened by the integrated autocorrelation time τ_int (see Phase E), reported alongside effective sample size N_eff.
- Cross-solvent (ΔSASA, ΔHydration, ΔH-bond lifetime, effect sizes) and cross-motif comparisons as before, on adequately-sampled ensembles only.

---

## Phase 4: Downstream computational analysis — COMPLETED

In-silico proteolysis (Step 4.1), bioactivity prediction (BIOPEP-UWM, PeptideRanker, AntiCP; CME confirmed as a catalogued antioxidative peptide and the only candidate above the PeptideRanker 0.5 threshold), molecular docking (CB-Dock2 / AutoDock Vina against Keap1, ACE, COX-2), and ADMET profiling (RDKit) are complete, with results integrated in `phase4_integrated_summary.csv`. Docking re-runs are conditional on new candidate fragments emerging from the extended analysis.

---

## Phase E: Sampling adequacy and extended production

This phase was added after the 100 ns campaign revealed that the slow conformational modes of these peptides are not adequately sampled at 100 ns. Its purpose is to size the sampling correctly before committing the full replicate set, and to do so with explicit convergence diagnostics rather than by assumption.

### E1.1 — Hardware benchmark
Apple M5 Max (128 GB, 40 GPU cores, OpenCL): ~693–710 ns/day cool-start per system. Back-to-back runs incur a measured ~21% second-run thermal penalty (second run ~480 ns/day), so a two-run overnight batch is ~9 h, not ~8 h. A full power-off cool-down restores the cool-start rate.

### E1.2 — Diagnostic first-replicate campaign
Independent 100 ns trajectories from the equilibrated frame-0 coordinates with deterministic per-run velocity seeds. Completed and verified (100,000 frames each): GGE/CME/YIY in water; GGE and CME in glyceline. These runs are **diagnostic** — under the extended design they are not retained as production data, but serve as the source of diverse starting conformations for Phase E2.

### E1.3 — Autocorrelation checkpoint
A sampling-adequacy gate. The integrated autocorrelation time τ_int of the reported observable (motif SASA) and of the radius of gyration (the slow conformational mode) is estimated by FFT autocorrelation with Sokal automated windowing (τ_int = 1 + 2 Σ ρ(k), windowed at W ≥ 5·τ_int(W)); the estimator was validated against AR(1) processes to ~1%. Reported per trajectory: τ_int (ns), N_eff = N/τ_int, CI inflation √τ_int, and a reliability flag (trustworthy only when T ≥ 50·τ_int).

**Motif-SASA τ_int (single replicate):**

| System | τ_int (ns) | N_eff | Reliable (T ≥ 50τ) |
|--------|-----------|-------|--------------------|
| CME_water | 1.51 | 59 | yes |
| YIY_water | 2.85 | 32 | no |
| GGE_water | 4.51 | 20 | no |
| CME_glyceline | 5.30 | 17 | no |
| GGE_glyceline | 9.99 | 9 | no (soft lower bound) |

**Observations.** Even water under-samples at 100 ns for the slow modes. Glyceline roughly doubles-to-triples the SASA autocorrelation time relative to water (CME 1.5→5.3 ns; GGE 4.5→10.0 ns), consistent with a more viscous, more structured solvent slowing peptide conformational dynamics — in keeping with an indirect, solvent-restructuring mechanism. GGE_glyceline's diagnostic plots show bidirectional slow wandering and a non-plateauing τ_int(W) curve; its 9.99 ns is a lower bound that 100 ns cannot resolve (N_eff = 9). Radius-of-gyration τ_int is longer throughout (4.2–10.7 ns even in water), evidence that the compact↔extended interconversion behind the kinetic-trapping concern is genuinely slow.

**Verdict: EXTEND.** Both glyceline systems exceed the 5 ns threshold; even the fastest motif (CME) is over it in glyceline. 100 ns is insufficient for the slow systems.

### E2 — Extended-production design (current operative plan)

The slow, low-N_eff regime is one in which ergodic coverage of conformational basins, not within-trajectory decorrelation, is the binding constraint. The design therefore favours replicate breadth from diverse starts over single-trajectory length.

- **Duration:** 200 ns per replicate (≈19 τ_int post-equilibration even for the worst case; ~overnight footprint). Single 500 ns trajectories are deprecated for this regime because one trajectory cannot test ergodicity.
- **Diverse starts:** replicates are seeded from conformations spanning the Rg/SASA range of the E1.2 diagnostic trajectories (compact, extended, and intermediate), not from velocity reseeding of a single frame 0. This is the productive reuse of the diagnostic runs.
- **Convergence as diagnostic:** extended-start and compact/diverse-start sub-ensembles are compared; agreement of their SASA/Rg distributions demonstrates ergodicity and directly answers the kinetic-trapping critique. Divergence indicates trapping that length will not fix → escalate to enhanced sampling (Paper 2 / methods).
- **Target:** pooled N_eff ≥ 50–100 per system; replicate count adaptive to reach it. Equilibration discard scaled to τ_int (~20 ns for glyceline).
- **Tiering (leaner version, in progress):** fast water systems kept light (existing data + a few diverse-start replicates — CME_water is already N_eff 59); glyceline given the full treatment; reline probed with one replicate to establish its τ_int before committing its set.
- **Sequencing:** GGE_glyceline first (worst case) — run its diverse-start replicates and check extended/compact convergence early, before committing compute to the rest.
- **Compute (flagged as estimate):** leaner version ~11–13 days; full nine-system version ~16–18 days of machine time.

Runs use `extension/E2_1_Extended_Production.ipynb`, writing to `trajectories_extended/{system}_{NS}ns_{start}_r{rep}/`, kept separate from the 100 ns diagnostic runs.

---

## Phase 5: Manuscript preparation

Structure as previously planned (introduction → methods → results → discussion → conclusion), with the mechanistic findings (coordination structure, H-bond lifetimes) as the analytical backbone and the accessibility comparison reported with autocorrelation-corrected CIs. Outstanding structural fixes to the shared draft: missing Table 5 (solvent coordination) and a Section 3.4 paragraph truncation.

**Key limitations to address (updated).**
- Short capped fragments rather than full proteins in cellular context.
- 20 wt% aqueous DES is below eutectic composition.
- Non-polarisable force field with empirical 0.8 charge scaling for DES.
- **Sampling of slow conformational modes:** addressed directly via τ_int analysis and the extended diverse-start replicate design; results reported with honest N_eff-based bounds rather than raw-frame statistics.
- Computational predictions require experimental validation.

---

## Tools & environment

- **Hardware:** Apple M5 Max MacBook Pro (128 GB, 40 GPU cores, OpenCL / Metal 4).
- **Environment:** `des-peptide` conda env — OpenMM 8.5.1, MDTraj 1.11.1, AmberTools (classic solver).
- **Repositories:** `des-peptide-study` (MD, private); `des-peptide-dft` (Paper 2, planned).
- **Force fields:** ff14SB (peptides), GAFF2 + AM1-BCC (DES components), choline +0.8 / chloride −0.8.

---

## File structure

```
des-peptide-study/
├── systems/                         # 9 equilibrated systems (.prmtop, _prod.dcd)
├── extension/
│   ├── E1_2_Trajectory_Extension.ipynb            # single 100 ns run
│   ├── E1_2_Trajectory_Extension_Batch.ipynb      # batched 100 ns runs
│   ├── E1_3_Autocorrelation_Checkpoint.ipynb      # tau_int sampling gate
│   ├── E2_1_Extended_Production.ipynb             # 200 ns extended-campaign runs
│   ├── trajectories/                              # 100 ns diagnostic runs
│   ├── trajectories_extended/                     # 200 ns diverse-start replicates
│   └── analysis/
│       └── autocorrelation_checkpoint/            # tau_int CSV + diagnostic plots
├── analysis/                        # Phase 3 notebooks + figures
├── downstream/                      # Phase 4 (proteolysis, bioactivity, docking, admet)
└── manuscript/
```
