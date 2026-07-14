# DES-PEP-CESD-001

Consolidated computational campaign — corrected MD, enhanced sampling and targeted DFT for
GGE, CME and YIY peptides in water, reline and glyceline.

**Protocol:** DES-PEP-CESD-001 v1.0 (14 July 2026). Study owner: Ross Gibson.

## Central design

Replicated, reweighted free-energy comparison (well-tempered metadynamics) from independently
prepared corrected systems. Conventional MD qualifies the model; DFT is targeted mechanistic
validation. A valid negative or inconclusive solvent comparison is a **successful** outcome
(§1.2). The campaign fails only if the model, sampling, reweighting, reproducibility or
provenance chain cannot support the prespecified inference.

## Execution arrangements (agreed at Phase 0)

- **Production MD/metadynamics** runs on the study owner's **M5 Max MacBook Pro** (OpenMM GPU;
  the §19 budget was written for this platform). The analysis/build environment here does **not**
  run production MD.
- **Gate sign-offs and the solvent-code key** are held by **real human reviewers / data steward**.
  This repo's notebooks prepare each gate's **evidence package** and enforce phase sequencing;
  they do **not** sign gates and do **not** hold the blinding key.

See `governance/roles_and_arrangements.yaml`.

## Per-phase workflow

Each phase is a self-contained notebook that loads the frozen inputs, performs the phase's
verification/analysis, and emits a hashed evidence package + a pre-populated (unsigned) gate
certificate. **Once the human reviewers sign the gate PASS, that phase is committed and tagged
`gate-NN`.** The signed release output of phase *n* is the required input to phase *n+1*
(§4 failproof rule).

```
phases/phaseNN_*/         phaseNN_*.py (jupytext source) + phaseNN_*.ipynb (executed) + outputs/
governance/               roles, blinding charter, reviewer checklists, deviation log, gate certificates
sap/                      statistical analysis plan (shell -> frozen at Gate 7)
registries/ runs/         identifier scheme, run registry
manifests/                per-gate SHA-256 evidence manifests
models/ trajectories/ enhanced_sampling/ qc/ dft/ analysis/   §18.1 artifact classes (later phases)
tools/run_nb.py           in-process notebook executor (see note below)
```

### Note on notebook execution

The Claude Science sandbox blocks socket binding (TCP and IPC), so a standard ZMQ Jupyter
kernel cannot start here. `tools/run_nb.py` executes a notebook's code cells in-process and
writes captured outputs back into the `.ipynb`, producing a genuine executed notebook without a
kernel socket. On the M5 Max, the notebooks also run normally under Jupyter/VS Code.

```
python tools/run_nb.py            # (invoked per-phase; see each phase dir)
```

## Status

| Gate | Phase | Session verification | Reviewer sign-off |
|------|-------|----------------------|-------------------|
| 0 | Protocol/estimand/analysis freeze | 31/31 checks PASS; evidence ready | **AWAITING** |

## Reproducibility (§18)

- SHA-256 manifests for every gate's evidence (`manifests/`).
- Large trajectories are **never** committed (`.gitignore`); they live under a checksummed
  archival manifest with ≥2 physical copies (§18.2).
- Released files are never overwritten; corrections create a new linked version (§18.2).
