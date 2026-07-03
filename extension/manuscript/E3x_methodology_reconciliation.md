# E3.x methodology reconciliation — executed work vs the Paper 1 extension plan

**Purpose.** `Paper1_Extension_Workflow.md` (June 2026 planning) specifies a sampling and statistical
design for Phases E2 and E3 that was subsequently superseded during execution. This note records where
the executed campaign departs from that plan so the divergences are documented for the Phase E5
Methods rewrite (Sections 2.5–2.6) rather than rediscovered. The plan's *scientific* structure — the
three E3 questions (replicate significance, indirect-restructuring, H-bond anomaly) — is retained and
is the organising spine of the E3.x notebooks. What changed is the *sampling and statistical
machinery* beneath those questions.

Framing throughout: these are observations about what was executed, not claims that the plan was
wrong. The design evolved for stated reasons; the Methods section must describe what ran.

---

## Summary of departures

| Element | Paper 1 plan | Executed (E2.x / E3.x) | Implication for Methods |
|---|---|---|---|
| Replicate design | 3 velocity replicates from **identical** equilibrated coordinates, differing seed only; extended start | **Diverse starts** drawn by Rg percentile (compact/mid/open/extended); 3–4 adaptive starts/system; 31 runs total | Rewrite 2.5: describe diverse-start construction and the adaptive 3-vs-4 criterion |
| Trajectory length | 100 ns target (checkpoint to 500 ns) | 200 ns per diverse start | Update 2.5 lengths |
| Primary inference | Inter-replicate mean ± SEM (n=3) | **Pooled ACF-corrected SE** with leave-one-out (`convergence_stats.py`) | Rewrite 2.6 primary-inference description |
| CI construction | Within-trajectory **block-bootstrap** (secondary tier) | **Parametric normal-approximation**: `N_eff = N/τ_int`, `se = std(ddof=1)/√N_eff`, `ci95 = 1.96·se`; FFT-ACF, Sokal windowing c=5.0, τ floored at 1 | Remove "block-bootstrap" from 2.6; describe the parametric estimator |
| Autocorrelation role | Secondary | **Primary** (decisive: raw τ_int 160–1201 ps, N_eff 4–28) | Foreground τ_int correction in 2.6 |
| Effect sizes | Recomputed from replicate means | From pooled corrected means; pairwise z and distribution overlap | Update effect-size description |
| Convergence adjudication | Not specified beyond SEM | **Leave-one-out** pooled-stability (7b); "robust" / "dominated" per observable | Add to 2.6 |
| Coordination cutoffs | Pair-specific RDF first-minimum, ±0.05 nm sensitivity | **Same** — retained as planned | No change (this part matches) |
| Cutoff sensitivity reporting | "robust to cutoff variation" (implied pass/fail) | Magnitude-sensitivity coefficient (graded) + qualitative shell-persistence; **not** pass/fail | Describe sensitivity as reported |
| Backbone H-bonds | Phase-3 defaults, resolution unspecified | **1 ps resolution required** (stride=10 aliases ~87% of episodes, KS D=0.843); backbone-only scope (donor N → acceptor O); per-start run-lengths | Add resolution requirement + scope to 2.6 |

---

## Substantive points for the Methods narrative

**1. The statistical change is upstream of the estimator.** The plan's SEM(n=3) is the natural summary
of its sampling design — three replicates from one equilibrated well are exchangeable draws from a
single ensemble, so their between-replicate scatter estimates a sampling error. The executed diverse-
start design deliberately breaks that exchangeability (starts are dispersed across the Rg
distribution). Once starts are not exchangeable draws from one basin, SEM(n=3) is no longer a coherent
summary: it would conflate between-basin spread with sampling error and report a falsely convergent
number for exactly the systems where convergence fails. The design was replaced first; the statistic
followed.

**2. Where the replacement is more robust: validity, not precision.** Promoting ACF correction to
primary is the correct ordering given the diagnostics — an uncorrected interval would have been
optimistic by 14–36×. More importantly, identical-coordinate velocity replicates cannot detect
conformational heterogeneity (three seeds from one well relax back into it and report tight, mutually
consistent means whether or not other basins exist). The diverse-start design plus leave-one-out is
what surfaced the multi-basin behaviour in CME_glyceline and GGE_water — results the planned n=3
design would likely have mislabelled as converged. This is a genuine gain in inferential validity.

**3. Where it is a lateral move: the CI construction itself.** The executed parametric normal-
approximation interval is efficient, deterministic, and free of block-length choices, but it is not
categorically more robust than the planned block-bootstrap. It assumes approximate normality of the
mean and a well-estimated τ_int, and its fixed 1.96 multiplier is a normal quantile, not a t quantile.
At the pooled N_eff actually reported (SASA ≈ 81, Rg ≈ 93) this is immaterial (t correction ~1–2%).
At per-start N_eff ≈ 27–31 a t-interval would be ~4–5% wider, and at the 10 ns diagnostic N_eff of
4–28 a fixed-z interval is materially anti-conservative. The plan's block-bootstrap would have been
most valuable in precisely that low-N_eff regime, and it was dropped rather than retained as a
cross-check.

**Recommendation for E5.** Restore the block-bootstrap not as the primary estimator but as a
validation layer: a moving-block bootstrap (block length ≈ 2τ_int) run as a cross-check on the
parametric CIs for the low-N_eff and multi-basin systems would directly test the normal-approximation
adequacy where it is weakest, letting the Methods state the parametric intervals are
bootstrap-corroborated rather than assumed. This is cheap relative to the trajectory generation
already sunk and recovers what the plan's secondary tier was reaching for, in a form consistent with
diverse-start sampling.

**4. Estimand caveat (must appear in Methods and Limitations).** For ergodic systems the diverse-start
pooled mean and a single-equilibrium mean coincide. For multi-basin systems (CME_glyceline, GGE_water)
the pooled mean is an equal-frame-weight average over non-interconverting basins — an
ensemble-of-starts quantity, a *different estimand* from the plan's implied single equilibrium mean.
These systems are reported as distributions, not as converged point values; true populations require
enhanced sampling (Paper 2).

**5. Small-n is not resolved by either framework.** Between-condition replication remains 3–4 starts.
The executed framework manages this more explicitly (the four-start default exists because leave-one-
out is near-degenerate at three) but does not eliminate it; effective n for cross-solvent
generalisation stays modest and should be stated as a limitation.

---

## Language to purge from the draft Methods (2.6)

Both descriptors below are inaccurate for the executed work and must not survive into the revision:

- "inter-replicate mean ± SEM (n=3)" — the executed inference is pooled ACF-corrected SE with
  leave-one-out over diverse starts.
- "block-bootstrap" (both the E2 objective and the E3.1 secondary tier) — the executed CI is a
  parametric ACF-corrected normal-approximation SE. (If the moving-block bootstrap validation layer
  above is added, "block-bootstrap" may return, but explicitly as corroboration, not as the primary CI.)
