"""
convergence_stats.py
====================
Canonical E2.x convergence statistics for the DES-peptide study (Paper 1).

This module is the single source of truth for the autocorrelation-corrected
convergence statistics used in the E2.x extended-sampling campaign. It is
extracted VERBATIM (no numerical changes) from the cells of
``E2_2_Convergence_Analysis.ipynb`` so that downstream analysis (E3.x) can
import the identical routine and obtain numerically identical results.

    from convergence_stats import prepare_series, compute_stat, pairwise, pooled_stability

--------------------------------------------------------------------------
IMPORTANT — what the CI actually is
--------------------------------------------------------------------------
There is NO bootstrap anywhere in this routine, despite any description of it
as a "block-bootstrap". The 95% CIs are *parametric, autocorrelation-corrected
normal-approximation standard errors*:

        N_eff = N / tau_int
        se    = std(ddof=1) / sqrt(N_eff)
        ci95  = 1.96 * se

A literal block-bootstrap re-implementation in E3.x would diverge from E2.x by
construction. Do not re-implement; import this module.

--------------------------------------------------------------------------
Conventions (must not be changed without re-baselining all E2.x results)
--------------------------------------------------------------------------
Series preparation (prepare_series):
  - Discard an equilibration window from the FRONT of the raw (1 ps) series.
    E2.x default: discard the first 20 ns.
  - Stride to 10 ps (STRIDE = 10 over SAVE_INTERVAL_PS = 1.0 ps frames).
    All statistics are computed on the strided 10 ps series; sub-10 ps
    correlation structure is intentionally not resolved.

Autocorrelation / integrated time (integrated_time):
  - ACF: FFT-based (autocorr_func_1d), normalised to 1.0 at lag 0.
  - tau_int: Sokal automatic windowing with c = 5.0.
  - tau_int is FLOORED at 1.0 (strided lag). Consequence: for well-mixed
    series N_eff is capped at N (the strided frame count). This floor is a
    real part of the convention — do not remove it.

Per-start statistics (compute_stat):
  - dt_ps = SAVE_INTERVAL_PS * STRIDE = 10.0  (tau_int_ns = tau_int * dt_ps/1000)
  - N_eff = N_strided / tau_int        (NOT N/(2 tau) or N/(1+2 tau))
  - se    = std(ddof=1) / sqrt(N_eff)
  - ci95  = 1.96 * se

Pairwise comparison (pairwise):
  - z       = |mean_A - mean_B| / sqrt(se_A**2 + se_B**2)
  - overlap = 50-bin histogram-intersection of the two distributions.

Pooled-ensemble stability (pooled_stability), the "cell 7b" routine:
  - Concatenate the (equilibrated, strided) series across all starts with
    equal frame weight.
  - pooled_mean = mean(concatenated)
  - pooled_ci   = 1.96 * std(concatenated, ddof=1) / sqrt(sum of per-start N_eff)
    (note: the pooled SD is taken on the raw concatenation, so between-start
     mean differences enter the variance and widen the CI for split systems.)
  - leave-one-out: drop each start in turn; flag shift_gt_ci when the
    re-pooled mean moves beyond pooled_ci. A start that does so DOMINATES.

PRECONDITION for callers: pass series already produced by prepare_series
(equilibrated window, strided to 10 ps). compute_stat / pairwise /
pooled_stability do NOT re-slice or re-stride.
"""

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- constants
SAVE_INTERVAL_PS = 1.0     # raw trajectory frame interval (E2.1 save cadence)
STRIDE           = 10      # analyse every 10th frame -> 10 ps spacing
DT_PS            = SAVE_INTERVAL_PS * STRIDE   # = 10.0 ; tau_int_ns conversion
DISCARD_NS       = 20.0    # equilibration window discarded from the front
DEFAULT_C        = 5.0     # Sokal windowing constant
TAU_FLOOR        = 1.0     # integrated-time floor (strided lags)


# ---------------------------------------------------------------- preparation
def prepare_series(raw, discard_ns=DISCARD_NS, stride=STRIDE,
                   save_interval_ps=SAVE_INTERVAL_PS):
    """Slice a raw 1 ps observable series to the equilibrated window and stride it.

    Parameters
    ----------
    raw : 1-D array-like
        The full per-frame observable series at the raw save cadence
        (SAVE_INTERVAL_PS, default 1 ps), i.e. one value per saved frame.
    discard_ns : float
        Equilibration window discarded from the FRONT, in ns (E2.x: 20).
    stride : int
        Keep every ``stride``-th frame after the discard (E2.x: 10 -> 10 ps).
    save_interval_ps : float
        Raw frame interval in ps (E2.x: 1.0).

    Returns
    -------
    np.ndarray
        The equilibrated, strided series that all downstream stats consume.
    """
    raw = np.asarray(raw, dtype=float)
    discard_frames = int(round(discard_ns * 1000.0 / save_interval_ps))
    return raw[discard_frames::stride]


# ---------------------------------------------------------------- autocorr
def autocorr_func_1d(x):
    """FFT autocorrelation function, normalised to 1.0 at lag 0. (cell 3)"""
    x = np.asarray(x, float)
    x = x - x.mean()
    n = len(x)
    nfft = 1
    while nfft < 2 * n:
        nfft *= 2
    f = np.fft.fft(x, n=nfft)
    acf = np.fft.ifft(f * np.conjugate(f))[:n].real
    acf /= acf[0]
    return acf


def integrated_time(x, c=DEFAULT_C):
    """Integrated autocorrelation time, Sokal windowing, floored at 1.0. (cell 3)

    Returns tau in units of (strided) lags. Multiply by DT_PS/1000 for ns.
    """
    taus = 2.0 * np.cumsum(autocorr_func_1d(x)) - 1.0
    m = np.arange(len(taus)) < c * taus
    w = int(np.argmin(m)) if np.any(m) else len(taus) - 1
    return max(taus[w], TAU_FLOOR)


# ---------------------------------------------------------------- per-start
def compute_stat(series, dt_ps=DT_PS):
    """Per-start statistics for one observable series. (cell 6 ``stat``)

    Expects ``series`` already equilibrated + strided (see prepare_series).
    Returns dict with mean, tau_int_ns, N_eff, se, ci95.
    """
    series = np.asarray(series, float)
    n = len(series)
    tau = integrated_time(series)
    n_eff = n / tau
    se = series.std(ddof=1) / np.sqrt(n_eff)
    return {"mean": float(series.mean()),
            "tau_int_ns": tau * dt_ps / 1000.0,
            "N_eff": n_eff,
            "se": float(se),
            "ci95": 1.96 * float(se)}


# ---------------------------------------------------------------- pairwise
def _overlap(a, b, bins=50):
    """50-bin histogram-intersection overlap coefficient. (cell 7)"""
    lo, hi = min(a.min(), b.min()), max(a.max(), b.max())
    edges = np.linspace(lo, hi, bins + 1)
    bw = edges[1] - edges[0]
    ha, _ = np.histogram(a, bins=edges, density=True)
    hb, _ = np.histogram(b, bins=edges, density=True)
    return float(np.sum(np.minimum(ha, hb)) * bw)


def pairwise(data, runs, observables=("rg", "sasa")):
    """Pairwise z and overlap across starts. (cell 7)

    ``data`` is the nested dict ``data[label][obs] -> series`` and
    ``data[label]['stat'][obs] -> compute_stat(...)``. ``runs`` is the ordered
    list of start labels. Returns a DataFrame.
    """
    labels = list(runs)
    crows = []
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            A, B = labels[i], labels[j]
            for obs in observables:
                sa, sb = data[A]["stat"][obs], data[B]["stat"][obs]
                z = abs(sa["mean"] - sb["mean"]) / np.sqrt(sa["se"]**2 + sb["se"]**2)
                ovl = _overlap(data[A][obs], data[B][obs])
                crows.append({"pair": f"{A} vs {B}", "observable": obs.upper(),
                              "delta_mean": round(sa["mean"] - sb["mean"], 4),
                              "z": round(z, 2), "overlap": round(ovl, 2)})
    return pd.DataFrame(crows)


# ---------------------------------------------------------------- pooled 7b
def pooled_stability(data, runs, observables=("rg", "sasa")):
    """Pooled-ensemble mean + leave-one-out stability. (cell 7b)

    Returns (pooled_stats_df, loo_df). Verdict strings can be derived as in
    E2.2: robust if max |leave-one-out shift| <= pooled ci95, else DOMINATED.
    """
    pool_rows, loo_rows = [], []
    for obs in observables:
        series_all = {lab: data[lab][obs] for lab in runs}
        pooled = np.concatenate(list(series_all.values()))
        neff_sum = sum(data[lab]["stat"][obs]["N_eff"] for lab in runs)
        pooled_mean = float(pooled.mean())
        pooled_ci = 1.96 * pooled.std(ddof=1) / np.sqrt(neff_sum)
        pool_rows.append({"observable": obs.upper(), "n_starts": len(runs),
                          "pooled_mean": round(pooled_mean, 4),
                          "N_eff_pooled": round(neff_sum, 1),
                          "ci95": round(float(pooled_ci), 4)})
        for drop in runs:
            rest = [series_all[lab] for lab in runs if lab != drop]
            m_rest = float(np.concatenate(rest).mean())
            loo_rows.append({"observable": obs.upper(), "dropped": drop,
                             "mean_without": round(m_rest, 4),
                             "shift": round(m_rest - pooled_mean, 4),
                             "shift_gt_ci": bool(abs(m_rest - pooled_mean) > pooled_ci)})
    return pd.DataFrame(pool_rows), pd.DataFrame(loo_rows)


def verdict(pooled_stats, loo, observable):
    """Convenience: reproduce the E2.2 robust/DOMINATED verdict line."""
    ci = float(pooled_stats.loc[pooled_stats.observable == observable, "ci95"].iloc[0])
    mx = float(loo.loc[loo.observable == observable, "shift"].abs().max())
    tag = "robust to any single start" if mx <= ci else "DOMINATED by one start"
    return ci, mx, tag


# ---------------------------------------------------------------- self-test
if __name__ == "__main__":
    # Minimal numerical self-test: a correlated series should give tau > 1,
    # and a white-noise series should hit the floor (tau == 1.0, N_eff == N).
    rng = np.random.default_rng(0)
    white = rng.standard_normal(5000)
    assert integrated_time(white) == 1.0, "white-noise floor not engaged"
    # AR(1) with phi=0.9 has theoretical tau_int ~ (1+phi)/(1-phi) = 19
    phi, x = 0.9, np.zeros(20000)
    for i in range(1, len(x)):
        x[i] = phi * x[i-1] + rng.standard_normal()
    tau = integrated_time(x)
    assert tau > 5.0, f"AR(1) tau unexpectedly small: {tau}"
    s = compute_stat(white)
    assert abs(s["N_eff"] - len(white)) < 1e-9, "white N_eff should equal N"
    print(f"self-test OK: white tau=1.0 (floored), AR(1) tau={tau:.1f}, "
          f"white N_eff={s['N_eff']:.0f}")
