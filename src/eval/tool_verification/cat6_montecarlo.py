"""
Category 6: statistical properties by Monte Carlo.
"""
from __future__ import annotations

import numpy as np
import xarray as xr
from scipy import stats

from src.eval.tool_verification.common import ar1, case, monthly_axis
from src.tools.basic_calculator import analyze_correlation, detect_outliers, fit_1D_trend, fit_trend_map

CAT = "cat6_montecarlo"
REPS = 4000


def _field(values, n):
    """(cells, n) array -> (time, cell) DataArray on a monthly axis."""
    return xr.DataArray(values.T, dims=("time", "cell"),
                        coords={"time": monthly_axis("1960-01", n), "cell": np.arange(values.shape[0])})


@case("fit_trend_map", "false-positive rate and 95% interval coverage under AR(1) noise")
def trend_map_mc(rec):
    rng = np.random.default_rng(20)
    for n in (480, 120):
        for phi in (0.0, 0.3, 0.6, 0.9):
            judged = n == 480 and phi <= 0.6
            null = fit_trend_map(_field(ar1(rng, n, phi, (REPS,)), n), per_decade=False)
            fpr = float(null["significant"].mean())
            naive = float((null["p_value_unadjusted"] < 0.05).mean())
            rec.metric(f"n={n}, phi={phi}: false-positive rate (nominal 0.05)", fpr,
                       threshold=(0.03, 0.08) if judged else None, op="in" if judged else "info", tier="T1",
                       note=f"without the correction {naive:.3f}")
            slope = 0.02
            trend = slope * (np.arange(n) / 12.0)
            fitted = fit_trend_map(_field(trend + ar1(rng, n, phi, (REPS,)), n), per_decade=False)
            cover = float(((fitted["ci95_low"] <= slope) & (fitted["ci95_high"] >= slope)).mean())
            rec.metric(f"n={n}, phi={phi}: 95% interval coverage", cover,
                       threshold=(0.92, 0.97) if judged else None, op="in" if judged else "info", tier="T1")


@case("fit_1D_trend", "false-positive rate and coverage, linear and Theil-Sen")
def trend_1d_mc(rec):
    rng = np.random.default_rng(21)
    n, phi, reps = 240, 0.5, 600
    time = monthly_axis("1990-01", n)
    for method in ("linear", "theil_sen"):
        hits, covered = 0, 0
        for _ in range(reps):
            null = fit_1D_trend(xr.DataArray(ar1(rng, n, phi), dims="time", coords={"time": time}), method=method)
            hits += null["trend"]["p_value"] < 0.05
            trended = fit_1D_trend(xr.DataArray(0.02 * np.arange(n) / 12 + ar1(rng, n, phi), dims="time",
                                                coords={"time": time}), method=method)
            covered += trended["trend"]["ci95_low"] <= 0.02 <= trended["trend"]["ci95_high"]
        judged = method == "linear"
        rec.metric(f"{method}, n={n}, phi={phi}: false-positive rate", hits / reps,
                   threshold=(0.025, 0.09) if judged else None, op="in" if judged else "info", tier="T1",
                   note=f"{reps} replicates (s.e. ~0.009)" + ("" if judged else "; documented: Theil-Sen and "
                                                              "Mann-Kendall assume independence, so expected high"))
        rec.metric(f"{method}, n={n}, phi={phi}: 95% interval coverage", covered / reps,
                   threshold=(0.91, 0.975) if judged else None, op="in" if judged else "info", tier="T1")


@case("fit_1D_trend", "AICc model selection picks the generating model")
def trend_model_selection(rec):
    rng = np.random.default_rng(22)
    years = np.arange(1900, 2020)
    x = years - years[0]
    truths = {
        "quadratic": lambda: 2e-4 * x ** 2 + rng.normal(0, 0.15, x.size),
    }
    for truth, make in truths.items():
        chosen = [fit_1D_trend(xr.DataArray(make() + 5.0, dims="year", coords={"year": years}))["method_used"]
                  for _ in range(200)]
        share = float(np.mean([c == truth for c in chosen]))
        rec.metric(f"{truth} truth, positive series: share selected as {truth}", share, threshold=0.7, op=">=",
                   tier="T1", note=str({m: chosen.count(m) for m in set(chosen)}))
    # Control: centred on zero the series changes sign, exponential is not offered,
    # and the competition is between linear and the models with more parameters.
    linear = lambda: 0.01 * x + rng.normal(0, 0.15, x.size)
    chosen = [fit_1D_trend(xr.DataArray(linear() - 0.6, dims="year", coords={"year": years}))["method_used"]
              for _ in range(200)]
    rec.metric("linear truth crossing zero (no exponential candidate): share selected as linear",
               float(np.mean([c == "linear" for c in chosen])), op="info", tier="T1",
               note=str({m: chosen.count(m) for m in set(chosen)}))


@case("fit_1D_trend", "AICc model selection on a positive series generated with a linear trend")
def trend_model_selection_linear(rec):
    # Its own generator, so that the draws of the case above are not disturbed.
    rng = np.random.default_rng(22)
    years = np.arange(1900, 2020)
    x = years - years[0]
    chosen = [fit_1D_trend(xr.DataArray(0.01 * x + rng.normal(0, 0.15, x.size) + 5.0, dims="year",
                                        coords={"year": years}))["method_used"]
              for _ in range(200)]
    rec.metric("linear truth, positive series: share selected as linear",
               float(np.mean([c == "linear" for c in chosen])), threshold=0.7, op=">=", tier="T1",
               note=str({m: chosen.count(m) for m in set(chosen)}))


@case("analyze_correlation", "rejection rate for independent AR(1) series")
def correlation_mc(rec):
    rng = np.random.default_rng(23)
    for n, phi in ((240, 0.0), (240, 0.6), (240, 0.8), (120, 0.9)):
        cells = 2000
        v1 = _field(ar1(rng, n, phi, (cells,)), n)
        v2 = _field(ar1(rng, n, phi, (cells,)), n)
        result = analyze_correlation(v1, v2)
        rate = float((result["pvalue"] < 0.05).mean())
        naive = float((result["pvalue_unadjusted"] < 0.05).mean())
        judged = phi <= 0.6
        rec.metric(f"n={n}, phi={phi}: rejection rate with n_eff (nominal 0.05)", rate,
                   threshold=(0.03, 0.08) if judged else None, op="in" if judged else "info", tier="T1",
                   note=f"without the correction {naive:.3f}")


@case("detect_outliers", "false-alarm rates on Gaussian noise vs their theoretical values")
def outliers_mc(rec):
    rng = np.random.default_rng(24)
    x = xr.DataArray(rng.standard_normal(1_000_000), dims="time")
    q3 = stats.norm.ppf(0.75)
    theory = {
        "zscore": 2 * stats.norm.sf(3.0),
        "iqr": 2 * stats.norm.sf(q3 + 1.5 * 2 * q3),
        "mad": 2 * stats.norm.sf(3.0),
        "percentile": 0.01,
    }
    for method, expected in theory.items():
        rate = float(detect_outliers(x, method=method).mean())
        rec.metric(f"{method}: flagged share / theoretical {expected:.5f}", rate / expected, threshold=(0.93, 1.07),
                   op="in", tier="T1", note=f"flagged {rate:.5f}")
