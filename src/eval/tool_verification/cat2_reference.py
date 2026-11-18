"""
Category 2: comparison with independent reference implementations.
"""
from __future__ import annotations

import functools
import glob
import os

import gsw
import numpy as np
import pandas as pd
import xarray as xr
from scipy import signal, stats

from src.eval.tool_verification import refs
from src.eval.tool_verification.common import (ar1, case, coord_like, lat_lon_like, max_abs, monthly_axis,
                                               refill, rel_err, same_nan, zone_area)
from src.tools.basic_calculator import (analyze_correlation, compute_basic_stat, compute_difference_metrics,
                                        detect_outliers, fit_1D_trend, fit_multivariate_relation, fit_trend_map,
                                        temporal_aggregate)
from src.tools.data_loader import (load_areas, load_basin_mask, load_climatology_data, load_data_with_geolocator,
                                   load_depth_levels, load_landsea_mask, load_regular_data, load_volumes_for_cmip6)
from src.tools.data_openers import get_dataset_cfg, woa23_files
from src.tools.preprocessing import regrid_horizontal, regrid_vertical
from src.tools.scientific_calculator import (compute_anomaly_vs_own_climatology, compute_enso_nino_index,
                                             compute_pdo_index, compute_potential_density_anomaly)

CAT = "cat2_reference"
GSW_CV = os.path.join(os.path.dirname(gsw.__file__), "tests", "gsw_cv_v3_0.npz")
CMIP6_ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"


# ---------------------------------------------------------------------------
# compute_potential_density_anomaly: TEOS-10 check values and pointwise GSW
# ---------------------------------------------------------------------------

@case("compute_potential_density_anomaly", "TEOS-10 check casts: sigma0 and in-situ rho from SP and t")
def density_check_values(rec):
    cv = np.load(GSW_CV, allow_pickle=True)
    g = lambda key: np.asarray(cv[key], dtype="float64")
    sp, t, p, lat, lon = g("SP_chck_cast"), g("t_chck_cast"), g("p_chck_cast"), g("lat_chck_cast"), g("long_chck_cast")
    sigma0_cv, rho_cv, pt0 = g("sigma0"), g("rho"), g("pt0_from_t")
    tol_sigma, tol_rho = float(cv["sigma0_ca"]), float(cv["rho_ca"])
    worst = {"in-situ sigma0": 0.0, "potential sigma0": 0.0, "in-situ rho": 0.0}
    for k in range(sp.shape[1]):
        ok = np.isfinite(sp[:, k])
        depth = -gsw.z_from_p(p[ok, k], lat[k])
        coords = {"depth": depth, "lat": lat[k], "lon": lon[k]}
        temp = xr.DataArray(t[ok, k], dims="depth", coords=coords, attrs={"units": "degC"})
        sal = xr.DataArray(sp[ok, k], dims="depth", coords=coords, attrs={"units": "1"})
        s0 = compute_potential_density_anomaly(temp, sal, temperature_type="in_situ")
        worst["in-situ sigma0"] = max(worst["in-situ sigma0"], max_abs(s0.values, sigma0_cv[ok, k]))
        theta = xr.DataArray(pt0[ok, k], dims="depth", coords=coords, attrs={"units": "degC"})
        s0p = compute_potential_density_anomaly(theta, sal, temperature_type="potential")
        worst["potential sigma0"] = max(worst["potential sigma0"], max_abs(s0p.values, sigma0_cv[ok, k]))
        rho = compute_potential_density_anomaly(temp, sal, ref_pressure=None, temperature_type="in_situ")
        worst["in-situ rho"] = max(worst["in-situ rho"], max_abs(rho.values + 1000, rho_cv[ok, k]))
    for label, value in worst.items():
        tol = tol_rho if "rho" in label else tol_sigma
        rec.metric(f"{label} vs TEOS-10 check value", value, threshold=10 * tol, tier="T0",
                   method="gsw_cv_v3_0 check casts; depth = -z_from_p(p)",
                   note=f"check-value tolerance {tol:.1e}; 10x allows the depth->pressure round trip")


def _pointwise_sigma(temp, sal, kind, ref_pressure=0.0):
    """sigma by a scalar loop over every cell, laying out nothing by broadcasting."""
    lat, lon = lat_lon_like(temp)
    depth = coord_like(temp, "depth")
    t, s = temp.values, sal.values
    out = np.full(t.shape, np.nan)
    for index in zip(*np.where(np.isfinite(t) & np.isfinite(s))):
        p = gsw.p_from_z(-depth[index], lat[index])
        sa = gsw.SA_from_SP(s[index], p, lon[index], lat[index])
        ct = gsw.CT_from_pt(sa, t[index]) if kind == "potential" else gsw.CT_from_t(sa, t[index], p)
        out[index] = gsw.rho(sa, ct, p if ref_pressure is None else ref_pressure) - 1000.0
    return out


@case("compute_potential_density_anomaly", "broadcast fields vs a pointwise GSW loop (WOA23, three models)")
def density_pointwise(rec):
    cases = [("woa23 (in-situ)",
              load_climatology_data("woa23", "temperature", depth_range=[0, 3000], lat_range=[30, 33], lon_range=[-40, -37]),
              load_climatology_data("woa23", "salinity", depth_range=[0, 3000], lat_range=[30, 33], lon_range=[-40, -37]),
              "in_situ")]
    for model in ("CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"):
        kw = dict(model=model, time_range=["2000-01", "2000-02"], depth_range=[0, 2000], lat_range=[30, 33],
                  lon_range=[320, 323])
        cases.append((f"{model} (potential)", load_regular_data("cmip6_historical", "thetao", **kw),
                      load_regular_data("cmip6_historical", "so", **kw), "potential"))
    for label, temp, sal, kind in cases:
        temp, sal = temp.load(), sal.load()
        tool = compute_potential_density_anomaly(temp, sal)
        rec.metric(f"{label}: temperature type detected", tool.attrs["source_temperature"].startswith(
            "potential" if kind == "potential" else "in-situ"), op="true", tier="T0",
            note=tool.attrs["source_temperature"])
        ref = _pointwise_sigma(temp, sal, kind)
        rec.metric(f"{label}: sigma0 vs pointwise loop", max_abs(tool.values, ref), threshold=1e-10, tier="T0",
                   note=f"{int(np.isfinite(ref).sum())} cells, dims {temp.dims}")
        rec.metric(f"{label}: same cells missing", same_nan(tool.values, ref), op="true", tier="T0")
        tool = compute_potential_density_anomaly(temp, sal, ref_pressure=2000.0)
        rec.metric(f"{label}: sigma2 vs pointwise loop", max_abs(tool.values, _pointwise_sigma(temp, sal, kind, 2000.0)),
                   threshold=1e-10, tier="T0")


@case("compute_potential_density_anomaly", "cost of the wrong temperature route grows with depth")
def density_wrong_route(rec):
    kw = dict(model="GFDL-ESM4", time_range=["2000-01", "2000-01"], lat_range=[-30, -28], lon_range=[200, 202])
    temp = load_regular_data("cmip6_historical", "thetao", **kw).load()
    sal = load_regular_data("cmip6_historical", "so", **kw).load()
    right = compute_potential_density_anomaly(temp, sal)
    wrong = compute_potential_density_anomaly(temp, sal, temperature_type="in_situ")
    diff = np.abs(wrong - right).mean([d for d in right.dims if d != "depth"])
    at = lambda z: float(diff.sel(depth=z, method="nearest"))
    rec.metric("|sigma0 error| from treating thetao as in-situ, at 2000 / 4000 m",
               f"{at(2000):.4f} / {at(4000):.4f} kg m-3", op="info", tier="T3",
               note="docstring states about 0.01 at 2000 m and 0.026 at 4000 m")
    rec.metric("the error at 4000 m exceeds that at 2000 m", bool(at(4000) > at(2000)), op="true", tier="T1",
               note=f"500 m: {at(500):.4f}, 1000 m: {at(1000):.4f}")


# ---------------------------------------------------------------------------
# Correlation
# ---------------------------------------------------------------------------

@case("analyze_correlation", "time correlation vs scipy, effective sample size vs Bretherton et al. (1999)")
def correlation_scipy(rec):
    rng = np.random.default_rng(10)
    n, cells = 240, 12
    common = ar1(rng, n, 0.7, (cells,))
    x = common + ar1(rng, n, 0.5, (cells,))
    y = 0.6 * common + ar1(rng, n, 0.8, (cells,))
    time = monthly_axis("2000-01", n)
    v1 = xr.DataArray(x.T, dims=("time", "cell"), coords={"time": time})
    v2 = xr.DataArray(y.T, dims=("time", "cell"), coords={"time": time})
    for method, func in (("pearson", stats.pearsonr), ("spearman", stats.spearmanr)):
        result = analyze_correlation(v1, v2, method=method)
        r_ref = np.array([func(x[c], y[c])[0] for c in range(cells)])
        p_ref = np.array([func(x[c], y[c])[1] for c in range(cells)])
        rec.metric(f"{method}: coefficient vs scipy", max_abs(result["correlation"].values, r_ref), threshold=1e-12,
                   tier="T0")
        rec.metric(f"{method}: unadjusted p vs scipy", rel_err(result["pvalue_unadjusted"].values, p_ref),
                   threshold=1e-8, tier="T0")
        neff = np.array([refs.n_eff_pair(x[c], y[c]) for c in range(cells)])
        rec.metric(f"{method}: n_eff vs n(1 - r1x r1y)/(1 + r1x r1y)", rel_err(result["n_effective"].values, neff),
                   threshold=1e-12, tier="T0")
        t = r_ref * np.sqrt(np.maximum(neff - 2, 1) / (1 - r_ref ** 2))
        p_adj = 2 * stats.t.sf(np.abs(t), np.maximum(neff - 2, 1))
        rec.metric(f"{method}: adjusted p vs t-test on n_eff - 2 dof", rel_err(result["pvalue"].values, p_adj),
                   threshold=1e-8, tier="T0")
    result = analyze_correlation(v1, v2, detrend=True)
    r_ref = np.array([stats.pearsonr(signal.detrend(x[c]), signal.detrend(y[c]))[0] for c in range(cells)])
    rec.metric("detrend=True vs scipy.signal.detrend then pearsonr", max_abs(result["correlation"].values, r_ref),
               threshold=1e-12, tier="T0")


@case("analyze_correlation", "lead-lag sign convention")
def correlation_lag(rec):
    rng = np.random.default_rng(11)
    n, k = 300, 3
    time = monthly_axis("1990-01", n)
    driver = rng.standard_normal(n)
    follower = np.roll(driver, k)            # follower(t) = driver(t - k): driver leads by k
    v1 = xr.DataArray(follower, dims="time", coords={"time": time})
    v2 = xr.DataArray(driver, dims="time", coords={"time": time})
    r = {lag: float(analyze_correlation(v1, v2, lag=lag)["correlation"]) for lag in (-k, 0, k)}
    rec.metric(f"var2 leading var1 by {k} steps peaks at lag=+{k}", r[k], threshold=1 - 1e-12, op=">=", tier="T0",
               note=f"r(lag=-{k})={r[-k]:.3f}, r(0)={r[0]:.3f}")


@case("analyze_correlation", "area-weighted pattern correlation vs a weighted Pearson in numpy")
def correlation_area(rec):
    a = load_climatology_data("woa23", "temperature", depth_range=[0, 0]).isel(depth=0).load()
    b = load_climatology_data("woa23", "temperature", depth_range=[100, 100]).isel(depth=0).load()
    b = b.drop_vars("depth").assign_coords(depth=a["depth"])
    result = analyze_correlation(a, b, dim=["lat", "lon"], weighting="area")
    ref = refs.weighted_pearson(a.values, b.values, load_areas(dataset="woa23").values)
    rec.metric("WOA23 SST vs 100 m temperature", abs(float(result["correlation"]) - ref), threshold=1e-12, tier="T0",
               note=f"r = {ref:.6f}")
    unweighted = analyze_correlation(a, b, dim=["lat", "lon"])
    av, bv = a.values.astype("float64"), b.values.astype("float64")
    ok = np.isfinite(av) & np.isfinite(bv)
    rec.metric("weighting='none' vs scipy over the same cells",
               abs(float(unweighted["correlation"]) - stats.pearsonr(av[ok], bv[ok])[0]),
               threshold=1e-12, tier="T0")


# ---------------------------------------------------------------------------
# Trends
# ---------------------------------------------------------------------------

@case("fit_1D_trend", "linear, Theil-Sen, Mann-Kendall and AICc vs scipy and first principles")
def trend_reference(rec):
    rng = np.random.default_rng(12)
    n = 360
    time = monthly_axis("1985-01", n)
    y = 0.015 * np.arange(n) / 12 + 0.3 * ar1(rng, n, 0.6)
    series = xr.DataArray(y, dims="time", coords={"time": time}, attrs={"units": "degC"})
    result = fit_1D_trend(series, method="linear")
    x = (time - time[0]).days.values / 365.2425
    lr = stats.linregress(x, y)
    t = result["trend"]
    # The tool fits even the straight line with scipy's curve_fit, an iterative
    # least-squares solver, so it agrees with the closed form only to the
    # solver's convergence tolerance, not to machine precision.
    fit_note = "linear model fitted with curve_fit (iterative), not closed-form OLS"
    rec.metric("slope vs scipy.linregress", rel_err(t["slope"], lr.slope), threshold=1e-6, tier="T1", note=fit_note)
    rec.metric("intercept vs scipy.linregress", rel_err(result["params"]["b"], lr.intercept), threshold=1e-4, tier="T1",
               note=fit_note)
    rec.metric("unadjusted standard error vs linregress", rel_err(t["slope_se_unadjusted"], lr.stderr), threshold=1e-6,
               tier="T1", note=fit_note)
    rec.metric("unadjusted p vs linregress", rel_err(t["p_value_unadjusted"], lr.pvalue), threshold=1e-3, tier="T1",
               note=fit_note)
    residuals = y - (lr.slope * x + lr.intercept)
    neff, r1 = refs.n_eff_single(residuals)
    rec.metric("n_eff vs n(1 - r1)/(1 + r1) of the residuals", rel_err(t["n_effective"], neff), threshold=1e-6,
               tier="T1")
    se = np.sqrt(np.sum(residuals ** 2) / (neff - 2) / np.sum((x - x.mean()) ** 2))
    half = stats.t.ppf(0.975, neff - 2) * se
    rec.metric("95% interval = slope +- t(0.975, n_eff-2) * se_eff",
               rel_err([t["ci95_low"], t["ci95_high"]], [lr.slope - half, lr.slope + half]), threshold=1e-6, tier="T1")
    rss = float(np.sum(residuals ** 2))
    rec.metric("AICc of the linear fit", abs(result["metrics"]["aicc"] - refs.aicc(rss, n, 2)), threshold=1e-6,
               tier="T1")

    ts = fit_1D_trend(series, method="theil_sen")
    rec.metric("Theil-Sen slope vs brute-force median of pairwise slopes",
               rel_err(ts["trend"]["slope"], refs.theil_sen_slope(x, y)), threshold=1e-12, tier="T0")
    lo, hi = stats.theilslopes(y, x, alpha=0.95)[2:]
    rec.metric("Theil-Sen 95% interval vs scipy.theilslopes", max_abs([ts["trend"]["ci95_low"], ts["trend"]["ci95_high"]],
                                                                    [lo, hi]), threshold=1e-12, tier="T0")
    s, tau, z, p = refs.mann_kendall(y, continuity=False)
    rec.metric("Mann-Kendall tau vs S / (n(n-1)/2)", abs(ts["mann_kendall_tau"] - tau), threshold=1e-12, tier="T0")
    rec.metric("Mann-Kendall p vs normal approximation without continuity correction",
               rel_err(ts["trend"]["p_value"], p), threshold=1e-6, tier="T0",
               note=f"with the textbook continuity correction p would be {refs.mann_kendall(y)[3]:.4g}")

    quad = fit_1D_trend(series, method="quadratic")
    coef = np.polyfit(x, y, 2)
    rec.metric("quadratic coefficients vs numpy.polyfit", rel_err([quad["params"][k] for k in "abc"], coef),
               threshold=1e-6, tier="T0")


@case("fit_trend_map", "every cell equals fit_1D_trend and numpy.polyfit")
def trend_map_reference(rec):
    rng = np.random.default_rng(13)
    n, ny, nx = 300, 4, 5
    time = monthly_axis("1990-01", n)
    slopes = rng.uniform(-0.05, 0.05, (ny, nx))
    noise = ar1(rng, n, 0.5, (ny, nx)) * 0.4
    values = slopes[..., None] * (np.arange(n) / 12.0) + noise
    field = xr.DataArray(np.moveaxis(values, -1, 0), dims=("time", "lat", "lon"),
                         coords={"time": time, "lat": np.arange(ny), "lon": np.arange(nx)}, attrs={"units": "degC"})
    field[::17, 0, 0] = np.nan
    result = fit_trend_map(field, per_decade=False)
    keys = {"slope": "slope", "slope_se": "slope_se", "ci95_low": "ci95_low", "ci95_high": "ci95_high",
            "p_value": "p_value", "p_value_unadjusted": "p_value_unadjusted", "n_effective": "n_effective"}
    worst = {k: 0.0 for k in keys}
    poly = np.zeros((ny, nx))
    x = (time - time[0]).days.values / 365.2425
    for j in range(ny):
        for i in range(nx):
            ok = np.isfinite(field[:, j, i].values)
            poly[j, i] = np.polyfit(x[ok], field[:, j, i].values[ok], 1)[0]
            if (j, i) == (0, 0):
                continue
            one = fit_1D_trend(field[:, j, i], method="linear")["trend"]
            for k, k1 in keys.items():
                worst[k] = max(worst[k], rel_err(float(result[k][j, i]), one[k1]))
    for k, v in worst.items():
        rec.metric(f"gap-free cells, {k}: max relative difference from fit_1D_trend", v, threshold=1e-5, tier="T1",
                   note="fit_1D_trend's curve_fit limits agreement to ~1e-7 in the slope")
    rec.metric("slope vs numpy.polyfit (closed form)", rel_err(result["slope"].values, poly), threshold=1e-10, tier="T0")
    one = fit_1D_trend(field[:, 0, 0], method="linear")["trend"]
    rec.metric("series with gaps: n_eff equals fit_1D_trend's, as the docstring states",
               rel_err(float(result["n_effective"][0, 0]), one["n_effective"]), threshold=1e-5, tier="T1",
               note=f"lag-1 {float(result['lag1_autocorr'][0, 0]):.4f} (map: pairs across a gap dropped) vs "
                    f"{one['lag1_autocorr']:.4f} (1-D: finite points concatenated); p {float(result['p_value'][0, 0]):.3g} "
                    f"vs {one['p_value']:.3g}")


# ---------------------------------------------------------------------------
# Regression
# ---------------------------------------------------------------------------

@case("fit_multivariate_relation", "OLS inference vs QR least squares, VIF vs inverse correlation, ridge closed form")
def multivariate_reference(rec):
    rng = np.random.default_rng(14)
    ny, nx = 40, 50
    lat = np.linspace(-60, 60, ny)
    lon = np.linspace(-170, 170, nx)
    po4 = rng.uniform(0.1, 2.5, (ny, nx))
    temp = rng.uniform(0, 28, (ny, nx))
    sil = 5 * po4 + rng.normal(0, 0.8, (ny, nx))           # collinear with phosphate by construction
    no3 = 16 * po4 + 0.5 - 0.05 * temp + rng.normal(0, 0.6, (ny, nx))
    coords = {"lat": lat, "lon": lon}
    da = lambda v, name: xr.DataArray(v, dims=("lat", "lon"), coords=coords, name=name)
    target, preds = da(no3, "no3"), [da(po4, "po4"), da(temp, "temp"), da(sil, "sil")]
    result = fit_multivariate_relation(target, preds, method="linear")
    X = np.column_stack([po4.ravel(), temp.ravel(), sil.ravel()])
    ref = refs.ols_qr(X, no3.ravel())
    names = ["intercept", "po4", "temp", "sil"]
    rec.metric("coefficients vs QR least squares", rel_err([result["coefficients"][k] for k in names], ref["beta"]),
               threshold=1e-8, tier="T0", note=f"N:P slope {result['coefficients']['po4']:.3f} (planted 16)")
    d = result["diagnostics"]
    rec.metric("coefficient standard errors", rel_err([d["coefficient_se"][k] for k in names], ref["se"]),
               threshold=1e-8, tier="T0")
    rec.metric("F statistic", rel_err(d["f_statistic"], ref["f"]), threshold=1e-8, tier="T0")
    rec.metric("R2 of the full-sample fit", abs(result["metrics"]["r2"] - ref["r2"]), threshold=1e-10, tier="T0")
    vif = refs.vif_from_correlation(X)
    rec.metric("VIF vs diag(inverse correlation matrix)", rel_err([d["vif"][k] for k in names[1:]], vif),
               threshold=1e-8, tier="T0", note=f"VIF {np.round(vif, 2).tolist()}")
    rec.metric("collinearity warning raised when VIF >= 5", ("collinearity_warning" in d) == bool(vif.max() >= 5),
               op="true", tier="T0")
    single = fit_multivariate_relation(target, [preds[0]], method="linear")
    lr = stats.linregress(po4.ravel(), no3.ravel())
    rec.metric("single predictor: slope and intercept vs scipy.linregress",
               rel_err([single["coefficients"]["po4"], single["coefficients"]["intercept"]], [lr.slope, lr.intercept]),
               threshold=1e-10, tier="T0")

    alpha = 3.0
    ridge = fit_multivariate_relation(target, preds, method="ridge", alpha=alpha)
    Z = (X - X.mean(0)) / X.std(0)
    beta = np.linalg.solve(Z.T @ Z + alpha * np.eye(3), Z.T @ (no3.ravel() - no3.mean()))
    rec.metric("ridge (standardised) vs (Z'Z + alpha I)^-1 Z'y",
               rel_err([ridge["coefficients"][k] for k in names[1:]], beta), threshold=1e-8, tier="T0")
    rec.metric("ridge intercept = mean of y", rel_err(ridge["coefficients"]["intercept"], no3.mean()), threshold=1e-10,
               tier="T0")


# ---------------------------------------------------------------------------
# Robust statistics and difference metrics
# ---------------------------------------------------------------------------

@case("detect_outliers", "masks vs scipy and numpy, per cell along time")
def outliers_reference(rec):
    rng = np.random.default_rng(15)
    x = rng.standard_t(4, (600, 6))
    data = xr.DataArray(x, dims=("time", "cell"))
    ref = {
        "zscore": np.abs(stats.zscore(x, axis=0)) > 3.0,
        "iqr": (x < np.percentile(x, 25, axis=0) - 1.5 * stats.iqr(x, axis=0))
               | (x > np.percentile(x, 75, axis=0) + 1.5 * stats.iqr(x, axis=0)),
        "mad": np.abs(x - np.median(x, axis=0)) / stats.median_abs_deviation(x, axis=0, scale="normal") > 3.0,
        "percentile": (x < np.percentile(x, 0.5, axis=0)) | (x > np.percentile(x, 99.5, axis=0)),
    }
    for method, expected in ref.items():
        mask = detect_outliers(data, method=method).values
        rec.metric(f"{method}: elements that differ from the reference mask", int((mask != expected).sum()),
                   threshold=0, op="==", tier="T0", note=f"{int(expected.sum())} flagged by the reference")


@case("compute_difference_metrics", "weighted metrics vs numpy on a WOA23 field")
def difference_reference(rec):
    ref = load_climatology_data("woa23", "temperature", depth_range=[0, 0]).isel(depth=0).load()
    lat, lon = lat_lon_like(ref)
    data = refill(ref, ref.values + 0.3 * np.sin(np.radians(3 * lon)) + 0.2 * np.cos(np.radians(lat)))
    m = compute_difference_metrics(data, ref)
    w = load_areas(dataset="woa23").values
    d = data.values - ref.values
    ok = np.isfinite(d)
    w, d, a, b = w[ok], d[ok], data.values[ok], ref.values[ok]
    wm = lambda v: np.sum(w * v) / np.sum(w)
    expected = {
        "bias": wm(d), "rmse": np.sqrt(wm(d ** 2)), "mae": wm(np.abs(d)),
        "centered_rmse": np.sqrt(wm((d - wm(d)) ** 2)),
        "std_ratio": np.sqrt(wm((a - wm(a)) ** 2)) / np.sqrt(wm((b - wm(b)) ** 2)),
        "max_abs_error": np.abs(d).max(),
    }
    for key, value in expected.items():
        rec.metric(f"{key} vs numpy with load_areas weights", rel_err(m[key], value), threshold=1e-10, tier="T0")
    rec.metric("rmse^2 = bias^2 + centred_rmse^2", abs(m["rmse"] ** 2 - m["bias"] ** 2 - m["centered_rmse"] ** 2),
               threshold=1e-12, tier="T0")


# ---------------------------------------------------------------------------
# Time aggregation and anomalies
# ---------------------------------------------------------------------------

@case("temporal_aggregate", "running means vs pandas, monthly climatology and WOA seasons vs numpy")
def temporal_reference(rec):
    rng = np.random.default_rng(16)
    n = 120
    time = monthly_axis("2000-01", n)
    x = rng.standard_normal(n).cumsum()
    series = xr.DataArray(x, dims="time", coords={"time": time})
    s = pd.Series(x, index=time)
    for window, center in ((3, True), (4, True), (5, False)):
        tool = temporal_aggregate(series, freq="rolling", window=window, center=center).values
        ref = s.rolling(window, center=center, min_periods=window).mean().values
        rec.metric(f"rolling window={window} center={center} vs pandas", max_abs(tool, ref), threshold=1e-12,
                   tier="T0", note="NaN pattern must match too" if same_nan(tool, ref) else "NaN pattern differs")
        rec.metric(f"rolling window={window} center={center}: same ends missing", same_nan(tool, ref), op="true",
                   tier="T0")
    days = time.days_in_month.values.astype(float)
    clim = temporal_aggregate(series, freq="climatology")
    ref = np.array([np.sum((x * days)[time.month == m]) / np.sum(days[time.month == m]) for m in range(1, 13)])
    rec.metric("climatology (day-weighted per calendar month) vs numpy", max_abs(clim.values, ref), threshold=1e-12,
               tier="T0")
    seasons = temporal_aggregate(series, freq="season", season_convention="woa")
    q = (time.month.values - 1) // 3
    ref = [np.sum((x * days)[(time.year == y) & (q == k)]) / np.sum(days[(time.year == y) & (q == k)])
           for y in range(2000, 2010) for k in range(4)]
    rec.metric("WOA calendar-quarter seasons vs numpy", max_abs(seasons.values, ref), threshold=1e-12, tier="T0",
               note=f"labels {sorted(set(seasons['season'].values.tolist()))}")


@case("compute_anomaly_vs_own_climatology", "month and season anomalies vs pandas groupby")
def own_anomaly_reference(rec):
    rng = np.random.default_rng(17)
    n = 240
    time = monthly_axis("1995-01", n)
    x = 10 + 4 * np.cos(2 * np.pi * time.month.values / 12) + rng.standard_normal(n)
    series = xr.DataArray(x, dims="time", coords={"time": time})
    s = pd.Series(x, index=time)
    tool = compute_anomaly_vs_own_climatology(series, base_period=["2000-01", "2009-12"]).values
    base = s["2000-01":"2009-12"]
    ref = x - base.groupby(base.index.month).mean().reindex(time.month).values
    rec.metric("monthly anomaly on a 2000-2009 base vs pandas", max_abs(tool, ref), threshold=1e-12, tier="T0")
    season_of = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM", 6: "JJA", 7: "JJA", 8: "JJA",
                 9: "SON", 10: "SON", 11: "SON"}
    labels = pd.Index([season_of[m] for m in time.month])
    tool = compute_anomaly_vs_own_climatology(series, freq="season").values
    ref = x - s.groupby(labels).mean().reindex(labels).values
    rec.metric("meteorological-season anomaly vs pandas", max_abs(tool, ref), threshold=1e-12, tier="T0")


# ---------------------------------------------------------------------------
# Indices
# ---------------------------------------------------------------------------

@functools.lru_cache(None)
def ersst_record(start="1950-01", end="2019-12"):
    return load_regular_data("ersst", "sst", time_range=[start, end]).load()


@case("compute_enso_nino_index", "four Nino indices vs an independent numpy computation on ERSSTv5")
def nino_reference(rec):
    sst = ersst_record("1981-01", "2020-12")
    lat, lon = sst["lat"].values, sst["lon"].values
    weights_1d = zone_area(lat - 1, lat + 1, 2.0)
    month = sst["time"].dt.month.values
    boxes = {"nino1+2": (-10, 0, -90, -80), "nino3": (-5, 5, -150, -90), "nino3.4": (-5, 5, -170, -120),
             "nino4": (-5, 5, 160, -150)}
    for region, (la0, la1, lo0, lo1) in boxes.items():
        rows = (lat >= la0) & (lat <= la1)
        cols = ((lon >= lo0) | (lon <= lo1)) if lo1 < lo0 else ((lon >= lo0) & (lon <= lo1))
        box = sst.values[:, rows][:, :, cols]
        clim = np.stack([np.nanmean(box[month == m], axis=0) for m in range(1, 13)])
        anomaly = box - clim[month - 1]
        w = np.broadcast_to(weights_1d[rows][:, None], anomaly.shape[1:])
        valid = np.isfinite(anomaly)
        ref = np.nansum(anomaly * w, axis=(1, 2)) / np.sum(np.where(valid, w, 0), axis=(1, 2))
        tool = compute_enso_nino_index(sst, nino_region=region, climatology="self").values
        rec.metric(f"{region}: max |difference| over 1981-2020", max_abs(tool, ref), threshold=1e-5, tier="T0",
                   method="closed-form zone areas, inclusive box, own monthly climatology",
                   note=f"{int(rows.sum())} rows x {int(cols.sum())} columns")


@case("compute_pdo_index", "EOF vs an independent SVD of the same preprocessing on ERSSTv5")
def pdo_reference(rec):
    sst = ersst_record("1950-01", "2019-12")
    index, eof, explained, regression = compute_pdo_index(sst)
    lat, lon = sst["lat"].values, sst["lon"].values
    values = sst.values
    month = sst["time"].dt.month.values
    clim = np.stack([np.nanmean(values[month == m], axis=0) for m in range(1, 13)])
    anomaly = values - clim[month - 1]
    area = np.broadcast_to(zone_area(lat - 1, lat + 1, 2.0)[:, None], anomaly.shape[1:])
    valid = np.isfinite(anomaly)
    global_mean = np.nansum(anomaly * area, axis=(1, 2)) / np.sum(np.where(valid, area, 0), axis=(1, 2))
    rows = (lat >= 20) & (lat <= 70)
    lon360 = lon % 360
    cols = (lon360 >= 110) & (lon360 <= 260)
    box = anomaly[:, rows][:, :, cols] - global_mean[:, None, None]
    w = np.sqrt(area[rows][:, cols])
    matrix = (box * w).reshape(box.shape[0], -1)
    keep = np.isfinite(matrix).all(axis=0)
    pc, pattern, frac = refs.eof_svd(matrix[:, keep])
    pc = (pc - pc.mean()) / pc.std()
    r = np.corrcoef(index.values, pc)[0, 1]
    rec.metric("|corr(index, SVD PC)|", abs(r), threshold=0.99999, op=">=", tier="T0",
               method="numpy SVD; closed-form areas; own monthly anomalies and global-mean removal")
    rec.metric("explained variance: tool minus SVD (percentage points)", abs(explained - 100 * frac), threshold=1e-3,
               tier="T0", note=f"tool {explained:.4f} %")
    tool_eof = eof.values.reshape(-1)
    # the tool divides its eigenvector by sqrt(area) to return it in field units
    ref_eof = np.full(rows.sum() * cols.sum(), np.nan)
    ref_eof[keep] = pattern / w.reshape(-1)[keep]
    ok = np.isfinite(tool_eof) & np.isfinite(ref_eof)
    rec.metric("|corr(EOF pattern, SVD eigenvector)|", abs(np.corrcoef(tool_eof[ok], ref_eof[ok])[0, 1]),
               threshold=0.99999, op=">=", tier="T0")
    lonbox = lon[cols]
    western = eof.values[:, (lonbox >= 150) | (lonbox <= -170)]
    rec.metric("sign convention: mean EOF west of 170W in the box is negative (cool)", float(np.nanmean(western)),
               threshold=0.0, op="<=", tier="T0")
    reg_ref = np.nanmean((box - box.mean(0)) * ((index.values - index.values.mean()) / index.values.std())[:, None, None],
                         axis=0)
    ok = np.isfinite(regression.values) & np.isfinite(reg_ref) & keep.reshape(reg_ref.shape)
    rec.metric("regression pattern vs numpy regression on the standardised index",
               max_abs(regression.values[ok], reg_ref[ok]), threshold=1e-4, tier="T0", note="degC per standard deviation")


@case("compute_basic_stat", "area- and volume-weighted means on model grids vs numpy with the tools' own weights")
def basic_stat_models(rec):
    for model in ("CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"):
        tos = load_regular_data("cmip6_historical", "tos", model=model, time_range=["2000-01", "2000-03"]).load()
        tool = compute_basic_stat(tos, "area_weighted_mean", dim=["lat", "lon"]).values
        areas = load_areas(dataset="cmip6_historical", model=model).values
        ref = []
        areas = areas.astype("float64")
        for k in range(tos.sizes["time"]):
            v = tos.values[k].astype("float64")
            ok = np.isfinite(v) & np.isfinite(areas)
            ref.append(np.sum(v[ok] * areas[ok]) / np.sum(areas[ok]))
        rec.metric(f"{model} tos: global mean vs float64 numpy with load_areas", rel_err(tool, ref), threshold=1e-6,
                   tier="T0", note="float32 inputs: the tool accumulates in float32")
        thetao = load_regular_data("cmip6_historical", "thetao", model=model, time_range=["2000-01", "2000-01"],
                                   depth_range=[0, 500], lat_range=[-30, 30], lon_range=[150, 250]).load()
        tool = float(compute_basic_stat(thetao, "volume_weighted_mean"))
        vol = load_volumes_for_cmip6("cmip6_historical", model, depth_range=[0, 500], lat_range=[-30, 30],
                                     lon_range=[150, 250],
                                     **({"time_range": ["2000-01", "2000-01"]} if model == "GFDL-ESM4" else {}))
        vol = vol.isel(time=0) if "time" in vol.dims else vol
        v = thetao.isel(time=0).values.astype("float64")
        vv = vol.transpose(*[d for d in thetao.isel(time=0).dims]).values.astype("float64")
        ok = np.isfinite(v) & np.isfinite(vv)
        rec.metric(f"{model} thetao 0-500 m, 30S-30N: volume-weighted mean vs numpy with load_volumes_for_cmip6",
                   rel_err(tool, np.sum(v[ok] * vv[ok]) / np.sum(vv[ok])), threshold=1e-6, tier="T0")


# ---------------------------------------------------------------------------
# Regridding
# ---------------------------------------------------------------------------

@case("regrid_horizontal", "1deg -> 2deg against the analytic parallel-bounded conservative weights")
def regrid_to_coarse(rec):
    template = load_climatology_data("woa23", "temperature", depth_range=[0, 0]).isel(depth=0).load()
    lat, lon = lat_lon_like(template)
    f = np.sin(np.radians(lat)) * np.cos(np.radians(2 * lon)) + 0.5 * np.cos(np.radians(lat)) ** 2
    out = regrid_horizontal(refill(template, f), target="ersst")
    expected, _ = refs.conservative_regular(f, template["lat"].values, template["lon"].values,
                                            out["lat"].values, out["lon"].values)
    diff = np.abs(out.values - expected)
    rec.metric("smooth field: max difference from the analytic weights", float(np.nanmax(diff)), threshold=1e-3,
               tier="T1", method="refs.conservative_regular (overlap in sin(lat) x overlap in lon)",
               note=f"median {np.nanmedian(diff):.1e}; the remainder is great-circle vs parallel edges")
    out = regrid_horizontal(template, target="ersst", na_thres=1.0)
    expected, _ = refs.conservative_regular(template.values, template["lat"].values, template["lon"].values,
                                            out["lat"].values, out["lon"].values)
    both = np.isfinite(out.values) & np.isfinite(expected)
    rec.metric("real WOA23 SST with land (na_thres=1): max difference",
               float(np.nanmax(np.abs(out.values - expected)[both])), threshold=0.05, tier="T1",
               note=f"{int(both.sum())} cells compared")
    rec.metric("real WOA23 SST with land: coastal cells filled only by the tool",
               int((np.isfinite(out.values) & ~np.isfinite(expected)).sum()), op="info", tier="T3",
               note="a great-circle edge reaches slivers of a neighbouring wet cell that parallel edges do not")


@case("regrid_vertical", "WOA23 profiles onto GLODAP levels vs numpy.interp column by column")
def regrid_vertical_reference(rec):
    temp = load_climatology_data("woa23", "temperature", lat_range=[-20, -10], lon_range=[60, 70]).load()
    out = regrid_vertical(temp, target="glodap", fill_surface=False)
    levels = np.asarray(get_dataset_cfg("glodap")["depth"]["levels"], dtype=float)
    z = temp["depth"].values
    v = temp.values.reshape(z.size, -1)
    ref = np.full((levels.size, v.shape[1]), np.nan)
    for c in range(v.shape[1]):
        ok = np.isfinite(v[:, c])
        if ok.sum() >= 1:
            zc, vc = z[ok], v[ok, c]
            inside = (levels >= zc.min()) & (levels <= zc.max())
            # a gap in the column leaves the straddling target missing, as numpy.interp would not
            ref[inside, c] = np.interp(levels[inside], zc, vc)
    tool = out.values.reshape(levels.size, -1)
    both = np.isfinite(tool) & np.isfinite(ref)
    rec.metric("values where both are defined", max_abs(tool[both], ref[both]), threshold=1e-5, tier="T0",
               note="WOA23 is float32")
    rec.metric("cells defined by numpy.interp but missing in the tool", int((np.isfinite(ref) & ~np.isfinite(tool)).sum()),
               op="info", tier="T3", note="nonzero only where a column has interior gaps")


# ---------------------------------------------------------------------------
# Loaders against raw reads
# ---------------------------------------------------------------------------

def _match_regular(tool2d, tool_lat, tool_lon, raw2d, raw_lat, raw_lon):
    """Put a raw 2-D field onto the tool's (lat, lon) labels by coordinate value."""
    ii = np.array([int(np.argmin(np.abs(raw_lat - la))) for la in tool_lat])
    jj = np.array([int(np.argmin(np.abs(((raw_lon - lo + 180) % 360) - 180))) for lo in tool_lon])
    lat_ok = np.allclose(raw_lat[ii], tool_lat, atol=1e-4)
    lon_ok = np.allclose(((raw_lon[jj] - tool_lon + 180) % 360) - 180, 0, atol=1e-4)
    return raw2d[np.ix_(ii, jj)], lat_ok and lon_ok


@case("load_regular_data", "ERSSTv5, GFDL-ESM4 and CESM2 against raw netCDF4 reads")
def loader_timeseries_raw(rec):
    import netCDF4
    path = os.path.join(get_dataset_cfg("ersst")["root"], get_dataset_cfg("ersst")["file"])
    tool = load_regular_data("ersst", "sst", time_range=["1998-01", "1998-01"]).isel(time=0).load()
    with netCDF4.Dataset(path) as ds:
        times = netCDF4.num2date(ds["time"][:], ds["time"].units)
        k = [i for i, t in enumerate(times) if t.year == 1998 and t.month == 1][0]
        raw = np.ma.filled(ds["sst"][k].astype("float64"), np.nan)
        raw_lat, raw_lon = ds["lat"][:].astype(float), ds["lon"][:].astype(float)
    matched, coords_ok = _match_regular(tool.values, tool["lat"].values, tool["lon"].values, raw, raw_lat, raw_lon)
    rec.metric("ERSSTv5 Jan 1998: every cell equals the raw file", max_abs(tool.values, matched), threshold=0.0,
               op="<=", tier="T0", note="after the loader's latitude flip and 0-360 -> -180..180 roll")
    rec.metric("ERSSTv5: coordinates map one to one", coords_ok and same_nan(tool.values, matched), op="true", tier="T0")

    for model, variable in (("GFDL-ESM4", "tos"), ("CESM2", "thetao")):
        files = sorted(glob.glob(os.path.join(CMIP6_ROOT, model, variable, "*.nc")))
        tool = load_regular_data("cmip6_historical", variable, model=model, time_range=["2000-06", "2000-06"]).isel(time=0)
        if "depth" in tool.dims:
            tool = tool.isel(depth=[0, 10, 30])
        tool = tool.load()
        raw, found = None, None
        for f in files:
            with netCDF4.Dataset(f) as ds:
                times = ds["time"]
                dates = netCDF4.num2date(times[:], times.units, getattr(times, "calendar", "standard"))
                hits = [i for i, t in enumerate(dates) if t.year == 2000 and t.month == 6]
                if hits:
                    var = ds[variable]
                    raw = np.ma.filled(var[hits[0]].astype("float64"), np.nan)
                    if "lev" in ds.variables:
                        lev = ds["lev"][:].astype(float)
                        scale = 0.01 if "cm" in getattr(ds["lev"], "units", "").lower() or "centi" in getattr(ds["lev"], "units", "").lower() else 1.0
                        found = lev * scale
                    raw_lat, raw_lon = ds["lat"][:].astype(float), ds["lon"][:].astype(float)
                    break
        if raw is None:
            rec.metric(f"{model} {variable}: raw month found", False, op="true", tier="T0")
            continue
        if model == "GFDL-ESM4":
            matched, coords_ok = _match_regular(tool.values, tool["lat"].values, tool["lon"].values, raw, raw_lat, raw_lon)
            rec.metric(f"{model} {variable} Jun 2000 vs raw", max_abs(tool.values, matched), threshold=0.0, tier="T0")
            rec.metric(f"{model}: coordinates map one to one and land matches", coords_ok and same_nan(tool.values, matched),
                       op="true", tier="T0")
        else:
            raw = raw[[0, 10, 30]]
            rec.metric(f"{model} {variable} Jun 2000, levels 0/10/30, index by index vs raw", max_abs(tool.values, raw),
                       threshold=0.0, tier="T0", note="curvilinear: the loader keeps the index dimensions")
            rec.metric(f"{model}: depth = lev (cm) / 100", max_abs(tool["depth"].values, found[[0, 10, 30]]),
                       threshold=1e-4, tier="T0", note=f"{tool['depth'].values}")


@case("load_climatology_data", "WOA23, GLODAP and SEANOE against raw netCDF4 reads")
def loader_climatology_raw(rec):
    for variable, time_slice, raw_var in (("temperature", "annual", "t_an"), ("temperature", "monthly", "t_an"),
                                          ("nitrate", "annual", "n_an")):
        tool = load_climatology_data("woa23", variable, time_slice=time_slice, depth_range=[0, 500],
                                     **({"months": [7]} if time_slice == "monthly" else {})).load()
        files = woa23_files(variable, time_slice)["files"]
        path = files[6][1] if time_slice == "monthly" else files[0][1]
        raw, _ = refs.raw_nc(path, raw_var)
        raw = raw[0]
        depth, _ = refs.raw_nc(path, "depth")
        raw = raw[depth <= 500]
        values = tool.isel(month=0).values if time_slice == "monthly" else tool.values
        rec.metric(f"WOA23 {variable} {time_slice}{' (July)' if time_slice == 'monthly' else ''} 0-500 m vs raw",
                   max_abs(values, raw), threshold=0.0, tier="T0", note=os.path.basename(path))
        rec.metric(f"WOA23 {variable} {time_slice}: land pattern identical", same_nan(values, raw), op="true", tier="T0")
    cfg = get_dataset_cfg("glodap")
    path = os.path.join(cfg["root"], cfg["filename_template"].format(varfile=cfg["variables"]["TCO2"]["varfile"])) \
        if "TCO2" in cfg["variables"] else None
    if path:
        tool = load_climatology_data("glodap", "TCO2").load()
        raw, _ = refs.raw_nc(path, "TCO2")
        raw_lat, _ = refs.raw_nc(path, "lat")
        raw_lon, _ = refs.raw_nc(path, "lon")
        k = 5
        matched, coords_ok = _match_regular(tool.isel(depth=k).values, tool["lat"].values, tool["lon"].values,
                                            raw[k], raw_lat, raw_lon)
        rec.metric("GLODAP TCO2 level 5 vs raw (20.5..379.5 longitudes rolled)",
                   max_abs(tool.isel(depth=k).values, matched), threshold=0.0, tier="T0")
        rec.metric("GLODAP: coordinates one to one", coords_ok, op="true", tier="T0")
    cfg = get_dataset_cfg("seanoe_mld")
    version = cfg["default_version"]
    path = os.path.join(cfg["root"], cfg["versions"][version]["file"])
    tool = load_climatology_data("seanoe_mld", "mld_dr003").load()
    raw, attrs = refs.raw_nc(path, "mld_dr003")
    raw_lat, _ = refs.raw_nc(path, "lat")
    raw_lon, _ = refs.raw_nc(path, "lon")
    month_dim = raw.shape.index(12) if 12 in raw.shape else 0
    raw = np.moveaxis(raw, month_dim, 0)
    if raw.shape[1] == raw_lon.size:
        raw = np.swapaxes(raw, 1, 2)
    worst = 0.0
    for m in range(12):
        matched, coords_ok = _match_regular(tool.isel(month=m).values, tool["lat"].values, tool["lon"].values,
                                            raw[m], raw_lat, raw_lon)
        worst = max(worst, max_abs(tool.isel(month=m).values, matched))
    rec.metric(f"SEANOE mld_dr003 ({version}) all months vs raw", worst, threshold=0.0, tier="T0")


@case("load_depth_levels", "level centres and bounds vs the raw files")
def depth_levels_raw(rec):
    for variable, time_slice in (("temperature", "annual"), ("temperature", "monthly"), ("nitrate", "monthly"),
                                 ("oxygen", "annual")):
        levels, bounds = load_depth_levels(dataset="woa23", variable=variable, time_slice=time_slice)
        path = woa23_files(variable, time_slice)["files"][0][1]
        raw_b, _ = refs.raw_nc(path, "depth_bnds")
        raw_z, _ = refs.raw_nc(path, "depth")
        rec.metric(f"WOA23 {variable} {time_slice}: levels and bounds equal the file",
                   max(max_abs(levels.values, raw_z), max_abs(bounds.values, raw_b)), threshold=0.0, tier="T0",
                   note=f"{raw_z.size} levels, deepest layer {raw_b[-1].tolist()}")
    import netCDF4
    for model in ("CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"):
        levels, bounds = load_depth_levels(dataset="cmip6_historical", model=model)
        f = sorted(glob.glob(os.path.join(CMIP6_ROOT, model, "thetao", "*.nc")))[0]
        with netCDF4.Dataset(f) as ds:
            lev = ds["lev"]
            bname = getattr(lev, "bounds", "lev_bnds")
            scale_z = 0.01 if "c" in str(getattr(lev, "units", "m")).lower()[:1] else 1.0
            scale_b = 0.01 if "c" in str(getattr(ds[bname], "units", getattr(lev, "units", "m"))).lower()[:1] else 1.0
            z = lev[:].astype(float) * scale_z
            b = ds[bname][:].astype(float) * scale_b
            units = (getattr(lev, "units", "?"), getattr(ds[bname], "units", "(none)"))
        rec.metric(f"{model}: levels and bounds equal the file, each scaled by its own units",
                   max(max_abs(levels.values, z), max_abs(bounds.values, b)), threshold=1e-4, tier="T0",
                   note=f"lev units {units[0]!r}, {bname} units {units[1]!r}; top layer {bounds.values[0].tolist()}")


@case("load_regular_data", "cpc_indices: Nino, ONI and PDO series vs hand-parsed ASCII")
def cpc_raw(rec):
    root = get_dataset_cfg("cpc_indices")["root"]
    for product, file in (("ersst5_91-20", "ersst5.nino.mth.91-20.ascii"), ("ersst5_81-10", "ersst5.nino.mth.81-10.ascii"),
                          ("oisst", "sstoi.indices")):
        raw = refs.raw_cpc_monthly(os.path.join(root, file))
        for region in ("NINO1+2", "NINO3", "NINO4", "NINO3.4"):
            worst = 0.0
            for field, column in (("anomaly", f"{region}_ANOM"), ("total", region)):
                tool = load_regular_data("cpc_indices", region, product=product, field=field)
                ref = np.array([raw[(t.year, t.month)][column] for t in pd.DatetimeIndex(tool["time"].values)])
                worst = max(worst, max_abs(tool.values, ref))
                n = tool.size
            rec.metric(f"{product} {region}: anomaly and total equal the file", worst, threshold=0.0, tier="T0",
                       note=f"{n} months; file rows {len(raw)}")
    oni = refs.raw_oni(os.path.join(root, "oni.ascii.txt"))
    centre = {"DJF": 1, "JFM": 2, "FMA": 3, "MAM": 4, "AMJ": 5, "MJJ": 6, "JJA": 7, "JAS": 8, "ASO": 9, "SON": 10,
              "OND": 11, "NDJ": 12}
    tool = load_regular_data("cpc_indices", product="oni")
    ref = pd.Series([a for _, _, _, a in oni], index=pd.to_datetime([f"{y}-{centre[s]:02d}-01" for s, y, _, _ in oni]))
    rec.metric("ONI: values equal the file, each season at its centre month",
               max_abs(tool.values, ref.reindex(pd.DatetimeIndex(tool["time"].values)).values), threshold=0.0, tier="T0",
               note="DJF of year Y is placed at January of Y")
    pdo = refs.raw_pdo(os.path.join(root, "ersst.v6.pdo.dat"))
    tool = load_regular_data("cpc_indices", product="pdo")
    ref = np.array([pdo[(t.year, t.month)] for t in pd.DatetimeIndex(tool["time"].values)])
    rec.metric("PDO: values equal the file, 99.99 sentinels dropped", max_abs(tool.values, ref), threshold=0.0, tier="T0",
               note=f"{tool.size} of {len(pdo)} valid months")


@case("load_basin_mask", "basin codes vs the ASCII mask pivoted on its own coordinate columns")
def basin_mask_raw(rec):
    path = os.path.join(get_dataset_cfg("woa23")["masks_root"], "basinmask_01.msk")
    tool = load_basin_mask()
    columns = [f"Basin_{d}m" for d in tool["depth"].values]
    grids, _ = refs.raw_woa_mask(path, columns)
    worst, nan_ok = 0.0, True
    for k, column in enumerate(columns):
        worst = max(worst, max_abs(tool.isel(depth=k).values, grids[column]))
        nan_ok &= same_nan(tool.isel(depth=k).values, grids[column])
    rec.metric("every level identical to the pivoted table", worst, threshold=0.0, tier="T0",
               note=f"{len(columns)} levels")
    rec.metric("missing cells identical", nan_ok, op="true", tier="T0")


@case("load_landsea_mask", "bottom level index vs the ASCII mask pivoted on its own coordinate columns")
def landsea_mask_raw(rec):
    path = os.path.join(get_dataset_cfg("woa23")["masks_root"], "landsea_01.msk")
    import pandas as _pd
    column = [c for c in _pd.read_csv(path, comment="#", nrows=1).columns if c.lower().startswith("bottom_standard")][0]
    grids, _ = refs.raw_woa_mask(path, [column])
    raw = np.nan_to_num(grids[column], nan=0.0)
    rec.metric("identical to the pivoted table (cells absent from the table are 0)",
               max_abs(load_landsea_mask().values, raw), threshold=0.0, tier="T0")


@case("load_data_with_geolocator", "named-region masks vs an independent shapely intersection")
def geolocator_shapely(rec):
    import geopandas as gpd
    import shapely
    from pyproj import Geod
    marine = gpd.read_file(os.path.join(refs_path("natural_earth"), "10m_physical", "ne_10m_geography_marine_polys.shp"))
    lat_e, lon_e = np.arange(-90, 91, 1.0), np.arange(-180, 181, 1.0)
    west, south = np.meshgrid(lon_e[:-1], lat_e[:-1])
    east, north = np.meshgrid(lon_e[1:], lat_e[1:])
    boxes = shapely.box(west, south, east, north)          # (180, 360), the WOA23 cells
    shifted = shapely.box(west + 0.5, south, east + 0.5, north)   # the same cells moved half a cell east
    geod = Geod(ellps="WGS84")
    areas = load_areas(dataset="woa23").values
    for name in ("Mediterranean Sea", "Gulf of Mexico", "Caribbean Sea", "Bay of Bengal", "Red Sea", "Persian Gulf",
                 "Barents Sea", "Bering Sea", "Gulf of Guinea", "Arabian Sea"):
        rows = marine[marine["name"] == name]
        if rows.empty:
            rec.metric(f"{name}: present in Natural Earth", False, op="true", tier="T0")
            continue
        geom = shapely.union_all(rows.geometry.values)
        inter = shapely.area(shapely.intersection(boxes, geom))
        ref_mask = inter > 0
        field = load_data_with_geolocator("woa23", "temperature", loc_name=name, depth_range=[0, 0]).isel(depth=0)
        ocean = np.isfinite(load_climatology_data("woa23", "temperature", depth_range=[0, 0]).isel(depth=0).values)
        tool_mask = np.isfinite(field.values)
        ref_wet = ref_mask & ocean
        union = (tool_mask | ref_wet).sum()
        jaccard = (tool_mask & ref_wet).sum() / union if union else 0.0
        polygon_area = abs(geod.geometry_area_perimeter(geom)[0])
        rec.metric(f"{name}: Jaccard index of wet cells, tool vs shapely", float(jaccard), threshold=0.99, op=">=",
                   tier="T0", note=f"tool {int(tool_mask.sum())}, shapely {int(ref_wet.sum())} cells; "
                                   f"wet-cell area / polygon area = {areas[tool_mask].sum() / polygon_area:.2f}")
        moved = (shapely.area(shapely.intersection(shifted, geom)) > 0) & ocean
        both = (tool_mask | moved).sum()
        rec.metric(f"{name}: diagnostic - Jaccard against cells displaced +0.5 deg in longitude",
                   float((tool_mask & moved).sum() / both) if both else 0.0, op="info", tier="T3",
                   note="1.0 here means the mask was built on [c, c+1] instead of [c-0.5, c+0.5]")


def refs_path(kind):
    import yaml
    with open(os.path.join(os.path.dirname(__file__), "..", "..", "..", "configs", "default_paths.yaml")) as handle:
        cfg = yaml.safe_load(handle)
    return cfg["natural_earth_path"] if kind == "natural_earth" else cfg["cache_dir"]
