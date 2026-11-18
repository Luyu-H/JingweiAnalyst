import numpy as np
import xarray as xr

models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
model_colors = {
    "CESM2": "#1f77b4",
    "GFDL-ESM4": "#ff7f0e",
    "MPI-ESM1-2-LR": "#2ca02c",
}

def to_psu(da):
    attrs = dict(da.attrs)
    u = str(attrs.get("units", "")).strip().lower()
    if u in ("psu", "1", "1.0", "1e-3", "0.001", ""):
        out = da.copy()
        out.attrs = attrs
        out.attrs["units"] = "psu"
        return out
    return transform_unit(data=da, from_unit=attrs.get("units", "1"), to_unit="psu")

def val(da):
    if hasattr(da, "compute"):
        da = da.compute()
    return float(da)

def ens_label(da):
    for key in ("variant_label", "member_id", "ensemble", "realization", "ripf"):
        if key in da.attrs:
            return str(da.attrs[key])
    return "unknown"

# WOA23 1971-2000 SSS
woa_sal = load_climatology_data(
    dataset="woa23", variable="salinity", time_slice="annual", period="1971-2000"
)
woa_levels, woa_bounds = load_depth_levels(
    dataset="woa23", variable="salinity", time_slice="annual", period="1971-2000"
)
woa_levels = np.asarray(woa_levels, dtype=float)
woa_lev_0_10 = [float(l) for l in woa_levels if 0.0 <= l <= 10.0]
woa_sal_lev = to_psu(woa_sal.sel(depth=woa_lev_0_10))
woa_sss = compute_basic_stat(data=woa_sal_lev, op="volume_weighted_mean", dim="depth")
woa_sss.attrs.update(woa_sal_lev.attrs)
woa_sss.attrs["units"] = "psu"

lsm = load_landsea_mask()
woa_sss = woa_sss.where(lsm > 1)

woa_global_mean = val(
    compute_basic_stat(data=woa_sss, op="area_weighted_mean", dim=["lat", "lon"])
)

sss_hist_monthly = {}
sss_ssp_monthly = {}
sss_hist_clim = {}
sss_hist_clim_regrid = {}
global_mean = {}
diff_to_woa = {}
rmse = {}
bias_scalar = {}
correlation = {}
bias_fields = {}
annual_hist_global = {}
annual_ssp_global = {}
hist_fit = {}
ssp_fit = {}
hist_slope_dec = {}
ssp_slope_dec = {}
hist_ci_dec = {}
ssp_ci_dec = {}
hist_r2 = {}
ssp_r2 = {}
hist_sig = {}
ssp_sig = {}
hist_stable = {}
slope_ratio = {}
ratio_meaningful = {}
trend_maps = {}
trend_slope_regrid = {}
trend_p_regrid = {}
delta_fields = {}
delta_global = {}
delta_basin = {}
ensemble_labels = {}
depth_levels_used = {}

for model in models:
    # Historical monthly SSS
    so_hist = load_regular_data(
        dataset="cmip6_historical", variable="so", model=model,
        time_range=["1950-01", "2014-12"]
    )
    so_hist = to_psu(so_hist)

    lev_h, _ = load_depth_levels(
        dataset="cmip6_historical", model=model, variable="so"
    )
    lev_h = np.asarray(lev_h, dtype=float)
    lev_h_0_10 = [float(l) for l in lev_h if 0.0 <= l <= 10.0]

    sss_hist = compute_basic_stat(
        data=so_hist.sel(depth=lev_h_0_10),
        op="volume_weighted_mean", dim="depth"
    )
    sss_hist.attrs.update(so_hist.attrs)
    sss_hist.attrs["units"] = "psu"
    sss_hist_monthly[model] = sss_hist

    # ssp245 monthly SSS
    so_ssp = load_regular_data(
        dataset="cmip6_ssp245", variable="so", model=model,
        time_range=["2071-01", "2100-12"]
    )
    so_ssp = to_psu(so_ssp)

    lev_s, _ = load_depth_levels(
        dataset="cmip6_ssp245", model=model, variable="so"
    )
    lev_s = np.asarray(lev_s, dtype=float)
    lev_s_0_10 = [float(l) for l in lev_s if 0.0 <= l <= 10.0]

    sss_ssp = compute_basic_stat(
        data=so_ssp.sel(depth=lev_s_0_10),
        op="volume_weighted_mean", dim="depth"
    )
    sss_ssp.attrs.update(so_ssp.attrs)
    sss_ssp.attrs["units"] = "psu"
    sss_ssp_monthly[model] = sss_ssp

    depth_levels_used[model] = {
        "historical": lev_h_0_10,
        "ssp245": lev_s_0_10,
    }

    # 1971-2000 model climatology and global mean
    sss_hist_7100 = sss_hist.sel(time=slice("1971-01", "2000-12"))
    sss_clim = temporal_aggregate(data=sss_hist_7100, freq="mean", op="mean")
    sss_clim.attrs.update(sss_hist.attrs)
    sss_clim.attrs["units"] = "psu"
    sss_hist_clim[model] = sss_clim

    global_mean[model] = val(
        compute_basic_stat(data=sss_clim, op="area_weighted_mean", dim=["lat", "lon"])
    )
    diff_to_woa[model] = global_mean[model] - woa_global_mean

    # Bias comparison on WOA23 grid
    sss_clim_regrid = regrid_horizontal(
        data=sss_clim, target="woa23", method="conservative"
    )
    sss_hist_clim_regrid[model] = sss_clim_regrid

    comp_mask = build_comparison_mask(
        fields=[woa_sss, sss_clim_regrid], ocean_only=True
    )
    model_on = sss_clim_regrid.where(comp_mask)
    woa_on = woa_sss.where(comp_mask)

    metrics = compute_difference_metrics(
        data=model_on, reference=woa_on, dim=["lat", "lon"]
    )
    bias_fields[model] = metrics["difference"]
    rmse[model] = val(metrics["rmse"])
    bias_scalar[model] = val(metrics["bias"])

    corr = analyze_correlation(
        var1=model_on, var2=woa_on, dim=["lat", "lon"], method="pearson"
    )
    correlation[model] = val(corr["correlation"])

    # Monthly global-mean SSS and annual series
    sss_hist_global_monthly = compute_basic_stat(
        data=sss_hist, op="area_weighted_mean", dim=["lat", "lon"]
    )
    annual_hist = temporal_aggregate(
        data=sss_hist_global_monthly, freq="year", op="mean"
    )
    annual_hist_global[model] = annual_hist
    hist_fit[model] = fit_1D_trend(data=annual_hist, method="linear")

    sss_ssp_global_monthly = compute_basic_stat(
        data=sss_ssp, op="area_weighted_mean", dim=["lat", "lon"]
    )
    annual_ssp = temporal_aggregate(
        data=sss_ssp_global_monthly, freq="year", op="mean"
    )
    annual_ssp_global[model] = annual_ssp
    ssp_fit[model] = fit_1D_trend(data=annual_ssp, method="linear")

    hist_tr = hist_fit[model]["trend"]
    ssp_tr = ssp_fit[model]["trend"]

    hist_slope_dec[model] = float(hist_tr["slope"]) * 10.0
    ssp_slope_dec[model] = float(ssp_tr["slope"]) * 10.0
    hist_ci_dec[model] = (
        float(hist_tr["ci95_low"]) * 10.0,
        float(hist_tr["ci95_high"]) * 10.0,
    )
    ssp_ci_dec[model] = (
        float(ssp_tr["ci95_low"]) * 10.0,
        float(ssp_tr["ci95_high"]) * 10.0,
    )
    hist_r2[model] = float(hist_fit[model]["metrics"]["r2"])
    ssp_r2[model] = float(ssp_fit[model]["metrics"]["r2"])
    hist_sig[model] = bool(hist_tr["significant"])
    ssp_sig[model] = bool(ssp_tr["significant"])

    hist_stable[model] = bool(
        hist_sig[model]
        and abs(hist_slope_dec[model]) > 0
        and (hist_ci_dec[model][1] - hist_ci_dec[model][0])
        < 2.0 * abs(hist_slope_dec[model])
    )
    slope_ratio[model] = (
        ssp_slope_dec[model] / hist_slope_dec[model]
        if hist_slope_dec[model] != 0 else np.nan
    )
    ratio_meaningful[model] = bool(hist_sig[model] and hist_stable[model])

    # Cell-by-cell historical trend map
    annual_sss_hist_map = temporal_aggregate(
        data=sss_hist, freq="year", op="mean"
    )
    tmap = fit_trend_map(
        data=annual_sss_hist_map, dim="year", per_decade=True
    )
    trend_maps[model] = tmap

    trend_slope_regrid[model] = regrid_horizontal(
        data=tmap["slope"], target="woa23", method="conservative"
    )
    trend_p_regrid[model] = regrid_horizontal(
        data=tmap["p_value"], target="woa23", method="nearest_s2d"
    )

    # Delta SSS: ssp245 2071-2100 minus historical 1995-2014
    sss_ssp_mean = temporal_aggregate(data=sss_ssp, freq="mean", op="mean")
    sss_hist_9514 = sss_hist.sel(time=slice("1995-01", "2014-12"))
    sss_hist_mean_9514 = temporal_aggregate(
        data=sss_hist_9514, freq="mean", op="mean"
    )
    delta = sss_ssp_mean - sss_hist_mean_9514
    delta = delta.assign_attrs(sss_hist.attrs)
    delta.attrs["units"] = "psu"
    delta_fields[model] = delta

    delta_global[model] = val(
        compute_basic_stat(data=delta, op="area_weighted_mean", dim=["lat", "lon"])
    )

    delta_basin[model] = {}
    for basin in ["Atlantic", "Indian", "Pacific"]:
        d_b = select_basin(
            data=delta, basin=basin, include_marginal_seas=True
        )
        delta_basin[model][basin] = val(
            compute_basic_stat(
                data=d_b, op="area_weighted_mean", dim=["lat", "lon"]
            )
        )

    ensemble_labels[model] = {
        "historical": ens_label(so_hist),
        "ssp245": ens_label(so_ssp),
    }

summary = {
    "woa_global_mean_psu": woa_global_mean,
    "models": {},
}
for model in models:
    summary["models"][model] = {
        "global_mean_1971_2000_psu": global_mean[model],
        "difference_from_woa_psu": diff_to_woa[model],
        "rmse_vs_woa_psu": rmse[model],
        "correlation_vs_woa": correlation[model],
        "historical_trend_psu_per_decade": hist_slope_dec[model],
        "historical_trend_ci95_psu_per_decade": hist_ci_dec[model],
        "historical_r2": hist_r2[model],
        "historical_significant": hist_sig[model],
        "ssp245_trend_psu_per_decade": ssp_slope_dec[model],
        "ssp245_trend_ci95_psu_per_decade": ssp_ci_dec[model],
        "ssp245_r2": ssp_r2[model],
        "ssp245_significant": ssp_sig[model],
        "ssp245_over_historical_slope_ratio": slope_ratio[model],
        "ratio_meaningful": ratio_meaningful[model],
        "delta_sss_global_psu": delta_global[model],
        "delta_sss_basins_psu": delta_basin[model],
        "ensemble_labels": ensemble_labels[model],
        "sss_levels_used": depth_levels_used[model],
    }

print("WOA23 1971-2000 global mean SSS: {:.2f} psu".format(woa_global_mean))
for model in models:
    print("\n{}".format(model))
    print("  1971-2000 global mean SSS: {:.2f} psu; difference from WOA23: {:+.2f} psu".format(
        global_mean[model], diff_to_woa[model]
    ))
    print("  RMSE vs WOA23: {:.2f} psu; spatial Pearson r: {:.3f}".format(
        rmse[model], correlation[model]
    ))
    print("  historical 1950-2014 trend: {:.3g} psu/decade (95% CI {:.3g} to {:.3g}), R2={:.3g}, significant={}".format(
        hist_slope_dec[model], hist_ci_dec[model][0], hist_ci_dec[model][1],
        hist_r2[model], hist_sig[model]
    ))
    print("  ssp245 2071-2100 trend: {:.3g} psu/decade (95% CI {:.3g} to {:.3g}), R2={:.3g}, significant={}".format(
        ssp_slope_dec[model], ssp_ci_dec[model][0], ssp_ci_dec[model][1],
        ssp_r2[model], ssp_sig[model]
    ))
    print("  ssp245/historical slope ratio: {:.1f} ({})".format(
        slope_ratio[model],
        "meaningful" if ratio_meaningful[model] else "not meaningful/unstable"
    ))
    print("  delta SSS 2071-2100 minus 1995-2014: global {:.2f} psu; Atlantic {:.2f}, Indian {:.2f}, Pacific {:.2f}".format(
        delta_global[model], delta_basin[model]["Atlantic"],
        delta_basin[model]["Indian"], delta_basin[model]["Pacific"]
    ))
    print("  ensembles: historical={}, ssp245={}".format(
        ensemble_labels[model]["historical"], ensemble_labels[model]["ssp245"]
    ))
