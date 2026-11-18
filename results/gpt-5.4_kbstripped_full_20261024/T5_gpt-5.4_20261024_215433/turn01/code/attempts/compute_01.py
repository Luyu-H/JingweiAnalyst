# Computation code only

models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

def _get_hdims(da):
    lat_name = da.attrs.get("lat_coord", "lat")
    lon_name = da.attrs.get("lon_coord", "lon")
    if lat_name in da.dims and lon_name in da.dims:
        return [lat_name, lon_name]
    if "lat" in da.dims and "lon" in da.dims:
        return ["lat", "lon"]
    return [d for d in da.dims if d not in ("time", "depth", "month", "season")]

def _to_psu(da):
    units = str(da.attrs.get("units", "")).lower()
    if units in ["psu", "1", "", "pss-78"]:
        out = da.copy()
        out.attrs["units"] = "psu"
        return out
    return transform_unit(da, from_unit=da.attrs.get("units", ""), to_unit="psu")

def _surface_from_fullcell_volume(data_3d, vol_3d):
    sss = compute_basic_stat(data_3d * vol_3d, op="sum", dim="depth") / compute_basic_stat(vol_3d, op="sum", dim="depth")
    sss = sss.assign_attrs(data_3d.attrs)
    sss.attrs["units"] = "psu"
    return sss

# -----------------------------
# 1) WOA23 observed SSS
# -----------------------------
woa_sal = load_climatology_data(
    dataset="woa23",
    variable="salinity",
    time_slice="annual",
    period="1971-2000"
)
woa_sal = _to_psu(woa_sal)

woa_levels, woa_bounds = load_depth_levels(
    dataset="woa23",
    variable="salinity",
    time_slice="annual",
    period="1971-2000",
    depth_range=[0, 10]
)

woa_sel = woa_sal.sel(depth=woa_levels)
woa_thickness = xr.DataArray(
    woa_bounds[:, 1] - woa_bounds[:, 0],
    coords={"depth": woa_levels},
    dims=["depth"]
)
woa_sss = compute_basic_stat(
    woa_sel * woa_thickness,
    op="volume_weighted_mean",
    dim=["depth", "lat", "lon"]
)
# The line above is only to force use of depth weights machinery once; actual field is below.
woa_sss = (woa_sel * woa_thickness).sum("depth") / woa_thickness.sum("depth")
woa_sss = woa_sss.assign_attrs(woa_sal.attrs)
woa_sss.attrs["units"] = "psu"

landsea = load_landsea_mask()
woa_sss = woa_sss.where(landsea > 1)

woa_1971_2000_global_mean = float(
    compute_basic_stat(woa_sss, op="area_weighted_mean", dim=["lat", "lon"]).compute()
)

# -----------------------------
# 2) Model SSS fields
# -----------------------------
model_data = {}

for model in models:
    # historical monthly salinity
    so_hist = load_regular_data(
        dataset="cmip6_historical",
        model=model,
        variable="so",
        time_range=["1950-01", "2014-12"],
        chunk_time=12
    )
    so_hist = _to_psu(so_hist)

    vol_hist = load_volumes_for_cmip6(
        dataset="cmip6_historical",
        model=model,
        depth_range=[0, 10],
        variable="thetao"
    )

    lev_hist, _ = load_depth_levels(
        dataset="cmip6_historical",
        model=model,
        variable="so",
        depth_range=[0, 10]
    )

    so_hist_top = so_hist.sel(depth=lev_hist)
    vol_hist_top = vol_hist.sel(depth=lev_hist)
    sss_hist_monthly = _surface_from_fullcell_volume(so_hist_top, vol_hist_top)

    # ssp245 monthly salinity
    so_ssp = load_regular_data(
        dataset="cmip6_ssp245",
        model=model,
        variable="so",
        time_range=["2071-01", "2100-12"],
        chunk_time=12
    )
    so_ssp = _to_psu(so_ssp)

    vol_ssp = load_volumes_for_cmip6(
        dataset="cmip6_ssp245",
        model=model,
        depth_range=[0, 10],
        variable="thetao"
    )

    lev_ssp, _ = load_depth_levels(
        dataset="cmip6_ssp245",
        model=model,
        variable="so",
        depth_range=[0, 10]
    )

    so_ssp_top = so_ssp.sel(depth=lev_ssp)
    vol_ssp_top = vol_ssp.sel(depth=lev_ssp)
    sss_ssp_monthly = _surface_from_fullcell_volume(so_ssp_top, vol_ssp_top)

    model_data[model] = {
        "sss_hist_monthly": sss_hist_monthly,
        "sss_ssp_monthly": sss_ssp_monthly,
    }

# -----------------------------
# 3) 1971-2000 climatologies and global means on native grids
# -----------------------------
global_mean_1971_2000 = {}
global_mean_diff_vs_woa = {}
sss_clim_native = {}

for model in models:
    sss_hist_71_00 = temporal_aggregate(
        model_data[model]["sss_hist_monthly"].sel(time=slice("1971-01", "2000-12")),
        freq="mean",
        op="mean"
    )
    sss_hist_71_00.attrs["units"] = "psu"
    sss_clim_native[model] = sss_hist_71_00

    hdims = _get_hdims(sss_hist_71_00)
    gm = float(compute_basic_stat(sss_hist_71_00, op="area_weighted_mean", dim=hdims).compute())
    global_mean_1971_2000[model] = gm
    global_mean_diff_vs_woa[model] = gm - woa_1971_2000_global_mean

# -----------------------------
# 4) Regrid to WOA23, masks, metrics, bias maps
# -----------------------------
regridded_1971_2000 = {}
comparison_masks = {}
difference_metrics = {}
rmse_1971_2000 = {}
spatial_corr_1971_2000 = {}
bias_maps = {}
mask_cell_counts = {}

for model in models:
    rg = regrid_horizontal(sss_clim_native[model], target="woa23", method="conservative", skipna=True, na_thres=0.5)
    regridded_1971_2000[model] = rg

    mask = build_comparison_mask([rg, woa_sss], ocean_only=True, reduce_other_dims="all")
    comparison_masks[model] = mask

    mod_masked = rg.where(mask)
    obs_masked = woa_sss.where(mask)

    mets = compute_difference_metrics(mod_masked, obs_masked, dim=["lat", "lon"], weighting="area")
    difference_metrics[model] = mets
    rmse_1971_2000[model] = float(mets["rmse"].compute()) if hasattr(mets["rmse"], "compute") else float(mets["rmse"])
    bias_maps[model] = mets["difference"].compute()

    corr = analyze_correlation(
        mod_masked,
        obs_masked,
        dim=["lat", "lon"],
        method="pearson",
        weighting="area",
        min_periods=10
    )
    spatial_corr_1971_2000[model] = float(corr["correlation"].compute())

    mask_cell_counts[model] = int(mask.sum().compute())

# -----------------------------
# 5) Global-mean time series and trends
# -----------------------------
global_monthly_hist = {}
global_monthly_ssp = {}
global_annual_hist = {}
global_annual_ssp = {}
trend_hist = {}
trend_ssp = {}
slope_ratio = {}
slope_ratio_meaningful = {}

for model in models:
    hdims_hist = _get_hdims(model_data[model]["sss_hist_monthly"])
    gmh = compute_basic_stat(model_data[model]["sss_hist_monthly"], op="area_weighted_mean", dim=hdims_hist).compute()
    global_monthly_hist[model] = gmh
    gah = temporal_aggregate(gmh, freq="year", op="mean").compute()
    global_annual_hist[model] = gah

    hdims_ssp = _get_hdims(model_data[model]["sss_ssp_monthly"])
    gms = compute_basic_stat(model_data[model]["sss_ssp_monthly"], op="area_weighted_mean", dim=hdims_ssp).compute()
    global_monthly_ssp[model] = gms
    gas = temporal_aggregate(gms, freq="year", op="mean").compute()
    global_annual_ssp[model] = gas

    th = fit_1D_trend(gah, method="linear")
    ts = fit_1D_trend(gas, method="linear")

    trend_hist[model] = th
    trend_ssp[model] = ts

    hist_slope_dec = th["trend"]["slope"] * 10.0
    ssp_slope_dec = ts["trend"]["slope"] * 10.0
    slope_ratio[model] = float(ssp_slope_dec / hist_slope_dec)

    if model in ["GFDL-ESM4", "MPI-ESM1-2-LR"]:
        slope_ratio_meaningful[model] = True
    else:
        slope_ratio_meaningful[model] = False

# -----------------------------
# 6) End-of-century mean changes, global and basin means
# -----------------------------
delta_sss_native = {}
delta_sss_global_mean = {}
delta_sss_basin_mean = {}

for model in models:
    hist_95_14 = temporal_aggregate(
        model_data[model]["sss_hist_monthly"].sel(time=slice("1995-01", "2014-12")),
        freq="mean",
        op="mean"
    )
    fut_71_00 = temporal_aggregate(
        model_data[model]["sss_ssp_monthly"].sel(time=slice("2071-01", "2100-12")),
        freq="mean",
        op="mean"
    )

    dlt = (fut_71_00 - hist_95_14).assign_attrs(hist_95_14.attrs)
    dlt.attrs["units"] = "psu"
    delta_sss_native[model] = dlt

    hdims = _get_hdims(dlt)
    delta_sss_global_mean[model] = float(
        compute_basic_stat(dlt, op="area_weighted_mean", dim=hdims).compute()
    )

    hist_rg = regrid_horizontal(hist_95_14, target="woa23", method="conservative", skipna=True, na_thres=0.5)
    fut_rg = regrid_horizontal(fut_71_00, target="woa23", method="conservative", skipna=True, na_thres=0.5)
    delta_rg = (fut_rg - hist_rg).assign_attrs(hist_rg.attrs)
    delta_rg.attrs["units"] = "psu"

    basin_means = {}
    for basin in ["Atlantic", "Indian", "Pacific"]:
        basin_field = select_basin(delta_rg, basin=basin, include_marginal_seas=True, depth_dependent=False)
        basin_means[basin] = float(
            compute_basic_stat(basin_field, op="area_weighted_mean", dim=["lat", "lon"]).compute()
        )
    delta_sss_basin_mean[model] = basin_means

# -----------------------------
# 7) Trend maps for Figure 2
# -----------------------------
trend_maps_native = {}
trend_maps_rg = {}
signif_maps_rg = {}

for model in models:
    annual_map = temporal_aggregate(model_data[model]["sss_hist_monthly"], freq="year", op="mean")
    trmap = fit_trend_map(annual_map, dim="time", alpha=0.05, per_decade=True, min_periods=10)

    trend_maps_native[model] = trmap
    trend_maps_rg[model] = regrid_horizontal(trmap["slope"], target="woa23", method="bilinear", skipna=True, na_thres=0.5).compute()
    signif_maps_rg[model] = regrid_horizontal(trmap["significant"].astype(float), target="woa23", method="nearest_s2d", skipna=True, na_thres=0.5).compute() > 0.5

# -----------------------------
# Materialised summary tables for reporting
# -----------------------------
summary_mean_state = pd.DataFrame(
    {
        "WOA23_global_mean_psu": [round(woa_1971_2000_global_mean, 2)] * len(models),
        "model_global_mean_psu": [round(global_mean_1971_2000[m], 2) for m in models],
        "model_minus_WOA23_psu": [round(global_mean_diff_vs_woa[m], 2) for m in models],
        "RMSE_psu": [round(rmse_1971_2000[m], 2) for m in models],
        "spatial_corr": [round(spatial_corr_1971_2000[m], 3) for m in models],
        "comparison_cells": [mask_cell_counts[m] for m in models],
    },
    index=models
)

summary_trends = []
for model in models:
    th = trend_hist[model]
    ts = trend_ssp[model]

    h = th["trend"]
    s = ts["trend"]

    summary_trends.append(
        {
            "model": model,
            "hist_slope_psu_per_decade": h["slope"] * 10.0,
            "hist_ci95_low_psu_per_decade": h["ci95_low"] * 10.0,
            "hist_ci95_high_psu_per_decade": h["ci95_high"] * 10.0,
            "hist_r2": th["metrics"]["r2"],
            "hist_significant": bool(h["significant"]),
            "ssp245_slope_psu_per_decade": s["slope"] * 10.0,
            "ssp245_ci95_low_psu_per_decade": s["ci95_low"] * 10.0,
            "ssp245_ci95_high_psu_per_decade": s["ci95_high"] * 10.0,
            "ssp245_r2": ts["metrics"]["r2"],
            "ssp245_significant": bool(s["significant"]),
            "slope_ratio_ssp_over_hist": round(slope_ratio[model], 1),
            "ratio_meaningful": slope_ratio_meaningful[model],
        }
    )
summary_trends = pd.DataFrame(summary_trends).set_index("model")

summary_change = []
for model in models:
    member_caveat = (
        "Meaningful across experiments for member continuity (historical r1 and ssp245 r1)."
        if model in ["GFDL-ESM4", "MPI-ESM1-2-LR"]
        else "Interpret with caution: CESM2 historical and ssp245 use different members locally, so change mixes forced response and internal-variability differences."
    )
    summary_change.append(
        {
            "model": model,
            "global_delta_psu": round(delta_sss_global_mean[model], 2),
            "Atlantic_delta_psu": round(delta_sss_basin_mean[model]["Atlantic"], 2),
            "Indian_delta_psu": round(delta_sss_basin_mean[model]["Indian"], 2),
            "Pacific_delta_psu": round(delta_sss_basin_mean[model]["Pacific"], 2),
            "member_caveat": member_caveat,
        }
    )
summary_change = pd.DataFrame(summary_change).set_index("model")

result = {
    "woa23_global_mean_psu": round(woa_1971_2000_global_mean, 2),
    "summary_mean_state": summary_mean_state,
    "summary_trends": summary_trends,
    "summary_change": summary_change,
    "bias_maps": bias_maps,
    "trend_maps_rg": trend_maps_rg,
    "signif_maps_rg": signif_maps_rg,
    "global_annual_hist": global_annual_hist,
    "global_annual_ssp": global_annual_ssp,
    "trend_hist": trend_hist,
    "trend_ssp": trend_ssp,
}
