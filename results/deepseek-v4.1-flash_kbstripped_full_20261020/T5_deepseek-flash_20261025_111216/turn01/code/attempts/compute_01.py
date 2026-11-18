# ============================================================
# SSS in three CMIP6 models: mean-state bias, historical trend,
# and end-of-century change under ssp245.
# ============================================================

models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

# ------------------------------------------------------------
# Helper: materialise a scalar
# ------------------------------------------------------------
def to_float(x):
    if hasattr(x, "compute"):
        x = x.compute()
    return float(x)

# ------------------------------------------------------------
# Helper: area-weighted global mean on the field's own grid
# ------------------------------------------------------------
def global_mean(field):
    return to_float(
        compute_basic_stat(field, op="area_weighted_mean", dim=["lat", "lon"])
    )

# ------------------------------------------------------------
# Helper: area-weighted basin mean using the WOA basin mask
# ------------------------------------------------------------
def basin_mean(field, basin):
    field_basin = select_basin(
        field, basin=basin, include_marginal_seas=True, depth_dependent=False
    )
    return to_float(
        compute_basic_stat(field_basin, op="area_weighted_mean", dim=["lat", "lon"])
    )

# ------------------------------------------------------------
# Helper: SSS from levels whose depths lie in 0-10 m,
# each level weighted by its full cell volume.
# ------------------------------------------------------------
def make_sss(so, vol, sel_levels):
    so_sel = so.sel(depth=sel_levels)
    vol_sel = vol.sel(depth=sel_levels)
    if "time" in vol_sel.dims:
        vol_sel = vol_sel.sel(time=so_sel.time)
    num = (so_sel * vol_sel).sum("depth", skipna=True)
    den = vol_sel.where(so_sel.notnull()).sum("depth", skipna=True)
    sss = num / den
    sss.attrs = so.attrs.copy()
    sss.attrs["units"] = "psu"
    sss.attrs["long_name"] = "Sea surface salinity"
    return sss

# ============================================================
# 1. WOA23 reference SSS, 1971-2000 annual climatology
# ============================================================
woa_sal = load_climatology_data(
    dataset="woa23",
    variable="salinity",
    time_slice="annual",
    period="1971-2000",
)

woa_levels, woa_bounds = load_depth_levels(
    dataset="woa23",
    variable="salinity",
    time_slice="annual",
    period="1971-2000",
)

woa_sel_idx = np.where((woa_levels >= 0.0) & (woa_levels <= 10.0))[0]
woa_sel_levels = woa_levels[woa_sel_idx]
woa_thick = np.abs(woa_bounds[woa_sel_idx, 1] - woa_bounds[woa_sel_idx, 0])

woa_area = load_areas(dataset="woa23")
woa_vol = woa_area * xr.DataArray(
    woa_thick, dims="depth", coords={"depth": woa_sel_levels}
)

woa_sal_sel = woa_sal.sel(depth=woa_sel_levels)
woa_num = (woa_sal_sel * woa_vol).sum("depth", skipna=True)
woa_den = woa_vol.where(woa_sal_sel.notnull()).sum("depth", skipna=True)
woa_sss = woa_num / woa_den
woa_sss.attrs = woa_sal.attrs.copy()
woa_sss.attrs["units"] = "psu"
woa_sss.attrs["long_name"] = "Sea surface salinity"

# Keep only ocean cells using the WOA 1-degree land-sea mask
landsea = load_landsea_mask()
woa_sss = woa_sss.where(landsea > 1)

woa_gm_1971_2000 = global_mean(woa_sss)

# Optional basin-mask reference (surface labels)
basin_mask_surface = load_basin_mask(depth_levels="0m")

# ============================================================
# 2-3. Model SSS: historical 1950-2014 and ssp245 2071-2100
# ============================================================
model_results = {}

for model in models:
    # ---- load salinity and volumes -------------------------
    so_hist = load_regular_data(
        dataset="cmip6_historical",
        variable="so",
        model=model,
        time_range=["1950-01", "2014-12"],
    )
    so_ssp = load_regular_data(
        dataset="cmip6_ssp245",
        variable="so",
        model=model,
        time_range=["2071-01", "2100-12"],
    )

    levels, bounds = load_depth_levels(
        dataset="cmip6_historical", model=model, variable="so"
    )
    sel_idx = np.where((levels >= 0.0) & (levels <= 10.0))[0]
    sel_levels = levels[sel_idx]

    vol_hist = load_volumes_for_cmip6(
        dataset="cmip6_historical", model=model, variable="so"
    )
    vol_ssp = load_volumes_for_cmip6(
        dataset="cmip6_ssp245", model=model, variable="so"
    )
    if "time" in vol_hist.dims:
        vol_hist = vol_hist.sel(time=so_hist.time)
    if "time" in vol_ssp.dims:
        vol_ssp = vol_ssp.sel(time=so_ssp.time)

    # area is not strictly needed for the volume weighting here,
    # but loading it keeps the model grid identity explicit.
    area_hist = load_areas(
        dataset="cmip6_historical", model=model, variable="so"
    )

    # ---- monthly SSS = volume-weighted mean over 0-10 m ----
    sss_hist_monthly = make_sss(so_hist, vol_hist, sel_levels)
    sss_ssp_monthly = make_sss(so_ssp, vol_ssp, sel_levels)

    # ---- annual means --------------------------------------
    sss_hist_annual = temporal_aggregate(
        sss_hist_monthly, freq="year", op="mean"
    )
    sss_ssp_annual = temporal_aggregate(
        sss_ssp_monthly, freq="year", op="mean"
    )

    # ---- period means --------------------------------------
    sss_hist_1971_2000 = temporal_aggregate(
        sss_hist_monthly.sel(time=slice("1971-01", "2000-12")),
        freq="mean",
        op="mean",
    )
    sss_hist_1995_2014 = temporal_aggregate(
        sss_hist_monthly.sel(time=slice("1995-01", "2014-12")),
        freq="mean",
        op="mean",
    )
    sss_ssp_2071_2100 = temporal_aggregate(
        sss_ssp_monthly.sel(time=slice("2071-01", "2100-12")),
        freq="mean",
        op="mean",
    )

    # ---- global means --------------------------------------
    gm_hist_1971_2000 = global_mean(sss_hist_1971_2000)
    gm_ssp_2071_2100 = global_mean(sss_ssp_2071_2100)
    gm_hist_1995_2014 = global_mean(sss_hist_1995_2014)

    # ---- end-of-century change -----------------------------
    delta = sss_ssp_2071_2100 - sss_hist_1995_2014
    delta.attrs = sss_ssp_2071_2100.attrs.copy()
    delta.attrs["units"] = "psu"
    delta.attrs["long_name"] = "ssp245 minus historical SSS"

    gm_delta = global_mean(delta)
    basin_delta = {
        b: basin_mean(delta, b) for b in ["Atlantic", "Indian", "Pacific"]
    }

    # ---- annual global-mean series -------------------------
    gm_hist_series = compute_basic_stat(
        sss_hist_annual, op="area_weighted_mean", dim=["lat", "lon"]
    ).compute()
    gm_ssp_series = compute_basic_stat(
        sss_ssp_annual, op="area_weighted_mean", dim=["lat", "lon"]
    ).compute()

    # ---- trends of the global-mean series ------------------
    trend_hist = fit_1D_trend(gm_hist_series, method="linear")
    trend_ssp = fit_1D_trend(gm_ssp_series, method="linear")

    # ---- cell-by-cell historical trend map -----------------
    time_dim = "year" if "year" in sss_hist_annual.dims else "time"
    trend_map = fit_trend_map(
        sss_hist_annual, dim=time_dim, per_decade=True
    )

    # regrid for display only
    slope_regrid = regrid_horizontal(
        trend_map["slope"], target="woa23", method="conservative"
    )
    sig_regrid = (
        regrid_horizontal(
            trend_map["significant"].astype(float),
            target="woa23",
            method="nearest_s2d",
        )
        > 0.5
    )

    model_results[model] = {
        "sss_hist_1971_2000": sss_hist_1971_2000,
        "sss_ssp_2071_2100": sss_ssp_2071_2100,
        "sss_hist_1995_2014": sss_hist_1995_2014,
        "gm_hist_1971_2000": gm_hist_1971_2000,
        "gm_ssp_2071_2100": gm_ssp_2071_2100,
        "gm_hist_1995_2014": gm_hist_1995_2014,
        "gm_delta": gm_delta,
        "basin_delta": basin_delta,
        "gm_hist_series": gm_hist_series,
        "gm_ssp_series": gm_ssp_series,
        "trend_hist": trend_hist,
        "trend_ssp": trend_ssp,
        "trend_map": trend_map,
        "slope_regrid": slope_regrid,
        "sig_regrid": sig_regrid,
        "delta": delta,
    }

# ============================================================
# 5. Bias, RMSE and pattern correlation against WOA23
# ============================================================
model_regrids = {
    m: regrid_horizontal(
        model_results[m]["sss_hist_1971_2000"],
        target="woa23",
        method="conservative",
    )
    for m in models
}

comparison_fields = [model_regrids[m] for m in models] + [woa_sss]
comparison_mask = build_comparison_mask(comparison_fields, ocean_only=True)
mask_count = int(comparison_mask.sum().compute())

bias_maps = {}
metrics_dict = {}
corr_dict = {}

for m in models:
    bias_maps[m] = (model_regrids[m] - woa_sss).where(comparison_mask).compute()

    met = compute_difference_metrics(
        data=model_regrids[m].where(comparison_mask),
        reference=woa_sss.where(comparison_mask),
        dim=["lat", "lon"],
    )
    metrics_dict[m] = {
        "bias": to_float(met["bias"]),
        "rmse": to_float(met["rmse"]),
        "mae": to_float(met["mae"]),
        "max_abs_error": to_float(met["max_abs_error"]),
    }

    corr = analyze_correlation(
        var1=model_regrids[m].where(comparison_mask),
        var2=woa_sss.where(comparison_mask),
        dim=["lat", "lon"],
        method="pearson",
        weighting="area",
    )
    corr_dict[m] = to_float(corr["correlation"])

# ============================================================
# Report the requested numbers
# ============================================================
print("=" * 72)
print("WOA23 reference")
print("=" * 72)
print(f"WOA23 global mean 1971-2000 SSS: {woa_gm_1971_2000:.2f} psu")
print(f"Common comparison mask has {mask_count} cells")

print("\n" + "=" * 72)
print("Mean state 1971-2000, RMSE and pattern correlation")
print("=" * 72)
for m in models:
    gm = model_results[m]["gm_hist_1971_2000"]
    diff = gm - woa_gm_1971_2000
    print(f"\n{m}")
    print(f"  global mean SSS 1971-2000 : {gm:.2f} psu")
    print(f"  difference from WOA23     : {diff:+.2f} psu")
    print(f"  RMSE vs WOA23            : {metrics_dict[m]['rmse']:.2f} psu")
    print(f"  pattern correlation      : {corr_dict[m]:.3f}")

print("\n" + "=" * 72)
print("Trends of the annual global-mean SSS")
print("=" * 72)
for m in models:
    th = model_results[m]["trend_hist"]
    ts = model_results[m]["trend_ssp"]

    slope_h = th["trend"]["slope"] * 10.0
    ci_h_low = th["trend"]["ci95_low"] * 10.0
    ci_h_high = th["trend"]["ci95_high"] * 10.0
    r2_h = th["metrics"]["r2"]
    sig_h = th["trend"]["significant"]

    slope_s = ts["trend"]["slope"] * 10.0
    ci_s_low = ts["trend"]["ci95_low"] * 10.0
    ci_s_high = ts["trend"]["ci95_high"] * 10.0
    r2_s = ts["metrics"]["r2"]
    sig_s = ts["trend"]["significant"]

    ratio = slope_s / slope_h

    print(f"\n{m}")
    print(
        f"  historical 1950-2014 : {slope_h:.3g} psu/decade "
        f"(95% CI {ci_h_low:.3g} to {ci_h_high:.3g}), "
        f"R2={r2_h:.3f}, CI excludes zero: {sig_h}"
    )
    print(
        f"  ssp245 2071-2100    : {slope_s:.3g} psu/decade "
        f"(95% CI {ci_s_low:.3g} to {ci_s_high:.3g}), "
        f"R2={r2_s:.3f}, CI excludes zero: {sig_s}"
    )
    print(f"  ratio ssp245/historical slope: {ratio:.1f}")

print("\n" + "=" * 72)
print("End-of-century change: ssp245 2071-2100 minus historical 1995-2014")
print("=" * 72)
for m in models:
    print(f"\n{m}")
    print(f"  global  : {model_results[m]['gm_delta']:+.2f} psu")
    for b in ["Atlantic", "Indian", "Pacific"]:
        print(f"  {b:8s}: {model_results[m]['basin_delta'][b]:+.2f} psu")

ensemble_caveat = (
    "Ensemble caveat: CESM2 ssp245 is r10i1p1f1, while only the historical r1 "
    "member is available locally. The CESM2 delta SSS therefore contains the "
    "forced change plus the internal-variability difference between two "
    "different members. GFDL-ESM4 and MPI-ESM1-2-LR use r1 for both experiments "
    "and are continuous. No trend is fitted across the 2014/2015 join."
)
print("\n" + ensemble_caveat)

print(
    "\nThe slope ratio compares a 65-year historical rate with a 30-year "
    "late-century scenario rate over different periods and forcings; for CESM2 "
    "it also compares different ensemble members. It is therefore indicative "
    "only, not a like-for-like acceleration factor."
)
