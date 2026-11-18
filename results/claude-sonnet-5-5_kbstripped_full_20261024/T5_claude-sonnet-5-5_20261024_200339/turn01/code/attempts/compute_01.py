models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
colors = {"CESM2": "tab:blue", "GFDL-ESM4": "tab:orange", "MPI-ESM1-2-LR": "tab:green"}

# ---------------------------------------------------------------- 1. WOA23 reference
woa_so = load_climatology_data(dataset="woa23", variable="salinity", period="1971-2000",
                               time_slice="annual", depth_range=[0, 10])
woa_levels, woa_bounds = load_depth_levels(dataset="woa23", variable="salinity",
                                           depth_range=[0, 10], period="1971-2000")
landsea = load_landsea_mask()
basin_mask = load_basin_mask(depth_levels="0m")
basin_names = basin_mask.attrs.get("basin_names")
woa_areas = load_areas()

# SSS = mean of the 0, 5, 10 m levels, each weighted by its full cell volume
woa_sss = compute_basic_stat(woa_so, op="volume_weighted_mean", dim="depth").compute()
woa_sss = woa_sss.assign_attrs(units="psu")

# ---------------------------------------------------------------- 2. Model data
sss_hist, sss_ssp = {}, {}
grid_info = {}
for m in models:
    levs, bnds = load_depth_levels(dataset="cmip6_historical", model=m, variable="so",
                                   depth_range=[0, 10])
    a = load_areas(dataset="cmip6_historical", model=m, variable="so")
    v = load_volumes_for_cmip6(dataset="cmip6_historical", model=m, variable="so",
                               depth_range=[0, 10])
    grid_info[m] = {
        "levels_m": np.asarray(levs).tolist(),
        "layer_bounds_m": np.asarray(bnds).tolist(),
        "area_source": a.attrs.get("area_source"),
        "volume_source": v.attrs.get("volume_source"),
        "volume_dims": list(v.dims),
    }
    so_h = load_regular_data(dataset="cmip6_historical", variable="so", model=m,
                             time_range=["1950-01", "2014-12"], depth_range=[0, 10])
    so_s = load_regular_data(dataset="cmip6_ssp245", variable="so", model=m,
                             time_range=["2071-01", "2100-12"], depth_range=[0, 10])
    # volume-weighted mean over the 0-10 m levels (full cell volumes)
    sss_hist[m] = compute_basic_stat(so_h, op="volume_weighted_mean", dim="depth").assign_attrs(units="psu")
    sss_ssp[m] = compute_basic_stat(so_s, op="volume_weighted_mean", dim="depth").assign_attrs(units="psu")

# ---------------------------------------------------------------- 3. Mean-state bias 1971-2000
clim_native, clim_1deg = {}, {}
for m in models:
    c = temporal_aggregate(sss_hist[m].sel(time=slice("1971-01-01", "2000-12-31")),
                           freq="mean").compute()
    clim_native[m] = c
    clim_1deg[m] = regrid_horizontal(c, target="woa23", method="conservative")

cmask = build_comparison_mask([clim_1deg[m] for m in models] + [woa_sss], ocean_only=True)
mask_attrs = dict(cmask.attrs)
n_shared_cells = int(cmask.sum())

woa_m = woa_sss.where(cmask)
clim_m = {m: clim_1deg[m].where(cmask) for m in models}

bias_rows = []
bias_maps = {}
for m in models:
    dm = compute_difference_metrics(clim_m[m], woa_m)
    corr = analyze_correlation(clim_m[m], woa_m, dim=["lat", "lon"], weighting="area")
    r = float(corr["correlation"])
    bias_maps[m] = dm["difference"].compute()
    bias_rows.append({
        "model": m,
        "model_mean_SSS": round(float(dm["mean_data"]), 2),
        "WOA23_mean_SSS": round(float(dm["mean_reference"]), 2),
        "bias_vs_WOA23": round(float(dm["bias"]), 2),
        "RMSE": round(float(dm["rmse"]), 2),
        "centered_RMSE": round(float(dm["centered_rmse"]), 2),
        "pearson_r": round(r, 3),
        "n_valid": int(dm["n_valid"]),
    })
bias_table = pd.DataFrame(bias_rows)
print("Shared mask: cells =", n_shared_cells, "| attrs:", mask_attrs)
print("Note: means are shared-sample means on the common 1-degree ocean mask, not each source's own global ocean mean.")
print(bias_table.to_string(index=False))

# ---------------------------------------------------------------- 4. Annual means and global-mean trends
ann_hist, ann_ssp, gm_hist, gm_ssp = {}, {}, {}, {}
trend_rows = []
fits = {}
for m in models:
    ah = temporal_aggregate(sss_hist[m], freq="year")
    asp = temporal_aggregate(sss_ssp[m], freq="year")
    print(m, "incomplete hist:", ah.attrs.get("incomplete_periods"),
          "| incomplete ssp:", asp.attrs.get("incomplete_periods"))
    ah = ah.compute()
    asp = asp.compute()
    ann_hist[m], ann_ssp[m] = ah, asp
    gm_hist[m] = compute_basic_stat(ah, op="area_weighted_mean", dim=["lat", "lon"]).compute()
    gm_ssp[m] = compute_basic_stat(asp, op="area_weighted_mean", dim=["lat", "lon"]).compute()

    fh = fit_1D_trend(gm_hist[m], method="linear")
    fs = fit_1D_trend(gm_ssp[m], method="linear")
    fits[m] = {"hist": fh, "ssp": fs}
    th, ts = fh["trend"], fs["trend"]
    ratio = ts["slope"] / th["slope"] if th["slope"] != 0 else np.nan
    # slopes are per YEAR -> multiply by 10 for per decade
    meaningful = bool(th["significant"] and ts["significant"]) and m != "CESM2"
    trend_rows.append({
        "model": m,
        "hist_slope_psu_per_dec": float(f"{th['slope']*10:.3g}"),
        "hist_CI95": (float(f"{th['ci95_low']*10:.3g}"), float(f"{th['ci95_high']*10:.3g}")),
        "hist_R2": round(float(fh["metrics"]["r2"]), 3),
        "hist_CI_excludes_zero": bool(th["significant"]),
        "ssp_slope_psu_per_dec": float(f"{ts['slope']*10:.3g}"),
        "ssp_CI95": (float(f"{ts['ci95_low']*10:.3g}"), float(f"{ts['ci95_high']*10:.3g}")),
        "ssp_R2": round(float(fs["metrics"]["r2"]), 3),
        "ssp_CI_excludes_zero": bool(ts["significant"]),
        "ratio_ssp_over_hist": round(float(ratio), 1),
        "ratio_meaningful": meaningful,
    })
trend_table = pd.DataFrame(trend_rows)
ratio_statement = ("The ssp245/historical slope ratio compares two separately fitted segments "
                   "(1950-2014 vs 2071-2100) that are not a continuous trajectory; it is only "
                   "indicative, and is not meaningful where either confidence interval includes zero "
                   "or where the members differ (CESM2: historical r1i1p1f1, ssp245 r10i1p1f1). "
                   "No fit crosses 2014/2015.")
print(trend_table.to_string(index=False))
print(ratio_statement)

# ---------------------------------------------------------------- 5. Trend maps (1950-2014)
trend_maps, sig_maps = {}, {}
for m in models:
    tm = fit_trend_map(ann_hist[m], dim="year")
    slope_r = regrid_horizontal(tm["slope"], target="woa23", method="conservative")
    sig_r = regrid_horizontal(tm["significant"].astype(float), target="woa23", method="conservative")
    trend_maps[m] = slope_r.compute()
    sig_maps[m] = (sig_r.compute() > 0.5)
    trend_units = tm["units"]

# ---------------------------------------------------------------- 7. Delta SSS
dsss_rows = []
for m in models:
    h95 = temporal_aggregate(sss_hist[m].sel(time=slice("1995-01-01", "2014-12-31")), freq="mean").compute()
    s71 = temporal_aggregate(sss_ssp[m], freq="mean").compute()
    d = (s71 - h95).assign_attrs(units="psu")
    row = {"model": m}
    row["Global"] = round(float(compute_basic_stat(d, op="area_weighted_mean", dim=["lat", "lon"]).compute()), 2)
    for b in ["Atlantic", "Indian", "Pacific"]:
        db = select_basin(d, basin=b, include_marginal_seas=True)
        row[b] = round(float(compute_basic_stat(db, op="area_weighted_mean", dim=["lat", "lon"]).compute()), 2)
    dsss_rows.append(row)
dsss_table = pd.DataFrame(dsss_rows)
member_caveat = ("CESM2: historical r1i1p1f1 but ssp245 r10i1p1f1 (parent r10 absent locally), so its "
                 "delta SSS = forced change + internal-variability difference between two members. "
                 "GFDL-ESM4 and MPI-ESM1-2-LR use r1 in both experiments (continuous join). "
                 "Basins: WOA23 basinmask_01 surface labels, marginal seas folded in, Southern and Arctic excluded; "
                 "Global = all valid ocean cells. Models lacking a marginal sea have no counterpart.")
print(dsss_table.to_string(index=False))
print(member_caveat)
print(grid_info)
