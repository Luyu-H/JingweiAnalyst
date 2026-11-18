import numpy as np
import xarray as xr

MODELS = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
HIST = ["1950-01", "2014-12"]
SSP = ["2071-01", "2100-12"]
BASINS = ["Atlantic", "Indian", "Pacific"]


def depth_name(da):
    for c in ("depth", "lev", "olevel", "zlev"):
        if c in da.dims:
            return c
    raise ValueError("no depth dimension found in " + str(da.dims))


def vol_weighted_sss(field, weight):
    """Within-cell mean over the retained levels, each level weighted by its FULL volume."""
    d = depth_name(field)
    valid = field.notnull()
    w = xr.where(valid, weight, 0.0).fillna(0.0)
    num = (field.fillna(0.0) * w).sum(dim=d)
    den = w.sum(dim=d)
    out = num / den
    attrs = dict(field.attrs)
    attrs.update(units="psu",
                 long_name="sea surface salinity (levels shallower than 10 m, volume weighted)")
    return out.assign_attrs(**attrs)


def with_prov(da, ref):
    return da.assign_attrs(**ref.attrs)


# ---------------------------------------------------------------- observations
woa_sal = load_climatology_data(dataset="woa23", variable="salinity",
                                time_slice="annual", period="1971-2000",
                                depth_range=[0, 10])
lev_w, bnds_w = load_depth_levels(dataset="woa23", variable="salinity",
                                  time_slice="annual", period="1971-2000",
                                  depth_range=[0, 10])
bnds_w = np.asarray(bnds_w)
thick_w = xr.DataArray(bnds_w[:, 1] - bnds_w[:, 0],
                       dims=[depth_name(woa_sal)],
                       coords={depth_name(woa_sal): woa_sal[depth_name(woa_sal)].values},
                       attrs={"units": "m"})

sss_woa = vol_weighted_sss(woa_sal, thick_w).compute()

landsea_mask = load_landsea_mask().compute()
basin_mask = load_basin_mask().compute()
area_woa = load_areas().compute()
ocean_area_woa = float(area_woa.where(sss_woa.notnull()).sum().compute())

print("WOA23 levels used (m):", np.asarray(lev_w))
print("WOA23 layer bounds:", bnds_w.tolist())

# ---------------------------------------------------------------- model fields
clim_hist, clim_ssp, delta, ann_hist, ann_ssp = {}, {}, {}, {}, {}
gm_hist_year, gm_ssp_year, mean_9514 = {}, {}, {}
area_hist, area_ssp = {}, {}
gm_1971_2000, ocean_area, area_src, vol_src = {}, {}, {}, {}

for m in MODELS:
    # ---- historical
    so_h = load_regular_data(dataset="cmip6_historical", variable="so", model=m,
                             time_range=HIST, depth_range=[0, 10])
    vol_h = load_volumes_for_cmip6(dataset="cmip6_historical", model=m,
                                   variable="so", depth_range=[0, 10])
    ar_h = load_areas(dataset="cmip6_historical", model=m, variable="so")
    sss_h = vol_weighted_sss(so_h, vol_h)
    area_src[m] = str(ar_h.attrs.get("area_source", "n/a"))
    vol_src[m] = str(vol_h.attrs.get("volume_source", "n/a"))
    print(m, "historical levels (m):", np.asarray(so_h[depth_name(so_h)].values))

    # ---- ssp245
    so_s = load_regular_data(dataset="cmip6_ssp245", variable="so", model=m,
                             time_range=SSP, depth_range=[0, 10])
    vol_s = load_volumes_for_cmip6(dataset="cmip6_ssp245", model=m,
                                   variable="so", depth_range=[0, 10])
    ar_s = load_areas(dataset="cmip6_ssp245", model=m, variable="so")
    sss_s = vol_weighted_sss(so_s, vol_s)

    # ---- climatologies / means
    clim_hist[m] = with_prov(
        temporal_aggregate(sss_h.sel(time=slice("1971-01", "2000-12")), freq="mean"),
        sss_h).compute()
    mean_9514[m] = with_prov(
        temporal_aggregate(sss_h.sel(time=slice("1995-01", "2014-12")), freq="mean"),
        sss_h).compute()
    clim_ssp[m] = with_prov(temporal_aggregate(sss_s, freq="mean"), sss_s).compute()

    # ---- annual fields and annual global-mean series
    ann_hist[m] = with_prov(temporal_aggregate(sss_h, freq="year"), sss_h)
    ann_ssp[m] = with_prov(temporal_aggregate(sss_s, freq="year"), sss_s)

    gm_hist_year[m] = with_prov(
        compute_basic_stat(ann_hist[m], op="area_weighted_mean", dim=["lat", "lon"]),
        ann_hist[m]).compute()
    gm_ssp_year[m] = with_prov(
        compute_basic_stat(ann_ssp[m], op="area_weighted_mean", dim=["lat", "lon"]),
        ann_ssp[m]).compute()

    # ---- end-of-century change on the native grid
    delta[m] = (clim_ssp[m] - mean_9514[m]).assign_attrs(
        **{**clim_ssp[m].attrs, "units": "psu",
           "long_name": "SSS change ssp245 2071-2100 minus historical 1995-2014"})

    # ---- global 1971-2000 mean on the model's own areas
    gm_1971_2000[m] = float(compute_basic_stat(clim_hist[m], op="area_weighted_mean",
                                               dim=["lat", "lon"]).compute())
    a_h = ar_h.isel(time=0) if "time" in ar_h.dims else ar_h
    area_hist[m] = a_h
    ocean_area[m] = float(a_h.where(clim_hist[m].notnull()).sum().compute())

gm_woa_1971_2000 = float(compute_basic_stat(sss_woa, op="area_weighted_mean",
                                            dim=["lat", "lon"]).compute())

print("\n--- Global mean SSS 1971-2000 (psu) ---")
print("WOA23 : %.2f   ocean area %.3e m2" % (gm_woa_1971_2000, ocean_area_woa))
for m in MODELS:
    print("%-14s : %.2f   minus WOA23 %+0.2f   ocean area %.3e m2  (area_source=%s, volume_source=%s)"
          % (m, gm_1971_2000[m], gm_1971_2000[m] - gm_woa_1971_2000, ocean_area[m],
             area_src[m], vol_src[m]))

# ------------------------------------------------- comparison on the WOA23 grid
reg_hist = {}
for m in MODELS:
    reg_hist[m] = regrid_horizontal(clim_hist[m], target="woa23",
                                    method="conservative").compute()

cmp_mask = build_comparison_mask([reg_hist[m] for m in MODELS] + [sss_woa],
                                 ocean_only=True).compute()
n_cmp = int(cmp_mask.sum().compute())
print("\ncomparison cells:", n_cmp)

diff_maps, bias_woa, rmse_woa, corr_woa = {}, {}, {}, {}
mean_data_woa, mean_ref_woa, mae_woa = {}, {}, {}

for m in MODELS:
    a = reg_hist[m].where(cmp_mask)
    b = sss_woa.where(cmp_mask)
    dm = compute_difference_metrics(a, b)
    bias_woa[m] = float(dm["bias"])
    rmse_woa[m] = float(dm["rmse"])
    mae_woa[m] = float(dm["mae"])
    mean_data_woa[m] = float(dm["mean_data"])
    mean_ref_woa[m] = float(dm["mean_reference"])
    diff_maps[m] = dm["difference"].assign_attrs(units="psu",
                                                 long_name="SSS bias model - WOA23")
    cm = analyze_correlation(a, b, dim=["lat", "lon"], method="pearson", weighting="area")
    corr_woa[m] = float(cm["correlation"])

print("\n--- Bias / RMSE / pattern correlation vs WOA23 (1971-2000) ---")
for m in MODELS:
    print("%-14s bias %+0.2f psu   RMSE %.2f psu   r = %.3f   mean(model) %.2f  mean(WOA) %.2f"
          % (m, bias_woa[m], rmse_woa[m], corr_woa[m], mean_data_woa[m], mean_ref_woa[m]))

# ------------------------------------------------------------ historical trends
trend_hist, trend_ssp, gm_at_2014 = {}, {}, {}
for m in MODELS:
    trend_hist[m] = fit_1D_trend(gm_hist_year[m], method="linear")
    trend_ssp[m] = fit_1D_trend(gm_ssp_year[m], method="linear")

print("\n--- Linear trend of annual global-mean SSS (psu per decade, 95% CI) ---")
for m in MODELS:
    th, ts = trend_hist[m]["trend"], trend_ssp[m]["trend"]
    rh, rs = trend_hist[m]["metrics"]["r2"], trend_ssp[m]["metrics"]["r2"]
    print("%-14s HIST 1950-2014 %+.3e [%+.3e, %+.3e]  R2=%.3f  significant=%s"
          % (m, th["slope"] * 10, th["ci95_low"] * 10, th["ci95_high"] * 10, rh, th["significant"]))
    print("%-14s SSP  2071-2100 %+.3e [%+.3e, %+.3e]  R2=%.3f  significant=%s"
          % (m, ts["slope"] * 10, ts["ci95_low"] * 10, ts["ci95_high"] * 10, rs, ts["significant"]))
    print("%-14s ratio ssp/hist = %.1f" % (m, ts["slope"] / th["slope"]))

# ---------------------------------------------------------- end-of-century ΔSSS
delta_global, delta_basin = {}, {}
for m in MODELS:
    delta_global[m] = float(compute_basic_stat(delta[m], op="area_weighted_mean",
                                               dim=["lat", "lon"]).compute())
    for b in BASINS:
        db = select_basin(delta[m], basin=b, include_marginal_seas=True,
                          depth_dependent=False)
        delta_basin[(m, b)] = float(compute_basic_stat(db, op="area_weighted_mean",
                                                       dim=["lat", "lon"]).compute())

print("\n--- ΔSSS (ssp245 2071-2100 minus historical 1995-2014), psu ---")
for m in MODELS:
    print("%-14s global %+0.2f | Atlantic %+0.2f | Indian %+0.2f | Pacific %+0.2f"
          % (m, delta_global[m], delta_basin[(m, "Atlantic")],
             delta_basin[(m, "Indian")], delta_basin[(m, "Pacific")]))

# ------------------------------------------------- cell-by-cell historical trend
trend_map_slope, trend_map_sig = {}, {}
for m in MODELS:
    tm = fit_trend_map(ann_hist[m], dim="year", per_decade=True)
    slope_r = regrid_horizontal(tm["slope"], target="woa23", method="bilinear").compute()
    sig = tm["significant"].astype("float32").assign_attrs(**tm["significant"].attrs)
    sig_r = regrid_horizontal(sig, target="woa23", method="bilinear").compute()
    trend_map_slope[m] = slope_r.assign_attrs(units="psu per decade",
                                              long_name="SSS trend 1950-2014")
    trend_map_sig[m] = sig_r
    print("%-14s trend-map slope range %.3e .. %.3e psu/decade"
          % (m, float(slope_r.min().compute()), float(slope_r.max().compute())))

print("\nWOA23 basin names:", basin_mask.attrs.get("basin_names"))
print("land-sea mask: ocean where value > 1 (value is an index, not 0/1)")

result = {
    "global_mean_1971_2000": {m: round(gm_1971_2000[m], 2) for m in MODELS},
    "global_mean_1971_2000_woa": round(gm_woa_1971_2000, 2),
    "bias_vs_woa": {m: round(gm_1971_2000[m] - gm_woa_1971_2000, 2) for m in MODELS},
    "rmse": {m: round(rmse_woa[m], 2) for m in MODELS},
    "correlation": {m: round(corr_woa[m], 3) for m in MODELS},
    "delta_global": {m: round(delta_global[m], 2) for m in MODELS},
    "delta_basin": {f"{m}:{b}": round(delta_basin[(m, b)], 2) for m in MODELS for b in BASINS},
}
print(result)
