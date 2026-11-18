import numpy as np
import xarray as xr

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _depth_dim(da):
    for d in ("depth", "lev", "olevel", "z", "st_ocean"):
        if d in da.dims:
            return d
    raise ValueError(f"No depth dimension found in {da.dims}")

def _preserve_attrs(out, src, units=None, long_name=None):
    out = out.assign_attrs(**src.attrs)
    if units is not None:
        out.attrs["units"] = units
    if long_name is not None:
        out.attrs["long_name"] = long_name
    return out

def _sss_from_levels(sal, vol, levels, depth_dim=None, units="psu"):
    """Full-cell-volume weighted mean salinity over levels in 0-10 m."""
    if depth_dim is None:
        depth_dim = _depth_dim(sal)
    vd = _depth_dim(vol)
    if vd != depth_dim:
        vol = vol.rename({vd: depth_dim})

    lev = np.asarray(levels, dtype=float)
    lev_sel = lev[(lev >= 0.0) & (lev <= 10.0)]
    if lev_sel.size == 0:
        raise ValueError("No levels in 0-10 m")

    sal_sel = sal.sel({depth_dim: lev_sel})
    vol_sel = vol.sel({depth_dim: lev_sel})

    # Keep only positive volume and valid salinity; use full cell volume.
    w = vol_sel.where(np.isfinite(sal_sel) & (vol_sel > 0))
    num = (sal_sel * w).sum(dim=depth_dim, skipna=True)
    den = w.sum(dim=depth_dim, skipna=True)
    sss = (num / den).where(den > 0)

    sss = _preserve_attrs(
        sss,
        sal,
        units=units,
        long_name="sea surface salinity, 0-10 m full-cell-volume weighted",
    )
    return sss

# ---------------------------------------------------------------------------
# WOA23 reference SSS, 1971-2000 annual
# ---------------------------------------------------------------------------
sal_woa = load_climatology_data(
    dataset="woa23",
    variable="salinity",
    time_slice="annual",
    period="1971-2000",
    depth_range=[0, 10],
)
lev_woa, bnds_woa = load_depth_levels(
    dataset="woa23",
    variable="salinity",
    time_slice="annual",
    period="1971-2000",
    depth_range=[0, 10],
)
area_woa = load_areas(dataset="woa23")
landsea = load_landsea_mask()
basinmask = load_basin_mask()

dd_woa = _depth_dim(sal_woa)
lev_woa_axis = sal_woa[dd_woa].values.astype(float)

# Full layer thickness from the layer boundaries, not the overlap with 0-10 m.
thick_woa = np.abs(bnds_woa[:, 1] - bnds_woa[:, 0])
thick_woa_da = xr.DataArray(
    thick_woa,
    coords={dd_woa: lev_woa_axis},
    dims=(dd_woa,),
)
vol_woa = area_woa * thick_woa_da

sss_woa = _sss_from_levels(sal_woa, vol_woa, lev_woa_axis, depth_dim=dd_woa)
sss_woa = sss_woa.where(landsea > 1)
sss_woa = _preserve_attrs(
    sss_woa,
    sal_woa,
    units="psu",
    long_name="WOA23 1971-2000 SSS, 0-10 m full-cell-volume weighted",
)

woa_global_1971_2000 = float(
    compute_basic_stat(sss_woa, op="area_weighted_mean", dim=["lat", "lon"]).compute()
)

# ---------------------------------------------------------------------------
# CMIP6 historical and ssp245 SSS, 0-10 m volume weighted on native grids
# ---------------------------------------------------------------------------
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

hist_sss = {}
ssp_sss = {}
area_model = {}
vol_hist_model = {}
vol_ssp_model = {}
levels_hist_model = {}
levels_ssp_model = {}

for model in models:
    so_h = load_regular_data(
        dataset="cmip6_historical",
        variable="so",
        model=model,
        time_range=["1950-01", "2014-12"],
    )
    so_s = load_regular_data(
        dataset="cmip6_ssp245",
        variable="so",
        model=model,
        time_range=["2071-01", "2100-12"],
    )

    area_h = load_areas(
        dataset="cmip6_historical",
        model=model,
        variable="so",
    )
    vol_h = load_volumes_for_cmip6(
        dataset="cmip6_historical",
        model=model,
        variable="so",
    )
    vol_s = load_volumes_for_cmip6(
        dataset="cmip6_ssp245",
        model=model,
        variable="so",
    )

    lev_h, bnds_h = load_depth_levels(
        dataset="cmip6_historical",
        model=model,
        variable="so",
    )
    lev_s, bnds_s = load_depth_levels(
        dataset="cmip6_ssp245",
        model=model,
        variable="so",
    )

    dd_h = _depth_dim(so_h)
    dd_s = _depth_dim(so_s)

    sss_h = _sss_from_levels(so_h, vol_h, lev_h, depth_dim=dd_h)
    sss_s = _sss_from_levels(so_s, vol_s, lev_s, depth_dim=dd_s)

    hist_sss[model] = sss_h
    ssp_sss[model] = sss_s
    area_model[model] = area_h
    vol_hist_model[model] = vol_h
    vol_ssp_model[model] = vol_s
    levels_hist_model[model] = lev_h
    levels_ssp_model[model] = lev_s

# ---------------------------------------------------------------------------
# Model 1971-2000 SSS climatology and global means
# ---------------------------------------------------------------------------
hist_clim = {}
model_global_1971_2000 = {}

for model in models:
    sub = hist_sss[model].sel(time=slice("1971-01", "2000-12"))
    clim = temporal_aggregate(sub, freq="mean", op="mean")
    clim = _preserve_attrs(clim, hist_sss[model], units="psu")
    hist_clim[model] = clim

    gm = compute_basic_stat(
        clim, op="area_weighted_mean", dim=["lat", "lon"]
    ).compute()
    model_global_1971_2000[model] = float(gm)

# ---------------------------------------------------------------------------
# Put model climatologies on WOA23 grid, common mask, bias/RMSE/correlation
# ---------------------------------------------------------------------------
regridded_clim = {}
for model in models:
    regridded_clim[model] = regrid_horizontal(
        hist_clim[model],
        target="woa23",
        method="conservative",
    )

fields_for_mask = [regridded_clim[m] for m in models] + [sss_woa]
comp_mask = build_comparison_mask(fields=fields_for_mask, ocean_only=True)
mask_n = int(comp_mask.sum().compute())

woa_masked = sss_woa.where(comp_mask).assign_attrs(**sss_woa.attrs)

bias_metrics = {}
bias_map = {}
corr_coeff = {}

for model in models:
    m_masked = regridded_clim[model].where(comp_mask)
    m_masked = m_masked.assign_attrs(**regridded_clim[model].attrs)

    met = compute_difference_metrics(
        data=m_masked,
        reference=woa_masked,
        dim=None,
        weighting="auto",
    )

    bias_val = met["bias"]
    if hasattr(bias_val, "compute"):
        bias_val = bias_val.compute()
    rmse_val = met["rmse"]
    if hasattr(rmse_val, "compute"):
        rmse_val = rmse_val.compute()

    nv = met["n_valid"]
    if hasattr(nv, "compute"):
        nv = int(nv.compute())
    else:
        nv = int(nv)

    bias_metrics[model] = {
        "bias": float(bias_val),
        "rmse": float(rmse_val),
        "n_valid": nv,
        "weighting": str(met["weighting"]),
    }
    bias_map[model] = met["difference"].compute()

    cor = analyze_correlation(
        var1=m_masked,
        var2=woa_masked,
        dim=["lat", "lon"],
        method="pearson",
        weighting="area",
    )
    cv = cor["correlation"]
    if hasattr(cv, "compute"):
        cv = cv.compute()
    corr_coeff[model] = float(cv)

# ---------------------------------------------------------------------------
# Historical and ssp245 trends of the annual global mean SSS
# ---------------------------------------------------------------------------
hist_global_annual = {}
ssp_global_annual = {}
trend_hist = {}
trend_ssp = {}

for model in models:
    annual_h = temporal_aggregate(hist_sss[model], freq="year", op="mean")
    annual_h = _preserve_attrs(annual_h, hist_sss[model], units="psu")
    gm_h = compute_basic_stat(
        annual_h, op="area_weighted_mean", dim=["lat", "lon"]
    ).compute()
    hist_global_annual[model] = gm_h
    trend_hist[model] = fit_1D_trend(gm_h, method="linear")

    annual_s = temporal_aggregate(ssp_sss[model], freq="year", op="mean")
    annual_s = _preserve_attrs(annual_s, ssp_sss[model], units="psu")
    gm_s = compute_basic_stat(
        annual_s, op="area_weighted_mean", dim=["lat", "lon"]
    ).compute()
    ssp_global_annual[model] = gm_s
    trend_ssp[model] = fit_1D_trend(gm_s, method="linear")

# ---------------------------------------------------------------------------
# Cell-by-cell historical SSS trend maps, and display versions on WOA23 grid
# ---------------------------------------------------------------------------
trend_map_hist = {}
slope_regrid_hist = {}
sig_regrid_hist = {}

for model in models:
    annual_h = temporal_aggregate(hist_sss[model], freq="year", op="mean")
    annual_h = _preserve_attrs(annual_h, hist_sss[model], units="psu")

    tm = fit_trend_map(annual_h, dim="year", alpha=0.05, per_decade=True)
    trend_map_hist[model] = tm

    slope = tm["slope"]
    sig = tm["significant"]

    slope_regrid_hist[model] = regrid_horizontal(
        slope,
        target="woa23",
        method="conservative",
    ).compute()

    sig_float = sig.astype(float).assign_attrs(**slope.attrs)
    sig_regrid = regrid_horizontal(
        sig_float,
        target="woa23",
        method="nearest_s2d",
    ).compute()
    sig_regrid_hist[model] = sig_regrid > 0.5

# ---------------------------------------------------------------------------
# Delta SSS: ssp245 2071-2100 minus historical 1995-2014, global and basins
# ---------------------------------------------------------------------------
hist_mean_1995_2014 = {}
ssp_mean_2071_2100 = {}
delta_sss = {}
delta_global = {}
delta_basins = {}

for model in models:
    hmean = temporal_aggregate(
        hist_sss[model].sel(time=slice("1995-01", "2014-12")),
        freq="mean",
        op="mean",
    )
    hmean = _preserve_attrs(hmean, hist_sss[model], units="psu")

    smean = temporal_aggregate(ssp_sss[model], freq="mean", op="mean")
    smean = _preserve_attrs(smean, ssp_sss[model], units="psu")

    hist_mean_1995_2014[model] = hmean
    ssp_mean_2071_2100[model] = smean

    delta = (smean - hmean).assign_attrs(**hist_sss[model].attrs)
    delta.attrs["units"] = "psu"
    delta.attrs["long_name"] = (
        "SSS change, ssp245 2071-2100 minus historical 1995-2014"
    )
    delta_sss[model] = delta

    dg = compute_basic_stat(
        delta, op="area_weighted_mean", dim=["lat", "lon"]
    ).compute()
    delta_global[model] = float(dg)

    delta_basins[model] = {}
    for basin in ["Atlantic", "Indian", "Pacific"]:
        db = select_basin(
            delta,
            basin=basin,
            include_marginal_seas=True,
            depth_dependent=False,
        )
        val = compute_basic_stat(
            db, op="area_weighted_mean", dim=["lat", "lon"]
        ).compute()
        delta_basins[model][basin] = float(val)

# ---------------------------------------------------------------------------
# Summary print
# ---------------------------------------------------------------------------
def _trend_vals(tr):
    t = tr["trend"]
    slope = float(t["slope"]) * 10.0          # per year -> per decade
    lo = float(t["ci95_low"]) * 10.0
    hi = float(t["ci95_high"]) * 10.0
    r2 = float(tr["metrics"]["r2"])
    sig = (lo > 0.0) or (hi < 0.0)
    return slope, lo, hi, r2, sig

print("WOA23 1971-2000 global mean SSS: {:.2f} psu".format(woa_global_1971_2000))
print("\nModel global mean SSS, 1971-2000:")
for m in models:
    print("  {:16s} {:6.2f} psu   model - WOA23 = {:+.2f} psu".format(
        m, model_global_1971_2000[m], model_global_1971_2000[m] - woa_global_1971_2000
    ))

print("\nBias/RMSE/correlation against WOA23 on the common WOA23 1-degree mask ({} cells):".format(mask_n))
for m in models:
    print("  {:16s} bias={:+.2f} psu  RMSE={:.2f} psu  r={:.3f}".format(
        m, bias_metrics[m]["bias"], bias_metrics[m]["rmse"], corr_coeff[m]
    ))

print("\nLinear trends of annual global mean SSS:")
for m in models:
    sh, loh, hih, r2h, sigh = _trend_vals(trend_hist[m])
    ss, los, his, r2s, sigs = _trend_vals(trend_ssp[m])
    ratio = ss / sh if sh != 0.0 else np.nan
    print("  {:16s} historical 1950-2014: {:+.3g} psu/dec (95%CI {:+.3g} to {:+.3g}), R2={:.3f}, significant={}".format(
        m, sh, loh, hih, r2h, sigh
    ))
    print("  {:16s} ssp245 2071-2100:     {:+.3g} psu/dec (95%CI {:+.3g} to {:+.3g}), R2={:.3f}, significant={}".format(
        m, ss, los, his, r2s, sigs
    ))
    print("  {:16s} ssp245/historical slope ratio = {:.1f}".format(m, ratio))

print("\nDelta SSS (ssp245 2071-2100 minus historical 1995-2014):")
for m in models:
    print("  {:16s} global={:+.2f} psu; Atlantic={:+.2f}, Indian={:+.2f}, Pacific={:+.2f} psu".format(
        m,
        delta_global[m],
        delta_basins[m]["Atlantic"],
        delta_basins[m]["Indian"],
        delta_basins[m]["Pacific"],
    ))

print("\nCaveat: CESM2 historical is r1i1p1f1 while its ssp245 is r10i1p1f1, so its")
print("ssp245 slope and delta SSS mix forced change with internal-variability differences.")
print("GFDL-ESM4 and MPI-ESM1-2-LR use r1 in both experiments.")

# Keep variables needed for plotting:
# bias_map, slope_regrid_hist, sig_regrid_hist,
# hist_global_annual, ssp_global_annual, trend_hist, trend_ssp, models
