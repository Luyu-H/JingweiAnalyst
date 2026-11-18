import numpy as np
import xarray as xr

DEPTH_RANGE = [200.0, 1000.0]
RHO = 1025.0                      # kg m-3, prescribed by the question
SOURCES = ["WOA23", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
MODELS = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

# ----------------------------------------------------------------------------
# small helpers
# ----------------------------------------------------------------------------
def hdims(da):
    """Horizontal dimension names actually present on this array."""
    return [d for d in da.dims if d != "depth"]

def patch_attrs(derived, template):
    """Carry the loader provenance onto a derived field (keeping new units)."""
    extra = {k: v for k, v in template.attrs.items()
             if k not in derived.attrs and k != "units"}
    derived.attrs.update(extra)
    return derived

def np_ratio(no3, po4):
    r = compute_ratio(no3, po4)
    r = patch_attrs(r, no3)
    r.attrs["units"] = "1"
    return r

def vw_mean(da, dims):
    r = compute_basic_stat(da, op="volume_weighted_mean", dim=dims)
    return float(np.asarray(r.compute()))

def frac_percent(da, threshold, comparison, reference):
    out = compute_exceedance_fraction(da, threshold=threshold, comparison=comparison,
                                      dim=None, weighting="volume", reference=reference)
    return float(out["fraction_percent"])

def flatten_pair(x, y):
    xv = np.asarray(x.compute(), dtype=float).ravel()
    yv = np.asarray(y.compute(), dtype=float).ravel()
    ok = np.isfinite(xv) & np.isfinite(yv)
    return xv[ok], yv[ok]

# ----------------------------------------------------------------------------
# 1. OBSERVATIONS - WOA23 nitrate / phosphate, annual analysed climatology
# ----------------------------------------------------------------------------
woa_no3 = load_climatology_data(dataset="woa23", variable="nitrate",
                                time_slice="annual", depth_range=DEPTH_RANGE)
woa_po4 = load_climatology_data(dataset="woa23", variable="phosphate",
                                time_slice="annual", depth_range=DEPTH_RANGE)
woa_levels, woa_bounds = load_depth_levels(dataset="woa23", variable="nitrate",
                                           depth_range=DEPTH_RANGE)
woa_areas = load_areas(dataset="woa23")

woa_dz = xr.DataArray(np.diff(np.asarray(woa_bounds), axis=1).ravel(),
                      dims="depth", coords={"depth": np.asarray(woa_levels)})
woa_layer_volume = float((woa_areas * woa_dz).where(woa_no3.notnull()).sum().compute())

# ----------------------------------------------------------------------------
# 2. MODELS - CMIP6 historical no3 / po4, 1971-2000, each model's own volcello
# ----------------------------------------------------------------------------
FIELDS = {"WOA23": dict(no3=woa_no3, po4=woa_po4, levels=woa_levels,
                        volume=woa_layer_volume, var_source="WOA23 annual climatology")}

for m in MODELS:
    no3 = load_regular_data(dataset="cmip6_historical", variable="no3", model=m,
                            time_range=["1971-01", "2000-12"], depth_range=DEPTH_RANGE)
    po4 = load_regular_data(dataset="cmip6_historical", variable="po4", model=m,
                            time_range=["1971-01", "2000-12"], depth_range=DEPTH_RANGE)

    # 3. day-weighted 1971-2000 climatology (one annual-mean field per model)
    no3c = temporal_aggregate(no3, freq="mean", op="mean")
    po4c = temporal_aggregate(po4, freq="mean", op="mean")

    vol = load_volumes_for_cmip6(dataset="cmip6_historical", model=m,
                                 depth_range=DEPTH_RANGE, variable="no3")
    if "time" in vol.dims:
        vol = vol.isel(time=0, drop=True)
    vol_total = float(vol.sum().compute())

    lev, bnds = load_depth_levels(dataset="cmip6_historical", model=m,
                                  variable="no3", depth_range=DEPTH_RANGE)

    FIELDS[m] = dict(no3=no3c, po4=po4c, levels=lev, bounds=bnds,
                     volume=vol_total,
                     area_source=vol.attrs.get("area_source"),
                     volume_source=vol.attrs.get("volume_source"))

# ----------------------------------------------------------------------------
# 4-6. cell-by-cell N:P on each source's OWN grid, own levels, own volumes
# ----------------------------------------------------------------------------
for name in SOURCES:
    d = FIELDS[name]
    no3, po4 = d["no3"], d["po4"]
    ratio = np_ratio(no3, po4)
    d["ratio"] = ratio

    dims_all = ["depth"] + hdims(no3)

    # (a) headline: volume-weighted mean of the cell-by-cell ratio
    d["headline"] = vw_mean(ratio, dims_all)
    # (b) ratio of the volume-weighted means
    d["mean_no3"] = vw_mean(no3, dims_all)
    d["mean_po4"] = vw_mean(po4, dims_all)
    d["ratio_of_means"] = d["mean_no3"] / d["mean_po4"]

    # depth-mean of N:P at each horizontal cell (volumes/thicknesses as weights)
    d["layer"] = compute_basic_stat(ratio, op="volume_weighted_mean",
                                    dim="depth").compute()

    # 9. volume fractions
    d["frac_14_18"] = frac_percent(ratio, threshold=2.0, comparison="abs_less", reference=16.0)
    d["frac_below12"] = frac_percent(ratio, threshold=12.0, comparison="less", reference=0.0)

# ----------------------------------------------------------------------------
# 7. concentrations in umol/kg at rho = 1025 kg m-3 (models only; WOA already mass)
# ----------------------------------------------------------------------------
for name in SOURCES:
    d = FIELDS[name]
    if name == "WOA23":
        d["no3_uk"] = d["no3"]
        d["po4_uk"] = d["po4"]
    else:
        d["no3_uk"] = transform_unit(d["no3"], from_unit="mol m-3",
                                     to_unit="umol/kg", density=1.025)
        d["po4_uk"] = transform_unit(d["po4"], from_unit="mol m-3",
                                     to_unit="umol/kg", density=1.025)

# ----------------------------------------------------------------------------
# 8. OLS regression NO3 = a*PO4 + b over every valid cell, equally weighted
# ----------------------------------------------------------------------------
for name in SOURCES:
    d = FIELDS[name]
    xs, ys = flatten_pair(d["po4_uk"], d["no3_uk"])        # x = PO4, y = NO3
    idx = np.arange(xs.size)
    target = xr.DataArray(ys, dims="sample", coords={"sample": idx},
                          attrs={"units": "umol kg-1"})
    predic = xr.DataArray(xs, dims="sample", coords={"sample": idx},
                          attrs={"units": "umol kg-1"})
    fit = fit_multivariate_relation(target=target, predictors=[predic],
                                    method="linear", validation="random")

    # OLS coefficients (identical to the tool's full-sample fit)
    a_ols, b_ols = np.polyfit(xs, ys, 1)
    yhat = a_ols * xs + b_ols
    ss_res = float(np.sum((ys - yhat) ** 2))
    ss_tot = float(np.sum((ys - ys.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot

    d["fit"] = fit
    d["slope"] = float(a_ols)
    d["intercept"] = float(b_ols)
    d["r2"] = float(r2)
    d["n_cells"] = int(xs.size)

    # store a plot subsample (display only)
    sel = np.random.default_rng(0).choice(xs.size, size=min(12000, xs.size), replace=False)
    d["scatter_x"] = xs[sel]
    d["scatter_y"] = ys[sel]

# ----------------------------------------------------------------------------
# 10. basin means of the headline N:P (WOA23 basinmask_01, surface labels)
# ----------------------------------------------------------------------------
basin_mask = load_basin_mask(depth_levels="0m")
BASINS = ["Pacific", "Atlantic", "Indian"]

for name in SOURCES:
    d = FIELDS[name]
    dims_all = ["depth"] + hdims(d["no3"])
    d["basin"] = {}
    for b in BASINS:
        sel = select_basin(d["ratio"], basin=b, include_marginal_seas=True,
                           depth_dependent=False)
        d["basin"][b] = vw_mean(sel, dims_all)

    # 9b. volume fraction below N:P = 12 within Pacific and Indian
    d["frac_below12_basin"] = {}
    for b in ["Pacific", "Indian"]:
        sel = select_basin(d["ratio"], basin=b, include_marginal_seas=True,
                           depth_dependent=False)
        d["frac_below12_basin"][b] = frac_percent(sel, threshold=12.0,
                                                  comparison="less", reference=0.0)

# ----------------------------------------------------------------------------
# 11. like-for-like comparison mask: models regridded to the WOA23 grid
#     (display / diagnosis only - per-source numbers above stay native)
# ----------------------------------------------------------------------------
DISPLAY_MAPS = {}
for name in SOURCES:
    layer = FIELDS[name]["layer"]
    if name == "WOA23":
        DISPLAY_MAPS[name] = layer
    else:
        DISPLAY_MAPS[name] = regrid_horizontal(layer, target="woa23",
                                               method="conservative").compute()

comp_mask = build_comparison_mask(
    [DISPLAY_MAPS[n].where(np.isfinite(DISPLAY_MAPS[n])) for n in SOURCES],
    ocean_only=True)
n_common = int(comp_mask.sum().compute())

# ----------------------------------------------------------------------------
# results table
# ----------------------------------------------------------------------------
rows = []
for name in SOURCES:
    d = FIELDS[name]
    rows.append(dict(source=name,
                     depth_levels=int(d["no3"].sizes["depth"]),
                     layer_volume_m3=d["volume"],
                     np_mean_of_ratio=round(d["headline"], 2),
                     np_ratio_of_means=round(d["ratio_of_means"], 2),
                     slope=round(d["slope"], 2),
                     intercept_umolkg=round(d["intercept"], 2),
                     r2=round(d["r2"], 3),
                     pct_14_18=round(d["frac_14_18"], 1),
                     pct_below12=round(d["frac_below12"], 1),
                     pacific=round(d["basin"]["Pacific"], 2),
                     atlantic=round(d["basin"]["Atlantic"], 2),
                     indian=round(d["basin"]["Indian"], 2),
                     pct_below12_pacific=round(d["frac_below12_basin"]["Pacific"], 1),
                     pct_below12_indian=round(d["frac_below12_basin"]["Indian"], 1),
                     n_cells=d["n_cells"],
                     redfield_position=("at Redfield (16)" if abs(d["headline"] - 16) < 0.25
                                        else ("above Redfield" if d["headline"] > 16
                                              else "below Redfield"))))
results = pd.DataFrame(rows).set_index("source")

results_print = results.copy()
results_print["layer_volume_m3"] = results_print["layer_volume_m3"].map(lambda v: f"{v:.3e}")
print(results_print.to_string())

print("\nshared comparison mask (models regridded to WOA23 grid):", n_common, "ocean cells")
print("WOA23 200-1000 m layer volume (m3):", f"{woa_layer_volume:.3e}")
for m in MODELS:
    print(f"{m}: volume_source={FIELDS[m]['volume_source']}, "
          f"levels in 200-1000 m={FIELDS[m]['no3'].sizes['depth']}, "
          f"regression n={FIELDS[m]['n_cells']}")
print("WOA23 levels in 200-1000 m:", len(woa_levels))
print("basin mask names available:", basin_mask.attrs.get("basin_names"))
