# ============================================================================
# Nitracline depth: 3 CMIP6 models (1971-2000) vs WOA23 observed climatology
# ============================================================================
import numpy as np
import pandas as pd
import xarray as xr

DEPTH_LOAD_MAX = 1300.0     # keep the level just below 1000 m for the gradient
NITRA_MAX      = 1000.0     # nitracline searched in the top 1000 m
MODELS         = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

# ----------------------------------------------------------------------------
# 1. WOA23 observed nitrate (all-period annual climatology) + vertical grid
# ----------------------------------------------------------------------------
woa_no3 = load_climatology_data(dataset="woa23", variable="nitrate", time_slice="annual")
woa_levels, woa_bounds = load_depth_levels(dataset="woa23", variable="nitrate",
                                           time_slice="annual")
print("WOA23 nitrate levels (m):", np.asarray(woa_levels).size,
      "top/bottom:", float(np.min(woa_levels)), float(np.max(woa_levels)))

woa_nit = compute_nutricline_depth(woa_no3, min_depth=0.0, max_depth=NITRA_MAX).compute()
woa_nit_global = float(compute_basic_stat(woa_nit, op="area_weighted_mean",
                                          dim=["lat", "lon"]).compute())

# ----------------------------------------------------------------------------
# 2.-3. CMIP6 hist. no3 -> 1971-2000 day-weighted mean -> unit conv. -> nitracline
# ----------------------------------------------------------------------------
model_no3_clim   = {}
model_nit_native = {}
model_nit_1deg   = {}
model_nit_global = {}

for m in MODELS:
    raw = load_regular_data(dataset="cmip6_historical", variable="no3", model=m,
                            time_range=["1971-01", "2000-12"],
                            depth_range=[0.0, DEPTH_LOAD_MAX])
    clim = temporal_aggregate(raw, freq="mean", op="mean")          # day-weighted
    clim = transform_unit(clim, from_unit="mol m-3", to_unit="umol/kg")
    model_no3_clim[m] = clim

    nit = compute_nutricline_depth(clim, min_depth=0.0, max_depth=NITRA_MAX).compute()
    model_nit_native[m] = nit

    model_nit_global[m] = float(compute_basic_stat(nit, op="area_weighted_mean",
                                                   dim=["lat", "lon"]).compute())
    # step 5: onto the WOA23 1-degree grid (native global field, so conservative ok)
    model_nit_1deg[m] = regrid_horizontal(nit, target="woa23",
                                          method="conservative").compute()

print("\n=== Global mean nitracline depth (m) ===")
print(f"WOA23        : {woa_nit_global:.0f}")
for m in MODELS:
    print(f"{m:13s}: {model_nit_global[m]:.0f}   (model - WOA23 = "
          f"{model_nit_global[m] - woa_nit_global:+.0f} m)")

# ----------------------------------------------------------------------------
# 4./6. Basin definitions (WOA23 basinmask) and basin means
# ----------------------------------------------------------------------------
basin_codes_da = load_basin_mask(depth_levels="0m")     # surface labels, 2-D
BASINS = ["Pacific", "Atlantic", "Indian"]

basin_mean = {}
for b in BASINS:
    basin_mean[("WOA23", b)] = float(compute_basic_stat(
        select_basin(woa_nit, basin=b, include_marginal_seas=True),
        op="area_weighted_mean", dim=["lat", "lon"]).compute())
    for m in MODELS:
        basin_mean[(m, b)] = float(compute_basic_stat(
            select_basin(model_nit_1deg[m], basin=b, include_marginal_seas=True),
            op="area_weighted_mean", dim=["lat", "lon"]).compute())

print("\n=== Mean nitracline depth by basin (m, WOA23 mask, marginal seas folded in) ===")
basin_rows = []
for b in BASINS:
    row = {"basin": b, "WOA23": basin_mean[("WOA23", b)]}
    for m in MODELS:
        row[m] = basin_mean[(m, b)]
        row[m + "_diff"] = basin_mean[(m, b)] - basin_mean[("WOA23", b)]
    basin_rows.append(row)
basin_table = pd.DataFrame(basin_rows).set_index("basin")
print(basin_table.round(1).to_string())

print("\nBasin ordering (deep -> shallow):")
for name in ["WOA23"] + MODELS:
    order = sorted(BASINS, key=lambda b: -basin_mean[(name, b)])
    print(f"  {name:13s}: " + " > ".join(order))

# ----------------------------------------------------------------------------
# 5. Common WOA23 ocean cell set -> spatial correlation and RMSE
# ----------------------------------------------------------------------------
common_mask = build_comparison_mask(
    fields=[woa_nit] + [model_nit_1deg[m] for m in MODELS], ocean_only=True)
print("\nCommon comparison cells:", int(common_mask.sum().compute()))

corr_val, rmse_val, bias_val = {}, {}, {}
for m in MODELS:
    a = model_nit_1deg[m].where(common_mask)
    b = woa_nit.where(common_mask)
    corr_val[m] = float(analyze_correlation(var1=a, var2=b, dim=["lat", "lon"],
                                            method="pearson")["correlation"].compute())
    met = compute_difference_metrics(data=a, reference=b, dim=["lat", "lon"])
    # compute_difference_metrics returns plain floats for scalar metrics over all
    # reduced dims -> no .compute() here
    rmse_val[m] = float(met["rmse"])
    bias_val[m] = float(met["bias"])

print("\n=== Pattern statistics on the common WOA23 ocean cells ===")
for m in MODELS:
    print(f"{m:13s}: r = {corr_val[m]:.3f}   RMSE = {rmse_val[m]:.1f} m   "
          f"bias = {bias_val[m]:+.1f} m")

best_model = min(MODELS, key=lambda m: rmse_val[m])
print("Best model by RMSE:", best_model,
      f"(r={corr_val[best_model]:.3f}, RMSE={rmse_val[best_model]:.1f} m)")

summary = pd.DataFrame({
    "global_mean_m":   {**{"WOA23": woa_nit_global}, **model_nit_global},
    "diff_vs_WOA23_m": {**{"WOA23": 0.0},
                        **{m: model_nit_global[m] - woa_nit_global for m in MODELS}},
    "r":               {m: corr_val[m] for m in MODELS},
    "RMSE_m":          {m: rmse_val[m] for m in MODELS},
})
print("\n", summary.round(2).to_string())

# ----------------------------------------------------------------------------
# 7./8. Equatorial Pacific section, 120E-280E, 2S-2N (continuous longitude)
# ----------------------------------------------------------------------------
woa_sec = load_climatology_data(dataset="woa23", variable="nitrate", time_slice="annual",
                                lat_range=[-2.0, 2.0], lon_range=[120.0, 280.0],
                                depth_range=[0.0, DEPTH_LOAD_MAX])
sec_lat = np.asarray(woa_sec["lat"].values, dtype=float)
woa_sec_mean = compute_basic_stat(woa_sec, op="mean", dim="lat").compute()
sec_lon = np.asarray(woa_sec_mean["lon"].values, dtype=float)
print("\nSection grid:", sec_lat.size, "latitudes", sec_lat[[0, -1]],
      "|", sec_lon.size, "longitudes", sec_lon[0], "->", sec_lon[-1])

woa_sec_nit = compute_nutricline_depth(woa_sec_mean, min_depth=0.0,
                                       max_depth=NITRA_MAX).compute()

# template on the WOA23 section grid for the curvilinear models.
# Build it from an existing WOA23 field (isel a level) so it keeps the
# 'source_dataset' provenance attribute; a bare xr.DataArray would drop it and
# break compute_nutricline_depth's grid lookup.
template = woa_sec.isel(depth=0, drop=True)

model_sec, model_sec_nit = {}, {}
for m in MODELS:
    raw = load_regular_data(dataset="cmip6_historical", variable="no3", model=m,
                            time_range=["1971-01", "2000-12"],
                            lat_range=[-2.0, 2.0], lon_range=[120.0, 280.0],
                            depth_range=[0.0, DEPTH_LOAD_MAX])
    clim = temporal_aggregate(raw, freq="mean", op="mean")
    clim = transform_unit(clim, from_unit="mol m-3", to_unit="umol/kg")
    onsec = regrid_horizontal(clim, target=template, method="bilinear")
    onsec = regrid_vertical(onsec, target="woa23", method="linear")
    onsec = compute_basic_stat(onsec, op="mean", dim="lat").compute()
    model_sec[m] = onsec
    model_sec_nit[m] = compute_nutricline_depth(onsec, min_depth=0.0,
                                                max_depth=NITRA_MAX).compute()

print("\nSection nitracline: mean over the basin (m)")
print(f"  WOA23        : {float(woa_sec_nit.mean()):.0f}")
for m in MODELS:
    print(f"  {m:13s}: {float(model_sec_nit[m].mean()):.0f}")
