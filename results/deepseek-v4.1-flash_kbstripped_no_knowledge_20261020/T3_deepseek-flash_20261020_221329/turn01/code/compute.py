# ---------------------------------------------------------------
# 0. small helper to materialise lazy results
# ---------------------------------------------------------------
def _f(x):
    if hasattr(x, "compute"):
        x = x.compute()
    return float(np.asarray(x))

def _mat(x):
    return x.compute() if hasattr(x, "compute") else x

models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

# ---------------------------------------------------------------
# 1. CMIP6 historical monthly SST, 1971-2000, native grids
# ---------------------------------------------------------------
tos = {m: load_regular_data(dataset="cmip6_historical", variable="tos", model=m,
                            time_range=["1971-01", "2000-12"])
       for m in models}

# ---------------------------------------------------------------
# 2. 1971-2000 day-weighted long-term mean SST climatology
# ---------------------------------------------------------------
clim = {}
for m in models:
    c = temporal_aggregate(tos[m], freq="mean", op="mean")
    clim[m] = c.assign_attrs(**tos[m].attrs)   # keep loader provenance for weighting

# ---------------------------------------------------------------
# 3. WOA23 0 m annual reference (1971-2000) + land-sea mask + areas
# ---------------------------------------------------------------
woa = load_climatology_data(dataset="woa23", variable="temperature",
                            time_slice="annual", period="1971-2000",
                            depth_range=[0.0, 0.0])
woa0 = woa.isel(depth=0, drop=True).assign_attrs(**woa.attrs)

landsea = _mat(load_landsea_mask())
ocean = landsea > 1

areas_1deg = load_areas()   # 1-degree cell areas (kept for reference)

# ---------------------------------------------------------------
# 4. Native-grid global means of the model climatologies
#    (each model's own areacello, own ocean cells)
# ---------------------------------------------------------------
native_mean = {}
for m in models:
    native_mean[m] = _f(compute_basic_stat(clim[m], op="area_weighted_mean",
                                           dim=["lat", "lon"]))

# ---------------------------------------------------------------
# 5. Global mean of the WOA23 0 m field over the whole 1-deg grid
#    (all valid cells, no land-sea mask applied)
# ---------------------------------------------------------------
woa_global_mean = _f(compute_basic_stat(woa0, op="area_weighted_mean",
                                        dim=["lat", "lon"]))

# ---------------------------------------------------------------
# 6. Regrid each model climatology conservatively to the WOA23 grid
# ---------------------------------------------------------------
reg = {m: regrid_horizontal(clim[m], target="woa23", method="conservative")
       for m in models}

# ---------------------------------------------------------------
# 7. Common comparison mask: all 3 models + WOA23 + landsea>1
# ---------------------------------------------------------------
common = build_comparison_mask([reg[m] for m in models] + [woa0],
                               ocean_only=True, extra_mask=ocean)

reg_c = {m: reg[m].where(common) for m in models}
woa_c = woa0.where(common)

# ---------------------------------------------------------------
# 8. Comparison statistics on the common cells
# ---------------------------------------------------------------
comp_mean = {}
bias, rmse, corr, mae, crmse = {}, {}, {}, {}, {}
bias_fields = {}
corr_p = {}

for m in models:
    comp_mean[m] = _f(compute_basic_stat(reg_c[m], op="area_weighted_mean",
                                         dim=["lat", "lon"]))
    met = compute_difference_metrics(reg_c[m], woa_c, dim=["lat", "lon"],
                                     weighting="area")
    bias[m]    = _f(met["bias"])
    rmse[m]    = _f(met["rmse"])
    mae[m]     = _f(met["mae"])
    crmse[m]   = _f(met["centered_rmse"])
    bias_fields[m] = _mat(met["difference"])

    cor = analyze_correlation(reg_c[m], woa_c, dim=["lat", "lon"],
                              method="pearson", weighting="area")
    corr[m]   = _f(cor["correlation"])
    corr_p[m] = _f(cor["pvalue"])

comp_mean_woa = _f(compute_basic_stat(woa_c, op="area_weighted_mean",
                                      dim=["lat", "lon"]))

# ranking by smallest |bias|
ranking = sorted(models, key=lambda m: abs(bias[m]))

# ---------------------------------------------------------------
# 9. Tidy summary
# ---------------------------------------------------------------
report = pd.DataFrame({
    "model": models,
    "native_global_mean_degC": [round(native_mean[m], 2) for m in models],
    "common_cell_mean_degC":   [round(comp_mean[m], 2) for m in models],
    "bias_vs_WOA23_degC":      [round(bias[m], 2) for m in models],
    "rmse_degC":               [round(rmse[m], 3) for m in models],
    "mae_degC":                [round(mae[m], 3) for m in models],
    "centered_rmse_degC":      [round(crmse[m], 3) for m in models],
    "corr":                    [round(corr[m], 3) for m in models],
    "corr_pvalue":             [round(corr_p[m], 6) for m in models],
}).set_index("model")

woa_row = pd.DataFrame({
    "native_global_mean_degC": [round(woa_global_mean, 2)],
    "common_cell_mean_degC":   [round(comp_mean_woa, 2)],
}, index=["WOA23"])

print("Native-grid global mean SST (1971-2000):")
for m in models:
    print(f"  {m:15s} {native_mean[m]:8.2f} degC")
print(f"  {'WOA23 0 m':15s} {woa_global_mean:8.2f} degC  (whole 1-deg grid, all valid cells)")
print()
print("Common-cell comparison with WOA23:")
print(report.to_string())
print(f"  WOA23 common-cell mean: {comp_mean_woa:.2f} degC")
print()
print("Ranking by |bias| (smallest -> largest):", " < ".join(ranking))

# handy variables for plotting
lon1, lat1 = woa0["lon"].values, woa0["lat"].values
woa_ref_field = _mat(woa_c)
