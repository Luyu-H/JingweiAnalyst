# ---------------------------------------------------------------
# 1. CMIP6 historical SST: native-grid 1971-2000 climatology
# ---------------------------------------------------------------
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
time_range = ["1971-01", "2000-12"]

model_clim = {}
native_mean = {}

for m in models:
    tos = load_regular_data(dataset="cmip6_historical", variable="tos",
                            model=m, time_range=time_range)
    tos = transform_unit(tos, from_unit=tos.attrs.get("units", "K"), to_unit="degC")
    clim = temporal_aggregate(tos, freq="mean")          # day-weighted long-term mean
    model_clim[m] = clim
    # native-ocean global mean, weighted by the model's own areacello
    native_mean[m] = float(
        compute_basic_stat(clim, op="area_weighted_mean", dim=["lat", "lon"]).compute()
    )

# ---------------------------------------------------------------
# 2. WOA23 0 m annual SST climatology, 1971-2000 normal
# ---------------------------------------------------------------
woa = load_climatology_data(dataset="woa23", variable="temperature",
                            period="1971-2000", time_slice="annual",
                            depth_range=[0, 10])
woa_sst = woa.isel(depth=0, drop=True) if "depth" in woa.dims else woa
woa_sst = woa_sst.assign_attrs(units="degC")

woa_native_mean = float(
    compute_basic_stat(woa_sst, op="area_weighted_mean", dim=["lat", "lon"]).compute()
)

# WOA23 1-degree land-sea mask: ocean is (value > 1)
lsm = load_landsea_mask()
ocean = lsm > 1

# ---------------------------------------------------------------
# 3. Regrid each model climatology conservatively to the WOA23 grid
# ---------------------------------------------------------------
regridded = {}
for m in models:
    regridded[m] = regrid_horizontal(model_clim[m], target="woa23", method="conservative")

# ---------------------------------------------------------------
# 4. Common comparison mask (WOA ocean cells valid for every field)
# ---------------------------------------------------------------
mask = build_comparison_mask([regridded[m] for m in models] + [woa_sst],
                             ocean_only=True, extra_mask=ocean)

# ---------------------------------------------------------------
# 5. Comparison statistics on the common cells
# ---------------------------------------------------------------
comp_mean = {}
delta = {}
rmse = {}
corr = {}
bias_field = {}
woa_comp_mean = None

for m in models:
    f = regridded[m].where(mask)
    w = woa_sst.where(mask)

    dm = compute_difference_metrics(data=f, reference=w, dim=["lat", "lon"])
    comp_mean[m] = float(dm["mean_data"])
    delta[m] = float(dm["bias"])
    rmse[m] = float(dm["rmse"])
    bias_field[m] = dm["difference"].assign_attrs(
        units="degC", long_name=f"{m} SST bias (model - WOA23)"
    )
    if woa_comp_mean is None:
        woa_comp_mean = float(dm["mean_reference"])

    c = analyze_correlation(var1=f, var2=w, dim=["lat", "lon"],
                            method="pearson", weighting="area")
    corr[m] = float(c["correlation"])

# ---------------------------------------------------------------
# 6. Ranking by absolute mean bias
# ---------------------------------------------------------------
ranking = sorted(models, key=lambda mm: abs(delta[mm]))

print("Native-grid global mean SST (degC):")
for m in models:
    print(f"  {m:15s} {native_mean[m]:.2f}")
print(f"  {'WOA23 (full grid)':15s} {woa_native_mean:.2f}")

print("\nOn the common WOA23 comparison cells:")
print(f"  WOA23 mean SST: {woa_comp_mean:.2f} degC")
for m in models:
    print(f"  {m:15s} mean={comp_mean[m]:.2f}  delta_SST={delta[m]:+.2f}  "
          f"RMSE={rmse[m]:.3f}  r={corr[m]:.3f}")

print("\nRanking (smallest -> largest |delta_SST|):")
for i, m in enumerate(ranking, 1):
    print(f"  {i}. {m}  (|delta| = {abs(delta[m]):.2f} degC)")

summary = pd.DataFrame({
    "native_global_mean_degC": [native_mean[m] for m in models] + [woa_native_mean],
    "comparison_mean_degC":    [comp_mean[m] for m in models] + [woa_comp_mean],
    "delta_SST_degC":          [delta[m] for m in models] + [np.nan],
    "RMSE_degC":               [rmse[m] for m in models] + [np.nan],
    "spatial_r":               [corr[m] for m in models] + [np.nan],
}, index=models + ["WOA23"])
