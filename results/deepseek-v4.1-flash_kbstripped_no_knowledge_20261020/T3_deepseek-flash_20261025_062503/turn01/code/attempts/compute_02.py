models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

# ---------------------------------------------------------------
# 1-3. CMIP6 historical tos, 1971-2000 day-weighted climatology,
#      native-grid area-weighted global mean (native areacello)
# ---------------------------------------------------------------
clim_native = {}
native_global_mean = {}

for m in models:
    tos = load_regular_data(
        dataset="cmip6_historical",
        variable="tos",
        model=m,
        time_range=["1971-01", "2000-12"],
    )
    clim = temporal_aggregate(data=tos, freq="mean", op="mean")   # day-weighted long-term mean
    clim_native[m] = clim
    gm = compute_basic_stat(data=clim, op="area_weighted_mean", dim=["lat", "lon"])
    native_global_mean[m] = float(gm.squeeze().compute())

# ---------------------------------------------------------------
# 4-5. WOA23 0 m annual 1971-2000 temperature + global mean
# ---------------------------------------------------------------
woa = load_climatology_data(
    dataset="woa23",
    variable="temperature",
    time_slice="annual",
    period="1971-2000",
    depth_range=[0.0, 0.0],
)
woa_sfc = woa.squeeze("depth", drop=True)          # (lat, lon) on the 1-degree grid

woa_global_mean = float(
    compute_basic_stat(data=woa, op="area_weighted_mean", dim=["lat", "lon"])
    .squeeze()
    .compute()
)

landsea = load_landsea_mask()                       # ocean where value > 1

# ---------------------------------------------------------------
# 6. Regrid each model climatology onto the WOA23 1-degree grid
# ---------------------------------------------------------------
reg = {}
for m in models:
    reg[m] = regrid_horizontal(data=clim_native[m], target="woa23", method="conservative")

# ---------------------------------------------------------------
# 7. Common comparison cells
# ---------------------------------------------------------------
cmp_mask = build_comparison_mask(
    fields=[reg[m] for m in models] + [woa_sfc],
    ocean_only=True,
)

woa_cmp = woa_sfc.where(cmp_mask)

# ---------------------------------------------------------------
# 8. Mean SST on the comparison cells
# ---------------------------------------------------------------
woa_cmp_mean = float(
    compute_basic_stat(data=woa_cmp, op="area_weighted_mean", dim=["lat", "lon"])
    .squeeze()
    .compute()
)

model_cmp_mean = {}
for m in models:
    model_cmp_mean[m] = float(
        compute_basic_stat(data=reg[m].where(cmp_mask), op="area_weighted_mean", dim=["lat", "lon"])
        .squeeze()
        .compute()
    )

# ---------------------------------------------------------------
# 9. Bias, RMSE, correlation on the comparison cells
# ---------------------------------------------------------------
bias_mean = {}
rmse = {}
corr = {}
diff_field = {}

for m in models:
    mm = reg[m].where(cmp_mask)

    dm = compute_difference_metrics(data=mm, reference=woa_cmp, dim=["lat", "lon"])
    rmse[m] = float(np.asarray(dm["rmse"]).item())

    diff_field[m] = dm["difference"]
    bias_mean[m] = float(
        compute_basic_stat(data=diff_field[m], op="area_weighted_mean", dim=["lat", "lon"])
        .squeeze()
        .compute()
    )

    cr = analyze_correlation(var1=mm, var2=woa_cmp, dim=["lat", "lon"], method="pearson")
    corr[m] = float(np.asarray(cr["correlation"]).item())

# ---------------------------------------------------------------
# 10. Ranking by |delta_SST|
# ---------------------------------------------------------------
ranking = sorted(models, key=lambda m: abs(bias_mean[m]))

summary = pd.DataFrame(
    {
        "native_global_mean_SST_C": [native_global_mean[m] for m in models] + [woa_global_mean],
        "comparison_mean_SST_C": [model_cmp_mean[m] for m in models] + [woa_cmp_mean],
        "bias_dSST_C": [bias_mean[m] for m in models] + [np.nan],
        "RMSE_C": [rmse[m] for m in models] + [np.nan],
        "spatial_corr": [corr[m] for m in models] + [np.nan],
    },
    index=models + ["WOA23"],
)

print(summary.round(3))
print("\nWOA23 global mean SST (native 1-degree grid): %.2f degC" % woa_global_mean)
print("Ocean cells used in comparison:", int(cmp_mask.sum().compute()))
print("Ranking (smallest to largest |delta_SST|):", " < ".join(ranking))

# materialise maps for plotting
bias_maps = {m: diff_field[m].compute() for m in models}
woa_map = woa_cmp.compute()
lon1d = woa_sfc["lon"].values
lat1d = woa_sfc["lat"].values
