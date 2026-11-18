models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
start, end = "1971-01", "2000-12"

def scalar(x):
    if hasattr(x, "compute"):
        return float(x.compute())
    return float(x)

# --- WOA23 reference at 0 m, 1971-2000 annual climatology ---
woa = load_climatology_data(
    dataset="woa23",
    variable="temperature",
    time_slice="annual",
    period="1971-2000",
    depth_range=[0, 0],
)
if "depth" in woa.dims:
    if woa.sizes["depth"] == 1:
        woa2d = woa.isel(depth=0, drop=True)
    else:
        woa2d = woa.sel(depth=0.0, method="nearest")
else:
    woa2d = woa

landsea = load_landsea_mask()
ocean = landsea > 1
woa_ocean = woa2d.where(ocean)

# WOA23 native-grid area-weighted global mean over valid ocean cells
woa_native_mean = scalar(
    compute_basic_stat(data=woa_ocean, op="area_weighted_mean", dim=None)
)
woa_area = load_areas(dataset="woa23")
woa_native_area = scalar(woa_area.where(woa_ocean.notnull()).sum())

# --- CMIP6 models on their native grids ---
native_means = {}
native_areas = {}
model_clims = {}

for model in models:
    tos = load_regular_data(
        dataset="cmip6_historical",
        variable="tos",
        model=model,
        time_range=[start, end],
    )
    units = tos.attrs.get("units", "K")
    if "k" in units.lower() and "degc" not in units.lower():
        tos_c = transform_unit(data=tos, from_unit=units, to_unit="degC")
    else:
        tos_c = tos

    # Day-weighted 1971-2000 long-term annual-mean SST climatology
    clim = temporal_aggregate(data=tos_c, freq="mean", op="mean")
    model_clims[model] = clim

    # Native-grid area-weighted global mean
    native_means[model] = scalar(
        compute_basic_stat(data=clim, op="area_weighted_mean", dim=None)
    )

    area = load_regular_data(
        dataset="cmip6_historical",
        variable="areacello",
        model=model,
    )
    native_areas[model] = scalar(area.where(clim.notnull()).sum())

# --- Regrid each model climatology to the WOA23 1-degree grid ---
regridded = {}
for model in models:
    regridded[model] = regrid_horizontal(
        data=model_clims[model],
        target="woa23",
        method="conservative",
    )

# --- Common comparison mask over all three models and WOA23 ---
mask = build_comparison_mask(
    fields=[regridded[m] for m in models] + [woa2d],
    ocean_only=True,
)
mask_ncells = int(mask.sum().compute())

woa_cmp = woa2d.where(mask)
regridded_masked = {m: regridded[m].where(mask) for m in models}
area_cmp = load_areas(dataset="woa23").where(mask)

# --- Bias, RMSE, shared-cell means, and spatial correlation ---
metrics = {}
for model in models:
    res = compute_difference_metrics(
        data=regridded_masked[model],
        reference=woa_cmp,
        dim=None,
        weighting="auto",
    )
    corr_res = analyze_correlation(
        var1=regridded_masked[model],
        var2=woa_cmp,
        dim=["lat", "lon"],
        method="pearson",
        weighting="area",
    )
    bias_val = scalar(res["bias"])
    metrics[model] = {
        "mean_data": scalar(res["mean_data"]),
        "mean_reference": scalar(res["mean_reference"]),
        "bias": bias_val,
        "delta_SST": bias_val,
        "rmse": scalar(res["rmse"]),
        "mae": scalar(res["mae"]),
        "correlation": scalar(corr_res["correlation"]),
    }

ranking = sorted(models, key=lambda m: abs(metrics[m]["delta_SST"]))

bias_fields = {}
for model in models:
    bias_fields[model] = (regridded_masked[model] - woa_cmp).assign_attrs(
        units="degC",
        long_name=f"{model} SST bias",
    )

print("WOA23 native-grid area-weighted mean SST (degC):", round(woa_native_mean, 2))
print("WOA23 native area covered (m2):", woa_native_area)
for model in models:
    print(f"{model} native-grid area-weighted mean SST (degC):", round(native_means[model], 2))
    print(f"{model} native area covered (m2):", native_areas[model])
print("Common comparison mask cells:", mask_ncells)
for model in models:
    print(model, "metrics:", {k: round(v, 3) if isinstance(v, float) else v for k, v in metrics[model].items()})
print("Ranking (smallest to largest |delta_SST|):", ranking)
