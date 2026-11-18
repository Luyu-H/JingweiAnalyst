# Helper for materialising scalar results
def _scalar(x):
    if hasattr(x, "compute"):
        x = x.compute()
    return float(x)

models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

model_clims = {}
model_means_native = {}
model_areas_native = {}
model_areas_da = {}

# ------------------------------------------------------------------
# 1-3. CMIP6 SST climatologies, native-grid means and native areas
# ------------------------------------------------------------------
for model in models:
    tos = load_regular_data(
        dataset="cmip6_historical",
        variable="tos",
        model=model,
        time_range=["1971-01", "2000-12"],
    )

    # Remove obvious fill values; keep valid temperature values.
    tos = tos.where(np.isfinite(tos))
    tos = tos.where(np.abs(tos) < 1.0e10)

    # CMIP6 tos is in K; convert to degrees Celsius.
    tos = transform_unit(
        tos,
        from_unit=tos.attrs.get("units", "K"),
        to_unit="degC",
    )

    # Day-weighted 1971-2000 mean.
    sst_clim = temporal_aggregate(tos, freq="mean", op="mean")
    sst_clim = sst_clim.where(np.isfinite(sst_clim))
    sst_clim = sst_clim.where((sst_clim > -5.0) & (sst_clim < 45.0))
    model_clims[model] = sst_clim

    # Native-grid cell areas.
    areas = load_areas(
        dataset="cmip6_historical",
        model=model,
        variable="tos",
    )
    model_areas_da[model] = areas

    # Native-grid area-weighted global mean.
    mean_native = compute_basic_stat(
        sst_clim,
        op="area_weighted_mean",
        dim=["lat", "lon"],
    )
    model_means_native[model] = _scalar(mean_native)

    # Ocean area represented by each native mean.
    model_areas_native[model] = float(
        areas.where(np.isfinite(sst_clim)).sum().compute()
    )

# ------------------------------------------------------------------
# 4. WOA23 0 m annual 1971-2000 SST and its native global mean
# ------------------------------------------------------------------
woa_temp = load_climatology_data(
    dataset="woa23",
    variable="temperature",
    time_slice="annual",
    period="1971-2000",
    depth_range=[0, 0],
)

if "depth" in woa_temp.dims:
    woa_sst = woa_temp.isel(depth=0, drop=True)
else:
    woa_sst = woa_temp

woa_mask = load_landsea_mask()
woa_ocean = woa_mask > 1

woa_sst = woa_sst.where(np.isfinite(woa_sst))
woa_sst = woa_sst.where(np.abs(woa_sst) < 1.0e10)
woa_sst = woa_sst.where(woa_ocean & np.isfinite(woa_sst))
woa_sst = woa_sst.where((woa_sst > -5.0) & (woa_sst < 45.0))

woa_areas = load_areas(dataset="woa23", variable="temperature")

woa_mean_native = _scalar(
    compute_basic_stat(
        woa_sst,
        op="area_weighted_mean",
        dim=["lat", "lon"],
    )
)
woa_area_native = float(
    woa_areas.where(np.isfinite(woa_sst)).sum().compute()
)

# ------------------------------------------------------------------
# 5. Regrid each model climatology to the WOA23 1-degree grid
# ------------------------------------------------------------------
model_regrid = {}
for model in models:
    model_regrid[model] = regrid_horizontal(
        model_clims[model],
        target="woa23",
        method="conservative",
    )

# ------------------------------------------------------------------
# 6. One common comparison mask and shared-cell fields
# ------------------------------------------------------------------
fields_for_mask = [woa_sst] + [model_regrid[m] for m in models]
comparison_mask = build_comparison_mask(
    fields=fields_for_mask,
    ocean_only=True,
)

n_compare_cells = int(comparison_mask.sum().compute())

woa_shared = woa_sst.where(comparison_mask)
model_shared = {
    model: model_regrid[model].where(comparison_mask)
    for model in models
}

# ------------------------------------------------------------------
# 7. Shared-cell means, biases, RMSE and spatial correlations
# ------------------------------------------------------------------
woa_mean_shared = _scalar(
    compute_basic_stat(
        woa_shared,
        op="area_weighted_mean",
        dim=["lat", "lon"],
    )
)

model_means_shared = {}
model_metrics = {}
model_corr = {}

for model in models:
    model_means_shared[model] = _scalar(
        compute_basic_stat(
            model_shared[model],
            op="area_weighted_mean",
            dim=["lat", "lon"],
        )
    )

    d = compute_difference_metrics(
        data=model_shared[model],
        reference=woa_shared,
        dim=None,
        weighting="auto",
    )

    model_metrics[model] = {
        "bias": _scalar(d["bias"]),
        "rmse": _scalar(d["rmse"]),
        "centered_rmse": _scalar(d["centered_rmse"]),
        "mean_data": _scalar(d["mean_data"]),
        "mean_reference": _scalar(d["mean_reference"]),
        "n_valid": _scalar(d["n_valid"]),
    }

    c = analyze_correlation(
        var1=model_shared[model],
        var2=woa_shared,
        dim=["lat", "lon"],
        method="pearson",
        weighting="area",
    )
    model_corr[model] = _scalar(c["correlation"])

# ------------------------------------------------------------------
# 8. Ranking from smallest to largest absolute bias
# ------------------------------------------------------------------
ranking = sorted(models, key=lambda m: abs(model_metrics[m]["bias"]))

# Materialise fields needed for plotting.
bias_maps = {
    model: (model_shared[model] - woa_shared)
    .compute()
    .assign_attrs(units="degC")
    for model in models
}
woa_shared_c = woa_shared.compute()

# Collect key results.
results = {
    "native_means_degC": model_means_native,
    "native_woa_mean_degC": woa_mean_native,
    "native_areas_m2": model_areas_native,
    "native_woa_area_m2": woa_area_native,
    "comparison_cells": n_compare_cells,
    "shared_means_degC": model_means_shared,
    "shared_woa_mean_degC": woa_mean_shared,
    "metrics": model_metrics,
    "correlations": model_corr,
    "ranking_smallest_to_largest_abs_bias": ranking,
}

print("Native-grid 1971-2000 global mean SST (degC):")
for model in models:
    print(
        f"  {model}: {model_means_native[model]:.2f} degC, "
        f"area = {model_areas_native[model]:.3e} m2"
    )
print(f"  WOA23: {woa_mean_native:.2f} degC, area = {woa_area_native:.3e} m2")

print(f"\nShared comparison cells: {n_compare_cells}")

print("\nShared-cell mean SST (degC):")
for model in models:
    print(f"  {model}: {model_means_shared[model]:.2f}")
print(f"  WOA23: {woa_mean_shared:.2f}")

print("\nShared-cell metrics against WOA23:")
for model in models:
    mm = model_metrics[model]
    print(
        f"  {model}: delta_SST = {mm['bias']:.2f} degC, "
        f"RMSE = {mm['rmse']:.3f} degC, "
        f"r = {model_corr[model]:.3f}"
    )

print("\nRanking (smallest to largest |delta_SST|):")
print("  " + " < ".join(ranking))
