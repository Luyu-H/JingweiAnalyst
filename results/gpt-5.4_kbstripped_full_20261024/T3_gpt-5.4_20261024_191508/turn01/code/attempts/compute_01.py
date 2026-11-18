# Computation code only

models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
time_range = ["1971-01", "2000-12"]

# Containers
tos_monthly = {}
tos_clim_native = {}
native_global_means = {}

# 1) Load CMIP6 monthly tos and form 1971-2000 day-weighted climatological mean on each native grid
for model in models:
    tos_monthly[model] = load_regular_data(
        dataset="cmip6_historical",
        variable="tos",
        model=model,
        time_range=time_range,
    )
    tos_clim_native[model] = temporal_aggregate(
        data=tos_monthly[model],
        freq="mean",
        op="mean",
    )
    native_global_means[model] = float(
        compute_basic_stat(
            data=tos_clim_native[model],
            op="area_weighted_mean",
            dim=["lat", "lon"],
        ).compute()
    )

# 2) Load WOA23 1971-2000 annual temperature climatology at 0 m
woa_surface = load_climatology_data(
    dataset="woa23",
    variable="temperature",
    period="1971-2000",
    time_slice="annual",
    depth_range=[0.0, 0.0],
)

# Drop the singleton depth dimension if present while preserving attrs
if "depth" in woa_surface.dims and woa_surface.sizes["depth"] == 1:
    woa_surface = woa_surface.isel(depth=0, drop=True)

# 3) Load WOA23 land-sea mask and keep ocean where landsea_01 > 1
landsea_01 = load_landsea_mask()
woa_ocean_mask = (landsea_01 > 1)

# Apply ocean mask and valid-data mask to WOA23 field
woa_surface_ocean = woa_surface.where(woa_ocean_mask)

# Native-grid/global mean WOA23 on the 1-degree observational grid
woa_native_global_mean = float(
    compute_basic_stat(
        data=woa_surface_ocean,
        op="area_weighted_mean",
        dim=["lat", "lon"],
    ).compute()
)

# 4) Regrid each model climatology to WOA23 1-degree grid using conservative remapping
tos_regridded = {}
for model in models:
    tos_regridded[model] = regrid_horizontal(
        data=tos_clim_native[model],
        target="woa23",
        method="conservative",
        skipna=True,
        na_thres=0.5,
    )

# 5) Build one shared comparison mask across all three regridded models and WOA23
comparison_mask = build_comparison_mask(
    fields=[tos_regridded[m] for m in models] + [woa_surface_ocean],
    ocean_only=True,
    reduce_other_dims="all",
)

# Apply the common mask
woa_comp = woa_surface_ocean.where(comparison_mask)
tos_comp = {model: tos_regridded[model].where(comparison_mask) for model in models}

# Count comparison cells
comparison_cell_count = int(comparison_mask.sum().compute())

# 6) Shared-grid means and metrics against WOA23
shared_grid_means = {}
difference_metrics = {}
spatial_correlations = {}

for model in models:
    shared_grid_means[model] = {
        "model_mean": float(
            compute_basic_stat(
                data=tos_comp[model],
                op="area_weighted_mean",
                dim=["lat", "lon"],
            ).compute()
        ),
        "woa_mean": float(
            compute_basic_stat(
                data=woa_comp,
                op="area_weighted_mean",
                dim=["lat", "lon"],
            ).compute()
        ),
    }

    difference_metrics[model] = compute_difference_metrics(
        data=tos_comp[model],
        reference=woa_comp,
        dim=["lat", "lon"],
        weighting="area",
    )

    # Materialize scalar outputs
    difference_metrics[model] = {
        k: (float(v.compute()) if hasattr(v, "compute") else v)
        for k, v in difference_metrics[model].items()
        if k != "difference"
    } | {
        "difference": compute_difference_metrics(
            data=tos_comp[model],
            reference=woa_comp,
            dim=["lat", "lon"],
            weighting="area",
        )["difference"]
    }

    spatial_correlations[model] = analyze_correlation(
        var1=tos_comp[model],
        var2=woa_comp,
        dim=["lat", "lon"],
        method="pearson",
        weighting="area",
    )

    spatial_correlations[model] = {
        k: (float(v.compute()) if hasattr(v, "compute") else v)
        for k, v in spatial_correlations[model].items()
    }

# 7) Ranking by absolute mean bias
abs_bias_ranking = sorted(
    models,
    key=lambda m: abs(difference_metrics[m]["bias"])
)

# 8) Compact report table with requested precision
summary_rows = []
for model in models:
    summary_rows.append({
        "model": model,
        "native_global_mean_SST_degC": round(native_global_means[model], 2),
        "shared_grid_model_mean_degC": round(shared_grid_means[model]["model_mean"], 2),
        "shared_grid_WOA23_mean_degC": round(shared_grid_means[model]["woa_mean"], 2),
        "delta_SST_degC": round(difference_metrics[model]["bias"], 2),
        "RMSE_degC": round(difference_metrics[model]["rmse"], 3),
        "spatial_r_degCgrid": round(spatial_correlations[model]["correlation"], 3),
        "n_comp_cells": int(difference_metrics[model]["n_valid"]),
    })

results_table = pd.DataFrame(summary_rows)
woa_native_global_mean_rounded = round(woa_native_global_mean, 2)
abs_bias_ranking_with_values = [
    (m, round(difference_metrics[m]["bias"], 2), round(abs(difference_metrics[m]["bias"]), 2))
    for m in abs_bias_ranking
]

# 9) Bias fields for plotting
bias_fields = {}
for model in models:
    bias_fields[model] = (tos_comp[model] - woa_comp).assign_attrs(units="degC")

# 10) Variables explicitly available for downstream reporting / plotting
cesm2_native_global_mean = native_global_means["CESM2"]
gfdl_native_global_mean = native_global_means["GFDL-ESM4"]
mpi_native_global_mean = native_global_means["MPI-ESM1-2-LR"]

cesm2_bias = difference_metrics["CESM2"]["bias"]
gfdl_bias = difference_metrics["GFDL-ESM4"]["bias"]
mpi_bias = difference_metrics["MPI-ESM1-2-LR"]["bias"]

cesm2_rmse = difference_metrics["CESM2"]["rmse"]
gfdl_rmse = difference_metrics["GFDL-ESM4"]["rmse"]
mpi_rmse = difference_metrics["MPI-ESM1-2-LR"]["rmse"]

cesm2_corr = spatial_correlations["CESM2"]["correlation"]
gfdl_corr = spatial_correlations["GFDL-ESM4"]["correlation"]
mpi_corr = spatial_correlations["MPI-ESM1-2-LR"]["correlation"]
