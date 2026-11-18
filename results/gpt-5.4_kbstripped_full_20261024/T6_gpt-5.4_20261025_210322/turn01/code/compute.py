# WOA23 observed nitrate climatology and nitracline depth
woa_no3 = load_climatology_data(
    dataset="woa23",
    variable="nitrate",
    time_slice="annual"
)

# Keep only ocean cells using the WOA23 land-sea mask definition
landsea_01 = load_landsea_mask()
woa_no3 = woa_no3.where(landsea_01 > 1)

# Observed nitracline depth: search in top 1000 m, allowing the gradient at the deepest
# searched level to use the level just below 1000 m by computing from the full loaded profile.
woa_nitracline = compute_nutricline_depth(woa_no3, min_depth=0, max_depth=1000).assign_attrs(units="m")
woa_nitracline = woa_nitracline.compute()

# Global mean nitracline depth for WOA23
woa_global_mean = float(
    compute_basic_stat(woa_nitracline, op="area_weighted_mean", dim=["lat", "lon"]).compute()
)

# Basin means for WOA23
basins = ["Pacific", "Atlantic", "Indian"]
woa_basin_means = {}
for basin in basins:
    woa_basin_field = select_basin(
        woa_nitracline,
        basin=basin,
        include_marginal_seas=True,
        depth_dependent=False
    )
    woa_basin_means[basin] = float(
        compute_basic_stat(woa_basin_field, op="area_weighted_mean", dim=["lat", "lon"]).compute()
    )

# CMIP6 models: 1971-2000 mean nitrate, density conversion to umol kg-1, nitracline depth
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
model_nitracline_native = {}
model_nitracline_regridded = {}
model_global_means = {}
model_global_diffs = {}
model_basin_means = {}
model_basin_diffs = {}
model_sections = {}

for model in models:
    # Nitrate
    no3 = load_regular_data(
        dataset="cmip6_historical",
        model=model,
        variable="no3",
        time_range=["1971-01", "2000-12"],
        chunk_time=12
    )

    # Temperature and salinity for in-situ density
    thetao = load_regular_data(
        dataset="cmip6_historical",
        model=model,
        variable="thetao",
        time_range=["1971-01", "2000-12"],
        chunk_time=12
    )
    so = load_regular_data(
        dataset="cmip6_historical",
        model=model,
        variable="so",
        time_range=["1971-01", "2000-12"],
        chunk_time=12
    )

    no3_mean = temporal_aggregate(no3, freq="mean", op="mean")
    thetao_mean = temporal_aggregate(thetao, freq="mean", op="mean")
    so_mean = temporal_aggregate(so, freq="mean", op="mean")

    sigma_insitu = compute_potential_density_anomaly(
        temp_data=thetao_mean,
        sal_data=so_mean,
        ref_pressure=None,
        temperature_type="auto"
    )

    no3_umolkg = transform_unit(
        data=no3_mean,
        from_unit=str(no3_mean.attrs.get("units", "mol m-3")),
        to_unit="umol/kg",
        density=sigma_insitu
    )

    model_ncl = compute_nutricline_depth(no3_umolkg, min_depth=0, max_depth=1000).assign_attrs(units="m")
    model_ncl = model_ncl.compute()
    model_nitracline_native[model] = model_ncl

    # Native-grid global mean
    model_global_means[model] = float(
        compute_basic_stat(model_ncl, op="area_weighted_mean", dim=["lat", "lon"]).compute()
    )
    model_global_diffs[model] = model_global_means[model] - woa_global_mean

    # Native-grid basin means using WOA23 basin mask projected to the model grid
    model_basin_means[model] = {}
    model_basin_diffs[model] = {}
    for basin in basins:
        basin_field = select_basin(
            model_ncl,
            basin=basin,
            include_marginal_seas=True,
            depth_dependent=False
        )
        model_basin_means[model][basin] = float(
            compute_basic_stat(basin_field, op="area_weighted_mean", dim=["lat", "lon"]).compute()
        )
        model_basin_diffs[model][basin] = model_basin_means[model][basin] - woa_basin_means[basin]

    # Regrid nitracline map to WOA23 grid for map metrics and Figure 1
    model_nitracline_regridded[model] = regrid_horizontal(
        model_ncl,
        target="woa23",
        method="bilinear",
        skipna=True,
        na_thres=0.5
    ).compute()

# Shared comparison mask on WOA23 grid across WOA23 and all regridded model maps
shared_mask = build_comparison_mask(
    fields=[woa_nitracline] + [model_nitracline_regridded[m] for m in models],
    ocean_only=True,
    reduce_other_dims="all"
).compute()

shared_mask_cell_count = int(shared_mask.sum().compute())

woa_nitracline_masked = woa_nitracline.where(shared_mask)

# Spatial metrics against WOA23 on the common mask
model_metrics = {}
for model in models:
    mod_masked = model_nitracline_regridded[model].where(shared_mask)

    corr = analyze_correlation(
        var1=mod_masked,
        var2=woa_nitracline_masked,
        dim=["lat", "lon"],
        method="pearson",
        weighting="area"
    )
    diffm = compute_difference_metrics(
        data=mod_masked,
        reference=woa_nitracline_masked,
        dim=["lat", "lon"],
        weighting="area"
    )

    model_metrics[model] = {
        "correlation": float(corr["correlation"].compute()),
        "rmse": float(diffm["rmse"].compute() if hasattr(diffm["rmse"], "compute") else diffm["rmse"]),
        "bias": float(diffm["bias"].compute() if hasattr(diffm["bias"], "compute") else diffm["bias"]),
        "n_valid": int(diffm["n_valid"]),
    }

# Rank best model: highest correlation, then lowest RMSE
best_model = sorted(
    models,
    key=lambda m: (-model_metrics[m]["correlation"], model_metrics[m]["rmse"], abs(model_global_diffs[m]))
)[0]

# Basin ordering consistency
woa_order = tuple(sorted(basins, key=lambda b: woa_basin_means[b]))
model_basin_order = {
    model: tuple(sorted(basins, key=lambda b: model_basin_means[model][b]))
    for model in models
}
model_matches_woa_order = {
    model: (model_basin_order[model] == woa_order)
    for model in models
}

# Equatorial Pacific sections for Figure 2, all in umol kg-1 and top 1000 m
# Region crosses the dateline: use loader lon_range directly.
woa_eqpac_no3 = load_climatology_data(
    dataset="woa23",
    variable="nitrate",
    time_slice="annual",
    lat_range=[-2, 2],
    lon_range=[120, 280],
    depth_range=[0, 1100]
).where(load_landsea_mask().sel(lat=slice(-2, 2), lon=slice(120, 280)) > 1)

# For a latitude-band section, collapse latitude only with an unweighted mean to avoid
# invalid partial area-weighting over just one horizontal axis.
woa_eqpac_section = compute_basic_stat(
    woa_eqpac_no3,
    op="mean",
    dim=["lat"]
).compute()

woa_eqpac_ncl = compute_basic_stat(
    woa_nitracline.where((woa_nitracline.lat >= -2) & (woa_nitracline.lat <= 2)).sel(lon=slice(120, 280)),
    op="mean",
    dim=["lat"]
).compute()

section_order = ["WOA23"] + models
section_fields = {"WOA23": woa_eqpac_section}
section_ncl_lines = {"WOA23": woa_eqpac_ncl}

for model in models:
    no3_eq = load_regular_data(
        dataset="cmip6_historical",
        model=model,
        variable="no3",
        time_range=["1971-01", "2000-12"],
        lat_range=[-2, 2],
        lon_range=[120, 280],
        depth_range=[0, 1100],
        chunk_time=12
    )
    thetao_eq = load_regular_data(
        dataset="cmip6_historical",
        model=model,
        variable="thetao",
        time_range=["1971-01", "2000-12"],
        lat_range=[-2, 2],
        lon_range=[120, 280],
        depth_range=[0, 1100],
        chunk_time=12
    )
    so_eq = load_regular_data(
        dataset="cmip6_historical",
        model=model,
        variable="so",
        time_range=["1971-01", "2000-12"],
        lat_range=[-2, 2],
        lon_range=[120, 280],
        depth_range=[0, 1100],
        chunk_time=12
    )

    no3_eq_mean = temporal_aggregate(no3_eq, freq="mean", op="mean")
    thetao_eq_mean = temporal_aggregate(thetao_eq, freq="mean", op="mean")
    so_eq_mean = temporal_aggregate(so_eq, freq="mean", op="mean")

    sigma_eq = compute_potential_density_anomaly(
        temp_data=thetao_eq_mean,
        sal_data=so_eq_mean,
        ref_pressure=None,
        temperature_type="auto"
    )
    no3_eq_umolkg = transform_unit(
        data=no3_eq_mean,
        from_unit=str(no3_eq_mean.attrs.get("units", "mol m-3")),
        to_unit="umol/kg",
        density=sigma_eq
    )

    lat_coord = no3_eq_umolkg.attrs.get("lat_coord", "lat")

    # Latitude-band mean section: use unweighted mean over latitude-like axis only
    if "lat" in no3_eq_umolkg.dims:
        sec = compute_basic_stat(no3_eq_umolkg, op="mean", dim=["lat"]).compute()
    else:
        lat_name = no3_eq_umolkg[lat_coord].dims[0]
        sec = compute_basic_stat(no3_eq_umolkg, op="mean", dim=[lat_name]).compute()

    # Robust nitracline line for plotting: use already computed regridded nitracline and
    # collapse the 2S-2N band with an unweighted mean over latitude.
    ncl_reg = model_nitracline_regridded[model].sel(lat=slice(-2, 2), lon=slice(120, 280))
    ncl_line = compute_basic_stat(ncl_reg, op="mean", dim=["lat"]).compute()

    section_fields[model] = sec
    section_ncl_lines[model] = ncl_line

# Round/report tables
global_mean_table = pd.DataFrame({
    "source": ["WOA23"] + models,
    "global_mean_m": [round(woa_global_mean)] + [round(model_global_means[m]) for m in models],
    "difference_from_WOA23_m": [0] + [round(model_global_diffs[m]) for m in models],
})

basin_rows = []
for source in ["WOA23"] + models:
    for basin in basins:
        if source == "WOA23":
            basin_rows.append({
                "source": source,
                "basin": basin,
                "mean_m": round(woa_basin_means[basin]),
                "difference_from_WOA23_m": 0
            })
        else:
            basin_rows.append({
                "source": source,
                "basin": basin,
                "mean_m": round(model_basin_means[source][basin]),
                "difference_from_WOA23_m": round(model_basin_diffs[source][basin])
            })
basin_mean_table = pd.DataFrame(basin_rows)

metrics_table = pd.DataFrame({
    "model": models,
    "pearson_r": [round(model_metrics[m]["correlation"], 3) for m in models],
    "rmse_m": [round(model_metrics[m]["rmse"], 1) for m in models],
    "n_cells_common_mask": [shared_mask_cell_count for _ in models],
    "matches_WOA23_basin_order": [model_matches_woa_order[m] for m in models],
})

summary = {
    "woa_global_mean_m": round(woa_global_mean),
    "model_global_means_m": {m: round(model_global_means[m]) for m in models},
    "model_minus_woa_global_m": {m: round(model_global_diffs[m]) for m in models},
    "woa_basin_means_m": {b: round(woa_basin_means[b]) for b in basins},
    "model_basin_means_m": {m: {b: round(model_basin_means[m][b]) for b in basins} for m in models},
    "model_minus_woa_basin_m": {m: {b: round(model_basin_diffs[m][b]) for b in basins} for m in models},
    "map_metrics": {m: {"pearson_r": round(model_metrics[m]["correlation"], 3), "rmse_m": round(model_metrics[m]["rmse"], 1)} for m in models},
    "shared_mask_cell_count": shared_mask_cell_count,
    "best_model": best_model,
    "woa_basin_order_shallow_to_deep": woa_order,
    "model_basin_order_shallow_to_deep": model_basin_order,
    "model_matches_woa_order": model_matches_woa_order,
}

print(global_mean_table.to_string(index=False))
print()
print(basin_mean_table.to_string(index=False))
print()
print(metrics_table.to_string(index=False))
print()
print(summary)
