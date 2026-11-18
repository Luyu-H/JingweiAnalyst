# Computation code only

models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
time_2002 = ["2002-01", "2002-12"]
time_hist = ["1971-01", "2000-12"]
depth_0_50 = [0.0, 50.0]
target_depths = [0.0, 100.0, 500.0, 1000.0, 2000.0]

model_fields_2002 = {}
model_fields_hist = {}
model_means_2002 = {}
model_means_hist = {}
model_density_2002 = {}
model_density_hist = {}
model_converted_2002 = {}
model_converted_hist = {}
model_scalar_means = {}
model_0_50_maps_native = {}
model_0_50_maps_regridded = {}
depth_bias_profiles = {}
shared_mask = None

# Load GLODAP climatology
glodap_tco2 = load_climatology_data(
    dataset="glodap",
    variable="TCO2",
    depth_range=depth_0_50,
)
glodap_talk = load_climatology_data(
    dataset="glodap",
    variable="TAlk",
    depth_range=depth_0_50,
)
glodap_ph = load_climatology_data(
    dataset="glodap",
    variable="pHtsinsitutp",
    depth_range=depth_0_50,
)

# Also load full-depth GLODAP TCO2 for depthwise comparison
glodap_tco2_full = load_climatology_data(
    dataset="glodap",
    variable="TCO2",
)

# 0-50 m GLODAP scalar means (native grid, volume-weighted)
glodap_dic_0_50_mean = compute_basic_stat(glodap_tco2, op="volume_weighted_mean").compute()
glodap_talk_0_50_mean = compute_basic_stat(glodap_talk, op="volume_weighted_mean").compute()
glodap_ph_0_50_mean = compute_basic_stat(glodap_ph, op="volume_weighted_mean").compute()

glodap_scalar_means = {
    "dic_0_50_mean_umolkg": float(glodap_dic_0_50_mean.values),
    "talk_0_50_mean_umolkg": float(glodap_talk_0_50_mean.values),
    "ph_0_50_mean_total": float(glodap_ph_0_50_mean.values),
}

# GLODAP 0-50 m DIC map: each cell mean over its own 0-50 m levels
glodap_dic_0_50_map = compute_basic_stat(
    glodap_tco2, op="mean", dim="depth"
).assign_attrs(units=glodap_tco2.attrs.get("units", "umol kg-1"))

# GLODAP depth slices at requested depths
glodap_tco2_depths = regrid_vertical(glodap_tco2_full, target=target_depths, method="linear", fill_surface=True)

# Load model data, aggregate in time, derive density, convert units
for model in models:
    # 2002 monthly fields
    dissic_2002 = load_regular_data(
        dataset="cmip6_historical",
        model=model,
        variable="dissic",
        time_range=time_2002,
    )
    talk_2002 = load_regular_data(
        dataset="cmip6_historical",
        model=model,
        variable="talk",
        time_range=time_2002,
    )
    ph_2002 = load_regular_data(
        dataset="cmip6_historical",
        model=model,
        variable="ph",
        time_range=time_2002,
    )
    thetao_2002 = load_regular_data(
        dataset="cmip6_historical",
        model=model,
        variable="thetao",
        time_range=time_2002,
    )
    so_2002 = load_regular_data(
        dataset="cmip6_historical",
        model=model,
        variable="so",
        time_range=time_2002,
    )

    # 1971-2000 monthly fields
    dissic_hist = load_regular_data(
        dataset="cmip6_historical",
        model=model,
        variable="dissic",
        time_range=time_hist,
    )
    thetao_hist = load_regular_data(
        dataset="cmip6_historical",
        model=model,
        variable="thetao",
        time_range=time_hist,
    )
    so_hist = load_regular_data(
        dataset="cmip6_historical",
        model=model,
        variable="so",
        time_range=time_hist,
    )

    # Cell volumes support field
    volcello_0_50 = load_volumes_for_cmip6(
        dataset="cmip6_historical",
        model=model,
        depth_range=depth_0_50,
        variable="thetao",
    )

    model_fields_2002[model] = {
        "dissic": dissic_2002,
        "talk": talk_2002,
        "ph": ph_2002,
        "thetao": thetao_2002,
        "so": so_2002,
        "volcello_0_50": volcello_0_50,
    }
    model_fields_hist[model] = {
        "dissic": dissic_hist,
        "thetao": thetao_hist,
        "so": so_hist,
    }

    # Time means
    model_means_2002[model] = {
        "dissic": temporal_aggregate(dissic_2002, freq="mean", op="mean"),
        "talk": temporal_aggregate(talk_2002, freq="mean", op="mean"),
        "ph": temporal_aggregate(ph_2002, freq="mean", op="mean"),
        "thetao": temporal_aggregate(thetao_2002, freq="mean", op="mean"),
        "so": temporal_aggregate(so_2002, freq="mean", op="mean"),
    }
    model_means_hist[model] = {
        "dissic": temporal_aggregate(dissic_hist, freq="mean", op="mean"),
        "thetao": temporal_aggregate(thetao_hist, freq="mean", op="mean"),
        "so": temporal_aggregate(so_hist, freq="mean", op="mean"),
    }

    # In-situ density anomaly for unit conversion
    model_density_2002[model] = compute_potential_density_anomaly(
        temp_data=model_means_2002[model]["thetao"],
        sal_data=model_means_2002[model]["so"],
        ref_pressure=None,
    )
    model_density_hist[model] = compute_potential_density_anomaly(
        temp_data=model_means_hist[model]["thetao"],
        sal_data=model_means_hist[model]["so"],
        ref_pressure=None,
    )

    # Convert to umol kg-1
    model_converted_2002[model] = {
        "dissic": transform_unit(
            model_means_2002[model]["dissic"],
            from_unit=model_means_2002[model]["dissic"].attrs.get("units", "mol m-3"),
            to_unit="umol kg-1",
            density=model_density_2002[model],
        ),
        "talk": transform_unit(
            model_means_2002[model]["talk"],
            from_unit=model_means_2002[model]["talk"].attrs.get("units", "mol m-3"),
            to_unit="umol kg-1",
            density=model_density_2002[model],
        ),
        "ph": model_means_2002[model]["ph"],
    }
    model_converted_hist[model] = {
        "dissic": transform_unit(
            model_means_hist[model]["dissic"],
            from_unit=model_means_hist[model]["dissic"].attrs.get("units", "mol m-3"),
            to_unit="umol kg-1",
            density=model_density_hist[model],
        ),
    }

    # Native-grid 0-50 m scalar means
    model_dic_0_50_mean = compute_basic_stat(
        model_converted_2002[model]["dissic"].sel(depth=slice(depth_0_50[0], depth_0_50[1])),
        op="volume_weighted_mean",
    ).compute()
    model_talk_0_50_mean = compute_basic_stat(
        model_converted_2002[model]["talk"].sel(depth=slice(depth_0_50[0], depth_0_50[1])),
        op="volume_weighted_mean",
    ).compute()
    model_ph_0_50_mean = compute_basic_stat(
        model_converted_2002[model]["ph"].sel(depth=slice(depth_0_50[0], depth_0_50[1])),
        op="volume_weighted_mean",
    ).compute()
    model_dic_hist_0_50_mean = compute_basic_stat(
        model_converted_hist[model]["dissic"].sel(depth=slice(depth_0_50[0], depth_0_50[1])),
        op="volume_weighted_mean",
    ).compute()

    dic_diff = float(model_dic_0_50_mean.values) - glodap_scalar_means["dic_0_50_mean_umolkg"]
    talk_diff = float(model_talk_0_50_mean.values) - glodap_scalar_means["talk_0_50_mean_umolkg"]
    ph_diff = float(model_ph_0_50_mean.values) - glodap_scalar_means["ph_0_50_mean_total"]
    dic_change = float(model_dic_0_50_mean.values) - float(model_dic_hist_0_50_mean.values)
    dic_change_ratio = dic_change / dic_diff

    model_scalar_means[model] = {
        "dic_0_50_mean_umolkg": float(model_dic_0_50_mean.values),
        "dic_minus_glodap_umolkg": float(dic_diff),
        "talk_0_50_mean_umolkg": float(model_talk_0_50_mean.values),
        "talk_minus_glodap_umolkg": float(talk_diff),
        "ph_0_50_mean_total": float(model_ph_0_50_mean.values),
        "ph_minus_glodap": float(ph_diff),
        "dic_1971_2000_0_50_mean_umolkg": float(model_dic_hist_0_50_mean.values),
        "dic_change_2002_minus_1971_2000_umolkg": float(dic_change),
        "dic_change_to_glodap_bias_ratio": float(dic_change_ratio),
    }

    # Native-grid 0-50 m DIC map
    model_0_50_maps_native[model] = compute_basic_stat(
        model_converted_2002[model]["dissic"].sel(depth=slice(depth_0_50[0], depth_0_50[1])),
        op="mean",
        dim="depth",
    ).assign_attrs(units="umol kg-1")

    # Regrid model 0-50 m DIC maps to GLODAP grid
    model_0_50_maps_regridded[model] = regrid_horizontal(
        model_0_50_maps_native[model],
        target=glodap_dic_0_50_map,
        method="conservative",
        skipna=True,
        na_thres=0.5,
    )

# Shared comparison mask across GLODAP + all three regridded models
shared_mask = build_comparison_mask(
    fields=[glodap_dic_0_50_map] + [model_0_50_maps_regridded[m] for m in models],
    ocean_only=False,
    reduce_other_dims="all",
).compute()

masked_glodap_dic_0_50_map = glodap_dic_0_50_map.where(shared_mask)

# Map metrics and bias maps
map_metrics = {}
dic_bias_maps = {}

for model in models:
    masked_model_map = model_0_50_maps_regridded[model].where(shared_mask)
    metrics = compute_difference_metrics(
        data=masked_model_map,
        reference=masked_glodap_dic_0_50_map,
        dim=None,
        weighting="auto",
    )
    corr = analyze_correlation(
        var1=masked_model_map,
        var2=masked_glodap_dic_0_50_map,
        dim=["lat", "lon"],
        method="pearson",
        weighting="area",
        min_periods=10,
    )

    mean_bias = float(metrics["bias"].compute().values if hasattr(metrics["bias"], "compute") else metrics["bias"])
    rmse = float(metrics["rmse"].compute().values if hasattr(metrics["rmse"], "compute") else metrics["rmse"])
    spatial_corr = float(corr["correlation"].compute().values)

    map_metrics[model] = {
        "mean_bias_umolkg": mean_bias,
        "rmse_umolkg": rmse,
        "spatial_correlation": spatial_corr,
        "n_valid": int(metrics["n_valid"]),
        "weighting": metrics["weighting"],
    }

    dic_bias_maps[model] = metrics["difference"].compute().assign_attrs(units="umol kg-1")

# Depthwise DIC bias comparison on GLODAP grid
for model in models:
    model_dic_depths = regrid_vertical(
        model_converted_2002[model]["dissic"],
        target=target_depths,
        method="linear",
        fill_surface=True,
    )
    model_dic_depths_rg = regrid_horizontal(
        model_dic_depths,
        target=glodap_tco2_depths,
        method="conservative",
        skipna=True,
        na_thres=0.5,
    )

    depth_bias = {}
    for depth_val in target_depths:
        model_slice = model_dic_depths_rg.sel(depth=depth_val)
        glodap_slice = glodap_tco2_depths.sel(depth=depth_val)

        metrics_depth = compute_difference_metrics(
            data=model_slice,
            reference=glodap_slice,
            dim=None,
            weighting="auto",
        )
        depth_bias[depth_val] = float(
            metrics_depth["bias"].compute().values
            if hasattr(metrics_depth["bias"], "compute")
            else metrics_depth["bias"]
        )

    depth_bias_profiles[model] = depth_bias

# Rounded summary tables for reporting
summary_table = pd.DataFrame(
    {
        model: {
            "DIC_0_50_mean_umolkg": round(model_scalar_means[model]["dic_0_50_mean_umolkg"]),
            "DIC_minus_GLODAP_umolkg": round(model_scalar_means[model]["dic_minus_glodap_umolkg"]),
            "TAlk_0_50_mean_umolkg": round(model_scalar_means[model]["talk_0_50_mean_umolkg"]),
            "TAlk_minus_GLODAP_umolkg": round(model_scalar_means[model]["talk_minus_glodap_umolkg"]),
            "pH_0_50_mean_total": round(model_scalar_means[model]["ph_0_50_mean_total"], 3),
            "pH_minus_GLODAP": round(model_scalar_means[model]["ph_minus_glodap"], 3),
            "DIC_change_2002_minus_1971_2000_umolkg": round(model_scalar_means[model]["dic_change_2002_minus_1971_2000_umolkg"]),
            "DIC_change_to_GLODAP_bias_ratio": round(model_scalar_means[model]["dic_change_to_glodap_bias_ratio"], 2),
        }
        for model in models
    }
).T

glodap_summary = {
    "DIC_0_50_mean_umolkg": round(glodap_scalar_means["dic_0_50_mean_umolkg"]),
    "TAlk_0_50_mean_umolkg": round(glodap_scalar_means["talk_0_50_mean_umolkg"]),
    "pH_0_50_mean_total": round(glodap_scalar_means["ph_0_50_mean_total"], 3),
}

map_metrics_table = pd.DataFrame(
    {
        model: {
            "mean_bias_umolkg": round(map_metrics[model]["mean_bias_umolkg"]),
            "rmse_umolkg": round(map_metrics[model]["rmse_umolkg"]),
            "spatial_correlation": round(map_metrics[model]["spatial_correlation"], 3),
            "n_valid": map_metrics[model]["n_valid"],
        }
        for model in models
    }
).T

depth_bias_table = pd.DataFrame(depth_bias_profiles).T[target_depths].rename(
    columns={d: f"{int(d)} m" for d in target_depths}
)
depth_bias_table = depth_bias_table.applymap(lambda x: round(x))

# Consistency of sign across DIC, alkalinity, and pH
sign_consistency = {}
for model in models:
    signs = np.sign([
        model_scalar_means[model]["dic_minus_glodap_umolkg"],
        model_scalar_means[model]["talk_minus_glodap_umolkg"],
        model_scalar_means[model]["ph_minus_glodap"],
    ])
    sign_consistency[model] = bool((signs[0] == signs[1]) and (signs[1] == signs[2]))

# Simple overall ranking diagnostic
overall_score_table = pd.DataFrame(index=models)
overall_score_table["abs_dic_bias_scalar"] = [abs(model_scalar_means[m]["dic_minus_glodap_umolkg"]) for m in models]
overall_score_table["abs_talk_bias_scalar"] = [abs(model_scalar_means[m]["talk_minus_glodap_umolkg"]) for m in models]
overall_score_table["abs_ph_bias_scalar"] = [abs(model_scalar_means[m]["ph_minus_glodap"]) for m in models]
overall_score_table["abs_dic_map_bias"] = [abs(map_metrics[m]["mean_bias_umolkg"]) for m in models]
overall_score_table["dic_map_rmse"] = [map_metrics[m]["rmse_umolkg"] for m in models]
overall_score_table["dic_map_corr_penalty"] = [1.0 - map_metrics[m]["spatial_correlation"] for m in models]
overall_score_table["mean_abs_depth_bias"] = [
    np.mean([abs(depth_bias_profiles[m][d]) for d in target_depths]) for m in models
]
overall_score_table["score"] = overall_score_table.sum(axis=1)
best_model = overall_score_table["score"].idxmin()

# Interpretation helper text
ratio_interpretation = {}
for model in models:
    ratio = model_scalar_means[model]["dic_change_to_glodap_bias_ratio"]
    if abs(ratio) < 0.25:
        txt = "small period-change contribution relative to apparent 2002 bias"
    elif abs(ratio) < 0.75:
        txt = "moderate period-change contribution relative to apparent 2002 bias"
    else:
        txt = "large period-change contribution relative to apparent 2002 bias"
    ratio_interpretation[model] = txt
