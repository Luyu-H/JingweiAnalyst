import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# 1. GLODAP reference climatology
# ----------------------------------------------------------------------
glodap_tco2 = load_climatology_data(dataset="glodap", variable="TCO2")
glodap_talk = load_climatology_data(dataset="glodap", variable="TAlk")
glodap_ph = load_climatology_data(dataset="glodap", variable="pHtsinsitutp")

glodap_levels, glodap_bounds = load_depth_levels(dataset="glodap", variable="TCO2")
glodap_mask_0_50 = (glodap_bounds[:, 0] >= 0) & (glodap_bounds[:, 1] <= 50)
glodap_levels_0_50 = glodap_levels[glodap_mask_0_50]

glodap_tco2_0_50 = glodap_tco2.sel(depth=glodap_levels_0_50)
glodap_talk_0_50 = glodap_talk.sel(depth=glodap_levels_0_50)
glodap_ph_0_50 = glodap_ph.sel(depth=glodap_levels_0_50)

glodap_dic_mean = float(
    compute_basic_stat(glodap_tco2_0_50, op="volume_weighted_mean", dim=None)
    .compute()
    .values
)
glodap_alk_mean = float(
    compute_basic_stat(glodap_talk_0_50, op="volume_weighted_mean", dim=None)
    .compute()
    .values
)
glodap_ph_mean = float(
    compute_basic_stat(glodap_ph_0_50, op="volume_weighted_mean", dim=None)
    .compute()
    .values
)

glodap_dic_map = compute_basic_stat(
    glodap_tco2_0_50, op="mean", dim="depth"
).compute()

# ----------------------------------------------------------------------
# 2. Model 2002 annual means, convert DIC/TAlk to umol kg-1, 0-50 m
# ----------------------------------------------------------------------
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
model_results = {}

for model in models:
    # Load 2002 monthly fields
    dissic_2002 = load_regular_data(
        dataset="cmip6_historical",
        variable="dissic",
        model=model,
        time_range=["2002-01", "2002-12"],
    )
    talk_2002 = load_regular_data(
        dataset="cmip6_historical",
        variable="talk",
        model=model,
        time_range=["2002-01", "2002-12"],
        depth_range=[0, 100],
    )
    ph_2002 = load_regular_data(
        dataset="cmip6_historical",
        variable="ph",
        model=model,
        time_range=["2002-01", "2002-12"],
        depth_range=[0, 100],
    )
    thetao_2002 = load_regular_data(
        dataset="cmip6_historical",
        variable="thetao",
        model=model,
        time_range=["2002-01", "2002-12"],
    )
    so_2002 = load_regular_data(
        dataset="cmip6_historical",
        variable="so",
        model=model,
        time_range=["2002-01", "2002-12"],
    )

    # Annual means (day-weighted by temporal_aggregate)
    dissic_2002_ann = temporal_aggregate(dissic_2002, freq="mean", op="mean").compute()
    talk_2002_ann = temporal_aggregate(talk_2002, freq="mean", op="mean").compute()
    ph_2002_ann = temporal_aggregate(ph_2002, freq="mean", op="mean").compute()
    thetao_2002_ann = temporal_aggregate(thetao_2002, freq="mean", op="mean").compute()
    so_2002_ann = temporal_aggregate(so_2002, freq="mean", op="mean").compute()

    # Density anomaly for conversion
    sigma_2002 = compute_potential_density_anomaly(
        temp_data=thetao_2002_ann,
        sal_data=so_2002_ann,
        ref_pressure=0,
    ).compute()

    # Convert DIC and TAlk from mol m-3 to umol kg-1
    dissic_2002_umol = transform_unit(
        dissic_2002_ann,
        from_unit="mol m-3",
        to_unit="umol kg-1",
        density=sigma_2002,
    ).compute()

    sigma_2002_talk = sigma_2002.sel(depth=talk_2002_ann.depth)
    talk_2002_umol = transform_unit(
        talk_2002_ann,
        from_unit="mol m-3",
        to_unit="umol kg-1",
        density=sigma_2002_talk,
    ).compute()

    # Model depth levels wholly within 0-50 m
    model_levels, model_bounds = load_depth_levels(
        dataset="cmip6_historical", model=model, variable="dissic"
    )
    model_mask_0_50 = (model_bounds[:, 0] >= 0) & (model_bounds[:, 1] <= 50)
    model_levels_0_50 = model_levels[model_mask_0_50]

    dissic_2002_0_50 = dissic_2002_umol.sel(depth=model_levels_0_50)
    talk_2002_0_50 = talk_2002_umol.sel(depth=model_levels_0_50)
    ph_2002_0_50 = ph_2002_ann.sel(depth=model_levels_0_50)

    # Global 0-50 m means (volume-weighted over depth and horizontal)
    dic_mean = float(
        compute_basic_stat(dissic_2002_0_50, op="volume_weighted_mean", dim=None)
        .compute()
        .values
    )
    alk_mean = float(
        compute_basic_stat(talk_2002_0_50, op="volume_weighted_mean", dim=None)
        .compute()
        .values
    )
    ph_mean = float(
        compute_basic_stat(ph_2002_0_50, op="volume_weighted_mean", dim=None)
        .compute()
        .values
    )

    # 2-D 0-50 m mean DIC map (arithmetic mean over each cell's selected levels)
    dic_map = compute_basic_stat(
        dissic_2002_0_50, op="mean", dim="depth"
    ).compute()

    model_results[model] = {
        "dic_mean": dic_mean,
        "alk_mean": alk_mean,
        "ph_mean": ph_mean,
        "dic_diff": dic_mean - glodap_dic_mean,
        "alk_diff": alk_mean - glodap_alk_mean,
        "ph_diff": ph_mean - glodap_ph_mean,
        "dic_2002_umol": dissic_2002_umol,
        "dic_map": dic_map,
    }

# ----------------------------------------------------------------------
# 3. Period change 2002 minus 1971-2000 mean global 0-50 m DIC
# ----------------------------------------------------------------------
for model in models:
    dissic_hist = load_regular_data(
        dataset="cmip6_historical",
        variable="dissic",
        model=model,
        time_range=["1971-01", "2000-12"],
        depth_range=[0, 100],
    )
    thetao_hist = load_regular_data(
        dataset="cmip6_historical",
        variable="thetao",
        model=model,
        time_range=["1971-01", "2000-12"],
        depth_range=[0, 100],
    )
    so_hist = load_regular_data(
        dataset="cmip6_historical",
        variable="so",
        model=model,
        time_range=["1971-01", "2000-12"],
        depth_range=[0, 100],
    )

    dissic_hist_mean = temporal_aggregate(
        dissic_hist, freq="mean", op="mean"
    ).compute()
    thetao_hist_mean = temporal_aggregate(
        thetao_hist, freq="mean", op="mean"
    ).compute()
    so_hist_mean = temporal_aggregate(so_hist, freq="mean", op="mean").compute()

    sigma_hist = compute_potential_density_anomaly(
        temp_data=thetao_hist_mean,
        sal_data=so_hist_mean,
        ref_pressure=0,
    ).compute()

    dissic_hist_umol = transform_unit(
        dissic_hist_mean,
        from_unit="mol m-3",
        to_unit="umol kg-1",
        density=sigma_hist,
    ).compute()

    model_levels, model_bounds = load_depth_levels(
        dataset="cmip6_historical", model=model, variable="dissic"
    )
    model_mask_0_50 = (model_bounds[:, 0] >= 0) & (model_bounds[:, 1] <= 50)
    model_levels_0_50 = model_levels[model_mask_0_50]

    dissic_hist_0_50 = dissic_hist_umol.sel(depth=model_levels_0_50)
    dic_hist_mean = float(
        compute_basic_stat(dissic_hist_0_50, op="volume_weighted_mean", dim=None)
        .compute()
        .values
    )

    period_change = model_results[model]["dic_mean"] - dic_hist_mean
    ratio = period_change / model_results[model]["dic_diff"]

    model_results[model]["dic_hist_mean"] = dic_hist_mean
    model_results[model]["period_change"] = period_change
    model_results[model]["ratio"] = ratio

# ----------------------------------------------------------------------
# 4. Map comparison of 0-50 m mean DIC on the GLODAP 1-degree grid
# ----------------------------------------------------------------------
for model in models:
    model_map = model_results[model]["dic_map"]
    model_map_regrid = regrid_horizontal(
        model_map, target="glodap", method="conservative"
    ).compute()

    mask = build_comparison_mask(
        fields=[model_map_regrid, glodap_dic_map], ocean_only=True
    ).compute()

    metrics = compute_difference_metrics(
        data=model_map_regrid.where(mask),
        reference=glodap_dic_map.where(mask),
        dim=["lat", "lon"],
        weighting="area",
    )

    corr_dict = analyze_correlation(
        var1=model_map_regrid.where(mask),
        var2=glodap_dic_map.where(mask),
        dim=["lat", "lon"],
        method="pearson",
        weighting="area",
    )

    bias_map = (model_map_regrid - glodap_dic_map).where(mask)

    model_results[model].update(
        {
            "dic_map_regrid": model_map_regrid,
            "map_mask": mask,
            "bias_map": bias_map,
            "map_bias": float(metrics["bias"].values),
            "map_rmse": float(metrics["rmse"].values),
            "map_corr": float(corr_dict["correlation"].values),
        }
    )

# ----------------------------------------------------------------------
# 5. Vertical profile of global mean DIC bias at 0, 100, 500, 1000, 2000 m
# ----------------------------------------------------------------------
target_depths = [0, 100, 500, 1000, 2000]
glodap_dic_vert = regrid_vertical(
    glodap_tco2, target=target_depths, method="linear"
).compute()

for model in models:
    model_dic_regrid = regrid_horizontal(
        model_results[model]["dic_2002_umol"],
        target="glodap",
        method="conservative",
    ).compute()

    model_dic_vert = regrid_vertical(
        model_dic_regrid, target=target_depths, method="linear"
    ).compute()

    vertical_bias = {}
    for d in target_depths:
        m = model_dic_vert.sel(depth=d)
        g = glodap_dic_vert.sel(depth=d)
        mask_d = build_comparison_mask(
            fields=[m, g], ocean_only=True
        ).compute()
        metrics_d = compute_difference_metrics(
            data=m.where(mask_d),
            reference=g.where(mask_d),
            dim=["lat", "lon"],
            weighting="area",
        )
        vertical_bias[d] = float(metrics_d["bias"].values)

    model_results[model]["vertical_bias"] = vertical_bias

# ----------------------------------------------------------------------
# 6. Summarise, rank models, and assess sign consistency
# ----------------------------------------------------------------------
print("GLODAP 0-50 m means:")
print(
    f"  DIC = {glodap_dic_mean:.0f} umol kg-1, "
    f"TAlk = {glodap_alk_mean:.0f} umol kg-1, "
    f"pH = {glodap_ph_mean:.3f}"
)
print()

for model in models:
    r = model_results[model]
    print(f"{model}:")
    print(
        f"  0-50 m DIC = {r['dic_mean']:.0f} umol kg-1 "
        f"(model - GLODAP = {r['dic_diff']:.0f})"
    )
    print(
        f"  0-50 m TAlk = {r['alk_mean']:.0f} umol kg-1 "
        f"(model - GLODAP = {r['alk_diff']:.0f})"
    )
    print(
        f"  0-50 m pH = {r['ph_mean']:.3f} "
        f"(model - GLODAP = {r['ph_diff']:.3f})"
    )
    print(
        f"  DIC map: bias = {r['map_bias']:.0f} umol kg-1, "
        f"RMSE = {r['map_rmse']:.0f} umol kg-1, "
        f"spatial correlation = {r['map_corr']:.3f}"
    )
    print(
        "  Vertical DIC bias (umol kg-1): "
        + ", ".join(
            f"{d} m: {r['vertical_bias'][d]:.0f}" for d in target_depths
        )
    )
    print(
        f"  2002 minus 1971-2000 DIC change = {r['period_change']:.0f} umol kg-1"
    )
    print(
        f"  Ratio of period change to 2002 DIC bias = {r['ratio']:.2f}"
    )
    print()

# Simple ranking: sum of ranks for absolute DIC/TAlk/pH differences,
# map RMSE, and vertical RMS bias. Lower is better.
metrics_for_rank = pd.DataFrame(
    {
        model: {
            "abs_dic_diff": abs(model_results[model]["dic_diff"]),
            "abs_alk_diff": abs(model_results[model]["alk_diff"]),
            "abs_ph_diff_x1000": abs(model_results[model]["ph_diff"]) * 1000.0,
            "map_rmse": model_results[model]["map_rmse"],
            "vert_rmse": np.sqrt(
                np.mean(
                    [
                        model_results[model]["vertical_bias"][d] ** 2
                        for d in target_depths
                    ]
                )
            ),
        }
        for model in models
    }
).T

rank_table = metrics_for_rank.rank(axis=0, method="average")
metrics_for_rank["total_rank"] = rank_table.sum(axis=1)
best_model = metrics_for_rank["total_rank"].idxmin()

print("Model ranking (lower total rank is better):")
print(metrics_for_rank[["total_rank"]].sort_values("total_rank"))
print(f"Best overall model by combined rank: {best_model}")
print()

# Sign consistency across models for each variable
dic_signs = [np.sign(model_results[m]["dic_diff"]) for m in models]
alk_signs = [np.sign(model_results[m]["alk_diff"]) for m in models]
ph_signs = [np.sign(model_results[m]["ph_diff"]) for m in models]

dic_sign_consistent = len(set(dic_signs)) == 1
alk_sign_consistent = len(set(alk_signs)) == 1
ph_sign_consistent = len(set(ph_signs)) == 1

print("Sign consistency of model - GLODAP differences across models:")
print(f"  DIC: {dic_sign_consistent}  (signs = {dic_signs})")
print(f"  TAlk: {alk_sign_consistent}  (signs = {alk_signs})")
print(f"  pH: {ph_sign_consistent}  (signs = {ph_signs})")
print()

for model in models:
    r = model_results[model]
    print(
        f"{model} within-model signs: "
        f"DIC {np.sign(r['dic_diff']):+.0f}, "
        f"TAlk {np.sign(r['alk_diff']):+.0f}, "
        f"pH {np.sign(r['ph_diff']):+.0f}"
    )
