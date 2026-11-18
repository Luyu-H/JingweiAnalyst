# Load WOA23 annual climatology in the 200-1000 m layer
woa_no3 = load_climatology_data(
    dataset="woa23",
    variable="nitrate",
    time_slice="annual",
    depth_range=[200.0, 1000.0],
)
woa_po4 = load_climatology_data(
    dataset="woa23",
    variable="phosphate",
    time_slice="annual",
    depth_range=[200.0, 1000.0],
)

# WOA23 basin mask requested in the plan; loaded here for provenance/checking
basinmask_3d = load_basin_mask()

# CMIP6 models to analyse
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
time_range = ["1971-01", "2000-12"]
depth_range = [200.0, 1000.0]

# Storage
sources = {}
global_stats = {}
regression_stats = {}
between_14_18 = {}
basin_headline_np = {}
below12_stats = {}

# Add WOA23 as one source
sources["WOA23"] = {
    "no3_native": woa_no3,
    "po4_native": woa_po4,
    "no3_umolkg": woa_no3.assign_attrs(units="umol kg-1"),
    "po4_umolkg": woa_po4.assign_attrs(units="umol kg-1"),
}

# Load and climatologise each CMIP6 model
for model in models:
    no3_ts = load_regular_data(
        dataset="cmip6_historical",
        model=model,
        variable="no3",
        time_range=time_range,
        depth_range=depth_range,
        chunk_time=12,
    )
    po4_ts = load_regular_data(
        dataset="cmip6_historical",
        model=model,
        variable="po4",
        time_range=time_range,
        depth_range=depth_range,
        chunk_time=12,
    )

    # Only GFDL-ESM4 publishes time-varying volcello; fixed-volume models must omit time_range
    if model == "GFDL-ESM4":
        volcello = load_volumes_for_cmip6(
            dataset="cmip6_historical",
            model=model,
            time_range=time_range,
            depth_range=depth_range,
            variable="thetao",
        )
    else:
        volcello = load_volumes_for_cmip6(
            dataset="cmip6_historical",
            model=model,
            depth_range=depth_range,
            variable="thetao",
        )

    no3_clim = temporal_aggregate(no3_ts, freq="climatology", op="mean").mean("month")
    po4_clim = temporal_aggregate(po4_ts, freq="climatology", op="mean").mean("month")

    no3_umolkg = transform_unit(
        no3_clim,
        from_unit=no3_clim.attrs.get("units", "mol m-3"),
        to_unit="umol kg-1",
        density=1.025,
    )
    po4_umolkg = transform_unit(
        po4_clim,
        from_unit=po4_clim.attrs.get("units", "mol m-3"),
        to_unit="umol kg-1",
        density=1.025,
    )

    sources[model] = {
        "no3_native": no3_clim,
        "po4_native": po4_clim,
        "no3_umolkg": no3_umolkg,
        "po4_umolkg": po4_umolkg,
        "volcello": volcello,
    }

# Compute N:P fields and all requested statistics
for name, src in sources.items():
    no3_native = src["no3_native"]
    po4_native = src["po4_native"]

    np_ratio = compute_ratio(no3_native, po4_native).assign_attrs(units="1")
    src["np_ratio"] = np_ratio

    # Global means over 200-1000 m
    mean_np_headline = compute_basic_stat(
        np_ratio,
        op="volume_weighted_mean",
        dim=None,
    )
    mean_no3 = compute_basic_stat(
        no3_native,
        op="volume_weighted_mean",
        dim=None,
    )
    mean_po4 = compute_basic_stat(
        po4_native,
        op="volume_weighted_mean",
        dim=None,
    )

    mean_np_headline_val = float(mean_np_headline.compute())
    mean_no3_val = float(mean_no3.compute())
    mean_po4_val = float(mean_po4.compute())
    ratio_of_means_val = mean_no3_val / mean_po4_val

    global_stats[name] = {
        "mean_nitrate": mean_no3_val,
        "mean_phosphate": mean_po4_val,
        "headline_mean_np": mean_np_headline_val,
        "ratio_of_mean_no3_to_mean_po4": ratio_of_means_val,
    }

    # Regression in umol kg-1 with equal cell weighting
    reg = fit_multivariate_relation(
        target=src["no3_umolkg"],
        predictors=[src["po4_umolkg"]],
        method="linear",
        validation="blocked",
    )

    coeffs = reg["coefficients"]
    if isinstance(coeffs, dict):
        slope_val = float(coeffs.get("predictor_0"))
        intercept_val = float(coeffs.get("intercept"))
    else:
        coeffs_arr = np.asarray(coeffs)
        slope_val = float(coeffs_arr[0])
        intercept_val = float(coeffs_arr[1])

    regression_stats[name] = {
        "slope": slope_val,
        "intercept": intercept_val,
        "r2": float(reg["metrics"]["r2"]),
        "full_result": reg,
    }

    # Percent volume with 14 < N:P < 18
    frac_14_18 = compute_exceedance_fraction(
        np_ratio,
        threshold=2.0,
        comparison="abs_less",
        reference=16.0,
        dim=None,
        weighting="auto",
    )
    between_14_18[name] = {
        "fraction_percent": float(np.asarray(frac_14_18["fraction_percent"]).item())
    }

    # Basin headline mean N:P
    basin_headline_np[name] = {}
    for basin in ["Pacific", "Atlantic", "Indian"]:
        basin_ratio = select_basin(
            np_ratio,
            basin=basin,
            include_marginal_seas=True,
            depth_dependent=False,
        )
        basin_mean = compute_basic_stat(
            basin_ratio,
            op="volume_weighted_mean",
            dim=None,
        )
        basin_headline_np[name][basin] = float(basin_mean.compute())

    # Percent volume with N:P < 12: global, Pacific, Indian
    below12_stats[name] = {}
    frac_below12_global = compute_exceedance_fraction(
        np_ratio,
        threshold=12.0,
        comparison="less",
        reference=0.0,
        dim=None,
        weighting="auto",
    )
    below12_stats[name]["global"] = float(np.asarray(frac_below12_global["fraction_percent"]).item())

    for basin in ["Pacific", "Indian"]:
        basin_ratio = select_basin(
            np_ratio,
            basin=basin,
            include_marginal_seas=True,
            depth_dependent=False,
        )
        frac_below12_basin = compute_exceedance_fraction(
            basin_ratio,
            threshold=12.0,
            comparison="less",
            reference=0.0,
            dim=None,
            weighting="auto",
        )
        below12_stats[name][basin] = float(np.asarray(frac_below12_basin["fraction_percent"]).item())

    # 200-1000 m mean N:P minus 16 for map display
    if name == "WOA23":
        mean_np_map = np_ratio.mean("depth")
    else:
        mean_np_map = compute_basic_stat(
            np_ratio,
            op="volume_weighted_mean",
            dim="depth",
        )
    src["np_minus_16_map_native"] = (mean_np_map - 16.0).assign_attrs(units="1")

# Regrid map fields to a common regular grid for display only
map_fields_display = {}
for name, src in sources.items():
    native_map = src["np_minus_16_map_native"]
    if name == "WOA23":
        map_fields_display[name] = native_map
    else:
        map_fields_display[name] = regrid_horizontal(
            native_map,
            target="woa23",
            method="bilinear",
            skipna=True,
            na_thres=0.5,
        )

# Materialise map arrays for robust common color scaling
map_arrays_for_scale = []
for name in ["WOA23", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]:
    arr = map_fields_display[name].compute()
    map_fields_display[name] = arr
    map_arrays_for_scale.append(np.ravel(arr.values))

pooled_abs = np.abs(np.concatenate([a[np.isfinite(a)] for a in map_arrays_for_scale]))
map_vmax = float(np.nanpercentile(pooled_abs, 98))

# Build printable summary tables
summary_rows = []
for name in ["WOA23", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]:
    headline = global_stats[name]["headline_mean_np"]
    if np.isclose(headline, 16.0, atol=0.1):
        status = "at Redfield"
    elif headline > 16.0:
        status = "above Redfield"
    else:
        status = "below Redfield"

    summary_rows.append({
        "source": name,
        "headline_mean_N_to_P": round(global_stats[name]["headline_mean_np"], 2),
        "ratio_of_mean_NO3_to_mean_PO4": round(global_stats[name]["ratio_of_mean_no3_to_mean_po4"], 2),
        "reg_slope_a": round(regression_stats[name]["slope"], 2),
        "reg_intercept_b_umolkg": round(regression_stats[name]["intercept"], 2),
        "reg_r2": round(regression_stats[name]["r2"], 3),
        "pct_volume_14_to_18": round(between_14_18[name]["fraction_percent"], 1),
        "Pacific_mean_NP": round(basin_headline_np[name]["Pacific"], 2),
        "Atlantic_mean_NP": round(basin_headline_np[name]["Atlantic"], 2),
        "Indian_mean_NP": round(basin_headline_np[name]["Indian"], 2),
        "pct_volume_below12_global": round(below12_stats[name]["global"], 1),
        "pct_volume_below12_Pacific": round(below12_stats[name]["Pacific"], 1),
        "pct_volume_below12_Indian": round(below12_stats[name]["Indian"], 1),
        "stoichiometry_vs_redfield": status,
    })

summary_table = pd.DataFrame(summary_rows).set_index("source")

# Compare each model with WOA23 for low-N:P reproduction
woa_low_global = below12_stats["WOA23"]["global"]
woa_low_pac = below12_stats["WOA23"]["Pacific"]
woa_low_ind = below12_stats["WOA23"]["Indian"]

model_assessment = {}
for model in models:
    dg = below12_stats[model]["global"] - woa_low_global
    dp = below12_stats[model]["Pacific"] - woa_low_pac
    di = below12_stats[model]["Indian"] - woa_low_ind

    basin_match = (abs(dp) <= 5.0) and (abs(di) <= 5.0)
    global_match = abs(dg) <= 5.0

    if global_match and basin_match:
        verdict = "reproduces observed low-N:P volume and basin distribution"
    elif global_match:
        verdict = "reproduces global low-N:P volume but not basin distribution"
    else:
        verdict = "does not reproduce observed low-N:P volume and basin distribution"

    model_assessment[model] = {
        "delta_global_pct_points": round(dg, 1),
        "delta_pacific_pct_points": round(dp, 1),
        "delta_indian_pct_points": round(di, 1),
        "assessment": verdict,
    }

model_assessment_table = pd.DataFrame(model_assessment).T

# Convenience objects for plotting Figure 1
scatter_order = ["WOA23", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
scatter_samples = {}
all_po4 = []
all_no3 = []

for name in scatter_order:
    po4 = sources[name]["po4_umolkg"].compute()
    no3 = sources[name]["no3_umolkg"].compute()
    valid = np.isfinite(po4.values) & np.isfinite(no3.values)
    x = po4.values[valid].ravel()
    y = no3.values[valid].ravel()
    scatter_samples[name] = {"x": x, "y": y}
    all_po4.append(x)
    all_no3.append(y)

all_po4_concat = np.concatenate(all_po4)
all_no3_concat = np.concatenate(all_no3)

xmax = float(np.nanpercentile(all_po4_concat, 99.5))
ymax = float(np.nanpercentile(all_no3_concat, 99.5))
xlim = (0.0, xmax * 1.02)
ylim = (0.0, ymax * 1.02)

# Text outputs
print("Summary table for intermediate-water nitrate:phosphate stoichiometry (200-1000 m)")
print(summary_table.to_string())

print("\nModel assessment against WOA23 low-N:P volume and basin distribution")
print(model_assessment_table.to_string())
