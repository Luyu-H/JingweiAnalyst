# Computation code only

dataset = selected_data
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
time_range = ["1971-01", "2000-12"]
target_months = [3, 9]
month_names = {3: "March", 9: "September"}

results = {}
diff_maps_native = {}
pair_fields_native = {}
scatter_inputs = {}
extremes = {}
display_maps = []

for model in models:
    # Load monthly model data on native grid
    thetao = load_regular_data(
        dataset=dataset,
        variable="thetao",
        model=model,
        time_range=time_range,
        chunk_time=12,
    )
    so = load_regular_data(
        dataset=dataset,
        variable="so",
        model=model,
        time_range=time_range,
        chunk_time=12,
    )
    mlotst = load_regular_data(
        dataset=dataset,
        variable="mlotst",
        model=model,
        time_range=time_range,
        chunk_time=12,
    )

    # Ensure temperature metadata reflects potential temperature for TEOS-10 route
    thetao = thetao.assign_attrs(**{**thetao.attrs, "temperature_type": "potential"})

    # Load grid metadata requested in the plan
    volumes = load_volumes_for_cmip6(
        dataset=dataset,
        model=model,
        depth_range=[float(thetao.depth.min()), float(thetao.depth.max())],
        variable="thetao",
    )
    levels, bounds = load_depth_levels(
        dataset=dataset,
        model=model,
        variable="thetao",
    )

    # Use shallowest model level as reference depth
    ref_depth = float(np.asarray(levels)[0])

    # Offline mixed layer depth from monthly thetao and so on native grid
    mld_offline = compute_mixed_layer_depth(
        temp_data=thetao,
        sal_data=so,
        ref_depth=ref_depth,
        method="density",
        density_threshold=0.03,
        temperature_type="potential",
    ).assign_attrs(units="m", long_name=f"{model} offline mixed layer depth")

    results[model] = {}
    diff_maps_native[model] = {}
    pair_fields_native[model] = {}
    scatter_inputs[model] = {}
    extremes[model] = {}

    # Total ocean area from published mlotst mask on native grid
    mlotst_any = mlotst.isel(time=0)
    total_ocean_area = compute_basic_stat(
        mlotst_any.where(np.isfinite(mlotst_any)),
        op="area_weighted_mean",
        dim=["lat", "lon"],
    )
    area_field = xr.ones_like(mlotst_any).where(np.isfinite(mlotst_any))
    total_ocean_area_amount = compute_exceedance_fraction(
        area_field,
        threshold=0.5,
        comparison="greater",
        dim=["lat", "lon"],
        weighting="area",
        reference=0.0,
    )["total_amount"]
    if hasattr(total_ocean_area_amount, "compute"):
        total_ocean_area_amount = total_ocean_area_amount.compute()
    total_ocean_area_amount = float(total_ocean_area_amount)

    for month in target_months:
        month_label = month_names[month]

        # Extract all instances of the target month
        off_m = mld_offline.where(mld_offline["time"].dt.month == month, drop=True)
        pub_m = mlotst.where(mlotst["time"].dt.month == month, drop=True)

        # Strict shared mask: keep only cells valid in both fields for all 30 same-calendar months
        shared_all30 = (off_m.notnull() & pub_m.notnull()).all(dim="time")

        off_m_shared = off_m.where(shared_all30)
        pub_m_shared = pub_m.where(shared_all30)

        # 1971-2000 climatologies as plain mean over the 30 monthly values
        off_clim = compute_basic_stat(off_m_shared, op="mean", dim="time").assign_attrs(units="m")
        pub_clim = compute_basic_stat(pub_m_shared, op="mean", dim="time").assign_attrs(units="m")

        # Native-grid difference map for all reported metrics
        diff = (off_clim - pub_clim).assign_attrs(
            units="m",
            long_name=f"{model} {month_label} offline MLD minus mlotst",
        )

        # Metrics on native grid, area weighted
        metrics = compute_difference_metrics(
            data=off_clim,
            reference=pub_clim,
            dim=["lat", "lon"],
            weighting="area",
        )

        bias = metrics["bias"]
        if hasattr(bias, "compute"):
            bias = bias.compute()
        bias = float(bias)

        rmsd = metrics["rmse"]
        if hasattr(rmsd, "compute"):
            rmsd = rmsd.compute()
        rmsd = float(rmsd)

        corr = analyze_correlation(
            var1=off_clim,
            var2=pub_clim,
            dim=["lat", "lon"],
            method="pearson",
            weighting="area",
        )
        spatial_corr = corr["correlation"]
        if hasattr(spatial_corr, "compute"):
            spatial_corr = spatial_corr.compute()
        spatial_corr = float(spatial_corr)

        exc = compute_exceedance_fraction(
            data=diff,
            threshold=20.0,
            comparison="abs_greater",
            dim=["lat", "lon"],
            weighting="area",
            reference=0.0,
        )
        frac_gt20 = exc["fraction_percent"]
        if hasattr(frac_gt20, "compute"):
            frac_gt20 = frac_gt20.compute()
        frac_gt20 = float(frac_gt20)

        valid_area_amount = exc["total_amount"]
        if hasattr(valid_area_amount, "compute"):
            valid_area_amount = valid_area_amount.compute()
        valid_area_amount = float(valid_area_amount)
        valid_area_pct = 100.0 * valid_area_amount / total_ocean_area_amount

        # Regional concentration of disagreement
        pos_ext = locate_extreme(diff, op="max")
        neg_ext = locate_extreme(diff, op="min")
        abs_ext = locate_extreme(diff, op="abs_max")

        # Materialize all quantities that will be reported or plotted
        off_clim_c = off_clim.compute() if hasattr(off_clim, "compute") else off_clim
        pub_clim_c = pub_clim.compute() if hasattr(pub_clim, "compute") else pub_clim
        diff_c = diff.compute() if hasattr(diff, "compute") else diff
        mask_c = shared_all30.compute() if hasattr(shared_all30, "compute") else shared_all30

        results[model][month] = {
            "month_name": month_label,
            "bias_m": round(bias, 1),
            "rmsd_m": round(rmsd, 1),
            "corr": round(spatial_corr, 3),
            "area_absdiff_gt20_percent": round(frac_gt20, 1),
            "valid_area_percent_of_ocean": round(valid_area_pct, 1),
        }

        diff_maps_native[model][month] = diff_c
        pair_fields_native[model][month] = {
            "offline_mld": off_clim_c,
            "mlotst": pub_clim_c,
            "shared_mask": mask_c,
        }
        scatter_inputs[model][month] = {
            "x_offline": off_clim_c.values.ravel(),
            "y_mlotst": pub_clim_c.values.ravel(),
        }
        extremes[model][month] = {
            "max_positive": pos_ext,
            "max_negative": neg_ext,
            "max_abs": abs_ext,
        }

        # Regrid for display only
        diff_display = regrid_horizontal(diff_c, target="woa23", method="bilinear")
        diff_maps_native[model][f"{month}_display"] = diff_display.compute() if hasattr(diff_display, "compute") else diff_display
        display_maps.append(diff_maps_native[model][f"{month}_display"].values.ravel())

# Ranking criterion: lowest mean of normalized RMSD rank, |bias| rank, highest correlation rank,
# and lowest exceedance-area rank across March and September.
ranking_scores = []
for model in models:
    rmsd_mean = np.mean([results[model][m]["rmsd_m"] for m in target_months])
    abias_mean = np.mean([abs(results[model][m]["bias_m"]) for m in target_months])
    corr_mean = np.mean([results[model][m]["corr"] for m in target_months])
    exc_mean = np.mean([results[model][m]["area_absdiff_gt20_percent"] for m in target_months])
    ranking_scores.append(
        {
            "model": model,
            "rmsd_mean": rmsd_mean,
            "abias_mean": abias_mean,
            "corr_mean": corr_mean,
            "exc_mean": exc_mean,
        }
    )

rank_df = pd.DataFrame(ranking_scores).set_index("model")
rank_df["r_rmsd"] = rank_df["rmsd_mean"].rank(ascending=True, method="min")
rank_df["r_abias"] = rank_df["abias_mean"].rank(ascending=True, method="min")
rank_df["r_corr"] = rank_df["corr_mean"].rank(ascending=False, method="min")
rank_df["r_exc"] = rank_df["exc_mean"].rank(ascending=True, method="min")
rank_df["consistency_score"] = rank_df[["r_rmsd", "r_abias", "r_corr", "r_exc"]].mean(axis=1)
ranking = rank_df.sort_values(["consistency_score", "r_rmsd", "r_corr"]).index.tolist()

# Compact summary table
summary_rows = []
for model in models:
    for month in target_months:
        rr = results[model][month]
        summary_rows.append(
            {
                "model": model,
                "month": rr["month_name"],
                "mean_diff_m": rr["bias_m"],
                "rmsd_m": rr["rmsd_m"],
                "spatial_r": rr["corr"],
                "area_|diff|>20m_%": rr["area_absdiff_gt20_percent"],
                "valid_area_%_of_ocean": rr["valid_area_percent_of_ocean"],
            }
        )
summary_table = pd.DataFrame(summary_rows)

ranking_criterion = (
    "Models ranked by the mean of four native-grid March/September ranks: lower RMSD, "
    "lower absolute bias, higher spatial Pearson correlation, and lower area fraction with "
    "|offline MLD - mlotst| > 20 m."
)

diagnostic_note = (
    "Offline MLD and published mlotst are not expected to agree exactly because published mlotst "
    "is based on the model top level and a displaced-parcel density test at local pressure, and "
    "is typically archived as a monthly mean of instantaneous MLDs, whereas the offline operator "
    "is computed from monthly mean thetao/so profiles using a sigma-threshold crossing depth "
    "interpolated between vertical levels. This commonly makes the offline estimate systematically "
    "shallower and can amplify disagreement in deep winter mixed-layer regions and frontal zones."
)
