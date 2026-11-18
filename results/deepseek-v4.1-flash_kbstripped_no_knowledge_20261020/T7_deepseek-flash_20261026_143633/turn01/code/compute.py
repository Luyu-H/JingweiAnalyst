# Helper functions
def _merge_attrs(da, ref_attrs):
    da.attrs = {**ref_attrs, **da.attrs}
    return da


def _make_ratio(no3, po4, valid):
    r = compute_ratio(no3.where(valid), po4.where(valid))
    r.attrs = {**no3.attrs, **r.attrs}
    r = r.assign_attrs(units="1", long_name="N:P ratio")
    return r


def _points_from_fields(no3, po4, valid, max_points=20000, seed=12345):
    n3 = np.asarray(no3.where(valid).values).ravel()
    p4 = np.asarray(po4.where(valid).values).ravel()
    m = np.isfinite(n3) & np.isfinite(p4)
    x = p4[m].astype(float)
    y = n3[m].astype(float)
    if x.size > max_points:
        rng = np.random.default_rng(seed)
        idx = rng.choice(x.size, size=max_points, replace=False)
        x = x[idx]
        y = y[idx]
    return x, y


def _ols_direct(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m].astype(float)
    y = y[m].astype(float)
    a, b = np.polyfit(x, y, 1)
    yhat = a * x + b
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return float(a), float(b), float(r2)


def process_source(name, no3, po4, valid):
    ratio = _make_ratio(no3, po4, valid)

    # Headline mean: mean of cell-by-cell N:P
    mean_ratio = float(
        compute_basic_stat(ratio, op="volume_weighted_mean", dim=None)
    )

    # Ratio of mean nitrate to mean phosphate
    mean_no3_da = compute_basic_stat(
        no3.where(valid), op="volume_weighted_mean", dim=None
    )
    mean_po4_da = compute_basic_stat(
        po4.where(valid), op="volume_weighted_mean", dim=None
    )
    ratio_means = float(compute_ratio(mean_no3_da, mean_po4_da))

    # OLS over all valid cells, each cell weighted equally
    x_full, y_full = _points_from_fields(no3, po4, valid, max_points=10**12)
    a_direct, b_direct, r2_direct = _ols_direct(x_full, y_full)

    # Use the provided regression tool for diagnostics / validation,
    # but report the unambiguous direct OLS values.
    fit_result = fit_multivariate_relation(
        target=no3.where(valid),
        predictors=[po4.where(valid)],
        method="linear",
        validation="random",
        test_size=0.2,
        standardize=False,
    )
    r2_fit = None
    if isinstance(fit_result, dict) and "metrics" in fit_result and "r2" in fit_result["metrics"]:
        r2_fit = float(fit_result["metrics"]["r2"])
    del fit_result
    r2_report = r2_fit if r2_fit is not None else r2_direct

    # Scatter points for the figure (subsampled for plotting only)
    scatter_x, scatter_y = _points_from_fields(no3, po4, valid, max_points=20000)

    # Volume fractions
    frac_14_18 = float(
        compute_exceedance_fraction(
            ratio,
            threshold=2.0,
            comparison="abs_less",
            reference=16.0,
            dim=None,
            weighting="volume",
        )["fraction_percent"]
    )
    frac_lt12_global = float(
        compute_exceedance_fraction(
            ratio,
            threshold=12.0,
            comparison="less",
            reference=0.0,
            dim=None,
            weighting="volume",
        )["fraction_percent"]
    )

    # Basin statistics
    basin_means = {}
    frac_lt12_basin = {}
    for basin in ["Pacific", "Atlantic", "Indian"]:
        basin_ratio = select_basin(
            ratio, basin=basin, include_marginal_seas=True
        )
        basin_means[basin] = float(
            compute_basic_stat(basin_ratio, op="volume_weighted_mean", dim=None)
        )
        if basin in ["Pacific", "Indian"]:
            frac_lt12_basin[basin] = float(
                compute_exceedance_fraction(
                    basin_ratio,
                    threshold=12.0,
                    comparison="less",
                    reference=0.0,
                    dim=None,
                    weighting="volume",
                )["fraction_percent"]
            )

    # Map of 200-1000 m mean N:P minus 16
    map_ratio = compute_basic_stat(
        ratio, op="volume_weighted_mean", dim="depth"
    )
    map_minus16 = map_ratio - 16.0
    map_minus16.attrs = {
        **map_ratio.attrs,
        "units": "1",
        "long_name": "N:P minus 16",
    }

    if mean_ratio > 16.005:
        classification = "above Redfield"
    elif mean_ratio < 15.995:
        classification = "below Redfield"
    else:
        classification = "at Redfield"

    return {
        "mean_ratio": mean_ratio,
        "ratio_means": ratio_means,
        "slope": a_direct,
        "intercept": b_direct,
        "r2": r2_report,
        "r2_direct": r2_direct,
        "frac_14_18": frac_14_18,
        "frac_lt12_global": frac_lt12_global,
        "frac_lt12_pacific": frac_lt12_basin.get("Pacific", np.nan),
        "frac_lt12_indian": frac_lt12_basin.get("Indian", np.nan),
        "basin_means": basin_means,
        "classification": classification,
        "map_minus16": map_minus16,
        "scatter_x": scatter_x,
        "scatter_y": scatter_y,
    }


sources_order = ["WOA23", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
results = {}

# ----------------------------------------------------------------------
# WOA23 observed climatology
# ----------------------------------------------------------------------
woa_no3 = load_climatology_data(
    dataset="woa23",
    variable="nitrate",
    time_slice="annual",
    depth_range=[200, 1000],
).compute()
woa_po4 = load_climatology_data(
    dataset="woa23",
    variable="phosphate",
    time_slice="annual",
    depth_range=[200, 1000],
).compute()

valid_woa = woa_no3.notnull() & woa_po4.notnull()
results["WOA23"] = process_source("WOA23", woa_no3, woa_po4, valid_woa)
results["WOA23"]["map_display"] = results["WOA23"]["map_minus16"]

# ----------------------------------------------------------------------
# CMIP6 historical models
# ----------------------------------------------------------------------
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

for model in models:
    try:
        no3_m = load_regular_data(
            dataset="cmip6_historical",
            variable="no3",
            model=model,
            time_range=["1971-01", "2000-12"],
            depth_range=[200, 1000],
        )
        po4_m = load_regular_data(
            dataset="cmip6_historical",
            variable="po4",
            model=model,
            time_range=["1971-01", "2000-12"],
            depth_range=[200, 1000],
        )
    except Exception as exc:
        results[model] = {"missing": str(exc)}
        continue

    # Day-weighted 1971-2000 climatological mean
    no3_mean = temporal_aggregate(no3_m, freq="mean", op="mean")
    po4_mean = temporal_aggregate(po4_m, freq="mean", op="mean")
    no3_mean = _merge_attrs(no3_mean, no3_m.attrs)
    po4_mean = _merge_attrs(po4_mean, po4_m.attrs)

    # Convert mol m-3 to umol kg-1 at 1025 kg m-3 = 1.025 kg L-1
    no3_umol = transform_unit(
        no3_mean,
        from_unit=no3_mean.attrs.get("units", "mol m-3"),
        to_unit="umol kg-1",
        density=1.025,
    )
    po4_umol = transform_unit(
        po4_mean,
        from_unit=po4_mean.attrs.get("units", "mol m-3"),
        to_unit="umol kg-1",
        density=1.025,
    )
    no3_umol = _merge_attrs(no3_umol, no3_mean.attrs)
    po4_umol = _merge_attrs(po4_umol, po4_mean.attrs)
    no3_umol = no3_umol.compute()
    po4_umol = po4_umol.compute()

    valid_m = no3_umol.notnull() & po4_umol.notnull()
    res = process_source(model, no3_umol, po4_umol, valid_m)

    # Regrid native model map to WOA23 regular grid for display only
    res["map_display"] = regrid_horizontal(
        res["map_minus16"], target="woa23", method="bilinear"
    )
    results[model] = res

# ----------------------------------------------------------------------
# Print summary
# ----------------------------------------------------------------------
print("N:P stoichiometry of intermediate water (200-1000 m)")
print("=" * 62)
for src in sources_order:
    if src not in results:
        continue
    r = results[src]
    if "missing" in r:
        print(f"\n{src}: MISSING ({r['missing']})")
        continue

    print(f"\n{src}:")
    print(f"  Headline mean N:P (mean of cell-by-cell ratio) = {r['mean_ratio']:.2f}")
    print(f"  Ratio of mean nitrate to mean phosphate         = {r['ratio_means']:.2f}")
    print(
        f"  OLS: NO3 = {r['slope']:.2f} * PO4 + {r['intercept']:.2f}, "
        f"R2 = {r['r2']:.3f}"
    )
    print(f"  Volume fraction 14 < N:P < 18 = {r['frac_14_18']:.1f}%")
    print(
        f"  Volume fraction N:P < 12: global = {r['frac_lt12_global']:.1f}%, "
        f"Pacific = {r['frac_lt12_pacific']:.1f}%, "
        f"Indian = {r['frac_lt12_indian']:.1f}%"
    )
    print(
        f"  Basin mean N:P: Pacific = {r['basin_means']['Pacific']:.2f}, "
        f"Atlantic = {r['basin_means']['Atlantic']:.2f}, "
        f"Indian = {r['basin_means']['Indian']:.2f}"
    )
    print(f"  Classification vs Redfield 16: {r['classification']}")
