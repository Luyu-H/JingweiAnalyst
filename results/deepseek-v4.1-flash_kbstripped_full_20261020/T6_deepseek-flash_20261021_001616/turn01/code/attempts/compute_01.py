# Helpers
FILL_LIMIT = 1e10

def mask_fill(da):
    """Remove the CMIP6 1e20 missing-value sentinel before any arithmetic."""
    return da.where(da < FILL_LIMIT)

def with_template_attrs(new, template, **extra):
    """Carry the loader/grid provenance attrs forward while setting new metadata."""
    attrs = dict(template.attrs)
    attrs.update(new.attrs)
    attrs.update(extra)
    return new.assign_attrs(attrs)

# ---------------------------------------------------------------------------
# 1. WOA23 reference nitrate, land-sea mask, basin mask and areas
# ---------------------------------------------------------------------------
woa_no3 = load_climatology_data(
    dataset="woa23",
    variable="nitrate",
    time_slice="annual",
    depth_range=[0, 1200],
)
landsea = load_landsea_mask()
ocean_mask = landsea > 1
woa_no3 = woa_no3.where(ocean_mask)
woa_no3 = with_template_attrs(
    woa_no3, woa_no3,
    units=woa_no3.attrs.get("units", "umol kg-1"),
    long_name="WOA23 nitrate",
)

woa_areas = load_areas(dataset="woa23")
woa_levels, woa_bounds = load_depth_levels(
    dataset="woa23", variable="nitrate", time_slice="annual", depth_range=[0, 1200]
)
basin_mask_surface = load_basin_mask(depth_levels="0m")

# WOA23 nitracline: max vertical gradient searched over the top 1000 m,
# with the deepest gradient using the level just below 1000 m (loaded above).
woa_nc = compute_nutricline_depth(woa_no3, min_depth=0, max_depth=1000)
woa_nc = with_template_attrs(
    woa_nc, woa_no3, units="m", long_name="Nitracline depth"
)
woa_global_mean = float(
    compute_basic_stat(woa_nc, op="area_weighted_mean", dim=None).compute()
)

basin_names = ["Pacific", "Atlantic", "Indian"]
woa_basin_means = {}
for b in basin_names:
    bmasked = select_basin(
        woa_nc, basin=b, include_marginal_seas=True, depth_dependent=False
    )
    woa_basin_means[b] = float(
        compute_basic_stat(bmasked, op="area_weighted_mean", dim=None).compute()
    )

# ---------------------------------------------------------------------------
# 2. CMIP6 historical no3, thetao and so, 1971-2000, native grids
#    Convert no3 mol m-3 -> umol kg-1 using in-situ TEOS-10 density.
# ---------------------------------------------------------------------------
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

model_no3_mean = {}
model_nc = {}
model_global_mean = {}
model_basin_means = {}
model_areas = {}
model_levels = {}

for model in models:
    no3 = load_regular_data(
        dataset="cmip6_historical",
        variable="no3",
        model=model,
        time_range=["1971-01", "2000-12"],
        depth_range=[0, 1200],
    )
    thetao = load_regular_data(
        dataset="cmip6_historical",
        variable="thetao",
        model=model,
        time_range=["1971-01", "2000-12"],
        depth_range=[0, 1200],
    )
    so = load_regular_data(
        dataset="cmip6_historical",
        variable="so",
        model=model,
        time_range=["1971-01", "2000-12"],
        depth_range=[0, 1200],
    )

    no3 = mask_fill(no3)
    thetao = mask_fill(thetao)
    so = mask_fill(so)

    # In-situ density anomaly (sigma), then absolute rho for manual division.
    sigma = compute_potential_density_anomaly(
        temp_data=thetao, sal_data=so, ref_pressure=None
    )
    rho = (sigma + 1000.0).assign_attrs(units="kg m-3", long_name="In-situ density")

    # mol m-3 -> umol m-3, then divide by rho (kg m-3) -> umol kg-1.
    no3_vol = transform_unit(
        data=no3, from_unit="mol m-3", to_unit="umol m-3"
    )
    no3_mass = compute_ratio(data1=no3_vol, data2=rho)
    no3_mass = with_template_attrs(
        no3_mass, no3,
        units="umol kg-1",
        density_assumption="in-situ TEOS-10 density (sigma + 1000)",
    )

    # 1971-2000 mean, day-weighted by month length.
    no3_mean = temporal_aggregate(data=no3_mass, freq="mean", op="mean")
    no3_mean = with_template_attrs(no3_mean, no3_mass)
    model_no3_mean[model] = no3_mean

    # Nitracline on the native grid.
    nc = compute_nutricline_depth(no3_mean, min_depth=0, max_depth=1000)
    nc = with_template_attrs(
        nc, no3_mean, units="m", long_name="Nitracline depth"
    )
    model_nc[model] = nc

    model_areas[model] = load_areas(
        dataset="cmip6_historical", model=model, variable="no3"
    )
    model_levels[model] = load_depth_levels(
        dataset="cmip6_historical", model=model, variable="no3",
        depth_range=[0, 1200]
    )

    # Global and basin means on each model's own native grid.
    model_global_mean[model] = float(
        compute_basic_stat(nc, op="area_weighted_mean", dim=None).compute()
    )
    model_basin_means[model] = {}
    for b in basin_names:
        bmasked = select_basin(
            nc, basin=b, include_marginal_seas=True, depth_dependent=False
        )
        model_basin_means[model][b] = float(
            compute_basic_stat(bmasked, op="area_weighted_mean", dim=None).compute()
        )

# ---------------------------------------------------------------------------
# 3. Regrid nitracline fields to WOA23 1-degree grid, shared mask, metrics
# ---------------------------------------------------------------------------
nc_model_reg = {}
for model in models:
    nc_model_reg[model] = regrid_horizontal(
        data=model_nc[model], target="woa23", method="conservative"
    )

woa_nc_reg = woa_nc  # already on the WOA23 1-degree grid

fields_for_mask = [woa_nc_reg] + [nc_model_reg[m] for m in models]
cmp_mask = build_comparison_mask(fields=fields_for_mask, ocean_only=True)
cmp_mask = cmp_mask.compute()
n_cmp = int(cmp_mask.sum().compute().item())

woa_nc_cmp = woa_nc_reg.where(cmp_mask)
model_nc_cmp = {m: nc_model_reg[m].where(cmp_mask) for m in models}

corr_stats = {}
diff_stats = {}

for model in models:
    corr_dict = analyze_correlation(
        var1=model_nc_cmp[model],
        var2=woa_nc_cmp,
        dim=["lat", "lon"],
        method="pearson",
        weighting="area",
    )
    corr_val = float(corr_dict["correlation"].compute())
    p_val = float(corr_dict["pvalue"].compute())

    metrics = compute_difference_metrics(
        data=model_nc_cmp[model],
        reference=woa_nc_cmp,
        dim=["lat", "lon"],
        weighting="auto",
    )
    rmse_val = float(metrics["rmse"].compute())
    bias_val = float(metrics["bias"].compute())
    mean_model = float(metrics["mean_data"].compute())
    mean_ref = float(metrics["mean_reference"].compute())

    corr_stats[model] = {"correlation": corr_val, "pvalue": p_val}
    diff_stats[model] = {
        "rmse": rmse_val,
        "bias": bias_val,
        "mean_data": mean_model,
        "mean_reference": mean_ref,
    }

best_model_by_rmse = min(models, key=lambda m: diff_stats[m]["rmse"])
best_model_by_corr = max(models, key=lambda m: corr_stats[m]["correlation"])

woa_basin_order = sorted(
    basin_names, key=lambda b: woa_basin_means[b], reverse=True
)
basin_order = {}
basin_order_matches = {}
for model in models:
    basin_order[model] = sorted(
        basin_names, key=lambda b: model_basin_means[model][b], reverse=True
    )
    basin_order_matches[model] = basin_order[model] == woa_basin_order

# ---------------------------------------------------------------------------
# 4. Equatorial Pacific sections (2S-2N, 120E-80W) for Figure 2
# ---------------------------------------------------------------------------
section_fields = {}
section_lon = {}
section_nc = {}

# WOA23 section
woa_eq = woa_no3.sel(lat=slice(-2, 2))
lon_vals = woa_eq["lon"]
woa_eq = woa_eq.where(
    ((lon_vals >= 120) & (lon_vals <= 180))
    | ((lon_vals >= -180) & (lon_vals <= -80))
)
woa_section = compute_basic_stat(woa_eq, op="area_weighted_mean", dim="lat")
woa_section = with_template_attrs(
    woa_section, woa_no3, units=woa_no3.attrs.get("units", "umol kg-1")
)
section_fields["WOA23"] = woa_section
section_lon["WOA23"] = woa_section["lon"]

woa_section_nc = compute_nutricline_depth(
    woa_section, min_depth=0, max_depth=1000
)
section_nc["WOA23"] = with_template_attrs(
    woa_section_nc, woa_section, units="m", long_name="Nitracline depth"
)

# Model sections
for model in models:
    da = model_no3_mean[model]
    lat_name = da.attrs.get("lat_coord", "lat")
    lon_name = da.attrs.get("lon_coord", "lon")
    lat2d = da.coords[lat_name]
    lon2d = da.coords[lon_name]

    lat_mask = (lat2d >= -2) & (lat2d <= 2)
    lon_mask = (lon2d >= 120) | (lon2d <= -80)
    eq_mask = lat_mask & lon_mask

    da_eq = da.where(eq_mask)

    if len(lat2d.dims) > 1:
        lat_dim = lat2d.dims[0]
        lon_1d = lon2d.where(lat_mask).mean(dim=lat_dim)
    else:
        lat_dim = lat_name
        lon_1d = lon2d

    section = compute_basic_stat(
        da_eq, op="area_weighted_mean", dim=lat_dim
    )
    section = with_template_attrs(section, da, units="umol kg-1")
    section_fields[model] = section
    section_lon[model] = lon_1d

    sec_nc = compute_nutricline_depth(section, min_depth=0, max_depth=1000)
    section_nc[model] = with_template_attrs(
        sec_nc, section, units="m", long_name="Nitracline depth"
    )

# ---------------------------------------------------------------------------
# 5. Report
# ---------------------------------------------------------------------------
print("Global mean nitracline depth (m):")
print("WOA23:", round(woa_global_mean))
for m in models:
    print(
        m,
        round(model_global_mean[m]),
        "diff from WOA23:",
        round(model_global_mean[m] - woa_global_mean),
    )

print("\nBasin mean nitracline depths (m):")
print("WOA23:", {b: round(woa_basin_means[b]) for b in basin_names})
for m in models:
    print(m, {b: round(model_basin_means[m][b]) for b in basin_names})
    print(
        "  diffs from WOA23:",
        {b: round(model_basin_means[m][b] - woa_basin_means[b]) for b in basin_names},
    )

print("\nShared-mask metrics on WOA23 1-degree grid:")
print("comparison cells:", n_cmp)
for m in models:
    print(
        m,
        "r =", round(corr_stats[m]["correlation"], 3),
        "RMSE =", round(diff_stats[m]["rmse"], 1),
        "bias =", round(diff_stats[m]["bias"], 1),
        "mean model =", round(diff_stats[m]["mean_data"], 1),
        "mean WOA23 =", round(diff_stats[m]["mean_reference"], 1),
    )
print("Best by RMSE:", best_model_by_rmse)
print("Best by correlation:", best_model_by_corr)

print("\nBasin ordering (deep to shallow):")
print("WOA23:", woa_basin_order)
for m in models:
    print(m, basin_order[m], "matches WOA23:", basin_order_matches[m])
