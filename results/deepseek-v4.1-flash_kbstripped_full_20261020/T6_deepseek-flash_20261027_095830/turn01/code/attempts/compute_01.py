# Helper for materialising scalar results
def _scalar(x):
    return float(x.compute()) if hasattr(x, "compute") else float(x)

models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
basins = ["Pacific", "Atlantic", "Indian"]

# ------------------------------------------------------------------
# WOA23 observed nitrate climatology and nitracline
# ------------------------------------------------------------------
woa_nit_raw = load_climatology_data(
    dataset="woa23",
    variable="nitrate",
    time_slice="annual",
    depth_range=[0, 1200],
)
woa_attrs = woa_nit_raw.attrs.copy()

landsea = load_landsea_mask()
woa_nit = woa_nit_raw.where(landsea > 1)   # ocean is mask > 1
woa_nit.attrs.update(woa_attrs)
woa_nit = woa_nit.assign_attrs(units=woa_nit.attrs.get("units", "umol kg-1"))

woa_levels, woa_bounds = load_depth_levels(dataset="woa23", variable="nitrate")

woa_nit_depth = compute_nutricline_depth(
    data=woa_nit,
    min_depth=0,
    max_depth=1000,
)
woa_nit_depth.attrs.update(woa_attrs)
woa_nit_depth = woa_nit_depth.assign_attrs(units="m").compute()

woa_global = _scalar(
    compute_basic_stat(
        data=woa_nit_depth,
        op="area_weighted_mean",
        dim=["lat", "lon"],
    )
)

woa_basin = {}
for b in basins:
    sel = select_basin(
        data=woa_nit_depth,
        basin=b,
        include_marginal_seas=True,
        depth_dependent=False,
    )
    woa_basin[b] = _scalar(
        compute_basic_stat(
            data=sel,
            op="area_weighted_mean",
            dim=["lat", "lon"],
        )
    )

# ------------------------------------------------------------------
# CMIP6 historical nitrate, 1971-2000 day-weighted mean, nitracline
# ------------------------------------------------------------------
model_global = {}
model_basin = {}
model_nit_depth = {}
model_regrid = {}
model_attrs = {}

for m in models:
    no3_raw = load_regular_data(
        dataset="cmip6_historical",
        variable="no3",
        model=m,
        time_range=["1971-01", "2000-12"],
        depth_range=[0, 1200],
    )
    model_attrs[m] = no3_raw.attrs.copy()

    no3_mean = temporal_aggregate(no3_raw, freq="mean")
    no3_mean.attrs.update(model_attrs[m])

    nit = compute_nutricline_depth(
        data=no3_mean,
        min_depth=0,
        max_depth=1000,
    )
    nit.attrs.update(model_attrs[m])
    nit = nit.assign_attrs(units="m").compute()
    model_nit_depth[m] = nit

    model_global[m] = _scalar(
        compute_basic_stat(
            data=nit,
            op="area_weighted_mean",
            dim=["lat", "lon"],
        )
    )

    model_basin[m] = {}
    for b in basins:
        sel = select_basin(
            data=nit,
            basin=b,
            include_marginal_seas=True,
            depth_dependent=False,
        )
        model_basin[m][b] = _scalar(
            compute_basic_stat(
                data=sel,
                op="area_weighted_mean",
                dim=["lat", "lon"],
            )
        )

    # Regrid only for map comparison and Figure 1. Reported means stay native.
    rg = regrid_horizontal(
        data=nit,
        target="woa23",
        method="conservative",
    )
    model_regrid[m] = rg.compute()

# ------------------------------------------------------------------
# Shared WOA 1-degree mask for map comparison
# ------------------------------------------------------------------
common_mask = build_comparison_mask(
    fields=[
        woa_nit_depth,
        model_regrid["CESM2"],
        model_regrid["GFDL-ESM4"],
        model_regrid["MPI-ESM1-2-LR"],
    ],
    ocean_only=True,
).compute()
n_common = int(common_mask.sum().values)

woa_map = woa_nit_depth.where(common_mask)
model_map = {m: model_regrid[m].where(common_mask) for m in models}

# ------------------------------------------------------------------
# Spatial correlation and RMSE on the common WOA 1-degree mask
# ------------------------------------------------------------------
model_corr = {}
model_rmse = {}
model_bias = {}

for m in models:
    corr = analyze_correlation(
        var1=model_map[m],
        var2=woa_map,
        dim=["lat", "lon"],
        method="pearson",
        weighting="area",
    )
    model_corr[m] = _scalar(corr["correlation"])

    met = compute_difference_metrics(
        data=model_map[m],
        reference=woa_map,
        dim=["lat", "lon"],
        weighting="auto",
    )
    model_rmse[m] = _scalar(met["rmse"])
    model_bias[m] = _scalar(met["bias"])

global_diff = {m: model_global[m] - woa_global for m in models}
basin_diff = {
    b: {m: model_basin[m][b] - woa_basin[b] for m in models}
    for b in basins
}

woa_basin_order = sorted(basins, key=lambda b: woa_basin[b], reverse=True)
model_basin_order = {
    m: sorted(basins, key=lambda b: model_basin[m][b], reverse=True)
    for m in models
}
reproduces_basin_order = {
    m: model_basin_order[m] == woa_basin_order for m in models
}
best_model = min(models, key=lambda m: (model_rmse[m], abs(global_diff[m])))

# ------------------------------------------------------------------
# Figure 2 preparation: equatorial Pacific nitrate section
# ------------------------------------------------------------------
eq_lat_range = [-2, 2]
eq_lon_range = [120, 280]
zmax = 1200

# WOA23 section
woa_eq_raw = load_climatology_data(
    dataset="woa23",
    variable="nitrate",
    time_slice="annual",
    depth_range=[0, zmax],
    lat_range=eq_lat_range,
    lon_range=eq_lon_range,
)
woa_eq_attrs = woa_eq_raw.attrs.copy()
woa_eq = woa_eq_raw

lat_name_woa = woa_eq.attrs.get("lat_coord", "lat")
lon_name_woa = woa_eq.attrs.get("lon_coord", "lon")
mask_woa = (
    (woa_eq[lat_name_woa] >= eq_lat_range[0])
    & (woa_eq[lat_name_woa] <= eq_lat_range[1])
    & (woa_eq[lon_name_woa] >= eq_lon_range[0])
    & (woa_eq[lon_name_woa] <= eq_lon_range[1])
)
woa_eq = woa_eq.where(mask_woa)
woa_eq.attrs.update(woa_eq_attrs)
woa_eq = woa_eq.assign_attrs(units=woa_eq.attrs.get("units", "umol kg-1"))

woa_sec = compute_basic_stat(
    data=woa_eq,
    op="area_weighted_mean",
    dim=["lat"],
)
woa_sec.attrs.update(woa_eq_attrs)
woa_sec = woa_sec.assign_attrs(units=woa_eq.attrs.get("units", "umol kg-1"))
woa_sec = woa_sec.compute()

woa_nit_sec = compute_nutricline_depth(
    data=woa_sec,
    min_depth=0,
    max_depth=1000,
)
woa_nit_sec.attrs.update(woa_sec.attrs)
woa_nit_sec = woa_nit_sec.assign_attrs(units="m").compute()

sections = {"WOA23": woa_sec}
nit_sections = {"WOA23": woa_nit_sec}
density_assumptions = {}

# CMIP6 sections, converted to umol kg-1 using in-situ density
for m in models:
    no3_eq_raw = load_regular_data(
        dataset="cmip6_historical",
        variable="no3",
        model=m,
        time_range=["1971-01", "2000-12"],
        depth_range=[0, zmax],
        lat_range=eq_lat_range,
        lon_range=eq_lon_range,
    )
    thetao_eq_raw = load_regular_data(
        dataset="cmip6_historical",
        variable="thetao",
        model=m,
        time_range=["1971-01", "2000-12"],
        depth_range=[0, zmax],
        lat_range=eq_lat_range,
        lon_range=eq_lon_range,
    )
    so_eq_raw = load_regular_data(
        dataset="cmip6_historical",
        variable="so",
        model=m,
        time_range=["1971-01", "2000-12"],
        depth_range=[0, zmax],
        lat_range=eq_lat_range,
        lon_range=eq_lon_range,
    )

    no3_eq_attrs = no3_eq_raw.attrs.copy()

    no3_eq_mean = temporal_aggregate(no3_eq_raw, freq="mean")
    thetao_eq_mean = temporal_aggregate(thetao_eq_raw, freq="mean")
    so_eq_mean = temporal_aggregate(so_eq_raw, freq="mean")
    no3_eq_mean.attrs.update(no3_eq_attrs)

    sigma_eq = compute_potential_density_anomaly(
        temp_data=thetao_eq_mean,
        sal_data=so_eq_mean,
        ref_pressure=None,
    )

    conc_eq = transform_unit(
        data=no3_eq_mean,
        from_unit="mol m-3",
        to_unit="umol/kg",
        density=sigma_eq,
    )
    density_assumptions[m] = conc_eq.attrs.get("density_assumption", None)

    lat_name = conc_eq.attrs.get("lat_coord", "lat")
    lon_name = conc_eq.attrs.get("lon_coord", "lon")
    lat2 = conc_eq[lat_name]
    lon2 = conc_eq[lon_name]

    mask_eq = (
        (lat2 >= eq_lat_range[0])
        & (lat2 <= eq_lat_range[1])
        & (lon2 >= eq_lon_range[0])
        & (lon2 <= eq_lon_range[1])
    )
    conc_eq = conc_eq.where(mask_eq)

    lat_dims = lat2.dims
    ydim = lat_dims[0] if len(lat_dims) > 0 else "lat"

    lon_rep = lon2.where(mask_eq).mean(dim=ydim)

    sec = compute_basic_stat(
        data=conc_eq,
        op="area_weighted_mean",
        dim=["lat"],
    )
    sec = sec.assign_coords({lon_name: lon_rep})

    sec_reg = regrid_vertical(
        data=sec,
        target=woa_sec,
        method="linear",
    )
    sections[m] = sec_reg.compute()

    nit_sec = compute_nutricline_depth(
        data=sec,
        min_depth=0,
        max_depth=1000,
    )
    nit_sec = nit_sec.assign_coords({lon_name: sec[lon_name]})
    nit_sections[m] = nit_sec.compute()

# ------------------------------------------------------------------
# Report the key numbers
# ------------------------------------------------------------------
print("Global mean nitracline depth (m, nearest m):")
print(f"  WOA23: {int(round(woa_global))}")
for m in models:
    print(
        f"  {m}: {int(round(model_global[m]))}  "
        f"diff(model-WOA): {int(round(global_diff[m]))}"
    )

print("\nBasin mean nitracline depth (m, nearest m):")
for b in basins:
    print(f"  {b}: WOA23 {int(round(woa_basin[b]))}")
    for m in models:
        print(
            f"    {m}: {int(round(model_basin[m][b]))}  "
            f"diff: {int(round(basin_diff[b][m]))}"
        )

print(f"\nCommon-mask cells: {n_common}")
print("\nSpatial metrics vs WOA23 on the common WOA 1-degree mask:")
for m in models:
    print(
        f"  {m}: r={model_corr[m]:.3f}, "
        f"RMSE={model_rmse[m]:.1f} m, bias={model_bias[m]:.1f} m"
    )

print("\nWOA23 basin order (deepest first):", woa_basin_order)
for m in models:
    print(
        f"  {m} order: {model_basin_order[m]}, "
        f"reproduces={reproduces_basin_order[m]}"
    )

print("\nBest model by RMSE (tie broken by global bias):", best_model)
print("Density assumptions used for CMIP6 nitrate conversion:", density_assumptions)
