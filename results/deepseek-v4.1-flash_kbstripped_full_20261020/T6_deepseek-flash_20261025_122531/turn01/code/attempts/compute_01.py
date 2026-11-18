models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
basins = ["Pacific", "Atlantic", "Indian"]
time_range = ["1971-01", "2000-12"]
depth_range = [0.0, 2000.0]


def _as_float(x):
    if hasattr(x, "compute"):
        return float(x.compute())
    return float(x)


def _tag_grid(out, src, units, long_name=None):
    attrs = dict(src.attrs)
    attrs.update(dict(out.attrs))
    attrs["units"] = units
    if long_name is not None:
        attrs["long_name"] = long_name
    return out.assign_attrs(**attrs)


# ----------------------------------------------------------------------
# WOA23 reference nitrate climatology, mask and vertical grid
# ----------------------------------------------------------------------
woa_no3 = load_climatology_data(
    dataset="woa23",
    variable="nitrate",
    time_slice="annual",
    depth_range=depth_range,
)

woa_landsea = load_landsea_mask()
woa_basin_surface = load_basin_mask(depth_levels="0m")
woa_levels, woa_bounds = load_depth_levels(
    dataset="woa23",
    variable="nitrate",
    time_slice="annual",
    depth_range=depth_range,
)

woa_nc = compute_nutricline_depth(woa_no3, min_depth=0.0, max_depth=1000.0)
woa_nc = _tag_grid(woa_nc, woa_no3, units="m", long_name="Nitracline depth")

woa_global = _as_float(
    compute_basic_stat(woa_nc, op="area_weighted_mean", dim=["lat", "lon"])
)

woa_basin = {}
for b in basins:
    fld = select_basin(
        woa_nc,
        basin=b,
        include_marginal_seas=True,
        depth_dependent=False,
    )
    woa_basin[b] = _as_float(
        compute_basic_stat(fld, op="area_weighted_mean", dim=["lat", "lon"])
    )


# ----------------------------------------------------------------------
# CMIP6 models: 1971-2000 nitrate climatology, nitracline, density, units
# ----------------------------------------------------------------------
model_no3_clim = {}
model_nc = {}
model_rho = {}
model_no3_umolkg = {}
model_global = {}
model_basin = {}

for m in models:
    no3 = load_regular_data(
        dataset="cmip6_historical",
        variable="no3",
        model=m,
        time_range=time_range,
        depth_range=depth_range,
    )
    thetao = load_regular_data(
        dataset="cmip6_historical",
        variable="thetao",
        model=m,
        time_range=time_range,
        depth_range=depth_range,
    )
    so = load_regular_data(
        dataset="cmip6_historical",
        variable="so",
        model=m,
        time_range=time_range,
        depth_range=depth_range,
    )

    no3_clim = temporal_aggregate(no3, freq="mean", op="mean")
    thetao_clim = temporal_aggregate(thetao, freq="mean", op="mean")
    so_clim = temporal_aggregate(so, freq="mean", op="mean")

    nc = compute_nutricline_depth(no3_clim, min_depth=0.0, max_depth=1000.0)
    nc = _tag_grid(nc, no3_clim, units="m", long_name="Nitracline depth")

    rho = compute_potential_density_anomaly(
        thetao_clim,
        so_clim,
        ref_pressure=None,
    )

    no3_umolkg = transform_unit(
        no3_clim,
        from_unit="mol m-3",
        to_unit="umol kg-1",
        density=rho,
    )
    no3_umolkg = _tag_grid(
        no3_umolkg,
        no3_clim,
        units="umol kg-1",
        long_name="Nitrate concentration",
    )

    model_no3_clim[m] = no3_clim
    model_nc[m] = nc
    model_rho[m] = rho
    model_no3_umolkg[m] = no3_umolkg

    model_global[m] = _as_float(
        compute_basic_stat(nc, op="area_weighted_mean", dim=["lat", "lon"])
    )

    model_basin[m] = {}
    for b in basins:
        fld = select_basin(
            nc,
            basin=b,
            include_marginal_seas=True,
            depth_dependent=False,
        )
        model_basin[m][b] = _as_float(
            compute_basic_stat(fld, op="area_weighted_mean", dim=["lat", "lon"])
        )


# ----------------------------------------------------------------------
# Summary tables
# ----------------------------------------------------------------------
summary_rows = []
summary_rows.append(
    {
        "source": "WOA23",
        "global_mean_m": woa_global,
        **{f"{b}_m": woa_basin[b] for b in basins},
        "corr_vs_woa": np.nan,
        "rmse_vs_woa_m": np.nan,
    }
)
for m in models:
    summary_rows.append(
        {
            "source": m,
            "global_mean_m": model_global[m],
            **{f"{b}_m": model_basin[m][b] for b in basins},
            "corr_vs_woa": np.nan,
            "rmse_vs_woa_m": np.nan,
        }
    )
summary_df = pd.DataFrame(summary_rows)

global_diff = {m: model_global[m] - woa_global for m in models}
basin_diff = {
    m: {b: model_basin[m][b] - woa_basin[b] for b in basins}
    for m in models
}

diff_rows = []
for m in models:
    diff_rows.append(
        {
            "source": m,
            "global_diff_m": global_diff[m],
            **{f"{b}_diff_m": basin_diff[m][b] for b in basins},
        }
    )
diff_df = pd.DataFrame(diff_rows)


def _shallow_to_deep(d):
    return sorted(d, key=d.get)


woa_basin_order = _shallow_to_deep(woa_basin)
model_basin_order = {m: _shallow_to_deep(model_basin[m]) for m in models}
model_matches_woa_order = {
    m: model_basin_order[m] == woa_basin_order for m in models
}


# ----------------------------------------------------------------------
# Put nitracline fields on the WOA23 1-degree grid and score on one mask
# ----------------------------------------------------------------------
model_nc_regrid = {}
for m in models:
    model_nc_regrid[m] = regrid_horizontal(
        model_nc[m],
        target="woa23",
        method="conservative",
    )

fields_for_mask = [woa_nc] + [model_nc_regrid[m] for m in models]
shared_mask = build_comparison_mask(
    fields_for_mask,
    ocean_only=True,
    reduce_other_dims="all",
)
n_shared = int(shared_mask.sum().compute())

woa_cmp = woa_nc.where(shared_mask)
woa_cmp = _tag_grid(woa_cmp, woa_nc, units="m", long_name="Nitracline depth")

model_cmp = {}
corr_vals = {}
rmse_vals = {}
bias_vals = {}

for m in models:
    cmp = model_nc_regrid[m].where(shared_mask)
    cmp = _tag_grid(cmp, model_nc_regrid[m], units="m", long_name="Nitracline depth")
    model_cmp[m] = cmp

    corr_dict = analyze_correlation(
        cmp,
        woa_cmp,
        dim=["lat", "lon"],
        method="pearson",
        weighting="area",
    )
    corr_vals[m] = _as_float(corr_dict["correlation"])

    diff_dict = compute_difference_metrics(
        data=cmp,
        reference=woa_cmp,
        dim=["lat", "lon"],
        weighting="auto",
    )
    rmse_vals[m] = _as_float(diff_dict["rmse"])
    bias_vals[m] = _as_float(diff_dict["bias"])

best_model = min(models, key=lambda m: (rmse_vals[m], -corr_vals[m]))

for row in summary_rows:
    if row["source"] in models:
        row["corr_vs_woa"] = corr_vals[row["source"]]
        row["rmse_vs_woa_m"] = rmse_vals[row["source"]]
summary_df = pd.DataFrame(summary_rows)

score_df = pd.DataFrame(
    {
        "source": models,
        "corr_vs_woa": [corr_vals[m] for m in models],
        "rmse_vs_woa_m": [rmse_vals[m] for m in models],
        "bias_vs_woa_m": [bias_vals[m] for m in models],
    }
)


# ----------------------------------------------------------------------
# Equatorial Pacific sections for Figure 2
# ----------------------------------------------------------------------
woa_levels_arr = np.asarray(woa_levels)
woa_lev_1000 = woa_levels_arr[woa_levels_arr <= 1000.0]
deeper_levels = woa_levels_arr[woa_levels_arr > 1000.0]
if deeper_levels.size > 0:
    woa_lev_nc = np.concatenate([woa_lev_1000, [deeper_levels[0]]])
else:
    woa_lev_nc = woa_lev_1000

woa_sec_raw = load_climatology_data(
    dataset="woa23",
    variable="nitrate",
    time_slice="annual",
    lat_range=[-2.0, 2.0],
    lon_range=[120.0, 280.0],
    depth_range=[0.0, 1200.0],
)
woa_sec_lat = compute_basic_stat(woa_sec_raw, op="mean", dim="lat")
woa_sec = regrid_vertical(woa_sec_lat, target=woa_lev_nc, method="linear")
woa_sec_nc = compute_nutricline_depth(woa_sec, min_depth=0.0, max_depth=1000.0)
woa_sec_nc = _tag_grid(woa_sec_nc, woa_sec, units="m", long_name="Nitracline depth")

model_sec = {}
model_sec_nc = {}
for m in models:
    no3_reg = regrid_horizontal(
        model_no3_umolkg[m],
        target="woa23",
        method="conservative",
    )
    no3_reg = no3_reg.assign_coords(lon=(no3_reg.lon % 360.0)).sortby("lon")
    sec_raw = no3_reg.sel(
        lat=slice(-2.0, 2.0),
        lon=slice(120.0, 280.0),
    )
    sec_lat = compute_basic_stat(sec_raw, op="mean", dim="lat")
    sec = regrid_vertical(sec_lat, target=woa_lev_nc, method="linear")

    model_sec[m] = sec
    sec_nc = compute_nutricline_depth(sec, min_depth=0.0, max_depth=1000.0)
    sec_nc = _tag_grid(sec_nc, sec, units="m", long_name="Nitracline depth")
    model_sec_nc[m] = sec_nc
