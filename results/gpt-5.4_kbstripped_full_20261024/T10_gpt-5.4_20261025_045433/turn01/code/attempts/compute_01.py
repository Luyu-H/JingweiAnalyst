# Load observational climatology and masks
woa_o2 = load_climatology_data(
    dataset="woa23",
    variable="oxygen",
    period="1971-2000",
    time_slice="annual",
)
woa_aou = load_climatology_data(
    dataset="woa23",
    variable="apparent_oxygen_utilization",
    period="1971-2000",
    time_slice="annual",
)
woa_o2sat = load_climatology_data(
    dataset="woa23",
    variable="percent_oxygen_saturation",
    period="1971-2000",
    time_slice="annual",
)
landsea_01 = load_landsea_mask()
basinmask_01 = load_basin_mask()

# Helper functions
def _infer_horiz_dims(da):
    lat_name = da.attrs.get("lat_coord", "lat")
    lon_name = da.attrs.get("lon_coord", "lon")
    if lat_name in da.dims and lon_name in da.dims:
        return [lat_name, lon_name]
    if "lat" in da.dims and "lon" in da.dims:
        return ["lat", "lon"]
    candidates = [d for d in da.dims if d != "depth" and d != "time"]
    if len(candidates) >= 2:
        return candidates[-2:]
    raise ValueError(f"Could not infer horizontal dims from {da.dims}")

def _sum_over_all_dims(da):
    dims = list(da.dims)
    return compute_basic_stat(da, op="sum", dim=dims) if len(dims) > 0 else da

def _make_native_thickness_from_volume(volumes, areacello):
    horiz_dims = _infer_horiz_dims(volumes)
    area_b = areacello.broadcast_like(volumes.isel(depth=0, drop=True))
    thick = volumes / area_b
    thick = thick.assign_attrs(units="m", long_name="layer thickness")
    return thick, horiz_dims

def _omz_metrics_from_o2_mass_basis(o2_umolkg, volumes, areacello, threshold):
    omz_thickness, omz_mask = identify_oxygen_minimum_zone(o2_umolkg, threshold=threshold)
    omz_volume_cells = xr.where(omz_mask, volumes, 0.0).assign_attrs(units="m3")
    total_volume = float(_sum_over_all_dims(omz_volume_cells).compute())
    extreme = locate_extreme(omz_thickness.compute(), op="max")
    return {
        "threshold": threshold,
        "omz_thickness": omz_thickness,
        "omz_mask": omz_mask,
        "omz_volume_cells": omz_volume_cells,
        "total_volume_m3": total_volume,
        "max_thickness_m": float(extreme["value"]),
        "max_coords": extreme["coordinates"],
        "extreme": extreme,
    }

def _global_inventory_moles(o2_native):
    content = convert_concentration_to_content(o2_native)
    inventory = float(_sum_over_all_dims(content).compute())
    return content, inventory

def _period_mean_model_fields(dataset, model, time_range):
    o2 = load_regular_data(
        dataset=dataset, model=model, variable="o2", time_range=time_range, chunk_time=12
    )
    thetao = load_regular_data(
        dataset=dataset, model=model, variable="thetao", time_range=time_range, chunk_time=12
    )
    so = load_regular_data(
        dataset=dataset, model=model, variable="so", time_range=time_range, chunk_time=12
    )
    o2_mean = temporal_aggregate(o2, freq="mean", op="mean")
    thetao_mean = temporal_aggregate(thetao, freq="mean", op="mean")
    so_mean = temporal_aggregate(so, freq="mean", op="mean")
    return o2, thetao, so, o2_mean, thetao_mean, so_mean

def _load_model_grid(dataset, model):
    # Native cell volumes
    vols = load_volumes_for_cmip6(dataset=dataset, model=model, variable="thetao")
    vols_mean = vols
    if "time" in vols.dims:
        vols_mean = temporal_aggregate(vols, freq="mean", op="mean")
    # Horizontal cell area and bathymetric depth
    areacello = load_regular_data(dataset=dataset, model=model, variable="areacello")
    deptho = load_regular_data(dataset=dataset, model=model, variable="deptho")
    return vols_mean, areacello, deptho

# WOA23 layer thickness from depth coordinate and basin/land masks on native 1-degree grid
woa_depth = woa_o2["depth"].values.astype(float)
woa_edges = np.empty(len(woa_depth) + 1, dtype=float)
woa_edges[1:-1] = 0.5 * (woa_depth[:-1] + woa_depth[1:])
woa_edges[0] = 0.0
woa_edges[-1] = woa_depth[-1] + 0.5 * (woa_depth[-1] - woa_depth[-2])
woa_dz = xr.DataArray(
    np.diff(woa_edges),
    coords={"depth": woa_o2.depth},
    dims=("depth",),
    attrs={"units": "m", "long_name": "WOA23 standard layer thickness"},
)
woa_bottom_index = (landsea_01 - 1).clip(min=0, max=woa_o2.sizes["depth"]).fillna(0).astype(int)
woa_level_index = xr.DataArray(
    np.arange(1, woa_o2.sizes["depth"] + 1),
    coords={"depth": woa_o2.depth},
    dims=("depth",),
)
woa_ocean3d = woa_level_index <= woa_bottom_index
R = 6371000.0
dlat = np.deg2rad(1.0)
dlon = np.deg2rad(1.0)
lat_rad = np.deg2rad(woa_o2["lat"])
woa_area_2d = xr.DataArray(
    (R**2) * dlat * dlon * np.cos(lat_rad),
    coords={"lat": woa_o2.lat},
    dims=("lat",),
).broadcast_like(woa_o2.isel(depth=0, drop=True))
woa_area_2d = woa_area_2d.where(landsea_01 > 1)
woa_volumes = (woa_dz * woa_area_2d).where(woa_ocean3d)
woa_volumes = woa_volumes.assign_attrs(units="m3", long_name="WOA23 cell volume approximation")

# WOA OMZ diagnostics
woa_omz20 = _omz_metrics_from_o2_mass_basis(woa_o2, woa_volumes, woa_area_2d, threshold=20.0)
woa_omz60 = _omz_metrics_from_o2_mass_basis(woa_o2, woa_volumes, woa_area_2d, threshold=60.0)

# WOA basin volumes at 60 umol kg-1
woa_basin_results_60 = {}
for basin in ["Pacific", "Indian", "Atlantic"]:
    mask_basin = select_basin(
        woa_omz60["omz_volume_cells"],
        basin=basin,
        include_marginal_seas=True,
        depth_dependent=False,
    )
    basin_vol = float(_sum_over_all_dims(mask_basin).compute())
    woa_basin_results_60[basin] = basin_vol

# WOA OMZ AOU and O2 saturation means inside OMZ60, volume-weighted
woa_aou_in_omz60 = woa_aou.where(woa_omz60["omz_mask"])
woa_o2sat_in_omz60 = woa_o2sat.where(woa_omz60["omz_mask"])
woa_aou_mean_omz60 = float(
    compute_basic_stat(woa_aou_in_omz60, op="volume_weighted_mean", dim=["depth", "lat", "lon"]).compute()
)
woa_o2sat_mean_omz60 = float(
    compute_basic_stat(woa_o2sat_in_omz60, op="volume_weighted_mean", dim=["depth", "lat", "lon"]).compute()
)

# Model/experiment definitions
cases = {
    ("cmip6_historical", "GFDL-ESM4"): ["1971-01", "2000-12"],
    ("cmip6_historical", "MPI-ESM1-2-LR"): ["1971-01", "2000-12"],
    ("cmip6_ssp245", "CESM2"): ["2071-01", "2100-12"],
    ("cmip6_ssp245", "GFDL-ESM4"): ["2071-01", "2100-12"],
    ("cmip6_ssp245", "MPI-ESM1-2-LR"): ["2071-01", "2100-12"],
}

native_results = {}
period_means = {}

for (dataset, model), time_range in cases.items():
    o2_native, thetao_native, so_native, o2_mean_native, thetao_mean_native, so_mean_native = _period_mean_model_fields(
        dataset, model, time_range
    )
    vols_native, areacello_native, deptho_native = _load_model_grid(dataset, model)
    sigma_insitu = compute_potential_density_anomaly(
        temp_data=thetao_mean_native,
        sal_data=so_mean_native,
        ref_pressure=None,
        temperature_type="auto",
    )
    o2_mean_umolkg = transform_unit(
        o2_mean_native,
        from_unit=str(o2_mean_native.attrs.get("units", "mol m-3")),
        to_unit="umol/kg",
        density=sigma_insitu,
    )
    omz20 = _omz_metrics_from_o2_mass_basis(o2_mean_umolkg, vols_native, areacello_native, threshold=20.0)
    omz60 = _omz_metrics_from_o2_mass_basis(o2_mean_umolkg, vols_native, areacello_native, threshold=60.0)
    o2_content, o2_inventory = _global_inventory_moles(o2_mean_native)

    native_results[(dataset, model)] = {
        "time_range": time_range,
        "o2_native": o2_native,
        "thetao_native": thetao_native,
        "so_native": so_native,
        "o2_mean_native": o2_mean_native,
        "thetao_mean_native": thetao_mean_native,
        "so_mean_native": so_mean_native,
        "sigma_insitu": sigma_insitu,
        "o2_mean_umolkg": o2_mean_umolkg,
        "volumes": vols_native,
        "areacello": areacello_native,
        "deptho": deptho_native,
        "inventory_content": o2_content,
        "inventory_moles": o2_inventory,
        "omz20": omz20,
        "omz60": omz60,
    }
    period_means[(dataset, model)] = o2_mean_umolkg

# Historical-to-ssp245 percentage changes for the two models with both experiments
scenario_change_pct = {}
for model in ["GFDL-ESM4", "MPI-ESM1-2-LR"]:
    hist20 = native_results[("cmip6_historical", model)]["omz20"]["total_volume_m3"]
    ssp20 = native_results[("cmip6_ssp245", model)]["omz20"]["total_volume_m3"]
    hist60 = native_results[("cmip6_historical", model)]["omz60"]["total_volume_m3"]
    ssp60 = native_results[("cmip6_ssp245", model)]["omz60"]["total_volume_m3"]
    scenario_change_pct[model] = {
        "lt20_pct": 100.0 * (ssp20 - hist20) / hist20,
        "lt60_pct": 100.0 * (ssp60 - hist60) / hist60,
    }

# Regrid historical OMZ thickness maps to WOA grid for display and basin-comparison support
gfdl_hist_omz60_thick_rg = regrid_horizontal(
    native_results[("cmip6_historical", "GFDL-ESM4")]["omz60"]["omz_thickness"],
    target="woa23",
    method="bilinear",
    skipna=True,
    na_thres=0.5,
)
mpi_hist_omz60_thick_rg = regrid_horizontal(
    native_results[("cmip6_historical", "MPI-ESM1-2-LR")]["omz60"]["omz_thickness"],
    target="woa23",
    method="bilinear",
    skipna=True,
    na_thres=0.5,
)
woa_omz60_thickness = woa_omz60["omz_thickness"]

# Basin OMZ volumes at 60 umol kg-1 for the two historical models, using WOA basin definitions after regridding
# Construct regridded support volumes on the WOA grid from native OMZ mask * native volume, then basin-mask on WOA grid
gfdl_hist_omz60_vol_rg = regrid_horizontal(
    native_results[("cmip6_historical", "GFDL-ESM4")]["omz60"]["omz_volume_cells"],
    target="woa23",
    method="conservative",
    skipna=True,
    na_thres=0.5,
)
mpi_hist_omz60_vol_rg = regrid_horizontal(
    native_results[("cmip6_historical", "MPI-ESM1-2-LR")]["omz60"]["omz_volume_cells"],
    target="woa23",
    method="conservative",
    skipna=True,
    na_thres=0.5,
)

historical_basin_results_60 = {"GFDL-ESM4": {}, "MPI-ESM1-2-LR": {}}
for basin in ["Pacific", "Indian", "Atlantic"]:
    gfdl_b = select_basin(gfdl_hist_omz60_vol_rg, basin=basin, include_marginal_seas=True, depth_dependent=False)
    mpi_b = select_basin(mpi_hist_omz60_vol_rg, basin=basin, include_marginal_seas=True, depth_dependent=False)
    historical_basin_results_60["GFDL-ESM4"][basin] = float(_sum_over_all_dims(gfdl_b).compute())
    historical_basin_results_60["MPI-ESM1-2-LR"][basin] = float(_sum_over_all_dims(mpi_b).compute())

# Regrid scenario-change thickness fields and CESM2 ssp245 thickness for Figure 2
gfdl_change_omz60_native = (
    native_results[("cmip6_ssp245", "GFDL-ESM4")]["omz60"]["omz_thickness"]
    - native_results[("cmip6_historical", "GFDL-ESM4")]["omz60"]["omz_thickness"]
).assign_attrs(units="m", long_name="OMZ thickness change")
mpi_change_omz60_native = (
    native_results[("cmip6_ssp245", "MPI-ESM1-2-LR")]["omz60"]["omz_thickness"]
    - native_results[("cmip6_historical", "MPI-ESM1-2-LR")]["omz60"]["omz_thickness"]
).assign_attrs(units="m", long_name="OMZ thickness change")
cesm2_ssp245_omz60_thick_native = native_results[("cmip6_ssp245", "CESM2")]["omz60"]["omz_thickness"]

gfdl_change_omz60_rg = regrid_horizontal(
    gfdl_change_omz60_native, target="woa23", method="bilinear", skipna=True, na_thres=0.5
)
mpi_change_omz60_rg = regrid_horizontal(
    mpi_change_omz60_native, target="woa23", method="bilinear", skipna=True, na_thres=0.5
)
cesm2_ssp245_omz60_thick_rg = regrid_horizontal(
    cesm2_ssp245_omz60_thick_native, target="woa23", method="bilinear", skipna=True, na_thres=0.5
)

# Materialise mapped fields used in figures
woa_omz60_thickness_plot = woa_omz60_thickness.compute()
gfdl_hist_omz60_thick_plot = gfdl_hist_omz60_thick_rg.compute()
mpi_hist_omz60_thick_plot = mpi_hist_omz60_thick_rg.compute()
gfdl_change_omz60_plot = gfdl_change_omz60_rg.compute()
mpi_change_omz60_plot = mpi_change_omz60_rg.compute()
cesm2_ssp245_omz60_thick_plot = cesm2_ssp245_omz60_thick_rg.compute()

# Assemble concise report tables
summary_rows = []

# WOA row
summary_rows.append({
    "source": "WOA23",
    "omz_lt20_1e15m3": woa_omz20["total_volume_m3"] / 1e15,
    "omz_lt60_1e15m3": woa_omz60["total_volume_m3"] / 1e15,
    "max_thickness60_m": woa_omz60["max_thickness_m"],
    "max_thickness60_lat": woa_omz60["max_coords"].get("lat", np.nan),
    "max_thickness60_lon": woa_omz60["max_coords"].get("lon", np.nan),
    "inventory_moles": np.nan,
})

for (dataset, model), out in native_results.items():
    label = f"{model} {'historical' if dataset == 'cmip6_historical' else 'ssp245'}"
    coords60 = out["omz60"]["max_coords"]
    summary_rows.append({
        "source": label,
        "omz_lt20_1e15m3": out["omz20"]["total_volume_m3"] / 1e15,
        "omz_lt60_1e15m3": out["omz60"]["total_volume_m3"] / 1e15,
        "max_thickness60_m": out["omz60"]["max_thickness_m"],
        "max_thickness60_lat": coords60.get("lat", np.nan),
        "max_thickness60_lon": coords60.get("lon", np.nan),
        "inventory_moles": out["inventory_moles"],
    })

summary_table = pd.DataFrame(summary_rows)

scenario_change_table = pd.DataFrame([
    {
        "model": m,
        "pct_change_lt20": scenario_change_pct[m]["lt20_pct"],
        "pct_change_lt60": scenario_change_pct[m]["lt60_pct"],
    }
    for m in ["GFDL-ESM4", "MPI-ESM1-2-LR"]
])

basin_table_60 = pd.DataFrame([
    {
        "source": "WOA23",
        "Pacific_1e15m3": woa_basin_results_60["Pacific"] / 1e15,
        "Indian_1e15m3": woa_basin_results_60["Indian"] / 1e15,
        "Atlantic_1e15m3": woa_basin_results_60["Atlantic"] / 1e15,
    },
    {
        "source": "GFDL-ESM4 historical",
        "Pacific_1e15m3": historical_basin_results_60["GFDL-ESM4"]["Pacific"] / 1e15,
        "Indian_1e15m3": historical_basin_results_60["GFDL-ESM4"]["Indian"] / 1e15,
        "Atlantic_1e15m3": historical_basin_results_60["GFDL-ESM4"]["Atlantic"] / 1e15,
    },
    {
        "source": "MPI-ESM1-2-LR historical",
        "Pacific_1e15m3": historical_basin_results_60["MPI-ESM1-2-LR"]["Pacific"] / 1e15,
        "Indian_1e15m3": historical_basin_results_60["MPI-ESM1-2-LR"]["Indian"] / 1e15,
        "Atlantic_1e15m3": historical_basin_results_60["MPI-ESM1-2-LR"]["Atlantic"] / 1e15,
    },
])

# Better-model summary against WOA23 global OMZ volume
model_fit_rows = []
for model in ["GFDL-ESM4", "MPI-ESM1-2-LR"]:
    err20 = abs(native_results[("cmip6_historical", model)]["omz20"]["total_volume_m3"] - woa_omz20["total_volume_m3"])
    err60 = abs(native_results[("cmip6_historical", model)]["omz60"]["total_volume_m3"] - woa_omz60["total_volume_m3"])
    model_fit_rows.append({
        "model": model,
        "abs_error_lt20_1e15m3": err20 / 1e15,
        "abs_error_lt60_1e15m3": err60 / 1e15,
        "combined_abs_error_1e15m3": (err20 + err60) / 1e15,
    })
model_fit_table = pd.DataFrame(model_fit_rows).sort_values("combined_abs_error_1e15m3").reset_index(drop=True)
better_model = model_fit_table.loc[0, "model"]

scenario_sign_agreement = (
    np.sign(scenario_change_pct["GFDL-ESM4"]["lt20_pct"]) == np.sign(scenario_change_pct["MPI-ESM1-2-LR"]["lt20_pct"])
    and np.sign(scenario_change_pct["GFDL-ESM4"]["lt60_pct"]) == np.sign(scenario_change_pct["MPI-ESM1-2-LR"]["lt60_pct"])
)
scenario_change_assessment = {
    "better_model_global_omz_volume": better_model,
    "models_agree_on_sign_of_change": bool(scenario_sign_agreement),
    "gfdl_pct_change_lt20": scenario_change_pct["GFDL-ESM4"]["lt20_pct"],
    "gfdl_pct_change_lt60": scenario_change_pct["GFDL-ESM4"]["lt60_pct"],
    "mpi_pct_change_lt20": scenario_change_pct["MPI-ESM1-2-LR"]["lt20_pct"],
    "mpi_pct_change_lt60": scenario_change_pct["MPI-ESM1-2-LR"]["lt60_pct"],
}

# Rounded versions for direct reporting
summary_table_report = summary_table.copy()
for col in ["omz_lt20_1e15m3", "omz_lt60_1e15m3"]:
    summary_table_report[col] = summary_table_report[col].map(lambda x: float(f"{x:.3g}") if pd.notnull(x) else np.nan)
summary_table_report["max_thickness60_m"] = summary_table_report["max_thickness60_m"].map(lambda x: round(float(x), 1))
summary_table_report["max_thickness60_lat"] = summary_table_report["max_thickness60_lat"].map(lambda x: round(float(x), 2) if pd.notnull(x) else np.nan)
summary_table_report["max_thickness60_lon"] = summary_table_report["max_thickness60_lon"].map(lambda x: round(float(x), 2) if pd.notnull(x) else np.nan)
summary_table_report["inventory_moles"] = summary_table_report["inventory_moles"].map(lambda x: float(f"{x:.3g}") if pd.notnull(x) else np.nan)

scenario_change_table_report = scenario_change_table.copy()
scenario_change_table_report["pct_change_lt20"] = scenario_change_table_report["pct_change_lt20"].map(lambda x: round(float(x), 2))
scenario_change_table_report["pct_change_lt60"] = scenario_change_table_report["pct_change_lt60"].map(lambda x: round(float(x), 2))

basin_table_60_report = basin_table_60.copy()
for col in ["Pacific_1e15m3", "Indian_1e15m3", "Atlantic_1e15m3"]:
    basin_table_60_report[col] = basin_table_60_report[col].map(lambda x: float(f"{x:.3g}"))

woa_omz60_context = {
    "mean_aou_umolkg": round(woa_aou_mean_omz60, 1),
    "mean_percent_o2_saturation": round(woa_o2sat_mean_omz60, 1),
    "interpretation": (
        "high apparent oxygen utilization and low percent oxygen saturation indicate an OMZ maintained primarily by respiration acting on poorly ventilated waters"
        if (woa_aou_mean_omz60 > 0 and woa_o2sat_mean_omz60 < 100)
        else "the AOU and oxygen saturation do not indicate a classic respiration-dominated, poorly ventilated OMZ"
    ),
}
