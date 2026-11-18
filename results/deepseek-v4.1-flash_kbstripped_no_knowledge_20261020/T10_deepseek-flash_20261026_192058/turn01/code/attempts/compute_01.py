import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt


def get_scalar(x):
    if hasattr(x, "compute"):
        return float(x.compute())
    return float(x)


def time_mean(da):
    return temporal_aggregate(da, freq="mean", op="mean")


# -----------------------------------------------------------------------------
# WOA23 observed climatology and basin mask
# -----------------------------------------------------------------------------
o2_woa = load_climatology_data(
    dataset="woa23",
    variable="oxygen",
    time_slice="annual",
    period="1971-2000",
)
aou_woa = load_climatology_data(
    dataset="woa23",
    variable="apparent_oxygen_utilization",
    time_slice="annual",
    period="1971-2000",
)
sat_woa = load_climatology_data(
    dataset="woa23",
    variable="percent_oxygen_saturation",
    time_slice="annual",
    period="1971-2000",
)
basin_mask = load_basin_mask()

o2_woa_umol = transform_unit(
    o2_woa,
    from_unit=o2_woa.attrs.get("units", "umol kg-1"),
    to_unit="umol/kg",
)
aou_woa_umol = aou_woa.assign_attrs(units="umol/kg")
sat_woa = sat_woa.assign_attrs(units="percent")

woa_thick_20, woa_mask_20 = identify_oxygen_minimum_zone(
    o2_woa_umol, threshold=20.0
)
woa_thick_60, woa_mask_60 = identify_oxygen_minimum_zone(
    o2_woa_umol, threshold=60.0
)

woa_res20 = compute_exceedance_fraction(
    o2_woa_umol,
    threshold=20.0,
    comparison="less",
    dim=None,
    weighting="volume",
)
woa_res60 = compute_exceedance_fraction(
    o2_woa_umol,
    threshold=60.0,
    comparison="less",
    dim=None,
    weighting="volume",
)
woa_vol_20 = get_scalar(woa_res20["exceeding_amount"]) / 1e15
woa_vol_60 = get_scalar(woa_res60["exceeding_amount"]) / 1e15
woa_max_60 = locate_extreme(woa_thick_60, op="max")

aou_omz_mean = get_scalar(
    compute_basic_stat(
        aou_woa_umol.where(woa_mask_60),
        op="volume_weighted_mean",
        dim=["depth", "lat", "lon"],
    )
)
sat_omz_mean = get_scalar(
    compute_basic_stat(
        sat_woa.where(woa_mask_60),
        op="volume_weighted_mean",
        dim=["depth", "lat", "lon"],
    )
)

# -----------------------------------------------------------------------------
# Decode WOA23 basin-mask codes and fold marginal seas into parent basins
# -----------------------------------------------------------------------------
basin_names = basin_mask.attrs.get("basin_names")
codes = np.unique(basin_mask.values[~np.isnan(basin_mask.values)]).astype(int)

if isinstance(basin_names, dict):
    code_to_name = {int(k): str(v) for k, v in basin_names.items()}
else:
    names = list(basin_names) if basin_names is not None else []
    code_to_name = {}
    for c in codes:
        if c < len(names):
            code_to_name[c] = str(names[c])
        elif 0 <= c - 1 < len(names):
            code_to_name[c] = str(names[c - 1])
        else:
            code_to_name[c] = str(c)


def centroid_for_code(c):
    mask_c = basin_mask == c
    if "depth" in mask_c.dims:
        mask2 = mask_c.any("depth")
    else:
        mask2 = mask_c
    lat_vals = basin_mask["lat"].values
    lon_vals = basin_mask["lon"].values
    lat_grid, lon_grid = np.meshgrid(lat_vals, lon_vals, indexing="ij")
    m = mask2.transpose("lat", "lon").values
    if not np.any(m):
        return np.nan, np.nan
    return float(np.nanmean(lat_grid[m])), float(np.nanmean(lon_grid[m]))


def centroid_group(lat_c, lon_c):
    if not np.isfinite(lat_c) or not np.isfinite(lon_c):
        return None
    if lat_c < -60.0 or lat_c > 66.0:
        return None
    # Mediterranean, Black Sea, Baltic, North Sea, Norwegian Sea
    if -10.0 <= lon_c <= 42.0 and 30.0 <= lat_c <= 70.0:
        return "Atlantic"
    # Gulf of Mexico, Caribbean, Hudson Bay, Labrador Sea, Baffin Bay
    if -100.0 <= lon_c <= -50.0 and 8.0 <= lat_c <= 70.0:
        return "Atlantic"
    # South China Sea, Philippine Sea, Indonesian seas
    if 100.0 <= lon_c <= 130.0 and -15.0 <= lat_c <= 25.0:
        return "Pacific"
    # Red Sea, Persian Gulf, Arabian Sea, Bay of Bengal
    if 30.0 <= lon_c <= 100.0 and -5.0 <= lat_c <= 30.0:
        return "Indian"
    if 20.0 < lon_c <= 120.0:
        return "Indian"
    if -70.0 <= lon_c <= 20.0:
        return "Atlantic"
    return "Pacific"


code_to_group = {}
for c in codes:
    name = code_to_name.get(c, str(c)).lower()
    if "southern" in name or "arctic" in name:
        group = None
    elif "pacific" in name:
        group = "Pacific"
    elif "atlantic" in name:
        group = "Atlantic"
    elif "indian" in name:
        group = "Indian"
    else:
        lat_c, lon_c = centroid_for_code(c)
        group = centroid_group(lat_c, lon_c)
    code_to_group[c] = group

# WOA23 basin-restricted OMZ volume at the 60 threshold
basin_woa = basin_mask.sel(
    depth=o2_woa_umol["depth"].values, method="nearest"
).assign_coords(depth=o2_woa_umol["depth"])

basin_totals = {}
for group in ["Pacific", "Indian", "Atlantic"]:
    codes_g = [c for c, g in code_to_group.items() if g == group]
    group_mask = basin_woa.isin(codes_g)
    res = compute_exceedance_fraction(
        o2_woa_umol.where(group_mask),
        threshold=60.0,
        comparison="less",
        dim=None,
        weighting="volume",
    )
    basin_totals[("WOA23", group)] = get_scalar(res["exceeding_amount"]) / 1e15


# -----------------------------------------------------------------------------
# CMIP6 model processing
# -----------------------------------------------------------------------------
def model_basin_mask_on_grid(model_da, bm):
    bm_lat = "lat" if "lat" in bm.dims else bm.attrs.get("lat_coord", "lat")
    bm_lon = "lon" if "lon" in bm.dims else bm.attrs.get("lon_coord", "lon")

    bm_d = bm.sel(depth=model_da["depth"].values, method="nearest")

    lat_name = model_da.attrs.get("lat_coord", "lat")
    lon_name = model_da.attrs.get("lon_coord", "lon")

    if lat_name in model_da.dims and lon_name in model_da.dims:
        bm_m = bm_d.sel(
            {
                bm_lat: model_da[lat_name].values,
                bm_lon: model_da[lon_name].values,
            },
            method="nearest",
        )
        bm_m = bm_m.rename({bm_lat: lat_name, bm_lon: lon_name})
        bm_m = bm_m.assign_coords(
            {lat_name: model_da[lat_name], lon_name: model_da[lon_name]}
        )
    else:
        lat2 = model_da[lat_name]
        lon2 = model_da[lon_name]
        bm_m = bm_d.sel(
            {
                bm_lat: xr.DataArray(lat2.values, dims=lat2.dims),
                bm_lon: xr.DataArray(lon2.values, dims=lon2.dims),
            },
            method="nearest",
        )
        bm_m = bm_m.assign_coords(
            {lat_name: model_da[lat_name], lon_name: model_da[lon_name]}
        )

    bm_m = bm_m.assign_coords(depth=model_da["depth"])
    dims_order = [d for d in model_da.dims if d in bm_m.dims]
    bm_m = bm_m.transpose(*dims_order)
    return bm_m


def basin_totals_for_model(o2_field, bm, source_label, threshold=60.0):
    bm_model = model_basin_mask_on_grid(o2_field, bm)
    out = {}
    for group in ["Pacific", "Indian", "Atlantic"]:
        codes_g = [c for c, g in code_to_group.items() if g == group]
        group_mask = bm_model.isin(codes_g)
        res = compute_exceedance_fraction(
            o2_field.where(group_mask),
            threshold=threshold,
            comparison="less",
            dim=None,
            weighting="volume",
        )
        out[(source_label, group)] = get_scalar(res["exceeding_amount"]) / 1e15
    return out


def process_model(dataset, model, exp_label, time_range):
    o2_raw = load_regular_data(
        dataset=dataset,
        variable="o2",
        model=model,
        time_range=time_range,
    )
    thetao_raw = load_regular_data(
        dataset=dataset,
        variable="thetao",
        model=model,
        time_range=time_range,
    )
    so_raw = load_regular_data(
        dataset=dataset,
        variable="so",
        model=model,
        time_range=time_range,
    )

    o2_mean = time_mean(o2_raw)
    thetao_mean = time_mean(thetao_raw)
    so_mean = time_mean(so_raw)

    sigma = compute_potential_density_anomaly(
        thetao_mean,
        so_mean,
        ref_pressure=None,
        temperature_type="potential",
    )

    o2_umol = transform_unit(
        o2_mean,
        from_unit=o2_mean.attrs.get("units", "mol m-3"),
        to_unit="umol/kg",
        density=sigma,
    )

    content = convert_concentration_to_content(o2_mean)
    inventory = get_scalar(content.sum(skipna=True))

    thick20, mask20 = identify_oxygen_minimum_zone(o2_umol, threshold=20.0)
    thick60, mask60 = identify_oxygen_minimum_zone(o2_umol, threshold=60.0)

    res20 = compute_exceedance_fraction(
        o2_umol,
        threshold=20.0,
        comparison="less",
        dim=None,
        weighting="volume",
    )
    res60 = compute_exceedance_fraction(
        o2_umol,
        threshold=60.0,
        comparison="less",
        dim=None,
        weighting="volume",
    )
    vol20 = get_scalar(res20["exceeding_amount"]) / 1e15
    vol60 = get_scalar(res60["exceeding_amount"]) / 1e15

    max60 = locate_extreme(thick60, op="max")

    return {
        "o2_umol": o2_umol,
        "thick20": thick20,
        "mask20": mask20,
        "thick60": thick60,
        "mask60": mask60,
        "vol20": vol20,
        "vol60": vol60,
        "max60": max60,
        "inventory": inventory,
    }


hist_models = ["GFDL-ESM4", "MPI-ESM1-2-LR"]
ssp_models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

hist_results = {}
for model in hist_models:
    hist_results[model] = process_model(
        "cmip6_historical",
        model,
        "historical",
        ["1971-01", "2000-12"],
    )

ssp_results = {}
for model in ssp_models:
    ssp_results[model] = process_model(
        "cmip6_ssp245",
        model,
        "ssp245",
        ["2071-01", "2100-12"],
    )


# -----------------------------------------------------------------------------
# Assemble result dictionaries
# -----------------------------------------------------------------------------
omz_volumes = {}
max_thickness = {}
inventories = {}
o2_fields = {}
thick60_fields = {}

omz_volumes[("WOA23", "climatology", 20)] = woa_vol_20
omz_volumes[("WOA23", "climatology", 60)] = woa_vol_60
max_thickness[("WOA23", "climatology", 60)] = woa_max_60
o2_fields[("WOA23", "climatology")] = o2_woa_umol
thick60_fields[("WOA23", "climatology")] = woa_thick_60

for model in hist_models:
    res = hist_results[model]
    omz_volumes[(model, "historical", 20)] = res["vol20"]
    omz_volumes[(model, "historical", 60)] = res["vol60"]
    max_thickness[(model, "historical", 60)] = res["max60"]
    inventories[(model, "historical")] = res["inventory"]
    o2_fields[(model, "historical")] = res["o2_umol"]
    thick60_fields[(model, "historical")] = res["thick60"]

for model in ssp_models:
    res = ssp_results[model]
    omz_volumes[(model, "ssp245", 20)] = res["vol20"]
    omz_volumes[(model, "ssp245", 60)] = res["vol60"]
    max_thickness[(model, "ssp245", 60)] = res["max60"]
    inventories[(model, "ssp245")] = res["inventory"]
    o2_fields[(model, "ssp245")] = res["o2_umol"]
    thick60_fields[(model, "ssp245")] = res["thick60"]

scenario_changes = {}
for model in hist_models:
    for thr in [20, 60]:
        h = omz_volumes[(model, "historical", thr)]
        s = omz_volumes[(model, "ssp245", thr)]
        scenario_changes[(model, thr)] = (s - h) / h * 100.0

for model in hist_models:
    basin_totals.update(
        basin_totals_for_model(
            hist_results[model]["o2_umol"],
            basin_mask,
            model,
            60.0,
        )
    )

thick_change = {}
for model in hist_models:
    thick_change[model] = (
        ssp_results[model]["thick60"] - hist_results[model]["thick60"]
    )


# -----------------------------------------------------------------------------
# Regrid thickness maps to the WOA23 plotting grid
# -----------------------------------------------------------------------------
def copy_grid_attrs(da, source_da):
    out = da.copy()
    out.attrs.update(
        {k: v for k, v in source_da.attrs.items() if k != "units"}
    )
    out.attrs["units"] = "m"
    return out


thick_hist_plot = {}
for model in hist_models:
    da = copy_grid_attrs(
        hist_results[model]["thick60"],
        hist_results[model]["o2_umol"],
    )
    thick_hist_plot[model] = regrid_horizontal(
        da, target="woa23", method="bilinear"
    ).compute()

thick_ssp_plot = {}
for model in ssp_models:
    da = copy_grid_attrs(
        ssp_results[model]["thick60"],
        ssp_results[model]["o2_umol"],
    )
    thick_ssp_plot[model] = regrid_horizontal(
        da, target="woa23", method="bilinear"
    ).compute()

thick_change_plot = {}
for model in hist_models:
    da = copy_grid_attrs(
        thick_change[model],
        ssp_results[model]["o2_umol"],
    )
    thick_change_plot[model] = regrid_horizontal(
        da, target="woa23", method="bilinear"
    ).compute()

thick_woa_plot = woa_thick_60.compute()


# -----------------------------------------------------------------------------
# Summary tables
# -----------------------------------------------------------------------------
volume_rows = []
for (source, exp, thr), val in sorted(omz_volumes.items()):
    volume_rows.append(
        {
            "source": source,
            "experiment": exp,
            "threshold_umol_kg": thr,
            "OMZ_volume_1e15_m3": round(val, 3),
        }
    )
omz_volume_table = pd.DataFrame(volume_rows)


def get_coord(coords, names):
    for n in names:
        if n in coords:
            return coords[n]
    return np.nan


max_rows = []
for (source, exp, thr), md in sorted(max_thickness.items()):
    coords = md.get("coordinates", {})
    lat = get_coord(coords, ["lat", "latitude"])
    lon = get_coord(coords, ["lon", "longitude"])
    max_rows.append(
        {
            "source": source,
            "experiment": exp,
            "threshold_umol_kg": thr,
            "max_thickness_m": round(float(md["value"]), 1),
            "lat": round(float(lat), 2),
            "lon": round(float(lon), 2),
        }
    )
max_thickness_table = pd.DataFrame(max_rows)

inv_rows = []
for (source, exp), val in sorted(inventories.items()):
    inv_rows.append(
        {
            "model": source,
            "experiment": exp,
            "oxygen_inventory_mol": float(f"{val:.3g}"),
        }
    )
inventory_table = pd.DataFrame(inv_rows)

basin_rows = []
for (source, group), val in sorted(basin_totals.items()):
    basin_rows.append(
        {
            "source": source,
            "basin": group,
            "OMZ_volume_1e15_m3": round(val, 3),
        }
    )
basin_table = pd.DataFrame(basin_rows)

change_rows = []
for (model, thr), pct in sorted(scenario_changes.items()):
    change_rows.append(
        {
            "model": model,
            "threshold_umol_kg": thr,
            "percent_change": round(pct, 2),
        }
    )
scenario_change_table = pd.DataFrame(change_rows)

obs20 = omz_volumes[("WOA23", "climatology", 20)]
obs60 = omz_volumes[("WOA23", "climatology", 60)]
compare_rows = []
for model in hist_models:
    for thr, obs in [(20, obs20), (60, obs60)]:
        m = omz_volumes[(model, "historical", thr)]
        compare_rows.append(
            {
                "model": model,
                "threshold_umol_kg": thr,
                "WOA23_1e15_m3": round(obs, 3),
                "model_1e15_m3": round(m, 3),
                "abs_diff_1e15_m3": round(abs(m - obs), 3),
                "abs_pct_diff": round(abs(m - obs) / obs * 100.0, 2),
            }
        )
model_comparison_table = pd.DataFrame(compare_rows)

closer_20 = min(
    hist_models,
    key=lambda m: abs(omz_volumes[(m, "historical", 20)] - obs20),
)
closer_60 = min(
    hist_models,
    key=lambda m: abs(omz_volumes[(m, "historical", 60)] - obs60),
)

if aou_omz_mean > 50.0 and sat_omz_mean < 30.0:
    omz_interpretation = (
        "High AOU with low O2 saturation indicates respiration-driven oxygen "
        "consumption dominates within the OMZ."
    )
elif aou_omz_mean < 30.0 and sat_omz_mean < 30.0:
    omz_interpretation = (
        "Low AOU with low O2 saturation suggests poor ventilation, rather "
        "than intense respiration, controls the OMZ."
    )
else:
    omz_interpretation = (
        "The AOU and saturation values indicate a mixture of respiration and "
        "ventilation influences on the OMZ."
    )

print("OMZ volumes (1e15 m3):")
print(omz_volume_table.to_string(index=False))
print("\nMaximum OMZ thickness at 60 umol/kg:")
print(max_thickness_table.to_string(index=False))
print("\nOxygen inventories (mol):")
print(inventory_table.to_string(index=False))
print("\nBasin OMZ volumes at 60 umol/kg (1e15 m3):")
print(basin_table.to_string(index=False))
print("\nScenario change (%):")
print(scenario_change_table.to_string(index=False))
print("\nModel-observation comparison:")
print(model_comparison_table.to_string(index=False))
print(f"\nCloser to WOA23 at 20 umol/kg: {closer_20}")
print(f"Closer to WOA23 at 60 umol/kg: {closer_60}")
print(
    f"WOA23 OMZ mean AOU: {aou_omz_mean:.1f} umol/kg; "
    f"mean O2 saturation: {sat_omz_mean:.1f}%"
)
print(omz_interpretation)
