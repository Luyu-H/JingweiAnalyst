# =====================================================================
# Ocean carbon system: 3 CMIP6 models (2002 annual mean) vs GLODAPv2
# =====================================================================
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
DS = "cmip6_historical"

def depth_dim(da):
    for n in ("depth", "lev", "olevel", "z", "st_ocean"):
        if n in da.dims:
            return n
    raise KeyError(f"no depth dimension found in {da.dims}")

def restore(da, attrs, units, long_name=None):
    """Re-attach the loader provenance (grid identity) after an operation changed the quantity."""
    d = {k: v for k, v in attrs.items() if k != "units"}
    d["units"] = units
    if long_name is not None:
        d["long_name"] = long_name
    return da.assign_attrs(**d)

def per_kg(field, rho_abs, attrs, long_name):
    """mol m-3 -> mmol m-3 -> umol kg-1 using the cell-by-cell in-situ density (absolute, kg/m3)."""
    mm = transform_unit(field, field.attrs.get("units", "mol m-3"), "mmol m-3")
    out = mm * 1000.0 / rho_abs
    out = restore(out, attrs, "umol kg-1", long_name)
    out.attrs["density_assumption"] = (
        "cell-by-cell in-situ density (sigma computed with TEOS-10, ref_pressure=None, +1000); "
        "umol kg-1 = mmol m-3 * 1000 / rho")
    return out

def global_mean_0_50(da):
    """Single volume-weighted global mean over all of the source's own 0-50 m cells."""
    dd = depth_dim(da)
    sub = da.sel({dd: slice(0.0, 50.0)})
    return float(compute_basic_stat(data=sub, op="volume_weighted_mean", dim=None).compute())

# ---------------------------------------------------------------------
# Step 1: load the CMIP6 fields (monthly 2002; monthly 1971-2000 for DIC/T/S)
# ---------------------------------------------------------------------
state2002, attrs2002 = {}, {}
state_hist, attrs_hist = {}, {}

for m in models:
    st, at = {}, {}
    for v in ["dissic", "talk", "ph", "thetao", "so"]:
        da = load_regular_data(dataset=DS, variable=v, model=m,
                               time_range=["2002-01", "2002-12"])
        at[v] = dict(da.attrs)
        st[v] = da.mean(dim="time")            # 2002 annual-mean state
    state2002[m], attrs2002[m] = st, at

    h, ha = {}, {}
    for v in ["dissic", "thetao", "so"]:
        da = load_regular_data(dataset=DS, variable=v, model=m,
                               time_range=["1971-01", "2000-12"],
                               depth_range=[0.0, 50.0])
        ha[v] = dict(da.attrs)
        h[v] = da.mean(dim="time")             # 1971-2000 mean state, already 0-50 m
    state_hist[m], attrs_hist[m] = h, ha

# grid-weight provenance (lazy, only attrs are read)
grid_provenance = {}
for m in models:
    a = load_areas(dataset=DS, model=m, variable="thetao")
    v = load_volumes_for_cmip6(dataset=DS, model=m)
    grid_provenance[m] = {
        "area_source": a.attrs.get("area_source"),
        "volume_source": v.attrs.get("volume_source"),
    }

# ---------------------------------------------------------------------
# Step 3: put the carbon fields on a per-kilogram basis
# ---------------------------------------------------------------------
dic_kg, talk_kg, ph_kg, dic_hist_kg = {}, {}, {}, {}
density_note = {}

for m in models:
    st = state2002[m]
    rho = compute_potential_density_anomaly(temp_data=st["thetao"], sal_data=st["so"],
                                            ref_pressure=None, temperature_type="auto")
    rho_abs = rho + 1000.0

    dic_kg[m]  = per_kg(st["dissic"], rho_abs, attrs2002[m]["dissic"], "DIC")
    talk_kg[m] = per_kg(st["talk"],   rho_abs, attrs2002[m]["talk"],   "Total alkalinity")
    ph_kg[m]   = restore(st["ph"], attrs2002[m]["ph"], "1", "pH (total scale, in situ)")
    density_note[m] = dic_kg[m].attrs["density_assumption"]

    h = state_hist[m]
    rho_h = compute_potential_density_anomaly(temp_data=h["thetao"], sal_data=h["so"],
                                              ref_pressure=None, temperature_type="auto")
    dic_hist_kg[m] = per_kg(h["dissic"], rho_h + 1000.0, attrs_hist[m]["dissic"], "DIC")

# ---------------------------------------------------------------------
# Step 2: GLODAP reference climatology (full depth, subset later)
# ---------------------------------------------------------------------
glodap_tco2 = load_climatology_data(dataset="glodap", variable="TCO2")
glodap_talk = load_climatology_data(dataset="glodap", variable="TAlk")
glodap_ph   = load_climatology_data(dataset="glodap", variable="pHtsinsitutp")

# ---------------------------------------------------------------------
# Step 5: one global 0-50 m volume-weighted mean per source
# ---------------------------------------------------------------------
DIC_mean, TALK_mean, PH_mean = {}, {}, {}
for m in models:
    DIC_mean[m]  = global_mean_0_50(dic_kg[m])
    TALK_mean[m] = global_mean_0_50(talk_kg[m])
    PH_mean[m]   = global_mean_0_50(ph_kg[m])

DIC_mean["GLODAP"]  = global_mean_0_50(glodap_tco2)
TALK_mean["GLODAP"] = global_mean_0_50(glodap_talk)
PH_mean["GLODAP"]   = global_mean_0_50(glodap_ph)

DIC_diff  = {m: DIC_mean[m]  - DIC_mean["GLODAP"]  for m in models}
TALK_diff = {m: TALK_mean[m] - TALK_mean["GLODAP"] for m in models}
PH_diff   = {m: PH_mean[m]   - PH_mean["GLODAP"]   for m in models}

# ---------------------------------------------------------------------
# Step 6: 0-50 m DIC maps, regridded to the GLODAP 1-degree grid
# ---------------------------------------------------------------------
dic_map = {}
for m in models:
    dd = depth_dim(dic_kg[m])
    col = compute_basic_stat(data=dic_kg[m].sel({dd: slice(0.0, 50.0)}),
                             op="volume_weighted_mean", dim=dd)
    col = restore(col, attrs2002[m]["dissic"], "umol kg-1", "DIC 0-50 m")
    dic_map[m] = regrid_horizontal(data=col, target="glodap", method="conservative")

glodap_map = compute_basic_stat(data=glodap_tco2.sel(depth=slice(0.0, 50.0)),
                                op="volume_weighted_mean", dim="depth")
glodap_map = restore(glodap_map, glodap_tco2.attrs, "umol kg-1", "DIC 0-50 m")

comparison_mask = build_comparison_mask(
    fields=[dic_map[m] for m in models] + [glodap_map],
    ocean_only=False)

map_metrics, map_corr = {}, {}
for m in models:
    a = dic_map[m].where(comparison_mask)
    b = glodap_map.where(comparison_mask)
    met = compute_difference_metrics(data=a, reference=b, dim=["lat", "lon"])
    cor = analyze_correlation(var1=a, var2=b, dim=["lat", "lon"],
                              method="pearson", weighting="area")
    map_metrics[m] = met
    map_corr[m] = cor

bias_maps = {m: dic_map[m].where(comparison_mask) - glodap_map.where(comparison_mask)
             for m in models}

# ---------------------------------------------------------------------
# Step 7: vertical profile bias at 0, 100, 500, 1000, 2000 m
# ---------------------------------------------------------------------
target_depths = [0.0, 100.0, 500.0, 1000.0, 2000.0]

glodap_prof = regrid_vertical(data=glodap_tco2, target=target_depths, method="linear")

dic_prof = {}
for m in models:
    prof = regrid_vertical(data=dic_kg[m], target=target_depths, method="linear")
    prof = restore(prof, attrs2002[m]["dissic"], "umol kg-1", "DIC")
    dic_prof[m] = regrid_horizontal(data=prof, target="glodap", method="conservative")

profile_bias, profile_ncells = {}, {}
for m in models:
    row, nrow = {}, {}
    for d in target_depths:
        a = dic_prof[m].sel(depth=d, method="nearest")
        b = glodap_prof.sel(depth=d, method="nearest")
        met = compute_difference_metrics(data=a, reference=b, dim=["lat", "lon"])
        row[d] = float(met["bias"])
        nrow[d] = int(met["n_valid"])
    profile_bias[m] = row
    profile_ncells[m] = nrow

# ---------------------------------------------------------------------
# Step 8: 2002 minus 1971-2000 change of the global 0-50 m DIC mean
# ---------------------------------------------------------------------
DIC_hist_mean = {m: global_mean_0_50(dic_hist_kg[m]) for m in models}
period_change = {m: DIC_mean[m] - DIC_hist_mean[m] for m in models}
change_to_bias = {m: period_change[m] / DIC_diff[m] for m in models}

# ---------------------------------------------------------------------
# Summary objects
# ---------------------------------------------------------------------
summary = {
    "DIC_mean": DIC_mean, "DIC_diff": DIC_diff,
    "TALK_mean": TALK_mean, "TALK_diff": TALK_diff,
    "PH_mean": PH_mean, "PH_diff": PH_diff,
    "map_bias": {m: float(map_metrics[m]["bias"]) for m in models},
    "map_rmse": {m: float(map_metrics[m]["rmse"]) for m in models},
    "map_corr": {m: float(map_corr[m]["correlation"]) for m in models},
    "map_n_valid": {m: int(map_metrics[m]["n_valid"]) for m in models},
    "map_mean_model": {m: float(map_metrics[m]["mean_data"]) for m in models},
    "map_mean_glodap": {m: float(map_metrics[m]["mean_reference"]) for m in models},
    "n_mask_cells": int(comparison_mask.sum().compute()),
    "profile_bias": profile_bias,
    "profile_ncells": profile_ncells,
    "period_change": period_change,
    "change_to_bias": change_to_bias,
    "density_note": density_note,
    "grid_provenance": grid_provenance,
}
print(summary)
