import numpy as np

HIST_TR = ["1971-01", "2000-12"]
SSP_TR = ["2071-01", "2100-12"]
BASINS = ["Pacific", "Indian", "Atlantic"]

# ---------------------------------------------------------------- helpers
def load_mean(dataset, var, model, tr):
    """30-year day-weighted long-term mean on the native grid; values near 1e20 are fill."""
    da = load_regular_data(dataset=dataset, variable=var, model=model, time_range=tr)
    da = da.where(da < 1e19)
    return temporal_aggregate(da, freq="mean").compute()

def omz_stats(o2u):
    """OMZ volumes (1e15 m3) at 20/60, thickness at 60 and its max location, basin volumes at 60."""
    out = {}
    for th in (20, 60):
        r = compute_exceedance_fraction(o2u, threshold=th, comparison="less", weighting="volume")
        out[f"vol{th}"] = float(r["exceeding_amount"]) / 1e15
        out[f"weighting{th}"] = r["weighting"]
        out[f"criterion{th}"] = r["criterion"]
    thk, _ = identify_oxygen_minimum_zone(o2u, threshold=60)
    thk = thk.compute()
    ext = locate_extreme(thk, op="max")
    out["thk60"] = thk
    out["max_thk60"] = float(ext["value"])
    out["max_thk60_coords"] = ext["coordinates"]
    out["max_thk60_ties"] = ext["ties"]
    for b in BASINS:
        sb = select_basin(o2u, b, include_marginal_seas=True, depth_dependent=False)
        rb = compute_exceedance_fraction(sb, threshold=60, comparison="less", weighting="volume")
        out[f"basin_{b}"] = float(rb["exceeding_amount"]) / 1e15
    return out

# ---------------------------------------------------------------- WOA23
woa_o2 = load_climatology_data(dataset="woa23", variable="oxygen").compute()
woa_aou = load_climatology_data(dataset="woa23", variable="apparent_oxygen_utilization").compute()
woa_psat = load_climatology_data(dataset="woa23", variable="percent_oxygen_saturation").compute()
landsea = load_landsea_mask()
basinmask = load_basin_mask()

woa_res = omz_stats(woa_o2)

# OMZ-mean AOU and percent saturation (cells with oxygen < 60 umol/kg), volume weighted
omz_cells = woa_o2 < 60
aou_omz = woa_aou.where(omz_cells)
psat_omz = woa_psat.where(omz_cells)
woa_aou_mean = float(compute_basic_stat(aou_omz, op="volume_weighted_mean", dim=["depth", "lat", "lon"]))
woa_psat_mean = float(compute_basic_stat(psat_omz, op="volume_weighted_mean", dim=["depth", "lat", "lon"]))
if woa_aou_mean > 100 and woa_psat_mean < 40:
    omz_interp = ("High AOU and low percent saturation: the water has accumulated a large respiration "
                  "(remineralisation) oxygen debt because it is poorly ventilated (long time since "
                  "last surface contact); the zone is set by respiration acting on poorly ventilated water.")
else:
    omz_interp = ("Low AOU and relatively high saturation: the zone is mainly set by ventilation rather "
                  "than by strong respiration.")
omz_interp += (" Caveat: AOU = O2sat - O2 and %sat = 100*O2/O2sat both depend on the solubility "
               "formulation used for O2sat (WOA uses Garcia & Gordon, umol/kg) and AOU assumes surface "
               "saturation at formation, so the two are not independent and preformed disequilibrium biases AOU.")

# ---------------------------------------------------------------- CMIP6 model processing
runs = [
    ("cmip6_historical", "GFDL-ESM4", HIST_TR),
    ("cmip6_historical", "MPI-ESM1-2-LR", HIST_TR),
    ("cmip6_ssp245", "CESM2", SSP_TR),
    ("cmip6_ssp245", "GFDL-ESM4", SSP_TR),
    ("cmip6_ssp245", "MPI-ESM1-2-LR", SSP_TR),
]

res = {}          # (dataset, model) -> stats
grid_info = {}
for ds, mdl, tr in runs:
    o2 = load_mean(ds, "o2", mdl, tr)                 # mol m-3
    thetao = load_mean(ds, "thetao", mdl, tr)
    so = load_mean(ds, "so", mdl, tr)

    # in-situ density anomaly (ref_pressure=None), passed unchanged to transform_unit
    sigma = compute_potential_density_anomaly(thetao, so, ref_pressure=None)
    o2u = transform_unit(o2, from_unit=o2.attrs.get("units", "mol m-3"), to_unit="umol/kg", density=sigma)
    o2u = o2u.where(o2u < 1e19)

    r = omz_stats(o2u)
    r["density_assumption"] = o2u.attrs.get("density_assumption")
    r["density_mean_insitu_kgm3"] = float(compute_basic_stat(sigma, op="volume_weighted_mean",
                                                             dim=["depth"] + [d for d in sigma.dims if d not in ("depth",)])) + 1000.0

    # global oxygen inventory (volume-basis native field, model's own volcello)
    content = convert_concentration_to_content(o2)
    r["inventory_mol"] = float(compute_basic_stat(content, op="sum"))
    r["inventory_volume_source"] = content.attrs.get("volume_source")

    # grid provenance
    vol = load_volumes_for_cmip6(ds, mdl)
    vol_static = vol.isel(time=0) if "time" in vol.dims else vol
    area = load_areas(dataset=ds, model=mdl)
    grid_info[(ds, mdl)] = dict(
        volume_source=vol.attrs.get("volume_source"),
        area_source=area.attrs.get("area_source"),
        total_ocean_volume_1e15=float(vol_static.sum()) / 1e15,
        total_ocean_area_1e14=float(area.sum()) / 1e14,
    )
    res[(ds, mdl)] = r

H = "cmip6_historical"; S = "cmip6_ssp245"
gh, mh = res[(H, "GFDL-ESM4")], res[(H, "MPI-ESM1-2-LR")]
gs, ms, cs = res[(S, "GFDL-ESM4")], res[(S, "MPI-ESM1-2-LR")], res[(S, "CESM2")]

# ---------------------------------------------------------------- scenario change (r1 windows, GFDL & MPI only)
pct_change = {}
for name, hh, ss in (("GFDL-ESM4", gh, gs), ("MPI-ESM1-2-LR", mh, ms)):
    pct_change[name] = {
        "vol20_pct": round((ss["vol20"] - hh["vol20"]) / hh["vol20"] * 100, 2),
        "vol60_pct": round((ss["vol60"] - hh["vol60"]) / hh["vol60"] * 100, 2),
    }
same_sign_20 = np.sign(pct_change["GFDL-ESM4"]["vol20_pct"]) == np.sign(pct_change["MPI-ESM1-2-LR"]["vol20_pct"])
same_sign_60 = np.sign(pct_change["GFDL-ESM4"]["vol60_pct"]) == np.sign(pct_change["MPI-ESM1-2-LR"]["vol60_pct"])

# ---------------------------------------------------------------- model vs WOA23 (volumes)
keys = ["vol20", "vol60"] + [f"basin_{b}" for b in BASINS]
vol_compare = {}
for name, hh in (("GFDL-ESM4", gh), ("MPI-ESM1-2-LR", mh)):
    bias = {k: hh[k] - woa_res[k] for k in keys}
    relerr = {k: (hh[k] - woa_res[k]) / woa_res[k] * 100 for k in keys}
    vol_compare[name] = dict(bias_1e15=bias, rel_error_pct=relerr,
                             mean_abs_rel_error_pct=float(np.mean([abs(v) for v in relerr.values()])))
closer_model = min(vol_compare, key=lambda n: vol_compare[n]["mean_abs_rel_error_pct"])
closer_model_vol60 = min(vol_compare, key=lambda n: abs(vol_compare[n]["rel_error_pct"]["vol60"]))

# ---------------------------------------------------------------- maps on the 1-degree grid (display / cell-by-cell only)
def to_1deg(da):
    return regrid_horizontal(da, target="woa23", method="conservative").compute()

g_hist_r = to_1deg(gh["thk60"])
m_hist_r = to_1deg(mh["thk60"])
woa_thk_r = woa_res["thk60"]  # already on the WOA23 1-degree grid

mask1 = build_comparison_mask([g_hist_r, m_hist_r, woa_thk_r], ocean_only=True)
fig1_cells = dict(mask1.attrs)
fig1_ncells = int(mask1.sum())
g_hist_m = g_hist_r.where(mask1)
m_hist_m = m_hist_r.where(mask1)
woa_thk_m = woa_thk_r.where(mask1)

map_metrics = {
    "GFDL-ESM4": compute_difference_metrics(g_hist_m, woa_thk_m),
    "MPI-ESM1-2-LR": compute_difference_metrics(m_hist_m, woa_thk_m),
}
map_summary = {k: {m: float(v[m]) for m in ("bias", "rmse", "centered_rmse", "std_ratio")}
               for k, v in map_metrics.items()}

g_dthk = to_1deg(gs["thk60"] - gh["thk60"])
m_dthk = to_1deg(ms["thk60"] - mh["thk60"])
c_thk_r = to_1deg(cs["thk60"])
mask2 = build_comparison_mask([g_dthk, m_dthk, c_thk_r], ocean_only=True)
fig2_cells = dict(mask2.attrs)
fig2_ncells = int(mask2.sum())
g_dthk_m = g_dthk.where(mask2)
m_dthk_m = m_dthk.where(mask2)
c_thk_m = c_thk_r.where(mask2)

# ---------------------------------------------------------------- summary table
summary = {
    "WOA23": {k: woa_res[k] for k in keys + ["max_thk60", "max_thk60_coords"]},
    "GFDL-ESM4 hist": {k: gh[k] for k in keys + ["max_thk60", "max_thk60_coords", "inventory_mol"]},
    "MPI-ESM1-2-LR hist": {k: mh[k] for k in keys + ["max_thk60", "max_thk60_coords", "inventory_mol"]},
    "CESM2 ssp245": {k: cs[k] for k in ["vol20", "vol60", "max_thk60", "max_thk60_coords", "inventory_mol"]},
    "GFDL-ESM4 ssp245": {k: gs[k] for k in ["vol20", "vol60", "max_thk60", "max_thk60_coords", "inventory_mol"]},
    "MPI-ESM1-2-LR ssp245": {k: ms[k] for k in ["vol20", "vol60", "max_thk60", "max_thk60_coords", "inventory_mol"]},
}
notes = ("Volumes, inventories, basin totals and OMZ thickness are on each source's native grid and levels; "
         "the 1-degree conservative regridding and shared mask are used only for the maps and cell-by-cell metrics. "
         f"Density for unit conversion: model in-situ density field (ref_pressure=None), mmol m-3 -> umol/kg = x*1000/rho; "
         f"density_assumption (GFDL hist): {gh['density_assumption']}. "
         "Basins: select_basin on basinmask_01 with include_marginal_seas=True, depth_dependent=False; Pacific, Indian, Atlantic only. "
         "CESM2 has no historical o2: no change can be formed for it.")
print(summary)
print(pct_change)
print("WOA OMZ mean AOU, %sat:", woa_aou_mean, woa_psat_mean)
print("closer model (mean abs rel err all volumes):", closer_model, "| at vol60:", closer_model_vol60)
print("same sign of change 20/60:", same_sign_20, same_sign_60)
print(grid_info)
print(notes)
