models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
basins = ["Pacific", "Atlantic", "Indian"]
t_range = ["1971-01", "2000-12"]

# Vertical axes: top-1000 m search window plus the first level below 1000 m
def depth_limit(dataset, model=None, variable=None):
    kw = dict(dataset=dataset, variable=variable)
    if model is not None:
        kw["model"] = model
    levels, bounds = load_depth_levels(**kw)
    levels = np.asarray(levels, dtype=float)
    below = levels[levels > 1000.0]
    lim = float(below.min()) if below.size else float(levels.max())
    return lim, levels, np.asarray(bounds)

lim = {}
levels_info = {}
for m in models:
    lim[m], lv, bd = depth_limit("cmip6_historical", model=m, variable="no3")
    levels_info[m] = lv
lim["WOA23"], lv, bd = depth_limit("woa23", variable="nitrate")
levels_info["WOA23"] = lv

# Masks (for reporting)
landsea = load_landsea_mask()
n_ocean_cells_woa = int((landsea > 1).sum())
basin_names = load_basin_mask(depth_levels="0m").attrs.get("basin_names")

# WOA23 nitrate (all-period composite), annual climatology, native 1-degree grid
woa_no3 = load_climatology_data(dataset="woa23", variable="nitrate", time_slice="annual",
                                depth_range=[0, lim["WOA23"]])
woa_no3 = woa_no3.where(woa_no3 < 1e10)

def with_attrs(nit, src):
    return nit.assign_attrs({**src.attrs, **nit.attrs, "units": "m", "long_name": "nitracline depth"})

def nitracline(field):
    nit = compute_nutricline_depth(field, min_depth=0, max_depth=1000)
    nvalid = field.notnull().sum("depth")
    nit = nit.where(nvalid >= 3)          # too-short / land columns stay NaN, never 0
    return with_attrs(nit, field)

# Order of operations: 1971-2000 day-weighted mean state -> vertical gradient -> nitracline depth
no3_mean = {}
nit = {}
nit["WOA23"] = nitracline(woa_no3).compute()
for m in models:
    d = load_regular_data(dataset="cmip6_historical", variable="no3", model=m,
                          time_range=t_range, depth_range=[0, lim[m]])
    d = d.where(d < 1e10)                  # unmasked 1e20 fill -> NaN
    dm = temporal_aggregate(d, freq="mean").compute()
    no3_mean[m] = dm
    nit[m] = nitracline(dm).compute()

sources = ["WOA23"] + models

# Global and basin means on native grids (area weighted over each source's valid cells)
global_mean = {}
basin_mean = {s: {} for s in sources}
for s in sources:
    global_mean[s] = float(compute_basic_stat(nit[s], op="area_weighted_mean", dim=["lat", "lon"]).compute())
    for b in basins:
        nb = select_basin(nit[s], basin=b, include_marginal_seas=True, depth_dependent=False)
        basin_mean[s][b] = float(compute_basic_stat(nb, op="area_weighted_mean", dim=["lat", "lon"]).compute())

global_mean_rounded = {s: int(round(v)) for s, v in global_mean.items()}
global_diff = {m: int(round(global_mean[m] - global_mean["WOA23"])) for m in models}
basin_mean_rounded = {s: {b: int(round(v)) for b, v in basin_mean[s].items()} for s in sources}
basin_diff = {m: {b: int(round(basin_mean[m][b] - basin_mean["WOA23"][b])) for b in basins} for m in models}

# Basin ordering (deepest first) and comparison with WOA23's ordering
basin_order = {s: sorted(basins, key=lambda b: -basin_mean[s][b]) for s in sources}
order_matches_woa = {m: basin_order[m] == basin_order["WOA23"] for m in models}

# Shared WOA23 1-degree grid comparison
nit_reg = {m: regrid_horizontal(nit[m], target="woa23", method="conservative") for m in models}
mask = build_comparison_mask([nit["WOA23"]] + [nit_reg[m] for m in models], ocean_only=True)
mask_attrs = dict(mask.attrs)
n_cells_shared = int(mask.sum())
woa_m = nit["WOA23"].where(mask)
reg_m = {m: nit_reg[m].where(mask) for m in models}

metrics = {}
corr = {}
corr_area = {}
for m in models:
    r = compute_difference_metrics(reg_m[m], woa_m)
    metrics[m] = {
        "rmse_m": round(float(r["rmse"]), 1),
        "bias_m": float(r["bias"]),
        "centered_rmse_m": float(r["centered_rmse"]),
        "mean_data": float(r["mean_data"]),
        "mean_reference": float(r["mean_reference"]),
        "std_ratio": float(r["std_ratio"]),
        "n_valid": int(r["n_valid"]),
    }
    c = analyze_correlation(reg_m[m], woa_m, dim=["lat", "lon"])
    corr[m] = round(float(c["correlation"]), 3)
    ca = analyze_correlation(reg_m[m], woa_m, dim=["lat", "lon"], weighting="area")
    corr_area[m] = round(float(ca["correlation"]), 3)

# Best model: combine RMSE rank and correlation rank (tie -> lower RMSE)
rmse_rank = {m: i for i, m in enumerate(sorted(models, key=lambda k: metrics[k]["rmse_m"]))}
corr_rank = {m: i for i, m in enumerate(sorted(models, key=lambda k: -corr[k]))}
best_model = sorted(models, key=lambda k: (rmse_rank[k] + corr_rank[k], metrics[k]["rmse_m"]))[0]

# ---- Equatorial Pacific section, 2S-2N, 120E-80W, top 1000 m, in mmol m-3 ----
def to360(da):
    da = da.assign_coords(lon=(da.lon % 360)).sortby("lon")
    return da.sel(lon=slice(120, 280))

def latmean(da):
    # 1-degree regular grid, +-2 deg band: cos(lat) is proportional to cell area
    da = da.sel(lat=slice(-2, 2))
    return da.weighted(np.cos(np.deg2rad(da.lat))).mean("lat")

woa_mmol = transform_unit(woa_no3, from_unit="umol/kg", to_unit="mmol/m3")
density_assumption = {"WOA23": woa_mmol.attrs.get("density_assumption")}
sec = {"WOA23": latmean(to360(woa_mmol)).sel(depth=slice(0, 1000)).compute()}
line = {"WOA23": latmean(to360(nit["WOA23"])).compute()}

for m in models:
    d = load_regular_data(dataset="cmip6_historical", variable="no3", model=m,
                          time_range=t_range, depth_range=[0, lim[m]],
                          lat_range=[-4, 4], lon_range=[118, 282])
    d = d.where(d < 1e10)
    dm = temporal_aggregate(d, freq="mean")
    dm_mmol = transform_unit(dm, from_unit="mol m-3", to_unit="mmol/m3")
    density_assumption[m] = dm_mmol.attrs.get("density_assumption")
    reg = regrid_horizontal(dm_mmol, target="woa23", method="bilinear")
    sec[m] = latmean(to360(reg)).sel(depth=slice(0, 1000)).compute()
    line[m] = latmean(to360(nit_reg[m])).compute()

pooled_sec = np.concatenate([sec[s].values.ravel() for s in sources])
sec_vmax = float(np.nanpercentile(pooled_sec, 99))
pooled_nit = np.concatenate([nit["WOA23"].values.ravel()] + [nit_reg[m].values.ravel() for m in models])
nit_vmax = float(np.nanpercentile(pooled_nit[np.isfinite(pooled_nit)], 98))

print("Order: 1971-2000 time mean -> vertical gradient -> nitracline depth (search top 1000 m).")
print("Models are per volume (mol m-3), WOA23 per mass (umol kg-1); nitracline depth is scale-invariant.")
print("GFDL-ESM4 is gr (1-degree archive grid); CESM2 and MPI-ESM1-2-LR are gn; MPI's coarse grid limits its equatorial Pacific.")
print("Basin means: native grids, WOA surface-label basin mask, marginal seas folded in; Southern and Arctic excluded (global mean includes all ocean cells).")
print("Level windows (limit incl. level below 1000 m):", lim)
print("Global mean nitracline depth (m):", global_mean_rounded, "diff vs WOA23:", global_diff)
print("Basin mean (m):", basin_mean_rounded)
print("Basin diff vs WOA23 (m):", basin_diff)
print("Basin ordering deepest->shallowest:", basin_order, "matches WOA23:", order_matches_woa)
print("Shared mask cells:", n_cells_shared, "attrs:", mask_attrs, "WOA ocean cells:", n_ocean_cells_woa)
print("Shared-mask metrics (not the same as native-grid basin/global mean differences):", metrics)
print("Spatial Pearson r (unweighted):", corr, "area-weighted:", corr_area)
print("Best model (RMSE + correlation + bias):", best_model)
print("Section units mmol m-3; density assumption:", density_assumption)
