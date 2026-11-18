# ============================================================================
# 0. Small helpers
# ============================================================================
def _val(x):
    """Materialise a scalar (float / 0-d dask / 0-d DataArray) into a Python float."""
    if hasattr(x, "compute"):
        return float(x.compute())
    return float(x)


def _retag(da, src, units, long_name=None):
    """Copy the provenance attrs of `src` onto `da` (keeps the grid identity the
    weighting / regridding tools need) and set the physical units explicitly."""
    a = dict(src.attrs)
    a["units"] = units
    if long_name is not None:
        a["long_name"] = long_name
    return da.assign_attrs(a)


def nitracline_from_gradient(grad, zmax=1000.0, depthname="depth"):
    """Depth of the maximum vertical gradient, searched over levels in the top zmax m
    (the gradient at the deepest of those levels still uses the level just below it,
    because it is computed on the full loaded column before this selection)."""
    g = grad.where(grad[depthname] <= zmax)
    return g.idxmax(dim=depthname).astype("float64")


def lat_dim_of(da):
    """Name of the dimension along which latitude varies (works on curvilinear grids)."""
    latname = da.attrs.get("lat_coord", "lat")
    latc = da.coords[latname]
    if latc.ndim == 1:
        return latc.dims[0]
    cand = [d for d in da.dims if d not in ("time", "depth") and d in latc.dims]
    return max(cand, key=lambda d: float(latc.std(d).mean()))


def lon_name_of(da):
    return da.attrs.get("lon_coord", "lon")


# ============================================================================
# 1. WOA23 nitrate reference: annual climatology (all time periods) + masks
# ============================================================================
woa_no3 = load_climatology_data(dataset="woa23", variable="nitrate",
                                time_slice="annual", depth_range=[0, 1200])
woa_levels, woa_bounds = load_depth_levels(dataset="woa23", variable="nitrate",
                                           time_slice="annual", depth_range=[0, 1200])
landsea = load_landsea_mask()          # 1 = land, ocean is (> 1)
woa_areas = load_areas()               # 1-degree cell areas

# common concentration unit for everything that follows: umol kg-1
woa_no3 = transform_unit(woa_no3, from_unit=woa_no3.attrs["units"],
                         to_unit="umol/kg").compute()

# vertical gradient of nitrate on the WOA grid, then the nitracline
grad_woa = compute_vertical_gradient(woa_no3).compute()
nc_woa = nitracline_from_gradient(grad_woa)
nc_woa = _retag(nc_woa, woa_no3, "m", "nitracline depth")
nc_woa = nc_woa.where(landsea > 1).compute()          # ocean mask

# ============================================================================
# 2. CMIP6 no3, 1971-2000, native grids -> 1971-2000 climatology -> nitracline
# ============================================================================
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
clim, nc_mod = {}, {}

for m in models:
    raw = load_regular_data(dataset="cmip6_historical", variable="no3", model=m,
                            time_range=["1971-01", "2000-12"],
                            depth_range=[0, 1200])
    c = temporal_aggregate(raw, freq="mean", op="mean")           # day-weighted 1971-2000 mean
    c = transform_unit(c, from_unit=c.attrs["units"], to_unit="umol/kg").compute()
    clim[m] = c

    g = compute_vertical_gradient(c).compute()
    nc = nitracline_from_gradient(g)
    # model ocean mask: full nitrate column present within the top 1000 m
    valid = c.sel(depth=slice(0, 1000)).notnull().all("depth")
    nc = _retag(nc, c, "m", "nitracline depth").where(valid).compute()
    nc_mod[m] = nc

# ============================================================================
# 3. Global area-weighted mean nitracline depth (native grids)
# ============================================================================
glob = {"WOA23": _val(compute_basic_stat(nc_woa, op="area_weighted_mean", dim=None))}
for m in models:
    glob[m] = _val(compute_basic_stat(nc_mod[m], op="area_weighted_mean", dim=None))

glob_diff = {m: glob[m] - glob["WOA23"] for m in models}

# ============================================================================
# 4. Regrid model nitraclines onto the WOA23 1-degree grid (global fields first)
# ============================================================================
reg = {}
for m in models:
    reg[m] = regrid_horizontal(nc_mod[m], target="woa23", method="conservative").compute()

# ============================================================================
# 5. Basin means (WOA basin mask, marginal seas folded in, Southern/Arctic excluded)
# ============================================================================
basin_names = ["Pacific", "Atlantic", "Indian"]

def basin_mean_1deg(field, basin):
    f = select_basin(field, basin=basin, include_marginal_seas=True, depth_dependent=False)
    return _val(compute_basic_stat(f, op="area_weighted_mean", dim=["lat", "lon"]))

basin = {}
for b in basin_names:
    basin[("WOA23", b)] = basin_mean_1deg(nc_woa, b)
    for m in models:
        basin[(m, b)] = basin_mean_1deg(reg[m], b)

basin_diff = {(m, b): basin[(m, b)] - basin[("WOA23", b)] for m in models for b in basin_names}

# basin ordering (shallowest -> deepest) for each source
def _order(src):
    vals = {b: basin[(src, b)] for b in basin_names}
    return sorted(basin_names, key=lambda b: vals[b])

basin_order = {src: _order(src) for src in ["WOA23"] + models}
order_match = {m: (basin_order[m] == basin_order["WOA23"]) for m in models}

# ============================================================================
# 6. Common valid ocean cells on the WOA grid -> spatial correlation & RMSE
# ============================================================================
cmp_mask = build_comparison_mask(
    fields=[nc_woa] + [reg[m] for m in models],
    ocean_only=True,
    extra_mask=(landsea > 1),
    reduce_other_dims="all",
)

skill = {}
for m in models:
    a = reg[m].where(cmp_mask)
    b = nc_woa.where(cmp_mask)
    cor = _val(analyze_correlation(var1=a, var2=b, dim=["lat", "lon"],
                                   method="pearson")["correlation"])
    met = compute_difference_metrics(data=a, reference=b, dim=["lat", "lon"])
    skill[m] = {
        "correlation": cor,
        "rmse": _val(met["rmse"]),
        "bias": _val(met["bias"]),
        "centered_rmse": _val(met["centered_rmse"]),
    }

# best model: smallest RMSE (ties broken by highest correlation)
best_model = sorted(models, key=lambda m: (round(skill[m]["rmse"], 6),
                                           -skill[m]["correlation"]))[0]

# ============================================================================
# 7. Report
# ============================================================================
print("=" * 78)
print("GLOBAL MEAN NITRACLINE DEPTH (m)")
print("=" * 78)
print(f"{'WOA23':<16}{glob['WOA23']:8.1f}  -> {round(glob['WOA23']):5d} m")
for m in models:
    print(f"{m:<16}{glob[m]:8.1f}  -> {round(glob[m]):5d} m   "
          f"diff vs WOA23 = {glob_diff[m]:+7.1f} m ({round(glob_diff[m]):+d} m)")

print("\n" + "=" * 78)
print("BASIN MEAN NITRACLINE DEPTH (m)   [marginal seas folded, S & Arctic excluded]")
print("=" * 78)
hdr = f"{'':<16}" + "".join(f"{b:>12}" for b in basin_names)
print(hdr)
print(f"{'WOA23':<16}" + "".join(f"{round(basin[('WOA23', b)]):>12d}" for b in basin_names))
for m in models:
    print(f"{m:<16}" + "".join(f"{round(basin[(m, b)]):>12d}" for b in basin_names))
print("-" * len(hdr))
print("Differences vs WOA23 (m):")
for m in models:
    print(f"{m:<16}" + "".join(f"{round(basin_diff[(m, b)]):>+12d}" for b in basin_names))

print("\nBasin ordering (shallowest -> deepest):")
for src in ["WOA23"] + models:
    print(f"  {src:<16}{' < '.join(basin_order[src])}   "
          f"({', '.join(f'{b}={basin[(src, b)]:.0f}m' for b in basin_order[src])})")
for m in models:
    print(f"  {m:<16}matches WOA23 ordering: {order_match[m]}")

print("\n" + "=" * 78)
print("SPATIAL SKILL vs WOA23 (WOA23 1-degree grid, common ocean cells)")
print("=" * 78)
for m in models:
    print(f"{m:<16} r = {skill[m]['correlation']:.3f}   "
          f"RMSE = {skill[m]['rmse']:.1f} m   bias = {skill[m]['bias']:+.1f} m   "
          f"centered RMSE = {skill[m]['centered_rmse']:.1f} m")
print(f"\nBest model reproducing the WOA23 nitracline: {best_model}")

# ============================================================================
# 8. Figure 2 data: equatorial Pacific (2S-2N, 120E-280E) longitude-depth sections
# ============================================================================
box_lat = [-2.0, 2.0]
box_lon = [120.0, 280.0]     # crosses the dateline; loader returns a continuous axis

woa_sec_raw = load_climatology_data(dataset="woa23", variable="nitrate",
                                    time_slice="annual", lat_range=box_lat,
                                    lon_range=box_lon, depth_range=[0, 1200])
woa_sec_raw = transform_unit(woa_sec_raw, from_unit=woa_sec_raw.attrs["units"],
                             to_unit="umol/kg").compute()

sections = {}
for m in models:
    d = load_regular_data(dataset="cmip6_historical", variable="no3", model=m,
                          time_range=["1971-01", "2000-12"], depth_range=[0, 1200],
                          lat_range=box_lat, lon_range=box_lon)
    d = temporal_aggregate(d, freq="mean", op="mean")
    d = transform_unit(d, from_unit=d.attrs["units"], to_unit="umol/kg").compute()
    sections[m] = d

sec_fields = {}
sec_nc = {}
for name, d in [("WOA23", woa_sec_raw)] + [(m, sections[m]) for m in models]:
    ldim = lat_dim_of(d)
    s = d.mean(dim=ldim)                       # average over 2S-2N -> longitude-depth section
    lname = lon_name_of(d)
    if s[lname].ndim == 2:
        s = s.assign_coords({lname: s[lname].mean(dim=ldim)})
    s = s.drop_vars([c for c in s.coords if c not in (lname, "depth") and s[c].ndim > 1],
                    errors="ignore")
    if lname not in s.dims:
        coldim = [x for x in s.dims if x != "depth"][0]
        s = s.swap_dims({coldim: lname})
    s = s.assign_coords({lname: s[lname] % 360}).sortby(lname)
    s = s.compute()
    sec_fields[name] = s

    g_sec = compute_vertical_gradient(s).compute()
    sec_nc[name] = nitracline_from_gradient(g_sec).compute()

panel_order = ["WOA23"] + models
print("\nSection nitracline (equatorial Pacific, 2S-2N), range over 120E-280E:")
for name in panel_order:
    v = sec_nc[name].values
    v = v[np.isfinite(v)]
    print(f"  {name:<16}min {v.min():6.0f} m   max {v.max():6.0f} m   mean {v.mean():6.0f} m")
