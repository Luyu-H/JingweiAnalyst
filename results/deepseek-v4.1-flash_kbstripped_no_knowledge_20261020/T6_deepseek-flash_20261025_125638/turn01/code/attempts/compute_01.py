# Load WOA23 nitrate annual climatology and auxiliary fields
woa_no3 = load_climatology_data(
    dataset="woa23",
    variable="nitrate",
    time_slice="annual",
    depth_range=[0, 1500],
)
woa_landsea = load_landsea_mask()
woa_basin = load_basin_mask(depth_levels="0m")
woa_levels, woa_bounds = load_depth_levels(
    dataset="woa23", variable="nitrate", time_slice="annual", depth_range=[0, 1500]
)
woa_areas = load_areas(dataset="woa23")

# WOA23 nitracline depth
woa_nitri = compute_nutricline_depth(woa_no3, min_depth=0, max_depth=1000)
woa_nitri.attrs.update(woa_no3.attrs)
woa_nitri.attrs["units"] = "m"
woa_nitri.attrs["long_name"] = "Nitracline depth"

models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

model_no3_clim = {}
model_nitri_native = {}
model_no3_regrid = {}
model_nitri_regrid = {}

for m in models:
    da = load_regular_data(
        dataset="cmip6_historical",
        variable="no3",
        model=m,
        time_range=["1971-01", "2000-12"],
        depth_range=[0, 1500],
    )

    # 1971-2000 day-weighted long-term mean
    clim = temporal_aggregate(da, freq="mean", op="mean")
    clim.attrs.update(da.attrs)
    if "units" in da.attrs:
        clim.attrs["units"] = da.attrs["units"]

    # Nitracline depth on native grid
    n = compute_nutricline_depth(clim, min_depth=0, max_depth=1000)
    n.attrs.update(clim.attrs)
    n.attrs["units"] = "m"
    n.attrs["long_name"] = "Nitracline depth"

    model_no3_clim[m] = clim
    model_nitri_native[m] = n

    # Put both the climatological nitrate field and the nitracline map on the WOA23 grid
    model_no3_regrid[m] = regrid_horizontal(clim, target="woa23", method="conservative")
    model_nitri_regrid[m] = regrid_horizontal(n, target="woa23", method="conservative")

# Common comparison mask
fields_for_mask = [woa_nitri] + [model_nitri_regrid[m] for m in models]
mask = build_comparison_mask(fields=fields_for_mask, ocean_only=True).compute()

# All nitracline fields on the common WOA23 1-degree grid
nitri_fields = {"WOA23": woa_nitri}
for m in models:
    nitri_fields[m] = model_nitri_regrid[m]


def val(x):
    if hasattr(x, "compute"):
        x = x.compute()
    return float(x)


# Global area-weighted mean nitracline depth
global_mean = {}
for src, da in nitri_fields.items():
    da_masked = da.where(mask)
    global_mean[src] = val(
        compute_basic_stat(
            da_masked, op="area_weighted_mean", dim=["lat", "lon"]
        )
    )

# Basin means, with marginal seas folded into their parent ocean.
# Southern and Arctic are not selected, so they are excluded.
basins = ["Pacific", "Atlantic", "Indian"]
basin_mean = {src: {} for src in nitri_fields}

for src, da in nitri_fields.items():
    for b in basins:
        bda = select_basin(
            da, basin=b, include_marginal_seas=True, depth_dependent=False
        )
        bda = bda.where(mask)
        basin_mean[src][b] = val(
            compute_basic_stat(bda, op="area_weighted_mean", dim=["lat", "lon"])
        )

# Differences from WOA23
global_diff = {m: global_mean[m] - global_mean["WOA23"] for m in models}
basin_diff = {
    m: {b: basin_mean[m][b] - basin_mean["WOA23"][b] for b in basins}
    for m in models
}

# Spatial correlation and RMSE against WOA23 on the common mask
corr = {}
rmse = {}

for m in models:
    a = nitri_fields[m].where(mask)
    b = nitri_fields["WOA23"].where(mask)

    corr_res = analyze_correlation(
        var1=a,
        var2=b,
        dim=["lat", "lon"],
        method="pearson",
        weighting="area",
    )
    corr[m] = val(corr_res["correlation"])

    met = compute_difference_metrics(
        data=a,
        reference=b,
        dim=["lat", "lon"],
        weighting="area",
    )
    rmse[m] = val(met["rmse"])

best_model = min(models, key=lambda m: rmse[m])

woa_order = sorted(basins, key=lambda b: basin_mean["WOA23"][b])
basin_order_match = {}
for m in models:
    model_order = sorted(basins, key=lambda b: basin_mean[m][b])
    basin_order_match[m] = model_order == woa_order

# Figure 2 data: common nitrate unit (mol m-3) and equatorial Pacific sections
def to_mol_m3(da):
    u = da.attrs.get("units", "")
    key = u.lower().replace(" ", "")
    if key in ["molm-3", "mol/m3", "molm^-3"]:
        return da
    return transform_unit(da, from_unit=u, to_unit="mol m-3", density=1.025)


woa_no3_mol = to_mol_m3(woa_no3)
model_no3_mol_regrid = {m: to_mol_m3(model_no3_regrid[m]) for m in models}

section_sources = ["WOA23"] + models
section_fields = {}
section_nitri = {}

for src in section_sources:
    da = woa_no3_mol if src == "WOA23" else model_no3_mol_regrid[src]

    # Average over 2S-2N, then take 120E-80W on a continuous 0-360 longitude axis
    sec = da.sel(lat=slice(-2, 2)).mean("lat")
    sec = sec.assign_coords(lon=(sec.lon % 360)).sortby("lon")
    sec = sec.sel(lon=slice(120, 280))
    sec.attrs.update(da.attrs)
    sec = sec.compute()

    line = compute_nutricline_depth(sec, min_depth=0, max_depth=1000)
    line.attrs.update(sec.attrs)
    line.attrs["units"] = "m"
    line = line.compute()

    section_fields[src] = sec
    section_nitri[src] = line

# Print the requested summary
print("Global mean nitracline depth (m):")
for src in ["WOA23"] + models:
    print(f"  {src}: {global_mean[src]:.0f}")

print("Model global difference from WOA23 (m):")
for m in models:
    print(f"  {m}: {global_diff[m]:+.0f}")

print("Basin mean nitracline depth (m):")
for b in basins:
    print(f"  {b}:")
    for src in ["WOA23"] + models:
        print(f"    {src}: {basin_mean[src][b]:.0f}")
    for m in models:
        print(f"    {m} diff: {basin_diff[m][b]:+.0f}")

print("Correlation and RMSE vs WOA23:")
for m in models:
    print(f"  {m}: r={corr[m]:.3f}, RMSE={rmse[m]:.1f} m")

print(f"Best model by RMSE: {best_model}")
print(f"WOA23 basin order: {woa_order}")
for m in models:
    print(f"  {m} order match: {basin_order_match[m]}")
