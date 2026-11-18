models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
time_range = ["1971-01", "2000-12"]

# ---- 1. Load and build day-weighted 1971-2000 SST climatologies ----
clim = {}
native_mean = {}
fill_check = {}
area_source = {}
for m in models:
    tos = load_regular_data(dataset="cmip6_historical", variable="tos", model=m,
                            time_range=time_range)
    c = temporal_aggregate(tos, freq="mean", op="mean").compute()
    fill_check[m] = float(np.nanmax(np.abs(c.values)))   # must be far below 1e20
    assert fill_check[m] < 1e10, f"fill values leaked into {m}"
    clim[m] = c
    gm = compute_basic_stat(c, op="area_weighted_mean", dim=["lat", "lon"])
    native_mean[m] = float(gm.compute())
    area_source[m] = gm.attrs.get("weights", gm.attrs.get("area_source", "areacello"))

# ---- WOA23 0 m annual, 1971-2000 normal ----
woa_raw = load_climatology_data(dataset="woa23", variable="temperature",
                                time_slice="annual", period="1971-2000",
                                depth_range=[0, 1])
woa = woa_raw.isel(depth=0, drop=True).compute()
assert float(np.nanmax(np.abs(woa.values))) < 1e10
native_mean["WOA23"] = float(compute_basic_stat(woa, op="area_weighted_mean",
                                                dim=["lat", "lon"]).compute())

# ---- 4. Conservative regrid to the WOA23 1-degree grid ----
reg = {}
for m in models:
    reg[m] = regrid_horizontal(clim[m], target="woa23", method="conservative").compute()

# ---- 5. One common mask ----
mask = build_comparison_mask([reg[m] for m in models] + [woa], ocean_only=True)
mask = mask.compute()
n_mask_cells = int(mask.sum())
mask_attrs = dict(mask.attrs)

reg_m = {m: reg[m].where(mask) for m in models}
woa_m = woa.where(mask)

# ---- 6. Scores on shared sample ----
scores = {}
diff_maps = {}
for m in models:
    met = compute_difference_metrics(reg_m[m], woa_m)
    corr = analyze_correlation(reg_m[m], woa_m, dim=["lat", "lon"], weighting="area")
    diff = met["difference"].compute()
    diff_maps[m] = diff
    ext_pos = locate_extreme(diff, op="max")
    ext_neg = locate_extreme(diff, op="min")
    scores[m] = {
        "mean_model": float(met["mean_data"]),
        "mean_woa23": float(met["mean_reference"]),
        "bias": float(met["bias"]),
        "rmse": float(met["rmse"]),
        "centered_rmse": float(met["centered_rmse"]),
        "corr": float(corr["correlation"].compute()),
        "max_abs_error": float(met["max_abs_error"]),
        "max_abs_error_at": met["max_abs_error_at"],
        "largest_positive": (ext_pos["value"], ext_pos["coordinates"]),
        "largest_negative": (ext_neg["value"], ext_neg["coordinates"]),
        "n_valid": met["n_valid"],
    }

ranking = sorted(models, key=lambda m: abs(scores[m]["bias"]))

# ---- Report ----
print("Native-grid global means (degC): GFDL-ESM4 is a 1-deg 'gr' grid; CESM2/MPI on 'gn' grids")
for k, v in native_mean.items():
    print(f"  {k}: {v:.2f}")
print("Area sources:", area_source)
print("Max |value| check (no 1e20):", fill_check)
print(f"Common mask: {n_mask_cells} cells; attrs: {mask_attrs}")
print("Conservative regridding to WOA23 1-deg grid")
for m in models:
    s = scores[m]
    print(f"{m}: model mean {s['mean_model']:.2f}, WOA23 mean {s['mean_woa23']:.2f}, "
          f"bias {s['bias']:.2f}, RMSE {s['rmse']:.3f}, r {s['corr']:.3f}")
    print(f"   largest +: {s['largest_positive']}, largest -: {s['largest_negative']}, "
          f"abs max at {s['max_abs_error_at']}")
print("Ranking by |bias| (smallest to largest):", ranking)
