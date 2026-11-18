ds = "cmip6_historical"
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
month_names = {3: "March", 9: "September"}
time_range = ["1971-01", "2000-12"]
basins = ["Atlantic", "Pacific", "Indian", "Southern", "Arctic"]
BIG = 1e10  # anything above this is an unmasked fill value (about 1e20)

results = {}        # (model, month name) -> dict of reported numbers
diff_native = {}    # native-grid difference fields (shared mask applied)
diff_map = {}       # difference regridded to a regular grid, for display only
scatter = {}        # (model, month name) -> (x, y) arrays, one point per shared cell
model_info = {}     # per-model reference depth and total ocean area
basin_tables = {}

for model in models:
    levels, bounds = load_depth_levels(dataset=ds, model=model, variable="thetao")
    ref_depth = float(np.asarray(levels)[0])  # shallowest level of this model

    thetao = load_regular_data(dataset=ds, variable="thetao", model=model, time_range=time_range)
    so = load_regular_data(dataset=ds, variable="so", model=model, time_range=time_range)
    mlotst = load_regular_data(dataset=ds, variable="mlotst", model=model, time_range=time_range)
    if "depth" in mlotst.dims:
        mlotst = mlotst.isel(depth=0, drop=True)

    # put so on thetao's vertical grid only if the two differ
    if (so.sizes["depth"] != thetao.sizes["depth"]) or (
            not np.allclose(so["depth"].values, thetao["depth"].values)):
        so = so.interp(depth=thetao["depth"])

    # the model's own areacello, on the model's own grid
    area = load_areas(dataset=ds, model=model, variable="thetao")
    total_area = float(area.sum().compute())
    model_info[model] = {"ref_depth_m": ref_depth, "total_ocean_area_m2": total_area,
                         "area_source": area.attrs.get("area_source")}

    for mnum, mname in month_names.items():
        idx = np.flatnonzero(thetao["time"].dt.month.values == mnum)
        t_m = thetao.isel(time=idx)
        s_m = so.isel(time=idx)
        idx2 = np.flatnonzero(mlotst["time"].dt.month.values == mnum)
        mlot_m = mlotst.isel(time=idx2)
        assert len(idx) == 30 and len(idx2) == 30

        # MLD of each monthly-mean profile (density criterion 0.03 kg m-3, reference depth = shallowest level)
        mld = compute_mixed_layer_depth(t_m, s_m, ref_depth=ref_depth, method="density",
                                        density_threshold=0.03, temperature_type="auto").compute()
        mld = mld.where(mld < BIG)
        mlot_m = mlot_m.where(mlot_m < BIG).compute()
        mlot_m = mlot_m.assign_coords(time=mld["time"])
        mld = mld.assign_attrs({**mlot_m.attrs, "units": "m", "long_name": "operator MLD (sigma0 +0.03)"})

        # area removed by excluding no-crossing columns (ocean at the surface in all 30 months, but no crossing found)
        ocean_surf = t_m.isel(depth=0).notnull().all("time").compute()
        nocross = ocean_surf & mld.isnull().any("time")
        nocross_area = float(area.where(nocross).sum().compute())

        # climatologies: MLD of each monthly profile first, then the mean over the 30 instances
        mld_clim = temporal_aggregate(mld, freq="mean").compute()
        mlot_clim = temporal_aggregate(mlot_m, freq="mean").compute()

        # shared mask: both fields valid in all 30 instances (same grid per model, so no regridding needed)
        mask = (mld.notnull() & mlot_m.notnull()).all("time")
        a = mld_clim.where(mask)
        b = mlot_clim.where(mask)
        mask_cells = int(mask.sum().compute())
        mask_area = float(area.where(mask).sum().compute())

        met = compute_difference_metrics(a, b)  # operator MLD minus mlotst, area-weighted by own areacello
        diff = met["difference"]
        corr = analyze_correlation(a, b, dim=["lat", "lon"], weighting="area")
        exc = compute_exceedance_fraction(diff, threshold=20.0, comparison="abs_greater", dim=["lat", "lon"])
        ext = locate_extreme(diff, op="abs_max")

        # basin breakdown of the disagreement
        rows = []
        for bname in basins:
            dB = select_basin(diff, basin=bname)
            if int(dB.notnull().sum().compute()) == 0:
                continue
            mB = float(compute_basic_stat(dB, op="area_weighted_mean", dim=["lat", "lon"]).compute())
            r2 = float(compute_basic_stat(dB ** 2, op="area_weighted_mean", dim=["lat", "lon"]).compute())
            eB = compute_exceedance_fraction(dB, threshold=20.0, comparison="abs_greater", dim=["lat", "lon"])
            rows.append({"basin": bname, "mean_diff_m": mB, "rmse_m": float(np.sqrt(r2)),
                         "pct_area_abs_diff_gt_20m": float(eB["fraction_percent"]),
                         "area_Mkm2": float(eB["total_amount"]) / 1e12})
        basin_tables[(model, mname)] = pd.DataFrame(rows)

        key = (model, mname)
        results[key] = {
            "mean_diff_m": float(met["bias"]),
            "rmse_m": float(met["rmse"]),
            "centered_rmse_m": float(met["centered_rmse"]),
            "mean_operator_m": float(met["mean_data"]),
            "mean_mlotst_m": float(met["mean_reference"]),
            "corr": float(corr["correlation"]),
            "pct_abs_diff_gt_20m": float(exc["fraction_percent"]),
            "mask_cells": mask_cells,
            "mask_area_m2": mask_area,
            "mask_pct_of_ocean_area": 100.0 * mask_area / total_area,
            "nocross_area_removed_m2": nocross_area,
            "nocross_area_removed_pct_of_ocean": 100.0 * nocross_area / total_area,
            "std_ratio": float(met["std_ratio"]),
            "max_abs_diff_m": ext["value"],
            "max_abs_diff_where": ext["coordinates"],
            "n_valid": met["n_valid"],
        }
        diff_native[key] = diff

        # regrid for display only (all numbers above remain on the native grid)
        diff_map[key] = regrid_horizontal(diff, target="woa23", method="conservative").compute()

        # scatter data: one point per shared-mask cell; log axes need positive values
        xv = a.values.ravel()
        yv = b.values.ravel()
        ok = np.isfinite(xv) & np.isfinite(yv)
        pos = ok & (xv > 0) & (yv > 0)
        scatter[key] = (xv[pos], yv[pos])
        results[key]["n_points_scatter"] = int(pos.sum())
        results[key]["n_points_nonpositive_dropped_from_scatter"] = int(ok.sum() - pos.sum())

    del thetao, so, mlotst

summary = pd.DataFrame(results).T
summary.index.names = ["model", "month"]

# ranking: mean RMSE over the two months (criterion), correlation and >20 m share as supporting criteria
rank_tbl = summary.reset_index().groupby("model")[["rmse_m", "corr", "pct_abs_diff_gt_20m", "mean_diff_m"]].mean()
rank_tbl = rank_tbl.sort_values("rmse_m")
rank_tbl["rank"] = np.arange(1, len(rank_tbl) + 1)
ranking_criterion = ("Lowest mean (March, September) area-weighted RMSE of operator MLD minus mlotst on the shared mask; "
                     "higher correlation and a smaller >20 m area share are supporting criteria.")

caveats = ("Operator MLD is computed from monthly-mean thetao and so, whereas mlotst is the time mean of instantaneous "
           "diagnostics; the MLD criterion is nonlinear, so exact agreement is not expected. Other causes: vertical "
           "interpolation and reference depth, equation-of-state differences, and differing encodings of no-crossing columns.")

print(pd.DataFrame(model_info).T)
print(summary.round(3).to_string())
print(rank_tbl.round(3))
print(ranking_criterion)
print(caveats)
for k, v in basin_tables.items():
    print(k)
    print(v.round(2).to_string(index=False))
