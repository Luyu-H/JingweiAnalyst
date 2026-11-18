import numpy as np
import xarray as xr

models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
months = {3: "March", 9: "September"}
time_range = ["1971-01", "2000-12"]


def loc_of_max(da):
    """Return (lat, lon) of the largest absolute value, on the native grid."""
    ad = np.abs(da).compute()
    latn = da.attrs.get("lat_coord", "lat")
    lonn = da.attrs.get("lon_coord", "lon")
    vals = ad.values.astype(float).ravel()
    k = int(np.nanargmax(vals))
    lat2 = da[latn].broadcast_like(ad).values.astype(float).ravel()
    lon2 = da[lonn].broadcast_like(ad).values.astype(float).ravel()
    return float(lat2[k]), float(lon2[k])


results = {}
fig1_maps = {}
scatter = {}

for m in models:
    thetao = load_regular_data(dataset="cmip6_historical", variable="thetao",
                               model=m, time_range=time_range)
    so = load_regular_data(dataset="cmip6_historical", variable="so",
                           model=m, time_range=time_range)
    mlotst = load_regular_data(dataset="cmip6_historical", variable="mlotst",
                               model=m, time_range=time_range)
    areas = load_areas(dataset="cmip6_historical", model=m, variable="thetao")
    total_ocean_area = float(areas.sum(skipna=True).compute())

    for mon, mon_name in months.items():
        t_sel = thetao.sel(time=(thetao.time.dt.month == mon))
        s_sel = so.sel(time=(so.time.dt.month == mon))
        ml_sel = mlotst.sel(time=(mlotst.time.dt.month == mon))

        # Offline operator MLD: 0.03 kg m-3 potential-density threshold, 10 m reference
        mld_op = compute_mixed_layer_depth(
            temp_data=t_sel, sal_data=s_sel,
            method="density", density_threshold=0.03, ref_depth=10.0,
            temperature_type="potential",
        )
        mld_op, ml_sel = xr.align(mld_op, ml_sel, join="inner")

        # Cell is kept only if both fields are defined at every one of the 30 instances
        mask = (mld_op.notnull() & ml_sel.notnull()).all("time").compute()

        op_clim = temporal_aggregate(mld_op.where(mask), freq="mean").compute()
        ml_clim = temporal_aggregate(ml_sel.where(mask), freq="mean").compute()
        # carry provenance attrs so the weighting tools can find the native grid
        op_clim = op_clim.assign_attrs(**{**op_clim.attrs, **t_sel.attrs})
        ml_clim = ml_clim.assign_attrs(**{**ml_clim.attrs, **ml_sel.attrs})

        h_dims = list(op_clim.dims)

        mets = compute_difference_metrics(data=op_clim, reference=ml_clim,
                                          dim=h_dims, weighting="area")
        corr_res = analyze_correlation(var1=op_clim, var2=ml_clim,
                                       dim=h_dims, method="pearson",
                                       weighting="area")

        diff = (op_clim - ml_clim).assign_attrs(units="m")
        ex = compute_exceedance_fraction(data=diff, threshold=20.0,
                                         comparison="abs_greater",
                                         dim=h_dims, weighting="area")

        valid_area = float((mask.astype(float) * areas).sum(skipna=True).compute())
        valid_pct = 100.0 * valid_area / total_ocean_area

        lat_max, lon_max = loc_of_max(diff)

        results[(m, mon_name)] = dict(
            bias=float(mets["bias"]),
            rmse=float(mets["rmse"]),
            corr=float(corr_res["correlation"]),
            exceed_pct=float(ex["fraction_percent"]),
            valid_pct=valid_pct,
            max_err=float(mets["max_abs_error"]),
            max_lat=lat_max,
            max_lon=lon_max,
            n_valid=int(mets["n_valid"]),
        )

        # display-only regrid of the native-grid difference map
        fig1_maps[(m, mon_name)] = regrid_horizontal(diff, target="woa23", method="bilinear")

        scatter[(m, mon_name)] = (op_clim.values.astype(float).ravel(),
                                  ml_clim.values.astype(float).ravel())

rmse_avg = {m: 0.5 * (results[(m, "March")]["rmse"] +
                      results[(m, "September")]["rmse"]) for m in models}
ranking = sorted(models, key=lambda mm: rmse_avg[mm])
ranking_table = [(mm, round(rmse_avg[mm], 1)) for mm in ranking]

summary_rows = []
for m in models:
    for mon_name in ("March", "September"):
        r = results[(m, mon_name)]
        summary_rows.append(dict(
            model=m, month=mon_name,
            bias_m=round(r["bias"], 1),
            rmse_m=round(r["rmse"], 1),
            corr=round(r["corr"], 3),
            exceed20_pct=round(r["exceed_pct"], 1),
            valid_area_pct=round(r["valid_pct"], 1),
            max_abs_err_m=round(r["max_err"], 1),
            max_err_latlon=(round(r["max_lat"], 1), round(r["max_lon"], 1)),
            n_valid_cells=r["n_valid"],
        ))
summary_table = pd.DataFrame(summary_rows)

print(summary_table.to_string(index=False))
print("\nRanking (mean of March & September RMSE, lower = more consistent):")
for i, (mm, v) in enumerate(ranking_table, 1):
    print(f"  {i}. {mm}: {v} m")
