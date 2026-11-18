import numpy as np
import xarray as xr

models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
months = [3, 9]
MONTH_NAME = {3: "March", 9: "September"}


def _val(x):
    """Materialise a (possibly dask-backed) scalar/0-d object as a float."""
    if hasattr(x, "compute"):
        x = x.compute()
    return float(np.asarray(x))


def _restamp(da, base_attrs, units, long_name):
    """Keep the loader's grid provenance while re-declaring the physical quantity."""
    a = dict(base_attrs)
    a.update(dict(da.attrs))
    a["units"] = units
    a["long_name"] = long_name
    da.attrs = a
    return da


results = {}        # (model, month) -> dict of metric scalars
op_clim = {}        # (model, month) -> operator MLD March/Sept climatology, native grid
ml_clim = {}        # (model, month) -> published mlotst March/Sept climatology, native grid
diff_clim = {}      # (model, month) -> operator - mlotst on the common mask, native grid
diff_map = {}       # (model, month) -> same field regridded to 1 deg for display only
extreme_info = {}   # (model, month) -> locate_extreme output on the difference
coverage = {}       # (model, month) -> coverage percent

for model in models:
    th = load_regular_data(dataset="cmip6_historical", variable="thetao", model=model,
                           time_range=["1971-01", "2000-12"], chunk_time=6)
    sa = load_regular_data(dataset="cmip6_historical", variable="so", model=model,
                           time_range=["1971-01", "2000-12"], chunk_time=6)
    ml = load_regular_data(dataset="cmip6_historical", variable="mlotst", model=model,
                           time_range=["1971-01", "2000-12"], chunk_time=6)
    ar = load_regular_data(dataset="cmip6_historical", variable="areacello", model=model)

    base_attrs = dict(th.attrs)

    # the model's own ocean domain: cells where the published MLD is ever defined
    ocean_cells = ml.notnull().any("time").compute()
    total_ocean_area = float(ar.where(ocean_cells).sum().compute())

    for mm in months:
        th_m = th.sel(time=th.time.dt.month == mm)
        sa_m = sa.sel(time=sa.time.dt.month == mm)
        ml_m = ml.sel(time=ml.time.dt.month == mm)

        # offline MLD with the density criterion that defines mlotst (sigma0 offset 0.03
        # kg m-3 from the 10 m reference level, crossing interpolated between levels)
        mld = compute_mixed_layer_depth(temp_data=th_m, sal_data=sa_m,
                                        ref_depth=10.0,
                                        method="density",
                                        density_threshold=0.03).compute()

        # Common mask built directly on the shared native grid: a cell is kept only if
        # both fields are defined at every one of the 30 instances of this month.
        # (build_comparison_mask is for fields on *different* grids and refuses
        # curvilinear index dimensions, which CESM2 and MPI-ESM1-2-LR use.)
        common = ((mld.notnull()) & (ml_m.notnull())).all("time").compute()

        op_mc = mld.where(common).compute()
        ml_mc = ml_m.where(common).compute()

        oc = _restamp(temporal_aggregate(op_mc, freq="mean").compute(),
                      base_attrs, "m", "operator mixed layer depth")
        mc = _restamp(temporal_aggregate(ml_mc, freq="mean").compute(),
                      base_attrs, "m", "published mlotst")

        d = _restamp((oc - mc).compute(), base_attrs, "m",
                     "operator MLD minus mlotst")

        hdims = list(oc.dims)

        met = compute_difference_metrics(data=oc, reference=mc, dim=hdims,
                                         weighting="area")
        exc = compute_exceedance_fraction(d, threshold=20.0, comparison="abs_greater",
                                          dim=hdims, weighting="area")
        cor = analyze_correlation(var1=oc, var2=mc, dim=hdims, method="pearson")
        loc = locate_extreme(d, op="abs_max")

        cov = float((ar.where(common).sum() / ar.where(ocean_cells).sum()).compute()) * 100.0

        results[(model, mm)] = {
            "bias": _val(met["bias"]),
            "rmse": _val(met["rmse"]),
            "mae": _val(met["mae"]),
            "centered_rmse": _val(met["centered_rmse"]),
            "max_abs_error": _val(met["max_abs_error"]),
            "max_abs_error_at": met["max_abs_error_at"],
            "correlation": _val(cor["correlation"]),
            "pvalue": _val(cor["pvalue"]),
            "frac_gt20": _val(exc["fraction_percent"]),
            "n_common": int(np.asarray(exc["n_valid"])),
            "coverage_percent": cov,
        }
        coverage[(model, mm)] = cov
        extreme_info[(model, mm)] = {
            "value": _val(loc["value"]),
            "coordinates": loc["coordinates"],
            "description": loc["description"],
        }

        op_clim[(model, mm)] = oc
        ml_clim[(model, mm)] = mc
        diff_clim[(model, mm)] = d
        diff_map[(model, mm)] = regrid_horizontal(d, target="woa23",
                                                  method="bilinear").compute()

        print(f"{model:16s} {MONTH_NAME[mm]:10s} "
              f"bias={results[(model, mm)]['bias']:7.1f} m  "
              f"RMSD={results[(model, mm)]['rmse']:6.1f} m  "
              f"r={results[(model, mm)]['correlation']:6.3f}  "
              f"|d|>20m={results[(model, mm)]['frac_gt20']:5.1f}%  "
              f"coverage={cov:5.1f}%  "
              f"max|d|={results[(model, mm)]['max_abs_error']:6.1f} m "
              f"at {extreme_info[(model, mm)]['coordinates']}")

# ---------------------------------------------------------------- ranking
rank_table = []
for model in models:
    rms_mean = float(np.mean([results[(model, mm)]["rmse"] for mm in months]))
    rank_table.append((model, rms_mean))
rank_table.sort(key=lambda t: t[1])
ranking = [name for name, _ in rank_table]

print("\nRanking (most -> least internally consistent), criterion = mean of the "
      "March and September area-weighted RMSD over the common mask:")
for i, (name, rms) in enumerate(rank_table, start=1):
    print(f"  {i}. {name:16s} mean RMSD = {rms:6.2f} m")

summary_table = pd.DataFrame(
    [{"model": m, "month": MONTH_NAME[mm],
      "bias_m": results[(m, mm)]["bias"],
      "rmsd_m": results[(m, mm)]["rmse"],
      "correlation": results[(m, mm)]["correlation"],
      "pct_gt_20m": results[(m, mm)]["frac_gt20"],
      "coverage_pct": results[(m, mm)]["coverage_percent"]}
     for m in models for mm in months])
print("\n" + summary_table.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
