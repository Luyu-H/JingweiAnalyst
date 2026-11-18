# ---------------------------------------------------------------------------
# Internal consistency test of CMIP6 mlotst against an offline (operator) MLD
# computed from each model's own thetao / so, 1971-2000, March & September.
# ---------------------------------------------------------------------------
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
time_range = ["1971-01", "2000-12"]
month_defs = [("March", 3), ("September", 9)]


def _hdims(da):
    """Horizontal dimensions of a field (native names: lat/lon or nlat/nlon or j/i)."""
    return [d for d in da.dims if d not in ("time", "depth", "month", "season", "year")]


results = {}   # (model, month_name) -> dict of headline numbers
panels  = {}   # (model, month_name) -> fields for plotting
ref_depth_info = {}

for model in models:
    print("=" * 78)
    print(model)

    thetao = load_regular_data(dataset="cmip6_historical", variable="thetao", model=model,
                               time_range=time_range, chunk_time=6)
    so     = load_regular_data(dataset="cmip6_historical", variable="so", model=model,
                               time_range=time_range, chunk_time=6)
    mlotst = load_regular_data(dataset="cmip6_historical", variable="mlotst", model=model,
                               time_range=time_range, chunk_time=6)
    area   = load_areas(dataset="cmip6_historical", model=model, variable="thetao")
    levels, bounds = load_depth_levels(dataset="cmip6_historical", model=model, variable="thetao")

    lev = np.asarray(levels).ravel()
    ref_depth = float(lev[0])
    hd = _hdims(mlotst)
    ref_depth_info[model] = dict(ref_depth=ref_depth, n_levels=len(lev),
                                 top=float(lev[0]), bottom=float(lev[-1]))
    print("  horizontal dims :", hd)
    print("  ref_depth (=shallowest level centre): %.3f m ; levels %d (%.2f .. %.2f m)"
          % (ref_depth, len(lev), lev[0], lev[-1]))

    # ---- offline mixed layer depth, one per month (never average T/S first) ----
    mld = compute_mixed_layer_depth(
        temp_data=thetao, sal_data=so,
        ref_depth=ref_depth,
        method="density",
        density_threshold=0.03,
        temperature_type="potential",
    )
    mld = mld.compute()
    mld.attrs.update(thetao.attrs)      # keep grid provenance for the weighting helpers
    mld.attrs["units"] = "m"
    mld.attrs["long_name"] = "offline mixed layer depth (sigma0 difference 0.03 kg m-3)"

    mlotst = mlotst.compute()
    area = area.compute()

    total_area = float(compute_basic_stat(data=area, op="sum", dim=hd).compute())
    print("  total ocean area (areacello): %.4e m2" % total_area)

    for mname, mnum in month_defs:
        mld_m = mld.sel(time=(mld.time.dt.month == mnum))
        mlo_m = mlotst.sel(time=(mlotst.time.dt.month == mnum))
        mld_m, mlo_m = xr.align(mld_m, mlo_m, join="inner")
        print("  %-9s instances: %d" % (mname, mld_m.sizes["time"]))

        # cells where BOTH fields are defined in every one of the 30 instances
        mask = (mld_m.notnull() & mlo_m.notnull()).all("time").compute()

        # plain 30-year calendar-month mean of that month's values
        mld_c = mld_m.where(mask).mean("time").compute()
        mlo_c = mlo_m.where(mask).mean("time").compute()
        mld_c.attrs = dict(thetao.attrs); mld_c.attrs["units"] = "m"
        mlo_c.attrs = dict(thetao.attrs); mlo_c.attrs["units"] = "m"

        diff = mld_c - mlo_c
        diff.attrs = dict(thetao.attrs)
        diff.attrs["units"] = "m"
        diff.attrs["long_name"] = "operator MLD minus mlotst"

        # ---- metrics on the native grid, common mask, areacello weights ----
        dm = compute_difference_metrics(data=mld_c, reference=mlo_c, dim=hd, weighting="area")
        bias = float(dm["bias"]); rmse = float(dm["rmse"])
        crmse = float(dm["centered_rmse"])
        maxerr = float(dm["max_abs_error"])

        cr = analyze_correlation(var1=mld_c, var2=mlo_c, dim=hd, method="pearson")
        r = float(cr["correlation"]); pval = float(cr["pvalue"])

        exc = compute_exceedance_fraction(data=diff, threshold=20.0,
                                          comparison="abs_greater", dim=hd, weighting="area")
        ex_pct = float(exc["fraction_percent"])
        ex_area = float(exc["exceeding_amount"])

        common_area = float(compute_basic_stat(data=area.where(mask), op="sum", dim=hd).compute())
        cov_pct = 100.0 * common_area / total_area

        # area dropped because the sigma0 threshold was never reached (operator MLD missing)
        never = (mlo_m.notnull().all("time") & mld_m.isnull().any("time")).compute()
        never_area = float(compute_basic_stat(data=area.where(never), op="sum", dim=hd).compute())

        results[(model, mname)] = dict(
            bias=bias, rmse=rmse, crmse=crmse, max_abs_error=maxerr,
            corr=r, pvalue=pval, exceed_pct=ex_pct, exceed_area=ex_area,
            common_area=common_area, total_area=total_area,
            coverage_pct=cov_pct, never_area=never_area,
            never_pct=100.0 * never_area / total_area,
            mean_mld=float(dm["mean_data"]), mean_mlotst=float(dm["mean_reference"]),
        )

        print("    bias %+7.2f m | RMSE %7.2f m | r %.3f | >20 m %5.1f %% | coverage %5.1f %% (%.3e m2) | never-reached %.2e m2"
              % (bias, rmse, r, ex_pct, cov_pct, common_area, never_area))

        # values for the scatter panel (native cells, common mask)
        xv = mld_c.values
        yv = mlo_c.values
        mv = mask.values
        panels[(model, mname)] = dict(
            diff_map=regrid_horizontal(diff, target="woa23", method="bilinear",
                                       skipna=True, na_thres=1.0).compute(),
            x=np.asarray(xv)[mv], y=np.asarray(yv)[mv],
            bias=bias, rmse=rmse, corr=r, exceed_pct=ex_pct, coverage_pct=cov_pct,
        )

# ---------------------------------------------------------------------------
# Ranking: mean of the March and September RMSDs, lower = more consistent
# ---------------------------------------------------------------------------
ranking = []
for model in models:
    rms = 0.5 * (results[(model, "March")]["rmse"] + results[(model, "September")]["rmse"])
    ranking.append((model, rms))
ranking.sort(key=lambda t: t[1])

print("\n" + "=" * 78)
print("RANKING by mean of March & September RMSD (lower = more internally consistent)")
for i, (model, rms) in enumerate(ranking, 1):
    b = 0.5 * (results[(model, "March")]["bias"] + results[(model, "September")]["bias"])
    rr = 0.5 * (results[(model, "March")]["corr"] + results[(model, "September")]["corr"])
    ex = 0.5 * (results[(model, "March")]["exceed_pct"] + results[(model, "September")]["exceed_pct"])
    cv = 0.5 * (results[(model, "March")]["coverage_pct"] + results[(model, "September")]["coverage_pct"])
    print("  %d. %-14s mean RMSD %6.2f m | mean bias %+6.2f m | mean r %.3f | mean >20 m %5.1f %% | mean coverage %5.1f %%"
          % (i, model, rms, b, rr, ex, cv))

print("\nFull headline table (operator MLD minus mlotst):")
print("%-14s %-10s %8s %8s %8s %8s %8s %10s" %
      ("model", "month", "bias", "RMSE", "cRMSE", "r", ">20m %", "cover %"))
for model in models:
    for mname, _ in month_defs:
        d = results[(model, mname)]
        print("%-14s %-10s %+8.2f %8.2f %8.2f %8.3f %8.1f %10.1f"
              % (model, mname, d["bias"], d["rmse"], d["crmse"], d["corr"],
                 d["exceed_pct"], d["coverage_pct"]))

print("\nReference (shallowest) level used per model:", ref_depth_info)

# non-positive MLD cells (cannot be shown on logarithmic axes)
npos = {}
for model in models:
    for mname, _ in month_defs:
        x = panels[(model, mname)]["x"]; y = panels[(model, mname)]["y"]
        bad = int(np.sum(~(np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0))))
        npos[(model, mname)] = bad
print("Non-positive / non-finite cells excluded from the log-log scatter (not from the metrics):", npos)
