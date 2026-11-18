# ---------------------------------------------------------------------------
# Internal consistency of CMIP6 mlotst vs an offline MLD diagnostic
# ---------------------------------------------------------------------------
models = ['CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']
month_names = ['March', 'September']
month_nums = {'March': 3, 'September': 9}
time_range = ['1971-01', '2000-12']


def _num(x):
    """Materialise a scalar (float / 0-d numpy- or dask-backed DataArray)."""
    if hasattr(x, 'compute'):
        x = x.compute()
    return float(x)


results = {}        # results[(model, month)] -> dict of diagnostics
map_fields = {}     # map_fields[(model, month)] -> regridded difference field (display only)
scatter_vals = {}   # scatter_vals[(model, month)] -> (operator MLD, mlotst) on common cells

for model in models:
    print('=' * 78)
    print('Model:', model)

    thetao = load_regular_data(dataset='cmip6_historical', variable='thetao', model=model,
                               field='total', time_range=time_range,
                               depth_range=[0.0, 1000.0])
    so = load_regular_data(dataset='cmip6_historical', variable='so', model=model,
                           field='total', time_range=time_range,
                           depth_range=[0.0, 1000.0])
    mlotst = load_regular_data(dataset='cmip6_historical', variable='mlotst', model=model,
                               field='total', time_range=time_range)
    areacello = load_regular_data(dataset='cmip6_historical', variable='areacello', model=model)

    # --- offline MLD from the model's own T and S, density criterion (0.03 kg m-3
    #     relative to the 10 m reference level), crossing depth interpolated ------
    mld_op = compute_mixed_layer_depth(temp_data=thetao, sal_data=so,
                                       ref_depth=10.0, method='density',
                                       density_threshold=0.03)

    total_area = _num(compute_basic_stat(data=areacello, op='sum', dim=None))

    for mname in month_names:
        mnum = month_nums[mname]
        op_m = mld_op.sel(time=(mld_op.time.dt.month == mnum))
        pub_m = mlotst.sel(time=(mlotst.time.dt.month == mnum))
        n_months = int(op_m.sizes['time'])

        # both fields must be defined in EVERY instance of this calendar month
        mask = (op_m.notnull() & pub_m.notnull()).all('time')

        op_m = op_m.where(mask)
        pub_m = pub_m.where(mask)

        clim_op = temporal_aggregate(data=op_m, freq='climatology').squeeze('month', drop=True)
        clim_pub = temporal_aggregate(data=pub_m, freq='climatology').squeeze('month', drop=True)

        # keep the loader's grid provenance on the climatologies
        clim_op.attrs = {**clim_op.attrs,
                         **{k: v for k, v in mld_op.attrs.items() if k != 'units'}}
        clim_op.attrs['units'] = 'm'
        clim_pub.attrs = {**clim_pub.attrs,
                          **{k: v for k, v in mlotst.attrs.items() if k != 'units'}}
        clim_pub.attrs['units'] = 'm'

        # --- bias / RMSE (area weighted over the common valid cells) -------------
        metrics = compute_difference_metrics(data=clim_op, reference=clim_pub)
        bias = _num(metrics['bias'])
        rmse = _num(metrics['rmse'])
        crmse = _num(metrics['centered_rmse'])
        n_valid = int(metrics['n_valid'])

        # --- spatial Pearson correlation over the common cells ------------------
        corr = analyze_correlation(var1=clim_op, var2=clim_pub, dim=['lat', 'lon'])
        r = _num(corr['correlation'])
        pval = _num(corr['pvalue'])

        # --- difference field, extreme, exceedance ------------------------------
        diff_field = metrics['difference'].compute()
        diff_field.attrs = dict(metrics['difference'].attrs)
        diff_field.attrs['units'] = 'm'

        extreme = locate_extreme(data=diff_field, op='abs_max')

        abs_diff = abs(diff_field)
        abs_diff.attrs = dict(diff_field.attrs)
        abs_diff.attrs['units'] = 'm'
        exc = compute_exceedance_fraction(data=abs_diff, threshold=20.0,
                                          comparison='greater', weighting='area', dim=None)
        frac20 = _num(exc['fraction_percent'])

        common_area = _num(compute_basic_stat(data=areacello.where(mask), op='sum', dim=None))
        area_pct = 100.0 * common_area / total_area

        results[(model, mname)] = dict(bias=bias, rmse=rmse, centered_rmse=crmse,
                                       corr=r, pvalue=pval, frac20=frac20,
                                       common_area=common_area, total_area=total_area,
                                       area_pct=area_pct, extreme=extreme,
                                       n_months=n_months, n_valid=n_valid)

        # --- native-grid values for the scatter panel ---------------------------
        scatter_vals[(model, mname)] = (clim_op.compute().values.ravel(),
                                        clim_pub.compute().values.ravel())

        # --- regrid the difference for DISPLAY only -----------------------------
        map_fields[(model, mname)] = regrid_horizontal(data=diff_field, target='1deg',
                                                       method='bilinear').compute()

        print(f'  {mname:9s} ({n_months} fields)  bias={bias:+7.1f} m  rmse={rmse:6.1f} m  '
              f'r={r:+.3f}  |d|>20 m: {frac20:5.1f} %  common area: {area_pct:5.1f} %')
        print('           extreme:', extreme['description'])

    del thetao, so, mlotst, areacello, mld_op

# ---------------------------------------------------------------------------
# Ranking: most -> least internally consistent (mean of March & September RMSE)
# ---------------------------------------------------------------------------
model_score = {}
for model in models:
    rmses = [results[(model, mn)]['rmse'] for mn in month_names]
    fracs = [results[(model, mn)]['frac20'] for mn in month_names]
    model_score[model] = (float(np.mean(rmses)), float(np.mean(fracs)))

ranking = sorted(models, key=lambda m: (model_score[m][0], model_score[m][1]))

print()
print('=' * 78)
print('RANKING (most -> least internally consistent)')
print('primary criterion: mean of the March and September RMSE (operator MLD - mlotst)')
for i, m in enumerate(ranking, start=1):
    mr, mf = model_score[m]
    print(f'  {i}. {m:16s} mean RMSE = {mr:6.1f} m | mean |d|>20 m fraction = {mf:5.1f} %')

print()
print('=' * 78)
print('SUMMARY TABLE')
for model in models:
    for mname in month_names:
        d = results[(model, mname)]
        print(f'{model:16s} {mname:9s} '
              f'bias={d["bias"]:+7.1f}  rmse={d["rmse"]:6.1f}  r={d["corr"]:+.3f}  '
              f'|d|>20m={d["frac20"]:5.1f}%  valid area={d["area_pct"]:5.1f}%')

summary_table = pd.DataFrame(
    [{'model': m, 'month': mn,
      'bias_m': round(results[(m, mn)]['bias'], 1),
      'rmse_m': round(results[(m, mn)]['rmse'], 1),
      'corr': round(results[(m, mn)]['corr'], 3),
      'pct_absdiff_gt20m': round(results[(m, mn)]['frac20'], 1),
      'pct_of_ocean_area': round(results[(m, mn)]['area_pct'], 1)}
     for m in models for mn in month_names])
print(summary_table.to_string(index=False))
