# Helper to materialise scalars/DataArrays from lazy results.
def _val(x):
    if hasattr(x, 'compute'):
        x = x.compute()
    if hasattr(x, 'item'):
        return float(x.item())
    return float(x)

models = ['CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']
grid_labels = {'CESM2': 'gn', 'GFDL-ESM4': 'gr', 'MPI-ESM1-2-LR': 'gn'}
time_range = ['1971-01', '2000-12']
months = [3, 9]
month_names = {3: 'March', 9: 'September'}

results = []
store = {}

for model in models:
    gl = grid_labels[model]

    thetao = load_regular_data(
        dataset='cmip6_historical', variable='thetao', model=model,
        time_range=time_range
    )
    so = load_regular_data(
        dataset='cmip6_historical', variable='so', model=model,
        time_range=time_range
    )
    mlotst = load_regular_data(
        dataset='cmip6_historical', variable='mlotst', model=model,
        time_range=time_range
    )

    thetao.attrs['grid_label'] = gl
    so.attrs['grid_label'] = gl
    mlotst.attrs['grid_label'] = gl

    levels, bounds = load_depth_levels(
        dataset='cmip6_historical', model=model, variable='thetao'
    )
    ref_depth = float(levels[0])

    area = load_areas(
        dataset='cmip6_historical', model=model, variable='thetao'
    )
    total_area = _val(area.sum(skipna=True))

    op = compute_mixed_layer_depth(
        temp_data=thetao,
        sal_data=so,
        ref_depth=ref_depth,
        method='density',
        density_threshold=0.03,
        temperature_type='auto'
    )
    op = op.assign_attrs({**thetao.attrs, **op.attrs, 'units': 'm', 'grid_label': gl})

    for month in months:
        op_m = op.sel(time=op.time.dt.month == month)
        pub_m = mlotst.sel(time=mlotst.time.dt.month == month)
        op_m, pub_m = xr.align(op_m, pub_m, join='inner')

        # Same-grid comparison mask: retain only cells valid in all 30 instances.
        mask = (op_m.notnull() & pub_m.notnull()).all(dim='time').compute()

        op_mm = op_m.where(mask)
        pub_mm = pub_m.where(mask)

        op_clim = (
            temporal_aggregate(op_mm, freq='climatology', op='mean')
            .sel(month=month).squeeze(drop=True).compute()
        )
        pub_clim = (
            temporal_aggregate(pub_mm, freq='climatology', op='mean')
            .sel(month=month).squeeze(drop=True).compute()
        )

        op_clim = op_clim.assign_attrs({**op.attrs, 'units': 'm'})
        pub_clim = pub_clim.assign_attrs({**mlotst.attrs, 'units': 'm'})

        metrics = compute_difference_metrics(
            data=op_clim, reference=pub_clim, dim=None
        )
        bias = _val(metrics['bias'])
        rmse = _val(metrics['rmse'])

        corr_dict = analyze_correlation(
            var1=op_clim,
            var2=pub_clim,
            dim=['lat', 'lon'],
            method='pearson',
            weighting='area'
        )
        corr = _val(corr_dict['correlation'])

        diff = (op_clim - pub_clim).assign_attrs({**op_clim.attrs, 'units': 'm'})

        exceed = compute_exceedance_fraction(
            data=diff,
            threshold=20.0,
            comparison='abs_greater',
            dim=None,
            weighting='area'
        )
        exceed_pct = _val(exceed['fraction_percent'])
        exceed_area = _val(exceed['exceeding_amount'])

        retained_area = _val(area.where(mask).sum(skipna=True))
        covered_pct = 100.0 * retained_area / total_area
        n_retained = int(_val(mask.sum()))

        missing_any = op_m.isnull().any(dim='time').compute()
        missing_all = op_m.isnull().all(dim='time').compute()
        excluded_area = _val(area.where(missing_any).sum(skipna=True))
        never_area = _val(area.where(missing_all).sum(skipna=True))

        results.append({
            'model': model,
            'grid_label': gl,
            'month': month_names[month],
            'month_num': month,
            'bias_m': bias,
            'rmse_m': rmse,
            'corr': corr,
            'exceed20_pct': exceed_pct,
            'exceed20_area_km2': exceed_area / 1e6,
            'covered_pct': covered_pct,
            'retained_cells': n_retained,
            'retained_area_km2': retained_area / 1e6,
            'total_ocean_area_km2': total_area / 1e6,
            'op_missing_any_area_km2': excluded_area / 1e6,
            'op_missing_all_area_km2': never_area / 1e6,
            'ref_depth_m': ref_depth,
        })

        store[(model, month)] = {
            'op': op_clim,
            'pub': pub_clim,
            'diff': diff,
            'mask': mask,
            'bias': bias,
            'rmse': rmse,
            'corr': corr,
            'grid_label': gl,
        }

results_df = pd.DataFrame(results)
results_df['model'] = pd.Categorical(
    results_df['model'], categories=models, ordered=True
)
results_df['month_num'] = pd.Categorical(
    results_df['month_num'], categories=months, ordered=True
)
results_df = (
    results_df.sort_values(['model', 'month_num'])
    .reset_index(drop=True)
)

ranking_df = (
    results_df.groupby('model', observed=True)
    .agg(
        mean_rmse_m=('rmse_m', 'mean'),
        mean_abs_bias_m=('bias_m', lambda x: float(np.mean(np.abs(x)))),
        mean_corr=('corr', 'mean'),
        mean_exceed20_pct=('exceed20_pct', 'mean'),
        mean_covered_pct=('covered_pct', 'mean'),
    )
    .reset_index()
)
ranking_df = ranking_df.sort_values('mean_rmse_m', ascending=True).reset_index(drop=True)
ranking_df['rank'] = np.arange(1, len(ranking_df) + 1)

# Figure 1 preparation: display-only regridding to a common 1-degree grid.
diff_regridded = {}
for model in models:
    for month in months:
        d = store[(model, month)]['diff']
        d_reg = regrid_horizontal(
            data=d,
            target='1deg',
            method='conservative',
            skipna=True
        ).compute()
        diff_regridded[(model, month)] = d_reg

pooled_vals = np.concatenate(
    [np.asarray(d.values).ravel() for d in diff_regridded.values()]
)
pooled_vals = pooled_vals[np.isfinite(pooled_vals)]
fig1_vmax = (
    float(np.nanpercentile(np.abs(pooled_vals), 98))
    if pooled_vals.size else 1.0
)
if not np.isfinite(fig1_vmax) or fig1_vmax <= 0:
    fig1_vmax = 1.0

print(results_df[['model', 'grid_label', 'month', 'bias_m', 'rmse_m', 'corr',
                  'exceed20_pct', 'covered_pct', 'op_missing_any_area_km2']]
      .to_string(index=False))
print()
print(ranking_df.to_string(index=False))
