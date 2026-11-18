import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

models = ['CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']


def _scalar(x):
    """Materialise a possibly-lazy scalar and return a Python float."""
    if hasattr(x, 'compute'):
        x = x.compute()
    return float(x)


# ---------------------------------------------------------------------------
# WOA23 reference SSS: mean salinity of levels within 0-10 m, weighted by
# full cell volume (area x full layer thickness), then area-weighted globally.
# ---------------------------------------------------------------------------
woa_sal = load_climatology_data(
    dataset='woa23',
    variable='salinity',
    time_slice='annual',
    period='1971-2000',
    depth_range=[0, 10],
)
woa_sss = compute_basic_stat(data=woa_sal, op='volume_weighted_mean', dim='depth')
woa_sss = woa_sss.assign_attrs({**woa_sal.attrs, **woa_sss.attrs}).compute()

woa_global_mean_val = _scalar(
    compute_basic_stat(data=woa_sss, op='area_weighted_mean', dim=['lat', 'lon'])
)

# Masks and vertical-grid information are loaded for provenance / basin work.
landsea = load_landsea_mask()
basinmask = load_basin_mask()
woa_areas = load_areas(dataset='woa23')
woa_levels, woa_bounds = load_depth_levels(
    dataset='woa23',
    variable='salinity',
    time_slice='annual',
    period='1971-2000',
    depth_range=[0, 10],
)

# Storage for results
bias_fields = {}
trend_fields = {}
sig_fields = {}
hist_series = {}
ssp_series = {}
hist_fits = {}
ssp_fits = {}

model_global_mean = {}
model_global_diff = {}
model_bias = {}
model_rmse = {}
model_corr = {}
model_hist_trend = {}
model_ssp_trend = {}
model_ratio = {}
model_ratio_meaningful = {}
model_delta_global = {}
model_delta_basins = {}

for model in models:
    # -----------------------------------------------------------------------
    # Historical experiment: monthly SSS at every native model cell
    # -----------------------------------------------------------------------
    hist_so = load_regular_data(
        dataset='cmip6_historical',
        variable='so',
        model=model,
        time_range=['1950-01', '2014-12'],
        depth_range=[0, 10],
    )
    hist_sss_monthly = compute_basic_stat(
        data=hist_so, op='volume_weighted_mean', dim='depth'
    )
    hist_sss_monthly = hist_sss_monthly.assign_attrs(
        {**hist_so.attrs, **hist_sss_monthly.attrs}
    )

    # 1971-2000 day-weighted mean SSS at each native cell
    hist_period = hist_sss_monthly.sel(time=slice('1971-01', '2000-12'))
    hist_clim = temporal_aggregate(data=hist_period, freq='mean', op='mean')
    hist_clim = hist_clim.assign_attrs(
        {**hist_sss_monthly.attrs, **hist_clim.attrs}
    ).compute()

    # Native-grid global mean, using the model's own areacello through provenance
    gm = _scalar(
        compute_basic_stat(data=hist_clim, op='area_weighted_mean', dim=['lat', 'lon'])
    )
    model_global_mean[model] = gm
    model_global_diff[model] = gm - woa_global_mean_val

    # -----------------------------------------------------------------------
    # Comparison against WOA23 on the 1-degree WOA grid
    # -----------------------------------------------------------------------
    model_regrid = regrid_horizontal(
        data=hist_clim, target='woa23', method='conservative'
    )
    mask = build_comparison_mask(
        fields=[model_regrid, woa_sss], ocean_only=True
    )
    model_masked = model_regrid.where(mask)
    woa_masked = woa_sss.where(mask)

    metrics = compute_difference_metrics(
        data=model_masked,
        reference=woa_masked,
        dim=['lat', 'lon'],
        weighting='area',
    )
    bias_fields[model] = metrics['difference'].compute()
    model_bias[model] = _scalar(metrics['bias'])
    model_rmse[model] = _scalar(metrics['rmse'])

    corr = analyze_correlation(
        var1=model_masked,
        var2=woa_masked,
        dim=['lat', 'lon'],
        method='pearson',
        weighting='area',
    )
    model_corr[model] = _scalar(corr['correlation'])

    # -----------------------------------------------------------------------
    # Historical annual global-mean SSS and its trend
    # -----------------------------------------------------------------------
    hist_annual = temporal_aggregate(data=hist_sss_monthly, freq='year', op='mean')
    hist_annual = hist_annual.assign_attrs(
        {**hist_sss_monthly.attrs, **hist_annual.attrs}
    )
    hist_gm_annual = compute_basic_stat(
        data=hist_annual, op='area_weighted_mean', dim=['lat', 'lon']
    ).compute()
    hist_series[model] = hist_gm_annual

    hist_fit = fit_1D_trend(data=hist_gm_annual, method='linear')
    hist_fits[model] = hist_fit

    # Cell-by-cell trend map on the native grid; regrid only for display.
    trend_dim = 'year' if 'year' in hist_annual.dims else 'time'
    trend_map = fit_trend_map(
        data=hist_annual, dim=trend_dim, per_decade=True
    )
    trend_fields[model] = regrid_horizontal(
        data=trend_map['slope'], target='woa23', method='conservative'
    ).compute()

    sig_float = trend_map['significant'].astype(float)
    sig_regrid = regrid_horizontal(
        data=sig_float, target='woa23', method='nearest_s2d'
    )
    sig_fields[model] = (sig_regrid > 0.5).compute()

    # -----------------------------------------------------------------------
    # ssp245 experiment: monthly SSS, annual global mean, trend
    # -----------------------------------------------------------------------
    ssp_so = load_regular_data(
        dataset='cmip6_ssp245',
        variable='so',
        model=model,
        time_range=['2071-01', '2100-12'],
        depth_range=[0, 10],
    )
    ssp_sss_monthly = compute_basic_stat(
        data=ssp_so, op='volume_weighted_mean', dim='depth'
    )
    ssp_sss_monthly = ssp_sss_monthly.assign_attrs(
        {**ssp_so.attrs, **ssp_sss_monthly.attrs}
    )

    ssp_annual = temporal_aggregate(data=ssp_sss_monthly, freq='year', op='mean')
    ssp_annual = ssp_annual.assign_attrs(
        {**ssp_sss_monthly.attrs, **ssp_annual.attrs}
    )
    ssp_gm_annual = compute_basic_stat(
        data=ssp_annual, op='area_weighted_mean', dim=['lat', 'lon']
    ).compute()
    ssp_series[model] = ssp_gm_annual

    ssp_fit = fit_1D_trend(data=ssp_gm_annual, method='linear')
    ssp_fits[model] = ssp_fit

    # Trend diagnostics. fit_1D_trend returns a slope per year on the year axis.
    h = hist_fit['trend']
    s = ssp_fit['trend']
    model_hist_trend[model] = {
        'slope_per_decade': h['slope'] * 10.0,
        'ci95_low_per_decade': h['ci95_low'] * 10.0,
        'ci95_high_per_decade': h['ci95_high'] * 10.0,
        'r2': hist_fit['metrics']['r2'],
        'p_value': h['p_value'],
        'significant': bool(h['significant']),
    }
    model_ssp_trend[model] = {
        'slope_per_decade': s['slope'] * 10.0,
        'ci95_low_per_decade': s['ci95_low'] * 10.0,
        'ci95_high_per_decade': s['ci95_high'] * 10.0,
        'r2': ssp_fit['metrics']['r2'],
        'p_value': s['p_value'],
        'significant': bool(s['significant']),
    }
    model_ratio[model] = (
        float(s['slope'] / h['slope']) if h['slope'] != 0 else np.nan
    )
    # The ratio is only useful if the historical slope is distinguishable from 0.
    model_ratio_meaningful[model] = bool(h['significant'])

    # -----------------------------------------------------------------------
    # Delta SSS: ssp245 2071-2100 mean minus historical 1995-2014 mean
    # -----------------------------------------------------------------------
    hist_delta_period = hist_sss_monthly.sel(time=slice('1995-01', '2014-12'))
    hist_delta_mean = temporal_aggregate(
        data=hist_delta_period, freq='mean', op='mean'
    )
    hist_delta_mean = hist_delta_mean.assign_attrs(
        {**hist_sss_monthly.attrs, **hist_delta_mean.attrs}
    )

    ssp_mean = temporal_aggregate(
        data=ssp_sss_monthly, freq='mean', op='mean'
    )
    ssp_mean = ssp_mean.assign_attrs(
        {**ssp_sss_monthly.attrs, **ssp_mean.attrs}
    )

    delta_native = ssp_mean - hist_delta_mean
    delta_native = delta_native.assign_attrs(
        {**ssp_mean.attrs, **delta_native.attrs, 'units': 'psu'}
    )

    delta_regrid = regrid_horizontal(
        data=delta_native, target='woa23', method='conservative'
    )
    delta_mask = build_comparison_mask(
        fields=[delta_regrid, woa_sss], ocean_only=True
    )
    delta_masked = delta_regrid.where(delta_mask)
    delta_masked = delta_masked.assign_attrs(
        {**delta_regrid.attrs, **delta_masked.attrs}
    )

    model_delta_global[model] = _scalar(
        compute_basic_stat(
            data=delta_masked, op='area_weighted_mean', dim=['lat', 'lon']
        )
    )

    model_delta_basins[model] = {}
    for basin_name in ['Atlantic', 'Indian', 'Pacific']:
        db = select_basin(
            data=delta_masked,
            basin=basin_name,
            include_marginal_seas=True,
        )
        model_delta_basins[model][basin_name] = _scalar(
            compute_basic_stat(
                data=db, op='area_weighted_mean', dim=['lat', 'lon']
            )
        )

ensemble_caveat = (
    "Delta SSS combines the ssp245 and historical experiments; the two may "
    "use different ensemble members or variant labels for a given model."
)

# Print a compact numerical summary
print('WOA23 global mean SSS (psu):', f'{woa_global_mean_val:.2f}')
for m in models:
    ht = model_hist_trend[m]
    st = model_ssp_trend[m]
    print(f'\n{m}')
    print(
        f'  1971-2000 global mean SSS = {model_global_mean[m]:.2f} psu; '
        f'difference vs WOA23 = {model_global_diff[m]:+.2f} psu'
    )
    print(
        f'  bias = {model_bias[m]:+.2f} psu; RMSE = {model_rmse[m]:.2f} psu; '
        f'spatial r = {model_corr[m]:.3f}'
    )
    print(
        f'  historical 1950-2014 trend = {ht["slope_per_decade"]:.3g} psu/decade '
        f'(95% CI {ht["ci95_low_per_decade"]:.3g} to {ht["ci95_high_per_decade"]:.3g}), '
        f'R2 = {ht["r2"]:.3f}, significant = {ht["significant"]}'
    )
    print(
        f'  ssp245 2071-2100 trend = {st["slope_per_decade"]:.3g} psu/decade '
        f'(95% CI {st["ci95_low_per_decade"]:.3g} to {st["ci95_high_per_decade"]:.3g}), '
        f'R2 = {st["r2"]:.3f}, significant = {st["significant"]}'
    )
    print(
        f'  ssp245/historical slope ratio = {model_ratio[m]:.1f} '
        f'(meaningful = {model_ratio_meaningful[m]})'
    )
    print(
        f'  delta SSS global = {model_delta_global[m]:+.2f} psu; '
        f'Atlantic = {model_delta_basins[m]["Atlantic"]:+.2f}, '
        f'Indian = {model_delta_basins[m]["Indian"]:+.2f}, '
        f'Pacific = {model_delta_basins[m]["Pacific"]:+.2f}'
    )
print('\nCaveat:', ensemble_caveat)
