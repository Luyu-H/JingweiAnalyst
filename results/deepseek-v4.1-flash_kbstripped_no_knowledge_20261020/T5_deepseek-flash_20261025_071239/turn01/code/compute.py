models = ['CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']
time_hist = ['1950-01', '2014-12']
time_ssp = ['2071-01', '2100-12']
depth_sss = [0, 10]


def to_psu(da):
    """Convert salinity to practical salinity units (psu)."""
    u = str(da.attrs.get('units', '')).strip().lower()
    if u in ('psu', 'pss', 'practical salinity unit'):
        out = da.copy()
        out.attrs = {**da.attrs, 'units': 'psu'}
        return out

    # Common CMIP6 / WOA encodings where values are ~0.035 and need *1000
    if u in ('1e-3', '0.001', 'g/kg', 'kg/kg') or ('1e-3' in u) or ('0.001' in u):
        out = da * 1000.0
        out.attrs = {**da.attrs, 'units': 'psu'}
        return out

    if u in ('1', ''):
        # Check a small slice to distinguish 0.035 from 35.
        isel = {d: 0 for d in da.dims if d in ('time', 'depth')}
        sample = float(da.isel(isel).max(skipna=True).compute()) if isel else float(da.max(skipna=True).compute())
        if sample < 1.0:
            out = da * 1000.0
        else:
            out = da.copy()
        out.attrs = {**da.attrs, 'units': 'psu'}
        return out

    return transform_unit(da, from_unit=da.attrs.get('units', 'psu'), to_unit='psu')


def copy_grid_attrs(out, src, units=None):
    out.attrs = {**src.attrs}
    if units is not None:
        out.attrs['units'] = units
    return out


def compute_model_sss(so, vol):
    """SSS = full-cell-volume-weighted mean of model levels within 0-10 m."""
    so_psu = to_psu(so)
    num = (so_psu * vol).sum(dim='depth', skipna=True)
    den = vol.where(so_psu.notnull()).sum(dim='depth', skipna=True)
    sss = num / den
    return copy_grid_attrs(sss, so_psu, units='psu')


# ---------------------------------------------------------------------------
# Load CMIP6 salinity, cell volumes, cell areas and level definitions
# ---------------------------------------------------------------------------
so_hist = {}
so_ssp = {}
vol_hist = {}
vol_ssp = {}
area_hist = {}
area_ssp = {}
lev_hist = {}
lev_ssp = {}
sss_hist = {}
sss_ssp = {}
variant_hist = {}
variant_ssp = {}

for model in models:
    # Historical
    so_h = load_regular_data(
        dataset='cmip6_historical', variable='so', model=model,
        time_range=time_hist, depth_range=depth_sss
    )
    vol_h = load_volumes_for_cmip6(
        dataset='cmip6_historical', model=model, variable='so',
        depth_range=depth_sss
    )
    area_h = load_areas(dataset='cmip6_historical', model=model, variable='so')
    lev_h = load_depth_levels(dataset='cmip6_historical', model=model, variable='so')

    so_hist[model] = so_h
    vol_hist[model] = vol_h
    area_hist[model] = area_h
    lev_hist[model] = lev_h
    variant_hist[model] = so_h.attrs.get(
        'variant_label',
        so_h.attrs.get('member_id', so_h.attrs.get('ensemble_member', 'unknown'))
    )
    sss_hist[model] = compute_model_sss(so_h, vol_h)

    # SSP2-4.5
    so_s = load_regular_data(
        dataset='cmip6_ssp245', variable='so', model=model,
        time_range=time_ssp, depth_range=depth_sss
    )
    vol_s = load_volumes_for_cmip6(
        dataset='cmip6_ssp245', model=model, variable='so',
        depth_range=depth_sss
    )
    area_s = load_areas(dataset='cmip6_ssp245', model=model, variable='so')
    lev_s = load_depth_levels(dataset='cmip6_ssp245', model=model, variable='so')

    so_ssp[model] = so_s
    vol_ssp[model] = vol_s
    area_ssp[model] = area_s
    lev_ssp[model] = lev_s
    variant_ssp[model] = so_s.attrs.get(
        'variant_label',
        so_s.attrs.get('member_id', so_s.attrs.get('ensemble_member', 'unknown'))
    )
    sss_ssp[model] = compute_model_sss(so_s, vol_s)


# ---------------------------------------------------------------------------
# Load WOA23 salinity, depth bounds, cell areas, land-sea and basin masks
# ---------------------------------------------------------------------------
woa_sal = load_climatology_data(
    dataset='woa23', variable='salinity', period='1971-2000',
    time_slice='annual', depth_range=depth_sss
)
woa_lev, woa_bnd = load_depth_levels(
    dataset='woa23', variable='salinity', time_slice='annual',
    period='1971-2000', depth_range=depth_sss
)
woa_area = load_areas(dataset='woa23')
landsea = load_landsea_mask()
basin_mask = load_basin_mask(depth_levels='0m')

woa_sal_psu = to_psu(woa_sal)
woa_thick = np.abs(woa_bnd[:, 1] - woa_bnd[:, 0])
woa_weights = woa_area * xr.DataArray(
    woa_thick, dims='depth', coords={'depth': woa_sal_psu.depth}
)
num_woa = (woa_sal_psu * woa_weights).sum(dim='depth', skipna=True)
den_woa = woa_weights.where(woa_sal_psu.notnull()).sum(dim='depth', skipna=True)
sss_woa = num_woa / den_woa
sss_woa = copy_grid_attrs(sss_woa, woa_sal_psu, units='psu')


# ---------------------------------------------------------------------------
# 1971-2000 mean state: model and WOA23 global ocean means
# ---------------------------------------------------------------------------
woa_global_7100 = float(
    compute_basic_stat(sss_woa, op='area_weighted_mean', dim=['lat', 'lon']).compute()
)

clim_7100 = {}
model_global_7100 = {}
model_bias_global = {}

for model in models:
    clim = temporal_aggregate(
        sss_hist[model].sel(time=slice('1971-01', '2000-12')),
        freq='mean', op='mean'
    )
    clim = copy_grid_attrs(clim, sss_hist[model], units='psu')
    clim_7100[model] = clim
    model_global_7100[model] = float(
        compute_basic_stat(clim, op='area_weighted_mean', dim=['lat', 'lon']).compute()
    )
    model_bias_global[model] = model_global_7100[model] - woa_global_7100


# ---------------------------------------------------------------------------
# Spatial bias metrics against WOA23 on the WOA23 1-degree grid
# ---------------------------------------------------------------------------
bias_maps = {}
rmse_vals = {}
corr_vals = {}

for model in models:
    model_reg = regrid_horizontal(
        clim_7100[model], target='woa23', method='conservative'
    )
    mask = build_comparison_mask([model_reg, sss_woa], ocean_only=True)
    data_m = model_reg.where(mask)
    ref_m = sss_woa.where(mask)

    metrics = compute_difference_metrics(
        data=data_m, reference=ref_m, dim=['lat', 'lon'], weighting='area'
    )
    bias_maps[model] = metrics['difference'].compute()
    # compute_difference_metrics returns scalars as plain floats, so no .compute()
    rmse_vals[model] = float(metrics['rmse'])

    corr = analyze_correlation(
        data_m, ref_m, dim=['lat', 'lon'], method='pearson', weighting='area'
    )
    corr_vals[model] = float(corr['correlation'].compute())


# ---------------------------------------------------------------------------
# Historical and ssp245 global-mean SSS trends
# ---------------------------------------------------------------------------
hist_monthly_global = {}
hist_annual_global = {}
hist_fit = {}
hist_trend = {}

for model in models:
    ser = compute_basic_stat(
        sss_hist[model], op='area_weighted_mean', dim=['lat', 'lon']
    ).compute()
    ann = temporal_aggregate(ser, freq='year', op='mean').compute()
    fit = fit_1D_trend(ann, method='linear')

    hist_monthly_global[model] = ser
    hist_annual_global[model] = ann
    hist_fit[model] = fit
    hist_trend[model] = {
        'slope_psu_per_decade': float(fit['trend']['slope'] * 10.0),
        'ci95_low': float(fit['trend']['ci95_low'] * 10.0),
        'ci95_high': float(fit['trend']['ci95_high'] * 10.0),
        'r2': float(fit['metrics']['r2']),
        'significant': bool(fit['trend']['significant']),
        'units': 'psu per decade',
    }

ssp_monthly_global = {}
ssp_annual_global = {}
ssp_fit = {}
ssp_trend = {}

for model in models:
    ser = compute_basic_stat(
        sss_ssp[model], op='area_weighted_mean', dim=['lat', 'lon']
    ).compute()
    ann = temporal_aggregate(ser, freq='year', op='mean').compute()
    fit = fit_1D_trend(ann, method='linear')

    ssp_monthly_global[model] = ser
    ssp_annual_global[model] = ann
    ssp_fit[model] = fit
    ssp_trend[model] = {
        'slope_psu_per_decade': float(fit['trend']['slope'] * 10.0),
        'ci95_low': float(fit['trend']['ci95_low'] * 10.0),
        'ci95_high': float(fit['trend']['ci95_high'] * 10.0),
        'r2': float(fit['metrics']['r2']),
        'significant': bool(fit['trend']['significant']),
        'units': 'psu per decade',
    }

slope_ratio = {
    model: ssp_trend[model]['slope_psu_per_decade'] /
           hist_trend[model]['slope_psu_per_decade']
    for model in models
}


# ---------------------------------------------------------------------------
# Cell-by-cell historical trend maps (1950-2014) and display regridding
# ---------------------------------------------------------------------------
trend_slope_reg = {}
trend_sig_reg = {}

for model in models:
    ann_map = temporal_aggregate(sss_hist[model], freq='year', op='mean')
    year_dim = 'year' if 'year' in ann_map.dims else 'time'
    tmap = fit_trend_map(
        ann_map, dim=year_dim, per_decade=True, alpha=0.05
    )

    slope = tmap['slope']
    sig = tmap['significant']

    slope_reg = regrid_horizontal(
        slope, target='woa23', method='conservative'
    ).compute()

    sig_float = sig.astype(float)
    sig_float.attrs = {**sig.attrs}
    sig_reg = regrid_horizontal(
        sig_float, target='woa23', method='nearest_s2d'
    ).compute() > 0.5

    trend_slope_reg[model] = slope_reg
    trend_sig_reg[model] = sig_reg


# ---------------------------------------------------------------------------
# Delta SSS: ssp245 2071-2100 minus historical 1995-2014
# ---------------------------------------------------------------------------
delta_global = {}
delta_basins = {model: {} for model in models}
basins = ['Atlantic', 'Indian', 'Pacific']

for model in models:
    hist_1995_2014 = temporal_aggregate(
        sss_hist[model].sel(time=slice('1995-01', '2014-12')),
        freq='mean', op='mean'
    )
    ssp_2071_2100 = temporal_aggregate(
        sss_ssp[model], freq='mean', op='mean'
    )

    hist_1995_2014 = copy_grid_attrs(hist_1995_2014, sss_hist[model], units='psu')
    ssp_2071_2100 = copy_grid_attrs(ssp_2071_2100, sss_ssp[model], units='psu')

    delta = ssp_2071_2100 - hist_1995_2014
    delta = copy_grid_attrs(delta, sss_ssp[model], units='psu')

    delta_global[model] = float(
        compute_basic_stat(
            delta, op='area_weighted_mean', dim=['lat', 'lon']
        ).compute()
    )

    ssp_reg = regrid_horizontal(
        ssp_2071_2100, target='woa23', method='conservative'
    )
    hist_reg = regrid_horizontal(
        hist_1995_2014, target='woa23', method='conservative'
    )
    delta_reg = ssp_reg - hist_reg
    delta_reg = copy_grid_attrs(delta_reg, ssp_reg, units='psu')

    for basin in basins:
        bfield = select_basin(
            delta_reg, basin=basin,
            include_marginal_seas=True, depth_dependent=False
        )
        delta_basins[model][basin] = float(
            compute_basic_stat(
                bfield, op='area_weighted_mean', dim=['lat', 'lon']
            ).compute()
        )


# ---------------------------------------------------------------------------
# Ensemble-member caveat
# ---------------------------------------------------------------------------
same_member = {
    model: (variant_hist[model] == variant_ssp[model])
    for model in models
}


# ---------------------------------------------------------------------------
# Print numerical summary
# ---------------------------------------------------------------------------
print('Global mean 1971-2000 SSS (psu)')
print(f'WOA23: {woa_global_7100:.2f}')
for model in models:
    print(f'{model}: {model_global_7100[model]:.2f} '
          f'(model minus WOA23: {model_bias_global[model]:+.2f})')

print('\nSpatial metrics against WOA23 on the WOA23 1-degree grid')
for model in models:
    print(f'{model}: RMSE = {rmse_vals[model]:.2f} psu, '
          f'r = {corr_vals[model]:.3f}')

print('\nHistorical 1950-2014 global-mean SSS trend')
for model in models:
    h = hist_trend[model]
    print(f"{model}: {h['slope_psu_per_decade']:.3g} psu/decade, "
          f"95% CI ({h['ci95_low']:.3g}, {h['ci95_high']:.3g}), "
          f"R2 = {h['r2']:.3g}, significant = {h['significant']}")

print('\nSSP2-4.5 2071-2100 global-mean SSS trend')
for model in models:
    s = ssp_trend[model]
    print(f"{model}: {s['slope_psu_per_decade']:.3g} psu/decade, "
          f"95% CI ({s['ci95_low']:.3g}, {s['ci95_high']:.3g}), "
          f"R2 = {s['r2']:.3g}, significant = {s['significant']}")

print('\nRatio of ssp245 slope to historical slope (descriptive only)')
for model in models:
    print(f'{model}: {slope_ratio[model]:.1f}')

print('\nDelta SSS, ssp245 2071-2100 minus historical 1995-2014 (psu)')
for model in models:
    print(f"{model}: global = {delta_global[model]:+.2f}, "
          f"Atlantic = {delta_basins[model]['Atlantic']:+.2f}, "
          f"Indian = {delta_basins[model]['Indian']:+.2f}, "
          f"Pacific = {delta_basins[model]['Pacific']:+.2f}")

print('\nEnsemble members used')
for model in models:
    print(f'{model}: historical = {variant_hist[model]}, '
          f'ssp245 = {variant_ssp[model]}, same member = {same_member[model]}')

print('\nCaveat: the ssp245/historical slope ratio is descriptive only. '
      'The two periods have different mean states, forcing and separate uncertainties, '
      'and the ratio is especially weak if either slope is near zero. '
      'If the variant labels differ between a model\'s two experiments, '
      'delta and slope-ratio differences also include member-to-member differences.')
