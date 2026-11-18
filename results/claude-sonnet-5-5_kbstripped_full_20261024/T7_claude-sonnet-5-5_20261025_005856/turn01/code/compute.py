models = ['CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']
sources = ['WOA23'] + models
trange = ['1971-01', '2000-12']
drange = [200, 1000]
RHO = 1025.0                       # kg m-3, used only for the regression/scatter conversion
hdims = ['depth', 'lat', 'lon']    # curvilinear index dims are mapped by compute_basic_stat
basins = ['Pacific', 'Atlantic', 'Indian']

# ---------------------------------------------------------------- 1. load
native = {}        # per source: time-mean no3, po4 on the native grid and levels
weights = {}       # per source: volume (models) or layer thickness (WOA23)
vol_source = {}
layer_info = {}

for m in models:
    no3 = load_regular_data(dataset='cmip6_historical', variable='no3', model=m,
                            time_range=trange, depth_range=drange)
    po4 = load_regular_data(dataset='cmip6_historical', variable='po4', model=m,
                            time_range=trange, depth_range=drange)
    # day-weighted 1971-2000 climatology on the native grid
    no3c = temporal_aggregate(no3, freq='mean').compute()
    po4c = temporal_aggregate(po4, freq='mean').compute()
    # one common valid sample for NO3 and PO4 (PO4 > 0 so the ratio is defined)
    po4c = po4c.where(po4c > 0)
    po4c = po4c.where(no3c.notnull())
    no3c = no3c.where(po4c.notnull())
    native[m] = dict(no3=no3c, po4=po4c)

    vol = load_volumes_for_cmip6(dataset='cmip6_historical', model=m,
                                 depth_range=drange, variable='no3')
    vol_source[m] = vol.attrs.get('volume_source')
    if 'time' in vol.dims:
        vol = vol.mean('time')
    weights[m] = vol.compute()

    lv, bd = load_depth_levels(dataset='cmip6_historical', model=m, variable='no3',
                               depth_range=drange)
    layer_info[m] = dict(n_levels=len(lv), top=float(np.min(bd)), bottom=float(np.max(bd)))

no3w = load_climatology_data(dataset='woa23', variable='nitrate', depth_range=drange).compute()
po4w = load_climatology_data(dataset='woa23', variable='phosphate', depth_range=drange).compute()
po4w = po4w.where(po4w > 0)
po4w = po4w.where(no3w.notnull())
no3w = no3w.where(po4w.notnull())
native['WOA23'] = dict(no3=no3w, po4=po4w)

lv, bd = load_depth_levels(dataset='woa23', variable='nitrate', depth_range=drange)
bd = np.asarray(bd)
assert len(lv) == no3w.sizes['depth']
dz = xr.DataArray(bd[:, 1] - bd[:, 0], dims='depth', coords={'depth': no3w['depth'].values})
weights['WOA23'] = dz
vol_source['WOA23'] = 'layer thickness from WOA23 depth-level bounds x 1-degree cell areas'
layer_info['WOA23'] = dict(n_levels=len(lv), top=float(bd.min()), bottom=float(bd.max()))


def column_mean(field, w):
    """Thickness/volume-weighted mean over the depth dimension only, NaN-aware."""
    if w.ndim > 1:
        w, field = xr.align(w, field, join='override')
    ww = w.where(field.notnull())
    s = ww.sum('depth')
    out = (field * ww).sum('depth') / s
    return out.where(s > 0)


# ---------------------------------------------------------------- 2-5. native-grid statistics
stats = {}
colmaps = {}      # native column-mean fields (for later regridding)
ratio3d = {}
for s in sources:
    no3n, po4n = native[s]['no3'], native[s]['po4']
    ratio = compute_ratio(no3n, po4n).assign_attrs(units='1')
    ratio3d[s] = ratio

    mean_no3 = float(compute_basic_stat(no3n, 'volume_weighted_mean', dim=hdims).compute())
    mean_po4 = float(compute_basic_stat(po4n, 'volume_weighted_mean', dim=hdims).compute())
    headline = float(compute_basic_stat(ratio, 'volume_weighted_mean', dim=hdims).compute())

    ex1418 = compute_exceedance_fraction(ratio, threshold=2.0, comparison='abs_less', reference=16.0)
    exlow = compute_exceedance_fraction(ratio, threshold=12.0, comparison='less')

    d = dict(
        mean_no3_native_units=mean_no3,
        mean_po4_native_units=mean_po4,
        np_headline_cell_by_cell=headline,
        np_ratio_of_means=mean_no3 / mean_po4,
        pct_14_18=float(ex1418['fraction_percent']),
        vol_14_18_m3=float(ex1418['exceeding_amount']),
        total_volume_m3=float(exlow['total_amount']),
        pct_below12_global=float(exlow['fraction_percent']),
        vol_below12_global_m3=float(exlow['exceeding_amount']),
        basin_mean={}, basin_pct_below12={}, basin_vol_below12_m3={},
    )
    for b in basins:
        rb = select_basin(ratio, b, include_marginal_seas=True, depth_dependent=False)
        d['basin_mean'][b] = float(compute_basic_stat(rb, 'volume_weighted_mean', dim=hdims).compute())
        if b in ('Pacific', 'Indian'):
            exb = compute_exceedance_fraction(rb, threshold=12.0, comparison='less')
            d['basin_pct_below12'][b] = float(exb['fraction_percent'])
            d['basin_vol_below12_m3'][b] = float(exb['exceeding_amount'])
    stats[s] = d

    colmaps[s] = dict(
        ratio=column_mean(ratio, weights[s]).assign_attrs(units='1').compute(),
    )

# ---------------------------------------------------------------- 6. regrid and common mask
ref_lat, ref_lon = native['WOA23']['no3']['lat'], native['WOA23']['no3']['lon']


def snap(x):
    return x.assign_coords(lat=ref_lat.values, lon=ref_lon.values)


reg = {}
for s in sources:
    if s == 'WOA23':
        reg[s] = dict(no3=native[s]['no3'], po4=native[s]['po4'], colratio=colmaps[s]['ratio'])
    else:
        reg[s] = dict(
            no3=snap(regrid_horizontal(native[s]['no3'], target='woa23', method='conservative')).compute(),
            po4=snap(regrid_horizontal(native[s]['po4'], target='woa23', method='conservative')).compute(),
            colratio=snap(regrid_horizontal(colmaps[s]['ratio'], target='woa23', method='conservative')).compute(),
        )
    # keep only 3-D cells where both nutrients are present
    ok = reg[s]['no3'].notnull() & reg[s]['po4'].notnull()
    reg[s]['no3'] = reg[s]['no3'].where(ok)
    reg[s]['po4'] = reg[s]['po4'].where(ok)

valid_maps = []
for s in sources:
    v = (reg[s]['no3'].notnull() & reg[s]['po4'].notnull()).any('depth')
    v = v.astype(float).where(v)
    v.attrs.update(reg[s]['no3'].attrs)
    valid_maps.append(v)

mask = build_comparison_mask(valid_maps, ocean_only=True)
mask_info = dict(attrs=dict(mask.attrs), n_cells_shared=int(mask.sum()),
                 n_cells_each={s: int(v.notnull().sum()) for s, v in zip(sources, valid_maps)})

# ---------------------------------------------------------------- 7. unit conversion and regression
fits = {}
scatter = {}
rng = np.random.default_rng(0)
for s in sources:
    no3r, po4r = reg[s]['no3'], reg[s]['po4']
    if s != 'WOA23':
        no3r = transform_unit(no3r, from_unit='mol m-3', to_unit='umol/kg', density=RHO / 1000.0)
        po4r = transform_unit(po4r, from_unit='mol m-3', to_unit='umol/kg', density=RHO / 1000.0)
    no3m = no3r.where(mask)
    po4m = po4r.where(mask)
    pair_ok = no3m.notnull() & po4m.notnull()
    no3m = no3m.where(pair_ok)
    po4m = po4m.where(pair_ok)

    fit = fit_multivariate_relation(target=no3m, predictors=[po4m], method='linear')

    x = np.asarray(po4m.values).ravel()
    y = np.asarray(no3m.values).ravel()
    good = np.isfinite(x) & np.isfinite(y)
    x, y = x[good], y[good]
    np_slope, np_int = np.polyfit(x, y, 1)
    np_r2 = float(np.corrcoef(x, y)[0, 1] ** 2)

    slope = intercept = None
    co = fit.get('coefficients')
    if isinstance(co, dict):
        for k, v in co.items():
            val = float(np.asarray(v).squeeze())
            if 'intercept' in str(k).lower():
                intercept = val
            elif slope is None:
                slope = val
    coef_source = 'fit_multivariate_relation'
    if slope is None or intercept is None:
        slope, intercept, coef_source = float(np_slope), float(np_int), 'numpy polyfit (tool coefficients not parsable)'
    r2 = float(fit['metrics']['r2'])

    fits[s] = dict(a=slope, b_umol_kg=intercept, r2=r2, n_cells=int(len(x)),
                   np_check=dict(a=float(np_slope), b=float(np_int), r2=np_r2),
                   coef_source=coef_source, vif=fit['diagnostics'].get('vif') if isinstance(fit.get('diagnostics'), dict) else None)
    idx = rng.choice(len(x), size=min(8000, len(x)), replace=False)
    scatter[s] = dict(x=x[idx], y=y[idx], x_p995=float(np.percentile(x, 99.5)),
                      y_p995=float(np.percentile(y, 99.5)))

# ---------------------------------------------------------------- 8. Redfield and comparison with WOA23
col_masked = {s: reg[s]['colratio'].where(mask) for s in sources}
diff_metrics = {}
for m in models:
    dm = compute_difference_metrics(col_masked[m], col_masked['WOA23'])
    diff_metrics[m] = dict(bias=float(dm['bias']), rmse=float(dm['rmse']),
                           centered_rmse=float(dm['centered_rmse']), n_valid=int(dm['n_valid']))

tol = 0.5
redfield = {}
for s in sources:
    h = stats[s]['np_headline_cell_by_cell']
    redfield[s] = 'above 16' if h > 16 + tol else ('below 16' if h < 16 - tol else 'at 16 (within +/-0.5)')

low_np = {}
w_glob = stats['WOA23']['pct_below12_global']
for m in models:
    rg = stats[m]['pct_below12_global'] / w_glob if w_glob > 0 else np.nan
    rp = stats[m]['basin_pct_below12']['Pacific'] / max(stats['WOA23']['basin_pct_below12']['Pacific'], 1e-9)
    ri = stats[m]['basin_pct_below12']['Indian'] / max(stats['WOA23']['basin_pct_below12']['Indian'], 1e-9)
    low_np[m] = dict(global_ratio_to_woa=rg, pacific_ratio_to_woa=rp, indian_ratio_to_woa=ri,
                     reproduces_volume=bool(0.5 <= rg <= 2.0),
                     reproduces_basins=bool(0.5 <= rp <= 2.0 and 0.5 <= ri <= 2.0))

for s in sources:
    d = stats[s]
    print(f"{s}: headline N:P={d['np_headline_cell_by_cell']:.2f}, ratio of means={d['np_ratio_of_means']:.2f}, "
          f"14<N:P<18: {d['pct_14_18']:.1f}%, N:P<12: {d['pct_below12_global']:.1f}% "
          f"({d['vol_below12_global_m3']:.3e} m3), basins {({k: round(v, 2) for k, v in d['basin_mean'].items()})}, "
          f"<12 Pac/Ind: {d['basin_pct_below12']}, Redfield: {redfield[s]}")
    f = fits[s]
    print(f"   fit a={f['a']:.2f}, b={f['b_umol_kg']:.2f} umol/kg, R2={f['r2']:.3f}, n={f['n_cells']}")
print('volume sources:', vol_source)
print('layers:', layer_info)
print('mask:', mask_info['n_cells_shared'], mask_info['n_cells_each'])
print('diff vs WOA23:', diff_metrics)
print('low N:P:', low_np)

map_data = {s: dict(field=np.asarray(col_masked[s].values) - 16.0,
                    lon=np.asarray(col_masked[s]['lon'].values),
                    lat=np.asarray(col_masked[s]['lat'].values)) for s in sources}
