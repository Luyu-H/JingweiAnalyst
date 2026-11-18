import numpy as np
import pandas as pd
import xarray as xr

TIME_WIN   = ['1971-01', '2000-12']
DEPTH_WIN  = [200.0, 1000.0]
MODELS     = ['CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']
SOURCES    = ['WOA23'] + MODELS


def _num(x):
    """Materialise a scalar out of a (possibly lazy) DataArray / numpy value."""
    if isinstance(x, xr.DataArray):
        x = x.compute()
        v = np.asarray(x.values, dtype=float).ravel()
        return float(v[0]) if v.size == 1 else v
    return float(x)


def sub_layer(da, dmin=DEPTH_WIN[0], dmax=DEPTH_WIN[1]):
    """Keep levels whose centre lies strictly inside the 200-1000 m layer."""
    d = np.asarray(da['depth'].values, dtype=float)
    idx = np.nonzero((d > dmin) & (d < dmax))[0]
    if idx.size == 0:
        raise ValueError(f"no depth levels strictly inside {dmin}-{dmax} m")
    return da.isel(depth=idx)


# --------------------------------------------------------------------------
# 1. CMIP6 nitrate / phosphate, 1971-2000, converted to umol kg-1
# --------------------------------------------------------------------------
def prep_cmip(model):
    raw_no3 = load_regular_data(dataset='cmip6_historical', variable='no3', model=model,
                                time_range=TIME_WIN, depth_range=DEPTH_WIN)
    raw_po4 = load_regular_data(dataset='cmip6_historical', variable='po4', model=model,
                                time_range=TIME_WIN, depth_range=DEPTH_WIN)

    no3 = transform_unit(raw_no3, from_unit='mol m-3', to_unit='umol kg-1', density=1.025)
    po4 = transform_unit(raw_po4, from_unit='mol m-3', to_unit='umol kg-1', density=1.025)

    prov_no3 = dict(no3.attrs)
    prov_po4 = dict(po4.attrs)

    # day-weighted 1971-2000 annual mean (NOT mean(dim='time'))
    no3a = temporal_aggregate(no3, freq='mean').compute().assign_attrs(**prov_no3)
    po4a = temporal_aggregate(po4, freq='mean').compute().assign_attrs(**prov_po4)

    no3a = sub_layer(no3a)
    po4a = sub_layer(po4a)
    no3a = no3a.assign_attrs(units='umol kg-1')
    po4a = po4a.assign_attrs(units='umol kg-1')
    return no3a, po4a


fields = {}
for m in MODELS:
    fields[m] = prep_cmip(m)

# --------------------------------------------------------------------------
# 2. WOA23 observed climatology (already umol kg-1)
# --------------------------------------------------------------------------
woa_no3 = load_climatology_data(dataset='woa23', variable='nitrate',
                                time_slice='annual', depth_range=DEPTH_WIN)
woa_po4 = load_climatology_data(dataset='woa23', variable='phosphate',
                                time_slice='annual', depth_range=DEPTH_WIN)
woa_no3 = sub_layer(woa_no3).assign_attrs(units='umol kg-1')
woa_po4 = sub_layer(woa_po4).assign_attrs(units='umol kg-1')
fields['WOA23'] = (woa_no3, woa_po4)

# provenance of the volume weights actually used downstream
vol_source = {}
for m in MODELS:
    v = load_volumes_for_cmip6(dataset='cmip6_historical', model=m, depth_range=DEPTH_WIN)
    vol_source[m] = v.attrs.get('volume_source', 'unknown')
vol_source['WOA23'] = 'area x layer thickness (WOA23 1 deg grid)'


# --------------------------------------------------------------------------
# 3. cell-by-cell N:P field for every source
# --------------------------------------------------------------------------
def np_ratio_field(no3, po4):
    r = compute_ratio(no3, po4)
    attrs = dict(no3.attrs)
    attrs['units'] = '1'
    attrs['long_name'] = 'Nitrate to phosphate ratio'
    return r.assign_attrs(**attrs)


np_fields = {s: np_ratio_field(*fields[s]) for s in SOURCES}


# --------------------------------------------------------------------------
# 4. OLS  NO3 = a*PO4 + b  (equal weight per cell) via fit_multivariate_relation
# --------------------------------------------------------------------------
def _flatten(da):
    return da.stack(cell=list(da.dims))


def _extract_line(res, xv, yv):
    """Return (slope, intercept) from the tool's coefficients, order-robust."""
    co = res.get('coefficients')
    arr = None
    if co is not None:
        if isinstance(co, dict):
            vals = [co[k] for k in co if k != 'intercept'] + ([co['intercept']] if 'intercept' in co else [])
            arr = np.asarray(vals, dtype=float).ravel()
        else:
            arr = np.asarray(co, dtype=float).ravel()
    sl_ref, ic_ref = np.polyfit(xv, yv, 1)
    if arr is not None and arr.size == 2:
        if abs(arr[0] - sl_ref) <= abs(arr[1] - sl_ref):
            return float(arr[0]), float(arr[1])
        return float(arr[1]), float(arr[0])
    return float(sl_ref), float(ic_ref)


fits = {}
for s in SOURCES:
    no3, po4 = fields[s]
    xda = _flatten(po4)
    yda = _flatten(no3)
    xv = np.asarray(xda.values, dtype=float)
    yv = np.asarray(yda.values, dtype=float)
    good = np.isfinite(xv) & np.isfinite(yv)
    xv, yv = xv[good], yv[good]

    xs = xr.DataArray(xv, dims='cell', coords={'cell': np.arange(xv.size)})
    ys = xr.DataArray(yv, dims='cell', coords={'cell': np.arange(yv.size)})
    res = fit_multivariate_relation(target=ys, predictors=[xs], method='linear',
                                    validation='random', test_size=0.2)

    a, b = _extract_line(res, xv, yv)
    r2 = float(res['metrics']['r2'])
    fits[s] = dict(a=a, b=b, r2=r2, n=int(xv.size),
                   r2_numpy=float(np.corrcoef(xv, yv)[0, 1] ** 2))


# --------------------------------------------------------------------------
# 5. headline stoichiometry, basin split, volume shares
# --------------------------------------------------------------------------
BASINS = ['Pacific', 'Atlantic', 'Indian']
results = {}

for s in SOURCES:
    no3, po4 = fields[s]
    r = np_fields[s]
    dims = list(r.dims)

    # (i) mean of the cell-by-cell ratio  -- headline
    mean_ratio = _num(compute_basic_stat(r, op='volume_weighted_mean', dim=dims))

    # (ii) ratio of the volume-weighted means
    m_no3 = _num(compute_basic_stat(no3, op='volume_weighted_mean', dim=list(no3.dims)))
    m_po4 = _num(compute_basic_stat(po4, op='volume_weighted_mean', dim=list(po4.dims)))
    ratio_of_means = m_no3 / m_po4

    # share of water volume with 14 < N:P < 18
    band = compute_exceedance_fraction(r, threshold=2.0, comparison='abs_less',
                                       reference=16.0, weighting='volume')
    pct_14_18 = _num(band['fraction_percent'])
    band_m3 = _num(band['exceeding_amount'])
    total_m3 = _num(band['total_amount'])

    # basin means of the headline (Southern / Arctic kept out of these three)
    basin_means = {}
    for b in BASINS:
        rb = select_basin(r, basin=b, include_marginal_seas=True, depth_dependent=False)
        basin_means[b] = _num(compute_basic_stat(rb, op='volume_weighted_mean', dim=list(rb.dims)))

    # volume fraction with N:P < 12
    low_g = compute_exceedance_fraction(r, threshold=12.0, comparison='less',
                                        reference=0.0, weighting='volume')
    low_p = compute_exceedance_fraction(select_basin(r, basin='Pacific', include_marginal_seas=True,
                                                     depth_dependent=False),
                                        threshold=12.0, comparison='less',
                                        reference=0.0, weighting='volume')
    low_i = compute_exceedance_fraction(select_basin(r, basin='Indian', include_marginal_seas=True,
                                                     depth_dependent=False),
                                        threshold=12.0, comparison='less',
                                        reference=0.0, weighting='volume')

    results[s] = dict(
        mean_ratio=mean_ratio,
        ratio_of_means=ratio_of_means,
        mean_no3=m_no3, mean_po4=m_po4,
        pct_14_18=pct_14_18, band_m3=band_m3, total_m3=total_m3,
        basins=basin_means,
        pct_lt12_global=_num(low_g['fraction_percent']),
        m3_lt12_global=_num(low_g['exceeding_amount']),
        vol_global_m3=_num(low_g['total_amount']),
        pct_lt12_pacific=_num(low_p['fraction_percent']),
        m3_lt12_pacific=_num(low_p['exceeding_amount']),
        vol_pacific_m3=_num(low_p['total_amount']),
        pct_lt12_indian=_num(low_i['fraction_percent']),
        m3_lt12_indian=_num(low_i['exceeding_amount']),
        vol_indian_m3=_num(low_i['total_amount']),
        a=fits[s]['a'], b=fits[s]['b'], r2=fits[s]['r2'], n_cells=fits[s]['n'],
    )


# --------------------------------------------------------------------------
# 6. 200-1000 m column N:P map, put on the WOA23 1 deg grid for comparison
# --------------------------------------------------------------------------
maps_native = {}
maps_reg = {}
for s in SOURCES:
    r = np_fields[s]
    maps_native[s] = compute_basic_stat(r, op='volume_weighted_mean', dim=['depth']).compute()
    if s == 'WOA23':
        maps_reg[s] = maps_native[s].compute()
    else:
        maps_reg[s] = regrid_horizontal(maps_native[s], target='woa23',
                                        method='conservative').compute()

cmp_mask = build_comparison_mask([maps_reg[s] for s in SOURCES], ocean_only=True)
n_mask_cells = int(np.asarray(cmp_mask.values).sum())
n_mask_total = int(np.asarray(cmp_mask.values).size)

diff_metrics = {}
for m in MODELS:
    diff_metrics[m] = compute_difference_metrics(maps_reg[m].where(cmp_mask),
                                                 maps_reg['WOA23'].where(cmp_mask))

for m in MODELS:
    dm = diff_metrics[m]
    results[m]['bias_vs_woa'] = _num(dm['bias'])
    results[m]['rmse_vs_woa'] = _num(dm['rmse'])
    results[m]['crmse_vs_woa'] = _num(dm['centered_rmse'])
    results[m]['max_abs_err_vs_woa'] = _num(dm['max_abs_error'])
    results[m]['max_abs_err_at'] = dm['max_abs_error_at']
    results[m]['mean_data_masked'] = _num(dm['mean_data'])
    results[m]['mean_woa_masked'] = _num(dm['mean_reference'])
    results[m]['n_valid_masked'] = int(_num(dm['n_valid']))


# --------------------------------------------------------------------------
# 7. summary table
# --------------------------------------------------------------------------
rows = []
for s in SOURCES:
    R = results[s]
    rows.append({
        'source': s,
        'mean_NP(mean of ratios)': round(R['mean_ratio'], 2),
        'ratio_of_means': round(R['ratio_of_means'], 2),
        'a (slope)': round(R['a'], 2),
        'b (intercept)': round(R['b'], 2),
        'R2': round(R['r2'], 3),
        'n_cells': R['n_cells'],
        '%vol 14<NP<18': round(R['pct_14_18'], 1),
        'Pacific': round(R['basins']['Pacific'], 2),
        'Atlantic': round(R['basins']['Atlantic'], 2),
        'Indian': round(R['basins']['Indian'], 2),
        '%vol NP<12 global': round(R['pct_lt12_global'], 1),
        '%vol NP<12 Pacific': round(R['pct_lt12_pacific'], 1),
        '%vol NP<12 Indian': round(R['pct_lt12_indian'], 1),
    })
summary = pd.DataFrame(rows).set_index('source')
print("200-1000 m nitrate:phosphate stoichiometry")
print(summary.to_string())
print()
for s in SOURCES:
    R = results[s]
    print(f"{s}: volume of layer = {R['total_m3']:.3e} m3 ; "
          f"14<NP<18 volume = {R['band_m3']:.3e} m3 ; "
          f"NP<12 global = {R['m3_lt12_global']:.3e} m3 ; "
          f"Pacific = {R['m3_lt12_pacific']:.3e} m3 of {R['vol_pacific_m3']:.3e} m3 ; "
          f"Indian = {R['m3_lt12_indian']:.3e} m3 of {R['vol_indian_m3']:.3e} m3")
    print(f"    mean NO3 = {R['mean_no3']:.3f}, mean PO4 = {R['mean_po4']:.3f} umol kg-1 ; "
          f"volume-source = {vol_source[s]}")
print()
print(f"comparison mask: {n_mask_cells} ocean cells retained of {n_mask_total}")
for m in MODELS:
    R = results[m]
    print(f"{m} vs WOA23: bias = {R['bias_vs_woa']:+.3f}, RMSE = {R['rmse_vs_woa']:.3f}, "
          f"centred RMSE = {R['crmse_vs_woa']:.3f}, max |err| = {R['max_abs_err_vs_woa']:.2f} "
          f"at {R['max_abs_err_at']}, n = {R['n_valid_masked']}")

# scatter samples for the figure (subsampled for drawing only)
rng = np.random.default_rng(0)
scatter_data = {}
for s in SOURCES:
    no3, po4 = fields[s]
    x = np.asarray(po4.values, dtype=float).ravel()
    y = np.asarray(no3.values, dtype=float).ravel()
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size > 20000:
        sel = rng.choice(x.size, 20000, replace=False)
        x, y = x[sel], y[sel]
    scatter_data[s] = (x, y)

print("\ndone")
