import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# Helper to preserve provenance attributes through derived quantities
def carry_attrs(new, old):
    attrs = dict(old.attrs)
    attrs.update(new.attrs)
    new.attrs = attrs
    return new

def get_coeffs(fit):
    co = fit['coefficients']
    if isinstance(co, dict):
        b = None
        a = None
        for k, v in co.items():
            kl = k.lower()
            if 'intercept' in kl or kl in ['const', 'b']:
                b = float(v)
        for k, v in co.items():
            kl = k.lower()
            if 'intercept' not in kl and kl not in ['const', 'b']:
                a = float(v)
                break
        if a is None:
            for k, v in co.items():
                if 'intercept' not in k.lower():
                    a = float(v)
                    break
        if b is None:
            b = 0.0
        return a, b
    else:
        arr = np.asarray(co).ravel()
        if arr.size == 2:
            b = float(arr[0])
            a = float(arr[1])
        elif arr.size == 1:
            a = float(arr[0])
            b = 0.0
        else:
            b = float(arr[0])
            a = float(arr[1])
        return a, b

sources = ['WOA23', 'CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']
models = ['CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']

stats = {}
regression = {}
scatter = {}
depth_maps = {}

# -----------------------------------------------------------------------------
# WOA23 observed climatology
# -----------------------------------------------------------------------------
woa_no3 = load_climatology_data(dataset='woa23', variable='nitrate',
                                time_slice='annual', depth_range=[200, 1000])
woa_po4 = load_climatology_data(dataset='woa23', variable='phosphate',
                                time_slice='annual', depth_range=[200, 1000])
woa_no3, woa_po4 = xr.align(woa_no3, woa_po4, join='inner')
woa_no3 = woa_no3.compute()
woa_po4 = woa_po4.compute()
woa_no3.name = 'NO3'
woa_po4.name = 'PO4'
woa_ratio = compute_ratio(woa_no3, woa_po4).compute()
woa_ratio = carry_attrs(woa_ratio, woa_no3)
woa_ratio.name = 'NP_ratio'

mean_ratio_woa = float(compute_basic_stat(woa_ratio, op='volume_weighted_mean', dim=None).compute())
mean_no3_woa = float(compute_basic_stat(woa_no3, op='volume_weighted_mean', dim=None).compute())
mean_po4_woa = float(compute_basic_stat(woa_po4, op='volume_weighted_mean', dim=None).compute())
stats['WOA23'] = {
    'mean_ratio': mean_ratio_woa,
    'ratio_of_means': mean_no3_woa / mean_po4_woa
}

mask = np.isfinite(woa_no3.values) & np.isfinite(woa_po4.values) & (woa_po4.values > 0)
x = woa_po4.values[mask]
y = woa_no3.values[mask]
scatter['WOA23'] = (x, y)
target = xr.DataArray(y, dims='sample', name='NO3')
pred = xr.DataArray(x, dims='sample', name='PO4')
fit_woa = fit_multivariate_relation(target=target, predictors=[pred],
                                    method='linear', validation='random')
a, b = get_coeffs(fit_woa)
regression['WOA23'] = {'a': a, 'b': b, 'r2': fit_woa['metrics']['r2'], 'fit': fit_woa}

band = compute_exceedance_fraction(woa_ratio, threshold=2, comparison='abs_less',
                                   reference=16, weighting='volume', dim=None)
stats['WOA23']['pct_14_18'] = float(band['fraction_percent'])
low = compute_exceedance_fraction(woa_ratio, threshold=12, comparison='less',
                                  weighting='volume', dim=None)
stats['WOA23']['pct_lt12_global'] = float(low['fraction_percent'])

for basin in ['Pacific', 'Atlantic', 'Indian']:
    basin_ratio = select_basin(woa_ratio, basin=basin,
                               include_marginal_seas=True, depth_dependent=False)
    bm = float(compute_basic_stat(basin_ratio, op='volume_weighted_mean', dim=None).compute())
    stats['WOA23'][f'mean_ratio_{basin}'] = bm
    if basin in ['Pacific', 'Indian']:
        bl = compute_exceedance_fraction(basin_ratio, threshold=12, comparison='less',
                                         weighting='volume', dim=None)
        stats['WOA23'][f'pct_lt12_{basin}'] = float(bl['fraction_percent'])

dm = compute_basic_stat(woa_ratio, op='volume_weighted_mean', dim='depth').compute()
dm = carry_attrs(dm, woa_ratio)
depth_maps['WOA23'] = dm

# -----------------------------------------------------------------------------
# CMIP6 historical models
# -----------------------------------------------------------------------------
for model in models:
    no3 = load_regular_data(dataset='cmip6_historical', variable='no3', model=model,
                            time_range=['1971-01', '2000-12'], depth_range=[200, 1000])
    po4 = load_regular_data(dataset='cmip6_historical', variable='po4', model=model,
                            time_range=['1971-01', '2000-12'], depth_range=[200, 1000])

    no3_mean = temporal_aggregate(no3, freq='mean', op='mean')
    po4_mean = temporal_aggregate(po4, freq='mean', op='mean')

    no3_umol = transform_unit(no3_mean, from_unit='mol m-3', to_unit='umol/kg',
                              density=1.025)
    po4_umol = transform_unit(po4_mean, from_unit='mol m-3', to_unit='umol/kg',
                              density=1.025)
    no3_umol = carry_attrs(no3_umol, no3_mean).compute()
    po4_umol = carry_attrs(po4_umol, po4_mean).compute()
    no3_umol, po4_umol = xr.align(no3_umol, po4_umol, join='inner')
    no3_umol.name = 'NO3'
    po4_umol.name = 'PO4'

    ratio = compute_ratio(no3_umol, po4_umol).compute()
    ratio = carry_attrs(ratio, no3_umol)
    ratio.name = 'NP_ratio'

    mean_ratio = float(compute_basic_stat(ratio, op='volume_weighted_mean', dim=None).compute())
    mean_no3 = float(compute_basic_stat(no3_umol, op='volume_weighted_mean', dim=None).compute())
    mean_po4 = float(compute_basic_stat(po4_umol, op='volume_weighted_mean', dim=None).compute())
    stats[model] = {
        'mean_ratio': mean_ratio,
        'ratio_of_means': mean_no3 / mean_po4
    }

    mask = np.isfinite(no3_umol.values) & np.isfinite(po4_umol.values) & (po4_umol.values > 0)
    x = po4_umol.values[mask]
    y = no3_umol.values[mask]
    scatter[model] = (x, y)
    target = xr.DataArray(y, dims='sample', name='NO3')
    pred = xr.DataArray(x, dims='sample', name='PO4')
    fit = fit_multivariate_relation(target=target, predictors=[pred],
                                    method='linear', validation='random')
    a, b = get_coeffs(fit)
    regression[model] = {'a': a, 'b': b, 'r2': fit['metrics']['r2'], 'fit': fit}

    band = compute_exceedance_fraction(ratio, threshold=2, comparison='abs_less',
                                       reference=16, weighting='volume', dim=None)
    stats[model]['pct_14_18'] = float(band['fraction_percent'])
    low = compute_exceedance_fraction(ratio, threshold=12, comparison='less',
                                      weighting='volume', dim=None)
    stats[model]['pct_lt12_global'] = float(low['fraction_percent'])

    for basin in ['Pacific', 'Atlantic', 'Indian']:
        basin_ratio = select_basin(ratio, basin=basin,
                                   include_marginal_seas=True, depth_dependent=False)
        bm = float(compute_basic_stat(basin_ratio, op='volume_weighted_mean', dim=None).compute())
        stats[model][f'mean_ratio_{basin}'] = bm
        if basin in ['Pacific', 'Indian']:
            bl = compute_exceedance_fraction(basin_ratio, threshold=12, comparison='less',
                                             weighting='volume', dim=None)
            stats[model][f'pct_lt12_{basin}'] = float(bl['fraction_percent'])

    dm = compute_basic_stat(ratio, op='volume_weighted_mean', dim='depth').compute()
    dm = carry_attrs(dm, ratio)
    depth_maps[model] = dm

# -----------------------------------------------------------------------------
# Regrid model depth-mean N:P fields for display only
# -----------------------------------------------------------------------------
display_maps = {}
display_maps['WOA23'] = depth_maps['WOA23']
for model in models:
    dm_reg = regrid_horizontal(depth_maps[model], target='woa23', method='bilinear')
    display_maps[model] = dm_reg.compute()

# -----------------------------------------------------------------------------
# Print summary of all reported values
# -----------------------------------------------------------------------------
print("=" * 80)
print("N:P stoichiometry in 200-1000 m")
print("=" * 80)
for src in sources:
    s = stats[src]
    r = regression[src]
    print(f"\n{src}:")
    print(f"  Mean N:P (mean of cell-by-cell ratio)          = {s['mean_ratio']:.2f}")
    print(f"  Mean N:P (ratio of volume-weighted means)      = {s['ratio_of_means']:.2f}")
    print(f"  OLS: NO3 = {r['a']:.2f} * PO4 + {r['b']:.2f}, R2 = {r['r2']:.3f}")
    print(f"  % volume 14 < N:P < 18                         = {s['pct_14_18']:.1f}%")
    print(f"  Basin mean N:P: Pacific = {s['mean_ratio_Pacific']:.2f}, "
          f"Atlantic = {s['mean_ratio_Atlantic']:.2f}, Indian = {s['mean_ratio_Indian']:.2f}")
    print(f"  % volume N:P < 12: global = {s['pct_lt12_global']:.1f}%, "
          f"Pacific = {s['pct_lt12_Pacific']:.1f}%, Indian = {s['pct_lt12_Indian']:.1f}%")
    diff = s['mean_ratio'] - 16.0
    if abs(diff) < 0.005:
        rel = "at"
    elif diff > 0:
        rel = "above"
    else:
        rel = "below"
    print(f"  Relative to Redfield (16)                      = {rel} ({diff:+.2f})")

print("\n" + "=" * 80)
print("Low-N:P volume comparison (<12) across sources")
print("=" * 80)
for src in sources:
    s = stats[src]
    print(f"  {src:14s}: global {s['pct_lt12_global']:5.1f}% | "
          f"Pacific {s['pct_lt12_Pacific']:5.1f}% | Indian {s['pct_lt12_Indian']:5.1f}%")

# -----------------------------------------------------------------------------
# Prepare plotting limits and shared colour scale
# -----------------------------------------------------------------------------
all_x = np.concatenate([scatter[src][0] for src in sources])
all_y = np.concatenate([scatter[src][1] for src in sources])
xmax = float(np.nanpercentile(all_x, 99.5))
ymax = float(np.nanpercentile(all_y, 99.5))
ymax = max(ymax, 16.0 * xmax * 1.05)
xmax = np.ceil(xmax * 10.0) / 10.0
ymax = np.ceil(ymax / 5.0) * 5.0
if xmax <= 0:
    xmax = 1.0
if ymax <= 0:
    ymax = 1.0

arrays = []
for src in sources:
    arr = (display_maps[src] - 16.0).values.flatten()
    arrays.append(arr[np.isfinite(arr)])
pooled = np.concatenate(arrays)
vmax = float(np.nanpercentile(np.abs(pooled), 98))
if vmax == 0:
    vmax = 1.0

figure1_data = scatter
figure1_reg = regression
figure2_data = display_maps
fig1_xmax = xmax
fig1_ymax = ymax
fig2_vmax = vmax
