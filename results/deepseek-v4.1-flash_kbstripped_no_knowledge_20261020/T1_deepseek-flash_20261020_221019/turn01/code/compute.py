import numpy as np
import pandas as pd
import xarray as xr

time_range = ['1950-01', '2014-12']
models = ['CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']
source_order = ['ERSSTv5'] + models

# Load SST fields on each source's native grid, subset to a box containing Nino 3.4.
sst_fields = {}
sst_fields['ERSSTv5'] = load_regular_data(
    dataset='ersst',
    variable='sst',
    field='total',
    time_range=time_range,
    lat_range=[-10, 10],
    lon_range=[-180, -110],
)

for model in models:
    sst_fields[model] = load_regular_data(
        dataset='cmip6_historical',
        variable='tos',
        model=model,
        field='total',
        time_range=time_range,
        lat_range=[-10, 10],
        lon_range=[-180, -110],
    )


def get_nino_stats(label, sst):
    """Compute Nino 3.4 index and associated amplitude statistics."""
    idx = compute_enso_nino_index(
        data=sst,
        nino_region='nino3.4',
        depth_level=None,
        anomaly_freq='month',
        climatology='self',
        base_period=time_range,
    ).compute()

    std_all = float(
        compute_basic_stat(data=idx, op='std', dim='time').compute()
    )

    month_stds = []
    for m in range(1, 13):
        sub = idx.where(idx.time.dt.month == m, drop=True)
        s = float(
            compute_basic_stat(data=sub, op='std', dim='time').compute()
        )
        month_stds.append(s)
    month_stds = np.asarray(month_stds, dtype=float)

    peak_i = int(np.nanargmax(month_stds))
    min_i = int(np.nanargmin(month_stds))
    peak_month = peak_i + 1
    min_month = min_i + 1
    peak_val = float(month_stds[peak_i])
    min_val = float(month_stds[min_i])
    ratio = peak_val / min_val if min_val != 0 else np.nan

    vals = np.asarray(idx.values, dtype=float)
    exceed_count = int(np.sum(np.abs(vals) > 1.0))
    exceed_pct = 100.0 * exceed_count / vals.size

    return {
        'idx': idx,
        'std_all': std_all,
        'month_stds': month_stds,
        'peak_month': peak_month,
        'peak_val': peak_val,
        'min_month': min_month,
        'min_val': min_val,
        'ratio': ratio,
        'exceed_count': exceed_count,
        'exceed_pct': exceed_pct,
        'n_months': vals.size,
    }


stats = {label: get_nino_stats(label, sst_fields[label]) for label in source_order}
idx_dict = {label: stats[label]['idx'] for label in source_order}
std_all_dict = {label: stats[label]['std_all'] for label in source_order}
month_std_dict = {label: stats[label]['month_stds'] for label in source_order}

ersst_std = stats['ERSSTv5']['std_all']
model_bias = {}
model_ratio = {}
for model in models:
    model_bias[model] = stats[model]['std_all'] - ersst_std
    model_ratio[model] = stats[model]['std_all'] / ersst_std

ranking = sorted(models, key=lambda m: abs(model_bias[m]))

obs_peak = stats['ERSSTv5']['peak_month']
obs_ratio = stats['ERSSTv5']['ratio']
phase_lock_comparison = {}
for model in models:
    phase_lock_comparison[model] = {
        'peak_month': stats[model]['peak_month'],
        'peak_match': stats[model]['peak_month'] == obs_peak,
        'ratio': stats[model]['ratio'],
        'ratio_difference_from_obs': stats[model]['ratio'] - obs_ratio,
        'ratio_ratio_to_obs': (
            stats[model]['ratio'] / obs_ratio if obs_ratio != 0 else np.nan
        ),
    }

months_abbr = [
    'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'
]

print('Nino 3.4 monthly anomaly, 1950-01 to 2014-12, own 1950-2014 monthly climatology')
print('Source         | N months | std_all (degC) | exceed > |1| count | exceed %')
for label in source_order:
    s = stats[label]
    print(
        f"{label:14s} | {s['n_months']:3d}      | "
        f"{s['std_all']:6.2f}         | {s['exceed_count']:3d}              | "
        f"{s['exceed_pct']:5.1f}%"
    )

print('\nModel amplitude bias and ratio vs ERSSTv5:')
for model in models:
    print(
        f"{model:14s} | bias = {model_bias[model]:+.2f} degC | "
        f"ratio = {model_ratio[model]:.2f}"
    )
print('Ranking by smallest absolute amplitude bias:', ranking)

print('\nSeasonal (calendar-month) standard deviations (degC):')
print('Source         ' + ' '.join(f'{m:>5s}' for m in months_abbr))
for label in source_order:
    print(
        f"{label:14s} " +
        ' '.join(f"{v:5.2f}" for v in stats[label]['month_stds'])
    )

print('\nSeasonal phase-locking summary:')
for label in source_order:
    s = stats[label]
    print(
        f"{label:14s} | peak {months_abbr[s['peak_month']-1]} "
        f"{s['peak_val']:.2f} degC | min {months_abbr[s['min_month']-1]} "
        f"{s['min_val']:.2f} degC | max/min ratio {s['ratio']:.2f}"
    )

print('\nModel phase-locking comparison with ERSSTv5:')
for model in models:
    c = phase_lock_comparison[model]
    print(
        f"{model:14s} | peak month {months_abbr[c['peak_month']-1]} | "
        f"reproduces observed peak month: {c['peak_match']} | "
        f"ratio {c['ratio']:.2f} vs obs {obs_ratio:.2f}"
    )

# Convenience variables for visualization.
time_axis = idx_dict['ERSSTv5'].time
all_anom_values = np.concatenate([
    np.asarray(idx_dict[m].values, dtype=float) for m in source_order
])
ymin = float(np.nanmin(all_anom_values))
ymax = float(np.nanmax(all_anom_values))
