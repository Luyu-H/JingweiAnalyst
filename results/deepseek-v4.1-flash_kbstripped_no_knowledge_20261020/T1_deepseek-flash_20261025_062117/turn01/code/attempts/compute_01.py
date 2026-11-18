# Computation code only

# Common settings
time_range = ['1950-01', '2014-12']
base_period = ['1950-01', '2014-12']
lat_box = [-10.0, 10.0]
lon_box = [-180.0, -100.0]
months = np.arange(1, 13)

sources_order = ['ERSSTv5', 'CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']
model_names = ['CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']

source_specs = {
    'ERSSTv5': dict(dataset='ersst', variable='sst'),
    'CESM2': dict(dataset='cmip6_historical', variable='tos', model='CESM2'),
    'GFDL-ESM4': dict(dataset='cmip6_historical', variable='tos', model='GFDL-ESM4'),
    'MPI-ESM1-2-LR': dict(dataset='cmip6_historical', variable='tos', model='MPI-ESM1-2-LR'),
}

nino_series = {}
std_all = {}
seasonal_std = {}
peak_month = {}
peak_value = {}
min_month = {}
min_value = {}
season_ratio = {}
exceed_n_pos = {}
exceed_n_neg = {}
exceed_n = {}
exceed_pct = {}

for src in sources_order:
    spec = source_specs[src]
    load_kwargs = dict(
        dataset=spec['dataset'],
        variable=spec['variable'],
        time_range=time_range,
        lat_range=lat_box,
        lon_range=lon_box,
    )
    if 'model' in spec:
        load_kwargs['model'] = spec['model']

    field = load_regular_data(**load_kwargs)

    nino = compute_enso_nino_index(
        data=field,
        nino_region='nino3.4',
        anomaly_freq='month',
        climatology='self',
        base_period=base_period,
    ).compute()

    if nino.sizes.get('time', 0) != 780:
        raise ValueError(f"{src}: expected 780 months, got {nino.sizes.get('time')}")

    nino_series[src] = nino

    # All-month population standard deviation
    sd_all = compute_basic_stat(nino, op='std', dim='time').compute()
    std_all[src] = float(sd_all)

    # Calendar-month population standard deviations
    month_std_values = []
    for m in months:
        month_idx = np.where(nino.time.dt.month.values == m)[0]
        if len(month_idx) != 65:
            raise ValueError(f"{src}, month {m}: expected 65 values, got {len(month_idx)}")
        sub = nino.isel(time=month_idx)
        sd_m = compute_basic_stat(sub, op='std', dim='time').compute()
        month_std_values.append(float(sd_m))

    seasonal_std[src] = pd.Series(
        month_std_values,
        index=pd.Index(months, name='month'),
        name=src,
    )

    # Peak and minimum of the seasonal-amplitude cycle
    seas_da = xr.DataArray(
        month_std_values,
        coords={'month': months},
        dims='month',
        name=f'{src}_seasonal_std',
        attrs={'units': 'degC'},
    )
    emax = locate_extreme(seas_da, op='max')
    emin = locate_extreme(seas_da, op='min')

    idx_max = int(emax['indices']['month'])
    idx_min = int(emin['indices']['month'])
    peak_month[src] = int(months[idx_max])
    peak_value[src] = float(emax['value'])
    min_month[src] = int(months[idx_min])
    min_value[src] = float(emin['value'])
    season_ratio[src] = peak_value[src] / min_value[src]

    # Strict exceedance: anomaly > +1 and anomaly < -1
    ep = compute_exceedance_fraction(
        nino,
        threshold=1.0,
        comparison='greater',
        reference=0.0,
        dim='time',
        weighting='none',
    )
    # x < -1 is equivalent to x - (-2) < 1
    en = compute_exceedance_fraction(
        nino,
        threshold=1.0,
        comparison='less',
        reference=-2.0,
        dim='time',
        weighting='none',
    )
    ea = compute_exceedance_fraction(
        nino,
        threshold=1.0,
        comparison='abs_greater',
        reference=0.0,
        dim='time',
        weighting='none',
    )

    npos = int(ep['n_exceeding'])
    nneg = int(en['n_exceeding'])
    nabs = int(ea['n_exceeding'])
    if npos + nneg != nabs:
        raise ValueError(f"{src}: exceedance count mismatch: {npos}+{nneg} != {nabs}")
    if int(ep['n_valid']) != 780 or int(en['n_valid']) != 780:
        raise ValueError(f"{src}: expected 780 valid values for exceedance")

    exceed_n_pos[src] = npos
    exceed_n_neg[src] = nneg
    exceed_n[src] = npos + nneg
    exceed_pct[src] = (npos + nneg) / 780.0 * 100.0

# Model amplitude comparison against ERSSTv5
obs_std = std_all['ERSSTv5']
bias = {m: std_all[m] - obs_std for m in model_names}
ratio = {m: std_all[m] / obs_std for m in model_names}
ranking = sorted(model_names, key=lambda m: abs(bias[m]))

report_all_std = pd.Series(
    {src: std_all[src] for src in sources_order},
    name='std_all_C',
)

report_model_comparison = pd.DataFrame({
    'model': model_names,
    'std_C': [std_all[m] for m in model_names],
    'bias_C': [bias[m] for m in model_names],
    'ratio': [ratio[m] for m in model_names],
})
report_model_comparison['abs_bias_C'] = report_model_comparison['bias_C'].abs()
report_model_comparison['rank_abs_bias'] = report_model_comparison['abs_bias_C'].rank(method='min').astype(int)
report_model_comparison = report_model_comparison.sort_values('abs_bias_C').reset_index(drop=True)
ranking = report_model_comparison['model'].tolist()

report_seasonal_std = pd.DataFrame({
    src: seasonal_std[src] for src in sources_order
})

report_peaks = pd.DataFrame({
    'source': sources_order,
    'peak_month': [peak_month[src] for src in sources_order],
    'peak_std_C': [peak_value[src] for src in sources_order],
    'min_month': [min_month[src] for src in sources_order],
    'min_std_C': [min_value[src] for src in sources_order],
    'max_min_ratio': [season_ratio[src] for src in sources_order],
})
report_peaks['same_peak_month_as_ERSSTv5'] = report_peaks['peak_month'] == peak_month['ERSSTv5']
report_peaks['ratio_vs_ERSSTv5'] = report_peaks['max_min_ratio'] / season_ratio['ERSSTv5']

report_exceedance = pd.DataFrame({
    'source': sources_order,
    'n_gt_plus1': [exceed_n_pos[src] for src in sources_order],
    'n_lt_minus1': [exceed_n_neg[src] for src in sources_order],
    'n_abs_gt_1': [exceed_n[src] for src in sources_order],
    'percent_abs_gt_1': [exceed_pct[src] for src in sources_order],
})

# Rounded reporting copies
report_all_std_rounded = report_all_std.round(2)

report_model_comparison_rounded = report_model_comparison.copy()
report_model_comparison_rounded[['std_C', 'bias_C', 'ratio', 'abs_bias_C']] = (
    report_model_comparison_rounded[['std_C', 'bias_C', 'ratio', 'abs_bias_C']].round(2)
)

report_seasonal_std_rounded = report_seasonal_std.round(2)

report_peaks_rounded = report_peaks.copy()
report_peaks_rounded[['peak_std_C', 'min_std_C', 'max_min_ratio', 'ratio_vs_ERSSTv5']] = (
    report_peaks_rounded[['peak_std_C', 'min_std_C', 'max_min_ratio', 'ratio_vs_ERSSTv5']].round(2)
)

report_exceedance_rounded = report_exceedance.copy()
report_exceedance_rounded['percent_abs_gt_1'] = (
    report_exceedance_rounded['percent_abs_gt_1'].round(1)
)

print("All-month Nino3.4 anomaly standard deviation (degC):")
print(report_all_std_rounded)
print("\nModel amplitude comparison:")
print(report_model_comparison_rounded)
print("\nRanking (smallest to largest absolute amplitude bias):")
print(ranking)
print("\nSeasonal (calendar-month) standard deviations (degC):")
print(report_seasonal_std_rounded)
print("\nSeasonal amplitude peaks:")
print(report_peaks_rounded)
print("\nStrict exceedance of |Nino3.4 anomaly| > 1 degC:")
print(report_exceedance_rounded)
