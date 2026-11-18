# Computation code only
models = ['CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']
vars_2002 = ['dissic', 'talk', 'ph', 'thetao', 'so']
vars_lt = ['dissic', 'thetao', 'so']

# Load CMIP6 historical fields for 2002 and 1971-2000
cmip_2002 = {m: {} for m in models}
cmip_lt = {m: {} for m in models}
for m in models:
    for v in vars_2002:
        cmip_2002[m][v] = load_regular_data(
            dataset='cmip6_historical', variable=v, model=m,
            time_range=['2002-01', '2002-12']
        )
    for v in vars_lt:
        cmip_lt[m][v] = load_regular_data(
            dataset='cmip6_historical', variable=v, model=m,
            time_range=['1971-01', '2000-12']
        )

# Annual 2002 means and day-weighted 1971-2000 means
cmip_2002_ann = {m: {} for m in models}
cmip_lt_mean = {m: {} for m in models}
for m in models:
    for v in vars_2002:
        cmip_2002_ann[m][v] = temporal_aggregate(
            data=cmip_2002[m][v], freq='mean', op='mean'
        )
    for v in vars_lt:
        cmip_lt_mean[m][v] = temporal_aggregate(
            data=cmip_lt[m][v], freq='mean', op='mean'
        )

# Convert model DIC and total alkalinity from mol m-3 to umol kg-1
model_dic_2002 = {}
model_talk_2002 = {}
model_ph_2002 = {}
model_dic_lt = {}
for m in models:
    dens_2002 = compute_potential_density_anomaly(
        temp_data=cmip_2002_ann[m]['thetao'],
        sal_data=cmip_2002_ann[m]['so'],
        ref_pressure=None
    )
    dens_lt = compute_potential_density_anomaly(
        temp_data=cmip_lt_mean[m]['thetao'],
        sal_data=cmip_lt_mean[m]['so'],
        ref_pressure=None
    )

    model_dic_2002[m] = transform_unit(
        data=cmip_2002_ann[m]['dissic'],
        from_unit=cmip_2002_ann[m]['dissic'].attrs.get('units', 'mol m-3'),
        to_unit='umol/kg',
        density=dens_2002
    )
    model_talk_2002[m] = transform_unit(
        data=cmip_2002_ann[m]['talk'],
        from_unit=cmip_2002_ann[m]['talk'].attrs.get('units', 'mol m-3'),
        to_unit='umol/kg',
        density=dens_2002
    )
    model_dic_lt[m] = transform_unit(
        data=cmip_lt_mean[m]['dissic'],
        from_unit=cmip_lt_mean[m]['dissic'].attrs.get('units', 'mol m-3'),
        to_unit='umol/kg',
        density=dens_lt
    )
    model_ph_2002[m] = cmip_2002_ann[m]['ph']

# Load GLODAP mapped climatology
glodap_tco2 = load_climatology_data(
    dataset='glodap', variable='TCO2', time_slice='annual'
)
glodap_talk = load_climatology_data(
    dataset='glodap', variable='TAlk', time_slice='annual'
)
glodap_ph = load_climatology_data(
    dataset='glodap', variable='pHtsinsitutp', time_slice='annual'
)

def val(x):
    if hasattr(x, 'compute'):
        x = x.compute()
    return float(x)

# Global 0-50 m volume-weighted means
dic_means = {}
talk_means = {}
ph_means = {}
for m in models:
    dic_means[m] = val(compute_basic_stat(
        data=model_dic_2002[m].sel(depth=slice(0, 50)),
        op='volume_weighted_mean', dim=None
    ))
    talk_means[m] = val(compute_basic_stat(
        data=model_talk_2002[m].sel(depth=slice(0, 50)),
        op='volume_weighted_mean', dim=None
    ))
    ph_means[m] = val(compute_basic_stat(
        data=model_ph_2002[m].sel(depth=slice(0, 50)),
        op='volume_weighted_mean', dim=None
    ))

dic_means['GLODAP'] = val(compute_basic_stat(
    data=glodap_tco2.sel(depth=slice(0, 50)),
    op='volume_weighted_mean', dim=None
))
talk_means['GLODAP'] = val(compute_basic_stat(
    data=glodap_talk.sel(depth=slice(0, 50)),
    op='volume_weighted_mean', dim=None
))
ph_means['GLODAP'] = val(compute_basic_stat(
    data=glodap_ph.sel(depth=slice(0, 50)),
    op='volume_weighted_mean', dim=None
))

dic_diff = {m: dic_means[m] - dic_means['GLODAP'] for m in models}
talk_diff = {m: talk_means[m] - talk_means['GLODAP'] for m in models}
ph_diff = {m: ph_means[m] - ph_means['GLODAP'] for m in models}

# 0-50 m DIC maps (each cell's own 0-50 m vertical mean)
dic_map_0_50 = {}
for m in models:
    dic_map_0_50[m] = compute_basic_stat(
        data=model_dic_2002[m].sel(depth=slice(0, 50)),
        op='volume_weighted_mean', dim='depth'
    )
dic_map_0_50['GLODAP'] = compute_basic_stat(
    data=glodap_tco2.sel(depth=slice(0, 50)),
    op='volume_weighted_mean', dim='depth'
)

# Regrid model DIC maps to the GLODAP 1-degree grid
dic_map_regrid = {}
for m in models:
    dic_map_regrid[m] = regrid_horizontal(
        data=dic_map_0_50[m],
        target=dic_map_0_50['GLODAP'],
        method='conservative'
    )

# Map comparison on common cells
mask_map = {}
dic_bias_map = {}
map_metrics = {}
for m in models:
    mask_map[m] = build_comparison_mask(
        fields=[dic_map_regrid[m], dic_map_0_50['GLODAP']],
        ocean_only=True
    )
    a = dic_map_regrid[m].where(mask_map[m])
    b = dic_map_0_50['GLODAP'].where(mask_map[m])

    met = compute_difference_metrics(
        data=a, reference=b, dim=['lat', 'lon'], weighting='area'
    )
    corr = analyze_correlation(
        var1=a, var2=b, dim=['lat', 'lon'],
        method='pearson', weighting='area'
    )

    map_metrics[m] = {
        'bias': val(met['bias']),
        'rmse': val(met['rmse']),
        'correlation': val(corr['correlation'])
    }
    dic_bias_map[m] = (
        dic_map_regrid[m] - dic_map_0_50['GLODAP']
    ).where(mask_map[m]).compute()

# Vertical DIC bias profiles at requested depths
depths = [0, 100, 500, 1000, 2000]
vertical_bias = {m: {} for m in models}
for m in models:
    for d in depths:
        model_d = regrid_vertical(
            data=model_dic_2002[m], target=[d],
            method='linear', fill_surface=True
        ).squeeze('depth', drop=True)
        glodap_d = regrid_vertical(
            data=glodap_tco2, target=[d],
            method='linear', fill_surface=True
        ).squeeze('depth', drop=True)

        model_d_reg = regrid_horizontal(
            data=model_d, target=dic_map_0_50['GLODAP'],
            method='conservative'
        )
        mask_d = build_comparison_mask(
            fields=[model_d_reg, glodap_d], ocean_only=True
        )
        a = model_d_reg.where(mask_d)
        b = glodap_d.where(mask_d)

        met = compute_difference_metrics(
            data=a, reference=b, dim=['lat', 'lon'], weighting='area'
        )
        vertical_bias[m][d] = val(met['bias'])

# 1971-2000 global 0-50 m mean DIC and 2002 minus 1971-2000 change
dic_lt_mean_0_50 = {}
for m in models:
    dic_lt_mean_0_50[m] = val(compute_basic_stat(
        data=model_dic_lt[m].sel(depth=slice(0, 50)),
        op='volume_weighted_mean', dim=None
    ))

period_change = {m: dic_means[m] - dic_lt_mean_0_50[m] for m in models}
ratio = {
    m: period_change[m] / dic_diff[m] if dic_diff[m] != 0 else np.nan
    for m in models
}

# Summary tables
summary = pd.DataFrame({
    'DIC_mean': [dic_means[m] for m in models] + [dic_means['GLODAP']],
    'DIC_diff': [dic_diff[m] for m in models] + [np.nan],
    'TA_mean': [talk_means[m] for m in models] + [talk_means['GLODAP']],
    'TA_diff': [talk_diff[m] for m in models] + [np.nan],
    'pH_mean': [ph_means[m] for m in models] + [ph_means['GLODAP']],
    'pH_diff': [ph_diff[m] for m in models] + [np.nan],
}, index=models + ['GLODAP'])

print("Global 0-50 m means (DIC/TA: umol kg-1; pH: total scale)")
print(summary.round(3).to_string())

map_df = pd.DataFrame({m: map_metrics[m] for m in models}).T
print("\n0-50 m DIC map comparison on GLODAP 1-degree grid")
print(map_df.round(3).to_string())

vert_df = pd.DataFrame({m: vertical_bias[m] for m in models})
print("\nGlobal mean DIC bias (model - GLODAP) by depth (umol kg-1)")
print(vert_df.round(0).to_string())

period_df = pd.DataFrame({
    'DIC_2002_minus_1971_2000': [period_change[m] for m in models],
    'Ratio_to_2002_GLODAP_DIC_bias': [ratio[m] for m in models]
}, index=models)
print("\nPeriod change of global 0-50 m DIC and ratio to 2002 DIC bias")
print(period_df.round(2).to_string())

# Composite ranking of overall carbon-system reproduction
def safe_den(vals):
    v = np.nanmean(vals)
    if not np.isfinite(v) or v == 0:
        return 1.0
    return v

score = {}
for m in models:
    den_dic = safe_den([abs(dic_diff[mm]) for mm in models])
    den_talk = safe_den([abs(talk_diff[mm]) for mm in models])
    den_ph = safe_den([abs(ph_diff[mm]) for mm in models])
    den_rmse = safe_den([map_metrics[mm]['rmse'] for mm in models])

    score[m] = (
        abs(dic_diff[m]) / den_dic
        + abs(talk_diff[m]) / den_talk
        + abs(ph_diff[m]) / den_ph
        + map_metrics[m]['rmse'] / den_rmse
        + (1.0 - map_metrics[m]['correlation'])
    )

best_model = min(models, key=lambda m: score[m])
print(f"\nBest overall carbon-system reproduction (lowest composite score): {best_model}")

print("\nSign consistency of each model's three biases:")
for m in models:
    same_sign = (np.sign(dic_diff[m]) == np.sign(talk_diff[m]) == np.sign(ph_diff[m]))
    print(
        f"  {m}: DIC {dic_diff[m]:+.1f}, TA {talk_diff[m]:+.1f}, "
        f"pH {ph_diff[m]:+.3f}; same sign = {same_sign}"
    )

print("\nAcross-model sign pattern:")
for name, diffs in [('DIC', dic_diff), ('TA', talk_diff), ('pH', ph_diff)]:
    print(
        f"  {name}: " + ", ".join(f"{m} {diffs[m]:+.3f}" for m in models)
    )
