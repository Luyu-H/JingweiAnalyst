models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
months = [3, 9]
results = {}
scatter_data = {}
map_data = {}

def get_lat_lon(field, loc):
    if loc is None:
        return None, None
    lat_name = field.attrs.get('lat_coord', 'lat')
    lon_name = field.attrs.get('lon_coord', 'lon')
    if isinstance(loc, dict):
        if lat_name in loc and lon_name in loc:
            return loc[lat_name], loc[lon_name]
        point = field.isel(loc)
        lat = point.coords[lat_name].values if lat_name in point.coords else point[lat_name].values
        lon = point.coords[lon_name].values if lon_name in point.coords else point[lon_name].values
        return float(lat), float(lon)
    if isinstance(loc, tuple):
        point = field.isel(dict(zip(field.dims, loc)))
        lat = point.coords[lat_name].values if lat_name in point.coords else point[lat_name].values
        lon = point.coords[lon_name].values if lon_name in point.coords else point[lon_name].values
        return float(lat), float(lon)
    return None, None

for model in models:
    thetao = load_regular_data(dataset="cmip6_historical", variable="thetao", model=model, time_range=["1971-01", "2000-12"])
    so = load_regular_data(dataset="cmip6_historical", variable="so", model=model, time_range=["1971-01", "2000-12"])
    mlotst = load_regular_data(dataset="cmip6_historical", variable="mlotst", model=model, time_range=["1971-01", "2000-12"])
    areacello = load_regular_data(dataset="cmip6_historical", variable="areacello", model=model)
    
    levels, bounds = load_depth_levels(dataset="cmip6_historical", model=model, variable="thetao")
    ref_depth = float(levels[0])
    
    mld_op = compute_mixed_layer_depth(
        temp_data=thetao, 
        sal_data=so, 
        ref_depth=ref_depth, 
        method='density', 
        density_threshold=0.03, 
        temperature_type='potential'
    )
    
    mld_op_clim = temporal_aggregate(data=mld_op, freq='climatology', op='mean')
    mlotst_clim = temporal_aggregate(data=mlotst, freq='climatology', op='mean')
    
    total_ocean_area = float(areacello.sum().compute())
    
    for month in months:
        mld_op_month = mld_op.where(mld_op.time.dt.month == month, drop=True)
        mlotst_month = mlotst.where(mlotst.time.dt.month == month, drop=True)
        valid_mask = mld_op_month.notnull().all('time') & mlotst_month.notnull().all('time')
        
        mld_op_masked = mld_op_clim.sel(month=month).where(valid_mask)
        mlotst_masked = mlotst_clim.sel(month=month).where(valid_mask)
        
        metrics = compute_difference_metrics(data=mld_op_masked, reference=mlotst_masked, dim=None, weighting='area')
        bias = float(metrics['bias'])
        rmse = float(metrics['rmse'])
        max_err = float(metrics['max_abs_error'])
        max_loc = metrics['max_abs_error_at']
        
        corr_res = analyze_correlation(var1=mld_op_masked, var2=mlotst_masked, dim=['lat','lon'], method='pearson')
        corr = float(corr_res['correlation'])
        
        diff = mld_op_masked - mlotst_masked
        diff = diff.assign_attrs(units='m', source_dataset='cmip6_historical', grid_type=mld_op_masked.attrs.get('grid_type'))
        exceed_res = compute_exceedance_fraction(data=diff, threshold=20, comparison='abs_greater', weighting='area')
        exceed_pct = float(exceed_res['fraction_percent'])
        
        mask_area = float(areacello.where(valid_mask).sum().compute())
        coverage_pct = 100.0 * mask_area / total_ocean_area
        
        lat_val, lon_val = get_lat_lon(mld_op_masked, max_loc)
        
        results[(model, month)] = {
            'bias': bias,
            'rmse': rmse,
            'corr': corr,
            'exceed_pct': exceed_pct,
            'coverage_pct': coverage_pct,
            'max_err': max_err,
            'max_lat': lat_val,
            'max_lon': lon_val,
            'max_loc_raw': max_loc
        }
        
        mld_op_np = mld_op_masked.compute().values
        mlotst_np = mlotst_masked.compute().values
        x = mld_op_np.ravel()
        y = mlotst_np.ravel()
        valid = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
        scatter_data[(model, month)] = (x[valid], y[valid])
        
        diff_regrid = regrid_horizontal(diff, target='woa23', method='conservative')
        map_data[(model, month)] = {
            'diff': diff_regrid.compute().values,
            'lon': diff_regrid.lon.values,
            'lat': diff_regrid.lat.values
        }

# Ranking
ranking = {}
for month in months:
    month_name = 'March' if month == 3 else 'September'
    models_sorted = sorted(models, key=lambda m: results[(m, month)]['rmse'])
    ranking[month_name] = models_sorted

# Results DataFrame
rows = []
for model in models:
    for month in months:
        r = results[(model, month)]
        rows.append({
            'Model': model,
            'Month': 'March' if month == 3 else 'September',
            'Bias (m)': r['bias'],
            'RMSE (m)': r['rmse'],
            'Correlation': r['corr'],
            '>20m (%)': r['exceed_pct'],
            'Coverage (%)': r['coverage_pct'],
            'Max Abs Error (m)': r['max_err'],
            'Max Error Lat': r['max_lat'],
            'Max Error Lon': r['max_lon']
        })
results_df = pd.DataFrame(rows)
