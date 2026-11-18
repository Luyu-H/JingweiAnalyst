import glob
import os
import re

CMIP6_ROOT = '/data1/lorraine-20419/benchmark_cmip6/historical'
ERSST_PATH = '/data1/lorraine-20419/ocean_data/ERSSTv5/sst.mnmean.nc'

model_cfg = {
    'CESM2': {
        'grid': 'gn',
        'y_dim': 'nlat',
        'x_dim': 'nlon',
        'lat_coord': 'lat',
        'lon_coord': 'lon',
        'depth_coord': 'lev',
        'depth_scale': 0.01,
        'calendar': 'noleap',
        'ocean_mask_from': 'sftof',
    },
    'GFDL-ESM4': {
        'grid': 'gr',
        'y_dim': 'lat',
        'x_dim': 'lon',
        'lat_coord': 'lat',
        'lon_coord': 'lon',
        'depth_coord': 'lev',
        'depth_scale': 1.0,
        'calendar': 'noleap',
    },
    'MPI-ESM1-2-LR': {
        'grid': 'gn',
        'y_dim': 'j',
        'x_dim': 'i',
        'lat_coord': 'latitude',
        'lon_coord': 'longitude',
        'depth_coord': 'lev',
        'depth_scale': 1.0,
        'calendar': 'proleptic_gregorian',
    },
}


def parse_variant(path):
    m = re.search(r'_(r\d+i\d+p\d+f\d+)_', os.path.basename(path))
    return m.group(1) if m else None


def find_files(model, var, table, grid, prefer_variant=None):
    pattern = f"{var}_{table}_{model}_historical_*_{grid}*.nc"
    files = sorted(glob.glob(os.path.join(CMIP6_ROOT, model, var, pattern)))
    if not files:
        raise FileNotFoundError(f"No files found for pattern {pattern}")
    groups = {}
    for f in files:
        v = parse_variant(f)
        groups.setdefault(v, []).append(f)
    if prefer_variant is not None and prefer_variant in groups:
        chosen = prefer_variant
    else:
        chosen = sorted(groups.keys(), key=lambda v: (-len(groups[v]), str(v)))[0]
    return groups[chosen], chosen


def open_cmip6_files(files):
    if len(files) == 1:
        return xr.open_dataset(files[0], use_cftime=True, chunks={'time': 120})
    return xr.open_mfdataset(files, combine='by_coords', use_cftime=True,
                             chunks={'time': 120})


def select_period(ds):
    tv = ds['time'].values
    ym = np.array([t.year * 12 + t.month for t in tv])
    mask = (ym >= 1950 * 12 + 1) & (ym <= 2014 * 12 + 12)
    idx = np.where(mask)[0]
    ds = ds.isel(time=idx)
    if ds.sizes['time'] != 780:
        raise ValueError(f"Expected 780 time steps, got {ds.sizes['time']}")
    return ds


def load_fixed(model, var, prefer_variant=None):
    cfg = model_cfg[model]
    files, _ = find_files(model, var, 'Ofx', cfg['grid'],
                          prefer_variant=prefer_variant)
    ds = xr.open_dataset(files[0], use_cftime=True)
    da = ds[var]
    for d in list(da.dims):
        if da.sizes[d] == 1:
            da = da.isel({d: 0})
    return da


# ------------------------------------------------------------------
# Regional mean Nino 3.4: 5S-5N, 170W-120W = 190E-240E
# ------------------------------------------------------------------
source_data = {}

# ERSSTv5: cosine-latitude area weights on the native 2-degree grid
ersst_ds = xr.open_dataset(ERSST_PATH, use_cftime=True)
ersst_ds = select_period(ersst_ds)
ersst_sst = ersst_ds['sst']

lat = ersst_ds['lat']
lon = ersst_ds['lon']
lon360 = lon % 360
region_mask = (lat >= -5) & (lat <= 5) & (lon360 >= 190) & (lon360 <= 240)
region_mask = region_mask.reset_coords(drop=True)

weights = np.cos(np.deg2rad(lat)).reset_coords(drop=True)
w = (weights * region_mask).fillna(0.0)

valid = ersst_sst.notnull()
num = (ersst_sst * w).sum(dim=['lat', 'lon'], skipna=True)
den = (w * valid).sum(dim=['lat', 'lon'], skipna=True)
ersst_nino34 = (num / den).compute()
source_data['ERSSTv5'] = (ersst_nino34, 'obs')

# CMIP6 tos: native grid and published areacello
for model in ['CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']:
    cfg = model_cfg[model]
    files, variant = find_files(model, 'tos', 'Omon', cfg['grid'])
    ds = open_cmip6_files(files)
    ds = ds.sortby('time')
    ds = select_period(ds)

    da = ds['tos']
    if 'lev' in da.dims:
        if da.sizes['lev'] > 1:
            da = da.isel(lev=0)
        else:
            da = da.squeeze('lev', drop=True)

    y_dim = cfg['y_dim']
    x_dim = cfg['x_dim']

    lat = ds[cfg['lat_coord']]
    lon = ds[cfg['lon_coord']]
    lon360 = lon % 360
    region_mask = (lat >= -5) & (lat <= 5) & (lon360 >= 190) & (lon360 <= 240)
    region_mask = region_mask.reset_coords(drop=True)

    weights = load_fixed(model, 'areacello',
                         prefer_variant=variant).fillna(0.0).reset_coords(drop=True)

    if model == 'CESM2':
        sftof = load_fixed(model, 'sftof',
                           prefer_variant=variant).fillna(0.0).reset_coords(drop=True)
        weights = weights.where(sftof > 0, 0.0)
    else:
        weights = weights.where(weights > 0, 0.0)

    w = (weights * region_mask).fillna(0.0)

    valid = da.notnull()
    num = (da * w).sum(dim=[y_dim, x_dim], skipna=True)
    den = (w * valid).sum(dim=[y_dim, x_dim], skipna=True)
    regional = (num / den).compute()
    source_data[model] = (regional, variant)


# ------------------------------------------------------------------
# Nino 3.4 anomalies and statistics
# ------------------------------------------------------------------
sources_order = ['ERSSTv5', 'CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']
month_names_full = ['January', 'February', 'March', 'April', 'May', 'June',
                    'July', 'August', 'September', 'October', 'November', 'December']
month_names_abbr = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

results = {}
for name in sources_order:
    da, variant = source_data[name]
    vals = np.asarray(da.values, dtype=float)
    tv = da['time'].values
    years = np.array([t.year for t in tv])
    months = np.array([t.month for t in tv])

    if len(vals) != 780:
        raise ValueError(f"{name}: expected 780 values, got {len(vals)}")
    if not np.all((years >= 1950) & (years <= 2014)):
        raise ValueError(f"{name}: time outside 1950-2014")

    # Plain multi-year mean for each calendar month over 1950-2014
    clim = np.array([vals[months == m].mean() for m in range(1, 13)])
    anom = vals - clim[months - 1]

    sd = float(np.std(anom, ddof=0))
    monthly_sd = np.array([np.std(anom[months == m], ddof=0) for m in range(1, 13)])

    peak_idx = int(np.argmax(monthly_sd))
    min_idx = int(np.argmin(monthly_sd))
    peak_month = peak_idx + 1
    min_month = min_idx + 1
    peak_sd = float(monthly_sd[peak_idx])
    min_sd = float(monthly_sd[min_idx])
    season_ratio = float(peak_sd / min_sd)

    exceed_count = int(np.sum(np.abs(anom) > 1.0))
    exceed_pct = float(exceed_count / 780.0 * 100.0)

    time_years = np.array([t.year + (t.month - 0.5) / 12.0 for t in tv])

    results[name] = {
        'anom': anom,
        'time': tv,
        'time_years': time_years,
        'months': months,
        'years': years,
        'clim': clim,
        'sd': sd,
        'monthly_sd': monthly_sd,
        'peak_month': peak_month,
        'peak_sd': peak_sd,
        'min_month': min_month,
        'min_sd': min_sd,
        'ratio': season_ratio,
        'exceed_count': exceed_count,
        'exceed_pct': exceed_pct,
        'variant': variant,
    }

obs_sd = results['ERSSTv5']['sd']
model_stats = {}
for model in ['CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']:
    sd = results[model]['sd']
    bias = sd - obs_sd
    ratio = sd / obs_sd
    model_stats[model] = {
        'sd': sd,
        'bias': bias,
        'ratio': ratio,
        'abs_bias': abs(bias),
    }

ranking = sorted(model_stats.keys(), key=lambda m: model_stats[m]['abs_bias'])

monthly_sd_df = pd.DataFrame(
    {name: results[name]['monthly_sd'] for name in sources_order},
    index=month_names_abbr
)
monthly_sd_df.index.name = 'Calendar month'

print("CMIP6 variant used for tos:")
for name in ['CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']:
    print(f"  {name}: {results[name]['variant']}")

print("\nArea weighting:")
print("  ERSSTv5: cosine(latitude) on its native 2-degree grid, valid ocean cells only")
print("  CESM2: published areacello masked by Ofx sftof > 0")
print("  GFDL-ESM4: published areacello on the 1-degree gr grid, areacello > 0")
print("  MPI-ESM1-2-LR: published areacello, areacello > 0")

print("\nNino 3.4 anomaly SD over 780 months (1950-2014, ddof=0):")
for name in sources_order:
    print(f"  {name}: {results[name]['sd']:.2f} degC")

print("\nModel amplitude vs ERSSTv5:")
for model in ['CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']:
    st = model_stats[model]
    print(f"  {model}: SD={st['sd']:.2f} degC, bias={st['bias']:+.2f} degC, ratio={st['ratio']:.2f}")
print("  ranking (smallest to largest |bias|): " + ", ".join(ranking))

print("\nMonthly Nino 3.4 anomaly SD by calendar month (ddof=0):")
print(monthly_sd_df.round(2).to_string())

print("\nPeak and minimum seasonal SD:")
for name in sources_order:
    r = results[name]
    print(f"  {name}: peak {month_names_abbr[r['peak_month']-1]} ({r['peak_sd']:.2f} degC), "
          f"min {month_names_abbr[r['min_month']-1]} ({r['min_sd']:.2f} degC), "
          f"peak/min ratio {r['ratio']:.2f}")

print("\nExceedance of |Nino 3.4 anomaly| > 1.0 degC:")
for name in sources_order:
    r = results[name]
    print(f"  {name}: {r['exceed_count']} months, {r['exceed_pct']:.1f}% of 780")
