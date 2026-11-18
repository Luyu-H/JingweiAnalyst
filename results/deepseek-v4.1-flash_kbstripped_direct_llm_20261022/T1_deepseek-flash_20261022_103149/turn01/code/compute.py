import glob
import os
import re

import numpy as np
import pandas as pd
import xarray as xr

# -----------------------------------------------------------------------------
# Constants and helpers
# -----------------------------------------------------------------------------
TARGET = pd.period_range('1950-01', '2014-12', freq='M')
NINO_LAT_MIN, NINO_LAT_MAX = -5.0, 5.0
NINO_LON_MIN, NINO_LON_MAX = 190.0, 240.0   # 170W to 120W in 0-360 convention


def _period_index_from_time(da):
    t = da.time
    yrs = np.asarray(t.dt.year.values, dtype=int)
    mons = np.asarray(t.dt.month.values, dtype=int)
    return pd.PeriodIndex([pd.Period(year=y, month=m, freq='M')
                           for y, m in zip(yrs, mons)])


def _lat_lon_2d(ds, lat_name, lon_name):
    lat = ds[lat_name]
    lon = ds[lon_name]
    if lat.ndim == 1 and lon.ndim == 1:
        lat2, lon2 = xr.broadcast(lat, lon)
    else:
        lat2, lon2 = lat, lon
    lon_norm = ((lon2 % 360.0) + 360.0) % 360.0
    return lat2, lon_norm


def _nino_mask(lat2, lon2):
    return ((lat2 >= NINO_LAT_MIN) & (lat2 <= NINO_LAT_MAX) &
            (lon2 >= NINO_LON_MIN) & (lon2 <= NINO_LON_MAX))


def _regional_mean(da, w, ydim, xdim):
    """Area-weighted regional mean, ignoring missing cells month by month."""
    valid = da.notnull()
    num = (da.fillna(0.0) * w).sum(dim=[ydim, xdim])
    den = (valid * w).sum(dim=[ydim, xdim])
    return num / den


def _subset_time(ds, y0=1950, y1=2014):
    t = ds.time
    years = np.asarray(t.dt.year.values)
    idx = np.where((years >= y0) & (years <= y1))[0]
    return ds.isel(time=idx)


def _variant_from_files(files):
    for f in files:
        m = re.search(r'_historical_([^_]+)_', os.path.basename(f))
        if m:
            return m.group(1)
    return None


def _load_fixed(model, field, variant=None):
    root = '/data1/lorraine-20419/benchmark_cmip6/historical'
    files = sorted(glob.glob(f'{root}/{model}/{field}/*.nc'))
    if not files:
        return None
    if variant is not None:
        pref = [f for f in files if variant in os.path.basename(f)]
        if pref:
            files = pref
    ds = xr.open_dataset(files[0], use_cftime=True)
    da = ds[field] if field in ds.variables else ds[list(ds.data_vars)[0]]
    if 'time' in da.dims:
        da = da.isel(time=0, drop=True)
    return da.squeeze()


def _series_from_regional(reg):
    pi = _period_index_from_time(reg)
    s = pd.Series(np.asarray(reg.compute().values).ravel(), index=pi)
    s = s[~s.index.duplicated(keep='first')].reindex(TARGET)
    return s.astype(float)


# -----------------------------------------------------------------------------
# ERSSTv5 observation: native 2-degree grid
# -----------------------------------------------------------------------------
ersst_path = '/data1/lorraine-20419/ocean_data/ERSSTv5/sst.mnmean.nc'
ds_obs = xr.open_dataset(ersst_path)
ds_obs = _subset_time(ds_obs, 1950, 2014)

obs_var = ds_obs['sst']
obs_lat2, obs_lon2 = _lat_lon_2d(ds_obs, 'lat', 'lon')
obs_mask = _nino_mask(obs_lat2, obs_lon2)
obs_w = (np.cos(np.deg2rad(obs_lat2)) * obs_mask).fillna(0.0)

obs_reg = _regional_mean(obs_var, obs_w, 'lat', 'lon')
obs_s = _series_from_regional(obs_reg)

# -----------------------------------------------------------------------------
# CMIP6 historical tos: CESM2, GFDL-ESM4, MPI-ESM1-2-LR
# -----------------------------------------------------------------------------
cmip6_root = '/data1/lorraine-20419/benchmark_cmip6/historical'
MODELS = ['CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']
MODEL_CFG = {
    'CESM2': dict(ydim='nlat', xdim='nlon', lat='lat', lon='lon'),
    'GFDL-ESM4': dict(ydim='lat', xdim='lon', lat='lat', lon='lon'),
    'MPI-ESM1-2-LR': dict(ydim='j', xdim='i', lat='latitude', lon='longitude'),
}

model_s = {}
for model in MODELS:
    cfg = MODEL_CFG[model]
    files = sorted(glob.glob(f'{cmip6_root}/{model}/tos/*.nc'))
    variant = _variant_from_files(files)

    ds = xr.open_mfdataset(files, combine='by_coords', use_cftime=True)
    ds = _subset_time(ds, 1950, 2014)
    tos = ds['tos']

    lat2, lon2 = _lat_lon_2d(ds, cfg['lat'], cfg['lon'])
    regmask = _nino_mask(lat2, lon2)

    area = _load_fixed(model, 'areacello', variant=variant)
    if area is None:
        area = xr.ones_like(lat2)

    # CESM2's raw areacello is global.  The other models publish ocean-masked
    # areas, but applying sftof where available is harmless and consistent.
    sftof = _load_fixed(model, 'sftof', variant=variant)
    if sftof is not None:
        area = area * (sftof > 0.0)

    w = (area * regmask).fillna(0.0)
    reg = _regional_mean(tos, w, cfg['ydim'], cfg['xdim'])
    model_s[model] = _series_from_regional(reg)

# -----------------------------------------------------------------------------
# Diagnostics: Nino 3.4 amplitude and seasonal phase locking
# -----------------------------------------------------------------------------
def _diagnose(s):
    clim = s.groupby(s.index.month).mean()
    month_num = pd.Series(s.index.month, index=s.index)
    anom = s - month_num.map(clim)

    std = float(anom.std(ddof=0))
    monthly_std = (anom.groupby(anom.index.month)
                       .std(ddof=0)
                       .reindex(range(1, 13)))

    peak_month = int(monthly_std.idxmax())
    min_month = int(monthly_std.idxmin())
    peak_value = float(monthly_std.max())
    min_value = float(monthly_std.min())
    peak_to_min_ratio = float(peak_value / min_value) if min_value != 0 else np.nan

    exceed_count = int((anom.abs() > 1.0).sum())
    exceed_pct = float(100.0 * exceed_count / len(anom))

    return dict(
        series=s,
        anomaly=anom,
        std=std,
        monthly_std=monthly_std,
        peak_month=peak_month,
        peak_value=peak_value,
        min_month=min_month,
        min_value=min_value,
        peak_to_min_ratio=peak_to_min_ratio,
        exceed_count=exceed_count,
        exceed_pct=exceed_pct,
    )


obs = _diagnose(obs_s)
models = {m: _diagnose(model_s[m]) for m in MODELS}
all_diag = {'ERSSTv5': obs, **models}

# Standard deviation over all 780 months.
std_table = pd.DataFrame([
    {'source': name, 'std_C': d['std']}
    for name, d in all_diag.items()
]).set_index('source')

# Standard deviation by calendar month.
monthly_std_table = pd.DataFrame(
    {name: d['monthly_std'] for name, d in all_diag.items()}
)
monthly_std_table.index.name = 'month'

# Peak/minimum seasonal amplitude and exceedance statistics.
peak_table = pd.DataFrame([
    {
        'source': name,
        'peak_month': d['peak_month'],
        'peak_std_C': d['peak_value'],
        'min_month': d['min_month'],
        'min_std_C': d['min_value'],
        'peak_to_min_ratio': d['peak_to_min_ratio'],
        'exceed_count': d['exceed_count'],
        'exceed_pct': d['exceed_pct'],
    }
    for name, d in all_diag.items()
]).set_index('source')

# Model amplitude bias and ratio relative to ERSSTv5.
obs_std = obs['std']
bias_table = pd.DataFrame([
    {
        'model': m,
        'std_C': models[m]['std'],
        'amplitude_bias_C': models[m]['std'] - obs_std,
        'amplitude_ratio': models[m]['std'] / obs_std,
    }
    for m in MODELS
]).set_index('model')
bias_table['abs_bias_C'] = bias_table['amplitude_bias_C'].abs()
ranking_by_abs_bias = bias_table.sort_values('abs_bias_C').index.tolist()

# Explicit phase-locking comparison against the observation.
phase_locking = pd.DataFrame([
    {
        'source': name,
        'peak_month': d['peak_month'],
        'peak_std_C': d['peak_value'],
        'same_peak_month_as_obs': bool(d['peak_month'] == obs['peak_month']),
        'peak_to_min_ratio': d['peak_to_min_ratio'],
        'ratio_relative_to_obs': d['peak_to_min_ratio'] / obs['peak_to_min_ratio'],
    }
    for name, d in all_diag.items()
]).set_index('source')

# Convenient objects for plotting and for writing the answer.
months = np.arange(1, 13)
source_order = ['ERSSTv5', 'CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']
monthly_std = {name: monthly_std_table[name].values for name in source_order}
anomalies = {name: all_diag[name]['anomaly'] for name in source_order}
std_values = {name: all_diag[name]['std'] for name in source_order}
exceed_pct_values = {name: all_diag[name]['exceed_pct'] for name in source_order}
exceed_counts = {name: all_diag[name]['exceed_count'] for name in source_order}
peak_months = {name: all_diag[name]['peak_month'] for name in source_order}
peak_values = {name: all_diag[name]['peak_value'] for name in source_order}
