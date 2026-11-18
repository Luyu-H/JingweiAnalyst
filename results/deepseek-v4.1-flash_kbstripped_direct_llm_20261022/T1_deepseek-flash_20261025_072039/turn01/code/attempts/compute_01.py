import os
import glob
import re
import numpy as np
import pandas as pd
import xarray as xr

HIST_ROOT = '/data1/lorraine-20419/benchmark_cmip6/historical'
ERSST_PATH = '/data1/lorraine-20419/ocean_data/ERSSTv5/sst.mnmean.nc'

START_YEAR = 1950
END_YEAR = 2014
N_MONTHS = (END_YEAR - START_YEAR + 1) * 12
TIME_AXIS = pd.date_range(f'{START_YEAR}-01-01', f'{END_YEAR}-12-01', freq='MS')

MODEL_CFG = {
    'CESM2': {
        'y': 'nlat', 'x': 'nlon', 'lat': 'lat', 'lon': 'lon',
        'apply_sftof': True
    },
    'GFDL-ESM4': {
        'y': 'lat', 'x': 'lon', 'lat': 'lat', 'lon': 'lon',
        'apply_sftof': False
    },
    'MPI-ESM1-2-LR': {
        'y': 'j', 'x': 'i', 'lat': 'latitude', 'lon': 'longitude',
        'apply_sftof': False
    },
}
MODELS = ['CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']
SOURCES = ['ERSSTv5'] + MODELS


def _files(root, model, var):
    return sorted(glob.glob(os.path.join(root, model, var, '*.nc')))


def _variant(path):
    m = re.search(r'_(r\d+i\d+p\d+f\d+)_', os.path.basename(path))
    return m.group(1) if m else None


def _open_cmip_files(files):
    if not files:
        raise FileNotFoundError('No CMIP6 files found')
    if len(files) == 1:
        return xr.open_dataset(files[0], use_cftime=True, decode_times=True)
    return xr.open_mfdataset(
        files, combine='by_coords', use_cftime=True, decode_times=True
    )


def _select_years(da, y0=START_YEAR, y1=END_YEAR):
    da = da.sortby('time')
    yrs = da.time.dt.year
    good = (yrs >= y0) & (yrs <= y1)
    return da.isel(time=np.where(good.values)[0])


def _ensure_celsius(da):
    u = str(da.attrs.get('units', '')).strip().lower()
    if u in ('k', 'kelvin') or u.startswith('kelvin'):
        da = da - 273.15
        da.attrs = dict(da.attrs)
        da.attrs['units'] = 'degC'
    return da


def _region_mask(lat2d, lon2d, lat0=-5.0, lat1=5.0, lon0=190.0, lon1=240.0):
    lonn = ((lon2d % 360.0) + 360.0) % 360.0
    return (lat2d >= lat0) & (lat2d <= lat1) & (lonn >= lon0) & (lonn <= lon1)


def _stats_from_index(idx):
    idx = idx.sortby('time')

    clim = idx.groupby('time.month').mean(dim='time')
    anom = (idx.groupby('time.month') - clim).sortby('time')

    std_all = float(anom.std(dim='time').values)

    seasonal = anom.groupby('time.month').std(dim='time')
    seasonal = seasonal.reindex({seasonal.dims[0]: np.arange(1, 13)})
    seasonal_vals = seasonal.values.astype(float)

    peak_month = int(seasonal.idxmax().values.item())
    peak_std = float(seasonal.max().values)
    min_month = int(seasonal.idxmin().values.item())
    min_std = float(seasonal.min().values)
    peak_min_ratio = peak_std / min_std if min_std != 0 else np.nan

    exceed_count = int((np.abs(anom) > 1.0).sum(dim='time').values)
    n = int(anom.sizes['time'])
    exceed_pct = round(100.0 * exceed_count / n, 1) if n else np.nan

    return {
        'std': std_all,
        'seasonal_std': seasonal_vals,
        'peak_month': peak_month,
        'peak_std': peak_std,
        'min_month': min_month,
        'min_std': min_std,
        'peak_min_ratio': peak_min_ratio,
        'exceed_count': exceed_count,
        'exceed_pct': exceed_pct,
        'anom_values': anom.values.astype(float),
        'n_months': n,
    }


# ---------------------------------------------------------------------------
# ERSSTv5
# ---------------------------------------------------------------------------
ds_obs = xr.open_dataset(ERSST_PATH)
sst_obs = ds_obs['sst'].sel(
    time=slice(f'{START_YEAR}-01-01', f'{END_YEAR}-12-31')
)

lat_o = sst_obs['lat']
lon_o = sst_obs['lon']
lonn_o = ((lon_o % 360.0) + 360.0) % 360.0

lat_o_sel = lat_o[(lat_o >= -5.0) & (lat_o <= 5.0)]
lon_o_sel = lon_o[(lonn_o >= 190.0) & (lonn_o <= 240.0)]
sub_o = sst_obs.sel(lat=lat_o_sel, lon=lon_o_sel)

w_o = np.cos(np.deg2rad(sub_o['lat']))
w_o = w_o.broadcast_like(sub_o.isel(time=0, drop=True))
valid_o = sub_o.notnull()
w_o_valid = w_o.where(valid_o)

num_o = (sub_o * w_o_valid).sum(dim=('lat', 'lon'), skipna=True)
den_o = w_o_valid.sum(dim=('lat', 'lon'), skipna=True)
idx_obs = (num_o / den_o).compute()

nino34_index = {'ERSSTv5': idx_obs}
nino34_stats = {'ERSSTv5': _stats_from_index(idx_obs)}


# ---------------------------------------------------------------------------
# CMIP6 historical tos
# ---------------------------------------------------------------------------
for model in MODELS:
    cfg = MODEL_CFG[model]

    tos_files = _files(HIST_ROOT, model, 'tos')
    ds_tos = _open_cmip_files(tos_files)
    tos = _select_years(ds_tos['tos'])
    tos = _ensure_celsius(tos)

    area_files = _files(HIST_ROOT, model, 'areacello')
    if not area_files:
        raise FileNotFoundError(f'No areacello file for {model}')

    variant = _variant(tos_files[0])
    if variant is not None:
        cand = [f for f in area_files if variant in os.path.basename(f)]
        if cand:
            area_files = cand

    ds_area = xr.open_dataset(area_files[0])
    area = ds_area['areacello']
    if 'time' in area.dims:
        area = area.isel(time=0)
    area = area.squeeze()

    if cfg['apply_sftof']:
        sftof_files = _files(HIST_ROOT, model, 'sftof')
        if sftof_files:
            if variant is not None:
                cand = [f for f in sftof_files if variant in os.path.basename(f)]
                if cand:
                    sftof_files = cand
            ds_sftof = xr.open_dataset(sftof_files[0])
            sftof = ds_sftof['sftof']
            if 'time' in sftof.dims:
                sftof = sftof.isel(time=0)
            sftof = sftof.squeeze()
            area = area * sftof

    lat2d = tos[cfg['lat']]
    lon2d = tos[cfg['lon']]
    if lat2d.ndim == 1:
        lat2d, lon2d = xr.broadcast(lat2d, lon2d)

    rmask = _region_mask(lat2d, lon2d)
    weights = area.where(rmask).fillna(0.0)

    valid = tos.notnull()
    w = weights.where(valid)

    num = (tos * w).sum(dim=(cfg['y'], cfg['x']), skipna=True)
    den = w.sum(dim=(cfg['y'], cfg['x']), skipna=True)
    idx_m = (num / den).compute()

    nino34_index[model] = idx_m
    nino34_stats[model] = _stats_from_index(idx_m)


for src in SOURCES:
    assert nino34_stats[src]['n_months'] == N_MONTHS, (
        src, nino34_stats[src]['n_months'], N_MONTHS
    )


# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------
nino34_std = {s: nino34_stats[s]['std'] for s in SOURCES}
nino34_std_2dp = {s: round(nino34_stats[s]['std'], 2) for s in SOURCES}

nino34_seasonal_std = {
    s: nino34_stats[s]['seasonal_std'] for s in SOURCES
}
nino34_seasonal_std_2dp = {
    s: np.round(nino34_stats[s]['seasonal_std'], 2) for s in SOURCES
}

nino34_peak_month = {s: nino34_stats[s]['peak_month'] for s in SOURCES}
nino34_peak_std = {s: nino34_stats[s]['peak_std'] for s in SOURCES}
nino34_min_month = {s: nino34_stats[s]['min_month'] for s in SOURCES}
nino34_min_std = {s: nino34_stats[s]['min_std'] for s in SOURCES}
nino34_peak_min_ratio = {s: nino34_stats[s]['peak_min_ratio'] for s in SOURCES}

nino34_exceed_count = {s: nino34_stats[s]['exceed_count'] for s in SOURCES}
nino34_exceed_pct = {s: nino34_stats[s]['exceed_pct'] for s in SOURCES}
nino34_anom_values = {s: nino34_stats[s]['anom_values'] for s in SOURCES}

obs_std = nino34_std['ERSSTv5']
amplitude_bias = {m: nino34_std[m] - obs_std for m in MODELS}
amplitude_ratio = {m: nino34_std[m] / obs_std for m in MODELS}
amplitude_bias_2dp = {m: round(amplitude_bias[m], 2) for m in MODELS}
amplitude_ratio_2dp = {m: round(amplitude_ratio[m], 2) for m in MODELS}
amplitude_bias_ranking = sorted(MODELS, key=lambda m: abs(amplitude_bias[m]))

seasonal_std_df = pd.DataFrame(
    {s: nino34_stats[s]['seasonal_std'] for s in SOURCES},
    index=np.arange(1, 13)
)
seasonal_std_rounded_df = seasonal_std_df.round(2)

summary_df = pd.DataFrame({
    'nino34_std_degC': [nino34_std[s] for s in SOURCES],
    'nino34_std_degC_2dp': [nino34_std_2dp[s] for s in SOURCES],
    'exceed_count_abs_gt_1C': [nino34_exceed_count[s] for s in SOURCES],
    'exceed_pct': [nino34_exceed_pct[s] for s in SOURCES],
    'amplitude_bias_vs_ERSSTv5_degC': [np.nan] + [
        amplitude_bias_2dp[m] for m in MODELS
    ],
    'amplitude_ratio_vs_ERSSTv5': [np.nan] + [
        amplitude_ratio_2dp[m] for m in MODELS
    ],
}, index=SOURCES)

phase_locking_df = pd.DataFrame({
    'peak_month': [nino34_peak_month[s] for s in SOURCES],
    'peak_std_degC': [round(nino34_peak_std[s], 2) for s in SOURCES],
    'min_month': [nino34_min_month[s] for s in SOURCES],
    'min_std_degC': [round(nino34_min_std[s], 2) for s in SOURCES],
    'peak_to_min_ratio': [round(nino34_peak_min_ratio[s], 2) for s in SOURCES],
}, index=SOURCES)

obs_peak_month = nino34_peak_month['ERSSTv5']
obs_peak_min_ratio = nino34_peak_min_ratio['ERSSTv5']

phase_locking_assessment = {
    m: {
        'model_peak_month': nino34_peak_month[m],
        'obs_peak_month': obs_peak_month,
        'peak_month_matches_obs': nino34_peak_month[m] == obs_peak_month,
        'model_peak_to_min_ratio': round(nino34_peak_min_ratio[m], 2),
        'obs_peak_to_min_ratio': round(obs_peak_min_ratio, 2),
        'ratio_difference_model_minus_obs': round(
            nino34_peak_min_ratio[m] - obs_peak_min_ratio, 2
        ),
    }
    for m in MODELS
}
