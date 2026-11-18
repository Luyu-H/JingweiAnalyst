import os, re, glob
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.interpolate import griddata

# -----------------------------------------------------------------------------
# Paths and constants
# -----------------------------------------------------------------------------
WOA_ROOT = Path('/data1/lorraine-20419/ocean_data/WOA23/climatology')
WOA_MASKS = Path('/data1/lorraine-20419/ocean_data/WOA23/masks')
HIST_ROOT = Path('/data1/lorraine-20419/benchmark_cmip6/historical')
SSP_ROOT = Path('/data1/lorraine-20419/benchmark_cmip6/ssp245')

R_EARTH = 6371000.0
THRESH_20 = 20.0
THRESH_60 = 60.0

CMIP_MODELS = {
    'CESM2': {
        'grid_type': 'curvilinear',
        'y_dim': 'nlat', 'x_dim': 'nlon',
        'lat_coord': 'lat', 'lon_coord': 'lon',
        'depth_coord': 'lev', 'depth_scale': 0.01,
    },
    'GFDL-ESM4': {
        'grid_type': 'regular',
        'y_dim': 'lat', 'x_dim': 'lon',
        'lat_coord': 'lat', 'lon_coord': 'lon',
        'depth_coord': 'lev', 'depth_scale': 1.0,
    },
    'MPI-ESM1-2-LR': {
        'grid_type': 'curvilinear',
        'y_dim': 'j', 'x_dim': 'i',
        'lat_coord': 'latitude', 'lon_coord': 'longitude',
        'depth_coord': 'lev', 'depth_scale': 1.0,
    },
}

# -----------------------------------------------------------------------------
# Small helpers
# -----------------------------------------------------------------------------
def read_msk_table(path):
    """Read a WOA .msk ASCII table into a list of numeric rows."""
    rows = []
    with open(path, 'r', errors='ignore') as f:
        for line in f:
            s = line.strip()
            if not s or s[0] in '#%!':
                continue
            s = s.replace(',', ' ')
            parts = s.split()
            vals = []
            ok = True
            for p in parts:
                try:
                    vals.append(float(p))
                except ValueError:
                    ok = False
                    break
            if ok and len(vals) >= 3:
                rows.append(vals)
    return rows


def woa_var(ds, prefix):
    """Pick the objectively analysed field for a WOA variable from an open file."""
    for cand in (f'{prefix}_an', f'{prefix}_mn', prefix):
        if cand in ds.data_vars:
            return ds[cand]
    for v in ds.data_vars:
        if v.startswith(prefix) and v.endswith('_an'):
            return ds[v]
    return ds[list(ds.data_vars)[0]]


def open_woa_field(var_dir, period_dir, tag, prefix):
    path = WOA_ROOT / var_dir / period_dir / f'woa23_{tag}_{prefix}00_01.nc'
    ds = xr.open_dataset(path, decode_times=False)
    da = woa_var(ds, prefix)
    if 'time' in da.dims:
        da = da.isel(time=0)
    return da.squeeze(drop=True)


def files_for(root, model, var, exp=None, strict=False):
    d = Path(root) / model / var
    if not d.exists():
        return []
    fs = sorted(d.glob('*.nc'))
    if exp is not None:
        fs_exp = [f for f in fs if exp in f.name]
        if fs_exp:
            fs = fs_exp
        elif strict:
            return []
    return fs


def choose_single_variant(files):
    """If several ensemble variants are present, choose the one with most files."""
    tagged = []
    for f in files:
        m = re.search(r'(r\d+i\d+p\d+f\d+)', f.name)
        if m:
            tagged.append((m.group(1), f))
    if not tagged:
        return files
    counts = {}
    for v, f in tagged:
        counts[v] = counts.get(v, 0) + 1
    best = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
    return [f for v, f in tagged if v == best]


def pick_var(ds, name):
    if name in ds.data_vars:
        return ds[name]
    for v in ds.data_vars:
        if name in v:
            return ds[v]
    return ds[list(ds.data_vars)[0]]


def open_cmip_field(root, model, var, exp, y0=None, y1=None, fallback_root=None):
    fs = files_for(root, model, var, exp, strict=True)
    if not fs and fallback_root is not None:
        fs = files_for(fallback_root, model, var, exp, strict=True)
        if not fs:
            fs = files_for(fallback_root, model, var, None, strict=False)
    if not fs:
        raise FileNotFoundError(f'No files found for {model=} {var=} {exp=}')
    fs = choose_single_variant(fs)
    ds = xr.open_mfdataset(
        [str(f) for f in fs],
        combine='by_coords',
        use_cftime=True,
        chunks='auto',
    )
    da = pick_var(ds, var)
    if y0 is not None and 'time' in da.dims:
        da = time_mean(da, y0, y1)
    return da


def time_mean(da, y0, y1):
    if 'time' not in da.dims:
        return da
    yrs = da.time.dt.year
    sel = (yrs >= y0) & (yrs <= y1)
    return da.where(sel, drop=True).mean('time')


def in_situ_density(thetao, salinity, depth_m):
    """
    UNESCO 1983 in-situ density from temperature (deg C), practical salinity
    and depth (m). Pressure is approximated hydrostatically in dbar.
    """
    T = thetao.clip(min=-2.0, max=40.0)
    S = salinity.clip(min=0.0, max=42.0)
    P = 1.013 + 1.0055 * depth_m       # dbar

    rho_w = (
        999.842594
        + 6.793952e-2 * T
        - 9.095290e-3 * T**2
        + 1.001685e-4 * T**3
        - 1.120083e-6 * T**4
        + 6.536332e-9 * T**5
    )
    A = (
        8.24493e-1
        - 4.0899e-3 * T
        + 7.6438e-5 * T**2
        - 8.2467e-7 * T**3
        + 5.3875e-9 * T**4
    )
    B = -5.72466e-3 + 1.0227e-4 * T - 1.6546e-6 * T**2
    C = 4.8314e-4

    rho0 = rho_w + A * S + B * S**1.5 + C * S**2

    K_w = (
        19652.21
        + 148.4206 * T
        - 2.327105 * T**2
        + 1.360477e-2 * T**3
        - 5.155288e-5 * T**4
    )
    K_0 = (
        K_w
        + S * (54.6746 - 0.603459 * T + 1.09987e-2 * T**2 - 6.1670e-5 * T**3)
        + S**1.5 * (7.944e-2 + 1.6483e-2 * T - 5.3009e-4 * T**2)
    )
    A_p = 3.239908 + 1.43713e-3 * T + 1.16092e-4 * T**2 - 5.77905e-7 * T**3
    B_p = 2.2838e-3 - 1.0981e-5 * T - 1.6078e-6 * T**2
    C_p = 1.91075e-4

    K = (
        K_0
        + P * (A_p + S * B_p + S**1.5 * C_p)
        + P**2 * (
            8.50935e-5
            - 6.12293e-6 * T
            + 5.2787e-8 * T**2
            + S * (-9.9348e-7 + 2.0816e-8 * T + 9.1697e-10 * T**2)
        )
    )
    return rho0 / (1.0 - P / K)


def threshold_volume(field, vol, thr):
    """Volume (in 1e15 m3) where field is strictly below thr, missing cells excluded."""
    field, vol = xr.align(field, vol, join='inner')
    vol_b = vol.broadcast_like(field)
    valid = field.notnull() & vol_b.notnull() & (vol_b > 0)
    mask = valid & (field < thr)
    amount = vol_b.where(mask, 0.0).sum()
    return float(amount.compute()) / 1.0e15


def oxygen_inventory(o2_native, vol):
    """Global dissolved oxygen inventory in mol: sum(o2 * volcello)."""
    o2_native, vol = xr.align(o2_native, vol, join='inner')
    vol_b = vol.broadcast_like(o2_native)
    valid = o2_native.notnull() & vol_b.notnull() & (vol_b > 0)
    amount = (o2_native.where(valid, 0.0) * vol_b.where(valid, 0.0)).sum()
    return float(amount.compute())


def omz_thickness(field, thickness, depth_dim, thr=60.0):
    """Column-summed thickness of water with field < thr."""
    thick_b = thickness.broadcast_like(field)
    valid = field.notnull() & thick_b.notnull()
    cond = valid & (field < thr)
    return thick_b.where(cond, 0.0).sum(dim=depth_dim)


def max_location(da, lat, lon):
    arr = np.asarray(da.values, dtype=float)
    if not np.isfinite(arr).any():
        return np.nan, np.nan, np.nan
    k = np.nanargmax(arr)
    iy, ix = np.unravel_index(k, arr.shape)
    lat_arr = np.asarray(lat.values)
    lon_arr = np.asarray(lon.values)
    if lat_arr.ndim == 2:
        lat0 = lat_arr[iy, ix]
        lon0 = lon_arr[iy, ix]
    else:
        lat0 = lat_arr[iy]
        lon0 = lon_arr[ix]
    return float(arr[iy, ix]), float(lat0), float(lon0)


def latlon_to_xyz(lat, lon):
    latr = np.deg2rad(lat)
    lonr = np.deg2rad(lon)
    x = np.cos(latr) * np.cos(lonr)
    y = np.cos(latr) * np.sin(lonr)
    z = np.sin(latr)
    return np.column_stack([x.ravel(), y.ravel(), z.ravel()])


# -----------------------------------------------------------------------------
# 1. WOA23 masks and reference fields
# -----------------------------------------------------------------------------
# Land-sea mask
land_rows = read_msk_table(WOA_MASKS / 'landsea_01.msk')
land_rows = [r for r in land_rows if len(r) >= 3 and -90.0 <= r[0] <= 90.0 and -180.0 <= r[1] <= 180.0]
lat_ll = np.array([r[0] for r in land_rows], dtype=float)
lon_ll = np.array([r[1] for r in land_rows], dtype=float)
bottom_ll = np.array([r[2] for r in land_rows], dtype=float)

lat_ll_u = np.sort(np.unique(np.round(lat_ll, 5)))
lon_ll_u = np.sort(np.unique(np.round(lon_ll, 5)))
lat_ll_idx = {float(v): i for i, v in enumerate(lat_ll_u)}
lon_ll_idx = {float(v): i for i, v in enumerate(lon_ll_u)}

bottom_grid = np.full((len(lat_ll_u), len(lon_ll_u)), np.nan)
for la, lo, b in zip(lat_ll, lon_ll, bottom_ll):
    i = lat_ll_idx[float(np.round(la, 5))]
    j = lon_ll_idx[float(np.round(lo, 5))]
    bottom_grid[i, j] = b

bottom_da = xr.DataArray(
    bottom_grid,
    coords={'lat': lat_ll_u, 'lon': lon_ll_u},
    dims=('lat', 'lon'),
)

# Basin mask
basin_rows = read_msk_table(WOA_MASKS / 'basinmask_01.msk')
basin_rows = [r for r in basin_rows if len(r) >= 4 and -90.0 <= r[0] <= 90.0 and -180.0 <= r[1] <= 180.0]
nlev_b = max(len(r) - 2 for r in basin_rows)

lat_b = np.array([r[0] for r in basin_rows], dtype=float)
lon_b = np.array([r[1] for r in basin_rows], dtype=float)
lat_b_u = np.sort(np.unique(np.round(lat_b, 5)))
lon_b_u = np.sort(np.unique(np.round(lon_b, 5)))
lat_b_idx = {float(v): i for i, v in enumerate(lat_b_u)}
lon_b_idx = {float(v): i for i, v in enumerate(lon_b_u)}

basin_arr = np.full((nlev_b, len(lat_b_u), len(lon_b_u)), np.nan)
for r in basin_rows:
    i = lat_b_idx[float(np.round(r[0], 5))]
    j = lon_b_idx[float(np.round(r[1], 5))]
    vals = r[2:2 + nlev_b]
    if len(vals) < nlev_b:
        vals = vals + [np.nan] * (nlev_b - len(vals))
    basin_arr[:len(vals), i, j] = vals

# depth axis from WOA nitrate file
ds_depth = xr.open_dataset(
    WOA_ROOT / 'nitrate/all_time_periods/woa23_all_n00_01.nc',
    decode_times=False,
)
depth_name = 'depth' if 'depth' in ds_depth.coords else list(ds_depth.coords)[0]
z = np.asarray(ds_depth[depth_name].values, dtype=float)
z = np.sort(z)
nlev = len(z)

depth_b = z if nlev_b == nlev else np.arange(nlev_b, dtype=float)
basin_da = xr.DataArray(
    basin_arr,
    coords={'depth': depth_b, 'lat': lat_b_u, 'lon': lon_b_u},
    dims=('depth', 'lat', 'lon'),
).where(lambda x: (x > 0) & (x < 100))

# WOA fields on the 1971-2000 normal
o2_woa = open_woa_field('oxygen', '1971-2000_climate_normal', 'decav71A0', 'o')
aou_woa = open_woa_field('apparent_oxygen_utilization', '1971-2000_climate_normal', 'decav71A0', 'A')
o2sat_woa = open_woa_field('percent_oxygen_saturation', '1971-2000_climate_normal', 'decav71A0', 'O')

o2_woa = o2_woa.sortby('lat').sortby('lon').sortby('depth')
aou_woa = aou_woa.reindex_like(o2_woa, method='nearest')
o2sat_woa = o2sat_woa.reindex_like(o2_woa, method='nearest')

# WOA cell geometry: layer thickness and spherical cell area
bounds = np.empty(nlev + 1, dtype=float)
bounds[0] = 0.0
bounds[1:-1] = 0.5 * (z[:-1] + z[1:])
bounds[-1] = z[-1] + 0.5 * (z[-1] - z[-2])
thick_woa = np.diff(bounds)
thick_woa_da = xr.DataArray(
    thick_woa,
    coords={'depth': o2_woa.depth},
    dims=('depth',),
)

lat_woa = np.asarray(o2_woa['lat'].values, dtype=float)
lon_woa = np.asarray(o2_woa['lon'].values, dtype=float)
dlat_rad = np.deg2rad(abs(np.mean(np.diff(lat_woa))))
dlon_rad = np.deg2rad(abs(np.mean(np.diff(lon_woa))))
area_woa_vals = R_EARTH**2 * np.cos(np.deg2rad(lat_woa)) * dlat_rad * dlon_rad
area_woa = xr.DataArray(
    area_woa_vals[:, None] * np.ones((1, len(lon_woa))),
    coords={'lat': o2_woa.lat, 'lon': o2_woa.lon},
    dims=('lat', 'lon'),
)
vol_woa = thick_woa_da * area_woa
vol_woa = vol_woa.reindex(depth=o2_woa.depth, lat=o2_woa.lat, lon=o2_woa.lon)

# WOA ocean level mask from landsea bottom levels
bottom_woa = bottom_da.reindex(lat=o2_woa.lat, lon=o2_woa.lon, method='nearest')
deepest_idx = np.where(
    bottom_woa.values > 1.0,
    np.floor(bottom_woa.values - 1.0),
    0.0,
).astype(int)
deepest_idx = np.clip(deepest_idx, 0, nlev)
lev_num = np.arange(nlev) + 1
ocean_level_vals = lev_num[:, None, None] <= deepest_idx[None, :, :]
ocean_level_woa = xr.DataArray(
    ocean_level_vals,
    coords={'depth': o2_woa.depth, 'lat': o2_woa.lat, 'lon': o2_woa.lon},
    dims=('depth', 'lat', 'lon'),
)

# Basin surface labels and marginal-sea folding
basin_da = basin_da.reindex(
    depth=o2_woa.depth,
    lat=o2_woa.lat,
    lon=o2_woa.lon,
    method='nearest',
)
basin_surf = basin_da.isel(depth=0)


def fold_basin_surface(surf):
    vals = surf.values.astype(float)
    out = np.full(vals.shape, np.nan)
    labels = np.unique(vals[np.isfinite(vals) & (vals > 0) & (vals < 100)])
    lat = np.asarray(surf['lat'].values, dtype=float)
    lon = np.asarray(surf['lon'].values, dtype=float)
    LON, LAT = np.meshgrid(lon, lat)
    for code in labels:
        m = vals == code
        mlat = np.nanmean(LAT[m])
        lonr = np.deg2rad(LON[m])
        mlon = np.rad2deg(
            np.arctan2(np.nanmean(np.sin(lonr)), np.nanmean(np.cos(lonr)))
        ) % 360.0
        code_int = int(round(float(code)))
        if mlat < -50.0:
            parent = 5      # Southern
        elif mlat > 60.0:
            parent = 4      # Arctic
        elif code_int in (1, 2, 3):
            parent = code_int
        elif code_int in (6, 7, 8):
            parent = 1      # Mediterranean / Black / Baltic -> Atlantic
        elif code_int in (9, 10):
            parent = 3      # Red Sea / Persian Gulf -> Indian
        elif 20.0 <= mlon <= 110.0 and mlat < 35.0:
            parent = 3
        elif 110.0 < mlon <= 290.0:
            parent = 2
        else:
            parent = 1
        out[m] = parent
    return xr.DataArray(out, coords=surf.coords, dims=surf.dims)


basin_surf_folded = fold_basin_surface(basin_surf)


def map_basin_to_grid(basin_surf, target_lat, target_lon):
    lat_s = np.asarray(basin_surf['lat'].values, dtype=float)
    lon_s = np.asarray(basin_surf['lon'].values, dtype=float)
    vals = basin_surf.values.astype(float)
    valid = np.isfinite(vals) & (vals > 0)
    LON_S, LAT_S = np.meshgrid(lon_s, lat_s)
    xyz_s = latlon_to_xyz(LAT_S[valid], LON_S[valid])
    tree = cKDTree(xyz_s)

    if target_lat.ndim == 1 and target_lon.ndim == 1:
        TLON, TLAT = np.meshgrid(target_lon.values, target_lat.values)
        dims = (target_lat.dims[0], target_lon.dims[0])
        coords = {target_lat.name: target_lat, target_lon.name: target_lon}
    else:
        TLAT = np.asarray(target_lat.values)
        TLON = np.asarray(target_lon.values)
        dims = target_lat.dims
        coords = {target_lat.name: target_lat, target_lon.name: target_lon}

    xyz_t = latlon_to_xyz(TLAT, TLON)
    _, idx = tree.query(xyz_t)
    mapped = vals[valid][idx].reshape(TLAT.shape)
    return xr.DataArray(mapped, coords=coords, dims=dims)


# WOA source bundle
woa_source = {
    'name': 'WOA23',
    'model': 'WOA23',
    'exp': '1971-2000',
    'o2_umol': o2_woa,
    'o2_native': None,
    'vol': vol_woa,
    'area': area_woa,
    'thickness': thick_woa_da,
    'lat': o2_woa['lat'],
    'lon': o2_woa['lon'],
    'depth_dim': 'depth',
    'lat_coord': 'lat',
    'lon_coord': 'lon',
    'grid_type': 'regular',
    'ocean_mask': ocean_level_woa,
}

# -----------------------------------------------------------------------------
# 2. CMIP6 loading
# -----------------------------------------------------------------------------
def load_cmip_source(model, exp, y0, y1):
    info = CMIP_MODELS[model]
    root = SSP_ROOT if exp == 'ssp245' else HIST_ROOT

    o2_da = open_cmip_field(root, model, 'o2', exp)
    thetao_da = open_cmip_field(root, model, 'thetao', exp)
    so_da = open_cmip_field(root, model, 'so', exp)

    o2_mean = time_mean(o2_da, y0, y1).squeeze(drop=True)
    thetao_mean = time_mean(thetao_da, y0, y1).squeeze(drop=True)
    so_mean = time_mean(so_da, y0, y1).squeeze(drop=True)

    o2_mean, thetao_mean, so_mean = xr.align(
        o2_mean, thetao_mean, so_mean, join='inner'
    )

    lev = thetao_mean[info['depth_coord']]
    lev_m = lev * info['depth_scale']
    rho = in_situ_density(thetao_mean, so_mean, lev_m)
    o2_umol = (o2_mean * 1000.0 / rho).assign_attrs(units='umol kg-1')

    fallback = HIST_ROOT if exp == 'ssp245' else None
    vol = open_cmip_field(
        root, model, 'volcello', exp,
        y0=y0, y1=y1, fallback_root=fallback,
    ).squeeze(drop=True)
    area = open_cmip_field(
        root, model, 'areacello', exp,
        y0=None, y1=None, fallback_root=fallback,
    ).squeeze(drop=True)
    deptho = open_cmip_field(
        root, model, 'deptho', exp,
        y0=None, y1=None, fallback_root=fallback,
    ).squeeze(drop=True)

    o2_umol, o2_mean, vol = xr.align(o2_umol, o2_mean, vol, join='inner')

    area = area.reindex({
        info['y_dim']: o2_mean[info['y_dim']],
        info['x_dim']: o2_mean[info['x_dim']],
    })
    thickness = (vol / area).where(area > 0)

    lat = o2_mean[info['lat_coord']]
    lon = o2_mean[info['lon_coord']]

    return {
        'name': f'{model}_{exp}',
        'model': model,
        'exp': exp,
        'o2_umol': o2_umol,
        'o2_native': o2_mean,
        'vol': vol,
        'area': area,
        'deptho': deptho,
        'thickness': thickness,
        'lat': lat,
        'lon': lon,
        'depth_dim': info['depth_coord'],
        'lat_coord': info['lat_coord'],
        'lon_coord': info['lon_coord'],
        'grid_type': info['grid_type'],
    }


sources = {'WOA23': woa_source}

for model in ['GFDL-ESM4', 'MPI-ESM1-2-LR']:
    sources[f'{model}_historical'] = load_cmip_source(model, 'historical', 1971, 2000)

for model in ['CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']:
    sources[f'{model}_ssp245'] = load_cmip_source(model, 'ssp245', 2071, 2100)

# -----------------------------------------------------------------------------
# 3. OMZ threshold volumes, inventories, thickness maxima
# -----------------------------------------------------------------------------
threshold_volumes = {}
for key, src in sources.items():
    threshold_volumes[key] = {
        'below20_1e15_m3': threshold_volume(src['o2_umol'], src['vol'], THRESH_20),
        'below60_1e15_m3': threshold_volume(src['o2_umol'], src['vol'], THRESH_60),
    }

inventories = {}
for key in [
    'GFDL-ESM4_historical',
    'MPI-ESM1-2-LR_historical',
    'CESM2_ssp245',
    'GFDL-ESM4_ssp245',
    'MPI-ESM1-2-LR_ssp245',
]:
    src = sources[key]
    inventories[key] = oxygen_inventory(src['o2_native'], src['vol'])

omz_thickness_fields = {}
max_thickness = {}
for key, src in sources.items():
    thr = omz_thickness(
        src['o2_umol'], src['thickness'], src['depth_dim'], THRESH_60
    ).compute()
    omz_thickness_fields[key] = thr
    maxv, lat0, lon0 = max_location(thr, src['lat'], src['lon'])
    max_thickness[key] = {
        'max_thickness_m': maxv,
        'latitude': lat0,
        'longitude': lon0,
    }

# -----------------------------------------------------------------------------
# 4. WOA23 OMZ water-mass character: volume-weighted AOU and O2 saturation
# -----------------------------------------------------------------------------
aou_woa, o2sat_woa, o2_woa_a, vol_woa_a = xr.align(
    aou_woa, o2sat_woa, o2_woa, vol_woa, join='inner'
)
omz_mask_woa = (o2_woa_a < THRESH_60) & o2_woa_a.notnull() & vol_woa_a.notnull()

w_aou = vol_woa_a.where(omz_mask_woa & aou_woa.notnull(), 0.0)
aou_mean_woa_omz = float(
    (aou_woa * w_aou).sum().compute() / w_aou.sum().compute()
)

w_sat = vol_woa_a.where(omz_mask_woa & o2sat_woa.notnull(), 0.0)
o2sat_mean_woa_omz = float(
    (o2sat_woa * w_sat).sum().compute() / w_sat.sum().compute()
)

# -----------------------------------------------------------------------------
# 5. Basin OMZ volumes at the 60 umol kg-1 threshold
# -----------------------------------------------------------------------------
def basin_omz_volume(src, basin_labels, thr=THRESH_60):
    out = {}
    for code, name in [(2, 'Pacific'), (3, 'Indian'), (1, 'Atlantic')]:
        mask = (
            (src['o2_umol'] < thr)
            & src['o2_umol'].notnull()
            & src['vol'].notnull()
            & (basin_labels == code)
        )
        amount = src['vol'].where(mask, 0.0).sum().compute()
        out[name] = float(amount) / 1.0e15
    return out


basin_volumes = {}
basin_volumes['WOA23'] = basin_omz_volume(
    woa_source, basin_surf_folded, THRESH_60
)

for model in ['GFDL-ESM4', 'MPI-ESM1-2-LR']:
    key = f'{model}_historical'
    src = sources[key]
    labels_model = map_basin_to_grid(basin_surf_folded, src['lat'], src['lon'])
    basin_volumes[key] = basin_omz_volume(src, labels_model, THRESH_60)

# -----------------------------------------------------------------------------
# 6. Scenario changes for GFDL-ESM4 and MPI-ESM1-2-LR
# -----------------------------------------------------------------------------
scenario_changes = {}
for model in ['GFDL-ESM4', 'MPI-ESM1-2-LR']:
    h = threshold_volumes[f'{model}_historical']
    s = threshold_volumes[f'{model}_ssp245']
    pct20 = 100.0 * (s['below20_1e15_m3'] - h['below20_1e15_m3']) / h['below20_1e15_m3']
    pct60 = 100.0 * (s['below60_1e15_m3'] - h['below60_1e15_m3']) / h['below60_1e15_m3']
    scenario_changes[model] = {
        'below20_percent_change': float(pct20),
        'below60_percent_change': float(pct60),
    }

# -----------------------------------------------------------------------------
# 7. Which historical model is closer to WOA23?
# -----------------------------------------------------------------------------
woa_v20 = threshold_volumes['WOA23']['below20_1e15_m3']
woa_v60 = threshold_volumes['WOA23']['below60_1e15_m3']

model_woa_differences = {}
for model in ['GFDL-ESM4', 'MPI-ESM1-2-LR']:
    key = f'{model}_historical'
    d20 = threshold_volumes[key]['below20_1e15_m3'] - woa_v20
    d60 = threshold_volumes[key]['below60_1e15_m3'] - woa_v60
    model_woa_differences[model] = {
        'diff20_1e15_m3': float(d20),
        'diff60_1e15_m3': float(d60),
        'sum_abs_diff_1e15_m3': float(abs(d20) + abs(d60)),
    }

closer_model = min(
    model_woa_differences,
    key=lambda m: model_woa_differences[m]['sum_abs_diff_1e15_m3'],
)

# -----------------------------------------------------------------------------
# 8. Collect results and print summary
# -----------------------------------------------------------------------------
result = {
    'threshold_volumes_1e15_m3': threshold_volumes,
    'oxygen_inventory_mol': inventories,
    'max_omz_thickness_60': max_thickness,
    'woa_omz_aou_umol_kg': aou_mean_woa_omz,
    'woa_omz_o2sat_percent': o2sat_mean_woa_omz,
    'basin_omz_volumes_60_1e15_m3': basin_volumes,
    'scenario_changes_percent': scenario_changes,
    'model_woa_differences_1e15_m3': model_woa_differences,
    'closer_historical_model_to_WOA23': closer_model,
}

print('=' * 78)
print('Threshold OMZ volumes (1e15 m3)')
print('=' * 78)
for k, v in threshold_volumes.items():
    print(
        f"{k:28s}  <20: {v['below20_1e15_m3']:12.3e}   "
        f"<60: {v['below60_1e15_m3']:12.3e}"
    )

print()
print('=' * 78)
print('Maximum OMZ column thickness at O2 < 60 umol kg-1')
print('=' * 78)
for k, v in max_thickness.items():
    print(
        f"{k:28s}  {v['max_thickness_m']:10.1f} m   "
        f"at ({v['latitude']:7.2f}, {v['longitude']:8.2f})"
    )

print()
print('=' * 78)
print('Dissolved oxygen inventory (mol)')
print('=' * 78)
for k, v in inventories.items():
    print(f"{k:28s}  {v:12.3e}")

print()
print('=' * 78)
print('WOA23 OMZ water-mass character at O2 < 60 umol kg-1')
print('=' * 78)
print(f'Volume-weighted AOU       : {aou_mean_woa_omz:8.1f} umol kg-1')
print(f'Volume-weighted O2 sat    : {o2sat_mean_woa_omz:8.1f} percent')

print()
print('=' * 78)
print('Basin OMZ volumes at O2 < 60 umol kg-1 (1e15 m3)')
print('=' * 78)
for k, v in basin_volumes.items():
    print(
        f"{k:28s}  Pacific: {v['Pacific']:10.3e}   "
        f"Indian: {v['Indian']:10.3e}   Atlantic: {v['Atlantic']:10.3e}"
    )

print()
print('=' * 78)
print('Historical -> ssp245 OMZ volume change (percent, positive = expansion)')
print('=' * 78)
for k, v in scenario_changes.items():
    print(
        f"{k:28s}  <20: {v['below20_percent_change']:8.2f} %   "
        f"<60: {v['below60_percent_change']:8.2f} %"
    )

print()
print('=' * 78)
print('Historical model differences from WOA23 (1e15 m3)')
print('=' * 78)
for k, v in model_woa_differences.items():
    print(
        f"{k:28s}  d20: {v['diff20_1e15_m3']:10.3e}   "
        f"d60: {v['diff60_1e15_m3']:10.3e}   "
        f"|d20|+|d60|: {v['sum_abs_diff_1e15_m3']:10.3e}"
    )
print(f'Closer historical model to WOA23: {closer_model}')
