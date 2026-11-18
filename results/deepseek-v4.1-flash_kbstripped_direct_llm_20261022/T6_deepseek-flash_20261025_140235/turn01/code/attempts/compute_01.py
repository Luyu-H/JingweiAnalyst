import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.stats import binned_statistic

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
def read_woa_msk(path):
    """Read a WOA .msk ASCII table, skipping header/comment lines."""
    rows = []
    with open(path, 'r', errors='ignore') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            parts = s.split()
            if len(parts) < 2:
                continue
            try:
                lat = float(parts[0])
                lon = float(parts[1])
            except ValueError:
                continue
            vals = []
            for p in parts[2:]:
                try:
                    vals.append(float(p))
                except ValueError:
                    vals.append(np.nan)
            rows.append((lat, lon, vals))
    return rows


def compute_nitracline(no3, depth):
    """
    Nitracline = depth of maximum vertical nitrate gradient in the top 1000 m.
    The gradient at the deepest level <= 1000 m uses the level just below.
    no3 : array with last dimension = depth
    depth : 1-D array of depths (positive down)
    Returns array of nitracline depths with shape no3.shape[:-1].
    """
    no3 = np.asarray(no3, dtype=float)
    depth = np.asarray(depth, dtype=float)
    if depth[0] > depth[-1]:
        depth = depth[::-1]
        no3 = no3[..., ::-1]

    dno3 = np.diff(no3, axis=-1)
    dz = np.diff(depth)
    grad = dno3 / dz
    mid_depth = 0.5 * (depth[1:] + depth[:-1])

    # only search gradients whose upper level is within the top 1000 m
    valid = depth[:-1] <= 1000.0
    grad_valid = grad[..., valid]
    mid_valid = mid_depth[valid]

    if grad_valid.shape[-1] == 0:
        return np.full(no3.shape[:-1], np.nan)

    grad_filled = np.where(np.isfinite(grad_valid), grad_valid, -np.inf)
    idx = np.argmax(grad_filled, axis=-1)
    out = np.take(mid_valid, idx)
    all_invalid = np.all(grad_filled == -np.inf, axis=-1)
    out = np.where(all_invalid, np.nan, out)
    return out


def regrid_to_woa(model_lat2d, model_lon2d, model_data, woa_lat, woa_lon, max_dist=2.5):
    """
    Nearest-neighbour regridding from a (possibly curvilinear) model grid to the
    WOA 1-degree grid.  Only finite model values are used.
    """
    valid = np.isfinite(model_data)
    if not np.any(valid):
        return np.full((len(woa_lat), len(woa_lon)), np.nan)

    pts = np.column_stack([model_lat2d[valid], (model_lon2d[valid] % 360)])
    vals = model_data[valid]

    # periodic copies in longitude
    pts_ext = np.vstack([pts, pts + [0, 360], pts - [0, 360]])
    vals_ext = np.tile(vals, 3)

    tree = cKDTree(pts_ext)
    woa_lon360 = woa_lon % 360
    woa_lat2d, woa_lon2d = np.meshgrid(woa_lat, woa_lon360, indexing='ij')
    query_pts = np.column_stack([woa_lat2d.ravel(), woa_lon2d.ravel()])
    dist, idx = tree.query(query_pts, k=1)

    out = vals_ext[idx].reshape(len(woa_lat), len(woa_lon)).astype(float)
    out[dist.reshape(len(woa_lat), len(woa_lon)) > max_dist] = np.nan
    return out


def area_weighted_mean(data, mask, lat):
    """Area-weighted mean over a 2-D field on a regular lat-lon grid."""
    w = np.cos(np.deg2rad(lat))[:, None] * np.ones((1, data.shape[1]))
    valid = mask & np.isfinite(data)
    if not np.any(valid):
        return np.nan
    return np.sum(data[valid] * w[valid]) / np.sum(w[valid])


def compute_section(no3_np, lat2d, lon2d, depth, lon_bins, lat_min=-2, lat_max=2):
    """
    Build a longitude-depth section averaged over a latitude band.
    no3_np : (ny, nx, nlev)
    Returns (nlev, nlon_bins)
    """
    nlev = len(depth)
    nbin = len(lon_bins) - 1
    section = np.full((nlev, nbin), np.nan)
    lats = lat2d.ravel()
    lons = lon2d.ravel() % 360
    lat_mask = (lats >= lat_min) & (lats <= lat_max)

    for k in range(nlev):
        vals = no3_np[..., k].ravel()
        mask = lat_mask & np.isfinite(vals)
        if np.any(mask):
            stat, _, _ = binned_statistic(
                lons[mask], vals[mask], statistic='mean', bins=lon_bins
            )
            section[k, :] = stat
    return section


# ---------------------------------------------------------------------------
# WOA23 nitrate and masks
# ---------------------------------------------------------------------------
woa_root = '/data1/lorraine-20419/ocean_data/WOA23'
woa_nc = os.path.join(woa_root, 'climatology/nitrate/all_time_periods/woa23_all_n00_01.nc')
woa_ds = xr.open_dataset(woa_nc, decode_times=False)

# locate the objective-analysis nitrate field
woa_var = None
for v in woa_ds.data_vars:
    if 'n' in v.lower() and 'an' in v.lower():
        woa_var = v
        break
if woa_var is None:
    for v in woa_ds.data_vars:
        if woa_ds[v].ndim == 4:
            woa_var = v
            break
if woa_var is None:
    woa_var = list(woa_ds.data_vars)[0]

woa_no3_da = woa_ds[woa_var]
if 'time' in woa_no3_da.dims:
    woa_no3_da = woa_no3_da.isel(time=0)
woa_no3_da = woa_no3_da.squeeze()

woa_depth = woa_ds['depth'].values.astype(float)
woa_lat = woa_ds['lat'].values.astype(float)
woa_lon = woa_ds['lon'].values.astype(float)

# (lat, lon, depth)
woa_no3 = woa_no3_da.transpose('lat', 'lon', 'depth').values.astype(float)
nitracline_woa = compute_nitracline(woa_no3, woa_depth)

# land-sea mask
landsea_path = os.path.join(woa_root, 'masks/landsea_01.msk')
landsea_rows = read_woa_msk(landsea_path)
landsea = np.full((180, 360), np.nan)
for lat, lon, vals in landsea_rows:
    i = int(round(lat + 89.5))
    j = int(round(lon + 179.5))
    if 0 <= i < 180 and 0 <= j < 360:
        landsea[i, j] = vals[0]
ocean_mask = landsea > 1

# basin mask
basin_path = os.path.join(woa_root, 'masks/basinmask_01.msk')
basin_rows = read_woa_msk(basin_path)
nlev_basin = 102
basin3d = np.full((nlev_basin, 180, 360), np.nan)
for lat, lon, vals in basin_rows:
    i = int(round(lat + 89.5))
    j = int(round(lon + 179.5))
    if 0 <= i < 180 and 0 <= j < 360:
        for k, v in enumerate(vals):
            if k < nlev_basin:
                basin3d[k, i, j] = v

# surface basin code (first valid level)
basin_map = np.full((180, 360), np.nan)
has_valid = np.any(np.isfinite(basin3d), axis=0)
first_valid = np.argmax(np.isfinite(basin3d), axis=0)
for i in range(180):
    for j in range(360):
        if has_valid[i, j]:
            basin_map[i, j] = basin3d[first_valid[i, j], i, j]

# Fold marginal seas into parent oceans.
# WOA13 codes: 1 Atlantic, 2 Pacific, 3 Indian, 4 Med, 5 Baltic, 6 Black Sea,
# 7 Red Sea, 8 Persian Gulf, 9 Hudson Bay, 10 Southern, 11 Arctic.
basin_group = np.zeros((180, 360), dtype=int)
for code in [1, 4, 5, 6, 9]:
    basin_group[basin_map == code] = 1      # Atlantic
basin_group[basin_map == 2] = 2             # Pacific
for code in [3, 7, 8]:
    basin_group[basin_map == code] = 3      # Indian
# 10 and 11 remain 0 (excluded)

# ---------------------------------------------------------------------------
# CMIP6 historical no3, 1971-2000
# ---------------------------------------------------------------------------
cmip_root = '/data1/lorraine-20419/benchmark_cmip6/historical'
models = ['CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']
model_meta = {
    'CESM2':        {'y_dim': 'nlat', 'x_dim': 'nlon', 'depth_coord': 'lev', 'depth_scale': 0.01},
    'GFDL-ESM4':    {'y_dim': 'lat',  'x_dim': 'lon',  'depth_coord': 'lev', 'depth_scale': 1.0},
    'MPI-ESM1-2-LR':{'y_dim': 'j',    'x_dim': 'i',    'depth_coord': 'lev', 'depth_scale': 1.0},
}

nitracline_models = {}
model_no3_np = {}
model_lats = {}
model_lons = {}
model_depths = {}

for model in models:
    model_dir = os.path.join(cmip_root, model, 'no3')
    files = sorted(glob.glob(os.path.join(model_dir, '*.nc')))
    if not files:
        raise FileNotFoundError(f'No no3 files for {model} in {model_dir}')

    # keep only the ensemble member that actually appears in the filenames
    base = os.path.basename(files[0])
    parts = base.split('_')
    variant = parts[4] if len(parts) > 4 else None
    if variant:
        files = [f for f in files if variant in os.path.basename(f)]

    ds = xr.open_mfdataset(
        files, combine='by_coords', use_cftime=True, decode_times=True,
        chunks={'time': 12}
    )
    years = ds.time.dt.year
    time_mask = (years >= 1971) & (years <= 2000)
    ds_sel = ds.isel(time=time_mask)

    no3_da = ds_sel['no3'].mean('time').compute()

    lev = ds_sel[model_meta[model]['depth_coord']].values.astype(float)
    depth = lev * model_meta[model]['depth_scale']
    model_depths[model] = depth

    y_dim = model_meta[model]['y_dim']
    x_dim = model_meta[model]['x_dim']
    lat_da = ds_sel[y_dim]
    lon_da = ds_sel[x_dim]

    if lat_da.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon_da.values, lat_da.values)
    else:
        lat2d = lat_da.values
        lon2d = lon_da.values

    # (ny, nx, nlev)
    no3_np = no3_da.transpose(y_dim, x_dim, 'lev').values.astype(float)
    model_no3_np[model] = no3_np
    model_lats[model] = lat2d
    model_lons[model] = lon2d

    nitracline_native = compute_nitracline(no3_np, depth)
    nitracline_models[model] = regrid_to_woa(
        lat2d, lon2d, nitracline_native, woa_lat, woa_lon
    )

    ds.close()

# ---------------------------------------------------------------------------
# Global and basin statistics
# ---------------------------------------------------------------------------
basin_names = {1: 'Atlantic', 2: 'Pacific', 3: 'Indian'}
results = {}

# WOA
results['WOA23'] = {'global': area_weighted_mean(nitracline_woa, ocean_mask, woa_lat)}
basin_means_woa = {}
for bg, name in basin_names.items():
    mask = (basin_group == bg) & ocean_mask
    basin_means_woa[name] = area_weighted_mean(nitracline_woa, mask, woa_lat)
results['WOA23']['basins'] = basin_means_woa

# models
for model in models:
    results[model] = {
        'global': area_weighted_mean(nitracline_models[model], ocean_mask, woa_lat)
    }
    basin_means = {}
    for bg, name in basin_names.items():
        mask = (basin_group == bg) & ocean_mask
        basin_means[name] = area_weighted_mean(nitracline_models[model], mask, woa_lat)
    results[model]['basins'] = basin_means

# differences from WOA
for model in models:
    results[model]['global_diff'] = results[model]['global'] - results['WOA23']['global']
    for name in basin_names.values():
        results[model][f'{name}_diff'] = (
            results[model]['basins'][name] - results['WOA23']['basins'][name]
        )

# spatial correlation and RMSE against WOA on the WOA 1-degree grid
corr_rmse = {}
for model in models:
    valid = ocean_mask & np.isfinite(nitracline_woa) & np.isfinite(nitracline_models[model])
    x = nitracline_woa[valid]
    y = nitracline_models[model][valid]
    if x.size > 1:
        corr = float(np.corrcoef(x, y)[0, 1])
        rmse = float(np.sqrt(np.mean((x - y) ** 2)))
    else:
        corr, rmse = np.nan, np.nan
    corr_rmse[model] = {'corr': corr, 'rmse': rmse}

# observed and modelled basin ordering (deepest first)
woa_basin_order = sorted(basin_names.values(), key=lambda n: results['WOA23']['basins'][n], reverse=True)
model_basin_orders = {}
for model in models:
    model_basin_orders[model] = sorted(
        basin_names.values(), key=lambda n: results[model]['basins'][n], reverse=True
    )

# ---------------------------------------------------------------------------
# Longitude-depth sections through the equatorial Pacific (2S-2N)
# ---------------------------------------------------------------------------
lon_bins = np.arange(120, 281, 1)          # 120E to 80W
lon_centers = 0.5 * (lon_bins[:-1] + lon_bins[1:])

sections = {}
# WOA
woa_lon2d, woa_lat2d = np.meshgrid(woa_lon, woa_lat)
sections['WOA23'] = compute_section(woa_no3, woa_lat2d, woa_lon2d, woa_depth, lon_bins)

# CMIP6
for model in models:
    sections[model] = compute_section(
        model_no3_np[model], model_lats[model], model_lons[model],
        model_depths[model], lon_bins
    )

# Convert CMIP6 mol m-3 to umol kg-1 (reference density 1025 kg m-3).
# WOA is already in umol kg-1.
sections_umol = {'WOA23': sections['WOA23']}
for model in models:
    sections_umol[model] = sections[model] * 1.0e6 / 1025.0

# Nitracline line for each section (computed from the section's own profile)
nitracline_lines = {}
for key in ['WOA23'] + models:
    depth = woa_depth if key == 'WOA23' else model_depths[key]
    nitracline_lines[key] = compute_nitracline(sections[key].T, depth)

# ---------------------------------------------------------------------------
# Report-ready rounded values
# ---------------------------------------------------------------------------
summary = {}
for key in ['WOA23'] + models:
    summary[key] = {
        'global_m': round(float(results[key]['global']), 0),
        'basins_m': {n: round(float(results[key]['basins'][n]), 0) for n in basin_names.values()},
    }
    if key != 'WOA23':
        summary[key]['global_diff_m'] = round(float(results[key]['global_diff']), 0)
        summary[key]['basin_diff_m'] = {
            n: round(float(results[key][f'{n}_diff']), 0) for n in basin_names.values()
        }
        summary[key]['corr'] = round(corr_rmse[key]['corr'], 3)
        summary[key]['rmse_m'] = round(corr_rmse[key]['rmse'], 1)
