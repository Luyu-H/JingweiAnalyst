import os, glob
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------
def sigma0(T, S):
    """UNESCO 1983 potential density anomaly (sigma0) from potential
    temperature (deg C) and practical salinity."""
    T = T.astype(np.float32)
    S = S.astype(np.float32)
    rho_w = (999.842594 + 6.793952e-2*T - 9.095290e-3*T**2
             + 1.001685e-4*T**3 - 1.120083e-6*T**4 + 6.536332e-9*T**5)
    A = (8.24493e-1 - 4.0899e-3*T + 7.6438e-5*T**2
         - 8.2467e-7*T**3 + 5.3875e-9*T**4)
    B = (-5.72466e-3 + 1.0227e-4*T - 1.6546e-6*T**2)
    C = 4.8314e-4
    rho = rho_w + A*S + B*S**1.5 + C*S**2
    return (rho - 1000.0).astype(np.float32)


def compute_mld(theta, so, lev):
    """
    Operator MLD: depth where sigma0 exceeds sigma0(10 m) + 0.03 kg/m3.
    theta, so : (nt, nlev, ny, nx)
    lev       : (nlev,) positive down, metres
    Returns   : (nt, ny, nx) metres, NaN where undefined.
    """
    nt, nlev, ny, nx = theta.shape
    mld = np.full((nt, ny, nx), np.nan, dtype=np.float32)

    # Level bracketing the 10 m reference depth.
    if lev[0] > 10.0:
        i10 = 0
        w10 = 0.0
    else:
        i10 = int(np.clip(np.searchsorted(lev, 10.0) - 1, 0, nlev - 2))
        w10 = (10.0 - lev[i10]) / (lev[i10 + 1] - lev[i10])

    y_idx, x_idx = np.indices((ny, nx))

    for t in range(nt):
        T = theta[t]
        S = so[t]
        sig = sigma0(T, S)                 # (nlev, ny, nx)

        # Reference sigma0 at 10 m.
        if lev[0] > 10.0:
            T10 = T[0]
            S10 = S[0]
        else:
            T10 = T[i10] * (1 - w10) + T[i10 + 1] * w10
            S10 = S[i10] * (1 - w10) + S[i10 + 1] * w10
        sig_ref = sigma0(T10, S10)
        target = sig_ref + 0.03

        # First level where sigma0 reaches the criterion.
        exceeds = ((sig >= target[None, :, :])
                   & np.isfinite(sig)
                   & np.isfinite(target)[None, :, :])
        k = np.argmax(exceeds, axis=0)
        found = np.any(exceeds, axis=0)

        k_prev = np.maximum(k - 1, 0)
        sig_k = sig[k, y_idx, x_idx]
        sig_prev = sig[k_prev, y_idx, x_idx]
        d_k = lev[k]
        d_prev = lev[k_prev]

        valid_prev = (k == 0) | np.isfinite(sig_prev)
        denom = sig_k - sig_prev
        valid = (found & np.isfinite(sig_k) & valid_prev
                 & ((k == 0) | (denom > 0)))

        frac = np.where((k == 0) | (denom == 0), 0.0,
                        (target - sig_prev) / denom)
        mld_t = np.where(k == 0, d_k, d_prev + frac * (d_k - d_prev))
        mld_t = np.where(valid, mld_t, np.nan)
        mld[t] = mld_t

    return mld


def get_files(root, model, var, table, grid_label):
    path = os.path.join(root, model, var)
    pattern = os.path.join(
        path, f'{var}_{table}_{model}_historical_*_{grid_label}*.nc'
    )
    files = glob.glob(pattern)
    if not files:
        pattern = os.path.join(
            path, f'{var}_{table}_{model}_historical_*.nc'
        )
        files = glob.glob(pattern)
    if not files:
        pattern = os.path.join(path, f'{var}_*.nc')
        files = glob.glob(pattern)
    return sorted(files)


def load_area(root, model, grid_label, info):
    """Load native areacello and, for CESM2, apply the sftof ocean mask."""
    area_files = get_files(root, model, 'areacello', 'Ofx', grid_label)
    ds_area = xr.open_dataset(area_files[0])
    area_da = ds_area['areacello']

    if info.get('ocean_mask_from') == 'sftof':
        sftof_files = get_files(root, model, 'sftof', 'Ofx', grid_label)
        ds_sftof = xr.open_dataset(sftof_files[0])
        sftof = ds_sftof['sftof']
        area_da = area_da.where(sftof > 0)

    area = area_da.values.astype(np.float32)
    area = np.where(area > 0, area, np.nan)
    total_area = np.nansum(area)

    ds_area.close()
    if info.get('ocean_mask_from') == 'sftof':
        ds_sftof.close()

    return area, total_area


def get_latlon_2d(lat, lon):
    if lat.ndim == 1 and lon.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon, lat)
    else:
        lat2d, lon2d = lat, lon
    return lat2d, lon2d


def regrid_to_regular(lat2d, lon2d, field, dlat=2.0, dlon=2.0):
    """
    Simple binning to a regular lon/lat grid for display only.
    Curvilinear fields are averaged into dlat x dlon boxes.
    """
    lat_edges = np.arange(-90, 90 + dlat, dlat)
    lon_edges = np.arange(-180, 180 + dlon, dlon)
    lat_centers = (lat_edges[:-1] + lat_edges[1:]) / 2
    lon_centers = (lon_edges[:-1] + lon_edges[1:]) / 2

    lon2d = np.where(lon2d > 180, lon2d - 360, lon2d)

    iy = np.digitize(lat2d.ravel(), lat_edges) - 1
    ix = np.digitize(lon2d.ravel(), lon_edges) - 1
    valid = ((iy >= 0) & (iy < len(lat_centers))
             & (ix >= 0) & (ix < len(lon_centers))
             & np.isfinite(field.ravel()))

    sum_arr = np.zeros((len(lat_centers), len(lon_centers)), dtype=np.float32)
    cnt = np.zeros_like(sum_arr)
    np.add.at(sum_arr, (iy[valid], ix[valid]),
              field.ravel()[valid].astype(np.float32))
    np.add.at(cnt, (iy[valid], ix[valid]), 1)
    grid = np.where(cnt > 0, sum_arr / cnt, np.nan)
    return lon_centers, lat_centers, grid


# ----------------------------------------------------------------------
# Main computation
# ----------------------------------------------------------------------
root = '/data1/lorraine-20419/benchmark_cmip6/historical'
models = ['CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']

model_info = {
    'CESM2': {
        'grid_label': 'gn', 'y_dim': 'nlat', 'x_dim': 'nlon',
        'lat_coord': 'lat', 'lon_coord': 'lon', 'depth_coord': 'lev',
        'depth_units': 'centimeters', 'depth_scale': 0.01,
        'ocean_mask_from': 'sftof'
    },
    'GFDL-ESM4': {
        'grid_label': 'gr', 'y_dim': 'lat', 'x_dim': 'lon',
        'lat_coord': 'lat', 'lon_coord': 'lon', 'depth_coord': 'lev',
        'depth_units': 'm', 'depth_scale': 1.0
    },
    'MPI-ESM1-2-LR': {
        'grid_label': 'gn', 'y_dim': 'j', 'x_dim': 'i',
        'lat_coord': 'latitude', 'lon_coord': 'longitude',
        'depth_coord': 'lev', 'depth_units': 'm', 'depth_scale': 1.0
    }
}

results = {}

for model in models:
    info = model_info[model]
    grid_label = info['grid_label']

    theta_files = get_files(root, model, 'thetao', 'Omon', grid_label)
    so_files = get_files(root, model, 'so', 'Omon', grid_label)
    mlotst_files = get_files(root, model, 'mlotst', 'Omon', grid_label)

    ds_t = xr.open_mfdataset(theta_files, use_cftime=True, combine='by_coords')
    ds_s = xr.open_mfdataset(so_files, use_cftime=True, combine='by_coords')
    ds_m = xr.open_mfdataset(mlotst_files, use_cftime=True, combine='by_coords')

    ds_t = ds_t.sel(time=slice('1971-01-01', '2000-12-31'))
    ds_s = ds_s.sel(time=slice('1971-01-01', '2000-12-31'))
    ds_m = ds_m.sel(time=slice('1971-01-01', '2000-12-31'))

    area, total_area = load_area(root, model, grid_label, info)

    for month in [3, 9]:
        idx_t = np.where((ds_t.time.dt.month == month).values)[0]
        idx_s = np.where((ds_s.time.dt.month == month).values)[0]
        idx_m = np.where((ds_m.time.dt.month == month).values)[0]

        theta_da = ds_t.thetao.isel(time=idx_t).astype(np.float32)
        so_da = ds_s.so.isel(time=idx_s).astype(np.float32)
        ml_da = ds_m.mlotst.isel(time=idx_m).astype(np.float32)

        theta_da, so_da, ml_da = xr.align(theta_da, so_da, ml_da,
                                          join='inner')

        theta = theta_da.load().values
        so = so_da.load().values
        ml = ml_da.load().values

        lev = theta_da.lev.values.astype(np.float32)
        units = theta_da.lev.attrs.get('units', '').lower()
        if units in ['cm', 'centimeters', 'centimeter']:
            lev = lev * 0.01

        lat = theta_da[info['lat_coord']].values
        lon = theta_da[info['lon_coord']].values

        # Operator MLD from the model's own thetao / so.
        op = compute_mld(theta, so, lev)

        # Cells where both fields are defined in all 30 instances.
        valid = np.all(np.isfinite(op) & np.isfinite(ml), axis=0)

        op_mean = np.nanmean(op, axis=0)
        ml_mean = np.nanmean(ml, axis=0)
        op_mean[~valid] = np.nan
        ml_mean[~valid] = np.nan

        diff = op_mean - ml_mean

        mask = (valid & np.isfinite(diff)
                & np.isfinite(area) & (area > 0))

        a = area[mask]
        d = diff[mask]
        w = a / a.sum()

        # Area-weighted statistics.
        bias = np.sum(w * d)
        rmse = np.sqrt(np.sum(w * d**2))

        op_v = op_mean[mask]
        ml_v = ml_mean[mask]
        mean_op = np.sum(w * op_v)
        mean_ml = np.sum(w * ml_v)
        cov = np.sum(w * (op_v - mean_op) * (ml_v - mean_ml))
        var_op = np.sum(w * (op_v - mean_op)**2)
        var_ml = np.sum(w * (ml_v - mean_ml)**2)
        corr = (cov / np.sqrt(var_op * var_ml)
                if var_op > 0 and var_ml > 0 else np.nan)

        pct_gt20 = 100 * np.sum(a[np.abs(d) > 20]) / np.sum(a)
        pct_area = 100 * np.sum(a) / total_area

        results[(model, month)] = dict(
            bias=bias, rmse=rmse, corr=corr,
            pct_gt20=pct_gt20, pct_area=pct_area,
            op_mean=op_mean, ml_mean=ml_mean, diff=diff,
            mask=mask, lat=lat, lon=lon, area=area
        )

        print(f'{model} {month:02d}: bias={bias:.1f} m, '
              f'RMSE={rmse:.1f} m, r={corr:.3f}, '
              f'>20m={pct_gt20:.1f}%, valid area={pct_area:.1f}%')

    ds_t.close()
    ds_s.close()
    ds_m.close()

# Ranking criterion: lowest area-weighted RMSE averaged over March and
# September.  The operator and mlotst are not expected to agree exactly
# because the model's mlotst may use a different reference depth, a
# different density criterion/algorithm, and because the vertical
# interpolation and remapping differ.  Regions of large disagreement are
# visible in the maps (Fig. 1), typically high-latitude deep-convection
# and marginal-sea regions where the mixed layer is deep and bracketing
# levels are sparse.
model_rmse = {
    m: np.mean([results[(m, 3)]['rmse'], results[(m, 9)]['rmse']])
    for m in models
}
ranking = sorted(models, key=lambda m: model_rmse[m])

print('\nRanking (most to least internally consistent, '
      'by mean RMSE of March and September):')
for m in ranking:
    print(f'  {m}: mean RMSE = {model_rmse[m]:.1f} m')
