import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# Paths and metadata
# -----------------------------------------------------------------------------
cmip6_root = "/data1/lorraine-20419/benchmark_cmip6/historical"
woa_root = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
woa_masks_root = "/data1/lorraine-20419/ocean_data/WOA23/masks"

models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

models_info = {
    "CESM2": {
        "lat_coord": "lat",
        "lon_coord": "lon",
        "depth_coord": "lev",
        "depth_scale": 0.01,
    },
    "GFDL-ESM4": {
        "lat_coord": "lat",
        "lon_coord": "lon",
        "depth_coord": "lev",
        "depth_scale": 1.0,
    },
    "MPI-ESM1-2-LR": {
        "lat_coord": "latitude",
        "lon_coord": "longitude",
        "depth_coord": "lev",
        "depth_scale": 1.0,
    },
}

basin_codes = {"Atlantic": 1, "Pacific": 2, "Indian": 3}
basin_names = ["Pacific", "Atlantic", "Indian"]

# -----------------------------------------------------------------------------
# Generic helpers
# -----------------------------------------------------------------------------
def pick_var(ds, var):
    if var in ds.variables:
        return ds[var]
    for name in ds.data_vars:
        if name.lower() == str(var).lower():
            return ds[name]
    for name in ds.data_vars:
        if "bnd" not in name.lower() and ds[name].ndim >= 2:
            return ds[name]
    raise KeyError(f"Could not find variable {var!r} in {list(ds.data_vars)}")


def find_cmip6_files(model, var):
    d = os.path.join(cmip6_root, model, var)
    files = sorted(glob.glob(os.path.join(d, f"{var}_*.nc")))
    if not files:
        files = sorted(glob.glob(os.path.join(d, "*.nc")))
    files = [
        f for f in files
        if model in os.path.basename(f) and "historical" in os.path.basename(f)
    ]
    if not files:
        files = sorted(glob.glob(os.path.join(d, "*.nc")))
    if not files:
        raise FileNotFoundError(f"No files for {model}/{var} under {d}")
    return files


def open_cmip6_var(model, var):
    files = find_cmip6_files(model, var)
    das = []
    for f in files:
        ds = xr.open_dataset(f, use_cftime=True, decode_times=True)
        das.append(pick_var(ds, var))
    if len(das) == 1:
        return das[0]
    if all("time" in da.dims for da in das):
        da = xr.concat(das, dim="time")
        da = da.sortby("time")
        return da
    return das[0]


def get_cmip6_climatology(model, var, y0=1971, y1=2000):
    da = open_cmip6_var(model, var)
    if "time" not in da.dims:
        return da.squeeze(drop=True)
    years = da.time.dt.year
    tmask = (years >= y0) & (years <= y1)
    da = da.isel(time=np.where(tmask.values)[0])
    w = da.time.dt.days_in_month
    w_b = w.broadcast_like(da)
    valid = da.notnull()
    num = (da.fillna(0.0) * w_b).sum("time")
    den = (valid.astype(float) * w_b).sum("time")
    clim = (num / den).where(den > 0)
    return clim.compute()


def read_msk_first3(path):
    rows = []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if (not s) or s.startswith("#"):
                continue
            parts = s.replace(",", " ").split()
            vals = []
            for p in parts[:3]:
                try:
                    vals.append(float(p))
                except ValueError:
                    break
            if len(vals) == 3:
                rows.append(vals)
    return np.asarray(rows, dtype=float)


def clean_latlon_rows(arr):
    if arr.size == 0:
        raise ValueError("No numeric rows found")
    a, b, c = arr[:, 0], arr[:, 1], arr[:, 2]
    ok1 = (np.abs(a) <= 90.1) & (np.abs(b) <= 180.1) & np.isfinite(c)
    ok2 = (np.abs(b) <= 90.1) & (np.abs(a) <= 180.1) & np.isfinite(c)
    if ok1.sum() >= ok2.sum():
        return a[ok1], b[ok1], c[ok1]
    return b[ok2], a[ok2], c[ok2]


def table_to_grid(lat_vals, lon_vals, vals, lat_coord, lon_coord):
    nlat = len(lat_coord)
    nlon = len(lon_coord)
    grid = np.full((nlat, nlon), np.nan, dtype=float)
    lat_idx = np.abs(lat_vals[:, None] - lat_coord[None, :]).argmin(axis=1)
    lon_idx = np.abs(lon_vals[:, None] - lon_coord[None, :]).argmin(axis=1)
    grid[lat_idx, lon_idx] = vals
    return grid


def fill_code_nearest(code_grid, max_iter=300):
    f = code_grid.copy()
    for _ in range(max_iter):
        nan = np.isnan(f)
        if not np.any(nan):
            break
        neighs = []
        up = np.full_like(f, np.nan)
        up[1:, :] = f[:-1, :]
        neighs.append(up)
        down = np.full_like(f, np.nan)
        down[:-1, :] = f[1:, :]
        neighs.append(down)
        left = np.roll(f, 1, axis=1)
        neighs.append(left)
        right = np.roll(f, -1, axis=1)
        neighs.append(right)
        for neigh in neighs:
            take = nan & np.isfinite(neigh)
            f[take] = neigh[take]
            nan = np.isnan(f)
    return f


def compute_cell_areas(lat, lon, R=6371000.0):
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    lat_edges = np.empty(len(lat) + 1, dtype=float)
    lat_edges[1:-1] = 0.5 * (lat[:-1] + lat[1:])
    lat_edges[0] = lat[0] - 0.5 * (lat[1] - lat[0])
    lat_edges[-1] = lat[-1] + 0.5 * (lat[-1] - lat[-2])
    lon_edges = np.empty(len(lon) + 1, dtype=float)
    lon_edges[1:-1] = 0.5 * (lon[:-1] + lon[1:])
    lon_edges[0] = lon[0] - 0.5 * (lon[1] - lon[0])
    lon_edges[-1] = lon[-1] + 0.5 * (lon[-1] - lon[-2])
    dlon = np.deg2rad(np.diff(lon_edges))
    lat_rad = np.deg2rad(lat_edges)
    dsin = np.sin(lat_rad[1:]) - np.sin(lat_rad[:-1])
    return (R ** 2) * dsin[:, None] * dlon[None, :]


def get_latlon2d(da, info):
    lat = da[info["lat_coord"]]
    lon = da[info["lon_coord"]]
    if lat.ndim == 1 and lon.ndim == 1:
        y_dim = lat.dims[0]
        x_dim = lon.dims[0]
        lon2d_np, lat2d_np = np.meshgrid(lon.values, lat.values)
        coords = {y_dim: lat, x_dim: lon}
        lat2d = xr.DataArray(lat2d_np, dims=(y_dim, x_dim), coords=coords)
        lon2d = xr.DataArray(lon2d_np, dims=(y_dim, x_dim), coords=coords)
        return lat2d, lon2d
    lat2d, lon2d = xr.broadcast(lat, lon)
    return lat2d, lon2d


def nearest_indices_flat(lat_vals, lon_vals, lat_target, lon_target, chunk=20000):
    shape = np.asarray(lat_vals).shape
    lat_flat = np.asarray(lat_vals, dtype=float).ravel()
    lon_flat = np.asarray(lon_vals, dtype=float).ravel()
    n = lat_flat.size
    lat_idx = np.empty(n, dtype=int)
    lon_idx = np.empty(n, dtype=int)
    for start in range(0, n, chunk):
        end = min(start + chunk, n)
        lat_idx[start:end] = np.abs(
            lat_flat[start:end, None] - lat_target[None, :]
        ).argmin(axis=1)
        lon_idx[start:end] = np.abs(
            lon_flat[start:end, None] - lon_target[None, :]
        ).argmin(axis=1)
    return lat_idx.reshape(shape), lon_idx.reshape(shape)


def nitracline_from_field(field, depth_values, depth_dim, top=1000.0):
    z = np.abs(np.asarray(depth_values, dtype=float))
    order = np.argsort(z)
    if not np.array_equal(order, np.arange(len(z))):
        field = field.isel({depth_dim: order})
        z = z[order]

    diff_field = field.diff(depth_dim).rename({depth_dim: "lev_mid"})
    dz = np.diff(z)
    grad = diff_field / xr.DataArray(dz, dims=["lev_mid"])
    mid = xr.DataArray(0.5 * (z[:-1] + z[1:]), dims=["lev_mid"])
    shallow = xr.DataArray(z[:-1], dims=["lev_mid"])
    grad = grad.assign_coords(lev_mid=mid)
    grad = grad.where(shallow <= top)

    idx = grad.fillna(-np.inf).argmax("lev_mid")
    nc = mid.isel(lev_mid=idx)
    nc = nc.where(grad.notnull().any("lev_mid"))
    return nc, grad


def area_mean(da, area, mask=None):
    if mask is not None:
        da = da.where(mask)
        area = area.where(mask)
    valid = da.notnull() & area.notnull() & (area > 0)
    num = (da * area).where(valid).sum()
    den = area.where(valid).sum()
    return float((num / den).compute())


def fill_nearest_2d(field, allowed_mask, max_iter=80):
    filled = field.copy()
    filled[~allowed_mask] = np.nan
    for _ in range(max_iter):
        nan = np.isnan(filled) & allowed_mask
        if not np.any(nan):
            break
        sum_neigh = np.zeros_like(filled)
        count = np.zeros_like(filled)

        neigh = np.full_like(filled, np.nan)
        neigh[1:, :] = filled[:-1, :]
        valid = np.isfinite(neigh)
        sum_neigh[valid] += neigh[valid]
        count[valid] += 1

        neigh = np.full_like(filled, np.nan)
        neigh[:-1, :] = filled[1:, :]
        valid = np.isfinite(neigh)
        sum_neigh[valid] += neigh[valid]
        count[valid] += 1

        neigh = np.roll(filled, 1, axis=1)
        valid = np.isfinite(neigh)
        sum_neigh[valid] += neigh[valid]
        count[valid] += 1

        neigh = np.roll(filled, -1, axis=1)
        valid = np.isfinite(neigh)
        sum_neigh[valid] += neigh[valid]
        count[valid] += 1

        take = nan & (count > 0)
        filled[take] = sum_neigh[take] / count[take]
    return filled


def regrid_conservative(field, area, lat_idx, lon_idx, allowed_mask):
    vals = np.asarray(field.values, dtype=float)
    areas = np.asarray(area.values, dtype=float)
    valid = np.isfinite(vals) & np.isfinite(areas) & (areas > 0)
    nlat, nlon = allowed_mask.shape
    ids = lat_idx * nlon + lon_idx
    num = np.bincount(
        ids[valid],
        weights=(vals[valid] * areas[valid]),
        minlength=nlat * nlon,
    )
    den = np.bincount(
        ids[valid],
        weights=areas[valid],
        minlength=nlat * nlon,
    )
    out = np.full(nlat * nlon, np.nan, dtype=float)
    good = den > 0
    out[good] = num[good] / den[good]
    out = out.reshape(nlat, nlon)
    out = fill_nearest_2d(out, allowed_mask)
    out[~allowed_mask] = np.nan
    return out


def weighted_corr_np(a, b, w, mask):
    a = np.asarray(a, dtype=float)[mask]
    b = np.asarray(b, dtype=float)[mask]
    w = np.asarray(w, dtype=float)[mask]
    valid = np.isfinite(a) & np.isfinite(b) & np.isfinite(w)
    a = a[valid]
    b = b[valid]
    w = w[valid]
    w = w / w.sum()
    ma = np.sum(w * a)
    mb = np.sum(w * b)
    cov = np.sum(w * (a - ma) * (b - mb))
    va = np.sum(w * (a - ma) ** 2)
    vb = np.sum(w * (b - mb) ** 2)
    return float(cov / np.sqrt(va * vb))


def weighted_rmse_np(a, b, w, mask):
    a = np.asarray(a, dtype=float)[mask]
    b = np.asarray(b, dtype=float)[mask]
    w = np.asarray(w, dtype=float)[mask]
    valid = np.isfinite(a) & np.isfinite(b) & np.isfinite(w)
    a = a[valid]
    b = b[valid]
    w = w[valid]
    return float(np.sqrt(np.sum(w * (a - b) ** 2) / np.sum(w)))


def simple_density(theta, salt, depth):
    # Lightweight in-situ density estimate in kg m-3, sufficient for the
    # volume-to-mass conversion in Figure 2.  theta in degC, salt ~35.
    rho = 1027.0 - 0.15 * theta + 0.78 * (salt - 35.0) + 0.0045 * depth
    return rho.clip(min=1000.0, max=1050.0)


def section_band(da, area, lat2d, lon2d_180, z, depth_dim,
                 lat_min=-2.0, lat_max=2.0,
                 lon_min=120.0, lon_max=280.0, dlon=1.0):
    ydim, xdim = lat2d.dims
    da_t = da.transpose(depth_dim, ydim, xdim)
    area_t = area.transpose(ydim, xdim)
    latv = lat2d.transpose(ydim, xdim).values
    lonv = ((lon2d_180 + 360.0) % 360.0).transpose(ydim, xdim).values
    areav = area_t.values
    dav = da_t.values

    sel = (
        (latv >= lat_min) & (latv <= lat_max) &
        (lonv >= lon_min) & (lonv <= lon_max) &
        np.isfinite(areav) & (areav > 0)
    )

    lon_edges = np.arange(lon_min - 0.5, lon_max + 0.5 + 1e-9, dlon)
    centers = 0.5 * (lon_edges[:-1] + lon_edges[1:])
    nlonb = len(centers)
    bin_idx = np.floor((lonv - (lon_min - 0.5)) / dlon).astype(int)
    bin_idx = np.clip(bin_idx, 0, nlonb - 1)

    idx_sel = bin_idx[sel]
    area_sel = areav[sel]
    out = np.full((dav.shape[0], nlonb), np.nan, dtype=float)

    for k in range(dav.shape[0]):
        vals = dav[k][sel]
        valid = np.isfinite(vals) & np.isfinite(area_sel)
        if not np.any(valid):
            continue
        idx = idx_sel[valid]
        w = area_sel[valid]
        num = np.bincount(idx, weights=vals[valid] * w, minlength=nlonb)
        den = np.bincount(idx, weights=w, minlength=nlonb)
        good = den > 0
        out[k, good] = num[good] / den[good]

    return {
        "data": out,
        "depth": np.asarray(z, dtype=float),
        "lon": centers,
        "lon_edges": lon_edges,
    }


def nitracline_section(data, depth, top=1000.0):
    z = np.abs(np.asarray(depth, dtype=float))
    order = np.argsort(z)
    if not np.array_equal(order, np.arange(len(z))):
        data = data[order, :]
        z = z[order]
    dz = np.diff(z)
    grad = np.diff(data, axis=0) / dz[:, None]
    mid = 0.5 * (z[:-1] + z[1:])
    allowed = z[:-1] <= top
    grad = np.where(allowed[:, None], grad, np.nan)
    nc = np.full(data.shape[1], np.nan, dtype=float)
    for j in range(data.shape[1]):
        col = grad[:, j]
        if np.any(np.isfinite(col)):
            nc[j] = mid[np.nanargmax(col)]
    return nc

# -----------------------------------------------------------------------------
# WOA23 nitrate, masks, basin labels and areas
# -----------------------------------------------------------------------------
woa_path = os.path.join(
    woa_root, "nitrate", "all_time_periods", "woa23_all_n00_01.nc"
)
woa_ds = xr.open_dataset(woa_path, decode_times=False, mask_and_scale=True)

woa = None
for name in woa_ds.data_vars:
    ln = name.lower()
    if ln in ("n_an", "nitrate", "no3") or (ln.startswith("n") and "_an" in ln):
        woa = woa_ds[name]
        break
if woa is None:
    for name in woa_ds.data_vars:
        if woa_ds[name].ndim >= 3:
            woa = woa_ds[name]
            break
if woa is None:
    woa = woa_ds[list(woa_ds.data_vars)[0]]

rename = {}
for d in list(woa.dims):
    dl = d.lower()
    if dl in ("lat", "latitude", "y"):
        new = "lat"
    elif dl in ("lon", "longitude", "x"):
        new = "lon"
    elif dl in ("depth", "lev", "z", "deptht"):
        new = "depth"
    else:
        new = d
    if d != new:
        rename[d] = new
if rename:
    woa = woa.rename(rename)

if "time" in woa.dims:
    woa = woa.isel(time=0)
woa = woa.squeeze(drop=True)
woa = woa.sortby("lat").sortby("lon")
woa = woa.where(woa < 1e20)

woa_lat = np.asarray(woa["lat"].values, dtype=float)
woa_lon = np.asarray(woa["lon"].values, dtype=float)

z_woa_raw = np.abs(np.asarray(woa["depth"].values, dtype=float))
order_z = np.argsort(z_woa_raw)
if not np.array_equal(order_z, np.arange(len(z_woa_raw))):
    woa = woa.isel(depth=order_z)
    z_woa_raw = z_woa_raw[order_z]
z_woa = z_woa_raw

landsea_path = os.path.join(woa_masks_root, "landsea_01.msk")
arr_ls = read_msk_first3(landsea_path)
lat_ls, lon_ls, bot_ls = clean_latlon_rows(arr_ls)
landsea_grid = table_to_grid(lat_ls, lon_ls, bot_ls, woa_lat, woa_lon)
woa_ocean_mask_np = landsea_grid > 1

basin_path = os.path.join(woa_masks_root, "basinmask_01.msk")
arr_b = read_msk_first3(basin_path)
lat_b, lon_b, code_b = clean_latlon_rows(arr_b)
basin_code_grid = table_to_grid(lat_b, lon_b, code_b, woa_lat, woa_lon)
code_filled = fill_code_nearest(basin_code_grid, max_iter=300)

lon2d_w, lat2d_w = np.meshgrid(woa_lon, woa_lat)
codes = np.unique(basin_code_grid[np.isfinite(basin_code_grid)])
code_to_parent = {}

for code in codes:
    m = basin_code_grid == code
    lat_c = float(np.nanmean(lat2d_w[m]))
    lon_rad = np.deg2rad(lon2d_w[m])
    lon_c = float(
        np.rad2deg(np.arctan2(np.nanmean(np.sin(lon_rad)),
                              np.nanmean(np.cos(lon_rad))))
    )

    if code <= 0:
        parent = np.nan
    elif lat_c < -50.0 or lat_c > 65.0:
        parent = np.nan
    else:
        if lat_c > 30.0 and lon_c < 45.0:
            parent = 1  # Atlantic / Mediterranean / Black Sea / Baltic
        elif lon_c > 110.0 and lat_c > 0.0:
            parent = 2  # Pacific marginal seas, e.g. South China Sea
        elif 20.0 <= lon_c <= 120.0:
            parent = 3  # Indian
        elif lon_c >= 120.0 or lon_c <= -80.0:
            parent = 2  # Pacific
        else:
            parent = 1  # Atlantic
    code_to_parent[code] = parent

parent_grid = np.full_like(code_filled, np.nan, dtype=float)
for code, parent in code_to_parent.items():
    parent_grid[code_filled == code] = parent

woa_ocean_mask = xr.DataArray(
    woa_ocean_mask_np, dims=["lat", "lon"],
    coords={"lat": woa_lat, "lon": woa_lon}
)
woa_basin = xr.DataArray(
    parent_grid, dims=["lat", "lon"],
    coords={"lat": woa_lat, "lon": woa_lon}
)
woa_area = xr.DataArray(
    compute_cell_areas(woa_lat, woa_lon), dims=["lat", "lon"],
    coords={"lat": woa_lat, "lon": woa_lon}
)

woa_nc, woa_grad = nitracline_from_field(woa, z_woa, "depth")
woa_nc = woa_nc.where(woa_ocean_mask)

global_means = {}
basin_means = {}
reg_fields = {}
metrics = {}
section_umol = {}
section_nc = {}

global_means["WOA23"] = area_mean(woa_nc, woa_area)
basin_means["WOA23"] = {}
for bname in basin_names:
    basin_means["WOA23"][bname] = area_mean(
        woa_nc, woa_area, woa_basin == basin_codes[bname]
    )

lat2d_w_da = xr.DataArray(
    lat2d_w, dims=("lat", "lon"),
    coords={"lat": woa_lat, "lon": woa_lon}
)
lon2d_w_da = xr.DataArray(
    lon2d_w, dims=("lat", "lon"),
    coords={"lat": woa_lat, "lon": woa_lon}
)
woa_sec = section_band(woa, woa_area, lat2d_w_da, lon2d_w_da, z_woa, "depth")
section_umol["WOA23"] = woa_sec
section_nc["WOA23"] = nitracline_section(
    woa_sec["data"], woa_sec["depth"]
)

# -----------------------------------------------------------------------------
# CMIP6 models
# -----------------------------------------------------------------------------
for model in models:
    info = models_info[model]

    no3_clim = get_cmip6_climatology(model, "no3")
    z_model = np.abs(
        np.asarray(no3_clim[info["depth_coord"]].values, dtype=float)
        * info["depth_scale"]
    )

    nc, grad = nitracline_from_field(
        no3_clim, z_model, info["depth_coord"]
    )

    area = open_cmip6_var(model, "areacello").squeeze(drop=True)
    area = area.where(area < 1e19)
    area_b = area.broadcast_like(nc)

    lat2d, lon2d = get_latlon2d(no3_clim, info)
    lon2d_180 = ((lon2d + 180.0) % 360.0) - 180.0

    lat_idx, lon_idx = nearest_indices_flat(
        lat2d.values, lon2d_180.values, woa_lat, woa_lon
    )
    model_basin_np = parent_grid[lat_idx, lon_idx]
    model_basin = xr.DataArray(
        model_basin_np,
        dims=nc.dims,
        coords={d: nc[d] for d in nc.dims if d in nc.coords},
    )

    global_means[model] = area_mean(nc, area_b)
    basin_means[model] = {}
    for bname in basin_names:
        basin_means[model][bname] = area_mean(
            nc, area_b, model_basin == basin_codes[bname]
        )

    reg_fields[model] = regrid_conservative(
        nc, area_b, lat_idx, lon_idx, woa_ocean_mask_np
    )

    # Figure 2 density conversion
    thetao = get_cmip6_climatology(model, "thetao")
    so = get_cmip6_climatology(model, "so")
    thetao, so = xr.align(thetao, so, join="exact")
    z_da = xr.DataArray(z_model, dims=[info["depth_coord"]])
    rho = simple_density(thetao, so, z_da)
    no3_umol = (no3_clim * 1.0e6 / rho).assign_attrs(
        units="micromoles_per_kilogram",
        density_assumption=(
            "In-situ density estimated from thetao/so with a linear EOS; "
            "nitracline depth itself is scale-invariant."
        ),
    )

    sec_native = section_band(
        no3_clim, area_b, lat2d, lon2d_180, z_model, info["depth_coord"]
    )
    sec_umol = section_band(
        no3_umol, area_b, lat2d, lon2d_180, z_model, info["depth_coord"]
    )

    section_umol[model] = sec_umol
    section_nc[model] = nitracline_section(
        sec_native["data"], sec_native["depth"]
    )

# -----------------------------------------------------------------------------
# Common mask and comparison metrics on the WOA 1-degree grid
# -----------------------------------------------------------------------------
woa_nc_grid = np.asarray(woa_nc.values, dtype=float)
woa_area_grid = np.asarray(woa_area.values, dtype=float)

shared = np.isfinite(woa_nc_grid) & woa_ocean_mask_np
for model in models:
    shared &= np.isfinite(reg_fields[model])

shared_cell_count = int(np.sum(shared))

for model in models:
    corr = weighted_corr_np(
        reg_fields[model], woa_nc_grid, woa_area_grid, shared
    )
    rmse = weighted_rmse_np(
        reg_fields[model], woa_nc_grid, woa_area_grid, shared
    )
    metrics[model] = {"corr": corr, "rmse": rmse}

best_model = min(metrics, key=lambda m: metrics[m]["rmse"])

woa_order = tuple(
    sorted(basin_names, key=lambda b: basin_means["WOA23"][b], reverse=True)
)
basin_order_match = {}
for model in models:
    model_order = tuple(
        sorted(basin_names, key=lambda b: basin_means[model][b], reverse=True)
    )
    basin_order_match[model] = (model_order == woa_order)

map_fields = {model: reg_fields[model] for model in models}
map_fields["WOA23"] = woa_nc_grid
section_order = ["WOA23", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

# -----------------------------------------------------------------------------
# Report
# -----------------------------------------------------------------------------
print("Global mean nitracline depth (m, nearest metre):")
for src in ["WOA23"] + models:
    print(f"  {src:16s} {global_means[src]:8.0f}")
print("Model minus WOA23 (m, nearest metre):")
for model in models:
    print(f"  {model:16s} {global_means[model] - global_means['WOA23']:+8.0f}")

print("\nBasin mean nitracline depth (m, nearest metre):")
for bname in basin_names:
    print(f"  {bname}:")
    for src in ["WOA23"] + models:
        print(f"    {src:16s} {basin_means[src][bname]:8.0f}")
    for model in models:
        d = basin_means[model][bname] - basin_means["WOA23"][bname]
        print(f"    {model:16s} diff {d:+8.0f}")

print("\nSpatial comparison on WOA 1-degree grid:")
print(f"  Shared valid cells: {shared_cell_count}")
for model in models:
    print(
        f"  {model:16s} corr = {metrics[model]['corr']:.3f}, "
        f"RMSE = {metrics[model]['rmse']:.1f} m"
    )
print(f"\nBest model by area-weighted RMSE: {best_model}")
print(f"WOA23 basin order (deepest to shallowest): {woa_order}")
for model in models:
    model_order = tuple(
        sorted(basin_names, key=lambda b: basin_means[model][b], reverse=True)
    )
    print(
        f"  {model:16s} order = {model_order}, "
        f"matches WOA23 = {basin_order_match[model]}"
    )
print(
    "\nNote: WOA23 nitrate is the all-period annual composite; model fields "
    "are 1971-2000 climatologies, so the comparison carries an unavoidable "
    "period mismatch."
)
