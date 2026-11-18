import os
import glob
import re
import warnings
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from scipy import stats
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore", category=RuntimeWarning)

# -----------------------------------------------------------------------------
# Paths and model metadata
# -----------------------------------------------------------------------------
woa_root = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
woa_mask_root = "/data1/lorraine-20419/ocean_data/WOA23/masks"
cmip_hist_root = "/data1/lorraine-20419/benchmark_cmip6/historical"
cmip_ssp_root = "/data1/lorraine-20419/benchmark_cmip6/ssp245"

models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
model_colors = {"CESM2": "C0", "GFDL-ESM4": "C1", "MPI-ESM1-2-LR": "C2"}

model_meta = {
    "CESM2": {
        "lat_coord": "lat",
        "lon_coord": "lon",
        "depth_scale": 0.01,  # lev in centimetres
    },
    "GFDL-ESM4": {
        "lat_coord": "lat",
        "lon_coord": "lon",
        "depth_scale": 1.0,
    },
    "MPI-ESM1-2-LR": {
        "lat_coord": "latitude",
        "lon_coord": "longitude",
        "depth_scale": 1.0,
    },
}

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def lonlat_to_xyz(lat, lon):
    lat = np.deg2rad(np.asarray(lat, dtype=float))
    lon = np.deg2rad(np.asarray(lon, dtype=float))
    return np.column_stack(
        [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)]
    )


def make_regrid_indices(src_lat, src_lon, tgt_lat, tgt_lon):
    src_lat = np.asarray(src_lat)
    src_lon = np.asarray(src_lon)
    if src_lat.ndim == 1 and src_lon.ndim == 1:
        src_lon2, src_lat2 = np.meshgrid(src_lon, src_lat)
    else:
        src_lat2, src_lon2 = src_lat, src_lon
    tgt_lat2, tgt_lon2 = np.meshgrid(tgt_lat, tgt_lon, indexing="ij")
    src_xyz = lonlat_to_xyz(src_lat2.ravel(), src_lon2.ravel())
    tgt_xyz = lonlat_to_xyz(tgt_lat2.ravel(), tgt_lon2.ravel())
    tree = cKDTree(src_xyz)
    _, idx = tree.query(tgt_xyz, k=1)
    return idx.reshape(len(tgt_lat), len(tgt_lon))


def regrid_array(arr, idx):
    arr = np.asarray(arr)
    lead = arr.shape[:-2]
    nsp = arr.shape[-2] * arr.shape[-1]
    flat = arr.reshape(-1, nsp)
    out = flat[:, idx.ravel()]
    return out.reshape(lead + idx.shape)


def weighted_corr(x, y, w):
    x = np.asarray(x).ravel()
    y = np.asarray(y).ravel()
    w = np.asarray(w).ravel()
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(w)
    x, y, w = x[valid], y[valid], w[valid]
    if len(x) < 2:
        return np.nan
    wsum = w.sum()
    mx = np.sum(w * x) / wsum
    my = np.sum(w * y) / wsum
    cov = np.sum(w * (x - mx) * (y - my)) / wsum
    vx = np.sum(w * (x - mx) ** 2) / wsum
    vy = np.sum(w * (y - my) ** 2) / wsum
    if vx <= 0 or vy <= 0:
        return np.nan
    return cov / np.sqrt(vx * vy)


def linear_trend_ess(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid], y[valid]
    n = len(x)
    if n < 3:
        return dict(
            slope=np.nan, intercept=np.nan, ci_low=np.nan, ci_high=np.nan,
            r2=np.nan, excludes_zero=False, p_value=np.nan, n=n,
            n_eff=np.nan, rho=np.nan,
        )
    x_mean = x.mean()
    y_mean = y.mean()
    Sxx = np.sum((x - x_mean) ** 2)
    Sxy = np.sum((x - x_mean) * (y - y_mean))
    slope = Sxy / Sxx
    intercept = y_mean - slope * x_mean
    resid = y - (intercept + slope * x)
    sse = np.sum(resid ** 2)

    denom = np.sqrt(np.sum(resid[:-1] ** 2) * np.sum(resid[1:] ** 2))
    rho = np.sum(resid[:-1] * resid[1:]) / denom if denom > 0 else 0.0
    rho = float(np.clip(rho, -0.99, 0.99))
    n_eff = n * (1 - rho) / (1 + rho)
    n_eff = float(np.clip(n_eff, 2.0, n))

    df_base = max(n - 2, 1)
    sigma2 = sse / df_base
    se = np.sqrt(sigma2 / Sxx * (n / n_eff))
    df = max(n_eff - 2.0, 1.0)
    tcrit = stats.t.ppf(0.975, df)
    ci = tcrit * se

    syy = np.sum((y - y_mean) ** 2)
    r2 = 1.0 - sse / syy if syy > 0 else np.nan
    pval = 2.0 * stats.t.sf(abs(slope / se), df) if se > 0 else np.nan
    excludes_zero = bool((slope - ci > 0) or (slope + ci < 0))

    return dict(
        slope=float(slope),
        intercept=float(intercept),
        ci_low=float(slope - ci),
        ci_high=float(slope + ci),
        r2=float(r2),
        excludes_zero=excludes_zero,
        p_value=float(pval),
        n=int(n),
        n_eff=float(n_eff),
        rho=float(rho),
    )


def trend_map_ess(annual_da):
    years = annual_da["year"].values.astype(float)
    spatial_dims = [d for d in annual_da.dims if d != "year"]
    arr = annual_da.transpose("year", *spatial_dims).values
    n = arr.shape[0]
    y = arr.reshape(n, -1)
    x = years - years.mean()
    Sxx = np.sum(x ** 2)

    valid = np.isfinite(y)
    count = valid.sum(axis=0)
    y_sum = np.nansum(y, axis=0)
    y_mean = np.divide(
        y_sum, count,
        out=np.full_like(y_sum, np.nan, dtype=float),
        where=count > 0,
    )

    y_centered = y - y_mean[None, :]
    Sxy = np.nansum(x[:, None] * y_centered, axis=0)
    slope = Sxy / Sxx
    intercept = y_mean - slope * years.mean()
    resid = y - (intercept[None, :] + slope[None, :] * years[:, None])
    sse = np.nansum(resid ** 2, axis=0)

    num = np.nansum(resid[:-1, :] * resid[1:, :], axis=0)
    den1 = np.nansum(resid[:-1, :] ** 2, axis=0)
    den2 = np.nansum(resid[1:, :] ** 2, axis=0)
    denom = np.sqrt(den1 * den2)
    rho = np.divide(num, denom, out=np.zeros_like(num), where=denom > 0)
    rho = np.clip(rho, -0.99, 0.99)
    n_eff = n * (1 - rho) / (1 + rho)
    n_eff = np.clip(n_eff, 2, n)

    df_base = np.maximum(count - 2, 1)
    sigma2 = sse / df_base
    sigma2 = np.where(count >= 3, sigma2, np.nan)
    se = np.sqrt(sigma2 / Sxx * (n / n_eff))
    df = np.maximum(n_eff - 2, 1)
    pval = 2.0 * stats.t.sf(np.abs(slope / se), df)
    slope = np.where(count >= 3, slope, np.nan)
    pval = np.where(count >= 3, pval, np.nan)

    return slope.reshape(arr.shape[1:]), pval.reshape(arr.shape[1:])


def open_so_files(model, experiment, root):
    d = os.path.join(root, model, "so")
    files = sorted(glob.glob(os.path.join(d, "so_*.nc")))
    if not files:
        raise FileNotFoundError(f"No so files found for {model} {experiment} under {d}")
    variants = {}
    for f in files:
        m = re.search(r"_(r\d+i\d+p\d+f\d+)_", os.path.basename(f))
        key = m.group(1) if m else "unknown"
        variants.setdefault(key, []).append(f)
    variant = max(variants, key=lambda k: len(variants[k]))
    files = variants[variant]
    ds = xr.open_mfdataset(
        files, combine="by_coords", use_cftime=True, decode_times=True
    )
    return ds, variant, files


def open_fixed_field(model, var, root, fallback_root=None, table="Ofx"):
    dirs = [os.path.join(root, model, var)]
    if fallback_root is not None:
        dirs.append(os.path.join(fallback_root, model, var))
    for d in dirs:
        if not os.path.isdir(d):
            continue
        files = sorted(glob.glob(os.path.join(d, f"{var}_*.nc")))
        if not files:
            continue
        tab_files = [f for f in files if f"_{table}_" in os.path.basename(f)]
        if tab_files:
            files = tab_files

        # GFDL volcello is Omon and split in time; one file already has the full grid.
        if var == "volcello" and table == "Omon":
            ds = xr.open_dataset(files[0], use_cftime=True)
        elif len(files) == 1:
            ds = xr.open_dataset(files[0], use_cftime=True)
        else:
            ds = xr.open_mfdataset(files, combine="by_coords", use_cftime=True)
        if var in ds:
            return ds[var]
    return None


def grid_from_points(lat_vals, lon_vals, vals, tgt_lat, tgt_lon):
    grid = np.full((len(tgt_lat), len(tgt_lon)), np.nan, dtype=float)
    dlat = tgt_lat[1] - tgt_lat[0]
    dlon = tgt_lon[1] - tgt_lon[0]
    ilat = np.round((np.asarray(lat_vals) - tgt_lat[0]) / dlat).astype(int)
    ilon = np.round((np.asarray(lon_vals) - tgt_lon[0]) / dlon).astype(int)
    ok = (
        (ilat >= 0) & (ilat < len(tgt_lat))
        & (ilon >= 0) & (ilon < len(tgt_lon))
    )
    grid[ilat[ok], ilon[ok]] = np.asarray(vals)[ok]
    return grid


def fold_basin_mask(basin_grid, lat, lon):
    major_codes = np.array([1, 2, 3, 4, 5])
    folded = np.full_like(basin_grid, np.nan, dtype=float)

    for c in major_codes:
        folded[basin_grid == c] = c

    lat2d, lon2d = np.meshgrid(lat, lon, indexing="ij")
    major_mask = np.isin(basin_grid, major_codes)
    if not major_mask.any():
        return folded

    major_lat = lat2d[major_mask]
    major_lon = lon2d[major_mask]
    major_vals = basin_grid[major_mask]
    tree = cKDTree(lonlat_to_xyz(major_lat, major_lon))

    high_codes = [
        c for c in np.unique(basin_grid)
        if np.isfinite(c) and c > 0 and c not in major_codes
    ]
    for c in high_codes:
        mask_c = basin_grid == c
        if not mask_c.any():
            continue
        clat = lat2d[mask_c]
        clon = lon2d[mask_c]
        _, nearest_idx = tree.query(lonlat_to_xyz(clat, clon), k=1)
        nearest_major = major_vals[nearest_idx]
        vals, counts = np.unique(nearest_major, return_counts=True)
        parent = vals[np.argmax(counts)]
        folded[mask_c] = parent
    return folded


# -----------------------------------------------------------------------------
# WOA23: surface salinity, area, land-sea mask, basin mask
# -----------------------------------------------------------------------------
sal_file = os.path.join(
    woa_root, "salinity", "1971-2000_climate_normal", "woa23_decav71A0_s00_01.nc"
)
ds_sal = xr.open_dataset(sal_file, decode_times=False)

candidates = [v for v in ds_sal.data_vars if v.startswith("s") and "an" in v]
if not candidates:
    candidates = [v for v in ds_sal.data_vars if "an" in v]
if not candidates:
    candidates = list(ds_sal.data_vars)
sal_name = candidates[0]
sal = ds_sal[sal_name]

if "time" in sal.dims:
    sal = sal.isel(time=0)
sal = sal.drop_vars("time", errors="ignore")

depth = ds_sal["depth"].values

if "depth_bnds" in ds_sal:
    bnds = ds_sal["depth_bnds"].values
    if bnds.ndim == 3:
        bnds = bnds[0]
    if bnds.shape[0] != len(depth) and bnds.shape[1] == len(depth):
        bnds = bnds.T
    thickness = np.abs(bnds[:, 1] - bnds[:, 0])
else:
    edges = np.empty(len(depth) + 1, dtype=float)
    edges[1:-1] = 0.5 * (depth[:-1] + depth[1:])
    edges[0] = depth[0] - (depth[1] - depth[0]) / 2.0
    edges[-1] = depth[-1] + (depth[-1] - depth[-2]) / 2.0
    thickness = np.diff(edges)

sel_idx_woa = np.where((depth >= 0.0) & (depth <= 10.0))[0]
if len(sel_idx_woa) == 0:
    raise ValueError("No WOA23 levels found in 0-10 m")

sal_sel = sal.isel(depth=sel_idx_woa).drop_vars("depth", errors="ignore")

lat_woa = ds_sal["lat"].values
lon_woa = ds_sal["lon"].values

dlon_rad = np.deg2rad(np.abs(lon_woa[1] - lon_woa[0])) if len(lon_woa) > 1 else np.deg2rad(1.0)
lat_edges = np.deg2rad(
    np.concatenate(
        [
            [lat_woa[0] - (lat_woa[1] - lat_woa[0]) / 2.0],
            lat_woa + (lat_woa[1] - lat_woa[0]) / 2.0,
        ]
    )
)
area_lat = 6371000.0 ** 2 * dlon_rad * (
    np.sin(lat_edges[1:]) - np.sin(lat_edges[:-1])
)
area_woa_np = np.tile(area_lat[:, None], (1, len(lon_woa)))
area_woa = xr.DataArray(
    area_woa_np, dims=("lat", "lon"),
    coords={"lat": lat_woa, "lon": lon_woa}, name="woa_area",
)

vol_woa_np = area_woa_np[None, :, :] * thickness[sel_idx_woa][:, None, None]
vol_woa = xr.DataArray(
    vol_woa_np, dims=("depth", "lat", "lon"),
    coords={"depth": depth[sel_idx_woa], "lat": lat_woa, "lon": lon_woa},
    name="woa_level_volume",
)

num_woa = (sal_sel * vol_woa).sum("depth")
den_woa = vol_woa.sum("depth")
woa_sss = (num_woa / den_woa).assign_attrs(units="psu")

# WOA land-sea mask
landsea_path = os.path.join(woa_mask_root, "landsea_01.msk")
arr_ls = np.loadtxt(landsea_path, skiprows=2, delimiter=",")
landsea_grid = grid_from_points(arr_ls[:, 0], arr_ls[:, 1], arr_ls[:, 2], lat_woa, lon_woa)
landsea_da = xr.DataArray(
    landsea_grid, dims=("lat", "lon"),
    coords={"lat": lat_woa, "lon": lon_woa}, name="landsea",
)

# WOA basin mask, surface level
basin_path = os.path.join(woa_mask_root, "basinmask_01.msk")
arr_b = np.loadtxt(basin_path, skiprows=2, delimiter=",")

if arr_b.shape[1] > 10:
    lat_b = arr_b[:, 0]
    lon_b = arr_b[:, 1]
    basin_surface_vals = arr_b[:, 2]
else:
    lat_b_all = arr_b[:, 0]
    lon_b_all = arr_b[:, 1]
    lev_b = arr_b[:, 2]
    basin_all = arr_b[:, 3]
    sel_b = lev_b == np.min(lev_b)
    lat_b = lat_b_all[sel_b]
    lon_b = lon_b_all[sel_b]
    basin_surface_vals = basin_all[sel_b]

basin_grid = grid_from_points(lat_b, lon_b, basin_surface_vals, lat_woa, lon_woa)
basin_folded_np = fold_basin_mask(basin_grid, lat_woa, lon_woa)
basin_folded_da = xr.DataArray(
    basin_folded_np, dims=("lat", "lon"),
    coords={"lat": lat_woa, "lon": lon_woa}, name="basin_folded",
)

# -----------------------------------------------------------------------------
# CMIP6 processing per model
# -----------------------------------------------------------------------------
model_climo = {}
model_delta = {}
model_area_regrid = {}
model_trend_regrid = {}
model_pval_regrid = {}
gm_hist = {}
gm_ssp = {}
member_labels = {}
trend_stats_hist = {}
trend_stats_ssp = {}

for model in models:
    # Fixed fields, opened once and used for both experiments.
    area_native = open_fixed_field(model, "areacello", cmip_hist_root, table="Ofx")
    if area_native is None:
        raise FileNotFoundError(f"areacello not found for {model}")
    if "time" in area_native.dims:
        area_native = area_native.isel(time=0)
    area_native = area_native.drop_vars("time", errors="ignore")

    sftof = None
    if model == "CESM2":
        sftof = open_fixed_field(model, "sftof", cmip_hist_root, table="Ofx")
        if sftof is not None:
            if "time" in sftof.dims:
                sftof = sftof.isel(time=0)
            sftof = sftof.drop_vars("time", errors="ignore")
            area_native = area_native.where(sftof > 0)

    area_native = area_native.where(area_native > 0)

    if model == "GFDL-ESM4":
        vol = open_fixed_field(model, "volcello", cmip_hist_root, table="Omon")
    else:
        vol = open_fixed_field(model, "volcello", cmip_hist_root, table="Ofx")
    if vol is None:
        raise FileNotFoundError(f"volcello not found for {model}")
    if "time" in vol.dims:
        vol = vol.isel(time=0)
    vol = vol.drop_vars("time", errors="ignore")
    vol = vol.drop_vars("lev", errors="ignore")

    if sftof is not None:
        vol = vol.where(sftof > 0)
    vol = vol.where(vol > 0)

    # Historical so, 1950-2014
    ds_hist, var_hist, _ = open_so_files(model, "historical", cmip_hist_root)
    member_labels[(model, "historical")] = var_hist
    time_mask = (ds_hist.time.dt.year >= 1950) & (ds_hist.time.dt.year <= 2014)
    ds_hist = ds_hist.isel(time=time_mask.values)

    so_hist = ds_hist["so"]
    lev_hist = ds_hist["lev"].values
    lev_m_hist = lev_hist * model_meta[model]["depth_scale"]
    lev_idx_hist = np.where((lev_m_hist >= 0.0) & (lev_m_hist <= 10.0))[0]
    if len(lev_idx_hist) == 0:
        raise ValueError(f"No {model} historical levels found in 0-10 m")

    so_hist_sel = so_hist.isel(lev=lev_idx_hist).drop_vars("lev", errors="ignore")
    vol_hist_sel = vol.isel(lev=lev_idx_hist)
    num = (so_hist_sel * vol_hist_sel).sum("lev")
    den = vol_hist_sel.sum("lev")
    sss_hist_monthly = (num / den.where(den > 0)).assign_attrs(units="psu")
    sss_hist_annual = sss_hist_monthly.groupby("time.year").mean("time").compute()

    # ssp245 so, 2071-2100
    ds_ssp, var_ssp, _ = open_so_files(model, "ssp245", cmip_ssp_root)
    member_labels[(model, "ssp245")] = var_ssp
    time_mask = (ds_ssp.time.dt.year >= 2071) & (ds_ssp.time.dt.year <= 2100)
    ds_ssp = ds_ssp.isel(time=time_mask.values)

    so_ssp = ds_ssp["so"]
    lev_ssp = ds_ssp["lev"].values
    lev_m_ssp = lev_ssp * model_meta[model]["depth_scale"]
    lev_idx_ssp = np.where((lev_m_ssp >= 0.0) & (lev_m_ssp <= 10.0))[0]
    if len(lev_idx_ssp) == 0:
        raise ValueError(f"No {model} ssp245 levels found in 0-10 m")

    so_ssp_sel = so_ssp.isel(lev=lev_idx_ssp).drop_vars("lev", errors="ignore")
    vol_ssp_sel = vol.isel(lev=lev_idx_ssp)
    num = (so_ssp_sel * vol_ssp_sel).sum("lev")
    den = vol_ssp_sel.sum("lev")
    sss_ssp_monthly = (num / den.where(den > 0)).assign_attrs(units="psu")
    sss_ssp_annual = sss_ssp_monthly.groupby("time.year").mean("time").compute()

    # Global means on native grid using native areacello
    spatial_dims = [d for d in sss_hist_annual.dims if d != "year"]
    num_gm = (sss_hist_annual * area_native).sum(spatial_dims)
    den_gm = area_native.where(sss_hist_annual.notnull()).sum(spatial_dims)
    gm_hist[model] = (num_gm / den_gm).compute()

    num_gm = (sss_ssp_annual * area_native).sum(spatial_dims)
    den_gm = area_native.where(sss_ssp_annual.notnull()).sum(spatial_dims)
    gm_ssp[model] = (num_gm / den_gm).compute()

    trend_stats_hist[model] = linear_trend_ess(
        gm_hist[model]["year"].values, gm_hist[model].values
    )
    trend_stats_ssp[model] = linear_trend_ess(
        gm_ssp[model]["year"].values, gm_ssp[model].values
    )

    # Cell-by-cell annual trend over 1950-2014 on the native grid
    slope_map, pval_map = trend_map_ess(sss_hist_annual)

    # Regrid indices from native model grid to the WOA23 1-degree grid
    lat_coord = model_meta[model]["lat_coord"]
    lon_coord = model_meta[model]["lon_coord"]
    src_lat = ds_hist[lat_coord].values
    src_lon = ds_hist[lon_coord].values
    idx = make_regrid_indices(src_lat, src_lon, lat_woa, lon_woa)

    # 1971-2000 climatology: annual cell means first, then period mean
    hist_1971_2000 = sss_hist_annual.sel(year=slice(1971, 2000)).values
    hist_1971_2000_regrid = regrid_array(hist_1971_2000, idx)
    climo_regrid = np.nanmean(hist_1971_2000_regrid, axis=0)
    model_climo[model] = xr.DataArray(
        climo_regrid, dims=("lat", "lon"),
        coords={"lat": lat_woa, "lon": lon_woa}, name=f"{model}_climo",
    )

    # Delta SSS: ssp245 2071-2100 minus historical 1995-2014
    hist_1995_2014 = sss_hist_annual.sel(year=slice(1995, 2014)).values
    ssp_2071_2100 = sss_ssp_annual.sel(year=slice(2071, 2100)).values
    hist_1995_2014_regrid = regrid_array(hist_1995_2014, idx)
    ssp_2071_2100_regrid = regrid_array(ssp_2071_2100, idx)
    hist_mean = np.nanmean(hist_1995_2014_regrid, axis=0)
    ssp_mean = np.nanmean(ssp_2071_2100_regrid, axis=0)
    delta_np = ssp_mean - hist_mean
    delta_np = np.where(
        np.isfinite(hist_mean) & np.isfinite(ssp_mean), delta_np, np.nan
    )
    model_delta[model] = xr.DataArray(
        delta_np, dims=("lat", "lon"),
        coords={"lat": lat_woa, "lon": lon_woa}, name=f"{model}_delta",
    )

    # Regridded areacello
    area_regrid_np = regrid_array(area_native.values, idx)
    model_area_regrid[model] = xr.DataArray(
        area_regrid_np, dims=("lat", "lon"),
        coords={"lat": lat_woa, "lon": lon_woa}, name=f"{model}_area",
    )

    # Regrid the display trend map to WOA23 1-degree grid
    trend_regrid = regrid_array(slope_map * 10.0, idx)  # psu per decade
    pval_regrid = regrid_array(pval_map, idx)
    model_trend_regrid[model] = xr.DataArray(
        trend_regrid, dims=("lat", "lon"),
        coords={"lat": lat_woa, "lon": lon_woa}, name=f"{model}_trend",
    )
    model_pval_regrid[model] = xr.DataArray(
        pval_regrid, dims=("lat", "lon"),
        coords={"lat": lat_woa, "lon": lon_woa}, name=f"{model}_pval",
    )

# -----------------------------------------------------------------------------
# Common comparison mask and mean-state metrics
# -----------------------------------------------------------------------------
mask_common = (landsea_da > 1) & np.isfinite(woa_sss)
for model in models:
    mask_common = mask_common & np.isfinite(model_climo[model])

common_mask_cells = int(mask_common.sum().values)

woa_valid = mask_common & np.isfinite(area_woa) & np.isfinite(woa_sss)
woa_gm = float(
    (woa_sss.where(woa_valid) * area_woa.where(woa_valid)).sum()
    / area_woa.where(woa_valid).sum()
)

metrics = {}
bias_maps = {}

for model in models:
    valid = (
        woa_valid
        & np.isfinite(model_climo[model])
        & np.isfinite(model_area_regrid[model])
    )
    woa_m = woa_sss.where(valid)
    model_m = model_climo[model].where(valid)
    w_model = model_area_regrid[model].where(valid)

    model_gm = float((model_m * w_model).sum() / w_model.sum())
    bias = (model_climo[model] - woa_sss).where(valid)
    rmse = float(np.sqrt((bias ** 2 * w_model).sum() / w_model.sum()))
    corr = float(weighted_corr(model_m.values, woa_m.values, w_model.values))

    bias_maps[model] = bias
    metrics[model] = dict(
        global_mean=model_gm,
        woa_mean=woa_gm,
        diff=model_gm - woa_gm,
        rmse=rmse,
        corr=corr,
        metric_cells=int(valid.sum().values),
    )

# -----------------------------------------------------------------------------
# Delta SSS global and basin means
# -----------------------------------------------------------------------------
delta_results = {}
basin_codes = {"Atlantic": 1, "Indian": 3, "Pacific": 2}

for model in models:
    delta = model_delta[model]
    area = model_area_regrid[model]
    valid = np.isfinite(delta) & np.isfinite(area)

    global_delta = float(
        (delta.where(valid) * area.where(valid)).sum()
        / area.where(valid).sum()
    )

    basins = {}
    for bname, code in basin_codes.items():
        bmask = valid & (basin_folded_da == code)
        if int(bmask.sum().values) == 0:
            basins[bname] = np.nan
        else:
            basins[bname] = float(
                (delta.where(bmask) * area.where(bmask)).sum()
                / area.where(bmask).sum()
            )

    delta_results[model] = dict(
        global_delta=global_delta,
        basins=basins,
        valid_cells=int(valid.sum().values),
    )

# -----------------------------------------------------------------------------
# Trend summaries, slope ratios and member caveats
# -----------------------------------------------------------------------------
trend_summary = {}
member_caveat = {}

for model in models:
    h = trend_stats_hist[model]
    s = trend_stats_ssp[model]

    h_slope_dec = h["slope"] * 10.0
    h_ci_low_dec = h["ci_low"] * 10.0
    h_ci_high_dec = h["ci_high"] * 10.0

    s_slope_dec = s["slope"] * 10.0
    s_ci_low_dec = s["ci_low"] * 10.0
    s_ci_high_dec = s["ci_high"] * 10.0

    ratio = s["slope"] / h["slope"] if abs(h["slope"]) > 0 else np.nan

    trend_summary[model] = dict(
        historical=dict(
            slope_per_decade=h_slope_dec,
            ci_low_per_decade=h_ci_low_dec,
            ci_high_per_decade=h_ci_high_dec,
            r2=h["r2"],
            excludes_zero=h["excludes_zero"],
            p_value=h["p_value"],
            n=h["n"],
            n_eff=h["n_eff"],
            rho=h["rho"],
            slope_per_year=h["slope"],
            intercept=h["intercept"],
        ),
        ssp245=dict(
            slope_per_decade=s_slope_dec,
            ci_low_per_decade=s_ci_low_dec,
            ci_high_per_decade=s_ci_high_dec,
            r2=s["r2"],
            excludes_zero=s["excludes_zero"],
            p_value=s["p_value"],
            n=s["n"],
            n_eff=s["n_eff"],
            rho=s["rho"],
            slope_per_year=s["slope"],
            intercept=s["intercept"],
        ),
        ratio_ssp_over_hist=float(ratio),
    )

    hm = member_labels[(model, "historical")]
    sm = member_labels[(model, "ssp245")]
    if hm != sm:
        member_caveat[model] = (
            f"{model}: historical uses {hm}, ssp245 uses {sm}. "
            "The two experiments are different ensemble members, so delta SSS and the "
            "ssp245/historical slope ratio include internal-variability differences between members."
        )
    else:
        member_caveat[model] = (
            f"{model}: historical and ssp245 both use {hm}. "
            "The experiments remain separate periods; the slope ratio is not a continuous "
            "1950-2100 trend or an acceleration estimate, and the 2015-2070 gap is not bridged."
        )

# -----------------------------------------------------------------------------
# Store plot-ready variables and a results dictionary
# -----------------------------------------------------------------------------
global_mean_hist = gm_hist
global_mean_ssp = gm_ssp
trend_maps = model_trend_regrid
pval_maps = model_pval_regrid

results = {
    "woa_global_mean": woa_gm,
    "common_mask_cells": common_mask_cells,
    "metrics": metrics,
    "trend_summary": trend_summary,
    "delta": delta_results,
    "member_labels": member_labels,
    "member_caveat": member_caveat,
}

# Print a compact summary
print("=== WOA23 1971-2000 surface salinity ===")
print(f"WOA23 global mean SSS: {woa_gm:.2f} psu")
print(f"Common comparison mask cells: {common_mask_cells}")
print()

print("=== Mean-state comparison, 1971-2000 ===")
for model in models:
    m = metrics[model]
    print(
        f"{model}: mean={m['global_mean']:.2f} psu, "
        f"model-WOA={m['diff']:.2f} psu, "
        f"RMSE={m['rmse']:.2f} psu, "
        f"spatial r={m['corr']:.3f}, "
        f"metric cells={m['metric_cells']}"
    )
print()

print("=== Global-mean annual SSS trends ===")
for model in models:
    t = trend_summary[model]
    h = t["historical"]
    s = t["ssp245"]
    print(f"{model}:")
    print(
        f"  historical 1950-2014: {h['slope_per_decade']:.3g} psu/decade "
        f"[{h['ci_low_per_decade']:.3g}, {h['ci_high_per_decade']:.3g}], "
        f"R2={h['r2']:.3f}, excludes zero={h['excludes_zero']}, "
        f"n_eff={h['n_eff']:.1f}"
    )
    print(
        f"  ssp245 2071-2100:     {s['slope_per_decade']:.3g} psu/decade "
        f"[{s['ci_low_per_decade']:.3g}, {s['ci_high_per_decade']:.3g}], "
        f"R2={s['r2']:.3f}, excludes zero={s['excludes_zero']}, "
        f"n_eff={s['n_eff']:.1f}"
    )
    print(f"  ssp245/historical slope ratio: {t['ratio_ssp_over_hist']:.1f}")
    print(f"  caveat: {member_caveat[model]}")
print()

print("=== Delta SSS, ssp245 2071-2100 minus historical 1995-2014 ===")
for model in models:
    d = delta_results[model]
    b = d["basins"]
    print(
        f"{model}: global={d['global_delta']:.2f} psu, "
        f"Atlantic={b['Atlantic']:.2f} psu, "
        f"Indian={b['Indian']:.2f} psu, "
        f"Pacific={b['Pacific']:.2f} psu "
        f"(valid cells={d['valid_cells']})"
    )
