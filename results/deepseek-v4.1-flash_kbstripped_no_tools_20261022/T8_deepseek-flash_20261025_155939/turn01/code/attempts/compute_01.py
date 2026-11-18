import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import gsw
from scipy.interpolate import interp1d

# ------------------------------------------------------------------
# Paths and model metadata
# ------------------------------------------------------------------
CMIP6_ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"
GLODAP_ROOT = "/data1/lorraine-20419/ocean_data/GLODAPv2.2016b_MappedClimatologies"

models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
model_info = {
    "CESM2": dict(y="nlat", x="nlon", lat="lat", lon="lon", depth_scale=0.01),
    "GFDL-ESM4": dict(y="lat", x="lon", lat="lat", lon="lon", depth_scale=1.0),
    "MPI-ESM1-2-LR": dict(y="j", x="i", lat="latitude", lon="longitude", depth_scale=1.0),
}

depths = np.array([0.0, 100.0, 500.0, 1000.0, 2000.0])

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
def find_cmip_files(model, var):
    d = os.path.join(CMIP6_ROOT, model, var)
    return sorted(glob.glob(os.path.join(d, "*.nc")))

def open_cmip_dataset(model, var):
    files = find_cmip_files(model, var)
    if len(files) == 0:
        return None
    if len(files) == 1:
        ds = xr.open_dataset(files[0], use_cftime=True, decode_times=True)
    else:
        ds = xr.open_mfdataset(
            files, combine="by_coords", use_cftime=True, decode_times=True
        )
    return ds

def annual_mean(da, year):
    if da is None:
        return None
    if "time" not in da.dims:
        return da
    yrs = da.time.dt.year.values
    idx = np.where(yrs == year)[0]
    if len(idx) == 0:
        raise ValueError(f"No time steps found for year {year}")
    return da.isel(time=idx).mean("time", keep_attrs=True)

def period_mean(da, y0, y1):
    if da is None:
        return None
    if "time" not in da.dims:
        return da
    yrs = da.time.dt.year.values
    idx = np.where((yrs >= y0) & (yrs <= y1))[0]
    if len(idx) == 0:
        raise ValueError(f"No time steps found for period {y0}-{y1}")
    return da.isel(time=idx).mean("time", keep_attrs=True)

def to_3d_np(da, mi):
    da = da.squeeze(drop=True)
    if "time" in da.dims:
        da = da.mean("time")
    for d in ["lev", mi["y"], mi["x"]]:
        if d not in da.dims:
            raise ValueError(f"Expected dimension {d!r} not found in {da.dims}")
    return da.transpose("lev", mi["y"], mi["x"]).values.astype(float)

def get_lat_lon(ds, mi):
    lat = ds[mi["lat"]] if mi["lat"] in ds else ds.coords[mi["lat"]]
    lon = ds[mi["lon"]] if mi["lon"] in ds else ds.coords[mi["lon"]]
    if lat.ndim == 1 and lon.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon.values, lat.values)
    elif lat.ndim == 2 and lon.ndim == 2:
        lat2d = lat.values
        lon2d = lon.values
    else:
        lat2d, lon2d = xr.broadcast(lat, lon)
        lat2d = lat2d.values
        lon2d = lon2d.values
    return lat2d.astype(float), lon2d.astype(float)

def compute_density(theta, so, lev_m, lat2d, lon2d):
    p = gsw.p_from_z(-lev_m[:, None, None], lat2d[None, :, :])
    SA = gsw.SA_from_SP(so, p, lon2d[None, :, :], lat2d[None, :, :])
    CT = gsw.CT_from_pt(SA, theta)
    rho = gsw.rho(SA, CT, p)
    return rho

def column_mean_0_50(var, vol, lev_m, max_depth=50.0):
    sel = (lev_m >= 0.0) & (lev_m <= max_depth)
    if not np.any(sel):
        return np.full(var.shape[1:], np.nan), np.zeros(var.shape[1:])
    v = var[sel]
    w = vol[sel]
    valid = np.isfinite(v) & np.isfinite(w) & (w > 0.0)
    num = np.sum(np.where(valid, v * w, 0.0), axis=0)
    den = np.sum(np.where(valid, w, 0.0), axis=0)
    out = np.full(num.shape, np.nan, dtype=float)
    ok = den > 0.0
    out[ok] = num[ok] / den[ok]
    return out, den

def weighted_global(col, colvol):
    valid = np.isfinite(col) & np.isfinite(colvol) & (colvol > 0.0)
    if not np.any(valid):
        return np.nan
    return float(np.sum(col[valid] * colvol[valid]) / np.sum(colvol[valid]))

def level_thickness(depth):
    depth = np.asarray(depth, dtype=float)
    if depth.size == 1:
        return np.array([1.0])
    bounds = np.empty(depth.size + 1, dtype=float)
    bounds[1:-1] = 0.5 * (depth[:-1] + depth[1:])
    bounds[0] = max(0.0, depth[0] - 0.5 * (depth[1] - depth[0]))
    bounds[-1] = depth[-1] + 0.5 * (depth[-1] - depth[-2])
    return np.diff(bounds)

def load_glodap_var(varfile):
    path = os.path.join(GLODAP_ROOT, f"GLODAPv2.2016b.{varfile}.nc")
    ds = xr.open_dataset(path)
    if varfile in ds.data_vars:
        name = varfile
    else:
        name = None
        for v in ds.data_vars:
            if v.lower() == varfile.lower():
                name = v
                break
        if name is None:
            for v in ds.data_vars:
                if "depth_surface" in ds[v].dims or "depth" in ds[v].dims:
                    name = v
                    break
        if name is None:
            name = list(ds.data_vars)[0]
    da = ds[name]
    depth_dim = "depth_surface" if "depth_surface" in da.dims else da.dims[0]
    y_name = "lat" if "lat" in da.dims else ("latitude" if "latitude" in da.dims else da.dims[1])
    x_name = "lon" if "lon" in da.dims else ("longitude" if "longitude" in da.dims else da.dims[2])
    da = da.transpose(depth_dim, y_name, x_name)
    if "Depth" in ds:
        depth = ds["Depth"].values.astype(float)
    elif depth_dim in ds:
        depth = ds[depth_dim].values.astype(float)
    else:
        depth = ds.coords[depth_dim].values.astype(float)
    lat = ds[y_name].values.astype(float)
    lon = ds[x_name].values.astype(float)

    order_d = np.argsort(depth)
    depth = depth[order_d]
    da = da.isel({depth_dim: order_d})

    lon360 = lon % 360.0
    order_lon = np.argsort(lon360)
    lon_sorted = lon360[order_lon]
    da = da.isel({x_name: order_lon})
    return da.values.astype(float), depth, lat, lon_sorted

def column_mean_0_50_glodap(var, depth, dz, area2d, max_depth=50.0):
    sel = (depth >= 0.0) & (depth <= max_depth)
    v = var[sel]
    w = dz[sel][:, None, None] * area2d[None, :, :]
    valid = np.isfinite(v) & np.isfinite(w) & (w > 0.0)
    num = np.sum(np.where(valid, v * w, 0.0), axis=0)
    den = np.sum(np.where(valid, w, 0.0), axis=0)
    out = np.full(num.shape, np.nan, dtype=float)
    ok = den > 0.0
    out[ok] = num[ok] / den[ok]
    return out, den

def make_bounds_1d(vals):
    vals = np.asarray(vals, dtype=float)
    b = np.empty(vals.size + 1, dtype=float)
    if vals.size == 1:
        b[0] = vals[0] - 0.5
        b[1] = vals[0] + 0.5
    else:
        b[1:-1] = 0.5 * (vals[:-1] + vals[1:])
        b[0] = vals[0] - 0.5 * (vals[1] - vals[0])
        b[-1] = vals[-1] + 0.5 * (vals[-1] - vals[-2])
    return b

def regrid_conservative(field, lat2d, lon2d, target_lat, target_lon):
    """First-order conservative remap of a 2-D field to a regular lat-lon grid.

    Source cells are approximated as lon/lat boxes from neighbouring cell centres.
    Overlap area is computed on the sphere and normalised by the covered area.
    """
    field = np.asarray(field, dtype=float)
    lat2d = np.asarray(lat2d, dtype=float)
    lon2d = np.asarray(lon2d, dtype=float)
    ny, nx = field.shape

    lon = lon2d % 360.0
    if nx > 1:
        dx = np.diff(lon, axis=1)
        dx = ((dx + 180.0) % 360.0) - 180.0
        left = np.empty((ny, nx), dtype=float)
        right = np.empty((ny, nx), dtype=float)
        left[:, 0] = lon[:, 0] - 0.5 * np.abs(dx[:, 0])
        left[:, 1:] = lon[:, 1:] - 0.5 * dx
        right[:, :-1] = lon[:, :-1] + 0.5 * dx
        right[:, -1] = lon[:, -1] + 0.5 * np.abs(dx[:, -1])
    else:
        left = lon - 0.5
        right = lon + 0.5
    lon_min = np.minimum(left, right)
    lon_max = np.maximum(left, right)

    if ny > 1:
        dy = np.diff(lat2d, axis=0)
        bottom = np.empty((ny, nx), dtype=float)
        top = np.empty((ny, nx), dtype=float)
        bottom[0, :] = lat2d[0, :] - 0.5 * np.abs(dy[0, :])
        bottom[1:, :] = 0.5 * (lat2d[:-1, :] + lat2d[1:, :])
        top[:-1, :] = 0.5 * (lat2d[:-1, :] + lat2d[1:, :])
        top[-1, :] = lat2d[-1, :] + 0.5 * np.abs(dy[-1, :])
    else:
        bottom = lat2d - 0.5
        top = lat2d + 0.5

    lat_min = np.minimum(bottom, top)
    lat_max = np.maximum(bottom, top)

    nlat = len(target_lat)
    nlon = len(target_lon)
    lat_b = make_bounds_1d(target_lat)
    lon_b = np.clip(make_bounds_1d(target_lon), 0.0, 360.0)

    num = np.zeros((nlat, nlon), dtype=float)
    den = np.zeros((nlat, nlon), dtype=float)

    for j in range(ny):
        for i in range(nx):
            val = field[j, i]
            if not np.isfinite(val):
                continue

            lo = lon_min[j, i]
            hi = lon_max[j, i]
            width = hi - lo
            if width <= 0.0:
                continue

            intervals = []
            if width >= 360.0:
                intervals.append((0.0, 360.0))
            else:
                lo0 = lo % 360.0
                hi0 = lo0 + width
                if hi0 <= 360.0:
                    intervals.append((lo0, hi0))
                else:
                    intervals.append((lo0, 360.0))
                    intervals.append((0.0, hi0 - 360.0))

            la0 = max(lat_min[j, i], -90.0)
            la1 = min(lat_max[j, i], 90.0)
            if la1 <= la0:
                continue

            i0 = np.searchsorted(lat_b, la0, side="right") - 1
            i1 = np.searchsorted(lat_b, la1, side="left") - 1
            i0 = max(0, min(nlat - 1, i0))
            i1 = max(0, min(nlat - 1, i1))

            for ii in range(i0, i1 + 1):
                olat0 = max(la0, lat_b[ii])
                olat1 = min(la1, lat_b[ii + 1])
                if olat1 <= olat0:
                    continue
                lat_fac = np.sin(np.deg2rad(olat1)) - np.sin(np.deg2rad(olat0))

                for l0, l1 in intervals:
                    if l1 <= l0:
                        continue
                    j0 = np.searchsorted(lon_b, l0, side="right") - 1
                    j1 = np.searchsorted(lon_b, l1, side="left") - 1
                    j0 = max(0, min(nlon - 1, j0))
                    j1 = max(0, min(nlon - 1, j1))

                    for jj in range(j0, j1 + 1):
                        olon0 = max(l0, lon_b[jj])
                        olon1 = min(l1, lon_b[jj + 1])
                        if olon1 <= olon0:
                            continue
                        area = (olon1 - olon0) * lat_fac
                        num[ii, jj] += val * area
                        den[ii, jj] += area

    out = np.full((nlat, nlon), np.nan, dtype=float)
    ok = den > 0.0
    out[ok] = num[ok] / den[ok]
    return out

def vertical_interp_np(field, lev_m, target_depths):
    nlev, ny, nx = field.shape
    f2 = field.reshape(nlev, -1)
    order = np.argsort(lev_m)
    lev = lev_m[order]
    f2 = f2[order, :]
    func = interp1d(
        lev, f2, axis=0, kind="linear",
        bounds_error=False, fill_value="extrapolate", assume_sorted=True
    )
    out = func(np.asarray(target_depths))
    return out.reshape(len(target_depths), ny, nx)

def weighted_bias_rmse_corr(model, ref, area2d, mask=None):
    if mask is not None:
        model = np.where(mask, model, np.nan)
        ref = np.where(mask, ref, np.nan)
    valid = np.isfinite(model) & np.isfinite(ref) & np.isfinite(area2d)
    if not np.any(valid):
        return np.nan, np.nan, np.nan, 0
    w = area2d[valid]
    x = model[valid]
    y = ref[valid]
    diff = x - y
    mean_bias = float(np.sum(diff * w) / np.sum(w))
    rmse = float(np.sqrt(np.sum((diff ** 2) * w) / np.sum(w)))
    xm = np.sum(x * w) / np.sum(w)
    ym = np.sum(y * w) / np.sum(w)
    cov = np.sum(w * (x - xm) * (y - ym))
    denom = np.sqrt(np.sum(w * (x - xm) ** 2) * np.sum(w * (y - ym) ** 2))
    corr = float(cov / denom) if denom > 0 else np.nan
    return mean_bias, rmse, corr, int(np.sum(valid))

# ------------------------------------------------------------------
# GLODAPv2.2016b mapped climatology
# ------------------------------------------------------------------
glodap_dic, glodap_depth, glodap_lat, glodap_lon = load_glodap_var("TCO2")
glodap_talk, _, _, _ = load_glodap_var("TAlk")
glodap_ph, _, _, _ = load_glodap_var("pHtsinsitutp")

dz_glodap = level_thickness(glodap_depth)
area_glodap = np.cos(np.deg2rad(glodap_lat))[:, None] * np.ones((1, len(glodap_lon)))

glodap_dic_col, glodap_dic_colvol = column_mean_0_50_glodap(
    glodap_dic, glodap_depth, dz_glodap, area_glodap
)
glodap_talk_col, glodap_talk_colvol = column_mean_0_50_glodap(
    glodap_talk, glodap_depth, dz_glodap, area_glodap
)
glodap_ph_col, glodap_ph_colvol = column_mean_0_50_glodap(
    glodap_ph, glodap_depth, dz_glodap, area_glodap
)

glodap_dic_global = weighted_global(glodap_dic_col, glodap_dic_colvol)
glodap_talk_global = weighted_global(glodap_talk_col, glodap_talk_colvol)
glodap_ph_global = weighted_global(glodap_ph_col, glodap_ph_colvol)

glodap_map_mask = np.isfinite(glodap_dic_col)

glodap_dic_profile = []
for d in depths:
    idx = int(np.argmin(np.abs(glodap_depth - d)))
    glodap_dic_profile.append(glodap_dic[idx])
glodap_dic_profile = np.asarray(glodap_dic_profile)

# ------------------------------------------------------------------
# CMIP6 models
# ------------------------------------------------------------------
results = {
    "dic_global": {},
    "dic_diff": {},
    "talk_global": {},
    "talk_diff": {},
    "ph_global": {},
    "ph_diff": {},
    "map_bias": {},
    "map_rmse": {},
    "map_corr": {},
    "map_cells": {},
    "profile_bias": {},
    "period_change": {},
    "ratio": {},
}

bias_maps = {}
profile_bias = {}

for model in models:
    mi = model_info[model]
    print(f"Processing {model} ...")

    # --- load 2002 annual means ---
    ds_dic = open_cmip_dataset(model, "dissic")
    ds_theta = open_cmip_dataset(model, "thetao")
    ds_so = open_cmip_dataset(model, "so")
    ds_talk = open_cmip_dataset(model, "talk")
    ds_ph = open_cmip_dataset(model, "ph")
    ds_vol = open_cmip_dataset(model, "volcello")
    ds_area = open_cmip_dataset(model, "areacello")

    if ds_dic is None:
        raise FileNotFoundError(f"dissic not found for {model}")
    if ds_theta is None or ds_so is None:
        raise FileNotFoundError(f"thetao/so not found for {model}")

    dic2002_da = annual_mean(ds_dic["dissic"], 2002)
    theta2002_da = annual_mean(ds_theta["thetao"], 2002)
    so2002_da = annual_mean(ds_so["so"], 2002)
    talk2002_da = annual_mean(ds_talk["talk"], 2002) if ds_talk is not None else None
    ph2002_da = annual_mean(ds_ph["ph"], 2002) if ds_ph is not None else None
    vol2002_da = annual_mean(ds_vol["volcello"], 2002) if ds_vol is not None else None

    lat2d, lon2d = get_lat_lon(ds_dic, mi)
    lev_da = ds_dic["lev"] if "lev" in ds_dic else ds_dic.coords["lev"]
    lev_m = lev_da.values.astype(float) * mi["depth_scale"]

    dic2002 = to_3d_np(dic2002_da, mi)
    theta2002 = to_3d_np(theta2002_da, mi)
    so2002 = to_3d_np(so2002_da, mi)

    if vol2002_da is None:
        raise FileNotFoundError(f"volcello not found for {model}")
    vol2002 = to_3d_np(vol2002_da, mi)

    rho2002 = compute_density(theta2002, so2002, lev_m, lat2d, lon2d)
    dic2002_umol = dic2002 * 1.0e6 / rho2002

    dic2002_col, dic2002_colvol = column_mean_0_50(dic2002_umol, vol2002, lev_m)
    dic2002_global = weighted_global(dic2002_col, dic2002_colvol)

    if talk2002_da is not None:
        talk2002 = to_3d_np(talk2002_da, mi)
        talk2002_umol = talk2002 * 1.0e6 / rho2002
        talk2002_col, talk2002_colvol = column_mean_0_50(talk2002_umol, vol2002, lev_m)
        talk2002_global = weighted_global(talk2002_col, talk2002_colvol)
    else:
        talk2002_global = np.nan
        talk2002_col = np.full(dic2002_col.shape, np.nan)

    if ph2002_da is not None:
        ph2002 = to_3d_np(ph2002_da, mi)
        ph2002_col, ph2002_colvol = column_mean_0_50(ph2002, vol2002, lev_m)
        ph2002_global = weighted_global(ph2002_col, ph2002_colvol)
    else:
        ph2002_global = np.nan
        ph2002_col = np.full(dic2002_col.shape, np.nan)

    # --- 1971-2000 period mean DIC ---
    dic_period_da = period_mean(ds_dic["dissic"], 1971, 2000)
    theta_period_da = period_mean(ds_theta["thetao"], 1971, 2000)
    so_period_da = period_mean(ds_so["so"], 1971, 2000)
    vol_period_da = period_mean(ds_vol["volcello"], 1971, 2000) if ds_vol is not None else None

    dic_period = to_3d_np(dic_period_da, mi)
    theta_period = to_3d_np(theta_period_da, mi)
    so_period = to_3d_np(so_period_da, mi)
    vol_period = to_3d_np(vol_period_da, mi)

    rho_period = compute_density(theta_period, so_period, lev_m, lat2d, lon2d)
    dic_period_umol = dic_period * 1.0e6 / rho_period
    dic_period_col, dic_period_colvol = column_mean_0_50(dic_period_umol, vol_period, lev_m)
    dic_period_global = weighted_global(dic_period_col, dic_period_colvol)

    period_change = dic2002_global - dic_period_global
    dic_diff = dic2002_global - glodap_dic_global
    ratio = period_change / dic_diff if dic_diff != 0 else np.nan

    # --- map comparison ---
    dic_model_regrid = regrid_conservative(
        dic2002_col, lat2d, lon2d, glodap_lat, glodap_lon
    )
    map_bias, map_rmse, map_corr, map_cells = weighted_bias_rmse_corr(
        dic_model_regrid, glodap_dic_col, area_glodap, mask=glodap_map_mask
    )
    bias_map = np.where(
        glodap_map_mask & np.isfinite(dic_model_regrid),
        dic_model_regrid - glodap_dic_col,
        np.nan,
    )

    # --- vertical profile comparison ---
    dic_prof = vertical_interp_np(dic2002_umol, lev_m, depths)
    prof_bias = []
    for k, d in enumerate(depths):
        model_regrid_k = regrid_conservative(
            dic_prof[k], lat2d, lon2d, glodap_lat, glodap_lon
        )
        ref_k = glodap_dic_profile[k]
        valid = np.isfinite(model_regrid_k) & np.isfinite(ref_k)
        if np.any(valid):
            w = area_glodap[valid]
            bias = model_regrid_k[valid] - ref_k[valid]
            prof_bias.append(float(np.sum(bias * w) / np.sum(w)))
        else:
            prof_bias.append(np.nan)
    prof_bias = np.asarray(prof_bias)

    # --- store ---
    results["dic_global"][model] = dic2002_global
    results["dic_diff"][model] = dic_diff
    results["talk_global"][model] = talk2002_global
    results["talk_diff"][model] = talk2002_global - glodap_talk_global
    results["ph_global"][model] = ph2002_global
    results["ph_diff"][model] = ph2002_global - glodap_ph_global
    results["map_bias"][model] = map_bias
    results["map_rmse"][model] = map_rmse
    results["map_corr"][model] = map_corr
    results["map_cells"][model] = map_cells
    results["profile_bias"][model] = prof_bias
    results["period_change"][model] = period_change
    results["ratio"][model] = ratio

    bias_maps[model] = bias_map
    profile_bias[model] = prof_bias

    print(
        f"  DIC 0-50 m: {dic2002_global:.2f} umol/kg, "
        f"diff vs GLODAP: {dic_diff:.2f}"
    )
    print(
        f"  TAlk 0-50 m: {talk2002_global:.2f} umol/kg, "
        f"diff: {talk2002_global - glodap_talk_global:.2f}"
    )
    print(
        f"  pH 0-50 m: {ph2002_global:.4f}, "
        f"diff: {ph2002_global - glodap_ph_global:.4f}"
    )
    print(
        f"  map bias/RMSE/corr: {map_bias:.2f} / {map_rmse:.2f} / {map_corr:.3f} "
        f"({map_cells} cells)"
    )
    print(
        f"  period change: {period_change:.2f} umol/kg, "
        f"ratio to offset: {ratio:.2f}"
    )

# ------------------------------------------------------------------
# Summary tables / best model / sign consistency
# ------------------------------------------------------------------
summary = pd.DataFrame(
    {
        "DIC_global": [results["dic_global"][m] for m in models],
        "DIC_diff": [results["dic_diff"][m] for m in models],
        "TAlk_global": [results["talk_global"][m] for m in models],
        "TAlk_diff": [results["talk_diff"][m] for m in models],
        "pH_global": [results["ph_global"][m] for m in models],
        "pH_diff": [results["ph_diff"][m] for m in models],
        "map_bias": [results["map_bias"][m] for m in models],
        "map_rmse": [results["map_rmse"][m] for m in models],
        "map_corr": [results["map_corr"][m] for m in models],
        "period_change": [results["period_change"][m] for m in models],
        "ratio": [results["ratio"][m] for m in models],
    },
    index=models,
)

score = {}
for m in models:
    score[m] = (
        abs(results["dic_diff"][m])
        + abs(results["talk_diff"][m])
        + 1000.0 * abs(results["ph_diff"][m])
        + results["map_rmse"][m]
        + np.nansum(np.abs(results["profile_bias"][m]))
    )
best_model = min(models, key=lambda m: score[m])

signs = {
    m: (
        np.sign(results["dic_diff"][m]),
        np.sign(results["talk_diff"][m]),
        np.sign(results["ph_diff"][m]),
    )
    for m in models
}

print("\nSummary (unrounded):")
print(summary)
print(f"\nBest model by combined error score: {best_model}")
print("Signs of (DIC, TAlk, pH) model-minus-GLODAP offsets:")
for m in models:
    print(f"  {m}: {signs[m]}")
print(f"\nGLODAP 0-50 m: DIC={glodap_dic_global:.2f}, TAlk={glodap_talk_global:.2f}, pH={glodap_ph_global:.4f}")
print(f"Comparison-mask cells: {int(np.sum(glodap_map_mask))}")

# Make key objects available to the visualization block.
dic_global = results["dic_global"]
dic_diff = results["dic_diff"]
talk_global = results["talk_global"]
talk_diff = results["talk_diff"]
ph_global = results["ph_global"]
ph_diff = results["ph_diff"]
map_bias_values = results["map_bias"]
map_rmse_values = results["map_rmse"]
map_corr_values = results["map_corr"]
map_cells = results["map_cells"]
period_change_values = results["period_change"]
ratio_values = results["ratio"]
