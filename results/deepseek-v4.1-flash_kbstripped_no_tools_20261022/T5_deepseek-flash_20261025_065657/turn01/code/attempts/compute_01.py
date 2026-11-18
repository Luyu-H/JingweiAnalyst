import os
import glob
import re
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy import stats

# ---------------------------------------------------------------------------
# Paths and model metadata
# ---------------------------------------------------------------------------
HIST_ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"
SSP_ROOT = "/data1/lorraine-20419/benchmark_cmip6/ssp245"
WOA_ROOT = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
WOA_MASKS = "/data1/lorraine-20419/ocean_data/WOA23/masks"

MODELS = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

MODEL_INFO = {
    "CESM2": dict(
        y_dim="nlat", x_dim="nlon", lat_coord="lat", lon_coord="lon",
        depth_coord="lev", depth_scale=0.01,
    ),
    "GFDL-ESM4": dict(
        y_dim="lat", x_dim="lon", lat_coord="lat", lon_coord="lon",
        depth_coord="lev", depth_scale=1.0,
    ),
    "MPI-ESM1-2-LR": dict(
        y_dim="j", x_dim="i", lat_coord="latitude", lon_coord="longitude",
        depth_coord="lev", depth_scale=1.0,
    ),
}

COLORS = {
    "CESM2": "tab:blue",
    "GFDL-ESM4": "tab:orange",
    "MPI-ESM1-2-LR": "tab:green",
}

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def read_numeric_table(path):
    """Read whitespace-separated numeric rows, skipping comments/headers."""
    rows = []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#") or s.startswith("%"):
                continue
            parts = s.split()
            try:
                vals = [float(p) for p in parts]
            except ValueError:
                continue
            rows.append(vals)
    if len(rows) == 0:
        return np.empty((0, 0))
    maxlen = max(len(r) for r in rows)
    out = np.full((len(rows), maxlen), np.nan, dtype=float)
    for i, r in enumerate(rows):
        out[i, :len(r)] = r
    return out


def regular_area(lat, lon, earth_radius=6371000.0):
    """Area of regular latitude-longitude cells, in m2."""
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)

    def edges(v):
        e = np.empty(v.size + 1, dtype=float)
        e[1:-1] = 0.5 * (v[:-1] + v[1:])
        e[0] = v[0] - 0.5 * (v[1] - v[0])
        e[-1] = v[-1] + 0.5 * (v[-1] - v[-2])
        return e

    lat_e = edges(lat)
    lon_e = edges(lon)
    dsin = np.sin(np.deg2rad(lat_e[1:])) - np.sin(np.deg2rad(lat_e[:-1]))
    dlon = np.deg2rad(np.diff(lon_e))
    return earth_radius**2 * dsin[:, None] * dlon[None, :]


def weighted_mean(da, w, dims=None):
    """NaN-aware weighted mean of an xarray DataArray."""
    if dims is None:
        dims = da.dims
    valid = np.isfinite(da) & np.isfinite(w) & (w > 0)
    den = w.where(valid).sum(dim=dims, skipna=True)
    num = (da * w).where(valid).sum(dim=dims, skipna=True)
    return num / den


def area_weighted_global_mean(da, area, ydim, xdim):
    """Monthly or static area-weighted global mean."""
    valid = np.isfinite(da) & np.isfinite(area) & (area > 0)
    num = (da * area).where(valid).sum(dim=(ydim, xdim), skipna=True)
    den = area.where(valid).sum(dim=(ydim, xdim), skipna=True)
    return num / den


def year_mean(da, y0, y1):
    """Mean over all times whose year is in [y0, y1]."""
    if "time" not in da.dims:
        return da
    m = (da.time.dt.year >= y0) & (da.time.dt.year <= y1)
    return da.where(m, drop=True).mean("time")


def subset_time(da, y0, y1):
    if "time" not in da.dims:
        return da
    m = (da.time.dt.year >= y0) & (da.time.dt.year <= y1)
    return da.where(m, drop=True)


def parse_variant(path):
    m = re.search(r"(r\d+i\d+p\d+f\d+)", os.path.basename(path))
    return m.group(1) if m else None


def pick_data_var(ds, preferred):
    if preferred in ds.data_vars:
        return ds[preferred]
    # Prefer a variable whose name contains the preferred name.
    for v in ds.data_vars:
        if preferred.lower() in v.lower() and ds[v].ndim >= 2:
            return ds[v]
    # Otherwise take the first substantial data variable.
    for v in ds.data_vars:
        if ds[v].ndim >= 2:
            return ds[v]
    raise KeyError(f"No suitable data variable for {preferred!r} in {list(ds.data_vars)}")


def load_cmip_var(root, model, var, fallback_root=None):
    """Open all files for root/model/var, with optional fallback root."""
    files = []
    for r in [root, fallback_root]:
        if r is None:
            continue
        d = os.path.join(r, model, var)
        fs = sorted(glob.glob(os.path.join(d, "*.nc")))
        if fs:
            files = fs
            break
    if not files:
        return None, None
    ds = xr.open_mfdataset(
        files,
        combine="by_coords",
        use_cftime=True,
        decode_times=True,
        chunks="auto",
        parallel=True,
    )
    da = pick_data_var(ds, var)
    return da, files


def get_latlon_2d(da, info):
    """Return 2-D lat and lon arrays with shape (y, x)."""
    lat_name = info["lat_coord"]
    lon_name = info["lon_coord"]
    lat = da[lat_name]
    lon = da[lon_name]
    ny = da.sizes[info["y_dim"]]
    nx = da.sizes[info["x_dim"]]

    if lat.ndim == 1 and lon.ndim == 1:
        lon2d, lat2d = np.meshgrid(np.asarray(lon.values), np.asarray(lat.values))
    else:
        lat2d = np.asarray(lat.values)
        lon2d = np.asarray(lon.values)
        if lat2d.shape == (nx, ny) and lon2d.shape == (nx, ny):
            lat2d = lat2d.T
            lon2d = lon2d.T
    return lat2d, lon2d


def compute_sss(so, vol, info, ocean_mask=None):
    """0-10 m SSS, weighted by full level volcello for selected levels."""
    lev_name = info["depth_coord"]
    scale = info["depth_scale"]

    lev_vals = np.asarray(so[lev_name].values, dtype=float)
    depth_m = lev_vals * scale
    sel_vals = lev_vals[(depth_m >= 0.0) & (depth_m <= 10.0)]
    if sel_vals.size == 0:
        # Defensive fallback: nearest single level to 5 m.
        sel_vals = np.array([lev_vals[np.argmin(np.abs(depth_m - 5.0))]])

    so_sel = so.sel({lev_name: sel_vals})
    vol_sel = vol.sel({lev_name: sel_vals})

    if "time" in vol_sel.dims and "time" in so_sel.dims:
        vol_sel = vol_sel.reindex(time=so_sel.time, method="nearest")

    valid = np.isfinite(so_sel) & np.isfinite(vol_sel) & (vol_sel > 0)
    num = (so_sel.where(valid) * vol_sel.where(valid)).sum(dim=lev_name, skipna=True)
    den = vol_sel.where(valid).sum(dim=lev_name, skipna=True)
    sss = num / den
    sss = sss.where(den > 0)
    if ocean_mask is not None:
        sss = sss.where(ocean_mask)
    sss = sss.assign_attrs(
        units="psu",
        long_name="0-10 m sea surface salinity, full-level-volume weighted",
    )
    return sss


def load_and_compute_sss(root, model, info, y0, y1, fallback_root=None):
    so, so_files = load_cmip_var(root, model, "so", fallback_root=fallback_root)
    if so is None:
        raise FileNotFoundError(f"so not found for {model} in {root}")
    so = subset_time(so, y0, y1)

    area, _ = load_cmip_var(root, model, "areacello", fallback_root=fallback_root)
    if area is None and fallback_root is not None:
        area, _ = load_cmip_var(fallback_root, model, "areacello")
    if area is None:
        raise FileNotFoundError(f"areacello not found for {model}")
    if "time" in area.dims:
        area = area.isel(time=0, drop=True)

    sftof, _ = load_cmip_var(root, model, "sftof", fallback_root=fallback_root)
    if sftof is not None:
        if "time" in sftof.dims:
            sftof = sftof.isel(time=0, drop=True)
        ocean_mask = sftof > 0
    else:
        ocean_mask = area > 0

    area = area.where(ocean_mask)

    vol, _ = load_cmip_var(root, model, "volcello", fallback_root=fallback_root)
    if vol is None:
        raise FileNotFoundError(f"volcello not found for {model}")
    if "time" in vol.dims:
        vol = subset_time(vol, y0, y1)
        if "time" in so.dims:
            vol = vol.reindex(time=so.time, method="nearest")

    sss = compute_sss(so, vol, info, ocean_mask=ocean_mask)
    return sss, area, so_files, vol


def regrid_2d_to_regular(da, lat2d, lon2d, target_lat, target_lon, method="linear"):
    """Regrid a 2-D curvilinear/regular field to a regular lon-lat grid."""
    vals = np.asarray(da.values, dtype=float).ravel()
    lats = np.asarray(lat2d, dtype=float).ravel()
    lons = np.asarray(lon2d, dtype=float).ravel() % 360.0

    tlon = np.asarray(target_lon, dtype=float) % 360.0
    tlat = np.asarray(target_lat, dtype=float)
    TLON, TLAT = np.meshgrid(tlon, tlat)

    # Duplicate points across the seam so interpolation is global.
    pts_lon = np.concatenate([lons, lons - 360.0, lons + 360.0])
    pts_lat = np.tile(lats, 3)
    pts_val = np.tile(vals, 3)

    valid = np.isfinite(pts_lon) & np.isfinite(pts_lat) & np.isfinite(pts_val)
    pts = np.column_stack([pts_lon[valid], pts_lat[valid]])
    v = pts_val[valid]

    out = griddata(pts, v, (TLON, TLAT), method=method)
    if method == "linear" and np.any(~np.isfinite(out)):
        out_near = griddata(pts, v, (TLON, TLAT), method="nearest")
        out = np.where(np.isfinite(out), out, out_near)

    return xr.DataArray(
        out,
        coords={"lat": target_lat, "lon": target_lon},
        dims=("lat", "lon"),
        name=getattr(da, "name", None),
    )


def compute_difference_metrics(data, reference, weights=None, mask=None):
    """Area-weighted metrics over the common valid sample."""
    if mask is not None:
        data = data.where(mask)
        reference = reference.where(mask)
        if weights is not None:
            weights = weights.where(mask)

    d = np.asarray(data.values, dtype=float).ravel()
    r = np.asarray(reference.values, dtype=float).ravel()
    if weights is None:
        w = np.ones_like(d, dtype=float)
    else:
        w = np.asarray(weights.values, dtype=float).ravel()

    valid = np.isfinite(d) & np.isfinite(r) & np.isfinite(w) & (w > 0)
    d, r, w = d[valid], r[valid], w[valid]
    if d.size == 0:
        return dict(
            mean_data=np.nan, mean_reference=np.nan, bias=np.nan,
            rmse=np.nan, correlation=np.nan, n=0,
        )

    wsum = w.sum()
    mean_d = np.sum(w * d) / wsum
    mean_r = np.sum(w * r) / wsum
    bias = mean_d - mean_r
    rmse = np.sqrt(np.sum(w * (d - r) ** 2) / wsum)

    vd = d - mean_d
    vr = r - mean_r
    cov = np.sum(w * vd * vr) / wsum
    vard = np.sum(w * vd**2) / wsum
    varr = np.sum(w * vr**2) / wsum
    corr = cov / np.sqrt(vard * varr) if vard > 0 and varr > 0 else np.nan

    return dict(
        mean_data=float(mean_d),
        mean_reference=float(mean_r),
        bias=float(bias),
        rmse=float(rmse),
        correlation=float(corr),
        n=int(d.size),
    )


def fit_1D_trend(y, x=None, method="linear", use_ess=True):
    """Linear trend with optional lag-1 effective-sample-size correction."""
    y = np.asarray(y, dtype=float)
    if x is None:
        x = np.arange(y.size, dtype=float)
    else:
        x = np.asarray(x, dtype=float)

    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid], y[valid]
    n = x.size

    out = dict(
        slope=np.nan, intercept=np.nan, r2=np.nan, p=np.nan,
        ci_low=np.nan, ci_high=np.nan, se=np.nan,
        n=n, n_eff=np.nan,
        slope_uncorrected=np.nan, se_uncorrected=np.nan,
        ci_low_uncorrected=np.nan, ci_high_uncorrected=np.nan,
        p_uncorrected=np.nan,
    )
    if n < 3:
        return out

    xmean = x.mean()
    ymean = y.mean()
    Sxx = np.sum((x - xmean) ** 2)
    if Sxx <= 0:
        return out
    Sxy = np.sum((x - xmean) * (y - ymean))
    slope = Sxy / Sxx
    intercept = ymean - slope * xmean
    resid = y - (intercept + slope * x)
    SSE = np.sum(resid**2)
    SST = np.sum((y - ymean) ** 2)
    r2 = 1.0 - SSE / SST if SST > 0 else np.nan

    if n > 2:
        se_ols = np.sqrt(SSE / (n - 2) / Sxx)
        tcrit_ols = stats.t.ppf(0.975, n - 2)
        p_ols = 2.0 * stats.t.sf(np.abs(slope / se_ols), n - 2) if se_ols > 0 else np.nan
        ci_ols = (slope - tcrit_ols * se_ols, slope + tcrit_ols * se_ols)
    else:
        se_ols = np.nan
        p_ols = np.nan
        ci_ols = (np.nan, np.nan)

    out.update(
        slope_uncorrected=slope,
        se_uncorrected=se_ols,
        p_uncorrected=p_ols,
        ci_low_uncorrected=ci_ols[0],
        ci_high_uncorrected=ci_ols[1],
    )

    if use_ess and n > 3 and np.std(resid) > 0:
        r1 = np.corrcoef(resid[:-1], resid[1:])[0, 1]
        if not np.isfinite(r1):
            r1 = 0.0
        r1 = np.clip(r1, -0.99, 0.99)
        n_eff = n * (1.0 - r1) / (1.0 + r1)
        n_eff = min(float(n), max(3.0, n_eff))
    else:
        n_eff = float(n)

    if n_eff > 2 and Sxx > 0:
        se = np.sqrt(SSE / (n_eff - 2) / Sxx)
        df = n_eff - 2.0
        tcrit = stats.t.ppf(0.975, df)
        tstat = slope / se if se > 0 else np.nan
        pval = 2.0 * stats.t.sf(np.abs(tstat), df) if np.isfinite(tstat) else np.nan
        ci = (slope - tcrit * se, slope + tcrit * se)
    else:
        se = np.nan
        pval = np.nan
        ci = (np.nan, np.nan)

    out.update(
        slope=slope,
        intercept=intercept,
        r2=r2,
        p=pval,
        ci_low=ci[0],
        ci_high=ci[1],
        se=se,
        n_eff=n_eff,
    )
    return out


def cellwise_trend(annual_da, use_ess=True):
    """Cell-by-cell linear trend of an annual field with ESS p-values."""
    y = np.asarray(annual_da.values, dtype=float)
    if "year" in annual_da.coords:
        years = np.asarray(annual_da["year"].values, dtype=float)
    elif "time" in annual_da.coords:
        years = np.asarray(annual_da["time"].dt.year.values, dtype=float)
    else:
        years = np.arange(y.shape[0], dtype=float)

    x = years[:, None, None]
    valid = np.isfinite(y)
    n = valid.sum(axis=0)

    xsum = np.where(valid, x, 0.0).sum(axis=0)
    ysum = np.where(valid, y, 0.0).sum(axis=0)
    xmean = xsum / n
    ymean = ysum / n

    Sxx = np.where(valid, (x - xmean) ** 2, 0.0).sum(axis=0)
    Sxy = np.where(valid, (x - xmean) * (y - ymean), 0.0).sum(axis=0)
    slope = Sxy / Sxx
    intercept = ymean - slope * xmean
    resid = y - (intercept + slope * x)

    SSE = np.where(valid, resid**2, 0.0).sum(axis=0)
    SST = np.where(valid, (y - ymean) ** 2, 0.0).sum(axis=0)
    r2 = 1.0 - SSE / SST

    good = (n > 2) & (Sxx > 0) & np.isfinite(slope)

    if use_ess:
        valid_pair = valid[1:] & valid[:-1]
        r1_num = np.where(valid_pair, resid[1:] * resid[:-1], 0.0).sum(axis=0)
        r1_den1 = np.where(valid_pair, resid[1:] ** 2, 0.0).sum(axis=0)
        r1_den2 = np.where(valid_pair, resid[:-1] ** 2, 0.0).sum(axis=0)
        denom = np.sqrt(r1_den1 * r1_den2)
        r1 = np.divide(r1_num, denom, out=np.zeros_like(r1_num), where=denom > 0)
        r1 = np.where(np.isfinite(r1), r1, 0.0)
        r1 = np.clip(r1, -0.99, 0.99)
        n_eff = n * (1.0 - r1) / (1.0 + r1)
        n_eff = np.minimum(n, np.maximum(3.0, n_eff))
    else:
        n_eff = n.astype(float)

    pval = np.full_like(slope, np.nan, dtype=float)
    good_eff = good & (n_eff > 2)
    se = np.full_like(slope, np.nan, dtype=float)
    se[good_eff] = np.sqrt(SSE[good_eff] / (n_eff[good_eff] - 2.0) / Sxx[good_eff])
    df = n_eff - 2.0
    tstat = np.full_like(slope, np.nan, dtype=float)
    tstat[good_eff] = slope[good_eff] / se[good_eff]
    pval[good_eff] = 2.0 * stats.t.sf(np.abs(tstat[good_eff]), df[good_eff])

    return slope, pval, r2, n_eff


def select_basin(basin_da, basin_name, include_marginal_seas=True):
    """Boolean mask for a WOA basin, optionally folding marginal seas."""
    # WOA13-style basin numbering; extra codes are harmless if absent.
    mapping = {
        "Atlantic": [1, 6, 7, 8, 11],
        "Pacific": [2],
        "Indian": [3, 9, 10],
        "Arctic": [4],
        "Southern": [5],
    }
    ids = mapping.get(basin_name, [])
    if not include_marginal_seas and len(ids) > 1:
        ids = ids[:1]
    return basin_da.isin(ids)


def load_landsea_mask(masks_root, target_lat, target_lon):
    path = os.path.join(masks_root, "landsea_01.msk")
    arr = read_numeric_table(path)
    lat = arr[:, 0]
    lon = arr[:, 1]
    val = arr[:, 2]
    lat_idx = np.abs(lat[:, None] - np.asarray(target_lat)[None, :]).argmin(axis=1)
    lon_idx = np.abs(lon[:, None] - np.asarray(target_lon)[None, :]).argmin(axis=1)
    out = np.full((len(target_lat), len(target_lon)), np.nan, dtype=float)
    out[lat_idx, lon_idx] = val
    return xr.DataArray(out, coords={"lat": target_lat, "lon": target_lon}, dims=("lat", "lon"))


def load_basin_mask(masks_root, target_lat, target_lon):
    path = os.path.join(masks_root, "basinmask_01.msk")
    arr = read_numeric_table(path)
    lat = arr[:, 0]
    lon = arr[:, 1]
    # Third column is the level-0 basin value; sparse rows are ocean points only.
    val = arr[:, 2]
    lat_idx = np.abs(lat[:, None] - np.asarray(target_lat)[None, :]).argmin(axis=1)
    lon_idx = np.abs(lon[:, None] - np.asarray(target_lon)[None, :]).argmin(axis=1)
    out = np.full((len(target_lat), len(target_lon)), np.nan, dtype=float)
    out[lat_idx, lon_idx] = val
    return xr.DataArray(out, coords={"lat": target_lat, "lon": target_lon}, dims=("lat", "lon"))


# ---------------------------------------------------------------------------
# WOA23 salinity, masks, and comparison mask
# ---------------------------------------------------------------------------
woa_file = os.path.join(
    WOA_ROOT,
    "salinity",
    "1971-2000_climate_normal",
    "woa23_decav71A0_s00_01.nc",
)
ds_woa = xr.open_dataset(woa_file, decode_times=False)

# Find the objectively analysed salinity field.
woa_sal = None
for cand in ["s_an", "s", "salinity", "so"]:
    if cand in ds_woa.data_vars:
        woa_sal = ds_woa[cand]
        break
if woa_sal is None:
    woa_sal = pick_data_var(ds_woa, "s")

# Drop the single annual time slice if present.
if "time" in woa_sal.dims:
    woa_sal = woa_sal.isel(time=0, drop=True)

# Standardise coordinate names.
rename = {}
if "latitude" in woa_sal.dims:
    rename["latitude"] = "lat"
if "longitude" in woa_sal.dims:
    rename["longitude"] = "lon"
if "lev" in woa_sal.dims:
    rename["lev"] = "depth"
if rename:
    woa_sal = woa_sal.rename(rename)

woa_lat = np.asarray(ds_woa["lat"].values, dtype=float)
woa_lon = np.asarray(ds_woa["lon"].values, dtype=float)
woa_area_2d = regular_area(woa_lat, woa_lon)
woa_area = xr.DataArray(
    woa_area_2d, coords={"lat": woa_lat, "lon": woa_lon}, dims=("lat", "lon")
)

depth_vals = np.asarray(woa_sal["depth"].values, dtype=float)
depth_abs = np.abs(depth_vals)
target_levels = np.array([0.0, 5.0, 10.0])
sel_idx = np.array([np.argmin(np.abs(depth_abs - d)) for d in target_levels])
sel_idx = np.unique(sel_idx)

# Full level thickness from consecutive standard depths.
order = np.argsort(depth_abs)
depth_sorted = depth_abs[order]
thick_sorted = np.empty_like(depth_sorted)
thick_sorted[:-1] = np.diff(depth_sorted)
thick_sorted[-1] = thick_sorted[-2] if thick_sorted.size > 1 else 1.0
thick = np.empty_like(depth_abs)
thick[order] = thick_sorted
thick_sel = thick[sel_idx]

woa_sel = woa_sal.isel(depth=sel_idx)
vol_woa_np = woa_area_2d[None, :, :] * thick_sel[:, None, None]
vol_woa = xr.DataArray(
    vol_woa_np,
    coords={"depth": woa_sel["depth"], "lat": woa_lat, "lon": woa_lon},
    dims=("depth", "lat", "lon"),
)

valid_woa = np.isfinite(woa_sel) & np.isfinite(vol_woa) & (vol_woa > 0)
num_woa = (woa_sel.where(valid_woa) * vol_woa.where(valid_woa)).sum("depth", skipna=True)
den_woa = vol_woa.where(valid_woa).sum("depth", skipna=True)
woa_sss = (num_woa / den_woa).where(den_woa > 0)
woa_sss = woa_sss.assign_attrs(units="psu", long_name="0-10 m WOA23 salinity")

landsea = load_landsea_mask(WOA_MASKS, woa_lat, woa_lon)
woa_ocean = landsea > 1
woa_sss = woa_sss.where(woa_ocean)

basin_mask = load_basin_mask(WOA_MASKS, woa_lat, woa_lon)

comparison_mask = woa_ocean & np.isfinite(woa_sss)
mask_cell_count = int(comparison_mask.sum().compute().item())
woa_sss_masked = woa_sss.where(comparison_mask)
woa_area_masked = woa_area.where(comparison_mask)

woa_global_mean = float(
    weighted_mean(woa_sss_masked, woa_area_masked, dims=("lat", "lon")).compute().item()
)

# ---------------------------------------------------------------------------
# Model processing
# ---------------------------------------------------------------------------
bias_fields = {}
trend_fields = {}
trend_sig_fields = {}
series_hist = {}
series_ssp = {}
fit_params_hist = {}
fit_params_ssp = {}
model_results = {}

for model in MODELS:
    info = MODEL_INFO[model]
    ydim, xdim = info["y_dim"], info["x_dim"]

    # Historical 1950-2014 and ssp245 2071-2100 monthly SSS.
    sss_hist, area_hist, so_files_hist, vol_hist = load_and_compute_sss(
        HIST_ROOT, model, info, 1950, 2014, fallback_root=None
    )
    sss_ssp, area_ssp, so_files_ssp, vol_ssp = load_and_compute_sss(
        SSP_ROOT, model, info, 2071, 2100, fallback_root=HIST_ROOT
    )

    var_hist = parse_variant(so_files_hist[0]) if so_files_hist else None
    var_ssp = parse_variant(so_files_ssp[0]) if so_files_ssp else None

    # 1971-2000 model climatology and native global mean.
    sss_clim_native = year_mean(sss_hist, 1971, 2000).compute()
    native_global_mean = float(
        area_weighted_global_mean(sss_clim_native, area_hist, ydim, xdim).compute().item()
    )

    lat2d, lon2d = get_latlon_2d(sss_clim_native, info)

    sss_clim_reg = regrid_2d_to_regular(
        sss_clim_native, lat2d, lon2d, woa_lat, woa_lon, method="linear"
    )
    area_reg = regrid_2d_to_regular(
        area_hist, lat2d, lon2d, woa_lat, woa_lon, method="linear"
    )

    sss_clim_masked = sss_clim_reg.where(comparison_mask)
    area_masked = area_reg.where(comparison_mask)

    metrics = compute_difference_metrics(
        sss_clim_masked, woa_sss_masked, weights=area_masked
    )

    bias_fields[model] = (sss_clim_masked - woa_sss_masked).assign_attrs(
        units="psu", long_name="SSS bias vs WOA23 1971-2000"
    )

    # Global mean monthly SSS, then annual means.
    gm_hist = area_weighted_global_mean(sss_hist, area_hist, ydim, xdim).compute()
    gm_ssp = area_weighted_global_mean(sss_ssp, area_ssp, ydim, xdim).compute()

    annual_hist = gm_hist.groupby("time.year").mean("time").compute()
    annual_ssp = gm_ssp.groupby("time.year").mean("time").compute()
    series_hist[model] = annual_hist
    series_ssp[model] = annual_ssp

    years_hist = np.asarray(annual_hist["year"].values, dtype=float)
    years_ssp = np.asarray(annual_ssp["year"].values, dtype=float)

    fit_h = fit_1D_trend(annual_hist.values, x=years_hist, method="linear", use_ess=True)
    fit_s = fit_1D_trend(annual_ssp.values, x=years_ssp, method="linear", use_ess=True)
    fit_params_hist[model] = fit_h
    fit_params_ssp[model] = fit_s

    def trend_summary(fit):
        return dict(
            slope_per_decade=fit["slope"] * 10.0,
            ci_low_per_decade=fit["ci_low"] * 10.0,
            ci_high_per_decade=fit["ci_high"] * 10.0,
            r2=fit["r2"],
            p=fit["p"],
            significant=bool(
                np.isfinite(fit["ci_low"])
                and np.isfinite(fit["ci_high"])
                and (fit["ci_low"] > 0 or fit["ci_high"] < 0)
            ),
            n=fit["n"],
            n_eff=fit["n_eff"],
            slope_uncorrected_per_decade=fit["slope_uncorrected"] * 10.0,
            ci_low_uncorrected_per_decade=fit["ci_low_uncorrected"] * 10.0,
            ci_high_uncorrected_per_decade=fit["ci_high_uncorrected"] * 10.0,
        )

    hist_trend = trend_summary(fit_h)
    ssp_trend = trend_summary(fit_s)

    if np.isfinite(fit_h["slope"]) and abs(fit_h["slope"]) > 0:
        slope_ratio = fit_s["slope"] / fit_h["slope"]
    else:
        slope_ratio = np.nan
    ratio_meaningful = bool(
        np.isfinite(fit_h["ci_low"])
        and np.isfinite(fit_h["ci_high"])
        and (fit_h["ci_low"] > 0 or fit_h["ci_high"] < 0)
    )

    # Cell-by-cell 1950-2014 annual-mean trend.
    annual_sss_hist = sss_hist.groupby("time.year").mean("time").compute()
    slope_cell, pval_cell, r2_cell, neff_cell = cellwise_trend(annual_sss_hist, use_ess=True)
    slope_cell_dec = slope_cell * 10.0
    sig_cell = pval_cell < 0.05

    template = annual_sss_hist.isel(year=0, drop=True)
    trend_da = xr.DataArray(
        slope_cell_dec, coords=template.coords, dims=template.dims,
        name="sss_trend_per_decade",
    )
    sig_da = xr.DataArray(
        sig_cell.astype(float), coords=template.coords, dims=template.dims,
        name="sss_trend_significant",
    )

    trend_reg = regrid_2d_to_regular(
        trend_da, lat2d, lon2d, woa_lat, woa_lon, method="linear"
    ).where(comparison_mask)
    sig_reg = regrid_2d_to_regular(
        sig_da, lat2d, lon2d, woa_lat, woa_lon, method="nearest"
    )
    trend_fields[model] = trend_reg
    trend_sig_fields[model] = (sig_reg > 0.5) & comparison_mask

    # Delta SSS: ssp245 2071-2100 minus historical 1995-2014.
    hist_base = year_mean(sss_hist, 1995, 2014).compute()
    ssp_mean = year_mean(sss_ssp, 2071, 2100).compute()
    delta_native = (ssp_mean - hist_base).assign_attrs(
        units="psu", long_name="ssp245 minus historical 1995-2014 SSS"
    )

    delta_reg = regrid_2d_to_regular(
        delta_native, lat2d, lon2d, woa_lat, woa_lon, method="linear"
    )
    area_reg_hist = regrid_2d_to_regular(
        area_hist, lat2d, lon2d, woa_lat, woa_lon, method="linear"
    )

    valid_delta = np.isfinite(delta_reg) & np.isfinite(area_reg_hist) & woa_ocean
    delta_global = float(
        weighted_mean(
            delta_reg.where(valid_delta),
            area_reg_hist.where(valid_delta),
            dims=("lat", "lon"),
        ).compute().item()
    )

    delta_basins = {}
    for basin_name in ["Atlantic", "Indian", "Pacific"]:
        bmask = select_basin(basin_mask, basin_name, include_marginal_seas=True)
        v = valid_delta & bmask
        delta_basins[basin_name] = float(
            weighted_mean(
                delta_reg.where(v),
                area_reg_hist.where(v),
                dims=("lat", "lon"),
            ).compute().item()
        )

    if var_hist != var_ssp:
        caveat = (
            f"{model}: historical is {var_hist} but ssp245 is {var_ssp}; "
            "delta SSS and the slope ratio include internal-variability differences "
            "between two different ensemble members."
        )
    else:
        caveat = (
            f"{model}: historical {var_hist} and ssp245 {var_ssp} are the same "
            "continuous ensemble member."
        )

    model_results[model] = dict(
        variant_historical=var_hist,
        variant_ssp245=var_ssp,
        member_caveat=caveat,
        native_global_mean_1971_2000_psu=native_global_mean,
        sss_global_mean_psu=metrics["mean_data"],
        woa_mean_over_model_mask_psu=metrics["mean_reference"],
        bias_psu=metrics["bias"],
        rmse_psu=metrics["rmse"],
        correlation=metrics["correlation"],
        comparison_cells=metrics["n"],
        historical_trend=hist_trend,
        ssp245_trend=ssp_trend,
        slope_ratio_ssp_to_hist=slope_ratio,
        slope_ratio_1dp=round(float(slope_ratio), 1) if np.isfinite(slope_ratio) else np.nan,
        ratio_meaningful=ratio_meaningful,
        delta_global_psu=delta_global,
        delta_basins_psu=delta_basins,
    )

# ---------------------------------------------------------------------------
# Compact report
# ---------------------------------------------------------------------------
results = {
    "woa_global_mean_1971_2000_psu": woa_global_mean,
    "comparison_mask_cells": mask_cell_count,
    "models": model_results,
}

print("WOA23 1971-2000 global mean SSS (psu):", round(woa_global_mean, 2))
print("Common comparison-mask cells:", mask_cell_count)
print()
for model in MODELS:
    r = model_results[model]
    print(model)
    print("  variants:", r["variant_historical"], r["variant_ssp245"])
    print("  model 1971-2000 mean SSS (psu):", round(r["sss_global_mean_psu"], 2))
    print("  WOA mean over same mask (psu):", round(r["woa_mean_over_model_mask_psu"], 2))
    print("  bias model-WOA (psu):", round(r["bias_psu"], 2))
    print("  RMSE (psu):", round(r["rmse_psu"], 2))
    print("  spatial correlation:", round(r["correlation"], 3))
    print("  historical 1950-2014 trend (psu/decade): "
          f"{r['historical_trend']['slope_per_decade']:.3g} "
          f"[{r['historical_trend']['ci_low_per_decade']:.3g}, "
          f"{r['historical_trend']['ci_high_per_decade']:.3g}], "
          f"R2={r['historical_trend']['r2']:.3g}, "
          f"significant={r['historical_trend']['significant']}")
    print("  ssp245 2071-2100 trend (psu/decade): "
          f"{r['ssp245_trend']['slope_per_decade']:.3g} "
          f"[{r['ssp245_trend']['ci_low_per_decade']:.3g}, "
          f"{r['ssp245_trend']['ci_high_per_decade']:.3g}], "
          f"R2={r['ssp245_trend']['r2']:.3g}, "
          f"significant={r['ssp245_trend']['significant']}")
    print("  ssp/historical slope ratio (1 dp):", r["slope_ratio_1dp"],
          "meaningful:", r["ratio_meaningful"])
    print("  delta SSS global (psu):", round(r["delta_global_psu"], 2))
    print("  delta SSS basins (psu):",
          {k: round(v, 2) for k, v in r["delta_basins_psu"].items()})
    print("  caveat:", r["member_caveat"])
    print()
