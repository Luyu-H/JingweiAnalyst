# Computation code only.
# Open raw CMIP6 and WOA23 files directly and compute the requested SSS diagnostics.

import os
import glob
import re
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from scipy import stats as st
from scipy.spatial import cKDTree

# -----------------------------------------------------------------------------
# Paths and model metadata
# -----------------------------------------------------------------------------
woa_root = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
woa_masks_root = "/data1/lorraine-20419/ocean_data/WOA23/masks"
hist_root = "/data1/lorraine-20419/benchmark_cmip6/historical"
ssp_root = "/data1/lorraine-20419/benchmark_cmip6/ssp245"

models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

model_info = {
    "CESM2": dict(
        y_dim="nlat", x_dim="nlon", lat_coord="lat", lon_coord="lon",
        depth_coord="lev", depth_scale=0.01, ocean_mask_from="sftof"
    ),
    "GFDL-ESM4": dict(
        y_dim="lat", x_dim="lon", lat_coord="lat", lon_coord="lon",
        depth_coord="lev", depth_scale=1.0
    ),
    "MPI-ESM1-2-LR": dict(
        y_dim="j", x_dim="i", lat_coord="latitude", lon_coord="longitude",
        depth_coord="lev", depth_scale=1.0
    ),
}

# -----------------------------------------------------------------------------
# Small generic helpers
# -----------------------------------------------------------------------------
def find_files(root, model, var):
    d = os.path.join(root, model, var)
    return sorted(glob.glob(os.path.join(d, "*.nc")))


def open_var(root, model, var, y0=None, y1=None, fallback_root=None, first_only=False):
    files = find_files(root, model, var)
    if len(files) == 0 and fallback_root is not None:
        files = find_files(fallback_root, model, var)
    if len(files) == 0:
        raise FileNotFoundError(f"No files found for {model}/{var} under {root}")
    if first_only and len(files) > 1:
        files = files[:1]
    if len(files) == 1:
        ds = xr.open_dataset(files[0], use_cftime=True, decode_times=True)
    else:
        ds = xr.open_mfdataset(files, combine="by_coords", use_cftime=True, decode_times=True)
    if "time" in ds.dims and y0 is not None and y1 is not None:
        yrs = ds.time.dt.year
        m = (yrs >= y0) & (yrs <= y1)
        ds = ds.isel(time=np.where(m.values)[0])
    return ds, files


def get_variant(files):
    for f in files:
        m = re.search(r"_(r\d+i\d+p\d+f\d+)_", os.path.basename(f))
        if m is not None:
            return m.group(1)
    return "unknown"


def read_numeric_table(path, min_cols=3):
    rows = []
    with open(path, "r", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.replace(",", " ").split()
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue
            if len(vals) >= min_cols:
                rows.append(vals)
    return rows


def read_woa_landsea(path):
    rows = read_numeric_table(path, 3)
    arr = np.array([r[:3] for r in rows], dtype=float)
    mask = np.ones((180, 360), dtype=float)
    lat, lon, val = arr[:, 0], arr[:, 1], arr[:, 2]
    iy = np.rint(lat + 89.5).astype(int)
    ix = np.rint(lon + 179.5).astype(int)
    good = (iy >= 0) & (iy < 180) & (ix >= 0) & (ix < 360)
    mask[iy[good], ix[good]] = val[good]
    return mask


def read_woa_basin(path):
    rows = read_numeric_table(path, 3)
    basin = np.zeros((180, 360), dtype=float)
    for r in rows:
        lat, lon = r[0], r[1]
        iy = int(np.rint(lat + 89.5))
        ix = int(np.rint(lon + 179.5))
        if not (0 <= iy < 180 and 0 <= ix < 360):
            continue
        surf = 0.0
        for v in r[2:]:
            if np.isfinite(v) and v > 0:
                surf = v
                break
        basin[iy, ix] = surf
    return basin


def fold_basin_array(basin, lat, lon):
    """Fold marginal-sea labels (>5) onto the nearest major ocean label."""
    folded = basin.copy()
    major_codes = [1, 2, 3]
    major_mask = np.isin(basin, major_codes)
    if not major_mask.any():
        return folded
    lat2d, lon2d = np.meshgrid(lat, lon, indexing="ij")
    latr = np.deg2rad(lat2d)
    lonr = np.deg2rad(lon2d % 360.0)
    xyz = np.stack(
        [np.cos(latr) * np.cos(lonr),
         np.cos(latr) * np.sin(lonr),
         np.sin(latr)],
        axis=-1,
    )
    tree = cKDTree(xyz[major_mask])
    major_vals = basin[major_mask]
    for code in np.unique(basin):
        if code <= 0 or code in [1, 2, 3, 4, 5]:
            continue
        m = basin == code
        if not m.any():
            continue
        dist, idx = tree.query(xyz[m])
        folded[m] = major_vals[idx]
    return folded


def regular_grid_area(lat, lon):
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    lat_edges = np.empty(len(lat) + 1)
    lat_edges[1:-1] = 0.5 * (lat[:-1] + lat[1:])
    lat_edges[0] = lat[0] - 0.5 * (lat[1] - lat[0])
    lat_edges[-1] = lat[-1] + 0.5 * (lat[-1] - lat[-2])
    lon_edges = np.empty(len(lon) + 1)
    lon_edges[1:-1] = 0.5 * (lon[:-1] + lon[1:])
    lon_edges[0] = lon[0] - 0.5 * (lon[1] - lon[0])
    lon_edges[-1] = lon[-1] + 0.5 * (lon[-1] - lon[-2])
    R = 6371000.0
    dlon = np.deg2rad(np.diff(lon_edges))
    dsin = np.diff(np.sin(np.deg2rad(lat_edges)))
    return R**2 * dsin[:, None] * dlon[None, :]


def regrid_conservative(field, lat2d, lon2d, area, woa_lat, woa_lon):
    """Area-weighted conservative binning to the WOA 1-degree grid."""
    f = np.asarray(field, dtype=float)
    a = np.asarray(area, dtype=float)
    lat = np.asarray(lat2d, dtype=float)
    lon = np.asarray(lon2d, dtype=float)
    lon = ((lon + 180.0) % 360.0) - 180.0

    ny = len(woa_lat)
    nx = len(woa_lon)
    lat_edges = np.linspace(-90.0, 90.0, ny + 1)
    lon_edges = np.linspace(-180.0, 180.0, nx + 1)

    iy = np.searchsorted(lat_edges, lat, side="right") - 1
    ix = np.searchsorted(lon_edges, lon, side="right") - 1
    valid = (
        np.isfinite(f) & np.isfinite(a) & (a > 0)
        & (iy >= 0) & (iy < ny) & (ix >= 0) & (ix < nx)
    )
    idx = iy[valid] * nx + ix[valid]
    num = np.bincount(idx, weights=(f[valid] * a[valid]), minlength=ny * nx)
    den = np.bincount(idx, weights=a[valid], minlength=ny * nx)
    out = np.full(ny * nx, np.nan, dtype=float)
    np.divide(num, den, out=out, where=den > 0)
    return out.reshape(ny, nx)


def basin_on_model_grid(lat2d, lon2d, basin_arr, woa_lat, woa_lon):
    lat = np.asarray(lat2d, dtype=float)
    lon = np.asarray(lon2d, dtype=float)
    lon = ((lon + 180.0) % 360.0) - 180.0
    iy = np.rint(lat - woa_lat[0]).astype(int)
    ix = np.rint(lon - woa_lon[0]).astype(int)
    iy = np.clip(iy, 0, len(woa_lat) - 1)
    ix = np.clip(ix, 0, len(woa_lon) - 1)
    return basin_arr[iy, ix]


def area_weighted_series(field, area, ydim, xdim):
    a = area.broadcast_like(field)
    valid = field.notnull() & a.notnull() & (a > 0)
    num = (field.where(valid, 0) * a.where(valid, 0)).sum(dim=(ydim, xdim))
    den = a.where(valid, 0).sum(dim=(ydim, xdim))
    return (num / den).assign_attrs(units="psu")


def fit_trend_1d(y, t):
    """OLS trend with effective-sample-size correction."""
    y = np.asarray(y, dtype=float)
    t = np.asarray(t, dtype=float)
    good = np.isfinite(y) & np.isfinite(t)
    y = y[good]
    t = t[good]
    n = len(y)
    if n < 3:
        return dict(slope_per_decade=np.nan, se_per_decade=np.nan, ci95=np.nan,
                    r2=np.nan, p=np.nan, neff=np.nan, n=n,
                    intercept=np.nan, slope_per_year=np.nan)
    xm = t.mean()
    ym = y.mean()
    Sxx = ((t - xm) ** 2).sum()
    if Sxx <= 0:
        return dict(slope_per_decade=np.nan, se_per_decade=np.nan, ci95=np.nan,
                    r2=np.nan, p=np.nan, neff=np.nan, n=n,
                    intercept=np.nan, slope_per_year=np.nan)
    slope = ((t - xm) * (y - ym)).sum() / Sxx
    intercept = ym - slope * xm
    resid = y - (intercept + slope * t)
    ss_res = (resid ** 2).sum()
    ss_tot = ((y - ym) ** 2).sum()
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    maxlag = min(n // 4, 10)
    sum_rho = 0.0
    for k in range(1, maxlag + 1):
        a = resid[:-k]
        b = resid[k:]
        ma = a.mean()
        mb = b.mean()
        num = ((a - ma) * (b - mb)).sum()
        den = np.sqrt(((a - ma) ** 2).sum() * ((b - mb) ** 2).sum())
        rho = num / den if den > 0 else 0.0
        sum_rho += (1.0 - k / n) * rho

    neff = n / (1.0 + 2.0 * sum_rho)
    neff = max(2.0, min(float(n), neff))
    sigma2 = ss_res / (n - 2)
    se = np.sqrt(sigma2 / Sxx) * np.sqrt(n / neff)
    df = neff - 2.0
    if df > 0:
        tcrit = st.t.ppf(0.975, df)
        p = 2.0 * st.t.sf(abs(slope / se), df) if se > 0 else 1.0
    else:
        tcrit = np.nan
        p = np.nan

    slope_dec = slope * 10.0
    se_dec = se * 10.0
    ci = tcrit * se_dec
    return dict(slope_per_decade=slope_dec, se_per_decade=se_dec, ci95=ci,
                r2=r2, p=p, neff=neff, n=n, intercept=intercept,
                slope_per_year=slope)


def trend_map_ess(annual, ydim, xdim):
    """Cell-by-cell OLS trend with ESS correction. Returns slope, p, r2, neff."""
    da = annual.transpose("year", ydim, xdim)
    years = da["year"].values.astype(float)
    Y = da.values.astype(float)
    T = Y.shape[0]
    shape = Y.shape[1:]
    N = int(np.prod(shape))
    Y2 = Y.reshape(T, N)

    valid = np.isfinite(Y2)
    n = valid.sum(axis=0)
    slope = np.full(N, np.nan)
    pval = np.full(N, np.nan)
    r2out = np.full(N, np.nan)
    neffout = np.full(N, np.nan)

    full = n == T
    if full.any():
        Yf = Y2[:, full]
        xm = years.mean()
        Sxx = ((years - xm) ** 2).sum()
        ym = Yf.mean(axis=0)
        Sxy = ((years - xm)[:, None] * (Yf - ym)).sum(axis=0)
        sl = Sxy / Sxx
        intercept = ym - sl * xm
        resid = Yf - (intercept + sl * years[:, None])
        ss_res = (resid ** 2).sum(axis=0)
        ss_tot = ((Yf - ym) ** 2).sum(axis=0)
        r2f = np.where(ss_tot > 0, 1.0 - ss_res / ss_tot, np.nan)

        maxlag = min(T // 4, 10)
        sum_rho = np.zeros(Yf.shape[1])
        for k in range(1, maxlag + 1):
            a = resid[:-k, :]
            b = resid[k:, :]
            ma = a.mean(axis=0)
            mb = b.mean(axis=0)
            num = ((a - ma) * (b - mb)).sum(axis=0)
            den = np.sqrt(((a - ma) ** 2).sum(axis=0) * ((b - mb) ** 2).sum(axis=0))
            rho = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
            sum_rho += (1.0 - k / T) * rho

        nefff = T / (1.0 + 2.0 * sum_rho)
        nefff = np.clip(nefff, 2.0, T)
        sigma2 = ss_res / (T - 2)
        se = np.sqrt(sigma2 / Sxx) * np.sqrt(T / nefff)
        df = nefff - 2.0
        with np.errstate(invalid="ignore", divide="ignore"):
            p = 2.0 * st.t.sf(np.abs(sl / se), df)
        p = np.where(np.isfinite(p), p, 1.0)

        slope[full] = sl * 10.0
        pval[full] = p
        r2out[full] = r2f
        neffout[full] = nefff

    partial = (n > 2) & (n < T)
    for j in np.where(partial)[0]:
        y = Y2[:, j]
        m = np.isfinite(y)
        res = fit_trend_1d(y[m], years[m])
        slope[j] = res["slope_per_decade"]
        pval[j] = res["p"]
        r2out[j] = res["r2"]
        neffout[j] = res["neff"]

    coords = {ydim: annual[ydim], xdim: annual[xdim]}
    dims = (ydim, xdim)
    return (
        xr.DataArray(slope.reshape(shape), coords=coords, dims=dims),
        xr.DataArray(pval.reshape(shape), coords=coords, dims=dims),
        xr.DataArray(r2out.reshape(shape), coords=coords, dims=dims),
        xr.DataArray(neffout.reshape(shape), coords=coords, dims=dims),
    )


def compute_sss_annual(root, model, info, y0, y1, area):
    """Monthly SSS averaged over levels in 0-10 m, then annual mean per cell."""
    ds_so, so_files = open_var(root, model, "so", y0, y1)
    so = ds_so["so"]

    ds_vol, _ = open_var(root, model, "volcello", None, fallback_root=hist_root)
    vol = ds_vol["volcello"]

    if "time" in vol.dims:
        yrs = vol.time.dt.year
        m = (yrs >= y0) & (yrs <= y1)
        if bool(m.any()):
            vol = vol.isel(time=np.where(m.values)[0])
        else:
            vol = vol.mean("time")

    if "time" in vol.dims and "time" in so.dims:
        if vol.sizes["time"] == so.sizes["time"]:
            vol = vol.assign_coords(time=so.time)
        else:
            vol = vol.sel(time=so.time, method="nearest")

    lev_name = info["depth_coord"]
    lev_vals = so[lev_name].values.astype(float) * info["depth_scale"]
    sel = np.where((lev_vals >= 0.0) & (lev_vals <= 10.0))[0]
    if len(sel) == 0:
        sel = np.array([0])

    so_sel = so.isel({lev_name: sel})
    vol_sel = vol.isel({lev_name: sel})
    vol_sel = vol_sel.where(area > 0)

    num = (so_sel * vol_sel).sum(dim=lev_name)
    den = vol_sel.sum(dim=lev_name)
    sss = (num / den.where(den > 0)).assign_attrs(units="psu")

    annual = sss.groupby("time.year").mean("time").compute()
    return annual, sss, so_files


# -----------------------------------------------------------------------------
# WOA23 salinity, mask, basin mask, and cell volumes
# -----------------------------------------------------------------------------
woa_file = glob.glob(os.path.join(
    woa_root, "salinity", "1971-2000_climate_normal", "woa23_*_s00_01.nc"
))[0]
woa_ds = xr.open_dataset(woa_file, decode_times=False)

lat_name = "lat" if "lat" in woa_ds.coords else "latitude"
lon_name = "lon" if "lon" in woa_ds.coords else "longitude"
depth_name = "depth" if "depth" in woa_ds.coords else "lev"

woa_ds = woa_ds.sortby(lat_name).sortby(lon_name)

woa_var = None
for v in woa_ds.data_vars:
    if v.lower().endswith("_an"):
        woa_var = v
        break
if woa_var is None:
    for v in woa_ds.data_vars:
        if woa_ds[v].ndim >= 3:
            woa_var = v
            break
if woa_var is None:
    raise RuntimeError("Could not identify WOA23 salinity analysed field")

sal = woa_ds[woa_var]
if "time" in sal.dims:
    sal = sal.isel(time=0)

depth_vals = woa_ds[depth_name].values.astype(float)
sel_depth = np.where((depth_vals >= 0.0) & (depth_vals <= 10.0))[0]
if len(sel_depth) == 0:
    sel_depth = np.arange(min(3, len(depth_vals)))
sal_sel = sal.isel({depth_name: sel_depth})

thick_full = np.abs(np.gradient(depth_vals))
thick_sel = thick_full[sel_depth]
thick_da = xr.DataArray(
    thick_sel, dims=(depth_name,), coords={depth_name: depth_vals[sel_depth]}
)

woa_lat = woa_ds[lat_name].values.astype(float)
woa_lon = woa_ds[lon_name].values.astype(float)
area_woa = regular_grid_area(woa_lat, woa_lon)
area_woa_da = xr.DataArray(
    area_woa, dims=(lat_name, lon_name), coords={lat_name: woa_lat, lon_name: woa_lon}
)
vol_woa_da = area_woa_da * thick_da
sss_woa = (
    (sal_sel * vol_woa_da).sum(dim=depth_name) / vol_woa_da.sum(dim=depth_name)
).assign_attrs(units="psu")

landsea = read_woa_landsea(os.path.join(woa_masks_root, "landsea_01.msk"))
woa_ocean = landsea > 1.0

basin_raw = read_woa_basin(os.path.join(woa_masks_root, "basinmask_01.msk"))
basin_folded = fold_basin_array(basin_raw, woa_lat, woa_lon)

# -----------------------------------------------------------------------------
# Process each CMIP6 model
# -----------------------------------------------------------------------------
model_results = {}

for model in models:
    info = model_info[model]
    ydim = info["y_dim"]
    xdim = info["x_dim"]

    # Area weights on the native model grid.
    ds_area, _ = open_var(hist_root, model, "areacello", None, fallback_root=hist_root, first_only=True)
    area = ds_area["areacello"].squeeze()

    if info.get("ocean_mask_from") == "sftof":
        ds_sft, _ = open_var(hist_root, model, "sftof", None, fallback_root=hist_root, first_only=True)
        sft = ds_sft["sftof"].squeeze()
        area = area.where(sft > 0)

    annual_h, sss_h_monthly, hist_so_files = compute_sss_annual(
        hist_root, model, info, 1950, 2014, area
    )
    annual_s, sss_s_monthly, ssp_so_files = compute_sss_annual(
        ssp_root, model, info, 2071, 2100, area
    )

    lat_da = sss_h_monthly[info["lat_coord"]]
    lon_da = sss_h_monthly[info["lon_coord"]]
    if lat_da.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon_da.values.astype(float), lat_da.values.astype(float))
    else:
        lat2d = lat_da.values.astype(float)
        lon2d = lon_da.values.astype(float)

    gm_h = area_weighted_series(annual_h, area, ydim, xdim)
    gm_s = area_weighted_series(annual_s, area, ydim, xdim)
    fit_h = fit_trend_1d(gm_h.values, gm_h["year"].values)
    fit_s = fit_trend_1d(gm_s.values, gm_s["year"].values)

    trend_h, p_h, r2_h, neff_h = trend_map_ess(annual_h, ydim, xdim)

    # 1971-2000 mean state, regridded to WOA grid for the comparison.
    clim = annual_h.sel(year=slice(1971, 2000)).mean("year")
    clim_woa = regrid_conservative(
        clim.values, lat2d, lon2d, area.values, woa_lat, woa_lon
    )

    # End-of-century change: ssp245 2071-2100 minus historical 1995-2014.
    hist_mean = annual_h.sel(year=slice(1995, 2014)).mean("year")
    ssp_mean = annual_s.sel(year=slice(2071, 2100)).mean("year")
    delta = ssp_mean - hist_mean

    basin_labels = basin_on_model_grid(lat2d, lon2d, basin_folded, woa_lat, woa_lon)

    valid_delta = (
        np.isfinite(delta.values)
        & np.isfinite(area.values)
        & (area.values > 0)
    )
    if valid_delta.any():
        global_delta = float(
            (delta.values[valid_delta] * area.values[valid_delta]).sum()
            / area.values[valid_delta].sum()
        )
    else:
        global_delta = np.nan

    basin_delta = {}
    basin_counts = {}
    for code in [1, 2, 3]:
        valid_b = valid_delta & (basin_labels == code)
        basin_counts[code] = int(valid_b.sum())
        if valid_b.any():
            basin_delta[code] = float(
                (delta.values[valid_b] * area.values[valid_b]).sum()
                / area.values[valid_b].sum()
            )
        else:
            basin_delta[code] = np.nan

    trend_woa = regrid_conservative(
        trend_h.values, lat2d, lon2d, area.values, woa_lat, woa_lon
    )
    p_woa = regrid_conservative(
        p_h.values, lat2d, lon2d, area.values, woa_lat, woa_lon
    )

    model_results[model] = dict(
        area=area,
        lat2d=lat2d,
        lon2d=lon2d,
        annual_h=annual_h,
        annual_s=annual_s,
        gm_h=gm_h,
        gm_s=gm_s,
        fit_h=fit_h,
        fit_s=fit_s,
        trend_h=trend_h,
        p_h=p_h,
        r2_h=r2_h,
        neff_h=neff_h,
        clim_woa=clim_woa,
        trend_woa=trend_woa,
        p_woa=p_woa,
        delta=delta,
        global_delta=global_delta,
        basin_delta=basin_delta,
        basin_counts=basin_counts,
        hist_variant=get_variant(hist_so_files),
        ssp_variant=get_variant(ssp_so_files),
        ydim=ydim,
        xdim=xdim,
    )

# -----------------------------------------------------------------------------
# Common comparison mask on the WOA23 1-degree grid
# -----------------------------------------------------------------------------
mask_common = np.isfinite(sss_woa.values) & woa_ocean
for model in models:
    mask_common &= np.isfinite(model_results[model]["clim_woa"])

woa_vals = sss_woa.values
metrics = {}
for model in models:
    m = model_results[model]["clim_woa"]
    valid = mask_common
    w = area_woa[valid]
    mean_m = float(np.sum(m[valid] * w) / np.sum(w))
    mean_o = float(np.sum(woa_vals[valid] * w) / np.sum(w))
    bias = m - woa_vals
    rmse = float(np.sqrt(np.sum(w * (bias[valid] ** 2)) / np.sum(w)))
    corr = float(np.corrcoef(m[valid], woa_vals[valid])[0, 1])
    metrics[model] = dict(
        mean_model=mean_m,
        mean_woa=mean_o,
        bias=mean_m - mean_o,
        rmse=rmse,
        corr=corr,
        n_cells=int(valid.sum()),
    )
    model_results[model]["bias_woa"] = np.where(valid, bias, np.nan)

# -----------------------------------------------------------------------------
# Store plotting arrays
# -----------------------------------------------------------------------------
bias_maps = {model: model_results[model]["bias_woa"] for model in models}
trend_maps_woa = {model: model_results[model]["trend_woa"] for model in models}
p_maps_woa = {model: model_results[model]["p_woa"] for model in models}

time_series = {}
for model in models:
    r = model_results[model]
    time_series[model] = dict(
        hist_years=r["gm_h"]["year"].values.astype(float),
        hist_vals=r["gm_h"].values.astype(float),
        ssp_years=r["gm_s"]["year"].values.astype(float),
        ssp_vals=r["gm_s"].values.astype(float),
        fit_h=r["fit_h"],
        fit_s=r["fit_s"],
    )

# -----------------------------------------------------------------------------
# Print headline results
# -----------------------------------------------------------------------------
print("WOA23 salinity file:", woa_file)
print("Common comparison mask cells:", int(mask_common.sum()), "of", int(mask_common.size))

print("\nMean-state comparison, 1971-2000, WOA23 1-degree grid, area-weighted:")
print(f"{'Model':15s} {'mean_model':>10s} {'mean_WOA':>10s} {'bias':>8s} {'RMSE':>8s} {'corr':>8s}")
for model in models:
    m = metrics[model]
    print(
        f"{model:15s} {m['mean_model']:10.4f} {m['mean_woa']:10.4f} "
        f"{m['bias']:8.4f} {m['rmse']:8.4f} {m['corr']:8.4f}"
    )

print("\nGlobal-mean SSS trends, annual, area-weighted, ESS-corrected:")
for model in models:
    fh = model_results[model]["fit_h"]
    fs = model_results[model]["fit_s"]
    sig_h = abs(fh["slope_per_decade"]) > fh["ci95"] if np.isfinite(fh["ci95"]) else False
    sig_s = abs(fs["slope_per_decade"]) > fs["ci95"] if np.isfinite(fs["ci95"]) else False
    ratio = (
        fs["slope_per_decade"] / fh["slope_per_decade"]
        if np.isfinite(fh["slope_per_decade"]) and fh["slope_per_decade"] != 0
        else np.nan
    )
    print(
        f"{model:15s} hist={fh['slope_per_decade']:.4g} psu/dec, "
        f"95%CI=[{fh['slope_per_decade']-fh['ci95']:.4g}, {fh['slope_per_decade']+fh['ci95']:.4g}], "
        f"R2={fh['r2']:.4f}, sig={sig_h}; "
        f"ssp={fs['slope_per_decade']:.4g} psu/dec, "
        f"95%CI=[{fs['slope_per_decade']-fs['ci95']:.4g}, {fs['slope_per_decade']+fs['ci95']:.4g}], "
        f"R2={fs['r2']:.4f}, sig={sig_s}; ratio(ssp/hist)={ratio:.1f}"
    )
    print(
        f"    variants: historical={model_results[model]['hist_variant']}, "
        f"ssp245={model_results[model]['ssp_variant']}"
    )

print("\nDelta SSS, ssp245 2071-2100 minus historical 1995-2014:")
basin_names = {1: "Atlantic", 2: "Pacific", 3: "Indian"}
for model in models:
    r = model_results[model]
    print(
        f"{model:15s} global={r['global_delta']:.4f} psu; "
        + "; ".join(
            f"{basin_names[c]}={r['basin_delta'][c]:.4f} psu (n={r['basin_counts'][c]})"
            for c in [1, 2, 3]
        )
    )
    if r["hist_variant"] != r["ssp_variant"]:
        print(
            f"    Caveat: {model} uses historical {r['hist_variant']} but ssp245 {r['ssp_variant']}; "
            "the delta includes an inter-member internal-variability difference and is not a "
            "single continuous realisation."
        )

print(
    "\nThe ssp245/historical slope ratio is descriptive only: it compares a 30-yr period "
    "with a 65-yr period, different forcing levels, and for CESM2 different realisations. "
    "It is not a like-for-like forced-amplification ratio."
)
