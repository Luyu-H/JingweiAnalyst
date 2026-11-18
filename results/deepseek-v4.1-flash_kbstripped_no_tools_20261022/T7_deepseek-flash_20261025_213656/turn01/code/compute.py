import os
import re
import glob
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from scipy.ndimage import distance_transform_edt

# -----------------------------------------------------------------------------
# Paths and constants
# -----------------------------------------------------------------------------
CMIP6_ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"
WOA_ROOT = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
WOA_MASK_ROOT = "/data1/lorraine-20419/ocean_data/WOA23/masks"

DENSITY = 1025.0  # kg m-3, prescribed constant density
# 1 mol m-3 in water of density 1025 kg m-3 = 1e6 umol / 1025 kg
CMIP6_TO_UMOL_KG = 1.0e6 / DENSITY

MODELS = {
    "CESM2": {
        "y_dim": "nlat",
        "x_dim": "nlon",
        "lat_coord": "lat",
        "lon_coord": "lon",
        "depth_coord": "lev",
        "depth_scale": 0.01,  # lev is in cm
    },
    "GFDL-ESM4": {
        "y_dim": "lat",
        "x_dim": "lon",
        "lat_coord": "lat",
        "lon_coord": "lon",
        "depth_coord": "lev",
        "depth_scale": 1.0,
    },
    "MPI-ESM1-2-LR": {
        "y_dim": "j",
        "x_dim": "i",
        "lat_coord": "latitude",
        "lon_coord": "longitude",
        "depth_coord": "lev",
        "depth_scale": 1.0,
    },
}

SOURCE_ORDER = ["WOA23", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]


# -----------------------------------------------------------------------------
# CMIP6 file discovery and opening
# -----------------------------------------------------------------------------
def cmip6_files(model, var, table):
    """Return CMIP6 files for a model/variable/table, filtering to 1971-2000."""
    d = os.path.join(CMIP6_ROOT, model, var)
    pattern = os.path.join(d, f"{var}_{table}_{model}_historical_*.nc")
    files = sorted(glob.glob(pattern))
    if not files:
        files = sorted(glob.glob(os.path.join(d, "*.nc")))

    selected = []
    for f in files:
        m = re.search(r"_(\d{6})-(\d{6})\.nc$", os.path.basename(f))
        if m:
            start = int(m.group(1))
            end = int(m.group(2))
            if end < 197101 or start > 200012:
                continue
        selected.append(f)
    if not selected:
        selected = files
    return selected


def open_cmip6_var(model, var, table="Omon"):
    """Open a CMIP6 monthly field and subset to 1971-2000."""
    files = cmip6_files(model, var, table)
    if not files:
        raise FileNotFoundError(f"No CMIP6 files for {model} {var} {table}")

    if len(files) > 1:
        ds = xr.open_mfdataset(
            files,
            combine="by_coords",
            use_cftime=True,
            decode_times=True,
            chunks={"time": 120},
        )
    else:
        ds = xr.open_dataset(
            files[0],
            use_cftime=True,
            decode_times=True,
            chunks={"time": 120},
        )

    if "time" in ds.dims:
        ds = ds.sortby("time")
        if table == "Omon":
            mask = (ds.time.dt.year >= 1971) & (ds.time.dt.year <= 2000)
            idx = np.where(mask.values)[0]
            ds = ds.isel(time=idx)
    return ds


def open_cmip6_fixed(model, var, info, lev_idx=None):
    """Open a CMIP6 fixed field, or GFDL volcello from Omon and time-average it."""
    table = "Ofx"
    if var == "volcello" and model == "GFDL-ESM4":
        table = "Omon"

    files = cmip6_files(model, var, table)
    if not files and table == "Ofx":
        table = "Omon"
        files = cmip6_files(model, var, table)
    if not files:
        raise FileNotFoundError(f"No CMIP6 fixed files for {model} {var}")

    if table == "Omon":
        if len(files) > 1:
            ds = xr.open_mfdataset(
                files,
                combine="by_coords",
                use_cftime=True,
                decode_times=True,
                chunks={"time": 120},
            )
        else:
            ds = xr.open_dataset(
                files[0],
                use_cftime=True,
                decode_times=True,
                chunks={"time": 120},
            )
        ds = ds.sortby("time")
        mask = (ds.time.dt.year >= 1971) & (ds.time.dt.year <= 2000)
        idx = np.where(mask.values)[0]
        ds = ds.isel(time=idx)
        da = ds[var].mean("time").compute()
    else:
        if len(files) > 1:
            ds = xr.open_mfdataset(files, combine="by_coords")
        else:
            ds = xr.open_dataset(files[0])
        da = ds[var]
        if "time" in da.dims:
            da = da.isel(time=0)
        da = da.compute()

    if lev_idx is not None and info["depth_coord"] in da.dims:
        da = da.isel({info["depth_coord"]: lev_idx})
    return da


def get_latlon_2d(ds, info):
    """Return 2-D latitude and longitude arrays for a model grid."""
    lat = ds[info["lat_coord"]].values
    lon = ds[info["lon_coord"]].values
    if lat.ndim == 1 and lon.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon, lat)
    else:
        lat2d, lon2d = np.broadcast_arrays(lat, lon)
    return lat2d, lon2d


# -----------------------------------------------------------------------------
# Time-weighted monthly climatology
# -----------------------------------------------------------------------------
def monthly_climatology(da, time_dim="time"):
    """
    Calendar-month-length weighted monthly climatology and annual time mean.
    Returns (monthly_clim[month, ...], annual_mean[...]).
    """
    w = da[time_dim].dt.days_in_month.astype(float)
    month = da[time_dim].dt.month

    num = (da * w).groupby(month).sum(dim=time_dim)
    den = w.groupby(month).sum(dim=time_dim)
    clim = (num / den).compute()
    ann = ((da * w).sum(dim=time_dim) / w.sum(dim=time_dim)).compute()
    return clim, ann


# -----------------------------------------------------------------------------
# WOA23 opening and geometry
# -----------------------------------------------------------------------------
WOA_VAR_INFO = {
    "nitrate": {"subdir": "nitrate/all_time_periods", "tag": "all", "prefix": "n"},
    "phosphate": {"subdir": "phosphate/all_time_periods", "tag": "all", "prefix": "p"},
}


def open_woa_var(variable):
    info = WOA_VAR_INFO[variable]
    fname = f"woa23_{info['tag']}_{info['prefix']}00_01.nc"
    path = os.path.join(WOA_ROOT, info["subdir"], fname)
    return xr.open_dataset(path, decode_times=False)


def get_woa_field(ds, prefix):
    """Find the objectively analysed field in a WOA file."""
    candidates = []
    for v in ds.data_vars:
        if v in ds.coords:
            continue
        dims = ds[v].dims
        if "depth" in dims and "lat" in dims and "lon" in dims:
            candidates.append(v)
    if not candidates:
        raise ValueError("No WOA 3-D field found")

    preferred = f"{prefix}_an"
    if preferred in candidates:
        return ds[preferred]
    for v in candidates:
        if v.lower().endswith("_an"):
            return ds[v]
    return ds[candidates[0]]


def woa_depth_bounds(ds, depth_name="depth"):
    """Return lower and upper depth bounds for WOA levels."""
    d = ds[depth_name].values
    n = len(d)

    if "depth_bnds" in ds.variables:
        b = ds["depth_bnds"].values
        if b.ndim == 2:
            if b.shape[0] == n and b.shape[1] == 2:
                lower, upper = b[:, 0], b[:, 1]
            elif b.shape[1] == n and b.shape[0] == 2:
                lower, upper = b[0, :], b[1, :]
            else:
                raise ValueError("Unexpected depth_bnds shape")
        elif b.ndim == 1 and len(b) == 2 * n:
            lower, upper = b[:n], b[n:]
        else:
            raise ValueError("Unexpected depth_bnds shape")
    else:
        lower = np.empty(n)
        upper = np.empty(n)
        lower[0] = 0.0
        for i in range(1, n):
            lower[i] = 0.5 * (d[i - 1] + d[i])
        for i in range(n - 1):
            upper[i] = 0.5 * (d[i] + d[i + 1])
        upper[-1] = d[-1] + 0.5 * (d[-1] - d[-2])
        if n == 1:
            upper[0] = d[0]

    return np.minimum(lower, upper), np.maximum(lower, upper)


def woa_area(lat, lon):
    """Approximate 1-degree cell area on a sphere."""
    R = 6371000.0
    dlon = np.deg2rad(1.0)
    area_lat = R**2 * dlon * (
        np.sin(np.deg2rad(lat + 0.5)) - np.sin(np.deg2rad(lat - 0.5))
    )
    return np.repeat(area_lat[:, None], len(lon), axis=1)


# -----------------------------------------------------------------------------
# WOA basin mask
# -----------------------------------------------------------------------------
def read_basinmask_file(path, lat_woa, lon_woa):
    """
    Read the WOA sparse basin mask.
    The file has one row per ocean grid point, with latitude, longitude and the
    basin labels for the standard levels.  Because depth_dependent=False we only
    need the surface label (third column).
    """
    lats, lons, vals = [], [], []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#") or s.startswith("!"):
                continue
            parts = s.split()
            if len(parts) < 3:
                continue
            try:
                lat = float(parts[0])
                lon = float(parts[1])
                val = float(parts[2])
            except ValueError:
                continue
            lats.append(lat)
            lons.append(lon)
            vals.append(val)

    lats = np.asarray(lats)
    lons = np.asarray(lons)
    vals = np.asarray(vals)

    basin = np.full((len(lat_woa), len(lon_woa)), np.nan, dtype=float)
    lat_idx = np.abs(lat_woa[:, None] - lats[None, :]).argmin(axis=0)
    lon_idx = np.abs(lon_woa[:, None] - lons[None, :]).argmin(axis=0)
    basin[lat_idx, lon_idx] = vals
    return basin


def map_basin_codes(basin):
    """
    Map WOA basin numbers to parent basins:
      1 Atlantic, 2 Pacific, 3 Indian, 4 Arctic, 5 Southern.
    Marginal-sea codes >5 are assigned to the nearest open-ocean parent
    (include_marginal_seas=True).
    """
    parent = np.full_like(basin, np.nan, dtype=float)

    parent[np.isclose(basin, 1.0)] = 1.0
    parent[np.isclose(basin, 2.0)] = 2.0
    parent[np.isclose(basin, 3.0)] = 3.0
    parent[np.isclose(basin, 4.0)] = 4.0
    parent[np.isclose(basin, 5.0)] = 5.0

    finite = np.isfinite(basin)
    rounded = np.zeros_like(basin, dtype=int)
    rounded[finite] = np.round(basin[finite]).astype(int)
    unknown = finite & (basin > 0.0) & ~np.isin(rounded, [1, 2, 3, 4, 5])

    if np.any(unknown):
        distances = []
        for code in (1.0, 2.0, 3.0):
            distances.append(distance_transform_edt(~(parent == code)))
        distances = np.stack(distances, axis=0)
        nearest = np.argmin(distances, axis=0).astype(float) + 1.0
        parent[unknown] = nearest[unknown]

    return parent


# -----------------------------------------------------------------------------
# Weighted metrics and mapping helpers
# -----------------------------------------------------------------------------
def metrics_from_arrays(np_arr, no3_arr, po4_arr, vol_arr, basin_arr=None):
    """Compute volume-weighted metrics and unweighted OLS regression."""
    total_vol = float(np.sum(vol_arr))

    mean_np = float(np.sum(np_arr * vol_arr) / total_vol)
    ratio_means = float(np.sum(no3_arr * vol_arr) / np.sum(po4_arr * vol_arr))

    mask_14_18 = (np_arr > 14.0) & (np_arr < 18.0)
    vol_14_18 = float(np.sum(vol_arr[mask_14_18]))
    pct_14_18 = 100.0 * vol_14_18 / total_vol

    mask_below12 = np_arr < 12.0
    vol_below12 = float(np.sum(vol_arr[mask_below12]))
    pct_below12 = 100.0 * vol_below12 / total_vol

    slope, intercept = np.polyfit(po4_arr, no3_arr, 1)
    pred = slope * po4_arr + intercept
    ss_res = np.sum((no3_arr - pred) ** 2)
    ss_tot = np.sum((no3_arr - np.mean(no3_arr)) ** 2)
    r2 = float(1.0 - ss_res / ss_tot)

    out = {
        "mean_NP": mean_np,
        "ratio_of_means": ratio_means,
        "slope_a": float(slope),
        "intercept_b": float(intercept),
        "r2": r2,
        "pct_14_18": float(pct_14_18),
        "vol_14_18_m3": vol_14_18,
        "pct_below12": float(pct_below12),
        "vol_below12_m3": vol_below12,
        "total_volume_m3": total_vol,
    }

    if basin_arr is not None:
        for code, name in [(1.0, "Pacific"), (2.0, "Atlantic"), (3.0, "Indian")]:
            m = basin_arr == code
            if np.any(m):
                out[f"mean_NP_{name}"] = float(
                    np.sum(np_arr[m] * vol_arr[m]) / np.sum(vol_arr[m])
                )
                if name in ("Pacific", "Indian"):
                    out[f"pct_below12_{name}"] = float(
                        100.0
                        * np.sum(vol_arr[m & (np_arr < 12.0)])
                        / np.sum(vol_arr[m])
                    )
                else:
                    out[f"pct_below12_{name}"] = np.nan
            else:
                out[f"mean_NP_{name}"] = np.nan
                out[f"pct_below12_{name}"] = np.nan

    return out


def vertical_weighted_mean(field, vol, valid):
    """Volume-weighted vertical mean of a (lev, y, x) field."""
    num = np.nansum(np.where(valid, field * vol, np.nan), axis=0)
    den = np.nansum(np.where(valid, vol, np.nan), axis=0)
    out = np.full(num.shape, np.nan, dtype=float)
    mask = den > 0.0
    out[mask] = num[mask] / den[mask]
    return out


def conservative_regrid(values2d, lat2d, lon2d, area2d, common_lat, common_lon):
    """
    Area-weighted conservative binning to a regular 1-degree grid.
    This is used for the display map only.
    """
    nlat = len(common_lat)
    nlon = len(common_lon)

    lon360 = lon2d % 360.0
    j = np.floor(lon360).astype(int) % nlon
    i = np.floor(lat2d + 90.0).astype(int)
    i = np.clip(i, 0, nlat - 1)

    valid = (
        np.isfinite(values2d)
        & np.isfinite(area2d)
        & (area2d > 0.0)
    )

    sum_area = np.zeros((nlat, nlon), dtype=float)
    sum_val = np.zeros((nlat, nlon), dtype=float)

    np.add.at(sum_area, (i[valid], j[valid]), area2d[valid])
    np.add.at(sum_val, (i[valid], j[valid]), values2d[valid] * area2d[valid])

    out = np.full((nlat, nlon), np.nan, dtype=float)
    mask = sum_area > 0.0
    out[mask] = sum_val[mask] / sum_area[mask]
    return out


# -----------------------------------------------------------------------------
# WOA basin mask set-up
# -----------------------------------------------------------------------------
ds_basin_ref = open_woa_var("nitrate")
lat_woa = ds_basin_ref["lat"].values
lon_woa = ds_basin_ref["lon"].values

basin_raw = read_basinmask_file(
    os.path.join(WOA_MASK_ROOT, "basinmask_01.msk"),
    lat_woa,
    lon_woa,
)
basin_parent_2d = map_basin_codes(basin_raw)

# For model cells, fill land in the basin mask with the nearest ocean basin.
valid_basin = np.isfinite(basin_parent_2d)
basin_indices = distance_transform_edt(
    ~valid_basin, return_distances=False, return_indices=True
)
basin_parent_filled = basin_parent_2d[tuple(basin_indices)]


# -----------------------------------------------------------------------------
# Process WOA23
# -----------------------------------------------------------------------------
def process_woa():
    ds_n = open_woa_var("nitrate")
    ds_p = open_woa_var("phosphate")

    da_n = get_woa_field(ds_n, "n")
    da_p = get_woa_field(ds_p, "p")

    if "time" in da_n.dims:
        da_n = da_n.isel(time=0)
    if "time" in da_p.dims:
        da_p = da_p.isel(time=0)

    depth = ds_n["depth"].values
    lev_idx = np.where((depth >= 200.0) & (depth <= 1000.0))[0]
    if lev_idx.size == 0:
        raise ValueError("No WOA levels in 200-1000 m")

    da_n = da_n.isel(depth=lev_idx)
    da_p = da_p.isel(depth=lev_idx)

    lat = ds_n["lat"].values
    lon = ds_n["lon"].values

    lower, upper = woa_depth_bounds(ds_n, "depth")
    thick = np.abs(upper[lev_idx] - lower[lev_idx])
    area2d = woa_area(lat, lon)
    vol_woa = thick[:, None, None] * area2d[None, :, :]

    no3_all = da_n.values
    po4_all = da_p.values
    np_all = np.where(po4_all > 0.0, no3_all / po4_all, np.nan)

    valid = (
        np.isfinite(np_all)
        & np.isfinite(no3_all)
        & np.isfinite(po4_all)
        & (vol_woa > 0.0)
    )

    np_vals = np_all[valid]
    no3_vals = no3_all[valid]
    po4_vals = po4_all[valid]
    vol_vals = vol_woa[valid]

    basin_b = np.broadcast_to(basin_parent_2d[None, :, :], np_all.shape)
    basin_vals = basin_b[valid]

    metrics = metrics_from_arrays(
        np_vals, no3_vals, po4_vals, vol_vals, basin_vals
    )
    metrics["po4_vals"] = po4_vals
    metrics["no3_vals"] = no3_vals
    metrics["volume_source"] = "thickness from WOA depth bounds × 1° cell area"
    metrics["density_assumption"] = "not applicable (WOA mass basis)"
    metrics["map_np"] = vertical_weighted_mean(np_all, vol_woa, valid)

    lon2d, lat2d = np.meshgrid(lon, lat)
    metrics["map_lat2d"] = lat2d
    metrics["map_lon2d"] = lon2d
    metrics["map_area2d"] = area2d
    metrics["basin2d"] = basin_parent_2d
    return metrics


# -----------------------------------------------------------------------------
# Process one CMIP6 model
# -----------------------------------------------------------------------------
def process_model(model):
    info = MODELS[model]
    depth_name = info["depth_coord"]
    y_dim = info["y_dim"]
    x_dim = info["x_dim"]

    ds_no3 = open_cmip6_var(model, "no3", "Omon")
    da_no3 = ds_no3["no3"]

    lev = ds_no3[depth_name].values
    depth_m = lev * info["depth_scale"]
    lev_idx = np.where((depth_m >= 200.0) & (depth_m <= 1000.0))[0]
    if lev_idx.size == 0:
        raise ValueError(f"No {model} levels in 200-1000 m")

    da_no3 = da_no3.isel({depth_name: lev_idx})

    ds_po4 = open_cmip6_var(model, "po4", "Omon")
    da_po4 = ds_po4["po4"].isel({depth_name: lev_idx})

    no3_monthly, no3_ann = monthly_climatology(da_no3)
    po4_monthly, po4_ann = monthly_climatology(da_po4)

    no3_ann_umol = (no3_ann * CMIP6_TO_UMOL_KG).assign_attrs(
        units="umol kg-1", density_assumption="1025 kg m-3"
    )
    po4_ann_umol = (po4_ann * CMIP6_TO_UMOL_KG).assign_attrs(
        units="umol kg-1", density_assumption="1025 kg m-3"
    )

    np_monthly = xr.where(
        po4_monthly > 0.0,
        no3_monthly / po4_monthly,
        np.nan,
    )
    np_field = np_monthly.mean(dim="month").assign_attrs(units="1")

    vol = open_cmip6_fixed(model, "volcello", info, lev_idx=lev_idx)
    vol = vol.transpose(*np_field.dims)

    area = open_cmip6_fixed(model, "areacello", info, lev_idx=None)
    area = area.transpose(y_dim, x_dim)
    area2d = area.values

    lat2d, lon2d = get_latlon_2d(ds_no3, info)

    np_vals_all = np_field.values
    no3_vals_all = no3_ann_umol.values
    po4_vals_all = po4_ann_umol.values
    vol_vals_all = vol.values

    valid = (
        np.isfinite(np_vals_all)
        & np.isfinite(no3_vals_all)
        & np.isfinite(po4_vals_all)
        & (vol_vals_all > 0.0)
    )

    np_vals = np_vals_all[valid]
    no3_vals = no3_vals_all[valid]
    po4_vals = po4_vals_all[valid]
    vol_vals = vol_vals_all[valid]

    # Map WOA surface basin labels onto the model grid by nearest 1-degree cell.
    lon_mod = ((lon2d + 180.0) % 360.0) - 180.0
    lat_idx = np.floor(lat2d + 90.0).astype(int)
    lon_idx = np.floor(lon_mod + 180.0).astype(int) % len(lon_woa)
    lat_idx = np.clip(lat_idx, 0, len(lat_woa) - 1)
    basin2d = basin_parent_filled[lat_idx, lon_idx]

    basin_b = np.broadcast_to(basin2d[None, :, :], np_vals_all.shape)
    basin_vals = basin_b[valid]

    metrics = metrics_from_arrays(
        np_vals, no3_vals, po4_vals, vol_vals, basin_vals
    )
    metrics["po4_vals"] = po4_vals
    metrics["no3_vals"] = no3_vals
    metrics["volume_source"] = (
        "volcello (Omon 1971-2000 time mean)"
        if model == "GFDL-ESM4"
        else "volcello (Ofx)"
    )
    metrics["density_assumption"] = "1025 kg m-3"
    metrics["map_np"] = vertical_weighted_mean(
        np_vals_all, vol_vals_all, valid
    )
    metrics["map_lat2d"] = lat2d
    metrics["map_lon2d"] = lon2d
    metrics["map_area2d"] = area2d
    metrics["basin2d"] = basin2d
    return metrics


# -----------------------------------------------------------------------------
# Compute all four sources
# -----------------------------------------------------------------------------
source_results = {}
source_results["WOA23"] = process_woa()
for model in ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]:
    source_results[model] = process_model(model)

# Regrid the 200-1000 m horizontal N:P field for Figure 2 only.
common_lat = np.arange(-89.5, 90.0, 1.0)
common_lon = np.arange(0.5, 360.0, 1.0)

for name in SOURCE_ORDER:
    r = source_results[name]
    r["map_common"] = conservative_regrid(
        r["map_np"],
        r["map_lat2d"],
        r["map_lon2d"],
        r["map_area2d"],
        common_lat,
        common_lon,
    )
    r["map_common_anom"] = r["map_common"] - 16.0

# Shared limits for Figure 1.
all_po4 = np.concatenate([source_results[s]["po4_vals"] for s in SOURCE_ORDER])
all_no3 = np.concatenate([source_results[s]["no3_vals"] for s in SOURCE_ORDER])
xmax = float(np.nanpercentile(all_po4, 99.8) * 1.08)
ymax = float(np.nanpercentile(all_no3, 99.8) * 1.08)
xmax = max(xmax, 4.0)
ymax = max(ymax, 50.0)
fig1_xlim = (0.0, xmax)
fig1_ylim = (0.0, ymax)

# Robust shared colour limit for Figure 2.
pooled_abs_anom = np.concatenate(
    [np.abs(source_results[s]["map_common_anom"]).ravel() for s in SOURCE_ORDER]
)
vmax_map = float(np.nanpercentile(pooled_abs_anom, 98.0))
if not np.isfinite(vmax_map) or vmax_map <= 0.0:
    vmax_map = 1.0

# -----------------------------------------------------------------------------
# Summary table
# -----------------------------------------------------------------------------
rows = []
for name in SOURCE_ORDER:
    r = source_results[name]
    rows.append(
        {
            "source": name,
            "mean_NP": r["mean_NP"],
            "ratio_of_means": r["ratio_of_means"],
            "slope_a": r["slope_a"],
            "intercept_b": r["intercept_b"],
            "r2": r["r2"],
            "pct_14_18": r["pct_14_18"],
            "vol_14_18_m3": r["vol_14_18_m3"],
            "pct_below12": r["pct_below12"],
            "pct_below12_Pacific": r.get("pct_below12_Pacific", np.nan),
            "pct_below12_Indian": r.get("pct_below12_Indian", np.nan),
            "mean_NP_Pacific": r.get("mean_NP_Pacific", np.nan),
            "mean_NP_Atlantic": r.get("mean_NP_Atlantic", np.nan),
            "mean_NP_Indian": r.get("mean_NP_Indian", np.nan),
            "total_volume_m3": r["total_volume_m3"],
            "volume_source": r["volume_source"],
        }
    )

results_table = pd.DataFrame(rows)
results_table_rounded = results_table.copy()

for col in [
    "mean_NP",
    "ratio_of_means",
    "slope_a",
    "intercept_b",
    "mean_NP_Pacific",
    "mean_NP_Atlantic",
    "mean_NP_Indian",
]:
    results_table_rounded[col] = results_table_rounded[col].round(2)

results_table_rounded["r2"] = results_table_rounded["r2"].round(3)

for col in [
    "pct_14_18",
    "pct_below12",
    "pct_below12_Pacific",
    "pct_below12_Indian",
]:
    results_table_rounded[col] = results_table_rounded[col].round(1)

print("Rounded results:")
print(results_table_rounded.to_string(index=False))
