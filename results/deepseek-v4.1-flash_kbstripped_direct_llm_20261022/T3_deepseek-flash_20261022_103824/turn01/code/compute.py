# ============================================================
# Computation block
# ============================================================

import os
import glob
import re
import numpy as np
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree

# -----------------------------
# Helper functions
# -----------------------------
def to_numpy_float(da):
    """Return a float numpy array from an xarray object / array."""
    v = da.values if hasattr(da, "values") else da
    if np.ma.isMaskedArray(v):
        v = v.filled(np.nan)
    return np.asarray(v, dtype=float)


def align_2d(arr, target_shape):
    """Align a 2-D array to target_shape, transposing or reshaping if needed."""
    arr = np.asarray(arr, dtype=float)
    arr = np.squeeze(arr)
    if arr.shape == target_shape:
        return arr
    if arr.ndim == 2 and arr.T.shape == target_shape:
        return arr.T
    if arr.size == int(np.prod(target_shape)):
        return arr.reshape(target_shape)
    raise ValueError(f"Cannot align array shape {arr.shape} to {target_shape}")


def lonlat_to_xyz(lat, lon):
    """Convert lon/lat in degrees to 3-D Cartesian unit vectors."""
    lat_r = np.deg2rad(np.asarray(lat, dtype=float))
    lon_r = np.deg2rad(np.mod(np.asarray(lon, dtype=float), 360.0))
    return np.column_stack((
        np.cos(lat_r) * np.cos(lon_r),
        np.cos(lat_r) * np.sin(lon_r),
        np.sin(lat_r)
    ))


def read_landsea_mask(path, lat_target, lon_target):
    """
    Read the WOA23 land/sea mask sparse/full table and rebuild a
    (lat, lon) array.  Ocean is mask > 1.

    The file has arbitrary header/comment rows, so instead of relying on
    pandas header sniffing we parse line by line and keep only the rows
    whose first two fields are numeric (the actual data records).
    """
    lat_list = []
    lon_list = []
    val_list = []

    with open(path, "r") as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 3:
                continue
            try:
                la = float(parts[0])
                lo = float(parts[1])
                va = float(parts[2])
            except ValueError:
                # header / title / comment-like line
                continue
            lat_list.append(la)
            lon_list.append(lo)
            val_list.append(va)

    lat_vals = np.asarray(lat_list, dtype=float)
    lon_vals = np.asarray(lon_list, dtype=float)
    val_vals = np.asarray(val_list, dtype=float)

    # If the file happened to use a 0..360 convention, fold it back.
    if lon_target[0] < 0 and lon_vals.size and lon_vals.min() >= 0 and lon_vals.max() > 180:
        lon_vals = np.where(lon_vals > 180, lon_vals - 360.0, lon_vals)

    nlat = len(lat_target)
    nlon = len(lon_target)
    out = np.full((nlat, nlon), np.nan)

    lat_step = lat_target[1] - lat_target[0]
    lon_step = lon_target[1] - lon_target[0]

    lat_idx = np.rint((lat_vals - lat_target[0]) / lat_step).astype(int)
    lon_idx = np.rint((lon_vals - lon_target[0]) / lon_step).astype(int)

    good = (
        (lat_idx >= 0) & (lat_idx < nlat) &
        (lon_idx >= 0) & (lon_idx < nlon)
    )
    out[lat_idx[good], lon_idx[good]] = val_vals[good]
    return out


# -----------------------------
# Paths and model definitions
# -----------------------------
woa_root = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
woa_masks_root = "/data1/lorraine-20419/ocean_data/WOA23/masks"
hist_root = "/data1/lorraine-20419/benchmark_cmip6/historical"

model_cfgs = {
    "CESM2": {
        "y_dim": "nlat", "x_dim": "nlon",
        "lat_coord": "lat", "lon_coord": "lon",
        "depth_coord": "lev",
        "ocean_mask_from": "sftof",
    },
    "GFDL-ESM4": {
        "y_dim": "lat", "x_dim": "lon",
        "lat_coord": "lat", "lon_coord": "lon",
        "depth_coord": "lev",
    },
    "MPI-ESM1-2-LR": {
        "y_dim": "j", "x_dim": "i",
        "lat_coord": "latitude", "lon_coord": "longitude",
        "depth_coord": "lev",
    },
}
models = list(model_cfgs.keys())

# -----------------------------
# Read WOA23 1971-2000 annual 0 m temperature
# -----------------------------
woa_temp_file = os.path.join(
    woa_root, "temperature", "1971-2000_climate_normal",
    "woa23_decav71A0_t00_01.nc"
)
woa_ds = xr.open_dataset(woa_temp_file, decode_times=False)

# Prefer the objectively analysed field t_an.
if "t_an" in woa_ds.variables:
    woa_da = woa_ds["t_an"]
else:
    woa_da = None
    for v in woa_ds.data_vars:
        if woa_ds[v].attrs.get("standard_name", "") == "sea_water_temperature":
            woa_da = woa_ds[v]
            break
    if woa_da is None:
        woa_da = woa_ds[list(woa_ds.data_vars)[0]]

# Annual time slice and surface level.
if "time" in woa_da.dims:
    woa_da = woa_da.isel(time=0)
for dname in ["depth", "lev", "z"]:
    if dname in woa_da.dims:
        woa_da = woa_da.isel({dname: 0})
        break
woa_da = woa_da.squeeze()

# Coordinates.
if "lat" in woa_ds:
    lat_name = "lat"
elif "latitude" in woa_ds:
    lat_name = "latitude"
else:
    lat_name = list(woa_da.coords)[0]

if "lon" in woa_ds:
    lon_name = "lon"
elif "longitude" in woa_ds:
    lon_name = "longitude"
else:
    lon_name = list(woa_da.coords)[1]

lat_woa = woa_ds[lat_name].values.astype(float)
lon_woa = woa_ds[lon_name].values.astype(float)

# Put the array into (lat, lon) order.
if lat_name in woa_da.dims and lon_name in woa_da.dims:
    woa_da = woa_da.transpose(lat_name, lon_name)

woa_sst = to_numpy_float(woa_da)

# Ensure ascending lat/lon for the mask and plotting.
if lat_woa[0] > lat_woa[-1]:
    lat_woa = lat_woa[::-1]
    woa_sst = woa_sst[::-1, :]
if lon_woa[0] > lon_woa[-1]:
    lon_woa = lon_woa[::-1]
    woa_sst = woa_sst[:, ::-1]

# -----------------------------
# Read WOA23 land/sea mask
# -----------------------------
landsea_file = os.path.join(woa_masks_root, "landsea_01.msk")
mask_arr = read_landsea_mask(landsea_file, lat_woa, lon_woa)

woa_ocean = (mask_arr > 1) & np.isfinite(woa_sst)

# Area weights on the 1-degree WOA grid.
lat_rad = np.deg2rad(lat_woa)
w_woa = np.cos(lat_rad)[:, None] * np.ones((1, len(lon_woa)), dtype=float)

# Global mean over the whole WOA 1-degree ocean grid.
woa_global_mean = (
    np.nansum(woa_sst[woa_ocean] * w_woa[woa_ocean])
    / np.sum(w_woa[woa_ocean])
)
woa_global_mean_unweighted = np.nanmean(woa_sst[woa_ocean])

# -----------------------------
# Read each CMIP6 historical tos, compute native global mean,
# and nearest-neighbour regrid to the WOA 1-degree grid.
# -----------------------------
global_mean_models = {}
model_on_woa = {}
model_tos_clim = {}
native_ocean_masks = {}
native_lat2d = {}
native_lon2d = {}

for model, cfg in model_cfgs.items():
    # Locate tos files.
    tos_files = sorted(glob.glob(os.path.join(
        hist_root, model, "tos",
        f"tos_Omon_{model}_historical_*_*.nc"
    )))
    if len(tos_files) == 0:
        tos_files = sorted(glob.glob(os.path.join(hist_root, model, "tos", "*.nc")))
    if len(tos_files) == 0:
        raise FileNotFoundError(f"No tos files found for {model}")

    # Parse the actual variant from the filename.
    mvar = re.search(
        r"_historical_(r\d+i\d+p\d+f\d+)_",
        os.path.basename(tos_files[0])
    )
    variant = mvar.group(1) if mvar else None

    ds = xr.open_mfdataset(
        tos_files,
        combine="by_coords",
        use_cftime=True,
    )

    if "tos" in ds:
        tos_da = ds["tos"]
    else:
        tos_da = None
        for v in ds.data_vars:
            if ds[v].attrs.get("standard_name", "") == "sea_surface_temperature":
                tos_da = ds[v]
                break
        if tos_da is None:
            tos_da = ds[list(ds.data_vars)[0]]

    # Select 1971-01 through 2000-12 using the decoded cftime years.
    time_index = ds.indexes["time"]
    years = np.array([t.year for t in time_index])
    sel = (years >= 1971) & (years <= 2000)
    if not np.any(sel):
        raise ValueError(f"No 1971-2000 data found for {model}")

    tos_sel = tos_da.isel(time=np.where(sel)[0])
    tos_clim = tos_sel.mean(dim="time").compute()

    # CMIP6 tos is in kelvin.
    units = tos_da.attrs.get("units", "")
    if units.lower() in ("k", "kelvin", "degrees_kelvin", "degk", "degrees_k"):
        tos_clim = tos_clim - 273.15

    # Native lat/lon arrays.
    latv = ds[cfg["lat_coord"]].values.astype(float)
    lonv = ds[cfg["lon_coord"]].values.astype(float)
    if latv.ndim == 1:
        lon2d, lat2d = np.meshgrid(lonv, latv)
    else:
        lat2d = latv
        lon2d = lonv

    # Native cell areas.
    area_files = sorted(glob.glob(os.path.join(
        hist_root, model, "areacello",
        f"areacello_Ofx_{model}_historical_*_*.nc"
    )))
    if len(area_files) == 0:
        area_files = sorted(glob.glob(os.path.join(
            hist_root, model, "areacello", "*.nc"
        )))
    if len(area_files) == 0:
        raise FileNotFoundError(f"No areacello files found for {model}")

    if variant is not None:
        matched = [f for f in area_files if variant in os.path.basename(f)]
        if matched:
            area_files = matched

    area_ds = xr.open_dataset(area_files[0])
    area_arr = to_numpy_float(area_ds["areacello"].squeeze())

    target_shape = tos_clim.shape
    area_arr = align_2d(area_arr, target_shape)

    # Ocean area.  CESM2's raw areacello is global; apply sftof.
    if cfg.get("ocean_mask_from") == "sftof":
        sftof_files = sorted(glob.glob(os.path.join(
            hist_root, model, "sftof",
            f"sftof_Ofx_{model}_historical_*_*.nc"
        )))
        if len(sftof_files) == 0:
            sftof_files = sorted(glob.glob(os.path.join(
                hist_root, model, "sftof", "*.nc"
            )))

        if variant is not None and len(sftof_files) > 0:
            matched = [f for f in sftof_files if variant in os.path.basename(f)]
            if matched:
                sftof_files = matched

        if len(sftof_files) > 0:
            sf_ds = xr.open_dataset(sftof_files[0])
            sf_arr = to_numpy_float(sf_ds["sftof"].squeeze())
            sf_arr = align_2d(sf_arr, target_shape)
            if np.nanmax(sf_arr) > 1.5:
                sf_arr = sf_arr / 100.0
            area_ocean = area_arr * sf_arr
        else:
            area_ocean = area_arr.copy()
    else:
        area_ocean = area_arr.copy()

    area_ocean = np.where(
        np.isfinite(area_ocean) & (area_ocean > 0),
        area_ocean,
        0.0
    )

    tos_vals = to_numpy_float(tos_clim)
    if tos_vals.shape != target_shape:
        tos_vals = align_2d(tos_vals, target_shape)

    valid_native = np.isfinite(tos_vals) & (area_ocean > 0)

    # Native-grid global mean.
    global_mean_models[model] = (
        np.sum(tos_vals[valid_native] * area_ocean[valid_native])
        / np.sum(area_ocean[valid_native])
    )

    model_tos_clim[model] = tos_vals
    native_ocean_masks[model] = valid_native
    native_lat2d[model] = lat2d
    native_lon2d[model] = lon2d

    # Nearest-neighbour remap of the model climatology to the WOA grid.
    lat_flat = lat2d[valid_native]
    lon_flat = lon2d[valid_native]
    tos_flat = tos_vals[valid_native]

    xyz_native = lonlat_to_xyz(lat_flat, lon_flat)
    tree = cKDTree(xyz_native)

    lon_woa_mesh, lat_woa_mesh = np.meshgrid(lon_woa, lat_woa)
    xyz_woa = lonlat_to_xyz(lat_woa_mesh.ravel(), lon_woa_mesh.ravel())

    _, idx = tree.query(xyz_woa, k=1)
    model_regrid = tos_flat[idx].reshape(len(lat_woa), len(lon_woa))
    model_on_woa[model] = model_regrid

# -----------------------------
# Common comparison cells and statistics
# -----------------------------
common_mask = woa_ocean.copy()
for model in models:
    common_mask &= np.isfinite(model_on_woa[model])

# Unweighted statistics on the common WOA-grid cells.
comp_stats = {}
for model in models:
    model_vals = model_on_woa[model][common_mask]
    woa_vals = woa_sst[common_mask]
    diff = model_vals - woa_vals

    comp_stats[model] = {
        "mean_model": float(np.mean(model_vals)),
        "mean_woa": float(np.mean(woa_vals)),
        "mean_bias": float(np.mean(diff)),
        "rmse": float(np.sqrt(np.mean(diff ** 2))),
        "corr": float(np.corrcoef(model_vals, woa_vals)[0, 1]),
    }

# Area-weighted statistics on the same common cells, for reference.
comp_stats_weighted = {}
for model in models:
    model_vals = model_on_woa[model][common_mask]
    woa_vals = woa_sst[common_mask]
    diff = model_vals - woa_vals
    w = w_woa[common_mask]

    mean_model_w = np.sum(model_vals * w) / np.sum(w)
    mean_woa_w = np.sum(woa_vals * w) / np.sum(w)
    mean_bias_w = np.sum(diff * w) / np.sum(w)
    rmse_w = np.sqrt(np.sum((diff ** 2) * w) / np.sum(w))

    cov = np.sum(w * (model_vals - mean_model_w) * (woa_vals - mean_woa_w)) / np.sum(w)
    var_m = np.sum(w * (model_vals - mean_model_w) ** 2) / np.sum(w)
    var_w = np.sum(w * (woa_vals - mean_woa_w) ** 2) / np.sum(w)

    comp_stats_weighted[model] = {
        "mean_model": float(mean_model_w),
        "mean_woa": float(mean_woa_w),
        "mean_bias": float(mean_bias_w),
        "rmse": float(rmse_w),
        "corr": float(cov / np.sqrt(var_m * var_w)),
    }

ranking = sorted(models, key=lambda m: abs(comp_stats[m]["mean_bias"]))

# -----------------------------
# Print the requested numerical results
# -----------------------------
print("Global mean SST, 1971-2000 (native ocean grids, area-weighted):")
for model in models:
    print(f"  {model:16s}: {global_mean_models[model]:.2f} degC")
print(f"  {'WOA23':16s}: {woa_global_mean:.2f} degC")
print(f"  (WOA23 unweighted simple mean over ocean cells: {woa_global_mean_unweighted:.2f} degC)")

print("\nComparison on common WOA23 1-degree grid cells:")
print(f"  Number of common cells: {int(common_mask.sum())}")
for model in models:
    s = comp_stats[model]
    print(
        f"  {model:16s}: mean={s['mean_model']:.2f} degC, "
        f"bias={s['mean_bias']:.2f} degC, "
        f"RMSE={s['rmse']:.3f} degC, "
        f"corr={s['corr']:.3f}"
    )
print(f"  WOA23 mean on common cells: {comp_stats[models[0]]['mean_woa']:.2f} degC")

print("\nArea-weighted comparison statistics on common cells (for reference):")
for model in models:
    s = comp_stats_weighted[model]
    print(
        f"  {model:16s}: mean={s['mean_model']:.2f} degC, "
        f"bias={s['mean_bias']:.2f} degC, "
        f"RMSE={s['rmse']:.3f} degC, "
        f"corr={s['corr']:.3f}"
    )

print("\nRanking by absolute mean bias, smallest to largest:")
print("  " + " < ".join(ranking))

# Store maps of bias on the WOA grid for plotting.
model_bias_on_woa = {
    model: model_on_woa[model] - woa_sst
    for model in models
}
