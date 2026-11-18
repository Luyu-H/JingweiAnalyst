import os
import glob
import numpy as np
import pandas as pd
import xarray as xr

# ----------------------------------------------------------------------
# WOA23 reference: 1971-2000 temperature, annual, objectively analysed, 0 m
# ----------------------------------------------------------------------
woa_root = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
woa_file = os.path.join(
    woa_root,
    "temperature",
    "1971-2000_climate_normal",
    "woa23_decav71A0_t00_01.nc",
)

with xr.open_dataset(woa_file, decode_times=False) as ds_woa:
    varname = "t_an"
    if varname not in ds_woa:
        candidates = [v for v in ds_woa.data_vars if v.endswith("_an")]
        if not candidates:
            raise KeyError(f"No objectively analysed field found in {woa_file}")
        varname = candidates[0]

    da = ds_woa[varname]
    if "time" in da.dims:
        da = da.isel(time=0)
    if "depth" in da.dims:
        da = da.sel(depth=0.0, method="nearest")
    woa_sst = da.squeeze().load()

# Ensure ascending lat/lon for simple floor-based binning below
woa_sst = woa_sst.sortby("lat").sortby("lon")
woa_lat = woa_sst["lat"].values
woa_lon = woa_sst["lon"].values
woa_sst_vals = woa_sst.values.astype(float)

# ----------------------------------------------------------------------
# WOA23 land-sea mask: ocean where mask value > 1
# ----------------------------------------------------------------------
def read_woa_landsea_mask(path, target_lat, target_lon):
    rows = []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
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
            rows.append((lat, lon, val))

    df = pd.DataFrame(rows, columns=["lat", "lon", "val"])
    lon_vals = df["lon"].values.astype(float)
    # Registry says -180..180; make robust if a 0..360 variant ever appears.
    if np.nanmax(lon_vals) > 180.0:
        lon_vals = ((lon_vals + 180.0) % 360.0) - 180.0

    ocean = df["val"].values > 1.0
    i = np.floor(df["lat"].values.astype(float) + 90.0).astype(int)
    j = np.floor(lon_vals + 180.0).astype(int)

    ok = (i >= 0) & (i < 180) & (j >= 0) & (j < 360)
    mask_arr = np.zeros((180, 360), dtype=bool)
    mask_arr[i[ok], j[ok]] = ocean[ok]

    lat_grid = np.arange(-89.5, 90.0, 1.0)
    lon_grid = np.arange(-179.5, 180.0, 1.0)
    mask_da = xr.DataArray(
        mask_arr,
        coords={"lat": lat_grid, "lon": lon_grid},
        dims=("lat", "lon"),
    )
    mask_da = mask_da.reindex(lat=target_lat, lon=target_lon)
    return mask_da.fillna(False).astype(bool)


mask_file = "/data1/lorraine-20419/ocean_data/WOA23/masks/landsea_01.msk"
woa_ocean_mask = read_woa_landsea_mask(mask_file, woa_lat, woa_lon)
woa_valid = woa_ocean_mask.values & np.isfinite(woa_sst_vals)

# "Global mean over the whole 1-degree grid": simple mean over valid ocean cells.
woa_global_mean = float(np.nanmean(woa_sst_vals[woa_valid]))

# Also keep an area-weighted value for context, though the requested bar chart
# uses the simple 1-degree-grid mean above.
coslat = np.cos(np.deg2rad(woa_lat))
woa_area_weights = np.broadcast_to(coslat[:, None], woa_sst_vals.shape)
woa_global_mean_area_weighted = float(
    np.nansum(woa_sst_vals[woa_valid] * woa_area_weights[woa_valid])
    / np.nansum(woa_area_weights[woa_valid])
)

# ----------------------------------------------------------------------
# CMIP6 historical tos, 1971-2000
# ----------------------------------------------------------------------
cmip_root = "/data1/lorraine-20419/benchmark_cmip6/historical"
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]


def load_cmip6_tos_1971_2000(model):
    files = sorted(glob.glob(os.path.join(cmip_root, model, "tos", "*.nc")))
    if not files:
        raise FileNotFoundError(f"No tos files found for {model}")

    ds = xr.open_mfdataset(
        files,
        combine="by_coords",
        use_cftime=True,
        decode_times=True,
    )

    years = ds.time.dt.year.values
    idx = np.where((years >= 1971) & (years <= 2000))[0]
    if len(idx) == 0:
        raise ValueError(f"No tos data in 1971-2000 for {model}")

    ds = ds.isel(time=idx)
    tos = ds["tos"]
    tos_clim = tos.mean(dim="time", skipna=True)

    # Published cell areas
    area_files = sorted(glob.glob(os.path.join(cmip_root, model, "areacello", "*.nc")))
    if not area_files:
        raise FileNotFoundError(f"No areacello files found for {model}")
    if len(area_files) == 1:
        ds_area = xr.open_dataset(area_files[0])
    else:
        ds_area = xr.open_mfdataset(area_files, combine="by_coords")

    area = ds_area["areacello"].squeeze()
    if set(area.dims) == set(tos_clim.dims):
        area = area.transpose(*tos_clim.dims)
    else:
        area = area.transpose(*tos_clim.dims)

    # Ocean mask. CESM2's raw areacello is global; use sftof for CESM2.
    # GFDL-ESM4 and MPI-ESM1-2-LR publish already ocean-masked areacello.
    if model == "CESM2":
        sftof_files = sorted(glob.glob(os.path.join(cmip_root, model, "sftof", "*.nc")))
        if sftof_files:
            ds_sftof = xr.open_dataset(sftof_files[0])
            sftof = ds_sftof["sftof"].squeeze()
            if set(sftof.dims) == set(tos_clim.dims):
                sftof = sftof.transpose(*tos_clim.dims)
            ocean = sftof.values > 0.0
        else:
            ocean = area.values > 0.0
    else:
        ocean = area.values > 0.0

    tos_vals = tos_clim.values.astype(float)
    area_vals = area.values.astype(float)

    if tos_vals.shape != area_vals.shape:
        raise ValueError(f"Shape mismatch for {model}: tos {tos_vals.shape}, area {area_vals.shape}")

    valid = np.isfinite(tos_vals) & ocean & (area_vals > 0.0)
    if valid.any():
        global_mean = float(
            np.nansum(tos_vals[valid] * area_vals[valid])
            / np.nansum(area_vals[valid])
        )
    else:
        global_mean = np.nan

    # Native-grid latitude/longitude
    if model == "MPI-ESM1-2-LR":
        lat_name, lon_name = "latitude", "longitude"
    else:
        lat_name, lon_name = "lat", "lon"

    lat_vals = ds[lat_name].values
    lon_vals = ds[lon_name].values

    if lat_vals.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon_vals, lat_vals)
    else:
        lat2d = lat_vals
        lon2d = lon_vals

    if lat2d.shape != tos_vals.shape:
        if lat2d.T.shape == tos_vals.shape:
            lat2d = lat2d.T
            lon2d = lon2d.T
        else:
            raise ValueError(
                f"Lat/lon shape mismatch for {model}: "
                f"lat2d {lat2d.shape}, tos {tos_vals.shape}"
            )

    return {
        "model": model,
        "tos_clim": tos_vals,
        "area": area_vals,
        "ocean": ocean,
        "lat2d": lat2d,
        "lon2d": lon2d,
        "global_mean": global_mean,
    }


model_data = {m: load_cmip6_tos_1971_2000(m) for m in models}

# ----------------------------------------------------------------------
# Bin each native model climatology onto the WOA23 1-degree grid.
# The binning is area-weighted over all model ocean cells falling in each
# WOA grid cell.  This is a simple conservative-remapping approximation.
# ----------------------------------------------------------------------
def regrid_model_to_woa(md, target_lat, target_lon):
    field = md["tos_clim"]
    area = md["area"]
    lat2d = md["lat2d"]
    lon2d = md["lon2d"]
    ocean = md["ocean"]

    valid = np.isfinite(field) & ocean & (area > 0.0)

    f = field[valid].astype(float)
    a = area[valid].astype(float)
    la = lat2d[valid].astype(float)
    lo = lon2d[valid].astype(float)

    lo = ((lo + 180.0) % 360.0) - 180.0

    i = np.floor(la + 90.0).astype(int)
    j = np.floor(lo + 180.0).astype(int)
    i = np.clip(i, 0, len(target_lat) - 1)
    j = np.clip(j, 0, len(target_lon) - 1)

    num = np.zeros((len(target_lat), len(target_lon)), dtype=float)
    den = np.zeros_like(num)
    np.add.at(num, (i, j), f * a)
    np.add.at(den, (i, j), a)

    out = np.full_like(num, np.nan)
    good = den > 0.0
    out[good] = num[good] / den[good]
    return out


model_on_woa = {m: regrid_model_to_woa(model_data[m], woa_lat, woa_lon) for m in models}

# Common comparison cells: WOA ocean with valid WOA SST and all three models valid.
common_mask = woa_valid.copy()
for m in models:
    common_mask &= np.isfinite(model_on_woa[m])

woa_vals_on_common = woa_sst_vals[common_mask]
woa_mean_on_common = float(np.nanmean(woa_vals_on_common))

comparison_results = {}
for m in models:
    model_vals = model_on_woa[m][common_mask]
    bias = model_vals - woa_vals_on_common

    comparison_results[m] = {
        "global_mean_native_C": float(model_data[m]["global_mean"]),
        "mean_on_WOA_common_C": float(np.nanmean(model_vals)),
        "bias_C": float(np.nanmean(bias)),
        "rmse_C": float(np.sqrt(np.nanmean(bias ** 2))),
        "corr": float(np.corrcoef(model_vals, woa_vals_on_common)[0, 1]),
    }

ranking = sorted(models, key=lambda m: abs(comparison_results[m]["bias_C"]))

summary = pd.DataFrame(
    {
        m: {
            "global_mean_native_C": comparison_results[m]["global_mean_native_C"],
            "mean_on_WOA_common_C": comparison_results[m]["mean_on_WOA_common_C"],
            "bias_C": comparison_results[m]["bias_C"],
            "rmse_C": comparison_results[m]["rmse_C"],
            "corr": comparison_results[m]["corr"],
        }
        for m in models
    }
).T

summary["abs_bias_C"] = summary["bias_C"].abs()
summary = summary.sort_values("abs_bias_C")

print("WOA23 global mean (simple 1-degree-grid mean):", round(woa_global_mean, 2))
print("WOA23 global mean (cos-lat area-weighted):", round(woa_global_mean_area_weighted, 2))
print("WOA23 mean on comparison cells:", round(woa_mean_on_common, 2))
print("\nModel comparison summary:")
print(summary)
print("\nRanking from smallest to largest |bias|:", ranking)
