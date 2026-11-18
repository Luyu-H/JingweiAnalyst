import os
import re
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
    for dname in ("depth", "lev", "deptht", "z", "st_ocean"):
        if dname in da.dims:
            da = da.sel({dname: 0.0}, method="nearest")
            break
    woa_sst = da.squeeze().load()

# normalise coordinate names to lat/lon
lat_name = "lat" if "lat" in woa_sst.dims else "latitude"
lon_name = "lon" if "lon" in woa_sst.dims else "longitude"
woa_sst = woa_sst.rename({lat_name: "lat", lon_name: "lon"})

# Ensure ascending lat/lon for simple floor-based binning below
woa_sst = woa_sst.sortby("lat").sortby("lon").transpose("lat", "lon")
woa_lat = np.asarray(woa_sst["lat"].values, dtype=float)
woa_lon = np.asarray(woa_sst["lon"].values, dtype=float)
woa_sst_vals = np.asarray(woa_sst.values, dtype=float)
print("WOA23 grid:", woa_sst_vals.shape, "lat", woa_lat[0], "->", woa_lat[-1],
      "lon", woa_lon[0], "->", woa_lon[-1])

# ----------------------------------------------------------------------
# WOA23 land-sea mask: ocean where mask value > 1
# ----------------------------------------------------------------------
_NUM_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?")


def _line_numbers(line):
    """Pull the numeric fields out of one mask-file line.

    Handles whitespace-, comma-, semicolon- or tab-separated columns as well
    as columns that are formatted run-together (Fortran fixed width), e.g.
    '-89.5-179.5    1'.  Returns [] when the line holds no usable numbers
    (header / comment / footer text).
    """
    s = line.strip()
    if not s:
        return []

    tokens = [t for t in re.split(r"[\s,;|]+", s) if t]
    if len(tokens) >= 3:
        try:
            nums = [float(t.replace("D", "E").replace("d", "e")) for t in tokens]
            if len(nums) >= 3:
                return nums
        except ValueError:
            pass

    # Fallback: extract every numeric literal from the line.
    return [float(v.replace("D", "E").replace("d", "e")) for v in _NUM_RE.findall(s)]


def read_woa_landsea_mask(path, target_lat, target_lon):
    """Rebuild the 180x360 WOA 1-degree mask from the .msk table.

    Returns an xr.DataArray of booleans on (target_lat, target_lon), or None if
    the table could not be parsed into a plausible mask.
    """
    rows = []
    with open(path, "r", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            nums = _line_numbers(s)
            if len(nums) < 3:
                continue
            rows.append((nums[0], nums[1], nums[2]))

    if not rows:
        print(f"WARNING: no land-sea mask records parsed from {path}")
        return None

    vals = np.asarray(rows, dtype=float)
    lat_m = vals[:, 0]
    lon_m = vals[:, 1]
    lev_m = vals[:, 2]

    # Registry says -180..180; make robust if a 0..360 variant ever appears.
    if np.nanmax(lon_m) > 180.0:
        lon_m = ((lon_m + 180.0) % 360.0) - 180.0

    ocean = lev_m > 1.0                      # 1 = land, (value-1) = deepest level
    i = np.floor(lat_m + 90.0).astype(int)   # -90..90  -> 0..179
    j = np.floor(lon_m + 180.0).astype(int)  # -180..180 -> 0..359

    ok = (i >= 0) & (i < 180) & (j >= 0) & (j < 360)
    mask_arr = np.zeros((180, 360), dtype=bool)
    mask_arr[i[ok], j[ok]] = ocean[ok]

    n_ocean = int(mask_arr.sum())
    print(f"WOA land-sea mask: {len(rows)} records, {n_ocean} ocean cells (expected ~41088)")
    if not (20000 < n_ocean < 60000):
        print("WARNING: implausible ocean-cell count; ignoring the .msk mask.")
        return None

    lat_grid = np.arange(-89.5, 90.0, 1.0)
    lon_grid = np.arange(-179.5, 180.0, 1.0)

    tlat = np.asarray(target_lat, dtype=float)
    tlon = np.asarray(target_lon, dtype=float)

    # Fast, exact path: target grid is the native WOA 1-degree grid.
    if (tlat.size == lat_grid.size and tlon.size == lon_grid.size
            and np.allclose(tlat, lat_grid) and np.allclose(tlon, lon_grid)):
        return xr.DataArray(mask_arr,
                            coords={"lat": tlat, "lon": tlon},
                            dims=("lat", "lon"))

    mask_da = xr.DataArray(mask_arr,
                           coords={"lat": lat_grid, "lon": lon_grid},
                           dims=("lat", "lon"))
    mask_da = mask_da.reindex(lat=tlat, lon=tlon)
    if bool(mask_da.isnull().all()):
        # Coordinate values did not match exactly (float32 vs float64 etc.)
        mask_da = xr.DataArray(mask_arr,
                               coords={"lat": tlat, "lon": tlon},
                               dims=("lat", "lon"))
    return mask_da.fillna(False).astype(bool)


mask_file = "/data1/lorraine-20419/ocean_data/WOA23/masks/landsea_01.msk"
woa_ocean_mask = read_woa_landsea_mask(mask_file, woa_lat, woa_lon)

woa_finite = np.isfinite(woa_sst_vals)
if woa_ocean_mask is None:
    woa_valid = woa_finite
else:
    woa_valid = np.asarray(woa_ocean_mask.values, dtype=bool) & woa_finite
print("WOA cells used for the global mean:", int(woa_valid.sum()))

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


def open_many(paths, use_cftime=True):
    """Open one or several NetCDF files, concatenating along time if needed."""
    if len(paths) == 1:
        return xr.open_dataset(paths[0], use_cftime=use_cftime, decode_times=True)
    return xr.open_mfdataset(
        paths,
        combine="by_coords",
        use_cftime=use_cftime,
        decode_times=True,
        data_vars="minimal",
        coords="minimal",
        compat="override",
    )


def align_like(da, target_dims, label):
    """Transpose da so its dimensions follow target_dims (must be the same set)."""
    if set(da.dims) != set(target_dims):
        raise ValueError(f"{label}: dims {da.dims} do not match tos dims {target_dims}")
    return da.transpose(*target_dims)


def load_cmip6_tos_1971_2000(model):
    files = sorted(glob.glob(os.path.join(cmip_root, model, "tos", "*.nc")))
    if not files:
        raise FileNotFoundError(f"No tos files found for {model}")

    ds = open_many(files)

    years = np.asarray(ds.time.dt.year.values, dtype=int)
    idx = np.where((years >= 1971) & (years <= 2000))[0]
    if len(idx) == 0:
        raise ValueError(f"No tos data in 1971-2000 for {model}")

    ds = ds.isel(time=idx)
    tos = ds["tos"]
    tos_clim = tos.mean(dim="time", skipna=True).squeeze()
    if tos_clim.ndim != 2:
        raise ValueError(f"{model}: unexpected tos climatology dims {tos_clim.dims}")

    # Published cell areas
    area_files = sorted(glob.glob(os.path.join(cmip_root, model, "areacello", "*.nc")))
    if not area_files:
        raise FileNotFoundError(f"No areacello files found for {model}")
    ds_area = open_many(area_files)
    area = ds_area["areacello"].squeeze()
    area = align_like(area, tos_clim.dims, f"{model} areacello")

    # Ocean mask. CESM2's raw areacello is global; use sftof for CESM2.
    # GFDL-ESM4 and MPI-ESM1-2-LR publish already ocean-masked areacello.
    if model == "CESM2":
        sftof_files = sorted(glob.glob(os.path.join(cmip_root, model, "sftof", "*.nc")))
        if sftof_files:
            ds_sftof = open_many(sftof_files)
            sftof = ds_sftof["sftof"].squeeze()
            sftof = align_like(sftof, tos_clim.dims, f"{model} sftof")
            ocean = np.asarray(sftof.values, dtype=float) > 0.0
        else:
            print(f"WARNING: no sftof for {model}; using areacello > 0 as ocean mask.")
            ocean = np.asarray(area.values, dtype=float) > 0.0
    else:
        ocean = np.asarray(area.values, dtype=float) > 0.0

    tos_vals = np.asarray(tos_clim.values, dtype=float)
    area_vals = np.asarray(area.values, dtype=float)

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

    lat_vals = np.asarray(ds[lat_name].values, dtype=float)
    lon_vals = np.asarray(ds[lon_name].values, dtype=float)

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
for m in models:
    print(f"{m}: native global mean SST 1971-2000 = {model_data[m]['global_mean']:.4f} degC")

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
print("Common comparison cells:", int(common_mask.sum()))

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

# ----------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    bias_on_common = {}
    for m in models:
        b = np.full(woa_sst_vals.shape, np.nan)
        b[common_mask] = model_on_woa[m][common_mask] - woa_sst_vals[common_mask]
        bias_on_common[m] = b

    ref_on_common = np.full(woa_sst_vals.shape, np.nan)
    ref_on_common[common_mask] = woa_sst_vals[common_mask]

    finite_biases = np.concatenate(
        [bias_on_common[m][np.isfinite(bias_on_common[m])] for m in models]
    )
    vmax = float(np.nanmax(np.abs(finite_biases))) if finite_biases.size else 1.0
    vmax = max(vmax, 0.1)

    # ---- Figure 1: 2x2 global maps -------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    axes = axes.ravel()

    for ax, m in zip(axes[:3], models):
        pm = ax.pcolormesh(woa_lon, woa_lat, bias_on_common[m], shading="auto",
                           cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_title(f"{m} - WOA23 SST bias (1971-2000)")
        cb = fig.colorbar(pm, ax=ax, extend="both")
        cb.set_label("SST bias (°C)")

    pm = axes[3].pcolormesh(woa_lon, woa_lat, ref_on_common, shading="auto",
                            cmap="viridis")
    axes[3].set_title("WOA23 0 m SST reference (comparison cells)")
    cb = fig.colorbar(pm, ax=axes[3])
    cb.set_label("SST (°C)")

    for ax in axes:
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)

    fig.suptitle("SST 1971-2000: CMIP6 model bias vs WOA23 and WOA23 reference")
    fig.savefig("fig1_sst_bias_maps.png", dpi=150)
    plt.close(fig)

    # ---- Figure 2: bar chart of global means ---------------------------
    labels = models + ["WOA23"]
    values = [comparison_results[m]["global_mean_native_C"] for m in models] + [woa_global_mean]

    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    bars = ax.bar(labels, values, color=["#4C72B0", "#55A868", "#C44E52", "#8172B2"])
    for bar, v in zip(bars, values):
        ax.annotate(f"{v:.2f}",
                    xy=(bar.get_x() + bar.get_width() / 2.0, bar.get_height()),
                    xytext=(0, 4), textcoords="offset points",
                    ha="center", va="bottom", fontsize=10)
    ax.set_ylabel("Global mean SST (°C)")
    ax.set_title("1971-2000 global mean SST (native ocean grids)")
    ax.set_ylim(min(values) - 3.0 if np.isfinite(min(values)) else 0.0,
                max(values) + 3.0 if np.isfinite(max(values)) else 1.0)
    fig.savefig("fig2_global_mean_sst.png", dpi=150)
    plt.close(fig)

    print("\nSaved figures: fig1_sst_bias_maps.png, fig2_global_mean_sst.png")
except Exception as exc:  # plotting must not destroy the numerical results
    print(f"WARNING: figure generation failed: {type(exc).__name__}: {exc}")

# ----------------------------------------------------------------------
# Rounded report
# ----------------------------------------------------------------------
print("\n--- Rounded results ---")
for m in models:
    r = comparison_results[m]
    print(f"{m}: native global mean = {r['global_mean_native_C']:.2f} °C; "
          f"mean on WOA cells = {r['mean_on_WOA_common_C']:.2f} °C; "
          f"bias = {r['bias_C']:.2f} °C; RMSE = {r['rmse_C']:.3f} °C; "
          f"corr = {r['corr']:.3f}")
print(f"WOA23: global mean (1-degree grid) = {woa_global_mean:.2f} °C; "
      f"mean on comparison cells = {woa_mean_on_common:.2f} °C")
print("Ranking (smallest to largest |bias|):", ranking)
