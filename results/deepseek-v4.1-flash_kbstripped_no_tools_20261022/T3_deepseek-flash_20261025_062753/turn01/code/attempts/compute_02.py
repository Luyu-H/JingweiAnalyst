import os
import glob
import numpy as np
import pandas as pd
import xarray as xr

# ============================================================
# Paths
# ============================================================
woa_root   = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
woa_masks  = "/data1/lorraine-20419/ocean_data/WOA23/masks"
cmip_root  = "/data1/lorraine-20419/benchmark_cmip6/historical"

# ============================================================
# 1. WOA23 temperature, 1971-2000 annual, objectively analysed
#    field, surface level.
# ============================================================
woa_file = os.path.join(
    woa_root, "temperature", "1971-2000_climate_normal",
    "woa23_decav71A0_t00_01.nc"
)
ds_woa = xr.open_dataset(woa_file, decode_times=False)

# Pick the "an" field (objectively analysed mean).
woa_vars = [v for v in ds_woa.data_vars if v.endswith("_an")]
if not woa_vars:
    woa_vars = [v for v in ds_woa.data_vars if ds_woa[v].ndim >= 3]
woa_var = woa_vars[0]
woa = ds_woa[woa_var]

for d in ("time", "depth"):
    if d in woa.dims:
        woa = woa.isel({d: 0})
woa = woa.squeeze()
if woa.dims != ("lat", "lon"):
    woa = woa.transpose("lat", "lon")

woa_lat = ds_woa["lat"].values.astype(float)
woa_lon = ds_woa["lon"].values.astype(float)
woa_sst = woa.values.astype(float)

woa_da = xr.DataArray(
    woa_sst, dims=("lat", "lon"),
    coords={"lat": woa_lat, "lon": woa_lon}
)

# ============================================================
# 2. WOA23 land-sea mask.
#    Ocean cells are mask > 1.  Bottom_Standard_Level is the
#    last column; 1 = land, value-1 = deepest standard level.
# ============================================================
mask_file = os.path.join(woa_masks, "landsea_01.msk")

with open(mask_file, "r") as f:
    lines = f.readlines()
header_idx = 0
for i, line in enumerate(lines):
    s = line.strip()
    if s and not s.startswith("#"):
        header_idx = i
        break

mask_df = pd.read_csv(
    mask_file,
    sep=",",
    skiprows=header_idx + 1,
    header=None
)
lat_m = mask_df.iloc[:, 0].values.astype(float)
lon_m = mask_df.iloc[:, 1].values.astype(float)
bottom = mask_df.iloc[:, -1].values.astype(float)

lat_u = np.sort(np.unique(lat_m))
lon_u = np.sort(np.unique(lon_m))
mask_arr = np.full((len(lat_u), len(lon_u)), np.nan)

lat_idx = np.searchsorted(lat_u, lat_m)
lon_idx = np.searchsorted(lon_u, lon_m)
mask_arr[lat_idx, lon_idx] = bottom

mask_da = xr.DataArray(
    mask_arr, dims=("lat", "lon"),
    coords={"lat": lat_u, "lon": lon_u}
)
mask_da = mask_da.reindex(lat=woa_lat, lon=woa_lon, method="nearest")
ocean_woa = (mask_da.values > 1)

# ============================================================
# 3. WOA23 global mean over its full 1-degree ocean domain,
#    using cos(latitude) cell-area weights.
# ============================================================
R_earth = 6371000.0  # m
cell_area_woa = (
    (R_earth * np.pi / 180.0) ** 2
    * np.cos(np.deg2rad(woa_lat))[:, None]
    * np.ones((1, len(woa_lon)))
)

valid_woa = ocean_woa & np.isfinite(woa_sst)
woa_gm = float(
    np.sum(woa_sst[valid_woa] * cell_area_woa[valid_woa])
    / np.sum(cell_area_woa[valid_woa])
)
woa_area_m2 = float(np.sum(cell_area_woa[valid_woa]))

# ============================================================
# 4. Helpers for CMIP6
# ============================================================
def get_time_indices(da):
    """Return a slice covering 1971-01 to 2000-12 inclusive."""
    time = da.indexes["time"]
    years = np.array([t.year for t in time])
    months = np.array([t.month for t in time])
    starts = np.where((years == 1971) & (months == 1))[0]
    stops = np.where((years == 2000) & (months == 12))[0]
    if len(starts) == 0 or len(stops) == 0:
        raise ValueError("1971-01 to 2000-12 not found in time axis")
    return slice(starts[0], stops[-1] + 1)


def conservative_remap_to_woa(field, area, lat2d, lon2d,
                              target_lat, target_lon, nsamp=8):
    """
    First-order conservative remap to a regular 1-degree grid by
    supersampling each source cell with nsamp x nsamp points and
    distributing its area to the target cell that contains each point.
    """
    field = np.asarray(field, dtype=float)
    area = np.asarray(area, dtype=float)
    ny, nx = field.shape

    if lat2d.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon2d, lat2d)

    # Make the longitude axis continuous before averaging to corners.
    lon_unwrapped = np.rad2deg(np.unwrap(np.deg2rad(lon2d), axis=1))

    lat_p = np.pad(lat2d, ((1, 1), (1, 1)), mode="edge")
    lon_p = np.pad(lon_unwrapped, ((1, 1), (1, 1)), mode="edge")

    lat_c = 0.25 * (
        lat_p[:-1, :-1] + lat_p[1:, :-1]
        + lat_p[:-1, 1:] + lat_p[1:, 1:]
    )
    lon_c = 0.25 * (
        lon_p[:-1, :-1] + lon_p[1:, :-1]
        + lon_p[:-1, 1:] + lon_p[1:, 1:]
    )

    nlat = len(target_lat)
    nlon = len(target_lon)
    target_sum = np.zeros((nlat, nlon), dtype=float)
    target_area = np.zeros((nlat, nlon), dtype=float)

    valid_src = np.isfinite(field) & np.isfinite(area) & (area > 0)

    us = (np.arange(nsamp) + 0.5) / nsamp
    vs = (np.arange(nsamp) + 0.5) / nsamp
    w_per_sample = area / (nsamp * nsamp)

    for u in us:
        for v in vs:
            slat = (
                (1 - u) * (1 - v) * lat_c[:-1, :-1]
                + u * (1 - v) * lat_c[:-1, 1:]
                + (1 - u) * v * lat_c[1:, :-1]
                + u * v * lat_c[1:, 1:]
            )
            slon = (
                (1 - u) * (1 - v) * lon_c[:-1, :-1]
                + u * (1 - v) * lon_c[:-1, 1:]
                + (1 - u) * v * lon_c[1:, :-1]
                + u * v * lon_c[1:, 1:]
            )

            slon = ((slon + 180.0) % 360.0) - 180.0
            slat = np.clip(slat, -90.0 + 1e-6, 90.0 - 1e-6)

            ilon = np.floor(slon - (target_lon[0] - 0.5)).astype(int)
            ilat = np.floor(slat - (target_lat[0] - 0.5)).astype(int)

            good = (
                valid_src
                & (ilat >= 0) & (ilat < nlat)
                & (ilon >= 0) & (ilon < nlon)
            )
            if not np.any(good):
                continue

            np.add.at(
                target_area,
                (ilat[good], ilon[good]),
                w_per_sample[good]
            )
            np.add.at(
                target_sum,
                (ilat[good], ilon[good]),
                field[good] * w_per_sample[good]
            )

    out = np.full((nlat, nlon), np.nan, dtype=float)
    filled = target_area > 0
    out[filled] = target_sum[filled] / target_area[filled]
    return out


def regrid_regular_to_woa(da, target_lat, target_lon):
    """Longitude-shift and interpolate a regular 1-degree field."""
    da = da.sortby("lat")
    da = da.assign_coords(
        lon=((da.lon + 180.0) % 360.0) - 180.0
    ).sortby("lon")
    return da.interp(lat=target_lat, lon=target_lon, method="linear")

# ============================================================
# 5. CMIP6 historical tos, 1971-2000 climatology, native means,
#    and conservative remapping to the WOA23 1-degree grid.
# ============================================================
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

native_means = {}
native_areas = {}
regridded = {}

for model in models:
    tos_dir = os.path.join(cmip_root, model, "tos")
    tos_files = sorted(
        glob.glob(os.path.join(tos_dir, f"tos_Omon_{model}_historical_*.nc"))
    )
    if len(tos_files) == 0:
        raise FileNotFoundError(f"No tos files for {model} under {tos_dir}")

    ds_tos = xr.open_mfdataset(
        tos_files, combine="by_coords", use_cftime=True
    )
    tos = ds_tos["tos"].sortby("time")
    tos = tos.isel(time=get_time_indices(tos))

    # Replace the CMIP6 missing value before any arithmetic.
    tos = tos.where(tos < 1e19)
    tos_clim = tos.mean("time").compute()

    # Native cell areas.
    area_dir = os.path.join(cmip_root, model, "areacello")
    area_files = sorted(
        glob.glob(os.path.join(area_dir, f"areacello_Ofx_{model}_historical_*.nc"))
    )
    if len(area_files) == 0:
        raise FileNotFoundError(f"No areacello files for {model} under {area_dir}")

    ds_area = xr.open_dataset(area_files[0])
    area = ds_area["areacello"].squeeze()
    if "time" in area.dims:
        area = area.isel(time=0)
    area = area.compute()

    f = np.asarray(tos_clim.values, dtype=float)
    a = np.asarray(area.values, dtype=float)

    if f.shape != a.shape:
        area = area.transpose(*tos_clim.dims)
        a = np.asarray(area.values, dtype=float)

    valid = np.isfinite(f) & np.isfinite(a) & (a > 0)
    native_means[model] = float(
        np.nansum(f[valid] * a[valid]) / np.nansum(a[valid])
    )
    native_areas[model] = float(np.nansum(a[valid]))

    # Native lat/lon coordinates.
    if model == "MPI-ESM1-2-LR":
        lat_name, lon_name = "latitude", "longitude"
    else:
        lat_name, lon_name = "lat", "lon"

    lat2d = ds_tos[lat_name].values
    lon2d = ds_tos[lon_name].values

    if model == "GFDL-ESM4":
        da = xr.DataArray(
            f, dims=("lat", "lon"),
            coords={
                "lat": ds_tos["lat"].values,
                "lon": ds_tos["lon"].values
            }
        )
        reg = regrid_regular_to_woa(da, woa_lat, woa_lon)
        regridded[model] = np.asarray(reg.values, dtype=float)
    else:
        regridded[model] = conservative_remap_to_woa(
            f, a, lat2d, lon2d, woa_lat, woa_lon, nsamp=8
        )

# ============================================================
# 6. Shared comparison mask on the WOA23 1-degree grid.
# ============================================================
shared = ocean_woa & np.isfinite(woa_sst)
for model in models:
    shared = shared & np.isfinite(regridded[model])

n_shared = int(shared.sum())

# cos(latitude) area weights on the regular 1-degree WOA grid.
w2d = (
    np.cos(np.deg2rad(woa_lat))[:, None]
    * np.ones((1, len(woa_lon)))
)
w_shared = w2d[shared]
w_sum = float(np.sum(w_shared))

woa_shared_mean = float(
    np.nansum(woa_sst[shared] * w_shared) / w_sum
)

shared_means = {"WOA23": woa_shared_mean}
biases = {}
rmse_vals = {}
corr_vals = {}

for model in models:
    x = regridded[model][shared]
    y = woa_sst[shared]

    mx = float(np.nansum(x * w_shared) / w_sum)
    my = woa_shared_mean
    shared_means[model] = mx
    biases[model] = mx - my

    diff = x - y
    rmse_vals[model] = float(
        np.sqrt(np.nansum(diff**2 * w_shared) / w_sum)
    )

    cov = float(np.nansum(w_shared * (x - mx) * (y - my)) / w_sum)
    vx = float(np.nansum(w_shared * (x - mx)**2) / w_sum)
    vy = float(np.nansum(w_shared * (y - my)**2) / w_sum)
    corr_vals[model] = cov / np.sqrt(vx * vy)

ranking = sorted(models, key=lambda m: abs(biases[m]))

# ============================================================
# 7. Summary tables (unrounded values kept in the dictionaries)
# ============================================================
native_summary = pd.DataFrame({
    "model": models + ["WOA23"],
    "native_global_mean_C": [native_means[m] for m in models] + [woa_gm],
    "native_ocean_area_m2": [native_areas[m] for m in models] + [woa_area_m2],
})

shared_summary = pd.DataFrame({
    "model": models,
    "shared_mean_C": [shared_means[m] for m in models],
    "bias_C": [biases[m] for m in models],
    "rmse_C": [rmse_vals[m] for m in models],
    "corr": [corr_vals[m] for m in models],
})

print("Native-grid global means (°C):")
for m in models:
    print(f"  {m:16s} {native_means[m]:8.4f}   area = {native_areas[m]:.4e} m2")
print(f"  {'WOA23':16s} {woa_gm:8.4f}   area = {woa_area_m2:.4e} m2")

print(f"\nShared comparison cells: {n_shared}")
print(f"WOA23 shared mean: {woa_shared_mean:.4f} °C")
for m in models:
    print(
        f"{m:16s} mean={shared_means[m]:8.4f}  "
        f"bias={biases[m]:8.4f}  RMSE={rmse_vals[m]:8.4f}  "
        f"corr={corr_vals[m]:8.4f}"
    )
print("\nRanking by |bias|:")
for m in ranking:
    print(f"  {m:16s} {abs(biases[m]):.4f} °C")
