import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

CMIP_ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"
WOA_ROOT = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
WOA_MASK_ROOT = "/data1/lorraine-20419/ocean_data/WOA23/masks"

models = {
    "CESM2": {
        "y_dim": "nlat",
        "x_dim": "nlon",
        "lat_coord": "lat",
        "lon_coord": "lon",
    },
    "GFDL-ESM4": {
        "y_dim": "lat",
        "x_dim": "lon",
        "lat_coord": "lat",
        "lon_coord": "lon",
    },
    "MPI-ESM1-2-LR": {
        "y_dim": "j",
        "x_dim": "i",
        "lat_coord": "latitude",
        "lon_coord": "longitude",
    },
}


def open_cmip_paths(paths):
    """Open one or more CMIP6 netCDF files and concatenate along time if needed."""
    dss = [xr.open_dataset(p, use_cftime=True, decode_times=True) for p in paths]
    if len(dss) == 1:
        return dss[0]
    return xr.concat(
        dss,
        dim="time",
        combine_attrs="override",
        coords="minimal",
        compat="override",
    )


def get_var(ds, name):
    """Return the requested data variable, with a small fallback for renamed files."""
    if name in ds.data_vars:
        return ds[name]
    for v in ds.data_vars:
        if name.lower() in v.lower():
            return ds[v]
    return ds[list(ds.data_vars)[0]]


def mask_fill(da):
    """Convert CF fill/missing values and the common 1e20 sentinel to NaN."""
    out = da
    for att in ("_FillValue", "missing_value"):
        if att in out.attrs and out.attrs[att] is not None:
            fv = np.asarray(out.attrs[att], dtype=float).ravel()[0]
            out = out.where(out != fv)
    out = out.where(out < 1.0e19)
    return out


def to_celsius(da):
    """Convert a temperature DataArray to degrees Celsius if needed."""
    units = str(da.attrs.get("units", "")).lower()
    if units in ("k", "kelvin", "kelvins", "degk", "degrees_k") or "kelvin" in units:
        da = da - 273.15
        da.attrs = dict(da.attrs)
        da.attrs["units"] = "degrees_Celsius"
    return da


def select_years(da, y0, y1):
    """Select time steps whose calendar year is between y0 and y1 inclusive."""
    yrs = da["time"].dt.year
    idx = np.where((yrs >= y0) & (yrs <= y1))[0]
    return da.isel(time=idx)


def get_latlon_2d(da, lat_name, lon_name):
    """Return 2-D latitude/longitude fields, broadcasting 1-D coordinates when needed."""
    lat = da[lat_name]
    lon = da[lon_name]
    if lat.ndim == 2 and lon.ndim == 2:
        return lat, lon
    lat2d, lon2d = xr.broadcast(lat, lon)
    return lat2d, lon2d


def conservative_bin_to_woa(field, area, lat2d, lon2d, target_lat, target_lon, name=None):
    """
    Area-weighted conservative coarsening/binning to the WOA 1-degree grid.

    Each native grid cell is assigned to the 1-degree target cell containing its
    centre; target values are the areacello-weighted mean of the contributing
    native cells.  This is the conservative coarsening used for curvilinear grids.
    """
    area_b = area.broadcast_like(field)
    lat_b = lat2d.broadcast_like(field)
    lon_b = lon2d.broadcast_like(field)

    f = np.asarray(field.values, dtype=float).ravel()
    a = np.asarray(area_b.values, dtype=float).ravel()
    la = np.asarray(lat_b.values, dtype=float).ravel()
    lo = np.asarray(lon_b.values, dtype=float).ravel()

    good = (
        np.isfinite(f)
        & np.isfinite(a)
        & (a > 0)
        & np.isfinite(la)
        & np.isfinite(lo)
    )
    f, a, la, lo = f[good], a[good], la[good], lo[good]

    nlat, nlon = len(target_lat), len(target_lon)

    iy = np.floor(la + 90.0).astype(int)
    lo180 = (lo + 180.0) % 360.0 - 180.0
    ix = np.floor(lo180 + 180.0).astype(int)

    iy = np.clip(iy, 0, nlat - 1)
    ix = np.clip(ix, 0, nlon - 1)
    idx = iy * nlon + ix

    sum_a = np.bincount(idx, weights=a, minlength=nlat * nlon)
    sum_fa = np.bincount(idx, weights=f * a, minlength=nlat * nlon)

    out = np.full(nlat * nlon, np.nan, dtype=float)
    nz = sum_a > 0
    out[nz] = sum_fa[nz] / sum_a[nz]

    return xr.DataArray(
        out.reshape(nlat, nlon),
        dims=("lat", "lon"),
        coords={"lat": target_lat, "lon": target_lon},
        name=name if name is not None else getattr(field, "name", "regridded"),
    )


def read_woa_landsea(path):
    """Read the WOA land-sea mask table and rebuild the 180x360 grid from lat/lon rows."""
    rows = []
    with open(path, "r") as fh:
        for line in fh:
            s = line.split("#", 1)[0].strip()
            if not s or s[0] in "%!":
                continue
            parts = s.replace(",", " ").replace(";", " ").split()
            if not parts:
                continue
            try:
                vals = [float(p.replace("D", "E").replace("d", "e")) for p in parts]
            except ValueError:
                continue
            if len(vals) >= 3:
                rows.append(vals)

    if not rows:
        raise ValueError(f"No numeric rows parsed from WOA land-sea mask: {path}")

    # Keep the most common row length; data rows should all have the same length.
    lengths = [len(r) for r in rows]
    ncol = max(set(lengths), key=lengths.count)
    rows = [r for r in rows if len(r) == ncol]
    arr = np.asarray(rows, dtype=float)

    if arr.ndim == 1:
        arr = arr.reshape(1, -1)

    ncols = arr.shape[1]
    col_min = np.nanmin(arr, axis=0)
    col_max = np.nanmax(arr, axis=0)
    col_range = col_max - col_min

    # Identify the mask value column first: positive integer-like values.
    val_col = None
    for j in range(ncols):
        col = arr[:, j]
        finite = col[np.isfinite(col)]
        if finite.size == 0:
            continue
        if (
            np.nanmin(finite) >= 1.0
            and np.nanmax(finite) <= 200.0
            and np.allclose(finite, np.round(finite))
        ):
            val_col = j
            break

    # Identify latitude.
    lat_col = None
    for j in range(ncols):
        if j == val_col:
            continue
        if col_min[j] >= -90.0 and col_max[j] <= 90.0 and col_range[j] > 1.0:
            lat_col = j
            break
    if lat_col is None:
        lat_col = 0

    # Identify longitude.
    lon_col = None
    for j in range(ncols):
        if j in (lat_col, val_col):
            continue
        if col_min[j] >= -180.0 and col_max[j] <= 180.0 and col_range[j] > 1.0:
            if lon_col is None or (col_max[j] > 90.0 or col_min[j] < -90.0):
                lon_col = j
    if lon_col is None:
        lon_col = next((j for j in range(ncols) if j not in (lat_col, val_col)), 1)

    if val_col is None:
        val_col = next((j for j in range(ncols) if j not in (lat_col, lon_col)), ncols - 1)

    lat = arr[:, lat_col]
    lon = arr[:, lon_col]
    val = arr[:, val_col]

    # Normalise longitude to [-180, 180).
    lon = (lon + 180.0) % 360.0 - 180.0

    nlat, nlon = 180, 360
    mask = np.full((nlat, nlon), np.nan, dtype=float)

    iy = np.rint(lat + 89.5).astype(int)
    ix = np.rint(lon + 179.5).astype(int)
    good = (
        np.isfinite(lat)
        & np.isfinite(lon)
        & np.isfinite(val)
        & (iy >= 0)
        & (iy < nlat)
        & (ix >= 0)
        & (ix < nlon)
    )
    mask[iy[good], ix[good]] = val[good]

    return xr.DataArray(
        mask,
        dims=("lat", "lon"),
        coords={
            "lat": np.arange(-89.5, 90.0, 1.0),
            "lon": np.arange(-179.5, 180.0, 1.0),
        },
        name="landsea",
    )


def open_woa_temperature_0m():
    """Open WOA23 temperature, 1971-2000 annual, objectively analysed, at 0 m."""
    subdir = "temperature/1971-2000_climate_normal"
    filename = "woa23_decav71A0_t00_01.nc"
    path = os.path.join(WOA_ROOT, subdir, filename)

    ds = xr.open_dataset(path, decode_times=False)

    data_vars = [
        v for v in ds.data_vars
        if "bnd" not in v.lower() and "bounds" not in v.lower()
    ]

    var = None
    for cand in ("t_an", "temperature", "temp"):
        if cand in data_vars:
            var = cand
            break
    if var is None:
        for v in data_vars:
            if v.lower().startswith("t") or "_an" in v.lower():
                var = v
                break
    if var is None:
        var = data_vars[0]

    da = ds[var]

    if "time" in da.dims:
        da = da.isel(time=0)

    depth_dim = None
    for cand in ("depth", "lev", "z", "deptht"):
        if cand in da.dims:
            depth_dim = cand
            break
    if depth_dim is not None:
        da = da.sel({depth_dim: 0.0}, method="nearest")

    lat_name = None
    for cand in ("lat", "latitude", "y"):
        if cand in da.coords or cand in da.dims:
            lat_name = cand
            break

    lon_name = None
    for cand in ("lon", "longitude", "x"):
        if cand in da.coords or cand in da.dims:
            lon_name = cand
            break

    if lat_name != "lat" or lon_name != "lon":
        da = da.rename({lat_name: "lat", lon_name: "lon"})

    da = da.squeeze(drop=True)
    da = da.where(da < 1.0e30)
    return da.compute(), ds


# ---------------------------------------------------------------------------
# WOA23 reference field and mask
# ---------------------------------------------------------------------------
woa_temp, woa_ds = open_woa_temperature_0m()
landsea = read_woa_landsea(os.path.join(WOA_MASK_ROOT, "landsea_01.msk"))

if not (landsea.lat.equals(woa_temp.lat) and landsea.lon.equals(woa_temp.lon)):
    landsea = landsea.interp(lat=woa_temp.lat, lon=woa_temp.lon, method="nearest")

ocean_mask = landsea > 1

# WOA23 global mean over the whole 1-degree grid, lat-area weighted over valid ocean cells
woa_valid = np.isfinite(woa_temp) & ocean_mask
wlat = xr.DataArray(
    np.cos(np.deg2rad(woa_temp.lat.values)),
    dims=("lat",),
    coords={"lat": woa_temp.lat},
)
w2d = wlat.broadcast_like(woa_temp)

woa_global_mean = float(
    (
        (woa_temp * w2d).where(woa_valid).sum()
        / w2d.where(woa_valid).sum()
    ).compute()
)

# ---------------------------------------------------------------------------
# CMIP6 tos 1971-2000 climatologies and native-grid means
# ---------------------------------------------------------------------------
model_native_means = {}
model_regridded = {}
tos_clims = {}

target_lat = woa_temp.lat.values
target_lon = woa_temp.lon.values

for model, cfg in models.items():
    tos_files = sorted(glob.glob(os.path.join(CMIP_ROOT, model, "tos", "*.nc")))
    if not tos_files:
        raise FileNotFoundError(f"No tos files found for {model}")

    ds_tos = open_cmip_paths(tos_files)
    tos = mask_fill(get_var(ds_tos, "tos"))
    tos = to_celsius(tos)
    tos_sel = select_years(tos, 1971, 2000)
    tos_clim = tos_sel.mean(dim="time", skipna=True)
    tos_clims[model] = tos_clim

    area_files = sorted(glob.glob(os.path.join(CMIP_ROOT, model, "areacello", "*.nc")))
    if not area_files:
        raise FileNotFoundError(f"No areacello files found for {model}")

    ds_area = xr.open_dataset(area_files[0], use_cftime=True, decode_times=True)
    area = mask_fill(get_var(ds_area, "areacello"))
    if "time" in area.dims:
        area = area.isel(time=0)

    valid = np.isfinite(tos_clim)
    area_use = area.where(valid & (area > 0))

    num = (tos_clim * area_use).sum(
        dim=(cfg["y_dim"], cfg["x_dim"]), skipna=True
    )
    den = area_use.sum(dim=(cfg["y_dim"], cfg["x_dim"]), skipna=True)
    model_native_means[model] = float((num / den).compute())

    lat2d, lon2d = get_latlon_2d(
        tos_clim, cfg["lat_coord"], cfg["lon_coord"]
    )
    model_regridded[model] = conservative_bin_to_woa(
        tos_clim,
        area,
        lat2d,
        lon2d,
        target_lat,
        target_lon,
        name=model,
    )

# ---------------------------------------------------------------------------
# Shared comparison mask and statistics on the WOA23 1-degree grid
# ---------------------------------------------------------------------------
shared = ocean_mask & np.isfinite(woa_temp)
for model in models:
    shared = shared & np.isfinite(model_regridded[model])

comparison_mask = shared
mask_cell_count = int(comparison_mask.values.sum())

w_shared = w2d.where(comparison_mask)
woa_shared = woa_temp.where(comparison_mask)
woa_shared_mean = float((woa_shared * w_shared).sum() / w_shared.sum())

shared_means = {}
bias_means = {}
rmse = {}
corr = {}
bias_maps = {}

for model in models:
    mod = model_regridded[model]
    mod_shared = mod.where(comparison_mask)

    shared_means[model] = float(
        (mod_shared * w_shared).sum() / w_shared.sum()
    )

    delta_da = (mod - woa_temp).where(comparison_mask)
    bias_means[model] = float(
        (delta_da * w_shared).sum() / w_shared.sum()
    )
    bias_maps[model] = delta_da

    rmse_da = ((mod - woa_temp) ** 2).where(comparison_mask)
    rmse[model] = float(
        np.sqrt((rmse_da * w_shared).sum() / w_shared.sum())
    )

    xv = mod_shared.values.ravel()
    yv = woa_shared.values.ravel()
    wv = w_shared.values.ravel()

    good = (
        np.isfinite(xv)
        & np.isfinite(yv)
        & np.isfinite(wv)
        & (wv > 0)
    )
    xv, yv, wv = xv[good], yv[good], wv[good]

    wsum = wv.sum()
    mx = (wv * xv).sum() / wsum
    my = (wv * yv).sum() / wsum
    cov = (wv * (xv - mx) * (yv - my)).sum() / wsum
    vx = (wv * (xv - mx) ** 2).sum() / wsum
    vy = (wv * (yv - my) ** 2).sum() / wsum
    corr[model] = float(cov / np.sqrt(vx * vy))

ranking = sorted(list(models.keys()), key=lambda m: abs(bias_means[m]))

print("Native-grid 1971-2000 mean SST (deg C):")
for m in models:
    print(f"  {m}: {model_native_means[m]:.2f}")

print(f"WOA23 0 m whole-grid mean SST (deg C): {woa_global_mean:.2f}")
print(f"Shared comparison mask cells: {mask_cell_count}")
print(f"WOA23 shared-mask mean SST (deg C): {woa_shared_mean:.2f}")

print("\nShared-mask comparison:")
for m in models:
    print(
        f"  {m}: mean={shared_means[m]:.2f} deg C, "
        f"bias={bias_means[m]:.2f} deg C, "
        f"RMSE={rmse[m]:.3f} deg C, "
        f"corr={corr[m]:.3f}"
    )

print("\nRanking by |delta_SST| (smallest to largest):")
for i, m in enumerate(ranking, 1):
    print(f"  {i}. {m}: {bias_means[m]:+.2f} deg C")

# ---------------------------------------------------------------------------
# Figure 1: bias maps and WOA23 reference
# ---------------------------------------------------------------------------
bias_absmax = max(
    float(np.nanmax(np.abs(bias_maps[m].values))) for m in models
)

fig1, axes = plt.subplots(2, 2, figsize=(14, 10))
ax_bias = axes.ravel()[:3]

for ax, m in zip(ax_bias, models.keys()):
    im = ax.pcolormesh(
        target_lon,
        target_lat,
        bias_maps[m].values,
        shading="auto",
        cmap="RdBu_r",
        vmin=-bias_absmax,
        vmax=bias_absmax,
    )
    ax.set_title(f"{m} SST bias (model - WOA23)")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    cb = fig1.colorbar(im, ax=ax, orientation="horizontal", pad=0.06)
    cb.set_label("SST bias (°C)")

ax = axes.ravel()[3]
im = ax.pcolormesh(
    target_lon,
    target_lat,
    woa_shared.values,
    shading="auto",
    cmap="viridis",
)
ax.set_title("WOA23 0 m SST reference")
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
cb = fig1.colorbar(im, ax=ax, orientation="horizontal", pad=0.06)
cb.set_label("SST (°C)")

fig1.suptitle("1971-2000 SST: CMIP6 biases and WOA23 reference on shared comparison cells")
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()

# ---------------------------------------------------------------------------
# Figure 2: native-grid global mean SST bar chart
# ---------------------------------------------------------------------------
bar_labels = list(models.keys()) + ["WOA23"]
bar_values = [model_native_means[m] for m in models] + [woa_global_mean]

fig2, ax2 = plt.subplots(figsize=(8, 5))
bars = ax2.bar(bar_labels, bar_values, color=["C0", "C1", "C2", "0.5"])
for bar, v in zip(bars, bar_values):
    ax2.text(
        bar.get_x() + bar.get_width() / 2.0,
        bar.get_height(),
        f"{v:.2f}",
        ha="center",
        va="bottom",
        fontsize=9,
    )
ax2.set_ylabel("SST (°C)")
ax2.set_title("1971-2000 mean SST over native ocean grids")
ax2.set_ylim(min(bar_values) - 1.0, max(bar_values) + 1.0)
plt.tight_layout()
plt.show()
