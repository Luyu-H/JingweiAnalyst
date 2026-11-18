import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import griddata

WOA_ROOT  = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
MASK_ROOT = "/data1/lorraine-20419/ocean_data/WOA23/masks"
CMIP_ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"
MODELS = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

# ----------------------------------------------------------------------
# 1) WOA23 temperature, 1971-2000 normal, annual slice, objectively analysed field, 0 m
# ----------------------------------------------------------------------
woa_dir  = os.path.join(WOA_ROOT, "temperature", "1971-2000_climate_normal")
woa_path = os.path.join(woa_dir, "woa23_decav71A0_t00_01.nc")
if not os.path.exists(woa_path):
    cand = sorted(glob.glob(os.path.join(woa_dir, "woa23_*_t00_01.nc")))
    woa_path = cand[0]

woa_ds = xr.open_dataset(woa_path, decode_times=False)

an_vars = [v for v in woa_ds.data_vars if v.endswith("_an")]
vname   = an_vars[0] if an_vars else list(woa_ds.data_vars)[0]
da      = woa_ds[vname]

def _find_dim(darray, keys):
    for d in darray.dims:
        if any(k in d.lower() for k in keys):
            return d
    return None

latn = _find_dim(da, ["lat"])
lonn = _find_dim(da, ["lon"])
depn = _find_dim(da, ["depth", "zlev", "lev", "zt"])
timn = _find_dim(da, ["time"])

isel = {}
if depn is not None:
    isel[depn] = 0            # 0 m surface level
if timn is not None:
    isel[timn] = 0            # annual time slice (nn = 00)
woa_surf = da.isel(isel).squeeze(drop=True)

woa_lat = np.asarray(woa_ds[latn].values, dtype=float)
woa_lon = np.asarray(woa_ds[lonn].values, dtype=float)
woa_arr = np.asarray(woa_surf.values, dtype=float).squeeze()

# guard against un-decoded fill values
woa_arr = np.where(np.abs(woa_arr) > 1e10, np.nan, woa_arr)

# enforce ascending coordinate order
if woa_lat[0] > woa_lat[-1]:
    woa_lat = woa_lat[::-1]; woa_arr = woa_arr[::-1, :]
if woa_lon[0] > woa_lon[-1]:
    woa_lon = woa_lon[::-1]; woa_arr = woa_arr[:, ::-1]

# ----------------------------------------------------------------------
# 2) WOA23 1-degree land-sea mask (ocean where value > 1)
# ----------------------------------------------------------------------
mask_path = os.path.join(MASK_ROOT, "landsea_01.msk")
rows = []
with open(mask_path) as fh:
    for line in fh:
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
        except ValueError:
            continue
mask_df = pd.DataFrame(rows, columns=["lat", "lon", "bottom"])

dlat = woa_lat[1] - woa_lat[0]
dlon = woa_lon[1] - woa_lon[0]
ilat = np.rint((mask_df["lat"].values - woa_lat[0]) / dlat).astype(int)
ilon = np.rint((mask_df["lon"].values - woa_lon[0]) / dlon).astype(int)
ok   = (ilat >= 0) & (ilat < woa_lat.size) & (ilon >= 0) & (ilon < woa_lon.size)

mask_arr = np.full((woa_lat.size, woa_lon.size), np.nan)
mask_arr[ilat[ok], ilon[ok]] = mask_df["bottom"].values[ok].astype(float)

woa_ocean = (mask_arr > 1) & np.isfinite(woa_arr)

# global mean of the WOA field over the whole 1-degree grid (all ocean cells)
woa_global_mean = float(np.nanmean(woa_arr[woa_ocean]))

# (also recorded: cosine-latitude area weighted version, for reference)
_wlat  = np.cos(np.deg2rad(woa_lat))[:, None] * np.ones((1, woa_lon.size))
woa_global_mean_areawt = float(np.nansum(woa_arr * _wlat * woa_ocean) /
                               np.nansum(_wlat * woa_ocean))

# ----------------------------------------------------------------------
# 3) CMIP6 historical tos, 1971-2000 climatology on each native grid
# ----------------------------------------------------------------------
def files_in(model, name):
    d = os.path.join(CMIP_ROOT, model, name)
    fs = sorted(glob.glob(os.path.join(d, name + "*.nc")))
    if not fs:
        fs = sorted(glob.glob(os.path.join(d, "*.nc")))
    return fs

def open_var(model, name):
    fs = files_in(model, name)
    if not fs:
        return None, None
    if len(fs) == 1:
        ds = xr.open_dataset(fs[0], use_cftime=True)
    else:
        ds = xr.open_mfdataset(fs, combine="by_coords", use_cftime=True)
    return ds, fs

native_global_mean = {}
tos_on_woa         = {}
regrid_info        = {}

for model in MODELS:
    ds, _ = open_var(model, "tos")
    tos = ds["tos"]
    if "height" in tos.dims:
        tos = tos.isel(height=0)

    years = tos.time.dt.year
    keep  = np.where(((years >= 1971) & (years <= 2000)).values)[0]
    tos   = tos.isel(time=keep)

    # climatological mean (all months 1971-2000)
    tos_clim = tos.mean("time").compute()

    # native cell areas
    ads, _ = open_var(model, "areacello")
    area = ads["areacello"].squeeze().compute()

    # ocean/land mask: CESM2 needs sftof, the others publish pre-masked areas
    sds, _ = open_var(model, "sftof")
    if sds is not None:
        sftof = sds["sftof"].squeeze().compute()
        sea   = (np.asarray(sftof.values, dtype=float) > 0) & np.isfinite(np.asarray(sftof.values, dtype=float))
        sds.close()
    else:
        sea = True

    tvals = np.asarray(tos_clim.values, dtype=float).squeeze()
    avals = np.asarray(area.values,     dtype=float).squeeze()
    sea   = np.broadcast_to(sea, tvals.shape)
    valid = np.isfinite(tvals) & np.isfinite(avals) & (avals > 0) & sea

    native_global_mean[model] = float(np.nansum(tvals * avals * valid) /
                                      np.nansum(avals * valid))

    # 2-D lat/lon of the native grid
    latname = ds["tos"].dims[1]
    lonname = ds["tos"].dims[2]
    lat_c = "lat" if "lat" in ds.variables else ("latitude" if "latitude" in ds.variables else latname)
    lon_c = "lon" if "lon" in ds.variables else ("longitude" if "longitude" in ds.variables else lonname)
    lat2  = np.asarray(ds[lat_c].values, dtype=float)
    lon2  = np.asarray(ds[lon_c].values, dtype=float)
    if lat2.ndim == 1:
        lon2g, lat2g = np.meshgrid(lon2, lat2)
    else:
        lat2g, lon2g = lat2, lon2

    regrid_info[model] = (lat2g, lon2g, tvals, avals, valid)
    ds.close()
    ads.close()

# ----------------------------------------------------------------------
# 4) Regrid each model climatology onto the WOA23 1-degree grid
# ----------------------------------------------------------------------
def regrid_to_target(lat2, lon2, data, tlat, tlon):
    lon2 = np.mod(lon2, 360.0)
    m = np.isfinite(data) & np.isfinite(lat2) & np.isfinite(lon2)
    x, y, z = lon2[m], lat2[m], data[m]
    # longitude halo so the Delaunay triangulation is periodic
    xs = np.concatenate([x - 360.0, x, x + 360.0])
    ys = np.concatenate([y, y, y])
    zs = np.concatenate([z, z, z])
    TLON, TLAT = np.meshgrid(np.mod(tlon, 360.0), tlat)
    out = griddata((xs, ys), zs, (TLON, TLAT), method="linear")
    need = ~np.isfinite(out)
    if need.any():
        out[need] = griddata((xs, ys), zs, (TLON[need], TLAT[need]), method="nearest")
    return out

for model in MODELS:
    lat2g, lon2g, tvals, avals, valid = regrid_info[model]
    regridded = regrid_to_target(lat2g, lon2g, tvals, woa_lat, woa_lon)
    # a target cell only counts if the nearest native ocean point is reasonably close
    dist, _ = __import__("scipy").spatial.cKDTree(
        np.column_stack([np.mod(lon2g[valid], 360.0), lat2g[valid]])
    ).query(np.column_stack([np.mod(np.meshgrid(woa_lon, woa_lat)[0], 360.0).ravel(),
                             np.meshgrid(woa_lon, woa_lat)[1].ravel()]), k=1)
    dist = dist.reshape(woa_lat.size, woa_lon.size)
    regridded[dist > 1.5] = np.nan          # no native ocean point nearby -> leave undefined
    tos_on_woa[model] = regridded

# ----------------------------------------------------------------------
# 5) Common comparison cells and statistics
# ----------------------------------------------------------------------
common = woa_ocean.copy()
for model in MODELS:
    common &= np.isfinite(tos_on_woa[model])

stats = {}
for model in MODELS:
    x = tos_on_woa[model][common]
    y = woa_arr[common]
    d = x - y
    stats[model] = dict(
        model_mean=float(np.mean(x)),
        woa_mean=float(np.mean(y)),
        bias=float(np.mean(d)),
        rmse=float(np.sqrt(np.mean(d ** 2))),
        corr=float(np.corrcoef(x, y)[0, 1]),
    )

woa_mean_on_common = float(np.mean(woa_arr[common]))
n_common = int(common.sum())

ranking = sorted(MODELS, key=lambda m: abs(stats[m]["bias"]))

# ----------------------------------------------------------------------
# 6) Report
# ----------------------------------------------------------------------
print("WOA23 grid:", woa_arr.shape, " ocean cells:", int(woa_ocean.sum()),
      " common comparison cells:", n_common)
print()
print("Native-grid global mean SST 1971-2000 (degC):")
for m in MODELS:
    print(f"  {m:15s} {native_global_mean[m]:7.2f}")
print(f"  {'WOA23':15s} {woa_global_mean:7.2f}   (whole 1-deg grid; area-weighted {woa_global_mean_areawt:.2f})")
print()
print("On the common WOA 1-deg comparison cells:")
print(f"  WOA23 mean SST         {woa_mean_on_common:7.2f}")
for m in MODELS:
    s = stats[m]
    print(f"  {m:15s} mean {s['model_mean']:6.2f}  bias {s['bias']:+6.2f}  "
          f"RMSE {s['rmse']:6.3f}  r {s['corr']:6.3f}")
print()
print("Ranking (smallest -> largest |bias|):", " < ".join(ranking))
