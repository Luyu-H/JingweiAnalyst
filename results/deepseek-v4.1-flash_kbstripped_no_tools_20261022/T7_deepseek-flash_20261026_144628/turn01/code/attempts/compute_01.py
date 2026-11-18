# ============================================================================
# Compute: WOA23 vs CMIP6 nitrate/phosphate stoichiometry, 200-1000 m
# ============================================================================

import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Constants and paths
# ---------------------------------------------------------------------------
RHO = 1025.0
FILL_LIMIT = 1.0e19

woa_root = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
woa_mask_root = "/data1/lorraine-20419/ocean_data/WOA23/masks"
cmip6_root = "/data1/lorraine-20419/benchmark_cmip6/historical"

source_order = ["WOA23", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

models_info = {
    "CESM2": {
        "y_dim": "nlat",
        "x_dim": "nlon",
        "lat_coord": "lat",
        "lon_coord": "lon",
        "depth_coord": "lev",
        "depth_units": "centimeters",
        "depth_scale": 0.01,
    },
    "GFDL-ESM4": {
        "y_dim": "lat",
        "x_dim": "lon",
        "lat_coord": "lat",
        "lon_coord": "lon",
        "depth_coord": "lev",
        "depth_units": "m",
        "depth_scale": 1.0,
    },
    "MPI-ESM1-2-LR": {
        "y_dim": "j",
        "x_dim": "i",
        "lat_coord": "latitude",
        "lon_coord": "longitude",
        "depth_coord": "lev",
        "depth_units": "m",
        "depth_scale": 1.0,
    },
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def find_coord(ds, candidates):
    for c in candidates:
        if c in ds.coords or c in ds.variables:
            return c
    raise KeyError(f"None of {candidates} found in dataset")


def to_m(values, units):
    """Convert a depth array to metres according to its declared units."""
    u = str(units).lower().replace(" ", "")
    if "centimeter" in u or "centimetre" in u or u == "cm" or u.endswith("cm"):
        return values * 0.01
    if "millimeter" in u or "millimetre" in u or u == "mm" or u.endswith("mm"):
        return values * 0.001
    if "kilometer" in u or "kilometre" in u or u == "km" or u.endswith("km"):
        return values * 1000.0
    return values


def get_depth_info(ds, lev_name):
    """Return level centres and layer bounds in metres."""
    lev = ds[lev_name]
    lev_m = to_m(np.asarray(lev.values, dtype=float), lev.attrs.get("units", ""))

    bvals = None
    bname = lev_name + "_bnds"
    if bname in ds:
        b = ds[bname]
        if "time" in b.dims:
            b = b.isel(time=0)
        bvals = to_m(np.asarray(b.values, dtype=float), b.attrs.get("units", ""))

    if bvals is None:
        n = len(lev_m)
        edges = np.empty((n, 2), dtype=float)
        if n > 1:
            mids = 0.5 * (lev_m[:-1] + lev_m[1:])
            edges[1:, 0] = mids
            edges[:-1, 1] = mids
            edges[0, 0] = lev_m[0] - 0.5 * (lev_m[1] - lev_m[0])
            edges[-1, 1] = lev_m[-1] + 0.5 * (lev_m[-1] - lev_m[-2])
        else:
            edges[0] = [lev_m[0] - 1.0, lev_m[0] + 1.0]
        bvals = edges

    if bvals.ndim != 2 or bvals.shape[0] != len(lev_m):
        if bvals.ndim == 2 and bvals.shape[1] == len(lev_m):
            bvals = bvals.T
        else:
            raise ValueError("Depth bounds have incompatible shape")

    lo = np.minimum(bvals[:, 0], bvals[:, 1])
    hi = np.maximum(bvals[:, 0], bvals[:, 1])
    return lev_m, lo, hi


def umol_kg_factor(units, rho=RHO):
    """Factor converting the declared concentration units to umol kg-1."""
    u = str(units).lower().replace(" ", "")
    if "umol" in u or "µmol" in u:
        if "kg" in u:
            return 1.0
        if "m-3" in u or "m^-3" in u or "m**-3" in u:
            return 1.0 / rho
    if "mmol" in u:
        return 1000.0 / rho
    if "mol" in u:
        return 1.0e6 / rho
    return 1.0e6 / rho


def get_lonlat_2d(ds, latname, lonname):
    lat = ds[latname]
    lon = ds[lonname]
    if lat.ndim == 1 and lon.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon.values, lat.values)
    else:
        lat2d = np.asarray(lat.values)
        lon2d = np.asarray(lon.values)
        if lat2d.ndim == 1 or lon2d.ndim == 1:
            lon2d, lat2d = np.meshgrid(np.ravel(lon.values), np.ravel(lat.values))
    return lat2d, lon2d


def classify_basin_centroid(lon, lat):
    """Fold marginal seas into their parent ocean using a code centroid."""
    # Excluded basins
    if lat < -60.0:
        return "Southern"
    if lat > 66.0:
        return "Arctic"

    # Atlantic marginal seas: Mediterranean, Black Sea, Baltic, North Sea
    if (-10.0 <= lon <= 50.0) and (30.0 <= lat <= 72.0):
        return "Atlantic"
    # Hudson Bay
    if (-105.0 <= lon <= -65.0) and (45.0 <= lat <= 70.0):
        return "Atlantic"

    # Indian marginal seas: Red Sea, Persian Gulf
    if (30.0 <= lon <= 65.0) and (10.0 <= lat <= 35.0):
        return "Indian"

    # Pacific marginal seas: South China Sea, Japan Sea, Okhotsk, Bering
    if (100.0 <= lon <= 180.0) and (0.0 <= lat <= 65.0):
        return "Pacific"
    if (-180.0 <= lon <= -150.0) and (45.0 <= lat <= 65.0):
        return "Pacific"

    # Major basins by longitude
    if -70.0 <= lon <= 20.0:
        return "Atlantic"
    if 20.0 < lon <= 100.0:
        return "Indian"
    return "Pacific"


def read_basin_mask(path):
    rows = []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 3:
                continue
            vals = []
            ok = True
            for p in parts:
                p2 = p.replace("*", "nan").replace("D", "E").replace("d", "e")
                try:
                    vals.append(float(p2))
                except ValueError:
                    ok = False
                    break
            if ok and len(vals) >= 3:
                rows.append(vals)
    if not rows:
        raise ValueError("No numeric rows in basin mask")
    maxlen = max(len(r) for r in rows)
    arr = np.full((len(rows), maxlen), np.nan, dtype=float)
    for i, r in enumerate(rows):
        arr[i, : len(r)] = r
    return arr


def detect_latlon_cols(arr):
    ncol = arr.shape[1]
    candidates = list(range(min(4, ncol)))
    for i in candidates:
        if not np.any(np.isfinite(arr[:, i])):
            continue
        if np.nanmin(arr[:, i]) >= -90.1 and np.nanmax(arr[:, i]) <= 90.1:
            for j in candidates:
                if j == i:
                    continue
                if not np.any(np.isfinite(arr[:, j])):
                    continue
                if np.nanmin(arr[:, j]) >= -180.1 and np.nanmax(arr[:, j]) <= 360.1:
                    if np.nanmax(np.abs(arr[:, i] - arr[:, j])) > 1.0:
                        return i, j
    return 0, 1


def get_value_cols(arr, lat_col, lon_col):
    cols = [c for c in range(arr.shape[1]) if c not in (lat_col, lon_col)]
    if lat_col > 0 and 0 in cols:
        col0 = arr[:, 0]
        if np.all(np.isfinite(col0)) and np.allclose(col0, np.round(col0)):
            if np.all(np.diff(col0) >= 0):
                cols.remove(0)
    return cols


def sample_basin_grid(lat2d, lon2d, basin_grid, lat_centers, lon_centers):
    lat = np.asarray(lat2d)
    lon = np.asarray(lon2d)
    lon = ((lon + 180.0) % 360.0) - 180.0
    dlat = np.median(np.diff(lat_centers))
    dlon = np.median(np.diff(lon_centers))
    i = np.rint((lat - lat_centers[0]) / dlat).astype(int)
    j = np.rint((lon - lon_centers[0]) / dlon).astype(int)
    i = np.clip(i, 0, basin_grid.shape[0] - 1)
    j = np.clip(j, 0, basin_grid.shape[1] - 1)
    return basin_grid[i, j]


def vol_exceedance(np_ratio, vol, cond, extra_mask=None):
    m = np.isfinite(np_ratio) & np.isfinite(vol) & (vol > 0.0)
    if extra_mask is not None:
        m = m & extra_mask
    total = np.nansum(vol[m])
    num = np.nansum(vol[m & cond])
    pct = 100.0 * num / total if total > 0.0 else np.nan
    return pct, total, num


def compute_stats(np_ratio, nitrate, phosphate, vol, basin=None):
    valid = (
        np.isfinite(np_ratio)
        & np.isfinite(nitrate)
        & np.isfinite(phosphate)
        & np.isfinite(vol)
        & (vol > 0.0)
    )
    total_vol = np.nansum(vol[valid])
    if total_vol <= 0.0:
        return {
            "mean_np": np.nan,
            "ratio_means": np.nan,
            "mean_nitrate": np.nan,
            "mean_phosphate": np.nan,
            "ols_slope": np.nan,
            "ols_intercept": np.nan,
            "r2": np.nan,
            "pct_14_18": np.nan,
            "pct_lt12": np.nan,
            "basin_means": {"Atlantic": np.nan, "Pacific": np.nan, "Indian": np.nan},
            "basin_exceedance": {},
            "valid_count": 0,
            "total_volume": 0.0,
        }

    mean_np = np.nansum(np_ratio[valid] * vol[valid]) / total_vol
    mean_n = np.nansum(nitrate[valid] * vol[valid]) / total_vol
    mean_p = np.nansum(phosphate[valid] * vol[valid]) / total_vol
    ratio_means = mean_n / mean_p if mean_p != 0.0 else np.nan

    x = phosphate[valid]
    y = nitrate[valid]
    a, b = np.polyfit(x, y, 1)
    ypred = a * x + b
    ss_res = np.sum((y - ypred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else np.nan

    pct_14_18, vol_14_18, num_14_18 = vol_exceedance(
        np_ratio, vol, (np_ratio > 14.0) & (np_ratio < 18.0)
    )
    pct_lt12, vol_lt12, num_lt12 = vol_exceedance(np_ratio, vol, np_ratio < 12.0)

    basin_means = {}
    for bname, bval in [("Atlantic", 1), ("Pacific", 2), ("Indian", 3)]:
        if basin is not None:
            bm = basin == bval
        else:
            bm = np.ones_like(vol, dtype=bool)
        m = valid & bm
        if np.any(m):
            basin_means[bname] = np.nansum(np_ratio[m] * vol[m]) / np.nansum(vol[m])
        else:
            basin_means[bname] = np.nan

    basin_exceedance = {}
    for bname, bval in [("Pacific", 2), ("Indian", 3)]:
        if basin is not None:
            bm = basin == bval
        else:
            bm = np.ones_like(vol, dtype=bool)
        pct_lt, vtot_lt, vnum_lt = vol_exceedance(
            np_ratio, vol, np_ratio < 12.0, extra_mask=bm
        )
        pct_mid, vtot_mid, vnum_mid = vol_exceedance(
            np_ratio, vol, (np_ratio > 14.0) & (np_ratio < 18.0), extra_mask=bm
        )
        basin_exceedance[bname] = {
            "lt12_pct": pct_lt,
            "lt12_volume": vtot_lt,
            "lt12_volume_below": vnum_lt,
            "14_18_pct": pct_mid,
            "14_18_volume": vtot_mid,
            "14_18_volume_between": vnum_mid,
        }

    return {
        "mean_np": mean_np,
        "ratio_means": ratio_means,
        "mean_nitrate": mean_n,
        "mean_phosphate": mean_p,
        "ols_slope": a,
        "ols_intercept": b,
        "r2": r2,
        "pct_14_18": pct_14_18,
        "vol_14_18": vol_14_18,
        "vol_between_14_18": num_14_18,
        "pct_lt12": pct_lt12,
        "vol_lt12": vol_lt12,
        "vol_below_12": num_lt12,
        "basin_means": basin_means,
        "basin_exceedance": basin_exceedance,
        "valid_count": int(np.sum(valid)),
        "total_volume": total_vol,
    }


# ---------------------------------------------------------------------------
# WOA23 nitrate and phosphate annual objectively-analysed fields
# ---------------------------------------------------------------------------
def open_woa(var):
    if var == "nitrate":
        subdir = "nitrate/all_time_periods"
        tag = "all"
        prefix = "n"
    elif var == "phosphate":
        subdir = "phosphate/all_time_periods"
        tag = "all"
        prefix = "p"
    else:
        raise ValueError(var)

    path = os.path.join(
        woa_root, subdir, f"woa23_{tag}_{prefix}00_01.nc"
    )
    ds = xr.open_dataset(path, decode_times=False)

    candidates = [f"{prefix}_an", var, prefix, f"{prefix}_mn"]
    vname = None
    for c in candidates:
        if c in ds.data_vars:
            vname = c
            break
    if vname is None:
        vname = list(ds.data_vars)[0]

    da = ds[vname]
    if "time" in da.dims:
        da = da.isel(time=0)
    da = da.where(np.abs(da) < FILL_LIMIT)
    return da, ds, vname


no3_da, no3_ds, no3_vname = open_woa("nitrate")
po4_da, po4_ds, po4_vname = open_woa("phosphate")

no3_da, po4_da = xr.align(no3_da, po4_da, join="inner")

zname = find_coord(no3_ds, ["depth", "lev", "deptht"])
latname = find_coord(no3_ds, ["lat", "latitude"])
lonname = find_coord(no3_ds, ["lon", "longitude"])

no3_da = no3_da.transpose(zname, latname, lonname)
po4_da = po4_da.transpose(zname, latname, lonname)

lat_woa = np.asarray(no3_ds[latname].values, dtype=float)
lon_woa = np.asarray(no3_ds[lonname].values, dtype=float)
lon2d_woa, lat2d_woa = np.meshgrid(lon_woa, lat_woa)

z_m, zlo, zhi = get_depth_info(no3_ds, zname)
dz = zhi - zlo

lev_in = (z_m >= 200.0) & (z_m <= 1000.0)
if not np.any(lev_in):
    lev_in = (zhi > 200.0) & (zlo < 1000.0)

lev_idx_woa = np.where(lev_in)[0]
dz_sel = dz[lev_idx_woa]

no3_sel = np.asarray(no3_da.values, dtype=float)[lev_idx_woa, :, :]
po4_sel = np.asarray(po4_da.values, dtype=float)[lev_idx_woa, :, :]

valid_n = np.isfinite(no3_sel)
valid_p = np.isfinite(po4_sel)

w_n = dz_sel[:, None, None] * valid_n
num_n = np.nansum(no3_sel * w_n, axis=0)
den_n = np.nansum(w_n, axis=0)
no3_layer_woa = np.where(den_n > 0.0, num_n / den_n, np.nan)

w_p = dz_sel[:, None, None] * valid_p
num_p = np.nansum(po4_sel * w_p, axis=0)
den_p = np.nansum(w_p, axis=0)
po4_layer_woa = np.where(den_p > 0.0, num_p / den_p, np.nan)

valid_both_woa = valid_n & valid_p
dz_both = np.nansum(dz_sel[:, None, None] * valid_both_woa, axis=0)

coslat = np.cos(np.deg2rad(lat_woa))
dlon_rad = np.deg2rad(np.abs(np.median(np.diff(lon_woa)))) if len(lon_woa) > 1 else np.deg2rad(1.0)
vol_woa = coslat[:, None] * dlon_rad * dz_both

with np.errstate(divide="ignore", invalid="ignore"):
    np_woa = no3_layer_woa / po4_layer_woa
np_woa = np.where(np.isfinite(np_woa), np_woa, np.nan)

# ---------------------------------------------------------------------------
# WOA23 basin mask and parent-ocean labels
# ---------------------------------------------------------------------------
basin_path = os.path.join(woa_mask_root, "basinmask_01.msk")
basin_arr = read_basin_mask(basin_path)

lat_col, lon_col = detect_latlon_cols(basin_arr)
value_cols = get_value_cols(basin_arr, lat_col, lon_col)

basin_codes = np.full((len(lat_woa), len(lon_woa)), np.nan, dtype=float)
dlat = np.median(np.diff(lat_woa)) if len(lat_woa) > 1 else 1.0
dlon = np.median(np.diff(lon_woa)) if len(lon_woa) > 1 else 1.0

for row in basin_arr:
    lat_val = row[lat_col]
    lon_val = row[lon_col]
    if not np.isfinite(lat_val) or not np.isfinite(lon_val):
        continue
    lon_val = ((lon_val + 180.0) % 360.0) - 180.0
    vals = row[value_cols]
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        continue
    code = vals[0]  # surface label; depth_dependent=False
    i = int(np.rint((lat_val - lat_woa[0]) / dlat))
    j = int(np.rint((lon_val - lon_woa[0]) / dlon))
    if 0 <= i < len(lat_woa) and 0 <= j < len(lon_woa):
        basin_codes[i, j] = code

codes = np.unique(basin_codes[np.isfinite(basin_codes)])
code_parent = {}
for code in codes:
    m = basin_codes == code
    lon_rad = np.deg2rad(lon2d_woa[m])
    lon_c = np.rad2deg(
        np.arctan2(np.nanmean(np.sin(lon_rad)), np.nanmean(np.cos(lon_rad)))
    )
    lat_c = np.nanmean(lat2d_woa[m])
    code_parent[code] = classify_basin_centroid(lon_c, lat_c)

name_to_int = {"Atlantic": 1, "Pacific": 2, "Indian": 3, "Southern": 4, "Arctic": 5}
basin_parent_grid = np.full_like(basin_codes, np.nan, dtype=float)
for code, parent in code_parent.items():
    basin_parent_grid[basin_codes == code] = name_to_int[parent]

# ---------------------------------------------------------------------------
# CMIP6 loaders
# ---------------------------------------------------------------------------
def load_cmip6_clim(model, var):
    files = sorted(glob.glob(os.path.join(cmip6_root, model, var, "*.nc")))
    if len(files) == 0:
        return None, None

    ds = xr.open_mfdataset(
        files,
        combine="by_coords",
        use_cftime=True,
        decode_times=True,
    )

    vname = var if var in ds.data_vars else None
    if vname is None:
        for v in ds.data_vars:
            if v.lower() == var.lower():
                vname = v
                break
    if vname is None:
        vname = list(ds.data_vars)[0]

    da = ds[vname]
    units = da.attrs.get("units", "")

    time = ds["time"]
    mask_t = ((time.dt.year >= 1971) & (time.dt.year <= 2000)).values
    da = da.isel(time=mask_t)
    da = da.where(np.abs(da) < FILL_LIMIT)

    clim = da.groupby("time.month").mean("time")
    annual = clim.mean("month")

    factor = umol_kg_factor(units, RHO)
    annual = annual * factor
    annual.attrs["units"] = "micromoles_per_kilogram"
    return annual, ds


def load_cmip6_fixed(model, var):
    files = sorted(glob.glob(os.path.join(cmip6_root, model, var, "*.nc")))
    if len(files) == 0:
        return None, None

    ds = xr.open_mfdataset(
        files,
        combine="by_coords",
        use_cftime=True,
        decode_times=True,
    )

    vname = var if var in ds.data_vars else list(ds.data_vars)[0]
    da = ds[vname]

    if "time" in da.dims:
        time = ds["time"]
        mt = ((time.dt.year >= 1971) & (time.dt.year <= 2000)).values
        if np.any(mt):
            da = da.isel(time=mt).mean("time")
        else:
            da = da.isel(time=0)

    da = da.where(np.abs(da) < FILL_LIMIT)
    return da, ds


# ---------------------------------------------------------------------------
# Aggregate results
# ---------------------------------------------------------------------------
sources = {}
results = {}

stats_woa = compute_stats(
    np_woa, no3_layer_woa, po4_layer_woa, vol_woa, basin_parent_grid
)
sources["WOA23"] = {
    "np_ratio": np_woa,
    "nitrate": no3_layer_woa,
    "phosphate": po4_layer_woa,
    "vol": vol_woa,
    "lat": lat2d_woa,
    "lon": lon2d_woa,
    "basin": basin_parent_grid,
    "stats": stats_woa,
}
results["WOA23"] = stats_woa


for model in ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]:
    info = models_info[model]

    no3_ann, ds_no3 = load_cmip6_clim(model, "no3")
    po4_ann, ds_po4 = load_cmip6_clim(model, "po4")

    if no3_ann is None or po4_ann is None:
        results[model] = {"missing": "no3 or po4"}
        sources[model] = None
        continue

    vol_da, ds_vol = load_cmip6_fixed(model, "volcello")
    deptho_da, ds_deptho = load_cmip6_fixed(model, "deptho")

    if vol_da is None or deptho_da is None:
        results[model] = {"missing": "volcello or deptho"}
        sources[model] = None
        continue

    no3_ann, po4_ann, vol_da = xr.align(
        no3_ann, po4_ann, vol_da, join="exact"
    )

    lev_name = info["depth_coord"]
    lev_m, zlo_m, zhi_m = get_depth_info(ds_no3, lev_name)

    lev_in = (lev_m >= 200.0) & (lev_m <= 1000.0)
    if not np.any(lev_in):
        lev_in = (zhi_m > 200.0) & (zlo_m < 1000.0)

    lev_idx = np.where(lev_in)[0]

    no3_sel = no3_ann.isel({lev_name: lev_idx})
    po4_sel = po4_ann.isel({lev_name: lev_idx})
    vol_sel = vol_da.isel({lev_name: lev_idx})

    w_n = vol_sel.where(np.isfinite(no3_sel))
    no3_layer = (no3_sel * w_n).sum(lev_name) / w_n.sum(lev_name)

    w_p = vol_sel.where(np.isfinite(po4_sel))
    po4_layer = (po4_sel * w_p).sum(lev_name) / w_p.sum(lev_name)

    valid_both = np.isfinite(no3_sel) & np.isfinite(po4_sel)
    vol_both = vol_sel.where(valid_both).sum(lev_name)

    np_ratio = no3_layer / po4_layer

    lat2d, lon2d = get_lonlat_2d(ds_no3, info["lat_coord"], info["lon_coord"])

    deptho = deptho_da.transpose(info["y_dim"], info["x_dim"]).values.astype(float)
    ocean = np.isfinite(deptho) & (deptho > 0.0)

    np_ratio_np = np.where(ocean, np_ratio.values.astype(float), np.nan)
    no3_np = np.where(ocean, no3_layer.values.astype(float), np.nan)
    po4_np = np.where(ocean, po4_layer.values.astype(float), np.nan)
    vol_np = np.where(ocean, vol_both.values.astype(float), np.nan)

    basin_model = sample_basin_grid(
        lat2d, lon2d, basin_parent_grid, lat_woa, lon_woa
    )

    stats_model = compute_stats(
        np_ratio_np, no3_np, po4_np, vol_np, basin_model
    )

    sources[model] = {
        "np_ratio": np_ratio_np,
        "nitrate": no3_np,
        "phosphate": po4_np,
        "vol": vol_np,
        "lat": lat2d,
        "lon": lon2d,
        "basin": basin_model,
        "stats": stats_model,
    }
    results[model] = stats_model


# ---------------------------------------------------------------------------
# Print summary table
# ---------------------------------------------------------------------------
summary_rows = []
for src in source_order:
    st = results.get(src)
    if st is None or "missing" in st:
        summary_rows.append(
            {
                "source": src,
                "status": "missing",
                "mean_NP": np.nan,
                "ratio_of_means": np.nan,
                "slope": np.nan,
                "intercept": np.nan,
                "R2": np.nan,
                "pct_14_18": np.nan,
                "pct_lt12": np.nan,
                "NP_Atlantic": np.nan,
                "NP_Pacific": np.nan,
                "NP_Indian": np.nan,
                "pct_lt12_Pacific": np.nan,
                "pct_lt12_Indian": np.nan,
            }
        )
        continue

    summary_rows.append(
        {
            "source": src,
            "status": "ok",
            "mean_NP": st["mean_np"],
            "ratio_of_means": st["ratio_means"],
            "slope": st["ols_slope"],
            "intercept": st["ols_intercept"],
            "R2": st["r2"],
            "pct_14_18": st["pct_14_18"],
            "pct_lt12": st["pct_lt12"],
            "NP_Atlantic": st["basin_means"].get("Atlantic", np.nan),
            "NP_Pacific": st["basin_means"].get("Pacific", np.nan),
            "NP_Indian": st["basin_means"].get("Indian", np.nan),
            "pct_lt12_Pacific": st["basin_exceedance"].get("Pacific", {}).get("lt12_pct", np.nan),
            "pct_lt12_Indian": st["basin_exceedance"].get("Indian", {}).get("lt12_pct", np.nan),
        }
    )

summary_df = pd.DataFrame(summary_rows)
print(summary_df.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
