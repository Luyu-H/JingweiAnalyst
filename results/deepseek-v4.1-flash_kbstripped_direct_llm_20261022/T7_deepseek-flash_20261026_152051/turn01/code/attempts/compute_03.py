import os
import glob
import re
import warnings
from collections import Counter

import numpy as np
import pandas as pd
import xarray as xr

warnings.filterwarnings("ignore")

# -----------------------------------------------------------------------------
# Paths and constants
# -----------------------------------------------------------------------------
WOA_ROOT = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
WOA_MASKS = "/data1/lorraine-20419/ocean_data/WOA23/masks"
CMIP_ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"

R_EARTH = 6371000.0
RHO = 1025.0  # kg m-3, for CMIP6 mol m-3 -> umol kg-1 conversion

SOURCE_ORDER = ["WOA23", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
CMIP_MODELS = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

MODEL_CONFIG = {
    "CESM2": {"lat": "lat", "lon": "lon"},
    "GFDL-ESM4": {"lat": "lat", "lon": "lon"},
    "MPI-ESM1-2-LR": {"lat": "latitude", "lon": "longitude"},
}

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def get_variant(path):
    m = re.search(r"_(r\d+i\d+p\d+f\d+)_", os.path.basename(path))
    return m.group(1) if m else None


def filter_files_by_variant(files, variant):
    if variant is None:
        return files
    out = [f for f in files if get_variant(f) == variant]
    return out if out else files


def find_var_files(model, var):
    d = os.path.join(CMIP_ROOT, model, var)
    patterns = [
        os.path.join(d, f"{var}_*_{model}_historical_*.nc"),
        os.path.join(d, f"*{var}*.nc"),
    ]
    for pat in patterns:
        files = sorted(glob.glob(pat))
        if files:
            return files
    return []


def open_cmip_dataset(files):
    if not files:
        raise FileNotFoundError("No CMIP6 files found")
    ds0 = xr.open_dataset(files[0], use_cftime=True)
    has_time = "time" in ds0.dims
    ds0.close()

    kwargs = {"use_cftime": True}
    if has_time:
        kwargs["chunks"] = {"time": 60}

    if len(files) == 1:
        return xr.open_dataset(files[0], **kwargs)
    return xr.open_mfdataset(files, combine="by_coords", **kwargs)


def period_mean(da, y0=1971, y1=2000):
    if "time" in da.dims:
        if da.sizes.get("time", 0) > 1:
            yr = da["time"].dt.year
            da = da.where((yr >= y0) & (yr <= y1), drop=True)
            da = da.mean("time")
        else:
            # Static / singleton time field (e.g. an Ofx volcello): keep it.
            da = da.isel(time=0)
    return da


def get_latlon2d(ds, model):
    cfg = MODEL_CONFIG[model]
    lat_da = ds[cfg["lat"]]
    lon_da = ds[cfg["lon"]]
    if lat_da.ndim == 1 and lon_da.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon_da.values, lat_da.values)
    else:
        lat2d = lat_da.values
        lon2d = lon_da.values
    return lat2d, lon2d


def bounds_overlap_from_centers(z, zmin, zmax):
    z = np.asarray(z, dtype=float)
    if z.size == 0:
        return np.array([])
    if z.size == 1:
        b0 = z - 0.5
        b1 = z + 0.5
    else:
        mid = (z[1:] + z[:-1]) / 2.0
        b0 = np.empty_like(z)
        b1 = np.empty_like(z)
        b0[0] = z[0] - (mid[0] - z[0])
        b0[1:] = mid
        b1[:-1] = mid
        b1[-1] = z[-1] + (z[-1] - mid[-1])
    overlap = np.clip(np.minimum(b1, zmax) - np.maximum(b0, zmin), 0.0, None)
    return overlap


def get_depth_frac(ds, depth_dim="lev", zmin=200.0, zmax=1000.0):
    """Fraction of each model level that lies within [zmin, zmax].

    The vertical coordinate is forced to be one-dimensional: some CMIP6 files
    store ``lev`` with a trailing singleton (or as an auxiliary 2-D field),
    which would otherwise make the returned fraction two-dimensional and
    break the subsequent ``xr.DataArray`` construction.
    """
    z = np.atleast_1d(np.squeeze(np.asarray(ds[depth_dim].values, dtype=float)))
    if z.ndim > 1:
        z = z.reshape(z.shape[0], -1)[:, 0]
    z_units = str(ds[depth_dim].attrs.get("units", "m")).lower()
    if ("cm" in z_units) or ("centimeter" in z_units):
        z = z * 0.01

    bname = None
    for cand in [f"{depth_dim}_bnds", f"{depth_dim}_bounds", "lev_bnds", "depth_bnds"]:
        if cand in ds:
            bname = cand
            break

    b0 = b1 = None
    if bname is not None:
        b = np.squeeze(np.asarray(ds[bname].values, dtype=float))
        b_units = str(ds[bname].attrs.get("units", "m")).lower()
        if ("cm" in b_units) or ("centimeter" in b_units):
            b = b * 0.01
        if b.ndim == 2:
            if b.shape[0] != z.size and b.shape[-1] == z.size:
                b = b.T
            if b.shape[0] == z.size and b.shape[1] >= 2:
                b0 = np.minimum(b[:, 0], b[:, 1])
                b1 = np.maximum(b[:, 0], b[:, 1])

    if b0 is None or b1 is None:
        if z.size == 1:
            b0 = z - 0.5
            b1 = z + 0.5
        else:
            mid = (z[1:] + z[:-1]) / 2.0
            b0 = np.empty_like(z)
            b1 = np.empty_like(z)
            b0[0] = z[0] - (mid[0] - z[0])
            b0[1:] = mid
            b1[:-1] = mid
            b1[-1] = z[-1] + (z[-1] - mid[-1])

    thick = b1 - b0
    overlap = np.clip(np.minimum(b1, zmax) - np.maximum(b0, zmin), 0.0, None)
    frac = np.where(thick > 0, overlap / thick, 0.0)
    return np.asarray(frac, dtype=float).reshape(-1)


def open_woa_var(var, prefix):
    subdir = f"{var}/all_time_periods"
    path = os.path.join(WOA_ROOT, subdir, f"woa23_all_{prefix}00_01.nc")
    ds = xr.open_dataset(path, decode_times=False)
    vname = f"{prefix}_an"
    if vname not in ds.data_vars:
        cands = [v for v in ds.data_vars if v.endswith("_an")]
        if not cands:
            cands = [v for v in ds.data_vars if "an" in v.lower()]
        vname = cands[0]
    da = ds[vname]
    if "time" in da.dims:
        da = da.isel(time=0)
    return da, ds


# -----------------------------------------------------------------------------
# Basin mask reader
# -----------------------------------------------------------------------------
def _read_numeric_rows(path):
    """
    Read a sparse text table, skipping comment/header lines.  Each retained
    row is the list of numeric tokens found on that line (leading/trailing
    non-numeric tokens such as a basin name are tolerated).
    """
    rows = []
    with open(path, "r") as fh:
        for line in fh:
            s = line.strip()
            if not s or s[0] in "#!%":
                continue
            # strip inline comments
            for mark in ("#", "!", "%", ";"):
                p = s.find(mark)
                if p != -1:
                    s = s[:p]
            s = s.strip()
            if not s:
                continue
            toks = s.replace(",", " ").split()
            vals = []
            for t in toks:
                try:
                    vals.append(float(t))
                except ValueError:
                    continue
            if len(vals) >= 3:
                rows.append(vals)
    return rows


def _rasterize_basin2d(lat, lon, basin, woa_lat, woa_lon):
    nlat, nlon = len(woa_lat), len(woa_lon)
    dlat = woa_lat[1] - woa_lat[0]
    dlon = woa_lon[1] - woa_lon[0]

    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float).copy()
    basin = np.asarray(basin, dtype=float)

    # Normalise longitude to the WOA longitude range.
    lon = np.where(lon < woa_lon.min(), lon + 360.0, lon)
    lon = np.where(lon > woa_lon.max(), lon - 360.0, lon)

    i = np.rint((lat - woa_lat[0]) / dlat).astype(int)
    j = np.rint((lon - woa_lon[0]) / dlon).astype(int)

    basin2d = np.full((nlat, nlon), np.nan, dtype=float)
    ok = ((i >= 0) & (i < nlat) & (j >= 0) & (j < nlon) & np.isfinite(basin))
    basin2d[i[ok], j[ok]] = basin[ok]
    return basin2d


def load_woa_basin_mask(mask_path, woa_depth, woa_lat, woa_lon):
    """
    Parse the WOA23 basin mask.  The WOA mask is a sparse table with one row
    per ocean grid point; depending on the exact product version rows are
    either ``Latitude Longitude Basin`` (3 columns) or
    ``Latitude Longitude Level Basin`` (4 columns), or one basin column per
    level (>4 columns).  All are handled here, and the geometry is rasterised
    onto the WOA 1-degree grid.
    """
    rows = _read_numeric_rows(mask_path)
    if not rows:
        raise RuntimeError("Could not parse WOA basin mask")

    # The most common row width defines the table layout.
    width_counts = Counter(len(r) for r in rows)
    ncol = width_counts.most_common(1)[0][0]
    rows = [r for r in rows if len(r) == ncol]
    arr = np.array(rows, dtype=float)

    nlev = len(woa_depth)
    nlat = len(woa_lat)
    nlon = len(woa_lon)
    dlat = woa_lat[1] - woa_lat[0]
    dlon = woa_lon[1] - woa_lon[0]

    # ------------------------------------------------------------------
    # 3 columns: Latitude Longitude Basin (a purely horizontal mask)
    # ------------------------------------------------------------------
    if ncol == 3:
        # Longitude is the column whose values exceed +/-90.
        lon_col = None
        for c in range(3):
            if np.nanmax(np.abs(arr[:, c])) > 90:
                lon_col = c
                break
        if lon_col is None:
            spans = [np.nanmax(arr[:, c]) - np.nanmin(arr[:, c]) for c in range(3)]
            lon_col = int(np.argmax(spans))

        rem = [c for c in range(3) if c != lon_col]
        lat_col = None
        for c in rem:
            col = arr[:, c]
            if np.nanmin(col) >= -90 and np.nanmax(col) <= 90:
                if not np.all(np.abs(col - np.round(col)) < 1e-6):
                    lat_col = c
                    break
        if lat_col is None:
            for c in rem:
                col = arr[:, c]
                if np.nanmin(col) >= -90 and np.nanmax(col) <= 90:
                    lat_col = c
                    break
        if lat_col is None:
            lat_col = rem[0]
        basin_col = [c for c in rem if c != lat_col][0]

        basin2d = _rasterize_basin2d(arr[:, lat_col], arr[:, lon_col],
                                     arr[:, basin_col], woa_lat, woa_lon)
        basin3d = np.repeat(basin2d[None, :, :], nlev, axis=0)
        return basin2d, basin3d

    # ------------------------------------------------------------------
    # More than 4 columns: Latitude Longitude Basin_at_level_1 ... _N
    # ------------------------------------------------------------------
    if ncol > 4:
        lon_col = None
        for c in range(ncol):
            if np.nanmax(np.abs(arr[:, c])) > 90:
                lon_col = c
                break
        if lon_col is None:
            lon_col = 1
        rem = [c for c in range(ncol) if c != lon_col]
        lat_col = None
        for c in rem:
            col = arr[:, c]
            if np.nanmin(col) >= -90 and np.nanmax(col) <= 90:
                lat_col = c
                break
        if lat_col is None:
            lat_col = rem[0]
        bcols = [c for c in range(ncol) if c not in (lat_col, lon_col)]
        # take the level column nearest 500 m (levels ordered surface -> deep)
        kmid = int(np.abs(np.asarray(woa_depth) - 500.0).argmin())
        bsel = bcols[min(kmid, len(bcols) - 1)]
        basin = arr[:, bsel]
        basin = np.where(np.isfinite(basin) & (basin >= 0) & (basin <= 30),
                         basin, np.nan)
        basin2d = _rasterize_basin2d(arr[:, lat_col], arr[:, lon_col],
                                     basin, woa_lat, woa_lon)
        basin3d = np.repeat(basin2d[None, :, :], nlev, axis=0)
        return basin2d, basin3d

    # ------------------------------------------------------------------
    # Exactly 4 columns: Latitude Longitude Level Basin
    # ------------------------------------------------------------------
    lon_col = None
    for c in range(4):
        if np.nanmax(np.abs(arr[:, c])) > 90:
            lon_col = c
            break
    if lon_col is None:
        lon_col = 1

    rem = [c for c in range(4) if c != lon_col]
    lat_col = None
    for c in rem:
        col = arr[:, c]
        if np.nanmin(col) >= -90 and np.nanmax(col) <= 90:
            if not np.all(np.abs(col - np.round(col)) < 1e-6):
                lat_col = c
                break
    if lat_col is None:
        for c in rem:
            col = arr[:, c]
            if np.nanmin(col) >= -90 and np.nanmax(col) <= 90:
                lat_col = c
                break
    if lat_col is None:
        lat_col = rem[0]

    rem2 = [c for c in range(4) if c not in (lon_col, lat_col)]
    basin_col = None
    for c in rem2:
        col = arr[:, c]
        if (np.nanmax(col) <= 20 and np.nanmin(col) >= 0
                and np.all(np.abs(col - np.round(col)) < 1e-6)):
            basin_col = c
            break
    if basin_col is None and rem2:
        basin_col = rem2[0]
    level_col = [c for c in rem2 if c != basin_col][0] if len(rem2) > 1 else None

    lat = arr[:, lat_col]
    lon = arr[:, lon_col]
    basin = arr[:, basin_col]
    level = arr[:, level_col].astype(float) if level_col is not None else None

    lon = np.where(lon < woa_lon.min(), lon + 360.0, lon)
    lon = np.where(lon > woa_lon.max(), lon - 360.0, lon)

    i = np.rint((lat - woa_lat[0]) / dlat).astype(int)
    j = np.rint((lon - woa_lon[0]) / dlon).astype(int)

    if level is not None:
        if np.nanmax(level) <= (nlev + 1) and np.nanmin(level) >= 0:
            if np.nanmin(level) >= 1:
                k = np.rint(level).astype(int) - 1
            else:
                k = np.rint(level).astype(int)
        else:
            k = np.abs(level[:, None] - np.asarray(woa_depth)[None, :]).argmin(axis=1)
    else:
        k = np.zeros(len(lat), dtype=int)

    basin3d = np.full((nlev, nlat, nlon), np.nan, dtype=float)
    ok = ((i >= 0) & (i < nlat) & (j >= 0) & (j < nlon)
          & (k >= 0) & (k < nlev) & np.isfinite(basin))
    basin3d[k[ok], i[ok], j[ok]] = basin[ok]

    # Collapse to a 2-D basin map by taking the modal basin at each cell.
    basin2d = np.full((nlat, nlon), np.nan, dtype=float)
    for ii in range(nlat):
        for jj in range(nlon):
            vals = basin3d[:, ii, jj]
            vals = vals[np.isfinite(vals)]
            if vals.size:
                u, c = np.unique(vals.astype(int), return_counts=True)
                basin2d[ii, jj] = u[c.argmax()]

    return basin2d, basin3d


def build_parent_codes(basin2d, lat2d, lon2d):
    """
    Map WOA basin codes to parent-ocean codes:
      1 Atlantic, 2 Pacific, 3 Indian, 4 Southern, 5 Arctic, 0 other/land.
    WOA13/18 marginal-sea codes are folded into their parent ocean.
    """
    known = {
        1: 1,   # Atlantic
        6: 1,   # Mediterranean -> Atlantic
        7: 1,   # Black Sea -> Atlantic
        8: 1,   # Baltic -> Atlantic
        11: 1,  # Hudson Bay -> Atlantic
        2: 2,   # Pacific
        3: 3,   # Indian
        9: 3,   # Red Sea -> Indian
        10: 3,  # Persian Gulf -> Indian
        4: 4,   # Southern Ocean
        5: 5,   # Arctic Ocean
    }

    out = np.zeros(basin2d.shape, dtype=int)
    finite = np.isfinite(basin2d)
    codes = np.rint(basin2d[finite]).astype(int)

    for ci in np.unique(codes):
        m = finite & (np.rint(basin2d) == ci)
        if ci == 0:
            p = 0
        elif ci in known:
            p = known[ci]
        else:
            lons = lon2d[m]
            lats = lat2d[m]
            if lons.size == 0:
                p = 0
            else:
                lon_mean = np.nanmean(np.mod(lons + 180.0, 360.0) - 180.0)
                lat_mean = np.nanmean(lats)
                if lat_mean > 60 and lon_mean > 0:
                    p = 5
                elif lon_mean < -70 or lon_mean > 150:
                    p = 2
                elif lon_mean > 30 and lat_mean < 35:
                    p = 3
                else:
                    p = 1
        out[m] = p
    return out


def map_to_woa(lat2d, lon2d, field_woa, woa_lat, woa_lon):
    lat = np.asarray(lat2d, dtype=float)
    lon = np.asarray(lon2d, dtype=float)
    lon = np.mod(lon + 180.0, 360.0) - 180.0

    dlat = woa_lat[1] - woa_lat[0]
    dlon = woa_lon[1] - woa_lon[0]

    i = np.rint((lat - woa_lat[0]) / dlat).astype(int)
    j = np.rint((lon - woa_lon[0]) / dlon).astype(int)
    i = np.clip(i, 0, len(woa_lat) - 1)
    j = np.clip(j, 0, len(woa_lon) - 1)
    return field_woa[i, j]


def compute_stats(no3, po4, vol, basin_code):
    no3 = np.asarray(no3, dtype=float)
    po4 = np.asarray(po4, dtype=float)
    vol = np.asarray(vol, dtype=float)
    basin_code = np.asarray(basin_code)

    valid = (
        np.isfinite(no3)
        & np.isfinite(po4)
        & np.isfinite(vol)
        & (vol > 0)
        & (po4 > 0)
    )

    ratio = np.full(no3.shape, np.nan, dtype=float)
    ratio[valid] = no3[valid] / po4[valid]

    v = vol[valid]
    if v.size == 0:
        nan = np.nan
        return {
            "mean_NP_ratio": nan,
            "mean_NP_ratio_volw": nan,
            "mean_no3": nan,
            "mean_po4": nan,
            "ratio_of_means": nan,
            "ratio_of_volw_means": nan,
            "slope": nan,
            "intercept": nan,
            "R2": nan,
            "pct_14_18": nan,
            "basin_means": {"Pacific": nan, "Atlantic": nan, "Indian": nan},
            "basin_means_volw": {"Pacific": nan, "Atlantic": nan, "Indian": nan},
            "pct_below12_global": nan,
            "pct_below12_Pacific": nan,
            "pct_below12_Indian": nan,
            "Redfield_status": "NA",
            "Redfield_status_volw": "NA",
        }

    # Headline: arithmetic mean of the cell-by-cell N:P ratio.
    mean_NP_ratio = np.nanmean(ratio[valid])
    mean_NP_ratio_volw = np.sum(ratio[valid] * v) / np.sum(v)

    mean_no3 = np.nanmean(no3[valid])
    mean_po4 = np.nanmean(po4[valid])
    ratio_of_means = mean_no3 / mean_po4 if mean_po4 > 0 else np.nan

    mean_no3_volw = np.sum(no3[valid] * v) / np.sum(v)
    mean_po4_volw = np.sum(po4[valid] * v) / np.sum(v)
    ratio_of_volw_means = (
        mean_no3_volw / mean_po4_volw if mean_po4_volw > 0 else np.nan
    )

    # Ordinary least squares, each valid horizontal cell weighted equally.
    x = po4[valid]
    y = no3[valid]
    if x.size >= 2:
        slope, intercept = np.polyfit(x, y, 1)
        yhat = slope * x + intercept
        ss_res = np.sum((y - yhat) ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    else:
        slope = intercept = R2 = np.nan

    # Volume fraction with 14 < N:P < 18.
    pct_14_18 = (
        np.sum(vol[valid & (ratio > 14) & (ratio < 18)]) / np.sum(vol[valid]) * 100.0
    )

    basin_means = {}
    basin_means_volw = {}
    basin_pct_below12 = {}

    for bname, bcode in [("Pacific", 2), ("Atlantic", 1), ("Indian", 3)]:
        mb = valid & (basin_code == bcode)
        if np.any(mb):
            basin_means[bname] = np.nanmean(ratio[mb])
            basin_means_volw[bname] = np.sum(ratio[mb] * vol[mb]) / np.sum(vol[mb])
            basin_pct_below12[bname] = (
                np.sum(vol[mb & (ratio < 12)]) / np.sum(vol[mb]) * 100.0
            )
        else:
            basin_means[bname] = np.nan
            basin_means_volw[bname] = np.nan
            basin_pct_below12[bname] = np.nan

    pct_below12_global = (
        np.sum(vol[valid & (ratio < 12)]) / np.sum(vol[valid]) * 100.0
    )

    def redfield_status(value):
        if not np.isfinite(value):
            return "NA"
        if value > 16.0:
            return "Above Redfield"
        if value < 16.0:
            return "Below Redfield"
        return "At Redfield"

    return {
        "mean_NP_ratio": mean_NP_ratio,
        "mean_NP_ratio_volw": mean_NP_ratio_volw,
        "mean_no3": mean_no3,
        "mean_po4": mean_po4,
        "ratio_of_means": ratio_of_means,
        "ratio_of_volw_means": ratio_of_volw_means,
        "slope": slope,
        "intercept": intercept,
        "R2": R2,
        "pct_14_18": pct_14_18,
        "basin_means": basin_means,
        "basin_means_volw": basin_means_volw,
        "pct_below12_global": pct_below12_global,
        "pct_below12_Pacific": basin_pct_below12["Pacific"],
        "pct_below12_Indian": basin_pct_below12["Indian"],
        "Redfield_status": redfield_status(mean_NP_ratio),
        "Redfield_status_volw": redfield_status(mean_NP_ratio_volw),
    }


# -----------------------------------------------------------------------------
# WOA23 nitrate, phosphate, and basin mask
# -----------------------------------------------------------------------------
no3_da, _ = open_woa_var("nitrate", "n")
po4_da, _ = open_woa_var("phosphate", "p")
no3_da, po4_da = xr.align(no3_da, po4_da, join="inner")

depth_dim = "depth" if "depth" in no3_da.dims else "lev"
lat_dim = "lat"
lon_dim = "lon"

woa_depth = np.asarray(no3_da[depth_dim].values, dtype=float)
woa_lat = np.asarray(no3_da[lat_dim].values, dtype=float)
woa_lon = np.asarray(no3_da[lon_dim].values, dtype=float)

no3_woa = no3_da.transpose(depth_dim, lat_dim, lon_dim).values.astype(float)
po4_woa = po4_da.transpose(depth_dim, lat_dim, lon_dim).values.astype(float)

# Layer weights from 200-1000 m using mid-point layer bounds.
woa_overlap = bounds_overlap_from_centers(woa_depth, 200.0, 1000.0)
woa_valid = np.isfinite(no3_woa) & np.isfinite(po4_woa)
woa_w = woa_overlap[:, None, None] * woa_valid
woa_sumw = np.sum(woa_w, axis=0)

with np.errstate(invalid="ignore", divide="ignore"):
    no3_woa_layer = np.sum(no3_woa * woa_w, axis=0) / woa_sumw
    po4_woa_layer = np.sum(po4_woa * woa_w, axis=0) / woa_sumw

no3_woa_layer[woa_sumw <= 0] = np.nan
po4_woa_layer[woa_sumw <= 0] = np.nan

# WOA grid-cell area.
dlon_deg = np.abs(woa_lon[1] - woa_lon[0]) if woa_lon.size > 1 else 1.0
dlat_deg = np.abs(woa_lat[1] - woa_lat[0]) if woa_lat.size > 1 else 1.0
dlon = np.deg2rad(dlon_deg)
lat_lower = np.deg2rad(woa_lat - dlat_deg / 2.0)
lat_upper = np.deg2rad(woa_lat + dlat_deg / 2.0)
area_lat = R_EARTH**2 * dlon * (np.sin(lat_upper) - np.sin(lat_lower))
area2d_woa = area_lat[:, None] * np.ones((1, woa_lon.size))

vol_woa_layer = area2d_woa * woa_sumw

# WOA basin mask.
basin_mask_path = os.path.join(WOA_MASKS, "basinmask_01.msk")
basin2d_woa_codes, basin3d_woa_codes = load_woa_basin_mask(
    basin_mask_path, woa_depth, woa_lat, woa_lon
)

lat2d_woa = np.repeat(woa_lat[:, None], woa_lon.size, axis=1)
lon2d_woa = np.repeat(woa_lon[None, :], woa_lat.size, axis=0)
basin_parent_woa = build_parent_codes(basin2d_woa_codes, lat2d_woa, lon2d_woa)

sources = {}
sources["WOA23"] = {
    "no3": no3_woa_layer,
    "po4": po4_woa_layer,
    "vol": vol_woa_layer,
    "basin": basin_parent_woa,
    "lat2d": lat2d_woa,
    "lon2d": lon2d_woa,
}

# -----------------------------------------------------------------------------
# CMIP6 historical no3, po4, and volcello
# -----------------------------------------------------------------------------
for model in CMIP_MODELS:
    files_no3_all = find_var_files(model, "no3")
    variants = sorted(set(get_variant(f) for f in files_no3_all))
    variant = next((v for v in variants if v is not None), None)
    files_no3 = filter_files_by_variant(files_no3_all, variant)

    ds_no3 = open_cmip_dataset(files_no3)
    da_no3 = ds_no3["no3"] if "no3" in ds_no3 else ds_no3[[v for v in ds_no3.data_vars if v.startswith("no3")][0]]
    da_no3 = period_mean(da_no3, 1971, 2000)

    files_po4_all = find_var_files(model, "po4")
    files_po4 = filter_files_by_variant(files_po4_all, variant)
    ds_po4 = open_cmip_dataset(files_po4)
    da_po4 = ds_po4["po4"] if "po4" in ds_po4 else ds_po4[[v for v in ds_po4.data_vars if v.startswith("po4")][0]]
    da_po4 = period_mean(da_po4, 1971, 2000)

    files_vol_all = find_var_files(model, "volcello")
    files_vol = filter_files_by_variant(files_vol_all, variant)
    ds_vol = open_cmip_dataset(files_vol)
    da_vol = ds_vol["volcello"]
    da_vol = period_mean(da_vol, 1971, 2000)

    da_no3, da_po4, da_vol = xr.align(da_no3, da_po4, da_vol, join="inner")

    # CMIP6: mol m-3 -> umol kg-1 at constant density.
    da_no3 = da_no3 * (1.0e6 / RHO)
    da_po4 = da_po4 * (1.0e6 / RHO)

    depth_dim_name = "lev"
    frac_vals = get_depth_frac(ds_no3, depth_dim=depth_dim_name,
                               zmin=200.0, zmax=1000.0)
    ds_lev = np.asarray(ds_no3[depth_dim_name].values, dtype=float).ravel()
    target_lev = np.asarray(da_no3[depth_dim_name].values, dtype=float).ravel()

    if frac_vals.size != target_lev.size and frac_vals.size == ds_lev.size:
        # Align the depth fractions onto the (possibly restricted) levels
        # of the co-aligned data arrays.
        frac_ser = pd.Series(frac_vals, index=ds_lev)
        frac_vals = frac_ser.reindex(target_lev).to_numpy(dtype=float)

    frac_da = xr.DataArray(frac_vals, dims=[depth_dim_name],
                           coords={depth_dim_name: target_lev})

    valid = (
        np.isfinite(da_no3)
        & np.isfinite(da_po4)
        & np.isfinite(da_vol)
        & (da_vol > 0)
    )
    w = da_vol * frac_da * valid
    sumw = w.sum(depth_dim_name)

    no3_layer = ((da_no3 * w).sum(depth_dim_name) / sumw).where(sumw > 0)
    po4_layer = ((da_po4 * w).sum(depth_dim_name) / sumw).where(sumw > 0)
    vol_layer = sumw

    no3_2d = no3_layer.compute().values
    po4_2d = po4_layer.compute().values
    vol_2d = vol_layer.compute().values

    lat2d, lon2d = get_latlon2d(ds_no3, model)
    basin_model = map_to_woa(lat2d, lon2d, basin_parent_woa, woa_lat, woa_lon)

    sources[model] = {
        "no3": no3_2d,
        "po4": po4_2d,
        "vol": vol_2d,
        "basin": basin_model,
        "lat2d": lat2d,
        "lon2d": lon2d,
    }

    ds_no3.close()
    ds_po4.close()
    ds_vol.close()

# -----------------------------------------------------------------------------
# Statistics
# -----------------------------------------------------------------------------
stats = {}
for name in SOURCE_ORDER:
    s = sources[name]
    stats[name] = compute_stats(s["no3"], s["po4"], s["vol"], s["basin"])

# Model-minus-WOA differences.
for model in CMIP_MODELS:
    for key in [
        "mean_NP_ratio",
        "mean_NP_ratio_volw",
        "ratio_of_means",
        "ratio_of_volw_means",
        "slope",
        "intercept",
        "R2",
        "pct_14_18",
        "pct_below12_global",
        "pct_below12_Pacific",
        "pct_below12_Indian",
    ]:
        stats[model][f"delta_{key}_vs_WOA"] = (
            stats[model][key] - stats["WOA23"][key]
        )

    for b in ["Pacific", "Atlantic", "Indian"]:
        stats[model][f"delta_{b}_mean_NP_vs_WOA"] = (
            stats[model]["basin_means"][b] - stats["WOA23"]["basin_means"][b]
        )
        stats[model][f"delta_{b}_mean_NP_volw_vs_WOA"] = (
            stats[model]["basin_means_volw"][b]
            - stats["WOA23"]["basin_means_volw"][b]
        )

# Summary table.
rows = []
for name in SOURCE_ORDER:
    s = stats[name]
    rows.append(
        {
            "source": name,
            "mean_NP_ratio": s["mean_NP_ratio"],
            "mean_NP_ratio_volw": s["mean_NP_ratio_volw"],
            "ratio_of_means": s["ratio_of_means"],
            "ratio_of_volw_means": s["ratio_of_volw_means"],
            "slope": s["slope"],
            "intercept": s["intercept"],
            "R2": s["R2"],
            "pct_14_18": s["pct_14_18"],
            "Pacific_mean_NP": s["basin_means"]["Pacific"],
            "Atlantic_mean_NP": s["basin_means"]["Atlantic"],
            "Indian_mean_NP": s["basin_means"]["Indian"],
            "Pacific_mean_NP_volw": s["basin_means_volw"]["Pacific"],
            "Atlantic_mean_NP_volw": s["basin_means_volw"]["Atlantic"],
            "Indian_mean_NP_volw": s["basin_means_volw"]["Indian"],
            "pct_below12_global": s["pct_below12_global"],
            "pct_below12_Pacific": s["pct_below12_Pacific"],
            "pct_below12_Indian": s["pct_below12_Indian"],
            "Redfield_status": s["Redfield_status"],
            "Redfield_status_volw": s["Redfield_status_volw"],
        }
    )

summary_df = pd.DataFrame(rows).set_index("source")

summary_df_rounded = summary_df.copy()
for col in summary_df_rounded.columns:
    if col in ("Redfield_status", "Redfield_status_volw"):
        continue
    if col == "R2":
        summary_df_rounded[col] = summary_df_rounded[col].round(3)
    elif "pct" in col:
        summary_df_rounded[col] = summary_df_rounded[col].round(1)
    else:
        summary_df_rounded[col] = summary_df_rounded[col].round(2)

print(summary_df_rounded.to_string())
