import os
import glob
import itertools
import numpy as np
import pandas as pd
import xarray as xr

# -----------------------------------------------------------------------------
# Constants and paths
# -----------------------------------------------------------------------------
WOA_ROOT = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
WOA_MASK_ROOT = "/data1/lorraine-20419/ocean_data/WOA23/masks"
CMIP6_ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"

WOA_NITRATE_PATH = os.path.join(WOA_ROOT, "nitrate", "all_time_periods", "woa23_all_n00_01.nc")
WOA_PHOSPHATE_PATH = os.path.join(WOA_ROOT, "phosphate", "all_time_periods", "woa23_all_p00_01.nc")
BASIN_MASK_PATH = os.path.join(WOA_MASK_ROOT, "basinmask_01.msk")

R_EARTH = 6371000.0          # m
RHO = 1025.0                # kg m-3, prescribed constant for conversion/regression

ORDER = ["WOA23", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
MODELS = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

MODEL_INFO = {
    "CESM2": {
        "y_dim": "nlat", "x_dim": "nlon",
        "lat_coord": "lat", "lon_coord": "lon",
        "depth_coord": "lev", "depth_scale": 0.01,
        "calendar": "noleap", "ocean_mask_from": "sftof",
    },
    "GFDL-ESM4": {
        "y_dim": "lat", "x_dim": "lon",
        "lat_coord": "lat", "lon_coord": "lon",
        "depth_coord": "lev", "depth_scale": 1.0,
        "calendar": "noleap",
    },
    "MPI-ESM1-2-LR": {
        "y_dim": "j", "x_dim": "i",
        "lat_coord": "latitude", "lon_coord": "longitude",
        "depth_coord": "lev", "depth_scale": 1.0,
        "calendar": "proleptic_gregorian",
    },
}

# -----------------------------------------------------------------------------
# General helpers
# -----------------------------------------------------------------------------
def find_data_var(ds, candidates):
    for c in candidates:
        if c in ds.data_vars:
            return c
    for v in ds.data_vars:
        for c in candidates:
            if c.lower() in v.lower():
                return v
    if len(ds.data_vars) == 1:
        return list(ds.data_vars)[0]
    raise KeyError(f"Could not find data variable among {candidates}; have {list(ds.data_vars)}")


def get_dim_name(da, candidates):
    for c in candidates:
        if c in da.dims:
            return c
    for c in candidates:
        if c in da.coords:
            return c
    raise KeyError(f"None of {candidates} found in dims {list(da.dims)} coords {list(da.coords)}")


def woa_to_umol_kg(da):
    """Return WOA field in micromoles per kilogram, converting only if metadata require it."""
    u = str(da.attrs.get("units", "")).lower()
    da_out = da
    if ("micro" in u) or ("umol" in u) or ("µmol" in u) or ("micromol" in u):
        pass
    elif ("mol" in u) and ("kg" in u):
        da_out = da * 1.0e6
    da_out.attrs["units"] = "umol kg-1"
    return da_out


def great_circle_dist(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(np.deg2rad, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0) ** 2
    return 2.0 * np.arcsin(np.sqrt(a))


def read_msk_rows(path):
    """Read whitespace/comma separated numeric rows, skipping non-numeric headers."""
    rows = []
    with open(path, "r", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#") or s.startswith("%"):
                continue
            s = s.replace(",", " ")
            parts = s.split()
            if len(parts) == 0:
                continue
            try:
                vals = [float(p) for p in parts]
            except ValueError:
                continue
            rows.append(vals)
    return rows


def map_basin_to_latlon(lat2d, lon2d, basin_name_grid, mask_lat, mask_lon):
    lon_mod = ((np.asarray(lon2d) + 180.0) % 360.0) - 180.0
    li = np.round(np.asarray(lat2d) - mask_lat[0]).astype(int)
    lj = np.round(lon_mod - mask_lon[0]).astype(int)
    li = np.clip(li, 0, len(mask_lat) - 1)
    lj = np.clip(lj, 0, len(mask_lon) - 1)
    return basin_name_grid[li, lj]


# -----------------------------------------------------------------------------
# Load and fold WOA23 basin mask
# -----------------------------------------------------------------------------
def load_basin_mask(path):
    rows = read_msk_rows(path)
    if len(rows) == 0:
        raise ValueError(f"No numeric rows in basin mask {path}")

    sample_rows = [r[:6] for r in rows if len(r) >= 6]
    if len(sample_rows) == 0:
        sample_rows = [r[:3] for r in rows if len(r) >= 3]
    sample = np.array(sample_rows, dtype=float)

    n_unique = []
    for j in range(sample.shape[1]):
        vals = sample[:, j][np.isfinite(sample[:, j])]
        n_unique.append(len(np.unique(vals)))

    lat_col = None
    lon_col = None
    for j in range(sample.shape[1]):
        vals = sample[:, j]
        if n_unique[j] > 100 and np.nanmax(np.abs(vals)) <= 90.5:
            lat_col = j
            break
    for j in range(sample.shape[1]):
        if j == lat_col:
            continue
        vals = sample[:, j]
        if n_unique[j] > 100 and np.nanmax(np.abs(vals)) <= 180.5:
            lon_col = j
            break

    if lat_col is None or lon_col is None:
        start = 1 if np.nanmax(sample[:, 0]) > 1000.0 else 0
        if np.nanmax(np.abs(sample[:, start])) <= 90.5:
            lat_col, lon_col = start, start + 1
        else:
            lon_col, lat_col = start, start + 1

    lev_start = max(lat_col, lon_col) + 1

    lat = []
    lon = []
    surf = []
    for r in rows:
        if len(r) > lev_start:
            lat.append(r[lat_col])
            lon.append(r[lon_col])
            surf.append(r[lev_start])  # surface level label

    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    surf = np.asarray(surf, dtype=float)

    valid = np.isfinite(lat) & np.isfinite(lon) & np.isfinite(surf)
    lat = lat[valid]
    lon = lon[valid]
    surf = surf[valid]

    if np.nanmax(lon) > 180.0:
        lon = ((lon + 180.0) % 360.0) - 180.0

    lat_centers = np.arange(-89.5, 90.0, 1.0)
    lon_centers = np.arange(-179.5, 180.0, 1.0)

    basin_grid = np.full((len(lat_centers), len(lon_centers)), np.nan)
    li = np.round(lat - lat_centers[0]).astype(int)
    lj = np.round(lon - lon_centers[0]).astype(int)
    ok = (li >= 0) & (li < len(lat_centers)) & (lj >= 0) & (lj < len(lon_centers))
    basin_grid[li[ok], lj[ok]] = surf[ok]

    lon2d_grid, lat2d_grid = np.meshgrid(lon_centers, lat_centers)

    codes = np.unique(basin_grid[np.isfinite(basin_grid)])
    if len(codes) == 0:
        raise ValueError("No basin codes found")

    code_stats = []
    for c in codes:
        m = basin_grid == c
        n = int(m.sum())
        if n == 0:
            continue
        lats = lat2d_grid[m]
        lons = lon2d_grid[m]
        ang = np.deg2rad(lons)
        meanlon = np.rad2deg(np.arctan2(np.mean(np.sin(ang)), np.mean(np.cos(ang))))
        code_stats.append({
            "code": float(c),
            "n": n,
            "medlat": float(np.median(lats)),
            "meanlat": float(np.mean(lats)),
            "meanlon": float(meanlon),
        })

    south_candidates = [s for s in code_stats if s["medlat"] < -55.0]
    south_code = max(south_candidates, key=lambda s: s["n"])["code"] if south_candidates else None

    remaining = [s for s in code_stats if s["code"] != south_code]
    arctic_candidates = [s for s in remaining if s["medlat"] > 65.0]
    arctic_code = max(arctic_candidates, key=lambda s: s["n"])["code"] if arctic_candidates else None

    remaining = [s for s in code_stats if s["code"] not in (south_code, arctic_code)]
    main_stats = sorted(remaining, key=lambda s: s["n"], reverse=True)[:3]

    refs = {
        "Atlantic": (-30.0, 20.0),
        "Pacific": (180.0, 0.0),
        "Indian": (80.0, 0.0),
    }
    main_names = ["Atlantic", "Pacific", "Indian"]

    best_perm = None
    best_cost = np.inf
    for perm in itertools.permutations(main_names):
        cost = 0.0
        for s, name in zip(main_stats, perm):
            rlon, rlat = refs[name]
            cost += great_circle_dist(s["meanlon"], s["meanlat"], rlon, rlat)
        if cost < best_cost:
            best_cost = cost
            best_perm = perm

    code_to_name = {}
    for s, name in zip(main_stats, best_perm):
        code_to_name[s["code"]] = name

    for s in code_stats:
        c = s["code"]
        if c in code_to_name:
            continue
        if c == south_code:
            code_to_name[c] = "Southern"
            continue
        if c == arctic_code:
            code_to_name[c] = "Arctic"
            continue

        lon = s["meanlon"]
        lat = s["meanlat"]

        # Fold common marginal seas into their parent ocean.
        if 95.0 <= lon <= 130.0 and 0.0 <= lat <= 45.0:
            name = "Pacific"
        elif 30.0 <= lon <= 60.0 and 5.0 <= lat <= 35.0:
            name = "Indian"
        elif -10.0 <= lon <= 45.0 and 30.0 <= lat <= 70.0:
            name = "Atlantic"
        elif -100.0 <= lon <= -60.0 and 10.0 <= lat <= 35.0:
            name = "Atlantic"
        else:
            dists = {name: great_circle_dist(lon, lat, *refs[name]) for name in main_names}
            name = min(dists, key=dists.get)
        code_to_name[c] = name

    basin_name_grid = np.full(basin_grid.shape, "Land", dtype=object)
    for c, name in code_to_name.items():
        basin_name_grid[basin_grid == c] = name

    return basin_name_grid, lat_centers, lon_centers


basin_name_grid, mask_lat, mask_lon = load_basin_mask(BASIN_MASK_PATH)

# -----------------------------------------------------------------------------
# Load WOA23 annual objectively analysed nitrate and phosphate
# -----------------------------------------------------------------------------
def load_woa_source(basin_name_grid, mask_lat, mask_lon):
    ds_n = xr.open_dataset(WOA_NITRATE_PATH, decode_times=False)
    ds_p = xr.open_dataset(WOA_PHOSPHATE_PATH, decode_times=False)

    var_n = find_data_var(ds_n, ["n_an", "nitrate_an", "nitrate", "n"])
    var_p = find_data_var(ds_p, ["p_an", "phosphate_an", "phosphate", "p"])

    da_n = ds_n[var_n]
    da_p = ds_p[var_p]

    if "time" in da_n.dims:
        da_n = da_n.isel(time=0)
    if "time" in da_p.dims:
        da_p = da_p.isel(time=0)

    da_n = da_n.squeeze(drop=True)
    da_p = da_p.squeeze(drop=True)

    depth_name = get_dim_name(da_n, ["depth", "deptht", "lev", "z"])
    lat_name = get_dim_name(da_n, ["lat", "latitude", "y"])
    lon_name = get_dim_name(da_n, ["lon", "longitude", "x"])

    da_n = da_n.transpose(depth_name, lat_name, lon_name)
    da_p = da_p.transpose(depth_name, lat_name, lon_name)

    lat = da_n[lat_name].values.astype(float)
    if lat[0] > lat[-1]:
        da_n = da_n.isel({lat_name: slice(None, None, -1)})
        da_p = da_p.isel({lat_name: slice(None, None, -1)})
        lat = lat[::-1]

    lon = da_n[lon_name].values.astype(float)
    if lon[0] > lon[-1]:
        da_n = da_n.isel({lon_name: slice(None, None, -1)})
        da_p = da_p.isel({lon_name: slice(None, None, -1)})
        lon = lon[::-1]

    depth = da_n[depth_name].values.astype(float)
    if depth[0] > depth[-1]:
        da_n = da_n.isel({depth_name: slice(None, None, -1)})
        da_p = da_p.isel({depth_name: slice(None, None, -1)})
        depth = depth[::-1]

    da_n, da_p = xr.align(da_n, da_p, join="inner")
    da_n = da_n.transpose(depth_name, lat_name, lon_name)
    da_p = da_p.transpose(depth_name, lat_name, lon_name)

    da_n = woa_to_umol_kg(da_n)
    da_p = woa_to_umol_kg(da_p)

    sel = np.where((depth >= 200.0) & (depth <= 1000.0))[0]
    if len(sel) == 0:
        raise ValueError("No WOA23 standard levels in 200-1000 m")

    bounds = np.empty(len(depth) + 1, dtype=float)
    bounds[1:-1] = 0.5 * (depth[:-1] + depth[1:])
    bounds[0] = depth[0] - 0.5 * (depth[1] - depth[0])
    bounds[-1] = depth[-1] + 0.5 * (depth[-1] - depth[-2])
    thick_full = np.diff(bounds)
    thick = thick_full[sel]

    da_n = da_n.isel({depth_name: sel})
    da_p = da_p.isel({depth_name: sel})

    no3 = da_n.values.astype(float)
    po4 = da_p.values.astype(float)

    no3[~np.isfinite(no3) | (np.abs(no3) > 1.0e30)] = np.nan
    po4[~np.isfinite(po4) | (np.abs(po4) > 1.0e30)] = np.nan

    dlat = np.deg2rad(np.abs(np.mean(np.diff(lat)))) if len(lat) > 1 else np.deg2rad(1.0)
    dlon = np.deg2rad(np.abs(np.mean(np.diff(lon)))) if len(lon) > 1 else np.deg2rad(1.0)

    area1d = (R_EARTH ** 2) * np.cos(np.deg2rad(lat)) * dlat * dlon
    area2d = area1d[:, None] * np.ones((1, len(lon)), dtype=float)
    vol = area2d[None, :, :] * thick[:, None, None]

    lon2d, lat2d = np.meshgrid(lon, lat)
    basin2d = map_basin_to_latlon(lat2d, lon2d, basin_name_grid, mask_lat, mask_lon)

    return {
        "no3": no3,
        "po4": po4,
        "vol": vol,
        "area2d": area2d,
        "lat2d": lat2d,
        "lon2d": lon2d,
        "lat1d": lat,
        "lon1d": lon,
        "basin2d": basin2d,
        "depth_levels_m": depth[sel],
        "volume_source": "WOA23 1-degree spherical cell area times depth-level bounds",
        "density_assumption": "WOA23 native umol kg-1",
    }


# -----------------------------------------------------------------------------
# Load CMIP6 historical no3 and po4
# -----------------------------------------------------------------------------
def parse_variant(fn):
    for tok in os.path.basename(fn).split("_"):
        if tok.startswith("r") and ("i" in tok) and ("p" in tok) and ("f" in tok):
            return tok
    return None


def open_cmip6_var(model, var):
    model_dir = os.path.join(CMIP6_ROOT, model, var)
    files = sorted(glob.glob(os.path.join(model_dir, f"{var}_*.nc")))
    if len(files) == 0:
        files = sorted(glob.glob(os.path.join(CMIP6_ROOT, model, f"{var}_*.nc")))
    if len(files) == 0:
        raise FileNotFoundError(f"No files found for {model} {var} under {model_dir}")

    variants = [parse_variant(f) for f in files]
    variant = next((v for v in variants if v is not None), None)
    if variant is not None:
        files = [f for f, v in zip(files, variants) if v == variant]

    if len(files) == 1:
        ds = xr.open_dataset(files[0], use_cftime=True, decode_times=True)
    else:
        ds = xr.open_mfdataset(
            files,
            combine="by_coords",
            use_cftime=True,
            decode_times=True,
        )

    vname = find_data_var(ds, [var])
    return ds[vname]


def time_mean_1971_2000(da):
    """Day-weighted 1971-2000 mean in the model calendar, or squeeze a static field."""
    if "time" not in da.dims:
        return da.squeeze(drop=True)

    yrs = da["time"].dt.year
    sel = (yrs >= 1971) & (yrs <= 2000)
    da = da.isel(time=sel)
    if da.sizes["time"] == 0:
        raise ValueError("No time steps in 1971-2000")

    w = da["time"].dt.days_in_month.astype(float)
    num = (da.fillna(0.0) * w).sum("time")
    den = (da.notnull() * w).sum("time")
    out = num / den
    out.attrs = da.attrs
    return out


def load_cmip6_source(model, basin_name_grid, mask_lat, mask_lon):
    info = MODEL_INFO[model]

    no3_raw = open_cmip6_var(model, "no3")
    po4_raw = open_cmip6_var(model, "po4")
    vol_raw = open_cmip6_var(model, "volcello")
    area_raw = open_cmip6_var(model, "areacello")

    no3_t = time_mean_1971_2000(no3_raw)
    po4_t = time_mean_1971_2000(po4_raw)
    vol_t = time_mean_1971_2000(vol_raw)
    area = area_raw.squeeze(drop=True)

    if model == "CESM2":
        sftof = open_cmip6_var(model, "sftof").squeeze(drop=True)
        area = area.where(sftof > 0.0)

    factor = 1.0e6 / RHO
    no3 = no3_t * factor
    po4 = po4_t * factor
    no3.attrs["units"] = "umol kg-1"
    po4.attrs["units"] = "umol kg-1"
    no3.attrs["density_assumption"] = "1025 kg m-3"
    po4.attrs["density_assumption"] = "1025 kg m-3"

    lev_name = info["depth_coord"]
    lev_vals = no3[lev_name].values.astype(float) * float(info["depth_scale"])
    sel_lev = np.where((lev_vals >= 200.0) & (lev_vals <= 1000.0))[0]
    if len(sel_lev) == 0:
        raise ValueError(f"No {model} levels in 200-1000 m")

    no3 = no3.isel({lev_name: sel_lev})
    po4 = po4.isel({lev_name: sel_lev})
    vol = vol_t.isel({lev_name: sel_lev})

    no3, po4, vol = xr.align(no3, po4, vol, join="inner")

    y_dim = no3.dims[1]
    x_dim = no3.dims[2]
    if area.dims != (y_dim, x_dim):
        area = area.transpose(y_dim, x_dim)
    area2d = area.values.astype(float)

    lat_vals = no3[info["lat_coord"]].values
    lon_vals = no3[info["lon_coord"]].values

    if lat_vals.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon_vals, lat_vals)
        lat1d = np.asarray(lat_vals, dtype=float)
        lon1d = np.asarray(lon_vals, dtype=float)
    else:
        lat2d = np.asarray(lat_vals, dtype=float)
        lon2d = np.asarray(lon_vals, dtype=float)
        lat1d = None
        lon1d = None

    basin2d = map_basin_to_latlon(lat2d, lon2d, basin_name_grid, mask_lat, mask_lon)

    volume_source = (
        "volcello (Omon monthly, day-weighted 1971-2000 mean)"
        if model == "GFDL-ESM4"
        else "volcello (Ofx fixed field)"
    )

    return {
        "no3": no3.values.astype(float),
        "po4": po4.values.astype(float),
        "vol": vol.values.astype(float),
        "area2d": area2d,
        "lat2d": lat2d,
        "lon2d": lon2d,
        "lat1d": lat1d,
        "lon1d": lon1d,
        "basin2d": basin2d,
        "depth_levels_m": lev_vals[sel_lev],
        "volume_source": volume_source,
        "density_assumption": "1025 kg m-3",
    }


# -----------------------------------------------------------------------------
# Build the four sources
# -----------------------------------------------------------------------------
sources = {}
sources["WOA23"] = load_woa_source(basin_name_grid, mask_lat, mask_lon)

for model in MODELS:
    sources[model] = load_cmip6_source(model, basin_name_grid, mask_lat, mask_lon)

# -----------------------------------------------------------------------------
# Diagnostics on each native grid
# -----------------------------------------------------------------------------
def compute_source_stats(src, label):
    no3 = src["no3"]
    po4 = src["po4"]
    vol = src["vol"]
    basin2d = src["basin2d"]

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = no3 / po4

    valid = (
        np.isfinite(no3)
        & np.isfinite(po4)
        & (po4 > 0.0)
        & np.isfinite(vol)
        & (vol > 0.0)
    )

    x = po4[valid].astype(float)
    y = no3[valid].astype(float)
    w = vol[valid].astype(float)
    r = ratio[valid].astype(float)

    if len(x) >= 2:
        a, b = np.polyfit(x, y, 1)
        yhat = a * x + b
        ss_res = np.sum((y - yhat) ** 2)
        ss_tot = np.sum((y - y.mean()) ** 2)
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else np.nan
    else:
        a, b, r2 = np.nan, np.nan, np.nan

    total_w = np.sum(w)
    mean_np = np.sum(r * w) / total_w if total_w > 0.0 else np.nan
    mean_no3 = np.sum(no3[valid] * w) / total_w if total_w > 0.0 else np.nan
    mean_po4 = np.sum(po4[valid] * w) / total_w if total_w > 0.0 else np.nan
    ratio_of_means = mean_no3 / mean_po4 if mean_po4 != 0.0 else np.nan

    pct_14_18 = np.sum(w[(r > 14.0) & (r < 18.0)]) / total_w * 100.0 if total_w > 0.0 else np.nan
    pct_below12 = np.sum(w[r < 12.0]) / total_w * 100.0 if total_w > 0.0 else np.nan

    basin3 = np.broadcast_to(basin2d, no3.shape)

    basin_results = {}
    for bname in ["Pacific", "Atlantic", "Indian"]:
        bm = valid & (basin3 == bname)
        if np.any(bm):
            wb = vol[bm].astype(float)
            rb = ratio[bm].astype(float)
            basin_results[bname] = np.sum(rb * wb) / np.sum(wb)
            basin_results[bname + "_pct_below12"] = (
                np.sum(wb[rb < 12.0]) / np.sum(wb) * 100.0
            )
        else:
            basin_results[bname] = np.nan
            basin_results[bname + "_pct_below12"] = np.nan

    if np.isfinite(mean_np):
        if mean_np > 16.0:
            redfield_status = "above Redfield"
        elif mean_np < 16.0:
            redfield_status = "below Redfield"
        else:
            redfield_status = "at Redfield"
    else:
        redfield_status = "undefined"

    return {
        "source": label,
        "mean_np": mean_np,
        "ratio_of_means": ratio_of_means,
        "a": a,
        "b": b,
        "r2": r2,
        "pct_14_18": pct_14_18,
        "pct_below12": pct_below12,
        "basins": basin_results,
        "redfield_status": redfield_status,
        "n_valid_cells": int(np.sum(valid)),
        "total_valid_volume_m3": float(total_w),
    }


results = {}
for name in ORDER:
    results[name] = compute_source_stats(sources[name], name)

summary_rows = []
for name in ORDER:
    res = results[name]
    summary_rows.append({
        "source": name,
        "mean_NP_cellwise": res["mean_np"],
        "ratio_of_volume_weighted_means": res["ratio_of_means"],
        "OLS_a": res["a"],
        "OLS_b_umol_kg": res["b"],
        "R2": res["r2"],
        "pct_volume_14_18": res["pct_14_18"],
        "pct_volume_below12": res["pct_below12"],
        "Pacific_mean_NP": res["basins"]["Pacific"],
        "Atlantic_mean_NP": res["basins"]["Atlantic"],
        "Indian_mean_NP": res["basins"]["Indian"],
        "Pacific_pct_below12": res["basins"]["Pacific_pct_below12"],
        "Indian_pct_below12": res["basins"]["Indian_pct_below12"],
        "redfield_status": res["redfield_status"],
        "volume_source": sources[name]["volume_source"],
        "density_assumption": sources[name]["density_assumption"],
    })

summary_df = pd.DataFrame(summary_rows).set_index("source")

# -----------------------------------------------------------------------------
# Scatter data for Figure 1
# -----------------------------------------------------------------------------
scatter_data = {}
for name in ORDER:
    src = sources[name]
    no3 = src["no3"]
    po4 = src["po4"]
    vol = src["vol"]
    valid = (
        np.isfinite(no3)
        & np.isfinite(po4)
        & (po4 > 0.0)
        & np.isfinite(vol)
        & (vol > 0.0)
    )
    scatter_data[name] = {
        "po4": po4[valid].ravel(),
        "no3": no3[valid].ravel(),
    }

# -----------------------------------------------------------------------------
# Column layer-mean N:P minus 16 for Figure 2, with display-only regridding
# -----------------------------------------------------------------------------
def regular_to_pm180(field, lat, lon):
    field = np.asarray(field, dtype=float).copy()
    lat = np.asarray(lat, dtype=float).copy()
    lon = np.asarray(lon, dtype=float).copy()

    if lat[0] > lat[-1]:
        field = field[::-1, :]
        lat = lat[::-1]

    if lon.max() > 180.0:
        lon_mod = ((lon + 180.0) % 360.0) - 180.0
        order = np.argsort(lon_mod)
        lon = lon_mod[order]
        field = field[:, order]

    return field, lat, lon


def bin_to_regular(field, lat2d, lon2d, area2d):
    lat_t = np.arange(-89.5, 90.0, 1.0)
    lon_t = np.arange(-179.5, 180.0, 1.0)

    lon_mod = ((np.asarray(lon2d) + 180.0) % 360.0) - 180.0
    li = np.round(np.asarray(lat2d) - lat_t[0]).astype(int)
    lj = np.round(lon_mod - lon_t[0]).astype(int)

    valid = (
        np.isfinite(field)
        & np.isfinite(area2d)
        & (area2d > 0.0)
        & (li >= 0) & (li < len(lat_t))
        & (lj >= 0) & (lj < len(lon_t))
    )

    sum_w = np.zeros((len(lat_t), len(lon_t)), dtype=float)
    sum_v = np.zeros_like(sum_w)

    np.add.at(sum_w, (li[valid], lj[valid]), area2d[valid])
    np.add.at(sum_v, (li[valid], lj[valid]), field[valid] * area2d[valid])

    out = np.full_like(sum_w, np.nan)
    m = sum_w > 0.0
    out[m] = sum_v[m] / sum_w[m]
    return out, lat_t, lon_t


map_fields = {}
for name in ORDER:
    src = sources[name]
    no3 = src["no3"]
    po4 = src["po4"]
    vol = src["vol"]

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = no3 / po4

    valid = (
        np.isfinite(no3)
        & np.isfinite(po4)
        & (po4 > 0.0)
        & np.isfinite(vol)
        & (vol > 0.0)
    )

    num = np.nansum(np.where(valid, ratio * vol, 0.0), axis=0)
    den = np.nansum(np.where(valid, vol, 0.0), axis=0)
    layer_mean = np.full_like(num, np.nan, dtype=float)
    m = den > 0.0
    layer_mean[m] = num[m] / den[m]
    layer_minus16 = layer_mean - 16.0

    if name in ("WOA23", "GFDL-ESM4"):
        field, lat_d, lon_d = regular_to_pm180(layer_minus16, src["lat1d"], src["lon1d"])
    else:
        field, lat_d, lon_d = bin_to_regular(
            layer_minus16, src["lat2d"], src["lon2d"], src["area2d"]
        )

    map_fields[name] = {"field": field, "lat": lat_d, "lon": lon_d}

# Print a compact numerical summary for inspection.
with pd.option_context("display.max_columns", None, "display.width", 200):
    print(summary_df)
