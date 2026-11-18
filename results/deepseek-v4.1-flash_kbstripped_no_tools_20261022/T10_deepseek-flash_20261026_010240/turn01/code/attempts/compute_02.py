import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree
from scipy.interpolate import griddata

# ---------------------------------------------------------------------------
# Paths and metadata
# ---------------------------------------------------------------------------
woa_root = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
woa_masks_root = "/data1/lorraine-20419/ocean_data/WOA23/masks"

cmip6_hist_root = "/data1/lorraine-20419/benchmark_cmip6/historical"
cmip6_ssp_root = "/data1/lorraine-20419/benchmark_cmip6/ssp245"

woa_vars = {
    "oxygen": {
        "prefix": "o",
        "periods": {
            "1971-2000": {"subdir": "oxygen/1971-2000_climate_normal", "tag": "decav71A0"},
            "all": {"subdir": "oxygen/all_time_periods", "tag": "all"},
        },
    },
    "apparent_oxygen_utilization": {
        "prefix": "A",
        "periods": {
            "1971-2000": {"subdir": "apparent_oxygen_utilization/1971-2000_climate_normal", "tag": "decav71A0"},
            "all": {"subdir": "apparent_oxygen_utilization/all_time_periods", "tag": "all"},
        },
    },
    "percent_oxygen_saturation": {
        "prefix": "O",
        "periods": {
            "1971-2000": {"subdir": "percent_oxygen_saturation/1971-2000_climate_normal", "tag": "decav71A0"},
            "all": {"subdir": "percent_oxygen_saturation/all_time_periods", "tag": "all"},
        },
    },
}

model_info = {
    "CESM2": {
        "y_dim": "nlat", "x_dim": "nlon",
        "lat_coord": "lat", "lon_coord": "lon",
        "depth_coord": "lev", "depth_scale": 0.01,
    },
    "GFDL-ESM4": {
        "y_dim": "lat", "x_dim": "lon",
        "lat_coord": "lat", "lon_coord": "lon",
        "depth_coord": "lev", "depth_scale": 1.0,
    },
    "MPI-ESM1-2-LR": {
        "y_dim": "j", "x_dim": "i",
        "lat_coord": "latitude", "lon_coord": "longitude",
        "depth_coord": "lev", "depth_scale": 1.0,
    },
}

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
def woa_path(var, period, nn):
    info = woa_vars[var]["periods"][period]
    subdir = info["subdir"]
    tag = info["tag"]
    prefix = woa_vars[var]["prefix"]
    fname = f"woa23_{tag}_{prefix}{nn:02d}_01.nc"
    return os.path.join(woa_root, subdir, fname)


def _woa_field(ds):
    """Return the objectively-analysed field from a WOA netCDF, robust to naming."""
    if "an" in ds.data_vars:
        return ds["an"]
    return ds[list(ds.data_vars)[0]]


def load_woa_monthly(var, preferred="1971-2000"):
    periods = woa_vars[var]["periods"]
    if preferred in periods:
        paths = [woa_path(var, preferred, nn) for nn in range(1, 13)]
        if all(os.path.exists(p) for p in paths):
            period = preferred
        elif "all" in periods:
            period = "all"
        else:
            period = preferred
    else:
        period = "all" if "all" in periods else list(periods.keys())[0]

    paths = [woa_path(var, period, nn) for nn in range(1, 13)]
    if not all(os.path.exists(p) for p in paths):
        raise FileNotFoundError(f"Missing WOA monthly files for {var} period {period}")

    ds_list = [xr.open_dataset(p, decode_times=False) for p in paths]
    da = xr.concat([_woa_field(d) for d in ds_list], dim="time")
    da = da.assign_coords(time=np.arange(1, 13))
    da = da.where(da < 1e19)

    days = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], dtype=float)
    w = xr.DataArray(days / days.sum(), dims="time", coords={"time": da.time})
    annual = (da * w).sum("time") / w.where(da.notnull()).sum("time")
    return annual, ds_list[0], period


def get_depth_bounds(ds):
    if "depth_bnds" in ds:
        b = ds["depth_bnds"].values
        if b.shape[0] == 2 and b.shape[1] == len(ds["depth"]):
            b = b.T
        return b
    elif "depth" in ds and "bounds" in ds["depth"].attrs:
        bname = ds["depth"].attrs["bounds"]
        b = ds[bname].values
        if b.shape[0] == 2 and b.shape[1] == len(ds["depth"]):
            b = b.T
        return b
    else:
        d = ds["depth"].values
        edges = np.zeros(len(d) + 1)
        edges[1:-1] = 0.5 * (d[:-1] + d[1:])
        edges[0] = 0.0
        edges[-1] = d[-1] + 0.5 * (d[-1] - d[-2])
        return np.column_stack([edges[:-1], edges[1:]])


def read_msk(path):
    """
    Robustly read a WOA .msk text table into a ragged list of float rows.

    The header of these files varies (comment lines, an optional '#', column
    names containing 'Latitude'/'Longitude').  Instead of trusting a fixed
    header offset, every line is scanned for numeric tokens; the first two
    tokens (assumed Latitude/Longitude) must be plausible geographic
    coordinates.  Empty fields (blanks, or levels below the sea floor) simply
    contribute no tokens, which is what we want since we take the first valid
    value downstream.
    """
    rows = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            vals = []
            for tok in s.split():
                try:
                    vals.append(float(tok.lstrip("+")))
                except ValueError:
                    continue
            if len(vals) < 3:
                continue
            la, lo = vals[0], vals[1]
            if not (-90.0 <= la <= 90.0 and -180.0 <= lo <= 180.0):
                continue
            rows.append(vals)
    return rows


def find_cmip6_files(root, model, var, exp):
    d = os.path.join(root, model, var)
    if not os.path.isdir(d):
        return []
    files = sorted(glob.glob(os.path.join(d, f"{var}_*.nc")))
    files = [f for f in files if f"_{exp}_" in os.path.basename(f)]
    return files


def get_period_mean(root, model, var, exp, start, end, allow_fallback=False):
    files = find_cmip6_files(root, model, var, exp)
    if len(files) == 0 and allow_fallback and exp != "historical":
        files = find_cmip6_files(cmip6_hist_root, model, var, "historical")
    if len(files) == 0:
        raise FileNotFoundError(f"No files for {model} {exp} {var}")

    ds0 = xr.open_dataset(files[0], use_cftime=True, decode_times=True)
    has_time = "time" in ds0.dims

    if has_time:
        ds = xr.open_mfdataset(files, combine="by_coords", use_cftime=True,
                               decode_times=True, chunks={"time": 120})
    else:
        ds = xr.open_dataset(files[0], use_cftime=True, decode_times=True)

    if var in ds:
        da = ds[var]
    else:
        candidates = [v for v in ds.data_vars if var in v.lower()]
        if candidates:
            da = ds[candidates[0]]
        else:
            da = ds[list(ds.data_vars)[0]]

    da = da.where(da < 1e19)
    if has_time:
        da = da.sel(time=slice(start, end))
        days = da.time.dt.days_in_month
        w = days / days.sum()
        mean_da = (da * w).sum("time") / w.where(da.notnull()).sum("time")
        return mean_da.compute()
    else:
        return da.compute()


def get_lat_lon_2d(ds, model):
    info = model_info[model]
    lat = ds[info["lat_coord"]].values
    lon = ds[info["lon_coord"]].values
    if lat.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon, lat)
    else:
        lat2d = lat
        lon2d = lon
    return lat2d, lon2d


def sph_xyz(lat, lon):
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)
    x = np.cos(lat) * np.cos(lon)
    y = np.cos(lat) * np.sin(lon)
    z = np.sin(lat)
    return np.column_stack([x, y, z])


def rho_insitu(T, S, P):
    """UNESCO 1983 (EOS-80) in-situ density (kg m-3). T degC, S PSU, P dbar."""
    T = np.asarray(T)
    S = np.asarray(S)
    P = np.asarray(P)

    rho_w = (999.842594 + 6.793952e-2 * T - 9.095290e-3 * T**2
             + 1.001685e-4 * T**3 - 1.120083e-6 * T**4 + 6.536332e-9 * T**5)
    A = (8.24493e-1 - 4.0899e-3 * T + 7.6438e-5 * T**2
         - 8.2467e-7 * T**3 + 5.3875e-9 * T**4)
    B = -5.72466e-3 + 1.0227e-4 * T - 1.6546e-6 * T**2
    C = 4.8314e-4
    rho_0 = rho_w + A * S + B * S**1.5 + C * S**2

    K_w = (19652.21 + 148.4206 * T - 2.327105 * T**2
           + 1.360477e-2 * T**3 - 5.155288e-5 * T**4)
    K_0 = (K_w
           + (54.6746 - 0.603459 * T + 1.09987e-2 * T**2 - 6.1670e-5 * T**3) * S
           + (7.944e-2 + 1.6483e-2 * T - 5.3009e-4 * T**2) * S**1.5)

    A_K = 3.239908 + 1.43713e-3 * T + 1.16092e-4 * T**2 - 5.77905e-7 * T**3
    B_K = 2.2838e-3 - 1.0981e-5 * T - 1.6078e-6 * T**2
    C_K = 1.91075e-4
    D_K = 8.50935e-5 - 6.12293e-6 * T + 5.2787e-8 * T**2
    E_K = -9.9348e-7 + 2.0816e-8 * T + 9.1697e-10 * T**2

    K = (K_0
         + P * (A_K + B_K * S + C_K * S**1.5)
         + P**2 * (D_K + E_K * S))
    return rho_0 / (1 - P / K)


def parent_basin(code, lat, lon):
    if code in [1, 2, 3, 4, 5]:
        return code
    # WOA13 marginal-sea groupings
    if code in [6, 7, 8, 11, 12, 13]:
        return 1  # Atlantic
    if code in [9, 10, 18, 19, 20]:
        return 3  # Indian
    if code in [14, 15, 16, 17]:
        return 2  # Pacific
    # Geographic fallback
    if (-10 <= lon <= 45) and (30 <= lat <= 70):
        return 1
    if (30 <= lon <= 100) and (-10 <= lat <= 30):
        return 3
    if ((100 <= lon <= 180) or (-180 <= lon <= -100)) and (0 <= lat <= 66):
        return 2
    if (-100 <= lon <= -50) and (10 <= lat <= 70):
        return 1
    if lat > 66:
        return 4
    if lat < -60:
        return 5
    if -100 <= lon <= 20:
        return 1
    if 20 < lon <= 120:
        return 3
    return 2


def regrid_to_woa(lat2d, lon2d, field):
    """Nearest/linear regrid of a curvilinear field to the WOA 1-degree grid."""
    mask = np.isfinite(field)
    if mask.sum() == 0:
        return np.full((len(lat_woa), len(lon_woa)), np.nan)

    lon360 = lon2d % 360
    lat_s = lat2d[mask]
    lon_s = lon360[mask]
    val_s = field[mask]

    lon_s_wrap = np.concatenate([lon_s, lon_s + 360, lon_s - 360])
    lat_s_wrap = np.concatenate([lat_s, lat_s, lat_s])
    val_s_wrap = np.concatenate([val_s, val_s, val_s])

    lon_t = lon_woa % 360
    lat_t = lat_woa
    lon_t_grid, lat_t_grid = np.meshgrid(lon_t, lat_t)

    out = griddata((lon_s_wrap, lat_s_wrap), val_s_wrap,
                   (lon_t_grid, lat_t_grid), method="linear")
    if np.isnan(out).any():
        out_near = griddata((lon_s_wrap, lat_s_wrap), val_s_wrap,
                            (lon_t_grid, lat_t_grid), method="nearest")
        out = np.where(np.isnan(out), out_near, out)
    return out


# ---------------------------------------------------------------------------
# WOA23 grid, masks and depth reference
# ---------------------------------------------------------------------------
ds_ref = xr.open_dataset(
    os.path.join(woa_root, "nitrate/all_time_periods/woa23_all_n00_01.nc"),
    decode_times=False,
)
depth_ref = ds_ref["depth"].values
lat_woa = ds_ref["lat"].values
lon_woa = ds_ref["lon"].values
nlat_woa = len(lat_woa)
nlon_woa = len(lon_woa)
lat_origin = float(lat_woa[0])
lon_origin = float(lon_woa[0])

# --- Land-sea mask (bottom standard-level per column) ----------------------
arr_landsea = read_msk(os.path.join(woa_masks_root, "landsea_01.msk"))
bottom_depth_2d = np.full((nlat_woa, nlon_woa), -1.0)
for vals in arr_landsea:
    la, lo, b = vals[0], vals[1], vals[2]
    if not np.isfinite(b):
        continue
    i = int(np.round(la - lat_origin))
    j = int(np.round(lo - lon_origin))
    if not (0 <= i < nlat_woa and 0 <= j < nlon_woa):
        continue
    if b <= 1:
        bottom_depth_2d[i, j] = -1.0
    elif b > 103:
        bottom_depth_2d[i, j] = 1e6
    else:
        idx = int(b) - 1
        bottom_depth_2d[i, j] = depth_ref[idx] if idx < len(depth_ref) else 1e6

# --- Basin mask ------------------------------------------------------------
arr_basin = read_msk(os.path.join(woa_masks_root, "basinmask_01.msk"))
basin_2d = np.zeros((nlat_woa, nlon_woa), dtype=int)
for vals in arr_basin:
    la, lo = vals[0], vals[1]
    candidates = [v for v in vals[2:] if np.isfinite(v) and v > 0]
    if not candidates:
        continue
    b = int(candidates[0])
    i = int(np.round(la - lat_origin))
    j = int(np.round(lo - lon_origin))
    if 0 <= i < nlat_woa and 0 <= j < nlon_woa:
        basin_2d[i, j] = b

basin_parent = np.zeros_like(basin_2d)
for i in range(nlat_woa):
    for j in range(nlon_woa):
        b = basin_2d[i, j]
        if b > 0:
            basin_parent[i, j] = parent_basin(b, lat_woa[i], lon_woa[j])

# --- WOA cell areas --------------------------------------------------------
R = 6371000.0
area_lat = R**2 * (np.pi / 180.0) * (
    np.sin(np.deg2rad(lat_woa + 0.5)) - np.sin(np.deg2rad(lat_woa - 0.5))
)
area2d = np.repeat(area_lat[:, None], nlon_woa, axis=1)

# ---------------------------------------------------------------------------
# WOA23 diagnostics
# ---------------------------------------------------------------------------
woa_o2, woa_ds, woa_o2_period = load_woa_monthly("oxygen", "1971-2000")
depth_bnds = get_depth_bounds(woa_ds)
layer_thick = depth_bnds[:, 1] - depth_bnds[:, 0]
vol_woa = area2d[None, :, :] * layer_thick[:, None, None]

valid_woa = (
    woa_o2.notnull().values
    & (woa_o2.depth.values[:, None, None] <= bottom_depth_2d[None, :, :])
)
o2_woa = woa_o2.values

V20_woa = float(vol_woa[valid_woa & (o2_woa < 20)].sum()) / 1e15
V60_woa = float(vol_woa[valid_woa & (o2_woa < 60)].sum()) / 1e15

thick_woa = np.where(valid_woa & (o2_woa < 60),
                     layer_thick[:, None, None], 0.0).sum(axis=0)
max_thick_woa = float(np.nanmax(thick_woa))
idx = np.unravel_index(np.nanargmax(thick_woa), thick_woa.shape)
lat_max_woa = float(lat_woa[idx[0]])
lon_max_woa = float(lon_woa[idx[1]])

woa_aou, woa_aou_ds, _ = load_woa_monthly("apparent_oxygen_utilization", "1971-2000")
woa_po2, woa_po2_ds, _ = load_woa_monthly("percent_oxygen_saturation", "1971-2000")

valid_aou = (
    woa_aou.notnull().values
    & (woa_aou.depth.values[:, None, None] <= bottom_depth_2d[None, :, :])
)
valid_po2 = (
    woa_po2.notnull().values
    & (woa_po2.depth.values[:, None, None] <= bottom_depth_2d[None, :, :])
)
omz_mask = valid_woa & (o2_woa < 60)

m_aou = valid_aou & omz_mask
mean_aou = float(np.nansum(woa_aou.values[m_aou] * vol_woa[m_aou])
                 / np.nansum(vol_woa[m_aou]))
m_po2 = valid_po2 & omz_mask
mean_po2 = float(np.nansum(woa_po2.values[m_po2] * vol_woa[m_po2])
                 / np.nansum(vol_woa[m_po2]))

basin_vols_woa = {}
for code in [1, 2, 3]:
    mask = valid_woa & (o2_woa < 60) & (basin_parent[None, :, :] == code)
    basin_vols_woa[code] = float(vol_woa[mask].sum()) / 1e15

results = {
    "woa": {
        "period": woa_o2_period,
        "V20": V20_woa,
        "V60": V60_woa,
        "max_thick": max_thick_woa,
        "lat_max": lat_max_woa,
        "lon_max": lon_max_woa,
        "basin_vols": basin_vols_woa,
        "mean_aou": mean_aou,
        "mean_po2": mean_po2,
        "thick_field": thick_woa,
    },
    "density_assumption": (
        "In-situ density from UNESCO 1983 EOS-80 using period-mean thetao and so; "
        "pressure approximated by depth in metres (1 m = 1 dbar)."
    ),
}

# ---------------------------------------------------------------------------
# CMIP6 diagnostics
# ---------------------------------------------------------------------------
model_results = {}


def compute_model_diagnostics(model, exp, start, end):
    root = cmip6_hist_root if exp == "historical" else cmip6_ssp_root

    o2_native = get_period_mean(root, model, "o2", exp, start, end, allow_fallback=False)
    thetao = get_period_mean(root, model, "thetao", exp, start, end, allow_fallback=False)
    so = get_period_mean(root, model, "so", exp, start, end, allow_fallback=False)
    volcello = get_period_mean(root, model, "volcello", exp, start, end, allow_fallback=True)
    areacello = get_period_mean(root, model, "areacello", exp, start, end, allow_fallback=True)
    # deptho loaded for completeness (not used directly)
    _ = get_period_mean(root, model, "deptho", exp, start, end, allow_fallback=True)

    o2_native, thetao, so, volcello, areacello = xr.align(
        o2_native, thetao, so, volcello, areacello, join="inner"
    )

    depth_m = o2_native.lev.values * model_info[model]["depth_scale"]
    p = depth_m[:, None, None] * np.ones_like(o2_native.values)
    rho = rho_insitu(thetao.values, so.values, p)

    o2_umol = o2_native.values * 1e6 / rho
    valid = (
        np.isfinite(o2_native.values)
        & (o2_native.values < 1e19)
        & np.isfinite(thetao.values)
        & np.isfinite(so.values)
    )

    vol = volcello.values
    if vol.shape != o2_native.shape:
        vol = np.broadcast_to(vol, o2_native.shape)

    V20 = float(np.nansum(vol[valid & (o2_umol < 20)])) / 1e15
    V60 = float(np.nansum(vol[valid & (o2_umol < 60)])) / 1e15
    inventory = float(np.nansum(o2_native.values[valid] * vol[valid]))

    area = areacello.values
    area3 = np.broadcast_to(area[None, :, :], vol.shape)
    thick = np.where(area3 > 0, vol / area3, 0.0)
    thick_60 = np.where(valid & (o2_umol < 60), thick, 0.0).sum(axis=0)
    max_thick = float(np.nanmax(thick_60))

    lat2d, lon2d = get_lat_lon_2d(o2_native, model)
    idx = np.unravel_index(np.nanargmax(thick_60), thick_60.shape)
    lat_max = float(lat2d[idx])
    lon_max = float(lon2d[idx])

    # Basin mapping
    lat_woa_grid, lon_woa_grid = np.meshgrid(lat_woa, lon_woa, indexing="ij")
    mask_ocean = basin_parent > 0
    pts = sph_xyz(lat_woa_grid[mask_ocean], lon_woa_grid[mask_ocean])
    vals = basin_parent[mask_ocean]
    tree = cKDTree(pts)
    q = sph_xyz(lat2d.ravel(), lon2d.ravel())
    _, nn_idx = tree.query(q, k=1)
    model_basin = vals[nn_idx].reshape(lat2d.shape)

    basin_vols = {}
    for code in [1, 2, 3]:
        mask = valid & (o2_umol < 60) & (model_basin[None, :, :] == code)
        basin_vols[code] = float(np.nansum(vol[mask])) / 1e15

    return {
        "V20": V20,
        "V60": V60,
        "inventory": inventory,
        "max_thick": max_thick,
        "lat_max": lat_max,
        "lon_max": lon_max,
        "basin_vols": basin_vols,
        "thick_field": thick_60,
        "lat2d": lat2d,
        "lon2d": lon2d,
    }


# Historical period
hist_start, hist_end = "1971-01-01", "2000-12-31"
for model in ["GFDL-ESM4", "MPI-ESM1-2-LR"]:
    model_results[(model, "historical")] = compute_model_diagnostics(
        model, "historical", hist_start, hist_end
    )

# ssp245 period
ssp_start, ssp_end = "2071-01-01", "2100-12-31"
for model in ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]:
    model_results[(model, "ssp245")] = compute_model_diagnostics(
        model, "ssp245", ssp_start, ssp_end
    )

results["models"] = model_results

# Historical-to-ssp245 changes
changes = {}
for model in ["GFDL-ESM4", "MPI-ESM1-2-LR"]:
    h = model_results[(model, "historical")]
    s = model_results[(model, "ssp245")]
    changes[model] = {
        "V20_pct": 100.0 * (s["V20"] - h["V20"]) / h["V20"],
        "V60_pct": 100.0 * (s["V60"] - h["V60"]) / h["V60"],
    }
results["changes"] = changes

# Model vs WOA comparison
compare = {}
for model in ["GFDL-ESM4", "MPI-ESM1-2-LR"]:
    h = model_results[(model, "historical")]
    compare[model] = {
        "V20_diff": h["V20"] - V20_woa,
        "V60_diff": h["V60"] - V60_woa,
        "V20_abs_diff": abs(h["V20"] - V20_woa),
        "V60_abs_diff": abs(h["V60"] - V60_woa),
    }
results["compare"] = compare
results["best_model_V20"] = min(compare, key=lambda m: compare[m]["V20_abs_diff"])
results["best_model_V60"] = min(compare, key=lambda m: compare[m]["V60_abs_diff"])
