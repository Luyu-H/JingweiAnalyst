import os
import glob
import re
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
ROOT_CMIP = "/data1/lorraine-20419/benchmark_cmip6"
ROOT_WOA = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
MASKS = "/data1/lorraine-20419/ocean_data/WOA23/masks"

# -----------------------------------------------------------------------------
# Small helpers
# -----------------------------------------------------------------------------
def cmip_files(model, exp, var):
    """Return CMIP6 files for one model/experiment/variable, choosing one variant."""
    d = os.path.join(ROOT_CMIP, exp, model, var)
    files = sorted(glob.glob(os.path.join(d, "*.nc")))
    files = [f for f in files if f"_{exp}_" in os.path.basename(f)]
    if not files:
        raise FileNotFoundError(f"No files for {model} {exp} {var} in {d}")
    variants = {}
    for f in files:
        bn = os.path.basename(f)
        m = re.search(r"_(r\d+i\d+p\d+f\d+)_", bn)
        v = m.group(1) if m else "unknown"
        variants.setdefault(v, []).append(f)
    if len(variants) > 1:
        files = max(variants.values(), key=len)
    return files


def select_time(ds, start, end):
    """Select a time slice, robust to cftime calendars."""
    if "time" not in ds.dims:
        return ds
    try:
        return ds.sel(time=slice(start, end))
    except Exception:
        sy, sm, sd = map(int, start.split("-"))
        ey, em, ed = map(int, end.split("-"))
        times = ds.time.values
        mask = []
        for t in times:
            yy, mm, dd = t.year, t.month, t.day
            ok = (yy > sy or (yy == sy and mm >= sm)) and (
                yy < ey or (yy == ey and mm <= em)
            )
            mask.append(ok)
        return ds.isel(time=np.array(mask))


def get_model_var(model, exp, var, start, end, fallback_historical=False):
    """Load a CMIP6 variable, time-averaged over the requested period."""
    try:
        files = cmip_files(model, exp, var)
        used_exp = exp
    except FileNotFoundError:
        if fallback_historical and exp != "historical":
            files = cmip_files(model, "historical", var)
            used_exp = "historical"
        else:
            raise

    ds = xr.open_mfdataset(
        files,
        combine="by_coords",
        use_cftime=True,
        decode_times=True,
    )

    vname = var
    if vname not in ds:
        candidates = [v for v in ds.data_vars if v.lower() == var.lower()]
        if not candidates:
            candidates = [v for v in ds.data_vars if var.lower() in v.lower()]
        if not candidates:
            raise KeyError(f"{var!r} not in {model} {exp}; vars={list(ds.data_vars)}")
        vname = candidates[0]

    da = ds[vname]

    # If a fixed field had to be taken from the historical run, do not try
    # to slice it with the scenario-period dates.  For a time-varying fallback
    # field, use the full available historical mean.
    if used_exp != exp and "time" in da.dims:
        da = da.mean("time")
        return da.load()

    ds = select_time(ds, start, end)
    da = ds[vname]
    if "time" in da.dims:
        da = da.mean("time")
    return da.load()


def read_numeric_table(path, min_cols=1):
    """Read whitespace-separated numeric rows from a .msk file robustly."""
    rows = []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#") or s.startswith("!"):
                continue
            s = s.replace(",", " ")
            parts = s.split()
            vals = []
            for p in parts:
                p = p.replace("D", "E").replace("d", "e")
                try:
                    vals.append(float(p))
                except ValueError:
                    continue
            if len(vals) >= min_cols:
                rows.append(vals)

    if not rows:
        raise ValueError(f"No numeric rows found in {path}")

    maxlen = max(len(r) for r in rows)
    arr = np.full((len(rows), maxlen), np.nan, dtype=float)
    for i, r in enumerate(rows):
        arr[i, :len(r)] = r
    return arr


def read_basin_mask(path, nlev):
    """Build a [nlev, 180, 360] basin array from basinmask_01.msk."""
    arr = read_numeric_table(path, min_cols=4)
    if arr.ndim != 2 or arr.shape[1] < 4:
        raise ValueError(f"basin mask {path} does not have at least 4 numeric columns")
    arr = arr[np.all(np.isfinite(arr[:, :4]), axis=1)]
    arr = arr[
        (arr[:, 0] >= 0)
        & (arr[:, 1] >= -90.0)
        & (arr[:, 1] <= 90.0)
        & (arr[:, 2] >= -180.0)
        & (arr[:, 2] <= 180.0)
    ]

    level = arr[:, 0].astype(int)
    lat = arr[:, 1]
    lon = arr[:, 2]
    basin = arr[:, 3].astype(int)

    if level.min() >= 1:
        lev_idx = level - 1
    else:
        lev_idx = level

    lat_idx = np.round((lat - (-89.5)) / 1.0).astype(int)
    lon_idx = np.round((lon - (-179.5)) / 1.0).astype(int)
    lon_idx = np.mod(lon_idx, 360)

    out = np.full((nlev, 180, 360), np.nan)
    ok = (
        (lev_idx >= 0)
        & (lev_idx < nlev)
        & (lat_idx >= 0)
        & (lat_idx < 180)
        & (lon_idx >= 0)
        & (lon_idx < 360)
    )
    out[lev_idx[ok], lat_idx[ok], lon_idx[ok]] = basin[ok]
    return out


def read_landsea_mask(path, nlev):
    """Build a [nlev, 180, 360] ocean boolean mask from landsea_01.msk."""
    arr = read_numeric_table(path, min_cols=3)
    if arr.ndim != 2 or arr.shape[1] < 3:
        raise ValueError(f"landsea mask {path} does not have at least 3 numeric columns")
    arr = arr[np.all(np.isfinite(arr[:, :3]), axis=1)]
    arr = arr[
        (arr[:, 0] >= -90.0)
        & (arr[:, 0] <= 90.0)
        & (arr[:, 1] >= -180.0)
        & (arr[:, 1] <= 180.0)
        & (arr[:, 2] >= 1)
    ]

    lat = arr[:, 0]
    lon = arr[:, 1]
    bottom = arr[:, 2].astype(int)

    lat_idx = np.round((lat - (-89.5)) / 1.0).astype(int)
    lon_idx = np.round((lon - (-179.5)) / 1.0).astype(int)
    lon_idx = np.mod(lon_idx, 360)

    out = np.zeros((nlev, 180, 360), dtype=bool)
    deep = np.where(bottom > 1, bottom - 1, -1)
    deep = np.minimum(deep, nlev - 1)
    for lev in range(nlev):
        ok = (
            (lat_idx >= 0)
            & (lat_idx < 180)
            & (lon_idx >= 0)
            & (lon_idx < 360)
            & (deep >= lev)
        )
        out[lev, lat_idx[ok], lon_idx[ok]] = True
    return out


def fold_basin_codes(b):
    """Fold WOA13 marginal-sea codes into their parent ocean."""
    b = b.copy()
    b[np.isin(b, [6, 7, 8])] = 1
    b[np.isin(b, [9, 10])] = 3
    b[np.isin(b, [11])] = 2
    b[np.isin(b, [12])] = 1
    return b


def get_woa_thickness(ds):
    """Layer thickness in metres for WOA standard levels."""
    if "depth_bnds" in ds:
        b = ds["depth_bnds"].values
        return b[:, 1] - b[:, 0]
    d = ds["depth"].values.astype(float)
    if len(d) == 1:
        return np.array([2.0 * d[0]])
    edges = np.empty(len(d) + 1)
    edges[0] = 0.0
    edges[1:-1] = 0.5 * (d[:-1] + d[1:])
    edges[-1] = d[-1] + 0.5 * (d[-1] - d[-2])
    return np.diff(edges)


def load_woa(subdir, tag, prefix, field="an"):
    """Load the annual objectively-analysed WOA23 field."""
    path = os.path.join(ROOT_WOA, subdir, f"woa23_{tag}_{prefix}00_01.nc")
    ds = xr.open_dataset(path, decode_times=False)

    vname = None
    for name in ds.data_vars:
        if name.lower() == f"{prefix.lower()}_{field}".lower():
            vname = name
            break
    if vname is None:
        for name in ds.data_vars:
            if name.lower().startswith(prefix.lower()) and field in name.lower():
                vname = name
                break
    if vname is None:
        for name in ds.data_vars:
            if "depth" in ds[name].dims:
                vname = name
                break
    if vname is None:
        raise KeyError(f"Could not find WOA field {prefix}_{field} in {path}")

    da = ds[vname]
    if "time" in da.dims:
        da = da.isel(time=0)
    da = da.squeeze()
    return da.load(), ds


def da_woa_to_array(da):
    """Return WOA array ordered as depth, lat, lon."""
    dims = [d for d in ["depth", "lat", "lon"] if d in da.dims]
    if len(dims) == 3:
        da = da.transpose(*dims)
    return da.values.astype(float)


def seawater_density(S, T, depth_m):
    """UNESCO/EOS-80 in-situ density, approximately, in kg m-3."""
    S = np.asarray(S, dtype=float)
    T = np.asarray(T, dtype=float)
    P = np.asarray(depth_m, dtype=float) / 10.0

    Sf = np.where(np.isfinite(S), S, 35.0)
    Tf = np.where(np.isfinite(T), T, 0.0)

    if P.ndim == 1 and T.ndim > 1:
        P = P.reshape((-1,) + (1,) * (T.ndim - 1))

    rho_w = (
        999.842594
        + 6.793952e-2 * Tf
        - 9.095290e-3 * Tf**2
        + 1.001685e-4 * Tf**3
        - 1.120083e-6 * Tf**4
        + 6.536332e-9 * Tf**5
    )
    A = (
        8.24493e-1
        - 4.0899e-3 * Tf
        + 7.6438e-5 * Tf**2
        - 8.2467e-7 * Tf**3
        + 5.3875e-9 * Tf**4
    )
    B = -5.72466e-3 + 1.0227e-4 * Tf - 1.6546e-6 * Tf**2
    C = 4.8314e-4
    rho0 = rho_w + A * Sf + B * Sf**1.5 + C * Sf**2

    Kw = (
        19652.21
        + 148.4206 * Tf
        - 2.327105 * Tf**2
        + 1.360477e-2 * Tf**3
        - 5.155288e-5 * Tf**4
    )
    K0 = (
        Kw
        + Sf
        * (
            54.6746
            - 0.603459 * Tf
            + 1.09987e-2 * Tf**2
            - 6.1670e-5 * Tf**3
        )
        + Sf**1.5
        * (7.944e-2 + 1.6483e-2 * Tf - 5.3009e-4 * Tf**2)
    )
    A2 = (
        3.239908
        + 1.43713e-3 * Tf
        + 1.16092e-4 * Tf**2
        - 5.77905e-7 * Tf**3
    )
    B2 = 2.209460e-2 - 1.262e-4 * Tf - 2.1371e-6 * Tf**2
    K = K0 + A2 * P + B2 * P**2
    return rho0 / (1.0 - P / K)


def map_basin_to_model(basin3d, woa_lat, woa_lon, woa_depth,
                       model_lat2d, model_lon2d, model_depth):
    """Nearest-neighbour mapping of the WOA basin mask onto a model grid."""
    lon180 = ((model_lon2d + 180.0) % 360.0) - 180.0
    li = np.round((model_lat2d - woa_lat[0]) / (woa_lat[1] - woa_lat[0])).astype(int)
    oi = np.round((lon180 - woa_lon[0]) / (woa_lon[1] - woa_lon[0])).astype(int)
    li = np.clip(li, 0, len(woa_lat) - 1)
    oi = np.mod(oi, len(woa_lon))
    di = np.abs(model_depth[:, None] - woa_depth[None, :]).argmin(axis=1)
    return basin3d[di[:, None, None], li[None, :, :], oi[None, :, :]]


def omz_metrics(o2_umolkg, volcello, thickness, lat2d, lon2d,
                basin3d=None, basin_map=None, thresholds=(20, 60)):
    """Volume, column thickness, max thickness, and basin volumes for OMZ thresholds."""
    valid = (
        np.isfinite(o2_umolkg)
        & np.isfinite(volcello)
        & (volcello > 0)
        & np.isfinite(thickness)
    )
    out = {}
    for thr in thresholds:
        mask = valid & (o2_umolkg < thr)
        out[f"vol_{thr}"] = float(np.nansum(volcello[mask]))
        thick_map = np.nansum(np.where(mask, thickness, 0.0), axis=0)
        out[f"thick_map_{thr}"] = thick_map
        if np.any(np.isfinite(thick_map)) and np.nanmax(thick_map) > 0:
            idx = np.nanargmax(thick_map)
            out[f"max_thick_{thr}"] = float(thick_map.flat[idx])
            out[f"max_lat_{thr}"] = float(lat2d.flat[idx])
            out[f"max_lon_{thr}"] = float(lon2d.flat[idx])
        else:
            out[f"max_thick_{thr}"] = np.nan
            out[f"max_lat_{thr}"] = np.nan
            out[f"max_lon_{thr}"] = np.nan
        if basin3d is not None and basin_map is not None:
            for code, bname in basin_map.items():
                bmask = mask & (basin3d == code)
                out[f"vol_{thr}_{bname}"] = float(np.nansum(volcello[bmask]))
    return out


# -----------------------------------------------------------------------------
# WOA23 observations
# -----------------------------------------------------------------------------
woa_o2_da, woa_ds = load_woa(
    "oxygen/1971-2000_climate_normal", "decav71A0", "o"
)
woa_aou_da, _ = load_woa(
    "apparent_oxygen_utilization/1971-2000_climate_normal", "decav71A0", "A"
)
woa_sat_da, _ = load_woa(
    "percent_oxygen_saturation/1971-2000_climate_normal", "decav71A0", "O"
)

o2_woa = da_woa_to_array(woa_o2_da)
aou_woa = da_woa_to_array(woa_aou_da)
sat_woa = da_woa_to_array(woa_sat_da)

woa_lat = woa_ds["lat"].values.astype(float)
woa_lon = woa_ds["lon"].values.astype(float)
woa_depth = woa_ds["depth"].values.astype(float)
woa_thick = get_woa_thickness(woa_ds)

nlev_woa = len(woa_depth)
basin_woa = read_basin_mask(os.path.join(MASKS, "basinmask_01.msk"), nlev_woa)
basin_woa = fold_basin_codes(basin_woa)
ocean_woa = read_landsea_mask(os.path.join(MASKS, "landsea_01.msk"), nlev_woa)

# WOA cell areas and volumes
Earth_R = 6371000.0
dlat = np.deg2rad(1.0)
dlon = np.deg2rad(1.0)
area_1d = Earth_R**2 * np.cos(np.deg2rad(woa_lat)) * dlat * dlon
area2d_woa = np.repeat(area_1d[:, None], len(woa_lon), axis=1)

vol_woa = area2d_woa[None, :, :] * woa_thick[:, None, None]
valid_woa = ocean_woa & np.isfinite(o2_woa)
vol_woa = np.where(valid_woa, vol_woa, 0.0)

thick3d_woa = np.broadcast_to(woa_thick[:, None, None], o2_woa.shape).astype(float)
thick3d_woa = np.where(valid_woa, thick3d_woa, np.nan)

lon2d_woa, lat2d_woa = np.meshgrid(woa_lon, woa_lat)

BASIN_MAP = {
    1: "Atlantic",
    2: "Pacific",
    3: "Indian",
    4: "Southern",
    5: "Arctic",
}

metrics_woa = omz_metrics(
    o2_woa, vol_woa, thick3d_woa,
    lat2d_woa, lon2d_woa,
    basin_woa, BASIN_MAP
)

mask60_woa = valid_woa & (o2_woa < 60.0)
w60 = np.where(mask60_woa, vol_woa, 0.0)
aou_mean = float(np.nansum(aou_woa * w60) / np.nansum(w60))
sat_mean = float(np.nansum(sat_woa * w60) / np.nansum(w60))


# -----------------------------------------------------------------------------
# CMIP6 processing
# -----------------------------------------------------------------------------
def process_model_combo(model, exp, start, end):
    """Load one model/experiment, convert O2 to umol/kg, and compute OMZ metrics."""
    o2_da = get_model_var(model, exp, "o2", start, end)

    lev = o2_da["lev"].values.astype(float)
    lev_units = o2_da["lev"].attrs.get("units", "").lower()
    if "centimeter" in lev_units or lev_units == "cm":
        lev = lev * 0.01

    if model == "GFDL-ESM4":
        lat1 = o2_da["lat"].values.astype(float)
        lon1 = o2_da["lon"].values.astype(float)
        lon2d, lat2d = np.meshgrid(lon1, lat1)
    else:
        latname = "lat" if "lat" in o2_da.coords else "latitude"
        lonname = "lon" if "lon" in o2_da.coords else "longitude"
        lat2d = o2_da[latname].values.astype(float)
        lon2d = o2_da[lonname].values.astype(float)

    thetao_da = get_model_var(model, exp, "thetao", start, end)
    so_da = get_model_var(model, exp, "so", start, end)
    vol_da = get_model_var(model, exp, "volcello", start, end, fallback_historical=True)
    area_da = get_model_var(model, exp, "areacello", start, end, fallback_historical=True)

    try:
        deptho_da = get_model_var(model, exp, "deptho", start, end, fallback_historical=True)
        deptho = deptho_da.values.astype(float)
    except Exception:
        deptho = None

    def lev_yx(da):
        dims = list(da.dims)
        depth_dim = "lev" if "lev" in dims else ("depth" if "depth" in dims else dims[0])
        other = [d for d in dims if d != depth_dim]
        if len(other) > 2:
            other = other[-2:]
        return da.transpose(depth_dim, other[0], other[1]).values.astype(float)

    o2 = lev_yx(o2_da)
    thetao = lev_yx(thetao_da)
    so = lev_yx(so_da)
    vol = lev_yx(vol_da)
    area = area_da.values.astype(float)
    if area.shape != lat2d.shape:
        area = area.T

    # Normalise practical salinity if the file declares 0.001-style units.
    if np.nanmax(so) > 1000.0:
        so = so / 1000.0
    if np.nanmax(so) < 1.0:
        so = so * 1000.0

    rho = seawater_density(so, thetao, lev)
    o2_umolkg = o2 * 1e6 / rho
    thickness = vol / area[None, :, :]

    basin_model = map_basin_to_model(
        basin_woa, woa_lat, woa_lon, woa_depth,
        lat2d, lon2d, lev
    )

    met = omz_metrics(
        o2_umolkg, vol, thickness,
        lat2d, lon2d,
        basin_model, BASIN_MAP
    )
    met["inventory_mol"] = float(np.nansum(o2 * vol))
    met["o2_umolkg"] = o2_umolkg
    met["volcello"] = vol
    met["thickness"] = thickness
    met["lat2d"] = lat2d
    met["lon2d"] = lon2d
    met["basin"] = basin_model
    met["lev"] = lev
    met["deptho"] = deptho
    return met


combos = [
    ("GFDL-ESM4_historical", "GFDL-ESM4", "historical", "1971-01-01", "2000-12-31"),
    ("MPI-ESM1-2-LR_historical", "MPI-ESM1-2-LR", "historical", "1971-01-01", "2000-12-31"),
    ("CESM2_ssp245", "CESM2", "ssp245", "2071-01-01", "2100-12-31"),
    ("GFDL-ESM4_ssp245", "GFDL-ESM4", "ssp245", "2071-01-01", "2100-12-31"),
    ("MPI-ESM1-2-LR_ssp245", "MPI-ESM1-2-LR", "ssp245", "2071-01-01", "2100-12-31"),
]

model_results = {}
for name, model, exp, start, end in combos:
    print(f"Processing {name} ...")
    model_results[name] = process_model_combo(model, exp, start, end)


# -----------------------------------------------------------------------------
# Scenario changes and model-observation comparison
# -----------------------------------------------------------------------------
pct_change = {}
for model in ["GFDL-ESM4", "MPI-ESM1-2-LR"]:
    h = model_results[f"{model}_historical"]
    s = model_results[f"{model}_ssp245"]
    pct_change[model] = {}
    for thr in [20, 60]:
        pct_change[model][thr] = (
            (s[f"vol_{thr}"] - h[f"vol_{thr}"]) / h[f"vol_{thr}"] * 100.0
        )

rel_err = {}
for model in ["GFDL-ESM4", "MPI-ESM1-2-LR"]:
    h = model_results[f"{model}_historical"]
    rel_err[model] = {}
    for thr in [20, 60]:
        rel_err[model][thr] = (
            (h[f"vol_{thr}"] - metrics_woa[f"vol_{thr}"])
            / metrics_woa[f"vol_{thr}"]
            * 100.0
        )


# -----------------------------------------------------------------------------
# Tabular summary
# -----------------------------------------------------------------------------
summary_rows = []

def add_summary(source, experiment, met):
    for thr in [20, 60]:
        row = {
            "source": source,
            "experiment": experiment,
            "threshold_umol_kg": thr,
            "volume_1e15_m3": met[f"vol_{thr}"] / 1e15,
            "max_thickness_m": met[f"max_thick_{thr}"],
            "max_lat": met[f"max_lat_{thr}"],
            "max_lon": met[f"max_lon_{thr}"],
        }
        if "inventory_mol" in met:
            row["inventory_mol"] = met["inventory_mol"]
        for b in ["Pacific", "Indian", "Atlantic"]:
            row[f"vol_{b}_1e15_m3"] = met.get(f"vol_{thr}_{b}", np.nan) / 1e15
        summary_rows.append(row)

add_summary("WOA23", "1971-2000", metrics_woa)
add_summary("GFDL-ESM4", "historical", model_results["GFDL-ESM4_historical"])
add_summary("MPI-ESM1-2-LR", "historical", model_results["MPI-ESM1-2-LR_historical"])
add_summary("CESM2", "ssp245", model_results["CESM2_ssp245"])
add_summary("GFDL-ESM4", "ssp245", model_results["GFDL-ESM4_ssp245"])
add_summary("MPI-ESM1-2-LR", "ssp245", model_results["MPI-ESM1-2-LR_ssp245"])

summary_df = pd.DataFrame(summary_rows)

print("\n--- OMZ volume, max thickness, inventory ---")
for name, met in [("WOA23", metrics_woa)] + list(model_results.items()):
    print(name)
    for thr in [20, 60]:
        print(
            f"  <{thr}: {met[f'vol_{thr}']/1e15:.3g} x 1e15 m3; "
            f"max thickness {met[f'max_thick_{thr}']:.1f} m at "
            f"({met[f'max_lat_{thr}']:.2f}, {met[f'max_lon_{thr}']:.2f})"
        )
    if "inventory_mol" in met:
        print(f"  inventory: {met['inventory_mol']:.3g} mol")

print("\n--- percentage change historical -> ssp245 ---")
for model in pct_change:
    print(model, {thr: f"{pct_change[model][thr]:.2f}%" for thr in [20, 60]})

print("\n--- relative error in historical OMZ volume vs WOA23 ---")
for model in rel_err:
    print(model, {thr: f"{rel_err[model][thr]:.2f}%" for thr in [20, 60]})
for thr in [20, 60]:
    better = min(rel_err, key=lambda m: abs(rel_err[m][thr]))
    print(f"At < {thr} umol/kg, closer historical model: {better}")

print("\n--- basin volumes at 60 umol/kg (1e15 m3) ---")
for name, met in [("WOA23", metrics_woa)] + list(model_results.items()):
    vals = {
        b: met.get(f"vol_60_{b}", np.nan) / 1e15
        for b in ["Pacific", "Indian", "Atlantic"]
    }
    print(name, {k: f"{v:.3g}" for k, v in vals.items()})

print(f"\nWOA OMZ<60 volume-weighted AOU: {aou_mean:.1f} umol/kg")
print(f"WOA OMZ<60 volume-weighted O2 saturation: {sat_mean:.1f}%")


# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------
hist_maps = [
    model_results["GFDL-ESM4_historical"]["thick_map_60"],
    model_results["MPI-ESM1-2-LR_historical"]["thick_map_60"],
    metrics_woa["thick_map_60"],
]
vmax_hist = max(np.nanmax(m) for m in hist_maps)

fig1, axes1 = plt.subplots(1, 3, figsize=(18, 5))
panels1 = [
    (axes1[0], model_results["GFDL-ESM4_historical"], "GFDL-ESM4 historical"),
    (axes1[1], model_results["MPI-ESM1-2-LR_historical"], "MPI-ESM1-2-LR historical"),
    (axes1[2], None, "WOA23 1971-2000"),
]
for ax, met, title in panels1:
    if met is None:
        lon = lon2d_woa
        lat = lat2d_woa
        data = metrics_woa["thick_map_60"]
    else:
        lon = met["lon2d"]
        lat = met["lat2d"]
        data = met["thick_map_60"]
    ax.pcolormesh(lon, lat, data, shading="auto", cmap="viridis",
                  vmin=0, vmax=vmax_hist)
    ax.set_title(title)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

sm_hist = ScalarMappable(norm=Normalize(vmin=0, vmax=vmax_hist), cmap="viridis")
sm_hist.set_array([])
fig1.colorbar(sm_hist, ax=axes1.ravel().tolist(), orientation="horizontal",
              fraction=0.05, pad=0.05,
              label="OMZ thickness < 60 umol/kg (m)")
fig1.savefig("figure1_omz_thickness_historical.png", dpi=150, bbox_inches="tight")

# Figure 2
ch_gfdl = (model_results["GFDL-ESM4_ssp245"]["thick_map_60"]
           - model_results["GFDL-ESM4_historical"]["thick_map_60"])
ch_mpi = (model_results["MPI-ESM1-2-LR_ssp245"]["thick_map_60"]
          - model_results["MPI-ESM1-2-LR_historical"]["thick_map_60"])
vmax_ch = max(np.nanmax(np.abs(ch_gfdl)), np.nanmax(np.abs(ch_mpi)))
vmax_ch = max(vmax_ch, 1e-6)
vmin_ch = -vmax_ch

fig2, axes2 = plt.subplots(1, 3, figsize=(18, 5))

met = model_results["GFDL-ESM4_ssp245"]
axes2[0].pcolormesh(met["lon2d"], met["lat2d"], ch_gfdl, shading="auto",
                    cmap="RdBu_r", vmin=vmin_ch, vmax=vmax_ch)
axes2[0].set_title("GFDL-ESM4 ssp245 - historical")
axes2[0].set_xlabel("Longitude")
axes2[0].set_ylabel("Latitude")

met = model_results["MPI-ESM1-2-LR_ssp245"]
axes2[1].pcolormesh(met["lon2d"], met["lat2d"], ch_mpi, shading="auto",
                    cmap="RdBu_r", vmin=vmin_ch, vmax=vmax_ch)
axes2[1].set_title("MPI-ESM1-2-LR ssp245 - historical")
axes2[1].set_xlabel("Longitude")
axes2[1].set_ylabel("Latitude")

met = model_results["CESM2_ssp245"]
vmax_cesm = max(np.nanmax(met["thick_map_60"]), 1e-6)
axes2[2].pcolormesh(met["lon2d"], met["lat2d"], met["thick_map_60"],
                    shading="auto", cmap="viridis", vmin=0, vmax=vmax_cesm)
axes2[2].set_title("CESM2 ssp245\n(no historical o2; change not available)")
axes2[2].set_xlabel("Longitude")
axes2[2].set_ylabel("Latitude")

sm_ch = ScalarMappable(norm=Normalize(vmin=vmin_ch, vmax=vmax_ch), cmap="RdBu_r")
sm_ch.set_array([])
fig2.colorbar(sm_ch, ax=[axes2[0], axes2[1]], orientation="horizontal",
              fraction=0.05, pad=0.05,
              label="Change in OMZ thickness (m)")

sm_cesm = ScalarMappable(norm=Normalize(vmin=0, vmax=vmax_cesm), cmap="viridis")
sm_cesm.set_array([])
fig2.colorbar(sm_cesm, ax=axes2[2], orientation="horizontal",
              fraction=0.05, pad=0.05,
              label="CESM2 OMZ thickness (m)")
fig2.savefig("figure2_omz_thickness_change.png", dpi=150, bbox_inches="tight")

plt.close("all")
print("\nSaved figure1_omz_thickness_historical.png and figure2_omz_thickness_change.png")
