# ---------------------------------------------------------------------------
# Computation
# ---------------------------------------------------------------------------
from pathlib import Path
from io import StringIO
import numpy as np
import pandas as pd
import xarray as xr

# ----------------------------- constants -----------------------------------
R_EARTH = 6371000.0  # m
WOA_ROOT = Path("/data1/lorraine-20419/ocean_data/WOA23/climatology")
WOA_MASKS = Path("/data1/lorraine-20419/ocean_data/WOA23/masks")
HIST_ROOT = Path("/data1/lorraine-20419/benchmark_cmip6/historical")
SSP_ROOT = Path("/data1/lorraine-20419/benchmark_cmip6/ssp245")

# ----------------------------- helpers -------------------------------------
def read_woa_landsea(path):
    """Read WOA landsea_01.msk.  Returns (ocean_bool_2d, raw_bottom_level_2d)."""
    arr = np.zeros((180, 360), dtype=int)
    with open(path, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                lat = float(parts[0])
                lon = float(parts[1])
                bot = float(parts[2])
            except ValueError:
                continue
            j = int(np.floor(lat + 90.0))
            i = int(np.floor(lon + 180.0))
            j = np.clip(j, 0, 179)
            i = np.mod(i, 360)
            arr[j, i] = int(bot)
    return arr > 1, arr


def read_woa_basin_mask(path, nlev=102, nlat=180, nlon=360):
    """
    Read WOA basinmask_01.msk.

    The file is a sparse table of ocean columns: latitude, longitude, then up
    to 102 basin codes for the standard levels.  Levels below the sea floor are
    absent/empty at the end of the row, so rows can have fewer than 104 fields.
    """
    arr = np.zeros((nlev, nlat, nlon), dtype=np.int16)
    with open(path, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                lat = float(parts[0])
                lon = float(parts[1])
            except ValueError:
                continue
            j = int(np.floor(lat + 90.0))
            i = int(np.floor(lon + 180.0))
            j = np.clip(j, 0, nlat - 1)
            i = np.mod(i, nlon)
            vals = parts[2:2 + nlev]
            for k, v in enumerate(vals):
                try:
                    b = int(float(v))
                except ValueError:
                    continue
                if b > 0:
                    arr[k, j, i] = b
    return arr


def layer_thickness_from_depth(depth):
    """Approximate layer thickness from depth centres."""
    depth = np.asarray(depth, dtype=float)
    if depth.size == 1:
        return np.array([depth[0]], dtype=float)
    edges = np.empty(depth.size + 1, dtype=float)
    if abs(depth[0]) < 1e-6:
        edges[0] = 0.0
    else:
        edges[0] = depth[0] - 0.5 * (depth[1] - depth[0])
    edges[1:-1] = 0.5 * (depth[:-1] + depth[1:])
    edges[-1] = depth[-1] + 0.5 * (depth[-1] - depth[-2])
    return np.diff(edges)


def woa_layer_thickness(ds, depth):
    """Use depth_bnds if available; otherwise derive from depth centres."""
    for bname in ("depth_bnds", "depth_bounds", "bounds"):
        if bname in ds:
            b = np.asarray(ds[bname].values)
            if b.ndim == 2:
                if b.shape[0] == len(depth) and b.shape[1] == 2:
                    return np.abs(b[:, 1] - b[:, 0])
                if b.shape[1] == len(depth) and b.shape[0] == 2:
                    return np.abs(b[1, :] - b[0, :])
    return layer_thickness_from_depth(depth)


def open_woa(root, subdir, tag, prefix, field="an"):
    """Open a WOA23 annual field and return (dataset, dataarray)."""
    path = Path(root) / subdir / f"woa23_{tag}_{prefix}00_01.nc"
    ds = xr.open_dataset(path, decode_times=False)

    da = None
    for name in ds.data_vars:
        if name == f"{prefix}_{field}":
            da = ds[name]
            break
    if da is None:
        for name in ds.data_vars:
            if name.startswith(prefix) and field in name:
                da = ds[name]
                break
    if da is None:
        for name in ds.data_vars:
            if "depth" in ds[name].dims and ds[name].ndim >= 3:
                da = ds[name]
                break

    if "time" in da.dims:
        da = da.isel(time=0, drop=True)
    da = da.squeeze(drop=True)
    return ds, da


def seawater_density(S, T, p):
    """
    UNESCO 1983 in-situ density.
    S: practical salinity, T: degC, p: pressure in bar.
    """
    T = np.asarray(T, dtype=float)
    S = np.asarray(S, dtype=float)
    p = np.asarray(p, dtype=float)

    rho_w = (
        999.842594
        + 6.793952e-2 * T
        - 9.095290e-3 * T**2
        + 1.001685e-4 * T**3
        - 1.120083e-6 * T**4
        + 6.536332e-9 * T**5
    )
    A = (
        0.824493
        - 4.0899e-3 * T
        + 7.6438e-5 * T**2
        - 8.2467e-7 * T**3
        + 5.3875e-9 * T**4
    )
    B = -5.72466e-3 + 1.0227e-4 * T - 1.6546e-6 * T**2
    C = 4.8314e-4
    rho0 = rho_w + A * S + B * S**1.5 + C * S**2

    Kw = (
        19652.21
        + 148.4206 * T
        - 2.327105 * T**2
        + 1.360477e-2 * T**3
        - 5.155288e-5 * T**4
    )
    K0 = (
        Kw
        + (54.6746 - 0.603459 * T + 1.09987e-2 * T**2 - 6.1670e-5 * T**3) * S
        + (7.944e-2 + 1.6483e-2 * T - 5.3009e-4 * T**2) * S**1.5
    )
    A_k = (
        3.239908
        + 1.43713e-3 * T
        + 1.16092e-4 * T**2
        - 5.77905e-7 * T**3
        + (2.2838e-3 - 1.0981e-5 * T - 1.6078e-6 * T**2) * S
        + 1.91075e-4 * S**1.5
    )
    B_k = (
        8.50935e-5
        - 6.12293e-6 * T
        + 5.2787e-8 * T**2
        + (-9.9348e-7 + 2.0816e-8 * T + 9.1697e-10 * T**2) * S
    )
    K = K0 + A_k * p + B_k * p**2
    return rho0 / (1.0 - p / K)


def density_from_thetao_so(theta, so, depth_m):
    S = np.asarray(so, dtype=float)
    if np.nanmax(S) < 1.0:       # some files store salinity as ~0.035
        S = S * 1000.0
    T = np.asarray(theta, dtype=float)
    p = np.asarray(depth_m, dtype=float)[:, None, None] / 10.0  # bar
    return seawater_density(S, T, p)


MODEL_INFO = {
    "CESM2": dict(
        y="nlat", x="nlon", lev="lev", depth_scale=0.01,
        lat_coord="lat", lon_coord="lon", ocean_mask_from="sftof"
    ),
    "GFDL-ESM4": dict(
        y="lat", x="lon", lev="lev", depth_scale=1.0,
        lat_coord="lat", lon_coord="lon"
    ),
    "MPI-ESM1-2-LR": dict(
        y="j", x="i", lev="lev", depth_scale=1.0,
        lat_coord="latitude", lon_coord="longitude"
    ),
}


def open_cmip_var(root, model, var, y0=None, y1=None):
    """Open one CMIP6 variable for a model, selecting years if requested."""
    base = Path(root) / model / var
    files = sorted(base.glob("*.nc"))
    if not files:
        raise FileNotFoundError(f"No files for {root}/{model}/{var}")
    if len(files) == 1:
        ds = xr.open_dataset(
            files[0], use_cftime=True, decode_times=True, chunks="auto"
        )
    else:
        ds = xr.open_mfdataset(
            files, combine="by_coords", use_cftime=True, decode_times=True,
            chunks="auto", compat="override"
        )

    if y0 is not None and "time" in ds.dims:
        years = ds.time.dt.year
        ds = ds.isel(time=(years >= y0) & (years <= y1))

    if var in ds:
        da = ds[var]
    else:
        da = None
        for name in ds.data_vars:
            if name == var or name.startswith(var):
                da = ds[name]
                break
        if da is None:
            for name in ds.data_vars:
                if ds[name].ndim >= 2:
                    da = ds[name]
                    break
    return ds, da


def load_fixed(root, model, var, y0, y1, fallback_root=None):
    """Load a fixed field, falling back to historical if needed."""
    try:
        ds, da = open_cmip_var(root, model, var, y0, y1)
    except FileNotFoundError:
        if fallback_root is None:
            return None
        ds, da = open_cmip_var(fallback_root, model, var, None, None)

    if "time" in da.dims:
        da = da.mean("time")
    da = da.squeeze(drop=True)
    return da


def get_lat_lon(da, model):
    """Return 2-D lat/lon arrays on the model grid."""
    info = MODEL_INFO[model]
    y, x = info["y"], info["x"]

    lat_name = info["lat_coord"]
    lon_name = info["lon_coord"]
    if lat_name not in da.coords:
        lat_name = "lat" if "lat" in da.coords else "latitude"
    if lon_name not in da.coords:
        lon_name = "lon" if "lon" in da.coords else "longitude"

    lat = da[lat_name]
    lon = da[lon_name]

    if lat.ndim == 1:
        lat2d, lon2d = np.meshgrid(lat.values, lon.values, indexing="ij")
    else:
        lat2d = lat.transpose(y, x).values if set(lat.dims) == {y, x} else lat.values
        lon2d = lon.transpose(y, x).values if set(lon.dims) == {y, x} else lon.values
    return lat2d, lon2d


def process_cmip(root, model, y0, y1, basin_woa=None, woa_depth=None,
                 fallback_root=None):
    """Compute OMZ volumes, thickness map, inventory and basin volumes."""
    info = MODEL_INFO[model]
    lev_name = info["lev"]
    y, x = info["y"], info["x"]

    ds_o2, o2 = open_cmip_var(root, model, "o2", y0, y1)
    ds_t, thetao = open_cmip_var(root, model, "thetao", y0, y1)
    ds_s, so = open_cmip_var(root, model, "so", y0, y1)

    o2_mean = o2.mean("time") if "time" in o2.dims else o2
    theta_mean = thetao.mean("time") if "time" in thetao.dims else thetao
    so_mean = so.mean("time") if "time" in so.dims else so

    o2_mean, theta_mean, so_mean = xr.align(
        o2_mean, theta_mean, so_mean, join="inner"
    )

    o2_da = o2_mean.transpose(lev_name, y, x)
    theta_da = theta_mean.transpose(lev_name, y, x)
    so_da = so_mean.transpose(lev_name, y, x)

    depth_m = o2_da[lev_name].values.astype(float) * info["depth_scale"]

    rho_np = density_from_thetao_so(
        theta_da.values, so_da.values, depth_m
    )
    rho = xr.DataArray(rho_np, coords=theta_da.coords, dims=theta_da.dims)
    o2_umol = o2_da * 1e6 / rho

    vol = load_fixed(root, model, "volcello", y0, y1, fallback_root)
    area = load_fixed(root, model, "areacello", y0, y1, fallback_root)
    deptho = load_fixed(root, model, "deptho", y0, y1, fallback_root)
    sftof = None
    if model == "CESM2":
        sftof = load_fixed(root, model, "sftof", y0, y1, fallback_root)

    vol = vol.transpose(lev_name, y, x)
    area = area.transpose(y, x)

    valid = np.isfinite(o2_da) & (vol > 0)

    if deptho is not None:
        deptho = deptho.transpose(y, x)
        valid = valid & (deptho > 0)

    if model == "CESM2" and sftof is not None:
        sftof = sftof.transpose(y, x)
        sftof_vals = sftof.compute()
        if float(np.nanmax(sftof_vals)) > 2.0:
            sftof = sftof / 100.0      # CMIP6 sftof is commonly in percent
        valid = valid & (sftof > 0)
        area_eff = area * sftof
    else:
        area_eff = area

    valid = valid & np.isfinite(area_eff) & (area_eff > 0)

    o2_umol_np = o2_umol.compute().values
    vol_np = vol.compute().values
    valid_np = valid.compute().values
    area_eff_np = area_eff.compute().values

    thick_np = vol_np / area_eff_np[None, :, :]
    thick_np = np.where(valid_np, thick_np, np.nan)

    omz20 = float(
        np.nansum(vol_np[valid_np & (o2_umol_np < 20.0)]) / 1e15
    )
    omz60 = float(
        np.nansum(vol_np[valid_np & (o2_umol_np < 60.0)]) / 1e15
    )

    thick_omz = np.where(
        valid_np & (o2_umol_np < 60.0), thick_np, 0.0
    ).sum(axis=0)
    ocean2d = valid_np.any(axis=0)
    thick60_map = np.where(ocean2d, thick_omz, np.nan)

    lat2d, lon2d = get_lat_lon(o2_da, model)

    if np.all(np.isnan(thick60_map)):
        max_thick60 = np.nan
        max_lat = np.nan
        max_lon = np.nan
    else:
        idx = np.unravel_index(
            np.nanargmax(thick60_map), thick60_map.shape
        )
        max_thick60 = float(thick60_map[idx])
        max_lat = float(lat2d[idx])
        max_lon = float(lon2d[idx])
        max_lon = ((max_lon + 180.0) % 360.0) - 180.0

    inventory = float(
        np.nansum(
            o2_da.compute().values[valid_np] * vol_np[valid_np]
        )
    )

    basin_omz60 = None
    if basin_woa is not None and woa_depth is not None:
        lon_mod = ((lon2d + 180.0) % 360.0) - 180.0
        i_idx = np.floor(lon_mod + 180.0).astype(int)
        i_idx = np.clip(i_idx, 0, 359)
        j_idx = np.floor(lat2d + 90.0).astype(int)
        j_idx = np.clip(j_idx, 0, 179)

        depth_idx = np.array(
            [np.argmin(np.abs(woa_depth - d)) for d in depth_m]
        )
        basin_model = basin_woa[
            depth_idx[:, None, None],
            j_idx[None, :, :],
            i_idx[None, :, :],
        ]

        basin_omz60 = {}
        for code, name in {1: "Atlantic", 2: "Pacific", 3: "Indian"}.items():
            m = (
                valid_np
                & (o2_umol_np < 60.0)
                & (basin_model == code)
            )
            basin_omz60[name] = float(np.nansum(vol_np[m]) / 1e15)

    return dict(
        omz20=omz20,
        omz60=omz60,
        max_thick60=max_thick60,
        max_lat=max_lat,
        max_lon=max_lon,
        inventory=inventory,
        basin_omz60=basin_omz60,
        thick60_map=thick60_map,
        lat2d=lat2d,
        lon2d=lon2d,
        o2_umol=o2_umol_np,
        vol=vol_np,
        valid=valid_np,
    )

# ----------------------------- WOA23 masks ---------------------------------
landsea_ocean, landsea_bottom = read_woa_landsea(
    WOA_MASKS / "landsea_01.msk"
)
basin_woa = read_woa_basin_mask(
    WOA_MASKS / "basinmask_01.msk", nlev=102
)

# ----------------------------- WOA23 fields --------------------------------
woa_o2_ds, woa_o2 = open_woa(
    WOA_ROOT,
    "oxygen/1971-2000_climate_normal",
    "decav71A0",
    "o",
    "an",
)
_, woa_o2sat = open_woa(
    WOA_ROOT,
    "percent_oxygen_saturation/1971-2000_climate_normal",
    "decav71A0",
    "O",
    "an",
)
_, woa_aou = open_woa(
    WOA_ROOT,
    "apparent_oxygen_utilization/1971-2000_climate_normal",
    "decav71A0",
    "A",
    "an",
)

woa_o2 = woa_o2.transpose("depth", "lat", "lon")
woa_o2sat = woa_o2sat.transpose("depth", "lat", "lon")
woa_aou = woa_aou.transpose("depth", "lat", "lon")

woa_o2, woa_o2sat, woa_aou = xr.align(
    woa_o2, woa_o2sat, woa_aou, join="inner"
)

depth_woa = woa_o2["depth"].values.astype(float)
lat_woa = woa_o2["lat"].values.astype(float)
lon_woa = woa_o2["lon"].values.astype(float)

thick_woa = woa_layer_thickness(woa_o2_ds, depth_woa)

# trim/pad basin mask to the oxygen depth axis
if basin_woa.shape[0] < len(depth_woa):
    pad = np.zeros(
        (len(depth_woa) - basin_woa.shape[0],) + basin_woa.shape[1:],
        dtype=basin_woa.dtype,
    )
    basin_woa = np.concatenate([basin_woa, pad], axis=0)
elif basin_woa.shape[0] > len(depth_woa):
    basin_woa = basin_woa[:len(depth_woa)]

# WOA cell areas and volumes
dlat = np.deg2rad(np.abs(np.diff(lat_woa)).mean()) if len(lat_woa) > 1 else np.deg2rad(1.0)
dlon = np.deg2rad(np.abs(np.diff(lon_woa)).mean()) if len(lon_woa) > 1 else np.deg2rad(1.0)
area_woa_1d = R_EARTH**2 * np.cos(np.deg2rad(lat_woa)) * dlat * dlon
area_woa = np.repeat(area_woa_1d[:, None], len(lon_woa), axis=1)

o2_woa_np = woa_o2.values
o2sat_woa_np = woa_o2sat.values
aou_woa_np = woa_aou.values

if np.nanmax(o2sat_woa_np) <= 1.5:
    o2sat_woa_np = o2sat_woa_np * 100.0

valid_woa = np.isfinite(o2_woa_np) & landsea_ocean[None, :, :]
vol_woa_np = (
    area_woa[None, :, :]
    * thick_woa[:, None, None]
    * valid_woa
)

woa_omz20_1e15 = float(
    np.nansum(vol_woa_np[valid_woa & (o2_woa_np < 20.0)]) / 1e15
)
woa_omz60_1e15 = float(
    np.nansum(vol_woa_np[valid_woa & (o2_woa_np < 60.0)]) / 1e15
)

woa_thick60_map = np.where(
    valid_woa & (o2_woa_np < 60.0),
    thick_woa[:, None, None],
    0.0,
).sum(axis=0)
ocean2d_woa = valid_woa.any(axis=0)
woa_thick60_map = np.where(ocean2d_woa, woa_thick60_map, np.nan)

idx_woa = np.unravel_index(
    np.nanargmax(woa_thick60_map), woa_thick60_map.shape
)
woa_max_thick60 = float(woa_thick60_map[idx_woa])
woa_max_lat = float(lat_woa[idx_woa[0]])
woa_max_lon = float(lon_woa[idx_woa[1]])

omz_mask_woa = valid_woa & (o2_woa_np < 60.0)
w_omz = vol_woa_np * omz_mask_woa

w_aou = w_omz * np.isfinite(aou_woa_np)
woa_aou_omz60_mean = float(
    np.nansum(aou_woa_np * w_aou) / np.nansum(w_aou)
)

w_o2sat = w_omz * np.isfinite(o2sat_woa_np)
woa_o2sat_omz60_mean = float(
    np.nansum(o2sat_woa_np * w_o2sat) / np.nansum(w_o2sat)
)

woa_basin_omz60 = {}
for code, name in {1: "Atlantic", 2: "Pacific", 3: "Indian"}.items():
    m = omz_mask_woa & (basin_woa == code)
    woa_basin_omz60[name] = float(np.nansum(vol_woa_np[m]) / 1e15)

# ----------------------------- CMIP6 historical -----------------------------
hist = {}
for model in ["GFDL-ESM4", "MPI-ESM1-2-LR"]:
    hist[model] = process_cmip(
        HIST_ROOT, model, 1971, 2000,
        basin_woa=basin_woa,
        woa_depth=depth_woa,
    )

# ----------------------------- CMIP6 ssp245 ---------------------------------
ssp = {}
for model in ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]:
    ssp[model] = process_cmip(
        SSP_ROOT, model, 2071, 2100,
        basin_woa=basin_woa,
        woa_depth=depth_woa,
        fallback_root=HIST_ROOT,
    )

# ----------------------------- changes --------------------------------------
pct_change_20 = {}
pct_change_60 = {}
for model in ["GFDL-ESM4", "MPI-ESM1-2-LR"]:
    pct_change_20[model] = (
        (ssp[model]["omz20"] - hist[model]["omz20"])
        / hist[model]["omz20"]
        * 100.0
    )
    pct_change_60[model] = (
        (ssp[model]["omz60"] - hist[model]["omz60"])
        / hist[model]["omz60"]
        * 100.0
    )

# closeness to WOA at each threshold
err60 = {
    m: abs(hist[m]["omz60"] - woa_omz60_1e15)
    for m in hist
}
err20 = {
    m: abs(hist[m]["omz20"] - woa_omz20_1e15)
    for m in hist
}
better_model_60 = min(err60, key=err60.get)
better_model_20 = min(err20, key=err20.get)

agree_sign_20 = (
    np.sign(pct_change_20["GFDL-ESM4"])
    == np.sign(pct_change_20["MPI-ESM1-2-LR"])
)
agree_sign_60 = (
    np.sign(pct_change_60["GFDL-ESM4"])
    == np.sign(pct_change_60["MPI-ESM1-2-LR"])
)

# convenience aliases for the answer
omz20_woa = woa_omz20_1e15
omz60_woa = woa_omz60_1e15
omz20_hist_gfdl = hist["GFDL-ESM4"]["omz20"]
omz60_hist_gfdl = hist["GFDL-ESM4"]["omz60"]
omz20_hist_mpi = hist["MPI-ESM1-2-LR"]["omz20"]
omz60_hist_mpi = hist["MPI-ESM1-2-LR"]["omz60"]
omz20_ssp_cesm2 = ssp["CESM2"]["omz20"]
omz60_ssp_cesm2 = ssp["CESM2"]["omz60"]
omz20_ssp_gfdl = ssp["GFDL-ESM4"]["omz20"]
omz60_ssp_gfdl = ssp["GFDL-ESM4"]["omz60"]
omz20_ssp_mpi = ssp["MPI-ESM1-2-LR"]["omz20"]
omz60_ssp_mpi = ssp["MPI-ESM1-2-LR"]["omz60"]

max_thick60_woa = woa_max_thick60
max_thick60_hist_gfdl = hist["GFDL-ESM4"]["max_thick60"]
max_thick60_hist_mpi = hist["MPI-ESM1-2-LR"]["max_thick60"]
max_thick60_ssp_cesm2 = ssp["CESM2"]["max_thick60"]
max_thick60_ssp_gfdl = ssp["GFDL-ESM4"]["max_thick60"]
max_thick60_ssp_mpi = ssp["MPI-ESM1-2-LR"]["max_thick60"]

inventory_hist_gfdl = hist["GFDL-ESM4"]["inventory"]
inventory_hist_mpi = hist["MPI-ESM1-2-LR"]["inventory"]
inventory_ssp_cesm2 = ssp["CESM2"]["inventory"]
inventory_ssp_gfdl = ssp["GFDL-ESM4"]["inventory"]
inventory_ssp_mpi = ssp["MPI-ESM1-2-LR"]["inventory"]

basin_omz60_woa = woa_basin_omz60
basin_omz60_hist_gfdl = hist["GFDL-ESM4"]["basin_omz60"]
basin_omz60_hist_mpi = hist["MPI-ESM1-2-LR"]["basin_omz60"]

woa_aou_omz60 = woa_aou_omz60_mean
woa_o2sat_omz60 = woa_o2sat_omz60_mean
