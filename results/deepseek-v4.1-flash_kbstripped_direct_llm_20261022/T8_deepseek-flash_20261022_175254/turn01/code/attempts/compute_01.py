import os
import glob
import warnings
import numpy as np
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore")

HIST_ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"
GLODAP_ROOT = "/data1/lorraine-20419/ocean_data/GLODAPv2.2016b_MappedClimatologies"

MODELS = {
    "CESM2": dict(
        y_dim="nlat", x_dim="nlon", lat="lat", lon="lon",
        depth_dim="lev", depth_scale=0.01, grid="curvilinear",
    ),
    "GFDL-ESM4": dict(
        y_dim="lat", x_dim="lon", lat="lat", lon="lon",
        depth_dim="lev", depth_scale=1.0, grid="regular",
    ),
    "MPI-ESM1-2-LR": dict(
        y_dim="j", x_dim="i", lat="latitude", lon="longitude",
        depth_dim="lev", depth_scale=1.0, grid="curvilinear",
    ),
}

profile_depths = [0, 100, 500, 1000, 2000]


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------

def _find_nc_files(model, var, root=HIST_ROOT):
    files = sorted(glob.glob(os.path.join(root, model, var, "*.nc")))
    if not files:
        files = sorted(glob.glob(os.path.join(root, model, "**", f"{var}*.nc"), recursive=True))
    if not files:
        raise FileNotFoundError(f"No files for {model} / {var}")
    return files


def open_cmip(model, var, root=HIST_ROOT):
    files = _find_nc_files(model, var, root=root)
    with xr.open_dataset(files[0], use_cftime=True, decode_times=True) as tmp:
        has_time = "time" in tmp.dims
    chunks = {"time": 12} if has_time else {}
    if len(files) == 1:
        return xr.open_dataset(
            files[0], use_cftime=True, decode_times=True, chunks=chunks
        )
    return xr.open_mfdataset(
        files, combine="by_coords", use_cftime=True, decode_times=True,
        chunks=chunks,
    )


def _select_time(da, year=None, period=None):
    if "time" not in da.dims:
        return da
    yrs = np.asarray(da.time.dt.year.values)
    if year is not None:
        idx = np.where(yrs == year)[0]
    elif period is not None:
        idx = np.where((yrs >= period[0]) & (yrs <= period[1]))[0]
    else:
        return da
    return da.isel(time=idx)


def mean_field(da, year=None, period=None):
    if "time" in da.dims:
        da = _select_time(da, year=year, period=period)
        da = da.mean("time")
    return da


def to_numpy_3d(da, depth_dim):
    dims = list(da.dims)
    if depth_dim not in dims:
        raise ValueError(f"Depth dimension {depth_dim!r} not in {dims}")
    other = [d for d in dims if d != depth_dim]
    da = da.transpose(depth_dim, *other)
    return np.asarray(da.values, dtype=float)


def get_depth_m(ds, cfg):
    lev = ds[cfg["depth_dim"]]
    vals = np.asarray(lev.values, dtype=float)
    units = str(lev.attrs.get("units", "")).lower()
    scale = cfg["depth_scale"]
    if "cm" in units:
        scale = 0.01
    elif units in ("m", "meter", "meters", "metre", "metres"):
        scale = 1.0
    return vals * scale


def get_lat_lon_2d(ds, cfg):
    lat = ds[cfg["lat"]]
    lon = ds[cfg["lon"]]
    if lat.ndim == 2:
        lat2d = np.asarray(lat.values, dtype=float)
        lon2d = np.asarray(lon.values, dtype=float)
    else:
        lat2d, lon2d = np.meshgrid(
            np.asarray(lat.values, dtype=float),
            np.asarray(lon.values, dtype=float),
            indexing="ij",
        )
    return lat2d, lon2d % 360.0


def get_area(model, cfg):
    ds = open_cmip(model, "areacello")
    da = ds["areacello"]
    if "time" in da.dims:
        da = da.mean("time")
    y, x = cfg["y_dim"], cfg["x_dim"]
    if y in da.dims and x in da.dims:
        da = da.transpose(y, x)
    else:
        da = da.transpose(da.dims[-2], da.dims[-1])
    area = np.asarray(da.values, dtype=float)
    ds.close()
    return area


def get_volcello(model, cfg, year=None, period=None):
    ds = open_cmip(model, "volcello")
    da = ds["volcello"]
    if "time" in da.dims:
        if year is not None:
            da = _select_time(da, year=year)
        elif period is not None:
            da = _select_time(da, period=period)
        da = da.mean("time")
    out = to_numpy_3d(da, cfg["depth_dim"])
    ds.close()
    return out


# ---------------------------------------------------------------------------
# EOS-80 in-situ density, pressure in bars
# ---------------------------------------------------------------------------

def rho_eos80(temp, salt, depth_m):
    T = np.asarray(temp, dtype=float) * 1.00024
    S = np.asarray(salt, dtype=float)
    p = np.asarray(depth_m, dtype=float)
    if p.ndim == 1:
        p = p.reshape((-1,) + (1,) * (T.ndim - 1))
    p = p / 10.0  # m -> dbar -> bar

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
        + S * (54.6746 - 0.603459 * T + 1.09987e-2 * T**2 - 6.1670e-5 * T**3)
        + S**1.5 * (7.944e-2 + 1.6483e-2 * T - 5.3009e-4 * T**2)
    )
    K = (
        K0
        + p * (3.239908 + 1.43713e-3 * T + 1.16092e-4 * T**2 - 5.77905e-7 * T**3)
        + p * S * (2.209460e-3 - 2.205960e-3 * T - 2.181160e-3 * T**2)
        + p * S**1.5 * (-3.01700e-4 + 3.58170e-5 * T)
        + p**2 * (8.50935e-5 - 6.12293e-6 * T + 5.2787e-8 * T**2)
        + p**2 * S * (-9.9348e-7 + 2.0816e-8 * T + 9.1697e-10 * T**2)
    )
    return rho0 / (1.0 - p / K)


def carbon_to_umolkg(carbon, temp, salt, depth_m, units):
    c = np.asarray(carbon, dtype=float)
    u = str(units).lower().replace(" ", "")
    if "umol" in u or "micromol" in u:
        return c
    if "molkg-1" in u and "m-3" not in u:
        return c * 1e6
    # CMIP6 standard: mol m-3 -> umol kg-1
    rho = rho_eos80(temp, salt, depth_m)
    return c * 1e6 / rho


# ---------------------------------------------------------------------------
# weighting / averaging / regridding
# ---------------------------------------------------------------------------

def weighted_cell_mean(field, weights):
    field = np.asarray(field, dtype=float)
    weights = np.broadcast_to(np.asarray(weights, dtype=float), field.shape)
    valid = np.isfinite(field) & np.isfinite(weights) & (weights > 0)
    num = np.nansum(np.where(valid, field * weights, 0.0), axis=0)
    den = np.nansum(np.where(valid, weights, 0.0), axis=0)
    out = np.full(num.shape, np.nan)
    good = den > 0
    out[good] = num[good] / den[good]
    return out


def volume_global_mean(field, weights):
    field = np.asarray(field, dtype=float)
    weights = np.broadcast_to(np.asarray(weights, dtype=float), field.shape)
    valid = np.isfinite(field) & np.isfinite(weights) & (weights > 0)
    if not np.any(valid):
        return np.nan
    return np.nansum(field[valid] * weights[valid]) / np.nansum(weights[valid])


def area_weighted_mean(field2d, area2d):
    field2d = np.asarray(field2d, dtype=float)
    area2d = np.asarray(area2d, dtype=float)
    valid = np.isfinite(field2d) & np.isfinite(area2d)
    if not np.any(valid):
        return np.nan
    return np.nansum(field2d[valid] * area2d[valid]) / np.nansum(area2d[valid])


def layer_thickness(depth, max_depth=50.0):
    d = np.asarray(depth, dtype=float)
    order = np.argsort(d)
    ds = d[order]
    n = len(ds)
    bounds = np.empty(n + 1, dtype=float)
    if n == 1:
        bounds[:] = [0.0, ds[0]]
    else:
        bounds[0] = max(0.0, ds[0] - 0.5 * (ds[1] - ds[0]))
        bounds[1:-1] = 0.5 * (ds[:-1] + ds[1:])
        bounds[-1] = ds[-1] + 0.5 * (ds[-1] - ds[-2])
    upper = np.maximum(bounds[:-1], 0.0)
    lower = np.minimum(bounds[1:], max_depth)
    thick = np.maximum(0.0, lower - upper)
    out = np.empty(n, dtype=float)
    out[order] = thick
    return out


def _latlon_xyz(lat, lon):
    lat = np.deg2rad(np.asarray(lat, dtype=float))
    lon = np.deg2rad(np.asarray(lon, dtype=float) % 360.0)
    return np.stack(
        [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)],
        axis=-1,
    )


def regrid_nearest(src_field, src_lat, src_lon, tgt_lat, tgt_lon, tgt_valid=None):
    src_field = np.asarray(src_field, dtype=float)
    src_lat = np.asarray(src_lat, dtype=float)
    src_lon = np.asarray(src_lon, dtype=float)
    valid_src = (
        np.isfinite(src_field)
        & np.isfinite(src_lat)
        & np.isfinite(src_lon)
    )
    src_vals = src_field[valid_src]
    out = np.full(np.asarray(tgt_lat).shape, np.nan)

    if src_vals.size == 0:
        return out

    src_xyz = _latlon_xyz(src_lat[valid_src], src_lon[valid_src])
    tree = cKDTree(src_xyz)

    tgt_lat = np.asarray(tgt_lat, dtype=float)
    tgt_lon = np.asarray(tgt_lon, dtype=float)
    if tgt_valid is None:
        use = np.ones(tgt_lat.shape, dtype=bool)
    else:
        use = np.asarray(tgt_valid, dtype=bool)
    use = use & np.isfinite(tgt_lat) & np.isfinite(tgt_lon)

    q = _latlon_xyz(tgt_lat[use], tgt_lon[use])
    _, idx = tree.query(q, k=1)
    out[use] = src_vals[idx]
    return out


def interp_to_depth(field3d, depth, target_depth):
    field3d = np.asarray(field3d, dtype=float)
    depth = np.asarray(depth, dtype=float)
    order = np.argsort(depth)
    depth = depth[order]
    field3d = field3d[order]

    out = np.full(field3d.shape[1:], np.nan)
    if target_depth <= depth[0]:
        return field3d[0].copy()
    if target_depth > depth[-1]:
        return out

    i = np.searchsorted(depth, target_depth) - 1
    i = int(np.clip(i, 0, len(depth) - 2))
    d0, d1 = depth[i], depth[i + 1]
    w = (target_depth - d0) / (d1 - d0)
    return field3d[i] * (1.0 - w) + field3d[i + 1] * w


# ---------------------------------------------------------------------------
# GLODAP reader
# ---------------------------------------------------------------------------

def _find_dim(ds, candidates):
    for name in candidates:
        if name in ds.dims or name in ds.coords:
            return name
    for d in list(ds.dims):
        if any(c in d.lower() for c in candidates):
            return d
    return None


def open_glodap(varfile, varname):
    path = os.path.join(GLODAP_ROOT, f"GLODAPv2.2016b.{varfile}.nc")
    ds = xr.open_dataset(path)
    da = ds[varname]

    lat_name = _find_dim(ds, ["lat", "latitude"])
    lon_name = _find_dim(ds, ["lon", "longitude"])
    depth_name = _find_dim(ds, ["depth_surface", "depth"])

    if depth_name is None:
        depth_name = [d for d in da.dims if "depth" in d.lower()][0]
    if lat_name is None:
        lat_name = [d for d in da.dims if "lat" in d.lower()][0]
    if lon_name is None:
        lon_name = [d for d in da.dims if "lon" in d.lower()][0]

    if "Depth" in ds:
        depth = np.asarray(ds["Depth"].values, dtype=float)
    else:
        depth = np.asarray(ds[depth_name].values, dtype=float)
    lat = np.asarray(ds[lat_name].values, dtype=float)
    lon = np.asarray(ds[lon_name].values, dtype=float) % 360.0

    rename = {}
    if depth_name in da.dims and depth_name != "depth":
        rename[depth_name] = "depth"
    if lat_name in da.dims and lat_name != "lat":
        rename[lat_name] = "lat"
    if lon_name in da.dims and lon_name != "lon":
        rename[lon_name] = "lon"
    if rename:
        da = da.rename(rename)

    da = da.assign_coords(
        depth=("depth", depth),
        lat=("lat", lat),
        lon=("lon", lon),
    )
    da = da.sortby("depth").sortby("lat").sortby("lon")
    da = da.transpose("depth", "lat", "lon")
    ds.close()
    return da


glodap_dic_da = open_glodap("TCO2", "TCO2")
glodap_talk_da = open_glodap("TAlk", "TAlk")
glodap_ph_da = open_glodap("pHtsinsitutp", "pHtsinsitutp")

glodap_depth = np.asarray(glodap_dic_da["depth"].values, dtype=float)
glodap_lat = np.asarray(glodap_dic_da["lat"].values, dtype=float)
glodap_lon = np.asarray(glodap_dic_da["lon"].values, dtype=float)

glodap_dic = np.asarray(glodap_dic_da.values, dtype=float)
glodap_talk = np.asarray(glodap_talk_da.reindex_like(glodap_dic_da).values, dtype=float)
glodap_ph = np.asarray(glodap_ph_da.reindex_like(glodap_dic_da).values, dtype=float)

glodap_thick_full = layer_thickness(glodap_depth, max_depth=50.0)
lev_g50 = glodap_depth <= 50.0

glodap_dic_50 = glodap_dic[lev_g50]
glodap_talk_50 = glodap_talk[lev_g50]
glodap_ph_50 = glodap_ph[lev_g50]
glodap_thick_50 = glodap_thick_full[lev_g50]

glodap_dic_map = weighted_cell_mean(glodap_dic_50, glodap_thick_50[:, None, None])
glodap_talk_map = weighted_cell_mean(glodap_talk_50, glodap_thick_50[:, None, None])
glodap_ph_map = weighted_cell_mean(glodap_ph_50, glodap_thick_50[:, None, None])

glodap_area = (
    np.cos(np.deg2rad(glodap_lat))[:, None]
    * np.ones((1, len(glodap_lon)), dtype=float)
)
glodap_vol_50 = glodap_area[None, :, :] * glodap_thick_50[:, None, None]

glodap_dic_global = volume_global_mean(glodap_dic_50, glodap_vol_50)
glodap_talk_global = volume_global_mean(glodap_talk_50, glodap_vol_50)
glodap_ph_global = volume_global_mean(glodap_ph_50, glodap_vol_50)

glodap_lat2d, glodap_lon2d = np.meshgrid(
    glodap_lat, glodap_lon, indexing="ij"
)
glodap_valid_map = np.isfinite(glodap_dic_map)


# ---------------------------------------------------------------------------
# process each CMIP6 model
# ---------------------------------------------------------------------------

def process_model(model):
    cfg = MODELS[model]

    ds_d = open_cmip(model, "dissic")
    ds_t = open_cmip(model, "talk")
    ds_p = open_cmip(model, "ph")
    ds_th = open_cmip(model, "thetao")
    ds_s = open_cmip(model, "so")

    depth_m_raw = get_depth_m(ds_d, cfg)
    order = np.argsort(depth_m_raw)
    depth_m = depth_m_raw[order]

    lat2d, lon2d = get_lat_lon_2d(ds_d, cfg)
    area = get_area(model, cfg)
    if area.shape != lat2d.shape:
        area = area.T

    def get_3d(da, year=None, period=None):
        return to_numpy_3d(
            mean_field(da, year=year, period=period),
            cfg["depth_dim"],
        )[order]

    dissic_2002 = get_3d(ds_d["dissic"], year=2002)
    talk_2002 = get_3d(ds_t["talk"], year=2002)
    ph_2002 = get_3d(ds_p["ph"], year=2002)
    theta_2002 = get_3d(ds_th["thetao"], year=2002)
    so_2002 = get_3d(ds_s["so"], year=2002)

    dissic_period = get_3d(ds_d["dissic"], period=(1971, 2000))
    theta_period = get_3d(ds_th["thetao"], period=(1971, 2000))
    so_period = get_3d(ds_s["so"], period=(1971, 2000))

    dic_units = ds_d["dissic"].attrs.get("units", "")
    talk_units = ds_t["talk"].attrs.get("units", "")

    dic_2002_umol = carbon_to_umolkg(
        dissic_2002, theta_2002, so_2002, depth_m[:, None, None], dic_units
    )
    talk_2002_umol = carbon_to_umolkg(
        talk_2002, theta_2002, so_2002, depth_m[:, None, None], talk_units
    )
    dic_period_umol = carbon_to_umolkg(
        dissic_period, theta_period, so_period, depth_m[:, None, None], dic_units
    )

    vol_2002 = get_volcello(model, cfg, year=2002)[order]
    vol_period = get_volcello(model, cfg, period=(1971, 2000))[order]

    lev50 = depth_m <= 50.0

    dic_map = weighted_cell_mean(dic_2002_umol[lev50], vol_2002[lev50])
    talk_map = weighted_cell_mean(talk_2002_umol[lev50], vol_2002[lev50])
    ph_map = weighted_cell_mean(ph_2002[lev50], vol_2002[lev50])

    dic_global = volume_global_mean(dic_2002_umol[lev50], vol_2002[lev50])
    talk_global = volume_global_mean(talk_2002_umol[lev50], vol_2002[lev50])
    ph_global = volume_global_mean(ph_2002[lev50], vol_2002[lev50])

    dic_global_map = area_weighted_mean(dic_map, area)
    talk_global_map = area_weighted_mean(talk_map, area)
    ph_global_map = area_weighted_mean(ph_map, area)

    dic_period_global = volume_global_mean(
        dic_period_umol[lev50], vol_period[lev50]
    )
    dic_period_map = weighted_cell_mean(
        dic_period_umol[lev50], vol_period[lev50]
    )
    dic_period_global_map = area_weighted_mean(dic_period_map, area)

    for ds in (ds_d, ds_t, ds_p, ds_th, ds_s):
        ds.close()

    return dict(
        depth_m=depth_m,
        lat2d=lat2d,
        lon2d=lon2d,
        area=area,
        dic_2002_umol=dic_2002_umol,
        talk_2002_umol=talk_2002_umol,
        ph_2002=ph_2002,
        dic_period_umol=dic_period_umol,
        vol_2002=vol_2002,
        vol_period=vol_period,
        dic_map=dic_map,
        talk_map=talk_map,
        ph_map=ph_map,
        dic_global=dic_global,
        talk_global=talk_global,
        ph_global=ph_global,
        dic_global_map=dic_global_map,
        talk_global_map=talk_global_map,
        ph_global_map=ph_global_map,
        dic_period_global=dic_period_global,
        dic_period_global_map=dic_period_global_map,
    )


model_results = {m: process_model(m) for m in MODELS}

for m in MODELS:
    model_results[m]["dic_bias_global"] = (
        model_results[m]["dic_global"] - glodap_dic_global
    )
    model_results[m]["talk_bias_global"] = (
        model_results[m]["talk_global"] - glodap_talk_global
    )
    model_results[m]["ph_bias_global"] = (
        model_results[m]["ph_global"] - glodap_ph_global
    )


# ---------------------------------------------------------------------------
# map metrics and bias maps on the GLODAP 1-degree grid
# ---------------------------------------------------------------------------

map_metrics = {}
dic_bias_maps = {}

for m in MODELS:
    src = model_results[m]
    model_on_glodap = regrid_nearest(
        src["dic_map"],
        src["lat2d"],
        src["lon2d"],
        glodap_lat2d,
        glodap_lon2d,
        glodap_valid_map,
    )
    valid = glodap_valid_map & np.isfinite(model_on_glodap)
    diff = model_on_glodap - glodap_dic_map

    bias_map = np.full(glodap_dic_map.shape, np.nan)
    bias_map[valid] = diff[valid]
    dic_bias_maps[m] = bias_map

    map_metrics[m] = dict(
        mean_bias=float(np.mean(diff[valid])),
        rmse=float(np.sqrt(np.mean(diff[valid] ** 2))),
        corr=float(np.corrcoef(model_on_glodap[valid], glodap_dic_map[valid])[0, 1])
        if np.sum(valid) > 1 else np.nan,
    )


# ---------------------------------------------------------------------------
# vertical profile of DIC bias at selected GLODAP depths
# ---------------------------------------------------------------------------

vertical_bias = {}

for m in MODELS:
    src = model_results[m]
    vals = {}
    for d in profile_depths:
        idx = int(np.argmin(np.abs(glodap_depth - d)))
        gl_slice = glodap_dic[idx]
        tgt_valid = np.isfinite(gl_slice)

        model_slice = interp_to_depth(src["dic_2002_umol"], src["depth_m"], d)
        model_on_glodap = regrid_nearest(
            model_slice,
            src["lat2d"],
            src["lon2d"],
            glodap_lat2d,
            glodap_lon2d,
            tgt_valid,
        )
        valid = tgt_valid & np.isfinite(model_on_glodap)
        vals[d] = (
            float(np.mean(model_on_glodap[valid] - gl_slice[valid]))
            if np.any(valid) else np.nan
        )
    vertical_bias[m] = vals


# ---------------------------------------------------------------------------
# period change and ratio to the 2002 model-minus-GLODAP bias
# ---------------------------------------------------------------------------

period_table = {}
for m in MODELS:
    change = (
        model_results[m]["dic_global"]
        - model_results[m]["dic_period_global"]
    )
    bias = model_results[m]["dic_bias_global"]
    ratio = float(change / bias) if abs(bias) > 1e-12 else np.nan
    period_table[m] = dict(
        dic_2002=float(model_results[m]["dic_global"]),
        dic_1971_2000=float(model_results[m]["dic_period_global"]),
        change=float(change),
        dic_bias_vs_glodap=float(bias),
        ratio=ratio,
    )


# ---------------------------------------------------------------------------
# tidy result tables
# ---------------------------------------------------------------------------

source_rows = []
for m in MODELS:
    source_rows.append(
        dict(
            source=m,
            DIC=model_results[m]["dic_global"],
            TAlk=model_results[m]["talk_global"],
            pH=model_results[m]["ph_global"],
        )
    )
source_rows.append(
    dict(source="GLODAP", DIC=glodap_dic_global, TAlk=glodap_talk_global, pH=glodap_ph_global)
)
global_mean_df = pd.DataFrame(source_rows).set_index("source")

global_bias_df = pd.DataFrame(
    {
        m: {
            "DIC": model_results[m]["dic_bias_global"],
            "TAlk": model_results[m]["talk_bias_global"],
            "pH": model_results[m]["ph_bias_global"],
        }
        for m in MODELS
    }
).T

global_mean_rounded = global_mean_df.copy()
global_mean_rounded[["DIC", "TAlk"]] = global_mean_rounded[["DIC", "TAlk"]].round(0)
global_mean_rounded["pH"] = global_mean_rounded["pH"].round(3)

global_bias_rounded = global_bias_df.copy()
global_bias_rounded[["DIC", "TAlk"]] = global_bias_rounded[["DIC", "TAlk"]].round(0)
global_bias_rounded["pH"] = global_bias_rounded["pH"].round(3)

map_metrics_df = pd.DataFrame.from_dict(map_metrics, orient="index")
map_metrics_rounded = map_metrics_df.copy()
map_metrics_rounded[["mean_bias", "rmse"]] = map_metrics_rounded[["mean_bias", "rmse"]].round(0)
map_metrics_rounded["corr"] = map_metrics_rounded["corr"].round(3)

vertical_bias_df = pd.DataFrame.from_dict(vertical_bias, orient="index")
vertical_bias_df = vertical_bias_df[profile_depths]
vertical_bias_rounded = vertical_bias_df.round(0)

period_df = pd.DataFrame.from_dict(period_table, orient="index")
period_rounded = period_df.copy()
for c in ["dic_2002", "dic_1971_2000", "change", "dic_bias_vs_glodap"]:
    period_rounded[c] = period_rounded[c].round(0)
period_rounded["ratio"] = period_rounded["ratio"].round(2)

sign_df = np.sign(global_bias_df[["DIC", "TAlk", "pH"]]).astype(int)

# aggregate error score: DIC/TAlk in umol kg-1, pH scaled by 1000
model_scores = pd.Series(
    {
        m: (
            abs(global_bias_df.loc[m, "DIC"])
            + abs(global_bias_df.loc[m, "TAlk"])
            + 1000.0 * abs(global_bias_df.loc[m, "pH"])
            + map_metrics_df.loc[m, "rmse"]
        )
        for m in MODELS
    }
)
best_model = model_scores.idxmin()

# sign-consistency statements
same_sign_across_models = {
    c: bool(sign_df[c].nunique() == 1) for c in sign_df.columns
}
within_model_consistency = {}
for m in MODELS:
    s = sign_df.loc[m]
    within_model_consistency[m] = bool(
        (s["DIC"] == s["TAlk"]) and (s["pH"] == -s["DIC"])
    )

period_implication = (
    "The 2002-minus-1971-2000 DIC change is of the same order as the "
    "model-minus-GLODAP DIC offset for at least some models. Because GLODAP "
    "is a 1972-2013 climatology while the model fields are for 2002, the "
    "offset cannot be interpreted as a pure model carbon-system error; it "
    "also contains a time-mismatch contribution."
)
