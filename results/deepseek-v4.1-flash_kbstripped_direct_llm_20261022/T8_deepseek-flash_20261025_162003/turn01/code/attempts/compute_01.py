import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

CMIP_ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"
GLODAP_ROOT = "/data1/lorraine-20419/ocean_data/GLODAPv2.2016b_MappedClimatologies"

MODELS = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
TARGET_DEPTHS = [0, 100, 500, 1000, 2000]

MODEL_INFO = {
    "CESM2": {"lat_coord": "lat", "lon_coord": "lon", "depth_coord": "lev"},
    "GFDL-ESM4": {"lat_coord": "lat", "lon_coord": "lon", "depth_coord": "lev"},
    "MPI-ESM1-2-LR": {"lat_coord": "latitude", "lon_coord": "longitude", "depth_coord": "lev"},
}


def _open_cmip(model, var):
    files = sorted(glob.glob(os.path.join(CMIP_ROOT, model, var, "*.nc")))
    if not files:
        files = sorted(glob.glob(os.path.join(CMIP_ROOT, model, var.lower(), "*.nc")))
    if not files:
        raise FileNotFoundError(f"No CMIP6 files found for {model} / {var}")
    if len(files) == 1:
        return xr.open_dataset(files[0], use_cftime=True, decode_times=True)
    return xr.open_mfdataset(
        files, combine="by_coords", use_cftime=True, decode_times=True, parallel=False
    )


def _get_dataarray(ds, name):
    if name in ds.data_vars:
        return ds[name]
    for v in ds.data_vars:
        if v.lower() == name.lower():
            return ds[v]
    raise KeyError(f"{name} not found in dataset; variables are {list(ds.data_vars)}")


def _select_years(da, y0, y1):
    if "time" not in da.dims:
        return da
    yrs = da.time.dt.year
    return da.sel(time=(yrs >= y0) & (yrs <= y1))


def _lev_m(da, model):
    lev = da[MODEL_INFO[model]["depth_coord"]]
    units = str(lev.attrs.get("units", "")).lower()
    if "centi" in units or units.strip() == "cm":
        return lev * 0.01
    return lev


def _get_latlon2d(da, model):
    info = MODEL_INFO[model]
    lat = da[info["lat_coord"]]
    lon = da[info["lon_coord"]]
    if lat.ndim == 1 and lon.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon.values, lat.values)
    else:
        lat2d = lat.values
        lon2d = lon.values
    return lat2d, lon2d


def _seawater_density(T, S, p_bar):
    """UNESCO 1983 in-situ density approximation, pressure in bar, density in kg/m3."""
    rho_w = (
        999.842594
        + 6.793952e-2 * T
        - 9.095290e-3 * T**2
        + 1.001685e-4 * T**3
        - 1.120083e-6 * T**4
        + 6.536332e-9 * T**5
    )
    A = (
        8.24493e-1
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
    rho = rho0 / (1.0 - p_bar / K0)
    return rho


def _remap_to_target(lat2d, lon2d, values, target_lat, target_lon):
    """Nearest-neighbour remap of a 2-D field to a regular lat/lon target grid."""
    lat2d = np.asarray(lat2d)
    lon2d = np.asarray(lon2d)
    values = np.asarray(values)

    valid = np.isfinite(values) & np.isfinite(lat2d) & np.isfinite(lon2d)
    if not valid.any():
        return np.full((len(target_lat), len(target_lon)), np.nan)

    lat = np.deg2rad(lat2d[valid])
    lon = np.deg2rad(lon2d[valid])
    xyz = np.column_stack(
        [
            np.cos(lat) * np.cos(lon),
            np.cos(lat) * np.sin(lon),
            np.sin(lat),
        ]
    )
    tree = cKDTree(xyz)

    tlon, tlat = np.meshgrid(np.deg2rad(target_lon), np.deg2rad(target_lat))
    q = np.column_stack(
        [
            np.cos(tlat).ravel() * np.cos(tlon).ravel(),
            np.cos(tlat).ravel() * np.sin(tlon).ravel(),
            np.sin(tlat).ravel(),
        ]
    )
    _, idx = tree.query(q, k=1)
    out = values[valid][idx]
    return out.reshape(len(target_lat), len(target_lon))


def _load_glodap(varfile, varname):
    path = os.path.join(GLODAP_ROOT, f"GLODAPv2.2016b.{varfile}.nc")
    ds = xr.open_dataset(path)
    da = _get_dataarray(ds, varname)

    depth_name = "depth_surface"
    if depth_name not in da.dims:
        for d in da.dims:
            if "depth" in d.lower():
                depth_name = d
                break

    lat_name = None
    lon_name = None
    for c in da.coords:
        if c.lower() in ("lat", "latitude"):
            lat_name = c
        if c.lower() in ("lon", "longitude"):
            lon_name = c

    new_lon = ((da[lon_name] + 180) % 360) - 180
    da = da.assign_coords({lon_name: new_lon}).sortby(lon_name)
    return da, lat_name, lon_name, depth_name


def _glodap_surface_map(da, depth_name, max_depth=50.0):
    depths = np.asarray(da[depth_name].values)
    sel = depths <= max_depth
    dsel = depths[sel]

    bounds = np.zeros(len(dsel) + 1)
    bounds[0] = 0.0
    if len(dsel) > 1:
        bounds[1:-1] = 0.5 * (dsel[:-1] + dsel[1:])
    bounds[-1] = max_depth
    dz = np.diff(bounds)

    da_s = da.isel({depth_name: np.where(sel)[0]})
    da_s = da_s.assign_coords({depth_name: dsel})
    dz_da = xr.DataArray(dz, dims=depth_name, coords={depth_name: dsel})
    return (da_s * dz_da).sum(depth_name) / dz.sum()


def _glodap_global_mean(da, depth_name, lat_name):
    vert = _glodap_surface_map(da, depth_name)
    vals = vert.values
    lat = da[lat_name].values
    w = np.cos(np.deg2rad(lat))[:, None]
    valid = np.isfinite(vals)
    num = np.nansum(vals * w * valid)
    den = np.nansum(w * valid)
    return float(num / den)


def _process_model_2002(model, glodap_dic_da, depth_name, target_lat, target_lon):
    ds_dic = _open_cmip(model, "dissic")
    ds_theta = _open_cmip(model, "thetao")
    ds_so = _open_cmip(model, "so")
    ds_talk = _open_cmip(model, "talk")
    ds_ph = _open_cmip(model, "ph")
    ds_vol = _open_cmip(model, "volcello")
    ds_area = _open_cmip(model, "areacello")

    def ann_mean(ds, var):
        da = _get_dataarray(ds, var)
        da = _select_years(da, 2002, 2002)
        if "time" in da.dims:
            da = da.mean("time")
        return da

    dic = ann_mean(ds_dic, "dissic")
    theta = ann_mean(ds_theta, "thetao")
    so = ann_mean(ds_so, "so")
    talk = ann_mean(ds_talk, "talk")
    ph = ann_mean(ds_ph, "ph")

    vol = _get_dataarray(ds_vol, "volcello")
    vol = _select_years(vol, 2002, 2002)
    if "time" in vol.dims:
        vol = vol.mean("time")

    area = _get_dataarray(ds_area, "areacello")
    if "time" in area.dims:
        area = area.mean("time")

    if model == "CESM2":
        ds_sftof = _open_cmip(model, "sftof")
        sftof = _get_dataarray(ds_sftof, "sftof")
        if "time" in sftof.dims:
            sftof = sftof.mean("time")
        ocean_mask = sftof > 0
        vol = vol.where(ocean_mask)
        area = area.where(ocean_mask)

    lev_m = _lev_m(dic, model)
    depth = lev_m.broadcast_like(theta)
    rho = _seawater_density(theta, so, depth / 10.0)

    dic_umolkg = dic * 1e6 / rho
    talk_umolkg = talk * 1e6 / rho

    vol_pos = vol.where(vol > 0)
    lev_mask = lev_m <= 50.0

    def wglobal(da):
        num = (da.where(lev_mask) * vol_pos.where(lev_mask)).sum()
        den = vol_pos.where(lev_mask).sum()
        return float((num / den).compute())

    dic_global = wglobal(dic_umolkg)
    alk_global = wglobal(talk_umolkg)
    ph_global = wglobal(ph)

    dic_map_da = (
        (dic_umolkg.where(lev_mask) * vol_pos.where(lev_mask)).sum("lev")
        / vol_pos.where(lev_mask).sum("lev")
    )
    dic_map = dic_map_da.compute().values
    lat2d, lon2d = _get_latlon2d(dic, model)

    dic_umolkg_np = dic_umolkg.compute().values
    lev_vals = np.asarray(lev_m.values)
    profile_bias = {}
    for zd in TARGET_DEPTHS:
        idx = int(np.argmin(np.abs(lev_vals - zd)))
        level_vals = dic_umolkg_np[idx, ...]
        model_on_glodap = _remap_to_target(
            lat2d, lon2d, level_vals, target_lat, target_lon
        )
        glodap_at_z = glodap_dic_da.sel({depth_name: zd}).values
        valid = np.isfinite(model_on_glodap) & np.isfinite(glodap_at_z)
        if valid.any():
            profile_bias[zd] = float(np.mean(model_on_glodap[valid] - glodap_at_z[valid]))
        else:
            profile_bias[zd] = np.nan

    return {
        "dic_global": dic_global,
        "alk_global": alk_global,
        "ph_global": ph_global,
        "dic_map": dic_map,
        "lat2d": lat2d,
        "lon2d": lon2d,
        "profile_bias": profile_bias,
    }


def _model_period_dic(model, y0, y1):
    ds_dic = _open_cmip(model, "dissic")
    ds_theta = _open_cmip(model, "thetao")
    ds_so = _open_cmip(model, "so")
    ds_vol = _open_cmip(model, "volcello")

    def mean_period(ds, var):
        da = _get_dataarray(ds, var)
        da = _select_years(da, y0, y1)
        if "time" in da.dims:
            da = da.mean("time")
        return da

    dic = mean_period(ds_dic, "dissic")
    theta = mean_period(ds_theta, "thetao")
    so = mean_period(ds_so, "so")

    vol = _get_dataarray(ds_vol, "volcello")
    vol = _select_years(vol, y0, y1)
    if "time" in vol.dims:
        vol = vol.mean("time")

    if model == "CESM2":
        ds_sftof = _open_cmip(model, "sftof")
        sftof = _get_dataarray(ds_sftof, "sftof")
        if "time" in sftof.dims:
            sftof = sftof.mean("time")
        vol = vol.where(sftof > 0)

    lev_m = _lev_m(dic, model)
    depth = lev_m.broadcast_like(theta)
    rho = _seawater_density(theta, so, depth / 10.0)
    dic_umolkg = dic * 1e6 / rho

    vol_pos = vol.where(vol > 0)
    lev_mask = lev_m <= 50.0

    num = (dic_umolkg.where(lev_mask) * vol_pos.where(lev_mask)).sum()
    den = vol_pos.where(lev_mask).sum()
    return float((num / den).compute())


# --- Load GLODAP mapped climatology and put it on a -180..180 sorted grid ---
glodap_dic_da, lat_name, lon_name, depth_name = _load_glodap("TCO2", "TCO2")
glodap_alk_da, _, _, _ = _load_glodap("TAlk", "TAlk")
glodap_ph_da, _, _, _ = _load_glodap("pHtsinsitutp", "pHtsinsitutp")

target_lat = np.asarray(glodap_dic_da[lat_name].values)
target_lon = np.asarray(glodap_dic_da[lon_name].values)

glodap_dic_map = _glodap_surface_map(glodap_dic_da, depth_name).values
glodap_alk_map = _glodap_surface_map(glodap_alk_da, depth_name).values
glodap_ph_map = _glodap_surface_map(glodap_ph_da, depth_name).values

glodap_dic_global = _glodap_global_mean(glodap_dic_da, depth_name, lat_name)
glodap_alk_global = _glodap_global_mean(glodap_alk_da, depth_name, lat_name)
glodap_ph_global = _glodap_global_mean(glodap_ph_da, depth_name, lat_name)

global_means = {
    "GLODAP": {
        "DIC": glodap_dic_global,
        "TAlk": glodap_alk_global,
        "pH": glodap_ph_global,
    }
}

bias_maps = {}
model_on_glodap_maps = {}
map_metrics = {}
profile_bias = {}
profile_models = {}

for model in MODELS:
    res = _process_model_2002(
        model, glodap_dic_da, depth_name, target_lat, target_lon
    )
    global_means[model] = {
        "DIC": res["dic_global"],
        "TAlk": res["alk_global"],
        "pH": res["ph_global"],
    }

    model_on_glodap = _remap_to_target(
        res["lat2d"], res["lon2d"], res["dic_map"], target_lat, target_lon
    )
    bias_map = model_on_glodap - glodap_dic_map
    valid = np.isfinite(bias_map)
    map_metrics[model] = {
        "mean_bias": float(np.mean(bias_map[valid])),
        "rmse": float(np.sqrt(np.mean(bias_map[valid] ** 2))),
        "corr": float(np.corrcoef(model_on_glodap[valid], glodap_dic_map[valid])[0, 1]),
    }

    bias_maps[model] = bias_map
    model_on_glodap_maps[model] = model_on_glodap
    profile_bias[model] = res["profile_bias"]
    profile_models[model] = res

dic_bias_0_50 = {
    m: global_means[m]["DIC"] - global_means["GLODAP"]["DIC"] for m in MODELS
}
alk_bias_0_50 = {
    m: global_means[m]["TAlk"] - global_means["GLODAP"]["TAlk"] for m in MODELS
}
ph_bias_0_50 = {
    m: global_means[m]["pH"] - global_means["GLODAP"]["pH"] for m in MODELS
}

period_change = {}
period_ratio = {}
dic_1971_2000 = {}
for model in MODELS:
    dic_1971_2000[model] = _model_period_dic(model, 1971, 2000)
    period_change[model] = global_means[model]["DIC"] - dic_1971_2000[model]
    period_ratio[model] = (
        period_change[model] / dic_bias_0_50[model]
        if abs(dic_bias_0_50[model]) > 1e-12
        else np.nan
    )

error_signs = {
    m: {
        "DIC": float(np.sign(dic_bias_0_50[m])),
        "TAlk": float(np.sign(alk_bias_0_50[m])),
        "pH": float(np.sign(ph_bias_0_50[m])),
    }
    for m in MODELS
}

sign_consistency = {
    "DIC": len({np.sign(dic_bias_0_50[m]) for m in MODELS}) == 1,
    "TAlk": len({np.sign(alk_bias_0_50[m]) for m in MODELS}) == 1,
    "pH": len({np.sign(ph_bias_0_50[m]) for m in MODELS}) == 1,
}

# A simple combined score: DIC map RMSE plus absolute surface biases, pH scaled by 1000.
model_scores = {
    m: (
        map_metrics[m]["rmse"]
        + abs(dic_bias_0_50[m])
        + abs(alk_bias_0_50[m])
        + 1000.0 * abs(ph_bias_0_50[m])
    )
    for m in MODELS
}
best_model = min(model_scores, key=model_scores.get)

# Useful summaries for the answer text.
summary_table = pd.DataFrame(
    {
        "DIC_0_50": {m: global_means[m]["DIC"] for m in ["GLODAP"] + MODELS},
        "TAlk_0_50": {m: global_means[m]["TAlk"] for m in ["GLODAP"] + MODELS},
        "pH_0_50": {m: global_means[m]["pH"] for m in ["GLODAP"] + MODELS},
        "DIC_bias": {m: dic_bias_0_50[m] for m in MODELS},
        "TAlk_bias": {m: alk_bias_0_50[m] for m in MODELS},
        "pH_bias": {m: ph_bias_0_50[m] for m in MODELS},
        "DIC_RMSE": {m: map_metrics[m]["rmse"] for m in MODELS},
        "DIC_corr": {m: map_metrics[m]["corr"] for m in MODELS},
        "period_change_DIC": {m: period_change[m] for m in MODELS},
        "period_ratio": {m: period_ratio[m] for m in MODELS},
    }
)

print(summary_table)
print("Best model by combined score:", best_model, model_scores)
print("Sign consistency:", sign_consistency)
print("Profile biases (umol/kg):", profile_bias)
