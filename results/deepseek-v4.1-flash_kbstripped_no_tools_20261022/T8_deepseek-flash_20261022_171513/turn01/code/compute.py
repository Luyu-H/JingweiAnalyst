# Computation code only
import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Paths and metadata
# ---------------------------------------------------------------------------
CMIP6_ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"
GLODAP_ROOT = "/data1/lorraine-20419/ocean_data/GLODAPv2.2016b_MappedClimatologies"

MODELS = {
    "CESM2": {
        "institution": "NCAR",
        "grid_type": "curvilinear",
        "grid_label": "gn",
        "y_dim": "nlat",
        "x_dim": "nlon",
        "shape": [384, 320],
        "lat_coord": "lat",
        "lon_coord": "lon",
        "lon_convention": "0_360",
        "depth_coord": "lev",
        "depth_units": "centimeters",
        "depth_scale": 0.01,
        "n_levels": 60,
        "calendar": "noleap",
        "ocean_mask_from": "sftof",
    },
    "GFDL-ESM4": {
        "institution": "NOAA-GFDL",
        "grid_type": "regular",
        "grid_label": "gr",
        "y_dim": "lat",
        "x_dim": "lon",
        "shape": [180, 360],
        "lat_coord": "lat",
        "lon_coord": "lon",
        "lon_convention": "0_360",
        "depth_coord": "lev",
        "depth_units": "m",
        "depth_scale": 1.0,
        "n_levels": 35,
        "calendar": "noleap",
    },
    "MPI-ESM1-2-LR": {
        "institution": "MPI-M",
        "grid_type": "curvilinear",
        "grid_label": "gn",
        "y_dim": "j",
        "x_dim": "i",
        "shape": [220, 256],
        "lat_coord": "latitude",
        "lon_coord": "longitude",
        "lon_convention": "0_360",
        "depth_coord": "lev",
        "depth_units": "m",
        "depth_scale": 1.0,
        "n_levels": 40,
        "calendar": "proleptic_gregorian",
    },
}

MODEL_LIST = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
DEPTHS = [0.0, 100.0, 500.0, 1000.0, 2000.0]

# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------
def find_cmip6_files(model, var):
    return sorted(glob.glob(os.path.join(CMIP6_ROOT, model, var, "*.nc")))

def open_files(files, chunk_time=False):
    if not files:
        return None
    chunks = {"time": 12} if chunk_time else None
    if len(files) == 1:
        return xr.open_dataset(
            files[0], decode_times=True, use_cftime=True, chunks=chunks
        )
    return xr.open_mfdataset(
        files,
        combine="by_coords",
        decode_times=True,
        use_cftime=True,
        chunks=chunks,
        data_vars="minimal",
        coords="minimal",
        compat="override",
    )

def get_dataarray(ds, var):
    if var in ds:
        return ds[var]
    for name in ds.data_vars:
        if name.lower() == var.lower():
            return ds[name]
    raise KeyError(f"{var!r} not found. Available data variables: {list(ds.data_vars)}")

def open_cmip6(model, var):
    files = find_cmip6_files(model, var)
    if not files:
        return None
    time_vars = {"dissic", "talk", "ph", "thetao", "so"}
    chunk_time = (var in time_vars) or any("Omon" in os.path.basename(f) for f in files)
    ds = open_files(files, chunk_time=chunk_time)
    return get_dataarray(ds, var)

def squeeze_singleton_time(da):
    if da is not None and "time" in da.dims and da.sizes["time"] == 1:
        return da.isel(time=0, drop=True)
    return da

def select_year(da, year):
    if da is None:
        return None
    if "time" not in da.dims:
        return da
    return da.where(da.time.dt.year == year, drop=True)

def select_period(da, y0, y1):
    if da is None:
        return None
    if "time" not in da.dims:
        return da
    return da.where((da.time.dt.year >= y0) & (da.time.dt.year <= y1), drop=True)

def depth_in_m(da, model_meta):
    lev_name = model_meta["depth_coord"]
    lev = da[lev_name]
    scale = float(model_meta.get("depth_scale", 1.0))
    return (lev * scale).assign_attrs(units="m")

def compute_potential_density_anomaly(theta, sal, depth_m, ref_pressure=None):
    """In-situ density from EOS-80, equivalent to ref_pressure=None here.

    theta and sal are potential temperature / practical salinity; depth_m is
    positive-down depth in metres.  The returned value is rho in kg m-3.
    """
    T = theta
    S = sal
    # pressure in bar: 10 m of seawater is about 1 bar
    p = depth_m / 10.0

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
    Kp = (
        K0
        + p
        * (
            3.239908
            + 1.43713e-3 * T
            + 1.16092e-4 * T**2
            - 5.77905e-7 * T**3
            + S * (2.2838e-3 - 1.0981e-5 * T - 1.6078e-6 * T**2)
            + S**1.5 * 1.91075e-4
        )
        + p**2
        * (
            8.50935e-5
            - 6.12293e-6 * T
            + 5.2787e-8 * T**2
            + S * (-9.9348e-7 + 2.0816e-8 * T + 9.1697e-10 * T**2)
        )
    )
    return rho0 / (1.0 - p / Kp)

def convert_mol_m3_to_umol_kg(da, rho):
    return (da * 1.0e6 / rho).assign_attrs(units="umol kg-1")

def volume_weighted_mean(field, weight, depth_da, max_depth=50.0):
    lev_name = depth_da.dims[0]
    sel = (depth_da >= 0.0) & (depth_da <= max_depth)
    f = field.where(sel)
    w = weight.where(sel)
    valid = np.isfinite(f) & np.isfinite(w) & (w > 0)
    return (f.where(valid) * w.where(valid)).sum() / w.where(valid).sum()

def level_thickness(depth_values, max_depth=None):
    depth = np.asarray(depth_values, dtype=float)
    n = len(depth)
    bounds = np.zeros(n + 1, dtype=float)
    bounds[1:-1] = 0.5 * (depth[:-1] + depth[1:])
    bounds[0] = depth[0] - 0.5 * (depth[1] - depth[0])
    bounds[-1] = depth[-1] + 0.5 * (depth[-1] - depth[-2])
    if max_depth is not None:
        bounds = np.clip(bounds, 0.0, max_depth)
    return np.diff(bounds)

def regular_cell_area(lat_values, lon_values):
    R = 6371000.0
    lat = np.asarray(lat_values, dtype=float)
    lon = np.asarray(lon_values, dtype=float)
    dlon = np.deg2rad(abs(np.mean(np.diff(lon)))) if len(lon) > 1 else np.deg2rad(1.0)

    lat_b = np.empty(len(lat) + 1, dtype=float)
    lat_b[1:-1] = 0.5 * (lat[:-1] + lat[1:])
    lat_b[0] = lat[0] - 0.5 * (lat[1] - lat[0])
    lat_b[-1] = lat[-1] + 0.5 * (lat[-1] - lat[-2])
    lat_b = np.clip(lat_b, -90.0, 90.0)

    area_lat = R**2 * dlon * (
        np.sin(np.deg2rad(lat_b[1:])) - np.sin(np.deg2rad(lat_b[:-1]))
    )
    return area_lat[:, None] * np.ones((1, len(lon)), dtype=float)

def latlon_2d(da, model_meta):
    lat_name = model_meta["lat_coord"]
    lon_name = model_meta["lon_coord"]
    lat = da[lat_name]
    lon = da[lon_name]
    if lat.ndim == 1 and lon.ndim == 1:
        ydim, xdim = model_meta["y_dim"], model_meta["x_dim"]
        lon2d, lat2d = np.meshgrid(lon.values, lat.values)
        lat_da = xr.DataArray(
            lat2d,
            dims=(ydim, xdim),
            coords={ydim: lat.values, xdim: lon.values},
            name=lat_name,
        )
        lon_da = xr.DataArray(
            lon2d,
            dims=(ydim, xdim),
            coords={ydim: lat.values, xdim: lon.values},
            name=lon_name,
        )
        return lat_da, lon_da
    return lat, lon

def regrid_to_glodap_da(da2d, lat2d, lon2d, area2d, glat, glon):
    """First-order conservative area-overlap regrid to the GLODAP 1-degree grid."""
    dims2d = da2d.dims
    if lat2d.dims != dims2d:
        lat2d = lat2d.transpose(*dims2d)
    if lon2d.dims != dims2d:
        lon2d = lon2d.transpose(*dims2d)
    if area2d.dims != dims2d:
        area2d = area2d.transpose(*dims2d)

    vals = np.asarray(da2d.values, dtype=float)
    latv = np.asarray(lat2d.values, dtype=float)
    lonv = np.asarray(lon2d.values, dtype=float)
    areav = np.asarray(area2d.values, dtype=float)

    valid = (
        np.isfinite(vals)
        & np.isfinite(latv)
        & np.isfinite(lonv)
        & np.isfinite(areav)
        & (areav > 0)
    )

    nlat = len(glat)
    nlon = len(glon)
    sumv = np.zeros((nlat, nlon), dtype=float)
    wsum = np.zeros((nlat, nlon), dtype=float)

    if valid.any():
        glat0 = float(glat.values[0])
        glon0 = float(glon.values[0])

        lat_idx = np.rint(latv[valid] - glat0).astype(int)
        lat_idx = np.clip(lat_idx, 0, nlat - 1)

        lon360 = np.mod(lonv[valid], 360.0)
        lon_idx = np.rint(np.mod(lon360 - glon0, 360.0)).astype(int) % nlon

        np.add.at(sumv, (lat_idx, lon_idx), vals[valid] * areav[valid])
        np.add.at(wsum, (lat_idx, lon_idx), areav[valid])

    with np.errstate(invalid="ignore", divide="ignore"):
        out = sumv / wsum
    out[wsum == 0] = np.nan

    return xr.DataArray(
        out,
        coords={glat.name: glat, glon.name: glon},
        dims=(glat.name, glon.name),
    )

def regrid_vertical(da, depth_da, target_depths):
    lev_name = depth_da.dims[0]
    da2 = da.assign_coords({lev_name: depth_da.values})
    da2 = da2.sortby(lev_name)

    first_depth = float(da2[lev_name][0])
    if first_depth > float(target_depths[0]):
        top = da2.isel({lev_name: 0}).expand_dims(
            {lev_name: [float(target_depths[0])]}
        )
        da2 = xr.concat([top, da2], dim=lev_name)

    out = da2.interp({lev_name: list(target_depths)}, method="linear")
    out = out.rename({lev_name: "depth"})
    out = out.assign_coords(depth=list(target_depths))
    return out

def build_comparison_mask(ref_map, model_maps, ocean_only=False):
    mask = np.isfinite(ref_map)
    for m in model_maps.values():
        mask = mask & np.isfinite(m)
    if ocean_only:
        mask = mask & (ref_map > 5.0)
    return mask

def compute_difference_metrics(model_map, ref_map, area, mask):
    valid = (
        mask
        & np.isfinite(model_map)
        & np.isfinite(ref_map)
        & np.isfinite(area)
        & (area > 0)
    )
    x = model_map.where(valid)
    y = ref_map.where(valid)
    w = area.where(valid)
    bias = x - y
    mean_bias = float((bias * w).sum().compute() / w.sum().compute())
    rmse = float(
        np.sqrt((bias**2 * w).sum().compute() / w.sum().compute())
    )
    return mean_bias, rmse

def analyze_correlation(model_map, ref_map, area, mask, weighting="area"):
    valid = (
        mask
        & np.isfinite(model_map)
        & np.isfinite(ref_map)
        & np.isfinite(area)
        & (area > 0)
    )
    x = model_map.where(valid)
    y = ref_map.where(valid)
    w = area.where(valid)
    mx = (x * w).sum() / w.sum()
    my = (y * w).sum() / w.sum()
    cov = ((x - mx) * (y - my) * w).sum() / w.sum()
    vx = ((x - mx) ** 2 * w).sum() / w.sum()
    vy = ((y - my) ** 2 * w).sum() / w.sum()
    return float((cov / np.sqrt(vx * vy)).compute())

# ---------------------------------------------------------------------------
# Load GLODAP
# ---------------------------------------------------------------------------
gds_tco2 = xr.open_dataset(
    os.path.join(GLODAP_ROOT, "GLODAPv2.2016b.TCO2.nc"), decode_times=False
)
gds_talk = xr.open_dataset(
    os.path.join(GLODAP_ROOT, "GLODAPv2.2016b.TAlk.nc"), decode_times=False
)
gds_ph = xr.open_dataset(
    os.path.join(GLODAP_ROOT, "GLODAPv2.2016b.pHtsinsitutp.nc"), decode_times=False
)

GLAT_NAME = "lat" if "lat" in gds_tco2.coords else "latitude"
GLON_NAME = "lon" if "lon" in gds_tco2.coords else "longitude"
GDEPTH_NAME = "depth_surface"

g_lat = gds_tco2[GLAT_NAME]
g_lon = gds_tco2[GLON_NAME]
g_depth_raw = (
    gds_tco2["Depth"] if "Depth" in gds_tco2.variables else gds_tco2[GDEPTH_NAME]
)
g_depth_da = xr.DataArray(
    np.asarray(g_depth_raw.values, dtype=float),
    coords={GDEPTH_NAME: np.asarray(g_depth_raw.values, dtype=float)},
    dims=(GDEPTH_NAME,),
)

def load_glodap_var(ds, var):
    da = get_dataarray(ds, var)
    da = da.transpose(GDEPTH_NAME, GLAT_NAME, GLON_NAME)
    da = da.where(np.isfinite(da))
    da = da.where(da < 1.0e10)
    return da.load()

glodap_dic = load_glodap_var(gds_tco2, "TCO2")
glodap_talk = load_glodap_var(gds_talk, "TAlk")
glodap_ph = load_glodap_var(gds_ph, "pHtsinsitutp")
glodap_dic, glodap_talk, glodap_ph = xr.align(
    glodap_dic, glodap_talk, glodap_ph, join="inner"
)

area_g_np = regular_cell_area(g_lat.values, g_lon.values)
area_g = xr.DataArray(
    area_g_np,
    coords={GLAT_NAME: g_lat, GLON_NAME: g_lon},
    dims=(GLAT_NAME, GLON_NAME),
)

thick_g_np = level_thickness(g_depth_da.values, max_depth=50.0)
thick_g = xr.DataArray(
    thick_g_np,
    coords={GDEPTH_NAME: g_depth_da.values},
    dims=(GDEPTH_NAME,),
)
vol_g = area_g * thick_g

glodap_dic_mean = float(
    volume_weighted_mean(glodap_dic, vol_g, g_depth_da, 50.0).compute()
)
glodap_talk_mean = float(
    volume_weighted_mean(glodap_talk, vol_g, g_depth_da, 50.0).compute()
)
glodap_ph_mean = float(
    volume_weighted_mean(glodap_ph, vol_g, g_depth_da, 50.0).compute()
)

glodap_dic_map = (
    glodap_dic.where(g_depth_da <= 50.0).mean(GDEPTH_NAME, skipna=True).compute()
)
glodap_dic_profile = regrid_vertical(glodap_dic, g_depth_da, DEPTHS)

# ---------------------------------------------------------------------------
# Process each CMIP6 model: 2002 annual mean and 1971-2000 mean
# ---------------------------------------------------------------------------
model_results = {}
model_dic_ann = {}
model_depth_m = {}
model_area = {}
model_dic_map_native = {}
model_dic_profile_bias = {}

density_assumption = (
    "in-situ density from EOS-80 with pressure estimated from depth "
    "(ref_pressure=None); no constant-density approximation"
)

for model in MODEL_LIST:
    meta = MODELS[model]
    lev_name = meta["depth_coord"]

    area_da = squeeze_singleton_time(open_cmip6(model, "areacello"))
    if area_da is None:
        raise FileNotFoundError(f"areacello not found for {model}")
    area_da = area_da.squeeze(drop=True)

    sftof_da = squeeze_singleton_time(open_cmip6(model, "sftof"))
    if sftof_da is not None:
        sftof_da = sftof_da.squeeze(drop=True)

    vol_da = open_cmip6(model, "volcello")
    if vol_da is None:
        raise FileNotFoundError(f"volcello not found for {model}")
    if "time" in vol_da.dims and vol_da.sizes["time"] == 1:
        vol_da = vol_da.isel(time=0, drop=True)

    # ---------------- 2002 monthly fields ----------------
    dic = select_year(open_cmip6(model, "dissic"), 2002)
    talk = select_year(open_cmip6(model, "talk"), 2002)
    ph = select_year(open_cmip6(model, "ph"), 2002)
    theta = select_year(open_cmip6(model, "thetao"), 2002)
    so = select_year(open_cmip6(model, "so"), 2002)

    if "time" in vol_da.dims:
        vol_2002 = select_year(vol_da, 2002)
        vol_2002, dic, talk, ph, theta, so = xr.align(
            vol_2002, dic, talk, ph, theta, so, join="inner"
        )
    else:
        vol_2002 = vol_da
        dic, talk, ph, theta, so = xr.align(
            dic, talk, ph, theta, so, join="inner"
        )

    depth_m = depth_in_m(dic, meta)
    rho = compute_potential_density_anomaly(theta, so, depth_m, ref_pressure=None)

    dic_c = convert_mol_m3_to_umol_kg(dic, rho)
    talk_c = convert_mol_m3_to_umol_kg(talk, rho)

    valid = (
        np.isfinite(dic)
        & np.isfinite(talk)
        & np.isfinite(ph)
        & np.isfinite(theta)
        & np.isfinite(so)
        & np.isfinite(rho)
        & (vol_2002 > 0)
    )
    if sftof_da is not None:
        valid = valid & (sftof_da > 0)

    dic_c = dic_c.where(valid)
    talk_c = talk_c.where(valid)
    ph_m = ph.where(valid)

    dic_ann = dic_c.mean("time", skipna=True).compute()
    talk_ann = talk_c.mean("time", skipna=True).compute()
    ph_ann = ph_m.mean("time", skipna=True).compute()

    if "time" in vol_2002.dims:
        vol_ann = (
            vol_2002.where(vol_2002 > 0).mean("time", skipna=True).compute()
        )
    else:
        vol_ann = vol_2002.where(vol_2002 > 0).compute()

    dic_mean_2002 = float(
        volume_weighted_mean(dic_ann, vol_ann, depth_m, 50.0).compute()
    )
    talk_mean_2002 = float(
        volume_weighted_mean(talk_ann, vol_ann, depth_m, 50.0).compute()
    )
    ph_mean_2002 = float(
        volume_weighted_mean(ph_ann, vol_ann, depth_m, 50.0).compute()
    )

    dic_map_native = (
        dic_ann.where(depth_m <= 50.0)
        .where(vol_ann > 0)
        .mean(lev_name, skipna=True)
        .compute()
    )

    # ---------------- 1971-2000 DIC ----------------
    dic_p = select_period(open_cmip6(model, "dissic"), 1971, 2000)
    theta_p = select_period(open_cmip6(model, "thetao"), 1971, 2000)
    so_p = select_period(open_cmip6(model, "so"), 1971, 2000)

    if "time" in vol_da.dims:
        vol_p = select_period(vol_da, 1971, 2000)
        vol_p, dic_p, theta_p, so_p = xr.align(
            vol_p, dic_p, theta_p, so_p, join="inner"
        )
    else:
        vol_p = vol_da
        dic_p, theta_p, so_p = xr.align(
            dic_p, theta_p, so_p, join="inner"
        )

    depth_m_p = depth_in_m(dic_p, meta)
    rho_p = compute_potential_density_anomaly(
        theta_p, so_p, depth_m_p, ref_pressure=None
    )
    dic_p_c = convert_mol_m3_to_umol_kg(dic_p, rho_p)

    valid_p = (
        np.isfinite(dic_p)
        & np.isfinite(theta_p)
        & np.isfinite(so_p)
        & np.isfinite(rho_p)
        & (vol_p > 0)
    )
    if sftof_da is not None:
        valid_p = valid_p & (sftof_da > 0)

    dic_p_c = dic_p_c.where(valid_p)
    dic_p_ann = dic_p_c.mean("time", skipna=True).compute()

    if "time" in vol_p.dims:
        vol_p_ann = (
            vol_p.where(vol_p > 0).mean("time", skipna=True).compute()
        )
    else:
        vol_p_ann = vol_p.where(vol_p > 0).compute()

    dic_p_mean = float(
        volume_weighted_mean(dic_p_ann, vol_p_ann, depth_m_p, 50.0).compute()
    )

    model_results[model] = {
        "DIC_mean": dic_mean_2002,
        "TAlk_mean": talk_mean_2002,
        "pH_mean": ph_mean_2002,
        "DIC_mean_1971_2000": dic_p_mean,
        "DIC_change_2002_minus_1971_2000": dic_mean_2002 - dic_p_mean,
    }
    model_dic_ann[model] = dic_ann
    model_depth_m[model] = depth_m
    model_area[model] = area_da
    model_dic_map_native[model] = dic_map_native

# ---------------------------------------------------------------------------
# Regrid model 0-50 m DIC maps to the GLODAP grid and compare
# ---------------------------------------------------------------------------
model_dic_map_glodap = {}
for model in MODEL_LIST:
    meta = MODELS[model]
    da_map = model_dic_map_native[model]
    area_da = model_area[model]
    lat2d, lon2d = latlon_2d(da_map, meta)
    model_dic_map_glodap[model] = regrid_to_glodap_da(
        da_map, lat2d, lon2d, area_da, g_lat, g_lon
    )

comparison_mask = build_comparison_mask(
    glodap_dic_map, model_dic_map_glodap, ocean_only=False
)
comparison_mask_count = int(comparison_mask.sum().compute())

glodap_dic_map_masked = glodap_dic_map.where(comparison_mask)
model_dic_map_masked = {
    model: model_dic_map_glodap[model].where(comparison_mask)
    for model in MODEL_LIST
}
dic_bias_maps = {
    model: (model_dic_map_masked[model] - glodap_dic_map_masked).compute()
    for model in MODEL_LIST
}

map_metrics = {}
for model in MODEL_LIST:
    mean_bias, rmse = compute_difference_metrics(
        model_dic_map_masked[model],
        glodap_dic_map_masked,
        area_g,
        comparison_mask,
    )
    corr = analyze_correlation(
        model_dic_map_masked[model],
        glodap_dic_map_masked,
        area_g,
        comparison_mask,
        weighting="area",
    )
    map_metrics[model] = {
        "DIC_map_bias": mean_bias,
        "DIC_map_rmse": rmse,
        "DIC_map_corr": corr,
    }

# ---------------------------------------------------------------------------
# Vertical DIC-bias profile at 0, 100, 500, 1000, 2000 m
# ---------------------------------------------------------------------------
for model in MODEL_LIST:
    meta = MODELS[model]
    dic_ann = model_dic_ann[model]
    depth_m = model_depth_m[model]
    area_da = model_area[model]
    lat2d, lon2d = latlon_2d(dic_ann, meta)

    model_prof = regrid_vertical(dic_ann, depth_m, DEPTHS).compute()

    reg_prof_levels = []
    for d in DEPTHS:
        da2d = model_prof.sel(depth=d)
        reg_prof_levels.append(
            regrid_to_glodap_da(da2d, lat2d, lon2d, area_da, g_lat, g_lon)
        )
    model_prof_reg = xr.concat(
        reg_prof_levels, dim=xr.IndexVariable("depth", DEPTHS)
    )

    bias_at_depth = []
    for d in DEPTHS:
        m2 = model_prof_reg.sel(depth=d)
        g2 = glodap_dic_profile.sel(depth=d)
        valid = (
            np.isfinite(m2)
            & np.isfinite(g2)
            & np.isfinite(area_g)
            & (area_g > 0)
        )
        bias = (m2 - g2).where(valid)
        w = area_g.where(valid)
        bias_at_depth.append(
            float((bias * w).sum().compute() / w.sum().compute())
        )
    model_dic_profile_bias[model] = bias_at_depth

# ---------------------------------------------------------------------------
# Period change, ratio to GLODAP DIC bias, summary tables and best model
# ---------------------------------------------------------------------------
period_results = {}
for model in MODEL_LIST:
    dic2002 = model_results[model]["DIC_mean"]
    dicperiod = model_results[model]["DIC_mean_1971_2000"]
    change = dic2002 - dicperiod
    glodap_bias = dic2002 - glodap_dic_mean
    ratio = change / glodap_bias if glodap_bias != 0 else np.nan
    period_results[model] = {
        "DIC_change": change,
        "DIC_ratio_to_GLODAP_bias": ratio,
        "DIC_2002_minus_GLODAP": glodap_bias,
    }

summary_rows = [
    {
        "model": "GLODAP",
        "DIC_mean": glodap_dic_mean,
        "DIC_minus_GLODAP": 0.0,
        "TAlk_mean": glodap_talk_mean,
        "TAlk_minus_GLODAP": 0.0,
        "pH_mean": glodap_ph_mean,
        "pH_minus_GLODAP": 0.0,
    }
]
for model in MODEL_LIST:
    summary_rows.append(
        {
            "model": model,
            "DIC_mean": model_results[model]["DIC_mean"],
            "DIC_minus_GLODAP": model_results[model]["DIC_mean"]
            - glodap_dic_mean,
            "TAlk_mean": model_results[model]["TAlk_mean"],
            "TAlk_minus_GLODAP": model_results[model]["TAlk_mean"]
            - glodap_talk_mean,
            "pH_mean": model_results[model]["pH_mean"],
            "pH_minus_GLODAP": model_results[model]["pH_mean"]
            - glodap_ph_mean,
        }
    )
summary_df = pd.DataFrame(summary_rows)

map_rows = []
for model in MODEL_LIST:
    map_rows.append(
        {
            "model": model,
            "DIC_map_bias": map_metrics[model]["DIC_map_bias"],
            "DIC_map_rmse": map_metrics[model]["DIC_map_rmse"],
            "DIC_map_corr": map_metrics[model]["DIC_map_corr"],
        }
    )
map_metrics_df = pd.DataFrame(map_rows)

profile_rows = []
for model in MODEL_LIST:
    row = {"model": model}
    for d, b in zip(DEPTHS, model_dic_profile_bias[model]):
        row[f"DIC_bias_{int(d)}m"] = b
    profile_rows.append(row)
profile_df = pd.DataFrame(profile_rows)

period_rows = []
for model in MODEL_LIST:
    period_rows.append(
        {
            "model": model,
            "DIC_change_2002_minus_1971_2000": period_results[model][
                "DIC_change"
            ],
            "ratio_to_GLODAP_DIC_bias": period_results[model][
                "DIC_ratio_to_GLODAP_bias"
            ],
        }
    )
period_df = pd.DataFrame(period_rows)

# Sign consistency
signs = {}
for model in MODEL_LIST:
    signs[model] = {
        "DIC": int(np.sign(model_results[model]["DIC_mean"] - glodap_dic_mean)),
        "TAlk": int(np.sign(model_results[model]["TAlk_mean"] - glodap_talk_mean)),
        "pH": int(np.sign(model_results[model]["pH_mean"] - glodap_ph_mean)),
    }

sign_across_models = {
    var: len({signs[m][var] for m in MODEL_LIST}) == 1
    for var in ["DIC", "TAlk", "pH"]
}
sign_within_model = {
    model: len(set(signs[model].values())) == 1 for model in MODEL_LIST
}
sign_consistency = {
    "signs": signs,
    "across_models": sign_across_models,
    "within_model": sign_within_model,
}

# Best-model score: normalized absolute mean biases, map RMSE, map correlation,
# and mean absolute profile bias.
dic_abs = np.array(
    [abs(model_results[m]["DIC_mean"] - glodap_dic_mean) for m in MODEL_LIST]
)
talk_abs = np.array(
    [abs(model_results[m]["TAlk_mean"] - glodap_talk_mean) for m in MODEL_LIST]
)
ph_abs = np.array(
    [abs(model_results[m]["pH_mean"] - glodap_ph_mean) for m in MODEL_LIST]
)
rmse_arr = np.array([map_metrics[m]["DIC_map_rmse"] for m in MODEL_LIST])
corr_arr = np.array([map_metrics[m]["DIC_map_corr"] for m in MODEL_LIST])
profile_abs = np.array(
    [np.mean(np.abs(model_dic_profile_bias[m])) for m in MODEL_LIST]
)

def normalize_score(a):
    a = np.asarray(a, dtype=float)
    positive = a[np.isfinite(a)]
    scale = np.mean(np.abs(positive)) if positive.size else 1.0
    return a / scale if scale > 0 else a

best_model_score_array = (
    normalize_score(dic_abs)
    + normalize_score(talk_abs)
    + normalize_score(ph_abs)
    + normalize_score(rmse_arr)
    + (1.0 - corr_arr)
    + normalize_score(profile_abs)
)
best_model = MODEL_LIST[int(np.nanargmin(best_model_score_array))]
best_model_scores = {
    model: float(best_model_score_array[i])
    for i, model in enumerate(MODEL_LIST)
}

# Values useful for plotting
g_lat_vals = g_lat.values
g_lon_vals = g_lon.values
depth_list = DEPTHS

# Rounded reporting table
report_table = summary_df.copy()
for col in ["DIC_mean", "DIC_minus_GLODAP", "TAlk_mean", "TAlk_minus_GLODAP"]:
    report_table[col] = report_table[col].round(0)
for col in ["pH_mean", "pH_minus_GLODAP"]:
    report_table[col] = report_table[col].round(3)

map_report_table = map_metrics_df.copy()
map_report_table["DIC_map_bias"] = map_report_table["DIC_map_bias"].round(0)
map_report_table["DIC_map_rmse"] = map_report_table["DIC_map_rmse"].round(0)
map_report_table["DIC_map_corr"] = map_report_table["DIC_map_corr"].round(3)

profile_report_table = profile_df.copy()
for col in profile_report_table.columns:
    if col.startswith("DIC_bias_"):
        profile_report_table[col] = profile_report_table[col].round(0)

period_report_table = period_df.copy()
period_report_table["DIC_change_2002_minus_1971_2000"] = period_report_table[
    "DIC_change_2002_minus_1971_2000"
].round(0)
period_report_table["ratio_to_GLODAP_DIC_bias"] = period_report_table[
    "ratio_to_GLODAP_DIC_bias"
].round(2)
