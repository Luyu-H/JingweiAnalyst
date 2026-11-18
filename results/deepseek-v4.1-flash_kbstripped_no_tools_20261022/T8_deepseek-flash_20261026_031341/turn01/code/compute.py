import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
CMIP6_ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"
GLODAP_ROOT = "/data1/lorraine-20419/ocean_data/GLODAPv2.2016b_MappedClimatologies"

models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
model_meta = {
    "CESM2": {
        "lat": "lat", "lon": "lon", "depth": "lev", "depth_scale": 0.01,
        "grid": "curvilinear",
    },
    "GFDL-ESM4": {
        "lat": "lat", "lon": "lon", "depth": "lev", "depth_scale": 1.0,
        "grid": "regular",
    },
    "MPI-ESM1-2-LR": {
        "lat": "latitude", "lon": "longitude", "depth": "lev", "depth_scale": 1.0,
        "grid": "curvilinear",
    },
}

R_EARTH = 6371000.0
depth_targets = np.array([0, 100, 500, 1000, 2000], dtype=float)

# ---------------------------------------------------------------------------
# File discovery / opening helpers
# ---------------------------------------------------------------------------
def find_cmip6_files(model, var):
    patterns = []
    if var in ["areacello", "sftof"]:
        patterns.append(f"{var}_Ofx_{model}_historical_*_*.nc")
    elif var == "volcello":
        patterns.append(f"{var}_Ofx_{model}_historical_*_*.nc")
        patterns.append(f"{var}_Omon_{model}_historical_*_*.nc")
    else:
        patterns.append(f"{var}_Omon_{model}_historical_*_*.nc")

    files = []
    for pat in patterns:
        files.extend(glob.glob(os.path.join(CMIP6_ROOT, model, var, pat)))

    if not files:
        dirpath = os.path.join(CMIP6_ROOT, model, var)
        if os.path.isdir(dirpath):
            files = glob.glob(os.path.join(dirpath, "*.nc"))
    return sorted(set(files))


def open_cmip6_var(model, var):
    files = find_cmip6_files(model, var)
    if not files:
        return None, None

    fname = os.path.basename(files[0])
    parts = fname.split("_")
    variant = parts[4] if len(parts) > 4 else None

    kwargs = dict(use_cftime=True, decode_times=True)
    if len(files) == 1:
        ds = xr.open_dataset(files[0], **kwargs)
    else:
        if var in ["dissic", "talk", "ph", "thetao", "so"] or (
            var == "volcello" and "Omon" in files[0]
        ):
            ds = xr.open_mfdataset(files, combine="by_coords", **kwargs)
        else:
            ds = xr.open_dataset(files[0], **kwargs)
    return ds, variant


def load_glodap(varfile, varname):
    path = os.path.join(GLODAP_ROOT, f"GLODAPv2.2016b.{varfile}.nc")
    ds = xr.open_dataset(path)
    da = ds[varname]
    if "Depth" in ds:
        da = da.assign_coords(depth_surface=ds["Depth"])
    return da, ds


# ---------------------------------------------------------------------------
# Physical helper functions
# ---------------------------------------------------------------------------
def compute_density(T, S, depth_m):
    """
    In-situ density (kg m-3) from potential temperature T (degC), practical
    salinity S, and depth (m), using the UNESCO 1983 (EOS-80) equation of
    state.  Pressure is approximated from depth as P[dbar] ~= depth[m],
    then converted to bars for the EOS-80 coefficients.
    """
    T = T.astype(float)
    S = S.astype(float)
    P = (depth_m / 10.0).astype(float)  # dbar -> bar

    T_b, P_b = xr.broadcast(T, P)

    rho0 = (
        999.842594
        + 6.793952e-2 * T_b
        - 9.095290e-3 * T_b**2
        + 1.001685e-4 * T_b**3
        - 1.120083e-6 * T_b**4
        + 6.536332e-9 * T_b**5
        + (
            0.824493
            - 4.0899e-3 * T_b
            + 7.6438e-5 * T_b**2
            - 8.2467e-7 * T_b**3
            + 5.3875e-9 * T_b**4
        )
        * S
        + (-5.72466e-3 + 1.0227e-4 * T_b - 1.6546e-6 * T_b**2) * S**1.5
        + 4.8314e-4 * S**2
    )

    K0 = (
        19652.21
        + 148.4206 * T_b
        - 2.327105 * T_b**2
        + 1.360477e-2 * T_b**3
        - 5.155288e-5 * T_b**4
        + (
            54.6746
            - 0.603459 * T_b
            + 1.09987e-2 * T_b**2
            - 6.1670e-5 * T_b**3
        )
        * S
        + (7.944e-2 + 1.6483e-2 * T_b - 5.3009e-4 * T_b**2) * S**1.5
    )

    K1 = (
        3.239908
        + 1.43713e-3 * T_b
        + 1.16092e-4 * T_b**2
        - 5.77905e-7 * T_b**3
        + (
            2.2838e-3
            - 1.0981e-5 * T_b
            - 1.6078e-6 * T_b**2
        )
        * S
        + 1.91075e-4 * S**1.5
    )

    K2 = (
        8.50935e-5
        - 6.12293e-6 * T_b
        + 5.2787e-8 * T_b**2
        + (-9.9348e-7 + 2.0816e-8 * T_b + 9.1697e-10 * T_b**2) * S
    )

    K = K0 + P_b * K1 + P_b**2 * K2
    rho = rho0 / (1.0 - P_b / K)
    return rho


def regrid_vertical(da, depth_name, target_depths):
    da = da.compute()
    return da.interp(
        {depth_name: target_depths},
        kwargs={"fill_value": "extrapolate"},
    )


def regrid_regular_to_glodap(da, lat_name, lon_name, target_lat, target_lon):
    da = da.compute()
    lon = da[lon_name]
    lon_shift = lon % 360.0
    lon_shift = xr.where(lon_shift < 20.5, lon_shift + 360.0, lon_shift)
    da = da.assign_coords({lon_name: lon_shift}).sortby(lon_name)
    out = da.interp(
        {lat_name: target_lat, lon_name: target_lon},
        method="linear",
    )
    out = out.rename({lat_name: "lat", lon_name: "lon"})
    return out


def regrid_curvilinear_to_glodap(da, lat2d, lon2d, target_lat, target_lon):
    da = da.compute()
    lat = lat2d.values.ravel()
    lon = lon2d.values.ravel() % 360.0
    lon = np.where(lon < 20.5, lon + 360.0, lon)
    vals = da.values.ravel()

    valid = np.isfinite(vals) & np.isfinite(lat) & np.isfinite(lon)
    points = np.column_stack([lat[valid], lon[valid]])
    values = vals[valid]

    target_lat_mesh, target_lon_mesh = np.meshgrid(
        target_lat, target_lon, indexing="ij"
    )
    grid = griddata(
        points,
        values,
        (target_lat_mesh, target_lon_mesh),
        method="linear",
    )
    return xr.DataArray(
        grid,
        coords={"lat": target_lat, "lon": target_lon},
        dims=("lat", "lon"),
    )


def area_weighted_mean_0_50(da, depth_name, area, depth_max=50.0):
    da = da.compute()
    area = area.compute()
    da_sub = da.sel({depth_name: da[depth_name] <= depth_max})
    da_sub_b, area_b = xr.broadcast(da_sub, area)
    mask = np.isfinite(da_sub_b) & np.isfinite(area_b)
    prod = (da_sub_b * area_b).where(mask)
    num = float(prod.sum().compute())
    den = float(area_b.where(mask).sum().compute())
    return num / den


# ---------------------------------------------------------------------------
# Load GLODAP climatology
# ---------------------------------------------------------------------------
tco2_glodap, ds_tco2 = load_glodap("TCO2", "TCO2")
talk_glodap, ds_talk = load_glodap("TAlk", "TAlk")
ph_glodap, ds_ph = load_glodap("pHtsinsitutp", "pHtsinsitutp")
sal_glodap, _ = load_glodap("salinity", "salinity")
temp_glodap, _ = load_glodap("temperature", "temperature")

glodap_lat = tco2_glodap["lat"].values
glodap_lon = tco2_glodap["lon"].values
glodap_depth = tco2_glodap["depth_surface"].values

dlat = np.deg2rad(1.0)
dlon = np.deg2rad(1.0)
glodap_area = xr.DataArray(
    (R_EARTH**2)
    * np.cos(np.deg2rad(glodap_lat))[:, None]
    * dlat
    * dlon
    * np.ones((len(glodap_lat), len(glodap_lon))),
    coords={"lat": glodap_lat, "lon": glodap_lon},
    dims=("lat", "lon"),
)

tco2_glodap_mean = area_weighted_mean_0_50(
    tco2_glodap, "depth_surface", glodap_area
)
talk_glodap_mean = area_weighted_mean_0_50(
    talk_glodap, "depth_surface", glodap_area
)
ph_glodap_mean = area_weighted_mean_0_50(
    ph_glodap, "depth_surface", glodap_area
)

tco2_glodap_2d = (
    tco2_glodap.sel(depth_surface=tco2_glodap.depth_surface <= 50.0)
    .mean("depth_surface")
    .compute()
)

tco2_vert_glodap = regrid_vertical(
    tco2_glodap, "depth_surface", depth_targets
)

# ---------------------------------------------------------------------------
# Process each CMIP6 model
# ---------------------------------------------------------------------------
model_results = {}
missing_report = {}

for model in models:
    meta = model_meta[model]
    lat_name = meta["lat"]
    lon_name = meta["lon"]
    depth_name = meta["depth"]
    depth_scale = meta["depth_scale"]
    grid_type = meta["grid"]

    missing_report[model] = []

    # Load monthly Omon fields
    ds_vars = {}
    for var in ["dissic", "talk", "ph", "thetao", "so"]:
        ds, variant = open_cmip6_var(model, var)
        if ds is None:
            missing_report[model].append(var)
        else:
            ds_vars[var] = ds

    # Load areacello
    ds_area, _ = open_cmip6_var(model, "areacello")
    if ds_area is None:
        missing_report[model].append("areacello")
        model_results[model] = {
            "dic_0_50": np.nan,
            "talk_0_50": np.nan,
            "ph_0_50": np.nan,
            "dic_7100_0_50": np.nan,
            "dic_2d_glodap": None,
            "depth_biases": np.full(len(depth_targets), np.nan),
            "depth_counts": np.zeros(len(depth_targets), dtype=int),
        }
        continue

    areacello = ds_area["areacello"]

    # CESM2 raw areacello covers the whole globe; mask by sftof
    if model == "CESM2":
        ds_sftof, _ = open_cmip6_var(model, "sftof")
        if ds_sftof is not None:
            sftof = ds_sftof["sftof"]
            areacello = areacello.where(sftof > 0)

    core_missing = any(
        v in missing_report[model] for v in ["dissic", "talk", "ph", "thetao", "so"]
    )
    if core_missing:
        model_results[model] = {
            "dic_0_50": np.nan,
            "talk_0_50": np.nan,
            "ph_0_50": np.nan,
            "dic_7100_0_50": np.nan,
            "dic_2d_glodap": None,
            "depth_biases": np.full(len(depth_targets), np.nan),
            "depth_counts": np.zeros(len(depth_targets), dtype=int),
        }
        continue

    # Convert depth coordinate to metres
    for var, ds in ds_vars.items():
        if depth_name in ds.coords:
            ds = ds.assign_coords(
                {depth_name: ds[depth_name] * depth_scale}
            )
            ds_vars[var] = ds

    def select_years(ds, var, years):
        time = ds["time"]
        if years == 2002:
            mask = time.dt.year == 2002
        else:
            mask = (time.dt.year >= 1971) & (time.dt.year <= 2000)
        return ds[var].isel(time=mask)

    dissic_2002 = select_years(ds_vars["dissic"], "dissic", 2002)
    talk_2002 = select_years(ds_vars["talk"], "talk", 2002)
    ph_2002 = select_years(ds_vars["ph"], "ph", 2002)
    thetao_2002 = select_years(ds_vars["thetao"], "thetao", 2002)
    so_2002 = select_years(ds_vars["so"], "so", 2002)

    dissic_7100 = select_years(ds_vars["dissic"], "dissic", 7100)
    thetao_7100 = select_years(ds_vars["thetao"], "thetao", 7100)
    so_7100 = select_years(ds_vars["so"], "so", 7100)

    # Align time axes
    (
        dissic_2002,
        talk_2002,
        ph_2002,
        thetao_2002,
        so_2002,
    ) = xr.align(
        dissic_2002,
        talk_2002,
        ph_2002,
        thetao_2002,
        so_2002,
        join="inner",
    )
    dissic_7100, thetao_7100, so_7100 = xr.align(
        dissic_7100,
        thetao_7100,
        so_7100,
        join="inner",
    )

    # Monthly in-situ density and per-mass concentrations
    lev_2002 = dissic_2002[depth_name]
    rho_2002 = compute_density(thetao_2002, so_2002, lev_2002)
    dissic_2002_umol = dissic_2002 * 1.0e6 / rho_2002
    talk_2002_umol = talk_2002 * 1.0e6 / rho_2002

    lev_7100 = dissic_7100[depth_name]
    rho_7100 = compute_density(thetao_7100, so_7100, lev_7100)
    dissic_7100_umol = dissic_7100 * 1.0e6 / rho_7100

    def temporal_mean(da):
        w = da.time.dt.days_in_month.astype(float)
        return da.weighted(w).mean("time")

    dissic_2002_mean = temporal_mean(dissic_2002_umol)
    talk_2002_mean = temporal_mean(talk_2002_umol)
    ph_2002_mean = temporal_mean(ph_2002)
    dissic_7100_mean = temporal_mean(dissic_7100_umol)

    # Native-grid 0-50 m global means
    dic_0_50 = area_weighted_mean_0_50(
        dissic_2002_mean, depth_name, areacello
    )
    talk_0_50 = area_weighted_mean_0_50(
        talk_2002_mean, depth_name, areacello
    )
    ph_0_50 = area_weighted_mean_0_50(
        ph_2002_mean, depth_name, areacello
    )
    dic_7100_0_50 = area_weighted_mean_0_50(
        dissic_7100_mean, depth_name, areacello
    )

    # 2-D 0-50 m DIC for map comparison
    dic_2d = (
        dissic_2002_mean.sel(
            {depth_name: dissic_2002_mean[depth_name] <= 50.0}
        )
        .mean(depth_name)
        .compute()
    )

    if grid_type == "regular":
        dic_2d_glodap = regrid_regular_to_glodap(
            dic_2d, lat_name, lon_name, glodap_lat, glodap_lon
        )
    else:
        lat2d = ds_vars["dissic"][lat_name]
        lon2d = ds_vars["dissic"][lon_name]
        dic_2d_glodap = regrid_curvilinear_to_glodap(
            dic_2d, lat2d, lon2d, glodap_lat, glodap_lon
        )

    # Fixed-depth DIC biases
    dic_vert = regrid_vertical(dissic_2002_mean, depth_name, depth_targets)

    depth_biases = []
    depth_counts = []
    for i, d in enumerate(depth_targets):
        model_2d = dic_vert.isel({depth_name: i})
        glodap_2d = tco2_vert_glodap.isel(depth_surface=i)

        if grid_type == "regular":
            model_2d_glodap = regrid_regular_to_glodap(
                model_2d, lat_name, lon_name, glodap_lat, glodap_lon
            )
        else:
            model_2d_glodap = regrid_curvilinear_to_glodap(
                model_2d, lat2d, lon2d, glodap_lat, glodap_lon
            )

        diff = model_2d_glodap - glodap_2d
        mask = np.isfinite(diff) & np.isfinite(glodap_area)
        w = glodap_area.where(mask)
        if int(mask.sum().compute()) > 0:
            bias = float((diff * w).sum().compute() / w.sum().compute())
            count = int(mask.sum().compute())
        else:
            bias = np.nan
            count = 0
        depth_biases.append(bias)
        depth_counts.append(count)

    model_results[model] = {
        "dic_0_50": dic_0_50,
        "talk_0_50": talk_0_50,
        "ph_0_50": ph_0_50,
        "dic_7100_0_50": dic_7100_0_50,
        "dic_2d_glodap": dic_2d_glodap,
        "depth_biases": np.array(depth_biases),
        "depth_counts": np.array(depth_counts),
    }

# ---------------------------------------------------------------------------
# Common-mask map comparison of 0-50 m DIC
# ---------------------------------------------------------------------------
fields = [tco2_glodap_2d] + [model_results[m]["dic_2d_glodap"] for m in models]
fields_aligned = xr.align(*fields, join="inner")

mask_map = xr.DataArray(
    np.ones_like(fields_aligned[0], dtype=bool),
    coords=fields_aligned[0].coords,
    dims=fields_aligned[0].dims,
)
for f in fields_aligned:
    mask_map = mask_map & np.isfinite(f)

masked_fields = [f.where(mask_map) for f in fields_aligned]
dic_glodap_masked = masked_fields[0]
model_masked = {m: masked_fields[i + 1] for i, m in enumerate(models)}

map_metrics = {}
for m in models:
    diff = model_masked[m] - dic_glodap_masked
    w = glodap_area.where(mask_map & np.isfinite(diff))
    bias = float((diff * w).sum().compute() / w.sum().compute())
    rmse = float(
        ((diff**2) * w).sum().compute() / w.sum().compute()
    ) ** 0.5

    x = model_masked[m].values
    y = dic_glodap_masked.values
    ww = glodap_area.where(mask_map).values
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(ww)
    xv = x[valid]
    yv = y[valid]
    wv = ww[valid]
    xm = np.sum(xv * wv) / np.sum(wv)
    ym = np.sum(yv * wv) / np.sum(wv)
    cov = np.sum((xv - xm) * (yv - ym) * wv) / np.sum(wv)
    vx = np.sum((xv - xm) ** 2 * wv) / np.sum(wv)
    vy = np.sum((yv - ym) ** 2 * wv) / np.sum(wv)
    corr = cov / np.sqrt(vx * vy)

    map_metrics[m] = {
        "bias": bias,
        "rmse": rmse,
        "corr": corr,
        "diff": diff,
    }

# ---------------------------------------------------------------------------
# Period change and ratio to GLODAP offset
# ---------------------------------------------------------------------------
for m in models:
    change = model_results[m]["dic_0_50"] - model_results[m]["dic_7100_0_50"]
    glodap_diff = model_results[m]["dic_0_50"] - tco2_glodap_mean
    ratio = change / glodap_diff if glodap_diff != 0 else np.nan
    model_results[m]["dic_change"] = change
    model_results[m]["glodap_diff"] = glodap_diff
    model_results[m]["ratio"] = ratio

# ---------------------------------------------------------------------------
# Summary tables (printed for the answer)
# ---------------------------------------------------------------------------
dic_talk_ph_summary = pd.DataFrame(
    {
        "model": models,
        "DIC_0_50": [model_results[m]["dic_0_50"] for m in models],
        "DIC_diff": [
            model_results[m]["dic_0_50"] - tco2_glodap_mean for m in models
        ],
        "TAlk_0_50": [model_results[m]["talk_0_50"] for m in models],
        "TAlk_diff": [
            model_results[m]["talk_0_50"] - talk_glodap_mean for m in models
        ],
        "pH_0_50": [model_results[m]["ph_0_50"] for m in models],
        "pH_diff": [
            model_results[m]["ph_0_50"] - ph_glodap_mean for m in models
        ],
    }
)
dic_talk_ph_summary["DIC_0_50"] = dic_talk_ph_summary["DIC_0_50"].round(0)
dic_talk_ph_summary["DIC_diff"] = dic_talk_ph_summary["DIC_diff"].round(0)
dic_talk_ph_summary["TAlk_0_50"] = dic_talk_ph_summary["TAlk_0_50"].round(0)
dic_talk_ph_summary["TAlk_diff"] = dic_talk_ph_summary["TAlk_diff"].round(0)
dic_talk_ph_summary["pH_0_50"] = dic_talk_ph_summary["pH_0_50"].round(3)
dic_talk_ph_summary["pH_diff"] = dic_talk_ph_summary["pH_diff"].round(3)

print("GLODAP 0-50 m means:")
print(f"  TCO2  = {tco2_glodap_mean:.0f} umol kg-1")
print(f"  TAlk  = {talk_glodap_mean:.0f} umol kg-1")
print(f"  pH    = {ph_glodap_mean:.3f}")
print()
print("Model 0-50 m means and differences from GLODAP:")
print(dic_talk_ph_summary.to_string(index=False))
print()

map_summary = pd.DataFrame(
    {
        "model": models,
        "bias": [map_metrics[m]["bias"] for m in models],
        "rmse": [map_metrics[m]["rmse"] for m in models],
        "corr": [map_metrics[m]["corr"] for m in models],
    }
)
map_summary["bias"] = map_summary["bias"].round(0)
map_summary["rmse"] = map_summary["rmse"].round(0)
map_summary["corr"] = map_summary["corr"].round(3)
print("0-50 m DIC map comparison on GLODAP 1-degree grid (common mask):")
print(map_summary.to_string(index=False))
print()

depth_summary = pd.DataFrame(
    {m: model_results[m]["depth_biases"] for m in models},
    index=depth_targets,
).T
depth_summary.columns = [f"{int(d)} m" for d in depth_targets]
print("Fixed-depth global mean DIC bias (model minus GLODAP, umol kg-1):")
print(depth_summary.round(0).to_string())
print()
print("Cell counts at the fixed depths:")
for m in models:
    print(m, dict(zip(depth_targets.astype(int), model_results[m]["depth_counts"])))
print()

period_summary = pd.DataFrame(
    {
        "model": models,
        "DIC_2002": [model_results[m]["dic_0_50"] for m in models],
        "DIC_1971_2000": [model_results[m]["dic_7100_0_50"] for m in models],
        "change": [model_results[m]["dic_change"] for m in models],
        "glodap_diff": [model_results[m]["glodap_diff"] for m in models],
        "ratio": [model_results[m]["ratio"] for m in models],
    }
)
period_summary["DIC_2002"] = period_summary["DIC_2002"].round(0)
period_summary["DIC_1971_2000"] = period_summary["DIC_1971_2000"].round(0)
period_summary["change"] = period_summary["change"].round(0)
period_summary["glodap_diff"] = period_summary["glodap_diff"].round(0)
period_summary["ratio"] = period_summary["ratio"].round(2)
print("Period change of global 0-50 m DIC (2002 minus 1971-2000):")
print(period_summary.to_string(index=False))
print()

print("Reported missing variables:")
for m in models:
    print(f"  {m}: {missing_report[m] if missing_report[m] else 'none'}")
