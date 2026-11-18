import os, glob, re
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator
from scipy import stats

cmip_root = "/data1/lorraine-20419/benchmark_cmip6/historical"
ersst_root = "/data1/lorraine-20419/ocean_data/ERSSTv5"
cpc_root = "/data1/lorraine-20419/ocean_data/CPC"

model_info = {
    "CESM2": {"y_dim": "nlat", "x_dim": "nlon", "lat_coord": "lat", "lon_coord": "lon"},
    "GFDL-ESM4": {"y_dim": "lat", "x_dim": "lon", "lat_coord": "lat", "lon_coord": "lon"},
    "MPI-ESM1-2-LR": {"y_dim": "j", "x_dim": "i", "lat_coord": "latitude", "lon_coord": "longitude"},
}
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
sources = ["ERSSTv5"] + models

# 1-degree cell-centre grid over the PDO domain 20N-70N, 110E-260E.
target_lat = np.arange(20.5, 70.0, 1.0)
target_lon = np.arange(110.5, 260.0, 1.0)


def find_files(model, var, table):
    d = os.path.join(cmip_root, model, var)
    if not os.path.isdir(d):
        return []
    files = sorted(glob.glob(os.path.join(d, f"{var}_{table}_{model}_historical_*.nc")))
    if not files:
        files = sorted(glob.glob(os.path.join(d, f"{var}_*_{model}_historical_*.nc")))
    return files


def parse_variant(files):
    for f in files:
        m = re.search(r"_(r\d+i\d+p\d+f\d+)_", os.path.basename(f))
        if m:
            return m.group(1)
    return None


def get_var_files(model, var, table):
    files = find_files(model, var, table)
    if not files:
        raise FileNotFoundError(f"No files for {model} {var} {table}")
    v = parse_variant(files)
    if v is not None:
        files = [f for f in files if f"_{v}_" in os.path.basename(f)]
    return files, v


def open_var(files, var):
    if len(files) == 1:
        ds = xr.open_dataset(files[0], use_cftime=True, decode_times=True)
        da = ds[var].load()
        ds.close()
    else:
        ds = xr.open_mfdataset(files, combine="by_coords", use_cftime=True, decode_times=True)
        da = ds[var].load()
        ds.close()
    return da


def select_years(da, start=1950, end=2014):
    da = da.sortby("time")
    yrs = da.time.dt.year
    m = (yrs >= start) & (yrs <= end)
    return da.isel(time=np.where(m.values)[0])


def regrid_regular_to_pdo(anom_da, lat_name, lon_name, target_lat, target_lon):
    da = anom_da.sortby(lat_name).sortby(lon_name)
    mask_da = xr.DataArray(
        np.isfinite(da.isel(time=0).values).astype(float),
        dims=(lat_name, lon_name),
        coords={lat_name: da[lat_name], lon_name: da[lon_name]},
    )
    anom_regrid = da.interp({lat_name: target_lat, lon_name: target_lon}, method="linear")
    mask_regrid = mask_da.interp({lat_name: target_lat, lon_name: target_lon}, method="linear") > 0.5
    anom_regrid = anom_regrid.where(mask_regrid)
    if lat_name != "lat" or lon_name != "lon":
        anom_regrid = anom_regrid.rename({lat_name: "lat", lon_name: "lon"})
    return anom_regrid


def regrid_curvilinear_to_pdo(anom_da, lat2d, lon2d, target_lat, target_lon, mask_native, buffer=5.0):
    ntime = anom_da.sizes["time"]
    lat_flat = np.asarray(lat2d).ravel()
    lon_flat = np.mod(np.asarray(lon2d).ravel(), 360.0)

    tlon, tlat = np.meshgrid(target_lon, target_lat)
    target_pts = np.column_stack([tlat.ravel(), tlon.ravel()])

    lat_min, lat_max = target_lat.min() - buffer, target_lat.max() + buffer
    lon_min, lon_max = target_lon.min() - buffer, target_lon.max() + buffer
    sel = (
        (lat_flat >= lat_min) & (lat_flat <= lat_max)
        & (lon_flat >= lon_min) & (lon_flat <= lon_max)
    )

    mask_flat = np.asarray(mask_native).ravel()

    # Interpolate a 0/1 ocean mask from all source points in the buffer.
    pts_all = np.column_stack([lat_flat[sel], lon_flat[sel]])
    vals_mask = mask_flat[sel].astype(float)
    interp_mask = LinearNDInterpolator(pts_all, vals_mask, fill_value=0.0)
    mask_target = interp_mask(target_pts).reshape(len(target_lat), len(target_lon)) > 0.5

    # Interpolate anomaly values from ocean points only, then apply the mask.
    ocean_sel = sel & mask_flat
    pts_ocean = np.column_stack([lat_flat[ocean_sel], lon_flat[ocean_sel]])
    vals = anom_da.values.reshape(ntime, -1)
    vals_ocean = vals[:, ocean_sel]
    interp_anom = LinearNDInterpolator(pts_ocean, vals_ocean.T, fill_value=np.nan)
    out = interp_anom(target_pts).T.reshape(ntime, len(target_lat), len(target_lon))
    out = np.where(mask_target[None, :, :], out, np.nan)

    return xr.DataArray(
        out,
        dims=("time", "lat", "lon"),
        coords={"time": anom_da.time, "lat": target_lat, "lon": target_lon},
    )


def process_cmip_model(model):
    info = model_info[model]

    tos_files, _ = get_var_files(model, "tos", "Omon")
    tos = open_var(tos_files, "tos")
    tos = select_years(tos)

    area_files, _ = get_var_files(model, "areacello", "Ofx")
    areacello = open_var(area_files, "areacello").squeeze()

    if model == "CESM2":
        sftof_files, _ = get_var_files(model, "sftof", "Ofx")
        sftof = open_var(sftof_files, "sftof").squeeze()
        areacello = areacello.where(sftof > 0)

    areacello = areacello.where(tos.isel(time=0).notnull())

    # Monthly climatology of each calendar month, then monthly anomalies.
    clim = tos.groupby("time.month").mean(dim="time")
    anom = tos.groupby("time.month") - clim

    y_dim, x_dim = info["y_dim"], info["x_dim"]
    den = areacello.sum(skipna=True)
    num = (anom * areacello).sum(dim=[y_dim, x_dim], skipna=True)
    gm = num / den
    anom_g = anom - gm

    lat_coord, lon_coord = info["lat_coord"], info["lon_coord"]
    if tos[lat_coord].ndim == 1:
        anom_pdo = regrid_regular_to_pdo(anom_g, lat_coord, lon_coord, target_lat, target_lon)
    else:
        lat2d = tos[lat_coord].values
        lon2d = tos[lon_coord].values
        mask_native = np.isfinite(tos.isel(time=0).values)
        anom_pdo = regrid_curvilinear_to_pdo(
            anom_g, lat2d, lon2d, target_lat, target_lon, mask_native=mask_native
        )
    return anom_pdo


def process_ersst():
    ds = xr.open_dataset(
        os.path.join(ersst_root, "sst.mnmean.nc"),
        use_cftime=True, decode_times=True
    )
    sst = ds["sst"].load()
    ds.close()
    sst = select_years(sst)

    clim = sst.groupby("time.month").mean(dim="time")
    anom = sst.groupby("time.month") - clim

    lat = sst["lat"]
    lon = sst["lon"]
    coslat = np.cos(np.deg2rad(lat.values))
    weights = xr.DataArray(
        np.broadcast_to(coslat[:, None], (len(lat), len(lon))),
        dims=("lat", "lon"),
        coords={"lat": lat, "lon": lon},
    )
    weights = weights.where(sst.isel(time=0).notnull())
    den = weights.sum(skipna=True)
    num = (anom * weights).sum(dim=["lat", "lon"], skipna=True)
    gm = num / den
    anom_g = anom - gm

    return regrid_regular_to_pdo(anom_g, "lat", "lon", target_lat, target_lon)


def compute_pdo(anom_da, source_name):
    lat = anom_da["lat"].values
    lon = anom_da["lon"].values
    X = anom_da.values
    ntime, nlat, nlon = X.shape

    X2 = X.reshape(ntime, -1)
    finite = np.isfinite(X2).all(axis=0)
    X2 = X2[:, finite]

    lat_grid = np.repeat(lat[:, None], nlon, axis=1).ravel()[finite]
    lon_grid = np.tile(lon, nlat).ravel()[finite]
    coslat = np.cos(np.deg2rad(lat_grid))
    w = np.sqrt(coslat)

    Xw = X2 * w[None, :]
    U, S, Vt = np.linalg.svd(Xw, full_matrices=False)
    pc = U[:, 0] * S[0]
    var_explained = float(S[0] ** 2 / np.sum(S ** 2) * 100.0)

    pc_std = (pc - pc.mean()) / pc.std(ddof=0)
    pattern_masked = (X2 * pc_std[:, None]).mean(axis=0)

    # Positive PDO index: negative SST anomalies in the central/western North Pacific.
    cw = (lat_grid >= 30) & (lat_grid <= 45) & (lon_grid >= 140) & (lon_grid <= 180)
    if cw.any() and np.nanmean(pattern_masked[cw]) > 0:
        pc_std = -pc_std
        pattern_masked = -pattern_masked

    lag12 = float(np.corrcoef(pc_std[:-12], pc_std[12:])[0, 1])
    amplitude = float(np.sqrt(np.sum(pattern_masked ** 2 * coslat) / np.sum(coslat)))

    pattern_full = np.full(nlat * nlon, np.nan)
    pattern_full[finite] = pattern_masked
    pattern_full = pattern_full.reshape(nlat, nlon)

    return {
        "source": source_name,
        "pc": pc_std,
        "time": anom_da["time"],
        "lat": lat,
        "lon": lon,
        "pattern": pattern_full,
        "var_explained": var_explained,
        "amplitude": amplitude,
        "lag12": lag12,
        "mask": finite.reshape(nlat, nlon),
    }


def weighted_pattern_corr(pa, pb, lat, lon):
    A = np.asarray(pa).ravel()
    B = np.asarray(pb).ravel()
    mask = np.isfinite(A) & np.isfinite(B)

    lat_grid = np.repeat(lat[:, None], len(lon), axis=1).ravel()[mask]
    w = np.cos(np.deg2rad(lat_grid))
    a = A[mask]
    b = B[mask]

    wsum = w.sum()
    am = np.sum(w * a) / wsum
    bm = np.sum(w * b) / wsum
    cov = np.sum(w * (a - am) * (b - bm))
    va = np.sum(w * (a - am) ** 2)
    vb = np.sum(w * (b - bm) ** 2)
    corr = float(cov / np.sqrt(va * vb))
    return corr, int(mask.sum())


pdo_results = {}
pdo_results["ERSSTv5"] = compute_pdo(process_ersst(), "ERSSTv5")
for model in models:
    pdo_results[model] = compute_pdo(process_cmip_model(model), model)

pattern_corrs = {}
for model in models:
    corr, n_shared = weighted_pattern_corr(
        pdo_results[model]["pattern"],
        pdo_results["ERSSTv5"]["pattern"],
        pdo_results["ERSSTv5"]["lat"],
        pdo_results["ERSSTv5"]["lon"],
    )
    pattern_corrs[model] = {"corr": corr, "n_shared": n_shared}

source_times = {}
for s in sources:
    t = pdo_results[s]["time"]
    source_times[s] = pd.to_datetime(pd.DataFrame({
        "year": t.dt.year.values,
        "month": t.dt.month.values,
        "day": 1,
    }))

ersst_index_series = pd.Series(
    pdo_results["ERSSTv5"]["pc"],
    index=source_times["ERSSTv5"]
)

# Published PDO index: year-row ASCII file.
pdo_path = os.path.join(cpc_root, "ersst.v6.pdo.dat")
pdo_rows = []
with open(pdo_path, "r") as f:
    for line in f:
        sline = line.strip()
        if not sline or sline.startswith("#"):
            continue
        parts = sline.replace(",", " ").split()
        if len(parts) >= 13 and re.fullmatch(r"\d{4}", parts[0]):
            year = int(parts[0])
            vals = []
            for x in parts[1:13]:
                try:
                    vals.append(float(x))
                except ValueError:
                    vals.append(np.nan)
            if len(vals) == 12:
                pdo_rows.append([year] + vals)

pdo_df = pd.DataFrame(pdo_rows, columns=["Year"] + list(range(1, 13)))
pdo_long = pdo_df.melt(id_vars="Year", var_name="month", value_name="pdo")
pdo_long["time"] = pd.to_datetime(pd.DataFrame({
    "year": pdo_long["Year"],
    "month": pdo_long["month"],
    "day": 1,
}))
published_pdo_series = pdo_long.set_index("time")["pdo"].sort_index()
published_pdo_series = published_pdo_series.replace(99.99, np.nan)
published_pdo_series = published_pdo_series.loc["1950-01-01":"2014-12-01"]

published_pdo_on_ersst_time = published_pdo_series.reindex(source_times["ERSSTv5"])
common = pd.concat([
    ersst_index_series.rename("ersst"),
    published_pdo_on_ersst_time.rename("pdo"),
], axis=1).dropna()

pdo_corr_ersst_published = float(common["ersst"].corr(common["pdo"]))

# Effective-sample-size-corrected significance for the two monthly series.
r1 = common["ersst"].autocorr(lag=1)
r2 = common["pdo"].autocorr(lag=1)
N = len(common)
N_eff = float(N * (1.0 - r1 * r2) / (1.0 + r1 * r2))
if N_eff > 2 and abs(pdo_corr_ersst_published) < 1.0:
    t_stat = pdo_corr_ersst_published * np.sqrt(
        (N_eff - 2.0) / (1.0 - pdo_corr_ersst_published ** 2)
    )
    pdo_corr_p_eff = float(2.0 * stats.t.sf(abs(t_stat), df=N_eff - 2.0))
else:
    pdo_corr_p_eff = np.nan

# PDO transition years from the ERSSTv5 standardized index.
smoothed_ersst = ersst_index_series.rolling(
    window=108, center=True, min_periods=108
).mean()

s_smooth = smoothed_ersst.dropna()
crossings = []
for i in range(1, len(s_smooth)):
    a = s_smooth.iloc[i - 1]
    b = s_smooth.iloc[i]
    if np.isfinite(a) and np.isfinite(b) and a * b < 0:
        crossings.append(s_smooth.index[i - 1].year)

transition_years = []
for yr in sorted(set(crossings)):
    if not transition_years or yr - transition_years[-1] >= 3:
        transition_years.append(int(yr))

documented_pacific_shifts = [1976, 1989, 1998]
documented_matches = {}
for yr in transition_years:
    matches = [d for d in documented_pacific_shifts if abs(yr - d) <= 2]
    if matches:
        documented_matches[yr] = matches

summary = pd.DataFrame([
    {
        "source": s,
        "explained_variance_pct": round(pdo_results[s]["var_explained"], 1),
        "amplitude_C_per_std": round(pdo_results[s]["amplitude"], 3),
        "lag12_autocorr": round(pdo_results[s]["lag12"], 3),
    }
    for s in sources
])
summary["pattern_corr_vs_ERSSTv5"] = np.nan
summary["shared_cells"] = np.nan
for model in models:
    idx = summary.index[summary["source"] == model][0]
    summary.loc[idx, "pattern_corr_vs_ERSSTv5"] = round(pattern_corrs[model]["corr"], 3)
    summary.loc[idx, "shared_cells"] = int(pattern_corrs[model]["n_shared"])

best_model = max(pattern_corrs, key=lambda m: pattern_corrs[m]["corr"])

ersst_amp = pdo_results["ERSSTv5"]["amplitude"]
ersst_lag12 = pdo_results["ERSSTv5"]["lag12"]
model_eval = pd.DataFrame([
    {
        "model": m,
        "pattern_corr": pattern_corrs[m]["corr"],
        "shared_cells": pattern_corrs[m]["n_shared"],
        "amplitude_ratio_to_ERSSTv5": pdo_results[m]["amplitude"] / ersst_amp,
        "lag12_diff_to_ERSSTv5": pdo_results[m]["lag12"] - ersst_lag12,
        "explained_variance_diff_to_ERSSTv5": (
            pdo_results[m]["var_explained"] - pdo_results["ERSSTv5"]["var_explained"]
        ),
    }
    for m in models
])

pdo_corr_interpretation = (
    "The Pearson correlation between the ERSSTv5-derived standardized leading-mode index "
    "and the published PDO monthly series measures temporal agreement in phase and variability "
    "over 1950-2014. It does not prove that the computed EOF is uniquely the PDO, does not "
    "validate the CMIP6 models, and does not establish a causal mechanism; no detrending was "
    "applied, so shared trend remains part of the agreement."
)

key_results = {
    "pattern_correlations": pattern_corrs,
    "pdo_correlation_with_published": round(pdo_corr_ersst_published, 3),
    "pdo_correlation_p_eff": pdo_corr_p_eff,
    "pdo_correlation_N_eff": N_eff,
    "transition_years": transition_years,
    "documented_matches": documented_matches,
    "best_model": best_model,
    "summary": summary,
    "model_eval": model_eval,
}
