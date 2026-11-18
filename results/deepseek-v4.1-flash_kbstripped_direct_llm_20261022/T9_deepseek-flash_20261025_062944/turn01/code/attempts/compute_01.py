import os
import re
import glob
import numpy as np
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree

# ---------------------------------------------------------------------------
# Paths and common PDO domain / target grid
# ---------------------------------------------------------------------------
ersst_path = "/data1/lorraine-20419/ocean_data/ERSSTv5/sst.mnmean.nc"
cpc_root = "/data1/lorraine-20419/ocean_data/CPC"
pdo_path = os.path.join(cpc_root, "ersst.v6.pdo.dat")
cmip_hist_root = "/data1/lorraine-20419/benchmark_cmip6/historical"

models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

# 1-degree centre grid, following the WOA/SEANOE half-degree-centre convention.
target_lat = np.arange(20.5, 70.0, 1.0)
target_lon = np.arange(110.5, 260.0, 1.0)   # 110.5E to 100.5W, 0-360 convention
target_lat2d, target_lon2d = np.meshgrid(target_lat, target_lon, indexing="ij")
ny_t, nx_t = target_lat2d.shape

dates = pd.date_range("1950-01-01", "2014-12-01", freq="MS")
months_expected = dates.month.values
nt_expected = len(dates)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def select_time_range(ds, y0=1950, y1=2014):
    """Select year window and return the dataset plus its month numbers."""
    ds = ds.sortby("time")
    yrs = ds.time.dt.year.values
    mons = ds.time.dt.month.values
    sel = (yrs >= y0) & (yrs <= y1)
    ds_sel = ds.isel(time=np.where(sel)[0])
    return ds_sel, mons[sel]

def to_celsius(da):
    """CMIP6 tos is normally in K; ERSST is already degC."""
    u = str(da.attrs.get("units", "")).strip().lower()
    if u in ("k", "kelvin") or u.endswith("k"):
        return da - 273.15
    return da

def _bilinear_solve(pts, T):
    """
    Invert the bilinear map of a curvilinear grid cell.
    pts = [P00, P10, P01, P11], each [lon, lat].
    Returns (u, v) if T is inside the cell, else None.
    """
    P00, P10, P01, P11 = pts
    a = P00
    b = P10 - P00
    c = P01 - P00
    d = P00 - P10 - P01 + P11

    best = None
    for u0, v0 in ((0.5, 0.5), (0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)):
        u, v = u0, v0
        ok = True
        for _ in range(25):
            F = a + b * u + c * v + d * u * v - T
            if np.linalg.norm(F) < 1e-12:
                break
            J = np.array(
                [
                    [b[0] + d[0] * v, c[0] + d[0] * u],
                    [b[1] + d[1] * v, c[1] + d[1] * u],
                ]
            )
            try:
                delta = np.linalg.solve(J, -F)
            except np.linalg.LinAlgError:
                ok = False
                break
            u += delta[0]
            v += delta[1]
            if not np.isfinite(u) or not np.isfinite(v):
                ok = False
                break
            if np.linalg.norm(delta) < 1e-12:
                break
        if ok and (-1e-5 <= u <= 1.0 + 1e-5) and (-1e-5 <= v <= 1.0 + 1e-5):
            P = a + b * u + c * v + d * u * v
            res = np.linalg.norm(P - T)
            if best is None or res < best[0]:
                best = (res, u, v)

    if best is not None and best[0] < 1e-4:
        return best[1], best[2]
    return None

def make_curvilinear_weights(lat2d, lon2d, target_lat, target_lon, max_candidates=12):
    """
    Precompute bilinear weights from a logically rectangular curvilinear grid
    to the target regular 1-degree grid.
    """
    lat2d = np.asarray(lat2d, dtype=float)
    lon2d = np.asarray(lon2d, dtype=float) % 360.0
    ny, nx = lat2d.shape

    lat_c = 0.25 * (
        lat2d[:-1, :-1] + lat2d[1:, :-1] + lat2d[:-1, 1:] + lat2d[1:, 1:]
    )
    lon_c = 0.25 * (
        lon2d[:-1, :-1] + lon2d[1:, :-1] + lon2d[:-1, 1:] + lon2d[1:, 1:]
    )
    lon_c = lon_c % 360.0

    tree = cKDTree(np.column_stack([lon_c.ravel(), lat_c.ravel()]))

    tlon = (target_lon.ravel() % 360.0)
    tlat = target_lat.ravel()
    target_pts = np.column_stack([tlon, tlat])

    dists, cell_ids = tree.query(target_pts, k=max_candidates)
    if max_candidates == 1:
        cell_ids = cell_ids[:, None]
        dists = dists[:, None]

    n_target = len(tlat)
    idx4 = np.zeros((n_target, 4), dtype=np.int64)
    w4 = np.zeros((n_target, 4), dtype=float)

    for it in range(n_target):
        best = None

        for kk in range(cell_ids.shape[1]):
            cell = cell_ids[it, kk]
            ci, cj = np.unravel_index(cell, (ny - 1, nx - 1))

            pts = np.array(
                [
                    [lon2d[ci, cj], lat2d[ci, cj]],
                    [lon2d[ci + 1, cj], lat2d[ci + 1, cj]],
                    [lon2d[ci, cj + 1], lat2d[ci, cj + 1]],
                    [lon2d[ci + 1, cj + 1], lat2d[ci + 1, cj + 1]],
                ],
                dtype=float,
            )

            clon = lon_c[ci, cj]
            pts[:, 0] = clon + ((pts[:, 0] - clon + 180.0) % 360.0 - 180.0)
            T = np.array(
                [
                    clon + ((tlon[it] - clon + 180.0) % 360.0 - 180.0),
                    tlat[it],
                ]
            )

            if not (
                T[0] >= pts[:, 0].min() - 1e-6
                and T[0] <= pts[:, 0].max() + 1e-6
                and T[1] >= pts[:, 1].min() - 1e-6
                and T[1] <= pts[:, 1].max() + 1e-6
            ):
                continue

            sol = _bilinear_solve(pts, T)
            if sol is None:
                continue

            u, v = sol
            P = (
                (1 - u) * (1 - v) * pts[0]
                + u * (1 - v) * pts[1]
                + (1 - u) * v * pts[2]
                + u * v * pts[3]
            )
            res = np.linalg.norm(P - T)

            if best is None or res < best[0]:
                best = (res, ci, cj, u, v)

        if best is not None:
            _, ci, cj, u, v = best
            w00 = (1 - u) * (1 - v)
            w10 = u * (1 - v)
            w01 = (1 - u) * v
            w11 = u * v
            idx4[it] = [
                ci * nx + cj,
                (ci + 1) * nx + cj,
                ci * nx + (cj + 1),
                (ci + 1) * nx + (cj + 1),
            ]
            w4[it] = [w00, w10, w01, w11]
        else:
            # Fallback: nearest corner of nearest cell.
            cell = cell_ids[it, 0]
            ci, cj = np.unravel_index(cell, (ny - 1, nx - 1))
            flat = np.array(
                [
                    ci * nx + cj,
                    (ci + 1) * nx + cj,
                    ci * nx + (cj + 1),
                    (ci + 1) * nx + (cj + 1),
                ]
            )
            pts = np.array(
                [
                    [lon2d[ci, cj], lat2d[ci, cj]],
                    [lon2d[ci + 1, cj], lat2d[ci + 1, cj]],
                    [lon2d[ci, cj + 1], lat2d[ci, cj + 1]],
                    [lon2d[ci + 1, cj + 1], lat2d[ci + 1, cj + 1]],
                ],
                dtype=float,
            )
            clon = lon_c[ci, cj]
            pts[:, 0] = clon + ((pts[:, 0] - clon + 180.0) % 360.0 - 180.0)
            T = np.array(
                [
                    clon + ((tlon[it] - clon + 180.0) % 360.0 - 180.0),
                    tlat[it],
                ]
            )
            jbest = np.argmin(np.linalg.norm(pts - T, axis=1))
            idx4[it] = [flat[jbest]] * 4
            w4[it] = [1.0, 0.0, 0.0, 0.0]

    return idx4, w4

def apply_curvilinear_weights(arr, idx4, w4):
    """Apply precomputed curvilinear bilinear weights to (time, y, x) data."""
    nt = arr.shape[0]
    flat = arr.reshape(nt, -1)
    out = np.empty((nt, idx4.shape[0]), dtype=float)
    for t in range(nt):
        out[t] = np.sum(flat[t, idx4] * w4, axis=1)
    return out.reshape(nt, target_lat.size, target_lon.size)

def compute_pdo(sst, months, target_lat2d, target_lon2d):
    """
    Leading EOF/PDO of monthly SST anomalies.
    Returns variance fraction, standardized index, regression map, amplitude,
    and lag-12-month autocorrelation.
    """
    sst = np.asarray(sst, dtype=float)
    months = np.asarray(months)

    # Monthly climatology over 1950-2014.
    clim = np.full((12, sst.shape[1], sst.shape[2]), np.nan, dtype=float)
    for m in range(1, 13):
        idxm = months == m
        if np.any(idxm):
            clim[m - 1] = np.nanmean(sst[idxm], axis=0)

    anom = sst - clim[months - 1]

    # Complete time series required at a grid point.
    valid = np.all(np.isfinite(anom), axis=0)
    X = anom[:, valid]

    lat_valid = target_lat2d[valid]
    lon_valid = target_lon2d[valid]

    # Standard PDO convention: area weighting by sqrt(cos(latitude)).
    w = np.sqrt(np.maximum(np.cos(np.deg2rad(lat_valid)), 0.0))
    Xw = X * w[None, :]

    U, S, Vt = np.linalg.svd(Xw, full_matrices=False)
    var_expl = float(S[0] ** 2 / np.sum(S ** 2))

    pc = U[:, 0] * S[0]
    idx_std = pc / np.std(pc, ddof=0)

    # Regression of original SST anomalies onto the standardized PC.
    amp = (X.T @ idx_std) / np.sum(idx_std ** 2)
    amap = np.full(sst.shape[1:], np.nan, dtype=float)
    amap[valid] = amp

    # Fix EOF sign to the PDO-positive convention:
    # warm western North Pacific and cool central/eastern North Pacific.
    west_mask = (
        (target_lat2d >= 20.0)
        & (target_lat2d <= 45.0)
        & (target_lon2d >= 120.0)
        & (target_lon2d <= 180.0)
    )
    east_mask = (
        (target_lat2d >= 20.0)
        & (target_lat2d <= 50.0)
        & (target_lon2d >= 180.0)
        & (target_lon2d <= 260.0)
    )
    west = np.nanmean(amap[west_mask])
    east = np.nanmean(amap[east_mask])
    if np.isfinite(west) and np.isfinite(east) and (west - east) < 0:
        amap = -amap
        idx_std = -idx_std

    amp_max = float(np.nanmax(np.abs(amap)))
    valid_map = np.isfinite(amap)
    w_map = np.sqrt(
        np.maximum(np.cos(np.deg2rad(target_lat2d[valid_map])), 0.0)
    )
    amp_rms = float(
        np.sqrt(np.average(amap[valid_map] ** 2, weights=w_map))
    )

    ac12 = float(np.corrcoef(idx_std[:-12], idx_std[12:])[0, 1])

    return {
        "var_expl": var_expl,
        "idx": idx_std,
        "amp_map": amap,
        "amp_max": amp_max,
        "amp_rms": amp_rms,
        "ac12": ac12,
        "valid": valid,
    }

def load_model_tos(model):
    """Load, select 1950-2014, convert to degC, and bilinearly regrid tos."""
    files = sorted(
        glob.glob(os.path.join(cmip_hist_root, model, "tos", "*.nc"))
    )
    if len(files) == 0:
        raise FileNotFoundError(f"No tos files for {model}")

    # Group by variant, in case more than one variant is present.
    groups = {}
    for f in files:
        parts = os.path.basename(f).split("_")
        variant = "unknown"
        for p in parts:
            if re.match(r"r\d+i\d+p\d+f\d+", p):
                variant = p
                break
        groups.setdefault(variant, []).append(f)

    if len(groups) > 1:
        if "r1i1p1f1" in groups:
            variant = "r1i1p1f1"
        else:
            variant = max(groups, key=lambda k: len(groups[k]))
    else:
        variant = list(groups.keys())[0]

    fs = sorted(groups[variant])
    if len(fs) == 1:
        ds = xr.open_dataset(fs[0], use_cftime=True)
    else:
        ds = xr.open_mfdataset(
            fs,
            combine="by_coords",
            use_cftime=True,
            data_vars="minimal",
            coords="minimal",
            compat="override",
        )

    ds, mons = select_time_range(ds, 1950, 2014)
    da = ds["tos"]
    da = to_celsius(da)

    ydim, xdim = {
        "CESM2": ("nlat", "nlon"),
        "GFDL-ESM4": ("lat", "lon"),
        "MPI-ESM1-2-LR": ("j", "i"),
    }[model]
    da = da.transpose("time", ydim, xdim)

    if model == "GFDL-ESM4":
        da = da.sortby("lat").sortby("lon")
        da_reg = da.interp(
            lat=target_lat, lon=target_lon, method="linear"
        )
        arr = da_reg.values.astype(float)
    else:
        latname = "lat" if model == "CESM2" else "latitude"
        lonname = "lon" if model == "CESM2" else "longitude"
        lat2d = ds[latname].values
        lon2d = ds[lonname].values
        if lat2d.ndim == 1:
            lon2d, lat2d = np.meshgrid(lon2d, lat2d)
        lon2d = lon2d % 360.0

        idx4, w4 = make_curvilinear_weights(
            lat2d, lon2d, target_lat, target_lon
        )
        arr = apply_curvilinear_weights(da.values.astype(float), idx4, w4)

    if len(mons) != nt_expected:
        raise ValueError(
            f"{model}: expected {nt_expected} months, got {len(mons)}"
        )
    return arr, mons

# ---------------------------------------------------------------------------
# Observations: ERSSTv5
# ---------------------------------------------------------------------------
ds_ersst = xr.open_dataset(ersst_path)
ds_ersst, months_obs = select_time_range(ds_ersst, 1950, 2014)
sst_obs = ds_ersst["sst"]
sst_obs = sst_obs.sortby("lat").sortby("lon")
sst_obs_reg = sst_obs.interp(
    lat=target_lat, lon=target_lon, method="linear"
)
arr_obs = sst_obs_reg.values.astype(float)

pdo_obs = compute_pdo(
    arr_obs, months_obs, target_lat2d, target_lon2d
)

# ---------------------------------------------------------------------------
# CMIP6 historical tos
# ---------------------------------------------------------------------------
pdo_models = {}
for model in models:
    arr_model, mons_model = load_model_tos(model)
    pdo_models[model] = compute_pdo(
        arr_model, mons_model, target_lat2d, target_lon2d
    )

# ---------------------------------------------------------------------------
# Published PDO index (NCEI, ERSSTv6)
# ---------------------------------------------------------------------------
pdo_df = pd.read_csv(pdo_path, sep=r"\s+", comment="#")
pdo_df.columns = [str(c) for c in pdo_df.columns]
pdo_df = pdo_df.rename(columns={pdo_df.columns[0]: "Year"})

month_cols = [
    "Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
]
pdo_long = pdo_df.melt(
    id_vars="Year",
    value_vars=[c for c in month_cols if c in pdo_df.columns],
    var_name="Month",
    value_name="pdo",
)
pdo_long["month"] = pd.to_datetime(
    pdo_long["Month"], format="%b"
).dt.month
pdo_long["time"] = pd.to_datetime(
    dict(
        year=pdo_long["Year"].astype(int),
        month=pdo_long["month"],
        day=1,
    )
)
pdo_series = (
    pdo_long.set_index("time")["pdo"].sort_index().replace(99.99, np.nan)
)
pdo_series = pdo_series.reindex(dates)
pub_idx = pdo_series.values.astype(float)

# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
sources = ["ERSSTv5"] + models

var_expl = {"ERSSTv5": pdo_obs["var_expl"]}
amp_max = {"ERSSTv5": pdo_obs["amp_max"]}
amp_rms = {"ERSSTv5": pdo_obs["amp_rms"]}
ac12 = {"ERSSTv5": pdo_obs["ac12"]}
idx_all = {"ERSSTv5": pdo_obs["idx"]}
amp_maps = {"ERSSTv5": pdo_obs["amp_map"]}

for model in models:
    var_expl[model] = pdo_models[model]["var_expl"]
    amp_max[model] = pdo_models[model]["amp_max"]
    amp_rms[model] = pdo_models[model]["amp_rms"]
    ac12[model] = pdo_models[model]["ac12"]
    idx_all[model] = pdo_models[model]["idx"]
    amp_maps[model] = pdo_models[model]["amp_map"]

pattern_corr = {}
for model in models:
    m = np.isfinite(pdo_obs["amp_map"]) & np.isfinite(
        pdo_models[model]["amp_map"]
    )
    pattern_corr[model] = float(
        np.corrcoef(
            pdo_obs["amp_map"][m],
            pdo_models[model]["amp_map"][m],
        )[0, 1]
    )

valid_pub = np.isfinite(pub_idx)
obs_pub_corr = float(
    np.corrcoef(
        pdo_obs["idx"][valid_pub], pub_idx[valid_pub]
    )[0, 1]
)

amp_ratio = {
    model: amp_max[model] / amp_max["ERSSTv5"] for model in models
}
ac12_diff = {model: ac12[model] - ac12["ERSSTv5"] for model in models}
best_model = max(pattern_corr, key=pattern_corr.get)

# ---------------------------------------------------------------------------
# Observed PDO phase-transition years
# ---------------------------------------------------------------------------
def compute_transition_years(idx, dates, window=108):
    """
    9-year (108-month) centred running mean; find sign crossings, then merge
    crossings less than 3 years after the last counted transition.
    """
    s = pd.Series(idx, index=pd.to_datetime(dates))
    smooth = s.rolling(window=window, center=True).mean()
    sv = smooth.dropna()

    crossings = []
    vals = sv.values
    signs = np.sign(vals)
    for i in range(1, len(sv)):
        if signs[i] == 0 or signs[i - 1] == 0:
            continue
        if signs[i] != signs[i - 1]:
            x0 = sv.index[i - 1]
            x1 = sv.index[i]
            y0 = vals[i - 1]
            y1 = vals[i]
            frac = 0.0 if y1 == y0 else (-y0 / (y1 - y0))
            cross_date = x0 + (x1 - x0) * frac
            crossings.append(cross_date.year)

    kept = []
    for yr in crossings:
        if not kept or (yr - kept[-1]) >= 3:
            kept.append(int(yr))
    return kept, smooth

transition_years, smooth_obs = compute_transition_years(
    pdo_obs["idx"], dates, window=108
)

documented_shifts = [1976, 1988, 1998, 2014]
