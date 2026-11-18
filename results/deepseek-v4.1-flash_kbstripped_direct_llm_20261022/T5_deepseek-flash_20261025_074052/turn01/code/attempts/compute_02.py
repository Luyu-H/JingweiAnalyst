import glob
import os
import re
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from scipy import stats
from scipy.spatial import cKDTree

# ---------------------------------------------------------------------------
# Paths and model metadata
# ---------------------------------------------------------------------------
WOA_ROOT = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
WOA_MASK_ROOT = "/data1/lorraine-20419/ocean_data/WOA23/masks"
HIST_ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"
SSP_ROOT = "/data1/lorraine-20419/benchmark_cmip6/ssp245"

MODELS = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

MODEL_SPECS = {
    "CESM2": {
        "y": "nlat", "x": "nlon",
        "lat": "lat", "lon": "lon",
        "lev": "lev", "depth_scale": 0.01,
    },
    "GFDL-ESM4": {
        "y": "lat", "x": "lon",
        "lat": "lat", "lon": "lon",
        "lev": "lev", "depth_scale": 1.0,
    },
    "MPI-ESM1-2-LR": {
        "y": "j", "x": "i",
        "lat": "latitude", "lon": "longitude",
        "lev": "lev", "depth_scale": 1.0,
    },
}

MODEL_COLORS = {
    "CESM2": "tab:blue",
    "GFDL-ESM4": "tab:orange",
    "MPI-ESM1-2-LR": "tab:green",
}

BASIN_NAMES = {
    1: "Atlantic",
    2: "Pacific",
    3: "Indian",
    4: "Arctic",
    5: "Southern",
}
REPORT_BASINS = ["Atlantic", "Indian", "Pacific"]

# ---------------------------------------------------------------------------
# Helpers for WOA mask files
# ---------------------------------------------------------------------------
def _read_woa_mask_rows(path):
    """Read whitespace-separated WOA mask file, skipping headers/comments."""
    rows = []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#") or s.startswith("%"):
                continue
            parts = s.split()
            if len(parts) < 3:
                continue
            try:
                lat = float(parts[0])
                lon = float(parts[1])
            except ValueError:
                continue
            rows.append((lat, lon, parts[2:]))
    return rows


def load_woa_landsea():
    path = os.path.join(WOA_MASK_ROOT, "landsea_01.msk")
    rows = _read_woa_mask_rows(path)
    vals = []
    for lat, lon, rest in rows:
        if not rest:
            continue
        try:
            v = float(rest[0])
        except ValueError:
            continue
        vals.append((lat, lon, v))
    df = pd.DataFrame(vals, columns=["lat", "lon", "bottom"])
    return (
        df.pivot_table(index="lat", columns="lon", values="bottom")
        .sort_index()
        .sort_index(axis=1)
    )


def load_woa_basin_surface():
    path = os.path.join(WOA_MASK_ROOT, "basinmask_01.msk")
    rows = _read_woa_mask_rows(path)
    vals = []
    for lat, lon, rest in rows:
        if not rest:
            continue
        basin = None
        for token in rest:
            try:
                basin = int(float(token))
                break
            except ValueError:
                continue
        if basin is None:
            continue
        vals.append((lat, lon, basin))
    df = pd.DataFrame(vals, columns=["lat", "lon", "basin"])
    return (
        df.pivot_table(index="lat", columns="lon", values="basin")
        .sort_index()
        .sort_index(axis=1)
    )


# ---------------------------------------------------------------------------
# Load WOA23 salinity 1971-2000 annual objective analysis
# ---------------------------------------------------------------------------
woa_file = os.path.join(
    WOA_ROOT,
    "salinity",
    "1971-2000_climate_normal",
    "woa23_decav71A0_s00_01.nc",
)
woa_ds = xr.open_dataset(woa_file, decode_times=False)

# salinity analysed field is normally s_an
if "s_an" in woa_ds:
    woa_s_var = "s_an"
else:
    candidates = [v for v in woa_ds.data_vars if v.endswith("_an")]
    woa_s_var = candidates[0]

woa_s_all = woa_ds[woa_s_var]
if "time" in woa_s_all.dims:
    woa_s_all = woa_s_all.isel(time=0)

depth_name = "depth" if "depth" in woa_s_all.dims else "lev"
depth_vals = woa_ds[depth_name].values.astype(float)
depth_idx = np.where(depth_vals <= 10.0)[0]
woa_s_top = woa_s_all.isel({depth_name: depth_idx})

# Full-level thickness weights for WOA.  Use depth_bnds if present; otherwise
# midpoint bounds.  Cell area cancels in the vertical mean at each point.
thick = None
if "depth_bnds" in woa_ds:
    b = woa_ds["depth_bnds"].values
    if b.ndim == 2 and b.shape[0] == len(depth_vals):
        b_use = b[depth_idx]
        thick = np.abs(b_use[:, 1] - b_use[:, 0])
    elif b.ndim == 2 and b.shape[1] == len(depth_vals):
        b_use = b[:, depth_idx].T
        thick = np.abs(b_use[:, 1] - b_use[:, 0])

if thick is None:
    d = depth_vals
    bounds = np.empty(len(d) + 1)
    bounds[0] = 0.0
    bounds[-1] = d[-1] + 0.5 * (d[-1] - d[-2])
    bounds[1:-1] = 0.5 * (d[:-1] + d[1:])
    thick_all = np.abs(np.diff(bounds))
    thick = thick_all[depth_idx]

thick_da = xr.DataArray(
    thick,
    dims=[depth_name],
    coords={depth_name: woa_s_top[depth_name]},
)
woa_sss = (woa_s_top * thick_da).sum(dim=depth_name) / thick_da.sum()

woa_lat = woa_ds["lat"].values
woa_lon = woa_ds["lon"].values

# Land-sea mask: 1 = land, >1 = ocean in this product.
# The mask coordinates and the WOA coordinates may differ at the level of
# floating-point round-off; use nearest-neighbour selection so the mask keeps
# its boolean dtype (reindex would inject NaN and upcast to object, which then
# breaks the later "& ocean_mask" operations).
landsea_da = xr.DataArray(
    (landsea_df.values > 1),
    dims=["lat", "lon"],
    coords={"lat": landsea_df.index.values, "lon": landsea_df.columns.values},
)
ocean_mask = (
    landsea_da
    .sel(lat=woa_lat, lon=woa_lon, method="nearest")
    .astype(bool)
)

woa_sss = woa_sss.where(ocean_mask)
woa_sss = woa_sss.assign_coords(lat=woa_lat, lon=woa_lon)

# Surface basin mask (same nearest-neighbour approach; land points are masked
# afterwards via ocean_mask)
basin_df = load_woa_basin_surface()
basin_da = xr.DataArray(
    basin_df.values,
    dims=["lat", "lon"],
    coords={"lat": basin_df.index.values, "lon": basin_df.columns.values},
)
basin = basin_da.sel(lat=woa_lat, lon=woa_lon, method="nearest")
basin = basin.where(ocean_mask)

# 1-degree cell-area proxy (cos(lat)); constant factors cancel
woa_area = xr.DataArray(
    np.cos(np.deg2rad(woa_lat))[:, None] * np.ones((1, len(woa_lon))),
    dims=["lat", "lon"],
    coords={"lat": woa_lat, "lon": woa_lon},
)
woa_area = woa_area.where(ocean_mask)

woa_global_mean_1971_2000 = float(
    (woa_sss * woa_area).sum() / woa_area.where(np.isfinite(woa_sss)).sum()
)


# ---------------------------------------------------------------------------
# CMIP6 helpers
# ---------------------------------------------------------------------------
def find_files(root, model, var):
    return sorted(glob.glob(os.path.join(root, model, var, "*.nc")))


def find_files_fallback(root, model, var, fallback_root=None):
    files = find_files(root, model, var)
    if not files and fallback_root is not None:
        files = find_files(fallback_root, model, var)
    return files


def open_cmip(root, model, var, fallback_root=None):
    files = find_files_fallback(root, model, var, fallback_root)
    if not files:
        raise FileNotFoundError(f"No files for {model}/{var}")
    if len(files) == 1:
        ds = xr.open_dataset(files[0], use_cftime=True)
    else:
        ds = xr.open_mfdataset(files, combine="by_coords", use_cftime=True)
    return ds, files


def parse_variant(files):
    for f in files:
        m = re.search(r"_(r\d+i\d+p\d+f\d+)_", os.path.basename(f))
        if m:
            return m.group(1)
    return None


def select_years(da, y0, y1):
    return da.sel(time=(da.time.dt.year >= y0) & (da.time.dt.year <= y1))


def load_model_sss(model, experiment, y0, y1):
    """Return monthly SSS (0-10 m, volume-weighted full levels), area, and variants."""
    spec = MODEL_SPECS[model]
    root = HIST_ROOT if experiment == "historical" else SSP_ROOT
    fallback = HIST_ROOT if experiment == "ssp245" else None

    ds_so, so_files = open_cmip(root, model, "so", None)
    so = ds_so["so"]
    so = select_years(so, y0, y1)

    lev_name = spec["lev"]
    lev_m = ds_so[lev_name].values.astype(float) * spec["depth_scale"]
    lev_idx = np.where((lev_m >= 0.0) & (lev_m <= 10.0))[0]
    so_top = so.isel({lev_name: lev_idx})

    ds_area, area_files = open_cmip(root, model, "areacello", fallback)
    area = ds_area["areacello"]
    if "time" in area.dims:
        area = area.isel(time=0)

    if model == "CESM2":
        ds_sftof, sftof_files = open_cmip(root, model, "sftof", fallback)
        sftof = ds_sftof["sftof"]
        if "time" in sftof.dims:
            sftof = sftof.isel(time=0)
        # apply ocean mask from sftof (do not scale by fraction)
        area = area.where(sftof > 0)

    ds_vol, vol_files = open_cmip(root, model, "volcello", fallback)
    vol = ds_vol["volcello"]
    if "time" in vol.dims:
        vol = select_years(vol, y0, y1)
    vol_top = vol.isel({lev_name: lev_idx})

    so_top, vol_top = xr.align(so_top, vol_top, join="inner")

    valid = np.isfinite(so_top) & np.isfinite(vol_top) & (vol_top > 0)
    num = (so_top * vol_top).where(valid).sum(dim=lev_name)
    den = vol_top.where(valid).sum(dim=lev_name)
    sss = (num / den).where(den > 0)

    return sss, area, so_files, vol_files, parse_variant(so_files), parse_variant(vol_files)


def weighted_global_mean(field, area, dims):
    """Area-weighted mean over selected spatial dims, skipping land/NaN."""
    w = area
    mask = np.isfinite(field) & np.isfinite(w) & (w > 0)
    num = (field.where(mask) * w).sum(dim=dims)
    den = w.where(mask).sum(dim=dims)
    return num / den


def coords_2d(da, spec):
    lat = da.coords[spec["lat"]]
    lon = da.coords[spec["lon"]]
    if lat.ndim == 1 and lon.ndim == 1:
        lat2d = np.repeat(lat.values[:, None], len(lon), axis=1)
        lon2d = np.repeat(lon.values[None, :], len(lat), axis=0)
    else:
        lat2d = lat.values
        lon2d = lon.values
    return lat2d, lon2d


def build_regridder(ref_2d, area_2d, spec, woa_lat, woa_lon):
    """Nearest-neighbour regridder from a model grid to the WOA 1-degree grid."""
    lat2d, lon2d = coords_2d(ref_2d, spec)
    ref = ref_2d.values
    area_vals = area_2d.values
    valid = (
        np.isfinite(ref)
        & np.isfinite(area_vals)
        & (area_vals > 0)
    )
    flat_idx = np.where(valid.ravel())[0]

    lat_f = lat2d.ravel()[flat_idx]
    lon_f = lon2d.ravel()[flat_idx] % 360.0

    target_lon, target_lat = np.meshgrid(woa_lon % 360.0, woa_lat)

    def xyz(lat, lon):
        lat_r = np.deg2rad(lat)
        lon_r = np.deg2rad(lon)
        return np.column_stack(
            (
                np.cos(lat_r) * np.cos(lon_r),
                np.cos(lat_r) * np.sin(lon_r),
                np.sin(lat_r),
            )
        )

    tree = cKDTree(xyz(lat_f, lon_f))
    _, ind = tree.query(xyz(target_lat.ravel(), target_lon.ravel()))
    ind = ind.reshape((len(woa_lat), len(woa_lon)))
    return flat_idx, ind


def regrid_field(field, flat_idx, ind, woa_lat, woa_lon):
    """Apply a nearest-neighbour regridder to a 2-D or leading-dim field."""
    if field.ndim == 2:
        arr = field.values.ravel()[flat_idx]
        out = arr[ind.ravel()].reshape((len(woa_lat), len(woa_lon)))
        return xr.DataArray(
            out, dims=["lat", "lon"], coords={"lat": woa_lat, "lon": woa_lon}
        )
    else:
        lead = field.dims[0]
        nlead = field.sizes[lead]
        arr = field.values.reshape(nlead, -1)[:, flat_idx]
        out = arr[:, ind.ravel()].reshape((nlead, len(woa_lat), len(woa_lon)))
        return xr.DataArray(
            out,
            dims=[lead, "lat", "lon"],
            coords={lead: field[lead], "lat": woa_lat, "lon": woa_lon},
        )


def trend_stats(years, values):
    x = np.asarray(years, dtype=float)
    y = np.asarray(values, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    n = len(x)
    res = stats.linregress(x, y)
    slope_dec = res.slope * 10.0
    se_dec = res.stderr * 10.0
    tcrit = stats.t.ppf(0.975, n - 2)
    ci = (slope_dec - tcrit * se_dec, slope_dec + tcrit * se_dec)
    r2 = res.rvalue ** 2
    excludes_zero = (ci[0] > 0) or (ci[1] < 0)
    return {
        "slope_psu_per_decade": float(slope_dec),
        "ci95_psu_per_decade": (float(ci[0]), float(ci[1])),
        "r2": float(r2),
        "p_value": float(res.pvalue),
        "se_psu_per_decade": float(se_dec),
        "n": int(n),
        "excludes_zero": bool(excludes_zero),
        "intercept": float(res.intercept),
        "slope_per_year": float(res.slope),
    }


def linear_trend_map(da):
    """Vectorised per-cell OLS slope and p-value for annual maps (year, lat, lon)."""
    arr = da.values
    nt = arr.shape[0]
    years = da.year.values.astype(float)
    x = years
    xm = x.mean()
    Sxx = ((x - xm) ** 2).sum()

    valid = np.isfinite(arr).all(axis=0)
    slope = np.full(valid.shape, np.nan)
    pval = np.full(valid.shape, np.nan)

    if not valid.any():
        return slope, pval

    yv = arr[:, valid]
    ym = yv.mean(axis=0)
    slope_v = ((x[:, None] - xm) * (yv - ym[None, :])).sum(axis=0) / Sxx
    intercept_v = ym - slope_v * xm
    fitted = intercept_v[None, :] + slope_v[None, :] * x[:, None]
    resid = yv - fitted
    dof = nt - 2
    sigma2 = (resid ** 2).sum(axis=0) / dof
    se_slope = np.sqrt(sigma2 / Sxx)
    tstat = slope_v / se_slope
    p_v = 2.0 * stats.t.sf(np.abs(tstat), dof)

    slope[valid] = slope_v
    pval[valid] = p_v
    return slope, pval


# ---------------------------------------------------------------------------
# Main model loop
# ---------------------------------------------------------------------------
model_global_mean_1971_2000 = {}
model_bias_global = {}
model_rmse = {}
model_corr = {}

trend_hist = {}
trend_ssp = {}
delta_global = {}
delta_basins = {}

bias_maps = {}
trend_maps = {}
pval_maps = {}

gm_hist_annual_all = {}
gm_ssp_annual_all = {}

variant_hist = {}
variant_ssp = {}
ensemble_caveat = {}

for model in MODELS:
    spec = MODEL_SPECS[model]
    y_dim, x_dim = spec["y"], spec["x"]

    # Historical monthly SSS, 1950-2014
    sss_hist, area, hist_so_files, hist_vol_files, vh_so, vh_vol = load_model_sss(
        model, "historical", 1950, 2014
    )
    variant_hist[model] = vh_so

    # SSP245 monthly SSS, 2071-2100
    sss_ssp, area_ssp, ssp_so_files, ssp_vol_files, vs_so, vs_vol = load_model_sss(
        model, "ssp245", 2071, 2100
    )
    variant_ssp[model] = vs_so

    if vh_so != vs_so:
        ensemble_caveat[model] = (
            f"{model} uses different ensemble members for historical ({vh_so}) "
            f"and ssp245 ({vs_so}); the delta mixes forced change with inter-member variability."
        )
    else:
        ensemble_caveat[model] = (
            f"{model} uses the same ensemble member ({vh_so}) for both experiments."
        )

    # ---- global annual-mean time series ----
    gm_hist_month = weighted_global_mean(sss_hist, area, (y_dim, x_dim)).compute()
    gm_hist_annual = gm_hist_month.groupby("time.year").mean("time")
    gm_hist_annual_all[model] = gm_hist_annual

    gm_ssp_month = weighted_global_mean(sss_ssp, area_ssp, (y_dim, x_dim)).compute()
    gm_ssp_annual = gm_ssp_month.groupby("time.year").mean("time")
    gm_ssp_annual_all[model] = gm_ssp_annual

    # ---- 1971-2000 model climatology and WOA comparison ----
    clim_native = sss_hist.sel(
        time=(sss_hist.time.dt.year >= 1971) & (sss_hist.time.dt.year <= 2000)
    ).mean("time").compute()

    flat_idx, ind = build_regridder(clim_native, area, spec, woa_lat, woa_lon)
    clim_woa = regrid_field(clim_native, flat_idx, ind, woa_lat, woa_lon)

    bias_map = (clim_woa - woa_sss).where(ocean_mask)
    bias_maps[model] = bias_map

    # area-weighted model global mean
    model_global_mean_1971_2000[model] = float(
        weighted_global_mean(
            sss_hist.sel(
                time=(sss_hist.time.dt.year >= 1971)
                & (sss_hist.time.dt.year <= 2000)
            ),
            area,
            (y_dim, x_dim),
        ).mean().compute()
    )
    model_bias_global[model] = (
        model_global_mean_1971_2000[model] - woa_global_mean_1971_2000
    )

    common = np.isfinite(bias_map) & np.isfinite(woa_sss) & ocean_mask
    b = bias_map.values[common.values]
    mv = clim_woa.values[common.values]
    ov = woa_sss.values[common.values]
    model_rmse[model] = float(np.sqrt(np.mean(b ** 2)))
    model_corr[model] = float(np.corrcoef(mv, ov)[0, 1])

    # ---- historical cell-by-cell trend, 1950-2014 ----
    annual_native = sss_hist.groupby("time.year").mean("time").compute()
    annual_woa = regrid_field(annual_native, flat_idx, ind, woa_lat, woa_lon)
    slope_map, p_map = linear_trend_map(annual_woa)
    trend_maps[model] = xr.DataArray(
        slope_map * 10.0,
        dims=["lat", "lon"],
        coords={"lat": woa_lat, "lon": woa_lon},
    ).where(ocean_mask)
    pval_maps[model] = xr.DataArray(
        p_map, dims=["lat", "lon"], coords={"lat": woa_lat, "lon": woa_lon}
    ).where(ocean_mask)

    # ---- global trend statistics ----
    trend_hist[model] = trend_stats(
        gm_hist_annual.year.values, gm_hist_annual.values
    )
    trend_ssp[model] = trend_stats(
        gm_ssp_annual.year.values, gm_ssp_annual.values
    )

    # ---- delta SSS: ssp245 2071-2100 minus historical 1995-2014 ----
    hist_mean_native = sss_hist.sel(
        time=(sss_hist.time.dt.year >= 1995) & (sss_hist.time.dt.year <= 2014)
    ).mean("time").compute()
    ssp_mean_native = sss_ssp.mean("time").compute()

    hist_mean_woa = regrid_field(hist_mean_native, flat_idx, ind, woa_lat, woa_lon)
    ssp_mean_woa = regrid_field(ssp_mean_native, flat_idx, ind, woa_lat, woa_lon)
    delta_map = (ssp_mean_woa - hist_mean_woa).where(ocean_mask)

    delta_global[model] = float(
        (delta_map * woa_area).sum()
        / woa_area.where(np.isfinite(delta_map)).sum()
    )

    delta_basins[model] = {}
    for bid, bname in BASIN_NAMES.items():
        if bname not in REPORT_BASINS:
            continue
        bmask = (basin == bid) & np.isfinite(delta_map) & ocean_mask
        delta_basins[model][bname] = float(
            (delta_map.where(bmask) * woa_area).sum()
            / woa_area.where(bmask).sum()
        )

# ---------------------------------------------------------------------------
# Derived ratios and summary table
# ---------------------------------------------------------------------------
trend_ratio = {}
ratio_meaningful = {}
ratio_note = {}

for model in MODELS:
    denom = trend_hist[model]["slope_psu_per_decade"]
    if abs(denom) > 1e-12:
        trend_ratio[model] = (
            trend_ssp[model]["slope_psu_per_decade"] / denom
        )
    else:
        trend_ratio[model] = np.nan
    ratio_meaningful[model] = False
    ratio_note[model] = (
        "The ssp245/historical slope ratio is not meaningful as a physical ratio: "
        "the two slopes are fitted to disjoint periods (1950-2014 and 2071-2100) "
        "under different forcing scenarios, and for CESM2 they also come from "
        "different ensemble members. It should not be read as a continuous "
        "acceleration or as a scenario-to-historical scaling."
    )

summary_rows = []
for model in MODELS:
    row = {
        "model": model,
        "global_mean_1971_2000_psu": model_global_mean_1971_2000[model],
        "bias_vs_WOA_psu": model_bias_global[model],
        "rmse_vs_WOA_psu": model_rmse[model],
        "corr_vs_WOA": model_corr[model],
        "hist_trend_psu_per_decade": trend_hist[model]["slope_psu_per_decade"],
        "hist_ci_low": trend_hist[model]["ci95_psu_per_decade"][0],
        "hist_ci_high": trend_hist[model]["ci95_psu_per_decade"][1],
        "hist_r2": trend_hist[model]["r2"],
        "hist_excludes_zero": trend_hist[model]["excludes_zero"],
        "ssp_trend_psu_per_decade": trend_ssp[model]["slope_psu_per_decade"],
        "ssp_ci_low": trend_ssp[model]["ci95_psu_per_decade"][0],
        "ssp_ci_high": trend_ssp[model]["ci95_psu_per_decade"][1],
        "ssp_r2": trend_ssp[model]["r2"],
        "ssp_excludes_zero": trend_ssp[model]["excludes_zero"],
        "ssp_over_hist_ratio": trend_ratio[model],
        "delta_global_psu": delta_global[model],
        "delta_Atlantic_psu": delta_basins[model]["Atlantic"],
        "delta_Indian_psu": delta_basins[model]["Indian"],
        "delta_Pacific_psu": delta_basins[model]["Pacific"],
        "historical_variant": variant_hist[model],
        "ssp245_variant": variant_ssp[model],
    }
    summary_rows.append(row)

summary_table = pd.DataFrame(summary_rows)
