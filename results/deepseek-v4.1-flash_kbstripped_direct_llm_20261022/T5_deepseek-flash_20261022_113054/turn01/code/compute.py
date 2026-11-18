# ============================================================
# Computation
# ============================================================
import os, glob, io
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from scipy import stats
from scipy.spatial import cKDTree

HIST_ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"
SSP_ROOT  = "/data1/lorraine-20419/benchmark_cmip6/ssp245"
WOA_ROOT  = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
MASK_ROOT = "/data1/lorraine-20419/ocean_data/WOA23/masks"

MODEL_CFG = {
    "CESM2": {
        "y": "nlat", "x": "nlon", "lat": "lat", "lon": "lon",
        "depth": "lev", "depth_scale": 0.01, "ocean_mask_from": "sftof",
    },
    "GFDL-ESM4": {
        "y": "lat", "x": "lon", "lat": "lat", "lon": "lon",
        "depth": "lev", "depth_scale": 1.0, "ocean_mask_from": None,
    },
    "MPI-ESM1-2-LR": {
        "y": "j", "x": "i", "lat": "latitude", "lon": "longitude",
        "depth": "lev", "depth_scale": 1.0, "ocean_mask_from": None,
    },
}
MODEL_NAMES = list(MODEL_CFG.keys())
MODEL_COLORS = {
    "CESM2": "tab:blue",
    "GFDL-ESM4": "tab:orange",
    "MPI-ESM1-2-LR": "tab:green",
}
BASIN_CODES = {1: "Atlantic", 2: "Pacific", 3: "Indian"}


# ------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------
def get_files(root, model, var, exp):
    d = os.path.join(root, model, var)
    files = sorted(glob.glob(os.path.join(d, "*.nc")))
    if exp is not None:
        files = [f for f in files if f"_{exp}_" in os.path.basename(f)]
    return files


def variant_from_file(f):
    parts = os.path.basename(f).split("_")
    return parts[4] if len(parts) >= 5 else None


def choose_variant_files(files):
    if not files:
        return [], None, []
    variants = {}
    for f in files:
        v = variant_from_file(f)
        variants.setdefault(v, []).append(f)
    best_v = max(variants.items(), key=lambda kv: len(kv[1]))[0]
    return variants[best_v], best_v, list(variants.keys())


def open_nc_files(files, use_cftime=True):
    if len(files) == 1:
        return xr.open_dataset(files[0], use_cftime=use_cftime)
    return xr.open_mfdataset(files, combine="by_coords", use_cftime=use_cftime)


def load_fixed(model, exp_root, exp_name, var):
    files = get_files(exp_root, model, var, exp_name)
    source = exp_name
    if not files:
        files = get_files(HIST_ROOT, model, var, "historical")
        source = "historical"
    if not files:
        return None, None, None, None
    files_sel, v, all_v = choose_variant_files(files)
    ds = xr.open_dataset(files_sel[0], use_cftime=True)
    if var not in ds:
        return None, v, source, all_v
    da = ds[var]
    if "time" in da.dims:
        da = da.isel(time=0)
    da = da.squeeze()
    return da, v, source, all_v


def mean_years(da, y0, y1):
    yrs = da.time.dt.year
    m = (yrs >= y0) & (yrs <= y1)
    return da.isel(time=np.where(m)[0]).mean("time")


def global_mean_native(data, area):
    data = np.asarray(data, dtype=float)
    area = np.asarray(area, dtype=float)
    valid = np.isfinite(data) & np.isfinite(area) & (area > 0)
    if not np.any(valid):
        return np.nan
    return float(np.nansum(data[valid] * area[valid]) / np.nansum(area[valid]))


def global_mean_series(sss, area):
    dims = [d for d in sss.dims if d != "time"]
    valid = np.isfinite(sss) & (area > 0)
    num = (sss.where(valid, 0) * area.where(valid, 0)).sum(dims)
    den = area.where(valid, 0).sum(dims)
    gm = num / den
    return gm.groupby("time.year").mean()


def linear_trend_stats(years, values):
    x = np.asarray(years, dtype=float)
    y = np.asarray(values, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]
    n = len(x)
    if n < 3:
        return {
            "slope_per_decade": np.nan, "ci95": (np.nan, np.nan),
            "r2": np.nan, "p": np.nan, "n": n,
            "excludes_zero": False,
        }
    xm = x.mean()
    ym = y.mean()
    sxx = np.sum((x - xm) ** 2)
    if sxx == 0:
        return {
            "slope_per_decade": np.nan, "ci95": (np.nan, np.nan),
            "r2": np.nan, "p": np.nan, "n": n,
            "excludes_zero": False,
        }
    slope = np.sum((x - xm) * (y - ym)) / sxx
    intercept = ym - slope * xm
    resid = y - (intercept + slope * x)
    sse = np.sum(resid ** 2)
    dof = n - 2
    se = np.sqrt(sse / dof / sxx) if sse > 0 else 0.0
    tcrit = stats.t.ppf(0.975, dof)
    slope_dec = slope * 10.0
    ci = (slope_dec - tcrit * se * 10.0, slope_dec + tcrit * se * 10.0)
    if se > 0:
        tstat = slope / se
        p = 2.0 * stats.t.sf(abs(tstat), dof)
    else:
        p = 1.0 if abs(slope) < 1e-15 else 0.0
    sst = np.sum((y - ym) ** 2)
    r2 = 1.0 - sse / sst if sst > 0 else np.nan
    return {
        "slope_per_decade": float(slope_dec),
        "ci95": (float(ci[0]), float(ci[1])),
        "r2": float(r2),
        "p": float(p),
        "n": int(n),
        "excludes_zero": bool(ci[0] > 0 or ci[1] < 0),
    }


def trend_map(annual_da):
    years = annual_da["year"].values.astype(float)
    y = annual_da.values.astype(float)
    n = len(years)
    t = years - years.mean()
    sxx = np.sum(t ** 2)
    y_mean = np.nanmean(y, axis=0)
    num = np.nansum(t[:, None, None] * (y - y_mean[None, :, :]), axis=0)
    slope = num / sxx
    intercept = y_mean - slope * years.mean()
    resid = y - (intercept[None, :, :] + slope[None, :, :] * years[:, None, None])
    sse = np.nansum(resid ** 2, axis=0)
    dof = n - 2
    with np.errstate(invalid="ignore", divide="ignore"):
        se = np.sqrt(sse / dof / sxx)
        tstat = slope / se
        p = 2.0 * stats.t.sf(np.abs(tstat), dof)
    return slope * 10.0, p


def regrid_to_woa(data_native, lat2d, lon2d, woa_lat, woa_lon, woa_ocean, max_dist_deg=10.0):
    data = np.asarray(data_native, dtype=float)
    valid = np.isfinite(data)
    out = np.full((len(woa_lat), len(woa_lon)), np.nan)
    if not np.any(valid):
        return out

    lat = np.asarray(lat2d, dtype=float)[valid]
    lon = np.asarray(lon2d, dtype=float)[valid]
    vals = data[valid]

    lon_rad = np.deg2rad(((lon + 180.0) % 360.0) - 180.0)
    lat_rad = np.deg2rad(lat)
    pts = np.column_stack([
        np.cos(lat_rad) * np.cos(lon_rad),
        np.cos(lat_rad) * np.sin(lon_rad),
        np.sin(lat_rad),
    ])
    tree = cKDTree(pts)

    lon_t, lat_t = np.meshgrid(woa_lon, woa_lat)
    lon_t_rad = np.deg2rad(((lon_t + 180.0) % 360.0) - 180.0)
    lat_t_rad = np.deg2rad(lat_t)
    tpts = np.column_stack([
        (np.cos(lat_t_rad) * np.cos(lon_t_rad)).ravel(),
        (np.cos(lat_t_rad) * np.sin(lon_t_rad)).ravel(),
        np.sin(lat_t_rad).ravel(),
    ])
    dist, idx = tree.query(tpts, k=1)
    vals_out = vals[idx].reshape(len(woa_lat), len(woa_lon))
    dist_deg = np.rad2deg(2.0 * np.arcsin(np.minimum(dist / 2.0, 1.0))).reshape(
        len(woa_lat), len(woa_lon)
    )
    out = np.where(dist_deg <= max_dist_deg, vals_out, np.nan)
    out = np.where(woa_ocean, out, np.nan)
    return out


def get_layer_thickness(ds, depth_name):
    depth = np.asarray(ds[depth_name].values, dtype=float)
    bname = None
    for v in ds.variables:
        low = v.lower()
        if ("bnd" in low or "bound" in low) and depth_name.lower() in low:
            bname = v
            break
    thick = None
    if bname is not None:
        b = np.asarray(ds[bname].values, dtype=float)
        if b.ndim == 2 and b.shape[1] == 2:
            thick = np.abs(b[:, 1] - b[:, 0])
        elif b.ndim == 1 and len(b) == len(depth) + 1:
            thick = np.abs(np.diff(b))
    if thick is None:
        thick = np.empty(len(depth), dtype=float)
        if len(depth) > 1:
            thick[:-1] = np.abs(np.diff(depth))
            thick[-1] = thick[-2]
        else:
            thick[:] = 1.0
    return thick


# ------------------------------------------------------------
# WOA23 mask readers
# ------------------------------------------------------------
def _read_msk_table(path):
    with open(path, "r") as f:
        lines = f.readlines()
    hidx = None
    for i, l in enumerate(lines):
        low = l.lower()
        if "lat" in low and "lon" in low:
            hidx = i
            break
    if hidx is None:
        hidx = 0
        while hidx < len(lines) and (
            lines[hidx].strip().startswith("#") or not lines[hidx].strip()
        ):
            hidx += 1
    header = lines[hidx].lstrip("#").strip()
    data_lines = [header]
    for l in lines[hidx + 1:]:
        s = l.strip()
        if s and not s.startswith("#"):
            data_lines.append(l)
    return pd.read_csv(io.StringIO("".join(data_lines)), delim_whitespace=True)


def load_landsea_mask():
    df = _read_msk_table(os.path.join(MASK_ROOT, "landsea_01.msk"))
    lat_col = [c for c in df.columns if "lat" in c.lower()][0]
    lon_col = [c for c in df.columns if "lon" in c.lower()][0]
    val_col = [
        c for c in df.columns
        if "bottom" in c.lower() or "level" in c.lower() or "sea" in c.lower()
    ][0]

    lat = df[lat_col].values.astype(float)
    lon = df[lon_col].values.astype(float)
    val = df[val_col].values.astype(float)

    mask_std = np.zeros((180, 360), dtype=bool)
    ilat = np.rint(lat - (-89.5)).astype(int)
    ilon = np.rint(lon - (-179.5)).astype(int) % 360
    good = (ilat >= 0) & (ilat < 180) & (ilon >= 0) & (ilon < 360)
    mask_std[ilat[good], ilon[good]] = val[good] > 1
    return mask_std


def load_basin_surface():
    df = _read_msk_table(os.path.join(MASK_ROOT, "basinmask_01.msk"))
    lat_col = [c for c in df.columns if "lat" in c.lower()][0]
    lon_col = [c for c in df.columns if "lon" in c.lower()][0]
    other = [c for c in df.columns if c not in (lat_col, lon_col)]

    lev_cols = [c for c in other if "lev" in c.lower() or "depth" in c.lower()]
    bas_cols = [c for c in other if "basin" in c.lower()]

    if lev_cols and bas_cols:
        lev_col = lev_cols[0]
        bas_col = bas_cols[0]
        lev_vals = pd.to_numeric(df[lev_col], errors="coerce").values.astype(float)
        if np.any(lev_vals == 0):
            surf = 0.0
        elif np.any(lev_vals == 1):
            surf = 1.0
        else:
            surf = np.nanmin(lev_vals)
        df_s = df[lev_vals == surf]
        vals = df_s[bas_col].values.astype(float)
        lat = df_s[lat_col].values.astype(float)
        lon = df_s[lon_col].values.astype(float)
    else:
        bas_col = None
        for c in other:
            if "basin" in c.lower() or "bas" in c.lower():
                bas_col = c
                break
        if bas_col is None:
            for c in other:
                if not str(c).isdigit():
                    bas_col = c
                    break
        vals = df[bas_col].values.astype(float)
        lat = df[lat_col].values.astype(float)
        lon = df[lon_col].values.astype(float)

    basin_std = np.full((180, 360), np.nan, dtype=float)
    ilat = np.rint(lat - (-89.5)).astype(int)
    ilon = np.rint(lon - (-179.5)).astype(int) % 360
    good = (ilat >= 0) & (ilat < 180) & (ilon >= 0) & (ilon < 360)
    basin_std[ilat[good], ilon[good]] = vals[good]
    return basin_std


# ------------------------------------------------------------
# Load WOA23 salinity and masks
# ------------------------------------------------------------
woa_candidates = glob.glob(
    os.path.join(WOA_ROOT, "salinity", "1971-2000_climate_normal",
                 "woa23_decav71A0_s00_01.nc")
)
if not woa_candidates:
    woa_candidates = glob.glob(
        os.path.join(WOA_ROOT, "salinity", "1971-2000_climate_normal", "*.nc")
    )
woa_file = woa_candidates[0]
woa_ds = xr.open_dataset(woa_file, decode_times=False)

# Find the objectively analysed salinity variable, e.g. s_an
vname = None
for v in woa_ds.data_vars:
    low = v.lower()
    if "an" in low and ("s" in low or "sal" in low):
        vname = v
        break
if vname is None:
    vname = [v for v in woa_ds.data_vars
             if any("depth" in d.lower() for d in woa_ds[v].dims)][0]

da = woa_ds[vname]
if "time" in da.dims:
    da = da.isel(time=0)

depth_name = [d for d in woa_ds.dims if "depth" in d.lower()][0]
depth_vals = woa_ds[depth_name].values.astype(float)
lev_idx_woa = np.where((depth_vals >= 0.0) & (depth_vals <= 10.0))[0]
da_top = da.isel({depth_name: lev_idx_woa})

thick_full = get_layer_thickness(woa_ds, depth_name)
thick = thick_full[lev_idx_woa]
thick_da = xr.DataArray(
    thick, dims=[depth_name], coords={depth_name: da_top[depth_name]}
)
valid_top = np.isfinite(da_top)
num_top = (da_top.where(valid_top, 0) * thick_da).sum(depth_name)
den_top = (valid_top * thick_da).sum(depth_name)
woa_sss_da = num_top / den_top

lat_name = [c for c in woa_ds.coords if "lat" in c.lower()][0]
lon_name = [c for c in woa_ds.coords if "lon" in c.lower()][0]
woa_lat = woa_ds[lat_name].values.astype(float)
woa_lon = woa_ds[lon_name].values.astype(float)
woa_sss = woa_sss_da.values.astype(float)

# WOA masks
try:
    landsea_std = load_landsea_mask()
    ilat_woa = np.rint(woa_lat - (-89.5)).astype(int)
    ilon_woa = np.rint(woa_lon - (-179.5)).astype(int) % 360
    woa_ocean = landsea_std[np.ix_(ilat_woa, ilon_woa)]
except Exception:
    woa_ocean = np.isfinite(woa_sss)

try:
    basin_std = load_basin_surface()
    basin_woa = basin_std[np.ix_(ilat_woa, ilon_woa)]
except Exception:
    basin_woa = np.full_like(woa_sss, np.nan, dtype=float)

woa_sss = np.where(woa_ocean, woa_sss, np.nan)
woa_area = np.cos(np.deg2rad(woa_lat))[:, None] * np.ones((1, len(woa_lon)))

valid_woa = woa_ocean & np.isfinite(woa_sss)
woa_global_mean = float(
    np.nansum(woa_sss[valid_woa] * woa_area[valid_woa])
    / np.nansum(woa_area[valid_woa])
)


# ------------------------------------------------------------
# Model SSS loader / segment processor
# ------------------------------------------------------------
def load_model_sss(model, root, exp_name, y0, y1):
    cfg = MODEL_CFG[model]
    files = get_files(root, model, "so", exp_name)
    files_sel, v_so, all_v_so = choose_variant_files(files)
    if not files_sel:
        raise RuntimeError(f"No 'so' files for {model} {exp_name}")

    ds = open_nc_files(files_sel, use_cftime=True)
    yrs = ds.time.dt.year
    tmask = (yrs >= y0) & (yrs <= y1)
    ds = ds.isel(time=np.where(tmask)[0])

    lev_vals = ds[cfg["depth"]].values.astype(float)
    lev_m = lev_vals * cfg["depth_scale"]
    lev_idx = np.where((lev_m >= 0.0) & (lev_m <= 10.0))[0]

    so = ds["so"].isel(lev=lev_idx)

    vol, v_vol, src_vol, _ = load_fixed(model, root, exp_name, "volcello")
    if vol is None:
        raise RuntimeError(f"No volcello for {model} {exp_name}")
    vol = vol.isel(lev=lev_idx)

    area, v_area, src_area, _ = load_fixed(model, root, exp_name, "areacello")
    if area is None:
        raise RuntimeError(f"No areacello for {model} {exp_name}")

    if cfg.get("ocean_mask_from") == "sftof":
        sftof, v_sf, src_sf, _ = load_fixed(model, root, exp_name, "sftof")
        if sftof is not None:
            area = area.where(sftof > 0)
            vol = vol.where(sftof > 0)
            so = so.where(sftof > 0)

    valid = np.isfinite(so) & (vol > 0)
    num = (so.where(valid, 0) * vol.where(valid, 0)).sum("lev")
    den = vol.where(valid, 0).sum("lev")
    sss = (num / den).where(den > 0)

    lat = ds[cfg["lat"]].values.astype(float)
    lon = ds[cfg["lon"]].values.astype(float)
    if lat.ndim == 1 and lon.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon, lat)
    else:
        lat2d = lat
        lon2d = lon

    return {
        "sss": sss,
        "area": area,
        "lat2d": lat2d,
        "lon2d": lon2d,
        "variant_so": v_so,
        "all_variants_so": all_v_so,
        "variant_vol": v_vol,
        "source_vol": src_vol,
        "variant_area": v_area,
        "source_area": src_area,
    }


def process_segment(model, root, exp_name, y0, y1, clim_periods, compute_trend=False):
    d = load_model_sss(model, root, exp_name, y0, y1)
    sss = d["sss"]
    area = d["area"]
    lat2d = d["lat2d"]
    lon2d = d["lon2d"]

    gm_annual = global_mean_series(sss, area)

    clim_maps = {}
    clim_native = {}
    clim_global_native = {}
    for name, (cy0, cy1) in clim_periods.items():
        cn = mean_years(sss, cy0, cy1)
        clim_native[name] = cn
        clim_global_native[name] = global_mean_native(cn, area)
        clim_maps[name] = regrid_to_woa(
            cn.values, lat2d, lon2d, woa_lat, woa_lon, woa_ocean
        )

    trend_woa = None
    p_woa = None
    if compute_trend:
        annual_cell = sss.groupby("time.year").mean("time")
        tr, p = trend_map(annual_cell)
        trend_woa = regrid_to_woa(
            tr, lat2d, lon2d, woa_lat, woa_lon, woa_ocean
        )
        p_woa = regrid_to_woa(
            p, lat2d, lon2d, woa_lat, woa_lon, woa_ocean
        )

    return {
        "gm_annual": gm_annual,
        "clim_maps": clim_maps,
        "clim_native": clim_native,
        "clim_global_native": clim_global_native,
        "trend_woa": trend_woa,
        "p_woa": p_woa,
        "variant_so": d["variant_so"],
        "all_variants_so": d["all_variants_so"],
        "variant_vol": d["variant_vol"],
        "source_vol": d["source_vol"],
        "variant_area": d["variant_area"],
        "source_area": d["source_area"],
    }


# ------------------------------------------------------------
# Process all models
# ------------------------------------------------------------
model_results = {}
bias_maps = {}
trend_maps = {}
p_maps = {}
gm_hist_dict = {}
gm_ssp_dict = {}

for model in MODEL_NAMES:
    hist = process_segment(
        model, HIST_ROOT, "historical", 1950, 2014,
        {"1971-2000": (1971, 2000), "1995-2014": (1995, 2014)},
        compute_trend=True,
    )
    ssp = process_segment(
        model, SSP_ROOT, "ssp245", 2071, 2100,
        {"2071-2100": (2071, 2100)},
        compute_trend=False,
    )

    gm_hist = hist["gm_annual"]
    gm_ssp = ssp["gm_annual"]
    gm_hist_dict[model] = gm_hist
    gm_ssp_dict[model] = gm_ssp

    model_global_1971_2000 = hist["clim_global_native"]["1971-2000"]
    bias_to_woa = model_global_1971_2000 - woa_global_mean

    model_clim_woa = hist["clim_maps"]["1971-2000"]
    bias_map = model_clim_woa - woa_sss

    valid = (
        woa_ocean
        & np.isfinite(model_clim_woa)
        & np.isfinite(woa_sss)
    )
    diff = model_clim_woa[valid] - woa_sss[valid]
    rmse = float(np.sqrt(np.mean(diff ** 2)))
    corr = float(np.corrcoef(model_clim_woa[valid], woa_sss[valid])[0, 1])

    hist_trend = linear_trend_stats(
        gm_hist["year"].values.astype(float), gm_hist.values
    )
    ssp_trend = linear_trend_stats(
        gm_ssp["year"].values.astype(float), gm_ssp.values
    )

    hist_slope = hist_trend["slope_per_decade"]
    ssp_slope = ssp_trend["slope_per_decade"]
    ratio = ssp_slope / hist_slope if abs(hist_slope) > 1e-15 else np.nan
    ratio_meaningful = bool(
        hist_trend["excludes_zero"] and abs(hist_slope) > 1e-4
    )

    hist_1995_2014_woa = hist["clim_maps"]["1995-2014"]
    ssp_2071_2100_woa = ssp["clim_maps"]["2071-2100"]
    delta_woa = ssp_2071_2100_woa - hist_1995_2014_woa

    delta_global_woa = global_mean_native(delta_woa, woa_area)
    delta_global_native = (
        ssp["clim_global_native"]["2071-2100"]
        - hist["clim_global_native"]["1995-2014"]
    )

    delta_basins = {}
    for code, bname in BASIN_CODES.items():
        mask_b = (
            woa_ocean
            & (basin_woa == code)
            & np.isfinite(delta_woa)
        )
        if np.any(mask_b):
            delta_basins[bname] = float(
                np.nansum(delta_woa[mask_b] * woa_area[mask_b])
                / np.nansum(woa_area[mask_b])
            )
        else:
            delta_basins[bname] = np.nan

    variants = {
        "historical_so": hist["variant_so"],
        "ssp245_so": ssp["variant_so"],
        "historical_volcello": hist["variant_vol"],
        "ssp245_volcello": ssp["variant_vol"],
        "historical_areacello": hist["variant_area"],
        "ssp245_areacello": ssp["variant_area"],
    }

    caveat = ""
    if hist["variant_so"] != ssp["variant_so"]:
        caveat = (
            f"historical so is {hist['variant_so']} but ssp245 so is "
            f"{ssp['variant_so']}; the delta change conflates the scenario "
            "change with an ensemble-member difference."
        )

    model_results[model] = {
        "global_mean_1971_2000": float(model_global_1971_2000),
        "bias_to_woa": float(bias_to_woa),
        "rmse_vs_woa": rmse,
        "corr_vs_woa": corr,
        "trend_historical": hist_trend,
        "trend_ssp245": ssp_trend,
        "ratio_ssp_to_hist": float(ratio) if np.isfinite(ratio) else np.nan,
        "ratio_meaningful": ratio_meaningful,
        "delta_global_woa": float(delta_global_woa),
        "delta_global_native": float(delta_global_native),
        "delta_basins": delta_basins,
        "variants": variants,
        "variant_caveat": caveat,
    }

    bias_maps[model] = bias_map
    trend_maps[model] = hist["trend_woa"]
    p_maps[model] = hist["p_woa"]


results = {
    "woa_global_mean_1971_2000": woa_global_mean,
    "models": model_results,
    "basin_codes": BASIN_CODES,
}

# Print a compact summary for convenience
print("WOA23 global mean SSS 1971-2000:", round(woa_global_mean, 2), "psu")
for m in MODEL_NAMES:
    r = model_results[m]
    print("\n", m)
    print("  global mean 1971-2000:", round(r["global_mean_1971_2000"], 2), "psu")
    print("  bias vs WOA23:", round(r["bias_to_woa"], 2), "psu")
    print("  RMSE vs WOA23:", round(r["rmse_vs_woa"], 2), "psu")
    print("  correlation vs WOA23:", round(r["corr_vs_woa"], 3))
    for seg in ["trend_historical", "trend_ssp245"]:
        t = r[seg]
        print(f"  {seg}: slope={t['slope_per_decade']:.3g} psu/dec, "
              f"95% CI=({t['ci95'][0]:.3g}, {t['ci95'][1]:.3g}), "
              f"R2={t['r2']:.3f}, excludes_zero={t['excludes_zero']}")
    print("  ratio ssp/hist:", round(r["ratio_ssp_to_hist"], 1),
          "meaningful:", r["ratio_meaningful"])
    print("  delta global (WOA grid):", round(r["delta_global_woa"], 2), "psu")
    print("  delta global (native):", round(r["delta_global_native"], 2), "psu")
    for b, v in r["delta_basins"].items():
        print(f"  delta {b}:", round(v, 2), "psu")
    if r["variant_caveat"]:
        print("  CAVEAT:", r["variant_caveat"])
