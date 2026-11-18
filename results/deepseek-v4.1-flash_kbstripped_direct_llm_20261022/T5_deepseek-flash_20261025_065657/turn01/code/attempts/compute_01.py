import os
import glob
import re
from scipy import stats

# -----------------------------------------------------------------------------
# Paths and model metadata
# -----------------------------------------------------------------------------
WOA_ROOT = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
WOA_MASKS_ROOT = "/data1/lorraine-20419/ocean_data/WOA23/masks"
HIST_ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"
SSP_ROOT = "/data1/lorraine-20419/benchmark_cmip6/ssp245"

MODELS = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
MODEL_INFO = {
    "CESM2": {"lat_coord": "lat", "lon_coord": "lon", "depth_scale": 0.01},
    "GFDL-ESM4": {"lat_coord": "lat", "lon_coord": "lon", "depth_scale": 1.0},
    "MPI-ESM1-2-LR": {"lat_coord": "latitude", "lon_coord": "longitude", "depth_scale": 1.0},
}

# -----------------------------------------------------------------------------
# Small helpers
# -----------------------------------------------------------------------------
def _depth_name(ds):
    for n in ["lev", "depth", "deptht", "z", "st_ocean"]:
        if n in ds.dims or n in ds.coords:
            return n
    return None


def _get_data_var(ds, name):
    if name in ds:
        return ds[name]
    for v in ds.data_vars:
        if v.lower() == name.lower():
            return ds[v]
    return ds[list(ds.data_vars)[0]]


def _depth_values_m(ds, name, default_scale=1.0):
    da = ds[name]
    vals = np.asarray(da.values, dtype=float)
    units = str(da.attrs.get("units", "")).strip().lower()
    if "cm" in units or "centimeter" in units or "centimetre" in units:
        vals = vals * 0.01
    elif units in ("m", "meter", "meters", "metre", "metres") or units.endswith(" m"):
        pass
    else:
        vals = vals * default_scale
    return vals


def _find_model_files(root, model, var, exp=None):
    d = os.path.join(root, model, var)
    files = sorted(glob.glob(os.path.join(d, f"{var}_*.nc")))
    if exp is not None:
        files_exp = [f for f in files if f"_{exp}_" in os.path.basename(f)]
        if files_exp:
            files = files_exp
    return files


def _variant_from_files(files):
    out = set()
    for f in files:
        m = re.search(r"_(r\d+i\d+p\d+f\d+)_", os.path.basename(f))
        if m:
            out.add(m.group(1))
    return sorted(out)


def _open_cmip_dataset(root, model, var, exp, time_slice=None,
                       fallback_root=None, fallback_exp=None):
    files = _find_model_files(root, model, var, exp)
    if not files and fallback_root is not None:
        files = _find_model_files(fallback_root, model, var, fallback_exp or exp)
    if not files:
        d = os.path.join(root, model, var)
        files = sorted(glob.glob(os.path.join(d, f"{var}_*.nc")))
    if not files:
        raise FileNotFoundError(f"No files found for {model} {var} under {root}")

    if len(files) == 1:
        ds = xr.open_dataset(files[0], decode_times=True, use_cftime=True)
    else:
        ds = xr.open_mfdataset(files, combine="by_coords", decode_times=True,
                               use_cftime=True)

    if time_slice is not None and "time" in ds.dims:
        ds = ds.sel(time=slice(time_slice[0], time_slice[1]))
    return ds


def _build_sss_monthly(ds_so, ds_vol, info):
    so = _get_data_var(ds_so, "so")
    dname = _depth_name(ds_so)
    dep = _depth_values_m(ds_so, dname, info.get("depth_scale", 1.0))
    lev_idx = np.where((dep >= 0.0) & (dep <= 10.0))[0]
    if len(lev_idx) == 0:
        lev_idx = np.array([0])

    so_sel = so.isel({dname: lev_idx})

    vol = _get_data_var(ds_vol, "volcello")
    dname_v = _depth_name(ds_vol)
    vol_sel = vol.isel({dname_v: lev_idx})
    if "time" in vol_sel.dims:
        vol_sel = vol_sel.isel(time=0, drop=True)
    if dname_v != dname:
        vol_sel = vol_sel.rename({dname_v: dname})

    num = (so_sel * vol_sel).sum(dim=dname)
    den = vol_sel.sum(dim=dname)
    return num / den


def _area_weighted_global_mean(da, area_np):
    ydim, xdim = da.dims[-2], da.dims[-1]
    area_np = np.asarray(area_np, dtype=float)
    area_np = np.where(np.isfinite(area_np) & (area_np > 0), area_np, np.nan)
    area_da = xr.DataArray(
        area_np,
        dims=(ydim, xdim),
        coords={ydim: da[ydim], xdim: da[xdim]},
    )
    mask = np.isfinite(da) & np.isfinite(area_da)
    num = (da.where(mask) * area_da).sum(dim=(ydim, xdim), skipna=True)
    den = area_da.where(mask).sum(dim=(ydim, xdim), skipna=True)
    return num / den


def _regrid_to_woa(field, lon2d, lat2d, area2d, woa_lat, woa_lon):
    f = np.asarray(field, dtype=float)
    lon = np.asarray(lon2d, dtype=float)
    lat = np.asarray(lat2d, dtype=float)
    area = np.asarray(area2d, dtype=float)
    if f.ndim != 2:
        raise ValueError("field must be 2-D")

    ny, nx = len(woa_lat), len(woa_lon)
    mask = (
        np.isfinite(f) & np.isfinite(lon) & np.isfinite(lat)
        & np.isfinite(area) & (area > 0)
    )
    out_flat = np.full(ny * nx, np.nan, dtype=float)
    if not np.any(mask):
        return out_flat.reshape(ny, nx)

    lonm = ((lon[mask] + 180.0) % 360.0) - 180.0
    latm = lat[mask]
    i = np.floor(lonm + 180.0).astype(int)
    j = np.floor(latm + 90.0).astype(int)
    i = np.clip(i, 0, nx - 1)
    j = np.clip(j, 0, ny - 1)
    flat = j * nx + i

    w = area[mask]
    num = np.bincount(flat, weights=f[mask] * w, minlength=ny * nx)
    den = np.bincount(flat, weights=w, minlength=ny * nx)
    good = den > 0
    out_flat[good] = num[good] / den[good]
    return out_flat.reshape(ny, nx)


def _linear_trend_stats(years, values):
    years = np.asarray(years, dtype=float)
    values = np.asarray(values, dtype=float)
    mask = np.isfinite(years) & np.isfinite(values)
    x = years[mask]
    y = values[mask]
    n = len(x)

    out = dict(
        n=n,
        slope_psu_per_decade=np.nan,
        ci_low_psu_per_decade=np.nan,
        ci_high_psu_per_decade=np.nan,
        r2=np.nan,
        p=np.nan,
        excludes_zero=False,
        intercept=np.nan,
        slope_per_year=np.nan,
    )
    if n < 3:
        return out

    xmean = x.mean()
    ymean = y.mean()
    Sxx = ((x - xmean) ** 2).sum()
    if Sxx <= 0:
        return out

    Sxy = ((x - xmean) * (y - ymean)).sum()
    slope = Sxy / Sxx
    intercept = ymean - slope * xmean
    resid = y - (intercept + slope * x)
    SSE = (resid ** 2).sum()
    SST = ((y - ymean) ** 2).sum()
    r2 = 1.0 - SSE / SST if SST > 0 else np.nan
    se = np.sqrt(SSE / (n - 2) / Sxx)
    tcrit = stats.t.ppf(0.975, n - 2)
    ci_low = slope - tcrit * se
    ci_high = slope + tcrit * se
    p = 2.0 * stats.t.sf(np.abs(slope / se), n - 2) if se > 0 else 1.0

    out.update(dict(
        slope_psu_per_decade=slope * 10.0,
        ci_low_psu_per_decade=ci_low * 10.0,
        ci_high_psu_per_decade=ci_high * 10.0,
        r2=r2,
        p=p,
        excludes_zero=bool(ci_low > 0 or ci_high < 0),
        intercept=intercept,
        slope_per_year=slope,
    ))
    return out


def _trend_map(annual_da):
    years = annual_da["year"].values.astype(float)
    y = np.asarray(annual_da.values, dtype=float)
    mask = np.isfinite(y)
    x = years[:, None, None]

    x0 = np.where(mask, x, 0.0)
    y0 = np.where(mask, y, 0.0)
    count = mask.sum(axis=0).astype(float)
    count_safe = np.where(count > 0, count, np.nan)

    xmean = x0.sum(axis=0) / count_safe
    ymean = y0.sum(axis=0) / count_safe

    Sxx = np.where(mask, (x - xmean[None, :, :]) ** 2, 0.0).sum(axis=0)
    Sxy = np.where(mask, (x - xmean[None, :, :]) * (y - ymean[None, :, :]), 0.0).sum(axis=0)
    slope = Sxy / Sxx
    intercept = ymean - slope * xmean

    resid = np.where(mask, y - (intercept[None, :, :] + slope[None, :, :] * x), 0.0)
    SSE = (resid ** 2).sum(axis=0)
    SST = np.where(mask, (y - ymean[None, :, :]) ** 2, 0.0).sum(axis=0)
    r2 = 1.0 - SSE / SST

    se = np.sqrt(SSE / (count - 2.0) / Sxx)
    tval = slope / se
    p = 2.0 * stats.t.sf(np.abs(tval), count - 2.0)

    bad = (count < 3) | (Sxx <= 0) | ~np.isfinite(slope)
    slope = np.where(bad, np.nan, slope)
    r2 = np.where(bad, np.nan, r2)
    p = np.where(bad, np.nan, p)
    return slope * 10.0, p, r2


def _weighted_mean_field(field, mask, area):
    m = mask & np.isfinite(field)
    if not np.any(m):
        return np.nan
    return np.nansum(field[m] * area[m]) / np.nansum(area[m])


def _read_msk_numeric(path):
    with open(path, "r", errors="ignore") as f:
        lines = f.readlines()

    start = 0
    for i, line in enumerate(lines):
        s = line.strip()
        if not s or s.startswith("#") or s.startswith("!"):
            continue
        parts = re.split(r"[\s,]+", s)
        if len(parts) >= 2:
            try:
                float(parts[0])
                float(parts[1])
                start = i
                break
            except ValueError:
                continue

    df = pd.read_csv(
        path,
        sep=r"[\s,]+",
        comment="#",
        header=None,
        skiprows=start,
        engine="python",
    )
    df = df.dropna(axis=1, how="all")
    for c in [0, 1]:
        df = df[pd.to_numeric(df[c], errors="coerce").notna()]
    return df


def _woa_layer_thickness(depth):
    z = np.asarray(depth, dtype=float)
    dz = np.empty_like(z)
    if len(z) > 1:
        dz[:-1] = z[1:] - z[:-1]
        dz[-1] = dz[-2]
    elif len(z) == 1:
        dz[0] = 1.0
    return np.abs(dz)


# -----------------------------------------------------------------------------
# WOA23 salinity, land-sea mask and basin mask
# -----------------------------------------------------------------------------
woa_file = os.path.join(
    WOA_ROOT,
    "salinity",
    "1971-2000_climate_normal",
    "woa23_decav71A0_s00_01.nc",
)
woa_ds = xr.open_dataset(woa_file, decode_times=False)

woa_var = None
for cand in ["s_an", "salinity", "s", "so"]:
    if cand in woa_ds:
        woa_var = cand
        break
if woa_var is None:
    for cand in woa_ds.data_vars:
        if cand.endswith("_an") or "sal" in cand.lower():
            woa_var = cand
            break
if woa_var is None:
    woa_var = list(woa_ds.data_vars)[0]

woa_da = woa_ds[woa_var]
if "time" in woa_da.dims:
    woa_da = woa_da.isel(time=0, drop=True)

woa_depth_name = _depth_name(woa_ds)
woa_depth = np.asarray(woa_ds[woa_depth_name].values, dtype=float)
woa_lat = np.asarray(woa_ds["lat"].values, dtype=float)
woa_lon = np.asarray(woa_ds["lon"].values, dtype=float)
ny_woa, nx_woa = len(woa_lat), len(woa_lon)

R_EARTH = 6371000.0
woa_area = (
    R_EARTH ** 2
    * np.deg2rad(woa_lat[1] - woa_lat[0])
    * np.deg2rad(woa_lon[1] - woa_lon[0])
    * np.cos(np.deg2rad(woa_lat))[:, None]
    * np.ones((1, nx_woa))
)

# Land-sea mask: ocean is mask > 1.
df_land = _read_msk_numeric(os.path.join(WOA_MASKS_ROOT, "landsea_01.msk"))
lat_land = df_land.iloc[:, 0].values.astype(float)
lon_land = df_land.iloc[:, 1].values.astype(float)
val_land = df_land.iloc[:, 2].values.astype(float)

landsea_arr = np.full((ny_woa, nx_woa), np.nan, dtype=float)
j_land = np.clip(
    np.round((lat_land - woa_lat[0]) / (woa_lat[1] - woa_lat[0])).astype(int),
    0, ny_woa - 1,
)
i_land = np.mod(
    np.round((lon_land - woa_lon[0]) / (woa_lon[1] - woa_lon[0])).astype(int),
    nx_woa,
)
landsea_arr[j_land, i_land] = val_land
woa_ocean = landsea_arr > 1.0

# Basin mask: take the surface basin value from the sparse table.
df_basin = _read_msk_numeric(os.path.join(WOA_MASKS_ROOT, "basinmask_01.msk"))
lat_basin = df_basin.iloc[:, 0].values.astype(float)
lon_basin = df_basin.iloc[:, 1].values.astype(float)
basin_surface_val = df_basin.iloc[:, 2].values.astype(float)

basin_surface = np.full((ny_woa, nx_woa), np.nan, dtype=float)
j_basin = np.clip(
    np.round((lat_basin - woa_lat[0]) / (woa_lat[1] - woa_lat[0])).astype(int),
    0, ny_woa - 1,
)
i_basin = np.mod(
    np.round((lon_basin - woa_lon[0]) / (woa_lon[1] - woa_lon[0])).astype(int),
    nx_woa,
)
basin_surface[j_basin, i_basin] = basin_surface_val

basin_codes = np.unique(basin_surface[np.isfinite(basin_surface)])
if basin_codes.size == 0:
    basin_fold = np.full_like(basin_surface, np.nan)
elif np.nanmax(basin_codes) <= 5:
    # Compact WOA13-style coding: 1 Atlantic, 2 Pacific, 3 Indian, 4 Arctic, 5 Southern.
    basin_fold = np.full_like(basin_surface, np.nan)
    basin_fold[np.isin(basin_surface, [1])] = 1  # Atlantic
    basin_fold[np.isin(basin_surface, [2])] = 2  # Pacific
    basin_fold[np.isin(basin_surface, [3])] = 3  # Indian
else:
    # Extended WOA13-style coding with marginal seas.
    basin_fold = np.full_like(basin_surface, np.nan)
    basin_fold[np.isin(basin_surface, [1, 4, 5, 6])] = 1  # Atlantic + Mediterranean/Black/Baltic
    basin_fold[np.isin(basin_surface, [2])] = 2           # Pacific
    basin_fold[np.isin(basin_surface, [3, 7, 8])] = 3     # Indian + Red Sea/Persian Gulf
    # 9 Arctic and 10 Southern remain NaN (excluded).

# WOA23 SSS over levels whose depths are in 0-10 m, weighted by full layer volume.
woa_lev_idx = np.where((woa_depth >= 0.0) & (woa_depth <= 10.0))[0]
if len(woa_lev_idx) == 0:
    woa_lev_idx = np.array([0])

woa_dz = _woa_layer_thickness(woa_depth)
woa_s_sel = np.asarray(
    woa_da.isel({woa_depth_name: woa_lev_idx}).values, dtype=float
)
woa_vol = woa_dz[woa_lev_idx][:, None, None] * woa_area[None, :, :]

valid_woa = np.isfinite(woa_s_sel)
num_woa = np.nansum(np.where(valid_woa, woa_s_sel * woa_vol, np.nan), axis=0)
den_woa = np.nansum(np.where(valid_woa, woa_vol, np.nan), axis=0)
with np.errstate(invalid="ignore", divide="ignore"):
    woa_sss_np = num_woa / den_woa
woa_sss_np[~woa_ocean] = np.nan

m_woa_global = woa_ocean & np.isfinite(woa_sss_np)
woa_sss_global_mean = (
    np.nansum(woa_sss_np[m_woa_global] * woa_area[m_woa_global])
    / np.nansum(woa_area[m_woa_global])
)

# -----------------------------------------------------------------------------
# CMIP6 processing
# -----------------------------------------------------------------------------
sss_global_mean = {}
sss_global_bias = {}
rmse_against_woa = {}
corr_against_woa = {}
bias_maps = {}

trend_hist_stats = {}
trend_ssp_stats = {}
trend_ratio = {}
trend_ratio_meaning = {}
trend_maps = {}
trend_p_maps = {}

delta_sss = {}
variant_caveats = {}
variant_hist = {}
variant_ssp = {}
global_series = {}

for model in MODELS:
    info = MODEL_INFO[model]

    hist_so_files = _find_model_files(HIST_ROOT, model, "so", "historical")
    ssp_so_files = _find_model_files(SSP_ROOT, model, "so", "ssp245")
    vh = _variant_from_files(hist_so_files)
    vs = _variant_from_files(ssp_so_files)
    variant_hist[model] = vh
    variant_ssp[model] = vs

    if vh and vs and vh != vs:
        variant_caveats[model] = (
            f"historical so variant {','.join(vh)} vs ssp245 so variant "
            f"{','.join(vs)}; the delta mixes different ensemble members."
        )
    elif vh and vs:
        variant_caveats[model] = (
            f"same so variant {','.join(vh)} in both experiments; no "
            "ensemble-member mismatch."
        )
    else:
        variant_caveats[model] = (
            "variant label could not be parsed for one experiment; check filenames."
        )

    ds_so_hist = _open_cmip_dataset(
        HIST_ROOT, model, "so", "historical",
        time_slice=("1950-01-01", "2014-12-31"),
    )
    ds_so_ssp = _open_cmip_dataset(
        SSP_ROOT, model, "so", "ssp245",
        time_slice=("2071-01-01", "2100-12-31"),
    )

    ds_vol_hist = _open_cmip_dataset(
        HIST_ROOT, model, "volcello", "historical",
    )
    ds_vol_ssp = _open_cmip_dataset(
        SSP_ROOT, model, "volcello", "ssp245",
        fallback_root=HIST_ROOT, fallback_exp="historical",
    )

    ds_area = _open_cmip_dataset(
        HIST_ROOT, model, "areacello", "historical",
        fallback_root=SSP_ROOT, fallback_exp="ssp245",
    )
    area_da = _get_data_var(ds_area, "areacello")
    if "time" in area_da.dims:
        area_da = area_da.isel(time=0, drop=True)
    area2d = np.asarray(area_da.values, dtype=float)
    area2d = np.where(np.isfinite(area2d) & (area2d > 0), area2d, np.nan)

    lon_da = ds_so_hist[info["lon_coord"]]
    lat_da = ds_so_hist[info["lat_coord"]]
    if lon_da.ndim == 1 and lat_da.ndim == 1:
        lon2d, lat2d = np.meshgrid(
            np.asarray(lon_da.values, dtype=float),
            np.asarray(lat_da.values, dtype=float),
        )
    else:
        lon2d = np.asarray(lon_da.values, dtype=float)
        lat2d = np.asarray(lat_da.values, dtype=float)

    sss_hist_monthly = _build_sss_monthly(ds_so_hist, ds_vol_hist, info)
    sss_ssp_monthly = _build_sss_monthly(ds_so_ssp, ds_vol_ssp, info)

    years_h = np.array([t.year for t in sss_hist_monthly.indexes["time"]])
    sss_hist_monthly = sss_hist_monthly.assign_coords(year=("time", years_h))
    annual_hist_sss = sss_hist_monthly.groupby("year").mean("time").compute()

    years_s = np.array([t.year for t in sss_ssp_monthly.indexes["time"]])
    sss_ssp_monthly = sss_ssp_monthly.assign_coords(year=("time", years_s))
    annual_ssp_sss = sss_ssp_monthly.groupby("year").mean("time").compute()

    annual_global_hist = _area_weighted_global_mean(
        annual_hist_sss, area2d
    ).compute()
    annual_global_ssp = _area_weighted_global_mean(
        annual_ssp_sss, area2d
    ).compute()

    sss_hist_clim = annual_hist_sss.sel(year=slice(1971, 2000)).mean("year")
    sss_hist_9514 = annual_hist_sss.sel(year=slice(1995, 2014)).mean("year")
    sss_ssp_mean = annual_ssp_sss.mean("year")

    sss_global_mean[model] = float(
        annual_global_hist.sel(year=slice(1971, 2000)).mean().values
    )
    sss_global_bias[model] = sss_global_mean[model] - woa_sss_global_mean

    model_clim_regrid = _regrid_to_woa(
        sss_hist_clim.values, lon2d, lat2d, area2d, woa_lat, woa_lon
    )
    bias = model_clim_regrid - woa_sss_np
    bias_maps[model] = bias

    m_cmp = (
        woa_ocean
        & np.isfinite(woa_sss_np)
        & np.isfinite(model_clim_regrid)
    )
    rmse_against_woa[model] = float(
        np.sqrt(np.mean((model_clim_regrid[m_cmp] - woa_sss_np[m_cmp]) ** 2))
    )
    corr_against_woa[model] = float(
        np.corrcoef(model_clim_regrid[m_cmp], woa_sss_np[m_cmp])[0, 1]
    )

    trend_hist_stats[model] = _linear_trend_stats(
        annual_global_hist["year"].values, annual_global_hist.values
    )
    trend_ssp_stats[model] = _linear_trend_stats(
        annual_global_ssp["year"].values, annual_global_ssp.values
    )

    sh = trend_hist_stats[model]["slope_psu_per_decade"]
    ss = trend_ssp_stats[model]["slope_psu_per_decade"]
    if np.isfinite(sh) and abs(sh) > 1e-12:
        trend_ratio[model] = ss / sh
    else:
        trend_ratio[model] = np.nan

    if (trend_hist_stats[model]["excludes_zero"]
            and trend_ssp_stats[model]["excludes_zero"]):
        trend_ratio_meaning[model] = (
            "ratio pairs slopes from disjoint periods and different forcings; "
            "it is not a continuation slope and should be read only as a "
            "descriptive comparison."
        )
    else:
        trend_ratio_meaning[model] = (
            "ratio is not meaningful because at least one segment's 95% CI "
            "includes zero (or the historical slope is near zero), so it is "
            "numerically unstable."
        )

    tr_map, p_map, _ = _trend_map(annual_hist_sss)
    trend_maps[model] = _regrid_to_woa(
        tr_map, lon2d, lat2d, area2d, woa_lat, woa_lon
    )
    trend_p_maps[model] = _regrid_to_woa(
        p_map, lon2d, lat2d, area2d, woa_lat, woa_lon
    )

    model_hist9514_regrid = _regrid_to_woa(
        sss_hist_9514.values, lon2d, lat2d, area2d, woa_lat, woa_lon
    )
    model_ssp_regrid = _regrid_to_woa(
        sss_ssp_mean.values, lon2d, lat2d, area2d, woa_lat, woa_lon
    )
    delta_regrid = model_ssp_regrid - model_hist9514_regrid

    delta_sss[model] = {
        "global": _weighted_mean_field(delta_regrid, woa_ocean, woa_area),
        "Atlantic": _weighted_mean_field(
            delta_regrid, woa_ocean & (basin_fold == 1), woa_area
        ),
        "Indian": _weighted_mean_field(
            delta_regrid, woa_ocean & (basin_fold == 3), woa_area
        ),
        "Pacific": _weighted_mean_field(
            delta_regrid, woa_ocean & (basin_fold == 2), woa_area
        ),
    }

    hist_years = annual_global_hist["year"].values
    hist_vals = annual_global_hist.values
    ssp_years = annual_global_ssp["year"].values
    ssp_vals = annual_global_ssp.values

    hist_fit = (
        trend_hist_stats[model]["intercept"]
        + trend_hist_stats[model]["slope_per_year"] * hist_years
    )
    ssp_fit = (
        trend_ssp_stats[model]["intercept"]
        + trend_ssp_stats[model]["slope_per_year"] * ssp_years
    )

    global_series[model] = dict(
        hist_years=hist_years,
        hist_sss=hist_vals,
        ssp_years=ssp_years,
        ssp_sss=ssp_vals,
        hist_fit=hist_fit,
        ssp_fit=ssp_fit,
    )
