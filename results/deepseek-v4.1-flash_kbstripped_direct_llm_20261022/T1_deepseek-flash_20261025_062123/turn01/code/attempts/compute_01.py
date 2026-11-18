import glob
import os

HIST_ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"
ERSST_PATH = "/data1/lorraine-20419/ocean_data/ERSSTv5/sst.mnmean.nc"

MODEL_CFG = {
    "CESM2": {
        "y": "nlat", "x": "nlon",
        "lat": "lat", "lon": "lon",
        "mask_from": "sftof",
    },
    "GFDL-ESM4": {
        "y": "lat", "x": "lon",
        "lat": "lat", "lon": "lon",
        "mask_from": None,
    },
    "MPI-ESM1-2-LR": {
        "y": "j", "x": "i",
        "lat": "latitude", "lon": "longitude",
        "mask_from": None,
    },
}

MODELS = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
SOURCES = ["ERSSTv5"] + MODELS


def _find_files(model, field):
    d = os.path.join(HIST_ROOT, model, field)
    patterns = [
        os.path.join(d, f"{field}_*.nc"),
        os.path.join(d, "*.nc"),
    ]
    files = []
    for pat in patterns:
        files = sorted(glob.glob(pat))
        if files:
            break
    if not files:
        files = sorted(glob.glob(os.path.join(HIST_ROOT, model, "**", f"{field}_*.nc"), recursive=True))
    return files


def open_cmip6_field(model, field):
    files = _find_files(model, field)
    if not files:
        raise FileNotFoundError(f"No files for {model}/{field} under {HIST_ROOT}")
    if len(files) == 1:
        ds = xr.open_dataset(files[0], use_cftime=True)
    else:
        ds = xr.open_mfdataset(files, combine="by_coords", use_cftime=True)
    if "time" in ds.dims:
        ds = ds.sortby("time")
    return ds


def select_window(ds):
    if "time" not in ds.dims:
        return ds
    yr = ds.time.dt.year
    m = (yr >= 1950) & (yr <= 2014)
    return ds.isel(time=np.where(m.values)[0]).sortby("time")


def region_mask(ds, lat_name, lon_name):
    lat = ds[lat_name]
    lon = ds[lon_name]
    lon360 = lon % 360.0
    return (lat >= -5.0) & (lat <= 5.0) & (lon360 >= 190.0) & (lon360 <= 240.0)


def weighted_nino34(da, weights, mask, spatial_dims):
    da_m = da.where(mask)
    w_m = weights.where(mask)
    num = (da_m * w_m).sum(dim=spatial_dims, skipna=True)
    den = w_m.where(da_m.notnull()).sum(dim=spatial_dims, skipna=True)
    return num / den


def make_df(idx, source):
    idx = idx.squeeze()
    t = idx.coords["time"]
    years = t.dt.year.values.astype(int)
    months = t.dt.month.values.astype(int)
    vals = np.asarray(idx.values, dtype=float).reshape(-1)

    df = pd.DataFrame({
        "source": source,
        "year": years,
        "month": months,
        "value": vals,
    })
    df = df[(df["year"] >= 1950) & (df["year"] <= 2014)].copy()
    df = df.sort_values(["year", "month"]).drop_duplicates(["year", "month"]).reset_index(drop=True)

    clim = df.groupby("month")["value"].mean()
    df["anom"] = df["value"] - df["month"].map(clim)
    return df


# ERSSTv5
ds_e = xr.open_dataset(ERSST_PATH, decode_times=True)
ds_e = select_window(ds_e)
da_e = ds_e["sst"]
lat_e = ds_e["lat"]
lon_e = ds_e["lon"]

w_e = xr.DataArray(
    np.cos(np.deg2rad(lat_e.values))[:, None] * np.ones((1, lon_e.size)),
    coords={"lat": lat_e, "lon": lon_e},
    dims=("lat", "lon"),
)
mask_e = region_mask(ds_e, "lat", "lon")
idx_e = weighted_nino34(da_e, w_e, mask_e, ("lat", "lon")).compute()

dfs = {"ERSSTv5": make_df(idx_e, "ERSSTv5")}


# CMIP6 historical models
for model in MODELS:
    cfg = MODEL_CFG[model]

    ds_tos = open_cmip6_field(model, "tos")
    ds_tos = select_window(ds_tos)
    da = ds_tos["tos"]

    units = str(da.attrs.get("units", "")).strip().lower()
    if units in ("k", "kelvin"):
        da = da - 273.15

    wds = open_cmip6_field(model, "areacello")
    wname = "areacello" if "areacello" in wds else list(wds.data_vars)[0]
    w = wds[wname]
    if "time" in w.dims:
        w = w.isel(time=0, drop=True)
    w = w.squeeze(drop=True)

    if cfg.get("mask_from"):
        try:
            mds = open_cmip6_field(model, cfg["mask_from"])
            mname = cfg["mask_from"] if cfg["mask_from"] in mds else list(mds.data_vars)[0]
            m = mds[mname]
            if "time" in m.dims:
                m = m.isel(time=0, drop=True)
            m = m.squeeze(drop=True)
            w = w.where(m > 0)
        except Exception:
            pass

    mask = region_mask(ds_tos, cfg["lat"], cfg["lon"])
    idx = weighted_nino34(da, w, mask, (cfg["y"], cfg["x"])).compute()

    dfs[model] = make_df(idx, model)


# Metrics
std_all = {}
seasonal_std = {}
peak_month = {}
peak_std = {}
min_month = {}
min_std = {}
peak_min_ratio = {}
exceed_count = {}
exceed_pct = {}
n_months = {}

for src in SOURCES:
    df = dfs[src]
    anom = df["anom"]
    n = len(anom)
    n_months[src] = n

    std_all[src] = float(anom.std(ddof=0))

    seas = anom.groupby(df["month"]).std(ddof=0).reindex(range(1, 13))
    seasonal_std[src] = seas.values.astype(float)

    peak_month[src] = int(seas.idxmax())
    peak_std[src] = float(seas.max())
    min_month[src] = int(seas.idxmin())
    min_std[src] = float(seas.min())
    peak_min_ratio[src] = float(seas.max() / seas.min())

    exc = anom.abs() > 1.0
    exceed_count[src] = int(exc.sum())
    exceed_pct[src] = float(100.0 * exceed_count[src] / n)

amplitude_bias = {m: float(std_all[m] - std_all["ERSSTv5"]) for m in MODELS}
amplitude_ratio = {m: float(std_all[m] / std_all["ERSSTv5"]) for m in MODELS}
bias_ranking = sorted(MODELS, key=lambda m: abs(amplitude_bias[m]))

obs_peak_month = peak_month["ERSSTv5"]
obs_peak_std = peak_std["ERSSTv5"]
obs_peak_min_ratio = peak_min_ratio["ERSSTv5"]

reproduces_peak_month = {m: bool(peak_month[m] == obs_peak_month) for m in MODELS}
peak_strength_ratio_to_obs = {m: float(peak_std[m] / obs_peak_std) for m in MODELS}
phase_locking_ratio_to_obs = {m: float(peak_min_ratio[m] / obs_peak_min_ratio) for m in MODELS}
