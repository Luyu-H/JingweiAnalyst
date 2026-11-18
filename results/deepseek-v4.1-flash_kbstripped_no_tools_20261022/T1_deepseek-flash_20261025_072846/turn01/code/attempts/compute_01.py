import os
import glob
import re
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

CMIP_ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"
ERSST_FILE = "/data1/lorraine-20419/ocean_data/ERSSTv5/sst.mnmean.nc"

MODEL_INFO = {
    "CESM2": {"lat_coord": "lat", "lon_coord": "lon"},
    "GFDL-ESM4": {"lat_coord": "lat", "lon_coord": "lon"},
    "MPI-ESM1-2-LR": {"lat_coord": "latitude", "lon_coord": "longitude"},
}
MODEL_ORDER = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
SOURCE_ORDER = ["ERSSTv5"] + MODEL_ORDER

MONTH_ABBR = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
MONTH_NAMES = ["January", "February", "March", "April", "May", "June",
               "July", "August", "September", "October", "November", "December"]

target_time = pd.date_range("1950-01-01", "2014-12-01", freq="MS")
target_months = target_time.month.values
N_MONTHS = len(target_time)
if N_MONTHS != 780:
    raise RuntimeError(f"Expected 780 months, got {N_MONTHS}")


def _variant_from_path(path):
    m = re.search(r"_(r\d+i\d+p\d+f\d+)_", os.path.basename(path))
    return m.group(1) if m else None


def _standardize_time(da):
    da = da.sortby("time")
    tvals = da["time"].values
    if np.issubdtype(tvals.dtype, np.datetime64):
        idx = pd.DatetimeIndex(tvals)
        years = idx.year.values
        months = idx.month.values
    else:
        years = np.array([v.year for v in tvals])
        months = np.array([v.month for v in tvals])

    mask = (
        ((years > 1950) | ((years == 1950) & (months >= 1)))
        & ((years < 2014) | ((years == 2014) & (months <= 12)))
    )
    da = da.isel(time=mask)
    if da.sizes["time"] != N_MONTHS:
        raise ValueError(
            f"Expected {N_MONTHS} months for 1950-01..2014-12, got {da.sizes['time']}"
        )
    da = da.assign_coords(time=target_time)
    return da


def _get_surface(da):
    for d in ("lev", "olevel", "depth", "deptht"):
        if d in da.dims:
            return da.isel({d: 0}, drop=True)
    return da


def _open_fixed(model, var):
    files = sorted(glob.glob(os.path.join(CMIP_ROOT, model, var, "*.nc")))
    if not files:
        return None, None
    ds = xr.open_dataset(files[0], decode_times=False)
    if var in ds:
        name = var
    else:
        candidates = [v for v in ds.data_vars if not v.endswith("_bnds")]
        if not candidates:
            raise KeyError(f"No data variable found in {files[0]}")
        name = candidates[0]
    return ds, ds[name]


def _regional_mean(field, area, lat_coord, lon_coord):
    lat = field[lat_coord]
    lon = field[lon_coord]
    lat2d, lon2d = xr.broadcast(lat, lon)
    lon360 = lon2d % 360.0

    # Nino 3.4: 5S-5N, 170W-120W = 190E-240E on a 0-360 longitude axis.
    box = (
        (lat2d >= -5.0)
        & (lat2d <= 5.0)
        & (lon360 >= 190.0)
        & (lon360 <= 240.0)
    )

    area2d = area.broadcast_like(lat2d)
    area_box = area2d.where(box)

    valid = np.isfinite(field)
    weights = area_box.where(valid)

    spatial_dims = [d for d in field.dims if d != "time"]
    num = (field * weights).sum(dim=spatial_dims, skipna=True)
    den = weights.sum(dim=spatial_dims, skipna=True)

    series = (num / den).compute()
    area_series = den.compute()
    return series.values, area_series.values


def load_cmip6_tos(model):
    info = MODEL_INFO[model]
    files = sorted(glob.glob(os.path.join(CMIP_ROOT, model, "tos", "*.nc")))
    if not files:
        raise FileNotFoundError(f"No tos files for {model}")

    variants = {_variant_from_path(f) for f in files if _variant_from_path(f)}
    if len(variants) != 1:
        raise ValueError(f"Expected one historical variant for {model}, found {variants}")
    variant = next(iter(variants))

    ds = xr.open_mfdataset(
        files, combine="by_coords", use_cftime=True, decode_times=True
    )
    field = _standardize_time(ds["tos"])
    field = _get_surface(field)

    area_ds, area = _open_fixed(model, "areacello")
    if area_ds is None:
        raise FileNotFoundError(f"No areacello for {model}")
    area = area.squeeze(drop=True)
    area_source = f"{model} areacello"

    sft_ds, sft = _open_fixed(model, "sftof")
    if sft_ds is not None:
        sft = sft.squeeze(drop=True)
        area = area.where(sft > 0)
        area_source = f"{model} areacello masked by sftof>0"

    series, area_series = _regional_mean(
        field, area, info["lat_coord"], info["lon_coord"]
    )

    ds.close()
    area_ds.close()
    if sft_ds is not None:
        sft_ds.close()

    return series, area_series, variant, area_source


def load_ersst():
    ds = xr.open_dataset(ERSST_FILE, use_cftime=False, decode_times=True)
    field = _standardize_time(ds["sst"])

    lat = ds["lat"]
    lon = ds["lon"]
    R = 6371000.0
    dlat = np.deg2rad(2.0)
    dlon = np.deg2rad(2.0)
    lat_rad = np.deg2rad(lat)
    area_lat = R**2 * dlon * (
        np.sin(lat_rad + dlat / 2.0) - np.sin(lat_rad - dlat / 2.0)
    )
    area = area_lat * xr.ones_like(lon)
    area_source = "ERSSTv5 regular 2-degree cell areas (R=6371 km)"

    series, area_series = _regional_mean(field, area, "lat", "lon")
    ds.close()
    return series, area_series, area_source


def _compute_stats(series, label):
    series = np.asarray(series, dtype=float)
    if series.size != N_MONTHS:
        raise ValueError(f"Expected {N_MONTHS} values, got {series.size}")

    clim = np.array(
        [np.nanmean(series[target_months == m]) for m in range(1, 13)]
    )
    anomaly = series - clim[target_months - 1]

    overall_std = float(np.nanstd(anomaly, ddof=0))
    seasonal_std = np.array(
        [np.nanstd(anomaly[target_months == m], ddof=0) for m in range(1, 13)]
    )

    peak_idx = int(np.nanargmax(seasonal_std))
    min_idx = int(np.nanargmin(seasonal_std))

    exceed_count = int(np.sum(np.abs(anomaly) > 1.0))
    exceed_pct = float(exceed_count / N_MONTHS * 100.0)

    return {
        "label": label,
        "series": series,
        "anomaly": anomaly,
        "climatology": clim,
        "overall_std": overall_std,
        "seasonal_std": seasonal_std,
        "peak_month": peak_idx + 1,
        "peak_std": float(seasonal_std[peak_idx]),
        "min_month": min_idx + 1,
        "min_std": float(seasonal_std[min_idx]),
        "peak_min_ratio": float(seasonal_std[peak_idx] / seasonal_std[min_idx]),
        "exceed_count": exceed_count,
        "exceed_pct": exceed_pct,
    }


def _add_area_stats(stats, area_series, area_source):
    area_m2 = np.asarray(area_series, dtype=float)
    stats["area_source"] = area_source
    stats["area_mean_1e6_km2"] = float(np.nanmean(area_m2) / 1e12)
    stats["area_min_1e6_km2"] = float(np.nanmin(area_m2) / 1e12)
    stats["area_max_1e6_km2"] = float(np.nanmax(area_m2) / 1e12)


results = {}

for model in MODEL_ORDER:
    series, area_series, variant, area_source = load_cmip6_tos(model)
    stats = _compute_stats(series, model)
    stats["variant"] = variant
    _add_area_stats(stats, area_series, area_source)
    results[model] = stats

series, area_series, area_source = load_ersst()
stats = _compute_stats(series, "ERSSTv5")
stats["variant"] = None
_add_area_stats(stats, area_series, area_source)
results["ERSSTv5"] = stats

ersst_std = results["ERSSTv5"]["overall_std"]
for model in MODEL_ORDER:
    results[model]["bias"] = results[model]["overall_std"] - ersst_std
    results[model]["ratio"] = results[model]["overall_std"] / ersst_std

amplitude_ranking = sorted(MODEL_ORDER, key=lambda m: abs(results[m]["bias"]))

ersst_peak_month = results["ERSSTv5"]["peak_month"]
ersst_peak_std = results["ERSSTv5"]["peak_std"]
ersst_peak_min_ratio = results["ERSSTv5"]["peak_min_ratio"]

for model in MODEL_ORDER:
    s = results[model]
    s["peak_month_match"] = bool(s["peak_month"] == ersst_peak_month)
    s["peak_strength_ratio"] = float(s["peak_std"] / ersst_peak_std)
    s["phase_locking_strength_ratio"] = float(
        s["peak_min_ratio"] / ersst_peak_min_ratio
    )

# Round only reported values.
for name, s in results.items():
    s["overall_std_2dp"] = round(float(s["overall_std"]), 2)
    s["seasonal_std_2dp"] = [round(float(v), 2) for v in s["seasonal_std"]]
    s["peak_std_2dp"] = round(float(s["peak_std"]), 2)
    s["min_std_2dp"] = round(float(s["min_std"]), 2)
    s["peak_min_ratio_2dp"] = round(float(s["peak_min_ratio"]), 2)
    s["exceed_pct_1dp"] = round(float(s["exceed_pct"]), 1)

    s["area_mean_1e6_km2_2dp"] = round(float(s["area_mean_1e6_km2"]), 2)
    s["area_min_1e6_km2_2dp"] = round(float(s["area_min_1e6_km2"]), 2)
    s["area_max_1e6_km2_2dp"] = round(float(s["area_max_1e6_km2"]), 2)

    if "bias" in s:
        s["bias_2dp"] = round(float(s["bias"]), 2)
        s["ratio_2dp"] = round(float(s["ratio"]), 2)
        s["peak_strength_ratio_2dp"] = round(float(s["peak_strength_ratio"]), 2)
        s["phase_locking_strength_ratio_2dp"] = round(
            float(s["phase_locking_strength_ratio"]), 2
        )

print("Nino 3.4 anomaly standard deviations, 1950-01..2014-12 (ddof=0)")
for name in SOURCE_ORDER:
    s = results[name]
    print(
        f"{name:14s} std={s['overall_std_2dp']:.2f} degC  "
        f"|anom|>1.0: {s['exceed_count']:3d}/{N_MONTHS} "
        f"({s['exceed_pct_1dp']:.1f}%)  "
        f"area={s['area_mean_1e6_km2_2dp']:.2f} x10^6 km^2"
    )
    print(f"               area source: {s['area_source']}")

print("\nAmplitude bias and ratio against ERSSTv5")
for model in MODEL_ORDER:
    s = results[model]
    print(
        f"{model:14s} bias={s['bias_2dp']:+.2f} degC  "
        f"ratio={s['ratio_2dp']:.2f}  "
        f"rank by |bias|={amplitude_ranking.index(model) + 1}"
    )

print("\nCalendar-month standard deviations (degC)")
print(f"{'source':14s}" + "".join(f"{m:>8s}" for m in MONTH_ABBR))
for name in SOURCE_ORDER:
    vals = results[name]["seasonal_std_2dp"]
    print(f"{name:14s}" + "".join(f"{v:8.2f}" for v in vals))

print("\nPeak / minimum seasonal amplitude")
for name in SOURCE_ORDER:
    s = results[name]
    print(
        f"{name:14s} peak {MONTH_NAMES[s['peak_month'] - 1]:9s} "
        f"{s['peak_std_2dp']:.2f} degC; "
        f"min {MONTH_NAMES[s['min_month'] - 1]:9s} "
        f"{s['min_std_2dp']:.2f} degC; "
        f"peak/min={s['peak_min_ratio_2dp']:.2f}"
    )

print("\nPhase-locking comparison with ERSSTv5")
for model in MODEL_ORDER:
    s = results[model]
    print(
        f"{model:14s} peak-month match={s['peak_month_match']}; "
        f"peak-strength ratio={s['peak_strength_ratio_2dp']:.2f}; "
        f"peak/min ratio relative to ERSST={s['phase_locking_strength_ratio_2dp']:.2f}"
    )

print("\nStrict exceedance counts: months with |Nino 3.4 anomaly| > 1.0 degC")
for name in SOURCE_ORDER:
    s = results[name]
    print(
        f"{name:14s} {s['exceed_count']:3d} / {N_MONTHS} = "
        f"{s['exceed_pct_1dp']:.1f}%"
    )
