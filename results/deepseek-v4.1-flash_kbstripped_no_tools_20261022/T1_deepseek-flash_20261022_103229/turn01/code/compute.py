import os
import glob
import re

root_cmip = "/data1/lorraine-20419/benchmark_cmip6/historical"
ersst_path = "/data1/lorraine-20419/ocean_data/ERSSTv5/sst.mnmean.nc"

model_info = {
    "CESM2": {
        "grid_label": "gn",
        "y_dim": "nlat",
        "x_dim": "nlon",
        "lat_coord": "lat",
        "lon_coord": "lon",
        "ocean_mask_from": "sftof",
    },
    "GFDL-ESM4": {
        "grid_label": "gr",
        "y_dim": "lat",
        "x_dim": "lon",
        "lat_coord": "lat",
        "lon_coord": "lon",
        "ocean_mask_from": None,
    },
    "MPI-ESM1-2-LR": {
        "grid_label": "gn",
        "y_dim": "j",
        "x_dim": "i",
        "lat_coord": "latitude",
        "lon_coord": "longitude",
        "ocean_mask_from": None,
    },
}

models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
sources = ["ERSSTv5"] + models


def find_cmip_files(model, variable, table, grid):
    d = os.path.join(root_cmip, model, variable)
    pattern = os.path.join(d, f"{variable}_{table}_{model}_historical_*_{grid}*.nc")
    return sorted(glob.glob(pattern))


def parse_variant(files):
    for f in files:
        m = re.search(r"_historical_(r\d+i\d+p\d+f\d+)_", os.path.basename(f))
        if m:
            return m.group(1)
    return None


def open_cmip_dataset(files, **kwargs):
    if not files:
        raise FileNotFoundError("No CMIP6 files found")
    dss = [xr.open_dataset(f, **kwargs) for f in files]
    if len(dss) == 1:
        return dss[0]
    ds = xr.concat(dss, dim="time")
    if "time" in ds.coords:
        ds = ds.sortby("time")
    return ds


def select_years_months(da, y0, y1):
    time = da["time"]
    ym = time.dt.year.values.astype(int) * 100 + time.dt.month.values.astype(int)
    idx = np.where((ym >= y0 * 100 + 1) & (ym <= y1 * 100 + 12))[0]
    return da.isel(time=idx)


def get_coord(ds, name):
    if name in ds.coords:
        return ds.coords[name]
    if name in ds.data_vars:
        return ds[name]
    raise KeyError(f"Coordinate {name!r} not found")


def squeeze_singletons(da, keep_dims):
    for d in list(da.dims):
        if d not in keep_dims and da.sizes[d] == 1:
            da = da.isel({d: 0}, drop=True)
    return da


def compute_nino34_cmip(model):
    info = model_info[model]
    grid = info["grid_label"]

    files_all = find_cmip_files(model, "tos", "Omon", grid)
    if not files_all:
        raise FileNotFoundError(f"No tos files for {model}")
    variant = parse_variant(files_all)
    files = [f for f in files_all if variant in os.path.basename(f)] if variant else files_all
    if not files:
        files = files_all

    ds = open_cmip_dataset(files, use_cftime=True, decode_times=True)
    tos = ds["tos"]
    tos = select_years_months(tos, 1950, 2014).sortby("time")
    tos = squeeze_singletons(tos, keep_dims=("time", info["y_dim"], info["x_dim"]))
    if tos.sizes["time"] != 780:
        raise ValueError(f"{model} tos has {tos.sizes['time']} months, expected 780")
    extra = [d for d in tos.dims if d not in ("time", info["y_dim"], info["x_dim"])]
    if extra:
        raise ValueError(f"Unexpected dimensions in {model} tos: {extra}")

    area_files_all = find_cmip_files(model, "areacello", "Ofx", grid)
    if not area_files_all:
        raise FileNotFoundError(f"No areacello files for {model}")
    area_files = [f for f in area_files_all if variant in os.path.basename(f)] if variant else area_files_all
    if not area_files:
        area_files = area_files_all[:1]
    area_ds = open_cmip_dataset(area_files, use_cftime=True, decode_times=True)
    area = area_ds["areacello"]
    if "time" in area.dims:
        area = area.isel(time=0, drop=True)
    area = squeeze_singletons(area, keep_dims=(info["y_dim"], info["x_dim"]))
    area = area.where(area > 0, 0.0).fillna(0.0)

    if info.get("ocean_mask_from") == "sftof":
        sftof_files_all = find_cmip_files(model, "sftof", "Ofx", grid)
        if not sftof_files_all:
            raise FileNotFoundError(f"sftof required for {model} but not found")
        sftof_files = [f for f in sftof_files_all if variant in os.path.basename(f)] if variant else sftof_files_all
        if not sftof_files:
            sftof_files = sftof_files_all[:1]
        sftof_ds = open_cmip_dataset(sftof_files, use_cftime=True, decode_times=True)
        sftof = sftof_ds["sftof"]
        if "time" in sftof.dims:
            sftof = sftof.isel(time=0, drop=True)
        sftof = squeeze_singletons(sftof, keep_dims=(info["y_dim"], info["x_dim"]))
        area = area.where((sftof > 0) & sftof.notnull(), 0.0)

    extra_area = [d for d in area.dims if d not in (info["y_dim"], info["x_dim"])]
    if extra_area:
        raise ValueError(f"Unexpected dimensions in {model} areacello: {extra_area}")

    lat = get_coord(ds, info["lat_coord"])
    lon = get_coord(ds, info["lon_coord"])
    lon360 = lon % 360.0
    region = (lat >= -5.0) & (lat <= 5.0) & (lon360 >= 190.0) & (lon360 <= 240.0)

    weights = area * region
    valid = tos.notnull()
    w = weights.where(valid, 0.0)
    num = (tos.fillna(0.0) * w).sum(dim=(info["y_dim"], info["x_dim"]))
    den = w.sum(dim=(info["y_dim"], info["x_dim"]))
    series = (num / den).where(den > 0)
    series = series.compute()
    series.name = "nino34_sst"
    return series, variant


def compute_nino34_ersst():
    ds = xr.open_dataset(ersst_path, use_cftime=True, decode_times=True)
    sst = ds["sst"]
    sst = select_years_months(sst, 1950, 2014).sortby("time")
    if sst.sizes["time"] != 780:
        raise ValueError(f"ERSSTv5 sst has {sst.sizes['time']} months, expected 780")

    lat = sst["lat"]
    lon = sst["lon"] % 360.0
    region = (lat >= -5.0) & (lat <= 5.0) & (lon >= 190.0) & (lon <= 240.0)

    w2d = xr.DataArray(
        np.cos(np.deg2rad(lat.values))[:, None] * np.ones((lat.size, lon.size)),
        dims=("lat", "lon"),
        coords={"lat": lat, "lon": sst["lon"]},
    )
    weights = w2d * region
    valid = sst.notnull()
    w = weights.where(valid, 0.0)
    num = (sst.fillna(0.0) * w).sum(dim=("lat", "lon"))
    den = w.sum(dim=("lat", "lon"))
    series = (num / den).where(den > 0)
    series = series.compute()
    series.name = "nino34_sst"
    return series, None


nino34_series = {}
variants = {}
nino34_series["ERSSTv5"], variants["ERSSTv5"] = compute_nino34_ersst()
for model in models:
    nino34_series[model], variants[model] = compute_nino34_cmip(model)


def compute_stats(series, source):
    vals = np.asarray(series.values, dtype=float)
    months = np.asarray(series.time.dt.month.values, dtype=int)
    if len(vals) != 780:
        raise ValueError(f"{source} has {len(vals)} values, expected 780")

    clim = np.array([np.nanmean(vals[months == m]) for m in range(1, 13)], dtype=float)
    anomaly = vals - clim[months - 1]

    std_all = float(np.nanstd(anomaly, ddof=0))
    std_month = np.array(
        [np.nanstd(anomaly[months == m], ddof=0) for m in range(1, 13)],
        dtype=float,
    )

    peak_idx = int(np.nanargmax(std_month))
    min_idx = int(np.nanargmin(std_month))
    peak_month = peak_idx + 1
    min_month = min_idx + 1
    peak_value = float(std_month[peak_idx])
    min_value = float(std_month[min_idx])
    ratio = float(peak_value / min_value) if min_value > 0 else np.nan

    n_exceed = int(np.sum(np.abs(anomaly) > 1.0))
    pct = float(100.0 * n_exceed / len(anomaly))

    decimal_year = (
        series.time.dt.year.values.astype(float)
        + (series.time.dt.month.values.astype(float) - 1.0) / 12.0
    )

    return {
        "source": source,
        "values": vals,
        "months": months,
        "climatology": clim,
        "anomaly": anomaly,
        "std_all": std_all,
        "std_month": std_month,
        "peak_month": peak_month,
        "peak_value": peak_value,
        "min_month": min_month,
        "min_value": min_value,
        "ratio": ratio,
        "n_exceed": n_exceed,
        "pct": pct,
        "decimal_year": decimal_year,
    }


stats = {src: compute_stats(nino34_series[src], src) for src in sources}

ersst_std = stats["ERSSTv5"]["std_all"]
amplitude = {}
for model in models:
    model_std = stats[model]["std_all"]
    amplitude[model] = {
        "model_std": model_std,
        "bias": model_std - ersst_std,
        "ratio": model_std / ersst_std,
        "abs_bias": abs(model_std - ersst_std),
    }

model_ranking = sorted(models, key=lambda m: amplitude[m]["abs_bias"])

observed_peak_month = stats["ERSSTv5"]["peak_month"]
observed_peak_value = stats["ERSSTv5"]["peak_value"]
observed_ratio = stats["ERSSTv5"]["ratio"]

peak_comparison = {}
for model in models:
    s = stats[model]
    peak_comparison[model] = {
        "matches_peak_month": bool(s["peak_month"] == observed_peak_month),
        "peak_month": s["peak_month"],
        "peak_value": s["peak_value"],
        "peak_value_diff": s["peak_value"] - observed_peak_value,
        "ratio": s["ratio"],
        "ratio_diff": s["ratio"] - observed_ratio,
    }

print("Nino 3.4 SST index: 5S-5N, 170W-120W (190E-240E), 1950-01 to 2014-12")
print("Anomaly baseline: each source's own 1950-2014 monthly climatology, calendar-month matched.")
print("Standard deviations use ddof=0 (population).")
print()

print(
    f"{'Source':<16} {'std_all':>8} {'peak_m':>6} {'peak_val':>9} "
    f"{'min_m':>6} {'min_val':>8} {'ratio':>7} {'n>|1|':>6} {'pct':>6}"
)
for src in sources:
    s = stats[src]
    print(
        f"{src:<16} {s['std_all']:8.2f} {s['peak_month']:6d} {s['peak_value']:9.2f} "
        f"{s['min_month']:6d} {s['min_value']:8.2f} {s['ratio']:7.2f} "
        f"{s['n_exceed']:6d} {s['pct']:6.1f}"
    )

print("\nMonthly standard deviations (degC):")
print(f"{'Source':<16}" + "".join([f"{m:>6d}" for m in range(1, 13)]))
for src in sources:
    print(f"{src:<16}" + "".join([f"{v:6.2f}" for v in stats[src]["std_month"]]))

print("\nMonthly climatology used as baseline (degC):")
print(f"{'Source':<16}" + "".join([f"{m:>6d}" for m in range(1, 13)]))
for src in sources:
    print(f"{src:<16}" + "".join([f"{v:6.2f}" for v in stats[src]["climatology"]]))

print("\nAmplitude versus ERSSTv5:")
for model in models:
    a = amplitude[model]
    print(
        f"{model:<16} std={a['model_std']:.2f}  bias={a['bias']:+.2f}  "
        f"ratio={a['ratio']:.2f}  |bias|={a['abs_bias']:.2f}"
    )
print("Ranking (smallest to largest |amplitude bias|):", " < ".join(model_ranking))

print("\nSeasonal phase-locking peak comparison with ERSSTv5:")
for model in models:
    p = peak_comparison[model]
    print(
        f"{model:<16} peak_month={p['peak_month']} "
        f"matches_obs_peak={p['matches_peak_month']} "
        f"peak_value={p['peak_value']:.2f} (diff {p['peak_value_diff']:+.2f}) "
        f"ratio={p['ratio']:.2f} (diff {p['ratio_diff']:+.2f})"
    )

print("\nVariants parsed from actual filenames:")
for model in models:
    print(f"{model:<16} {variants[model]}")
