import os
import glob
from collections import Counter

WOA_ROOT = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
WOA_MASKS = "/data1/lorraine-20419/ocean_data/WOA23/masks"
CMIP_ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"
RHO = 1025.0

MODELS = {
    "CESM2": dict(y="nlat", x="nlon", lat="lat", lon="lon", depth="lev",
                  depth_units="centimeters", depth_scale=0.01),
    "GFDL-ESM4": dict(y="lat", x="lon", lat="lat", lon="lon", depth="lev",
                      depth_units="m", depth_scale=1.0),
    "MPI-ESM1-2-LR": dict(y="j", x="i", lat="latitude", lon="longitude",
                          depth="lev", depth_units="m", depth_scale=1.0),
}
SOURCES = ["WOA23", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

REG_LAT = -89.5 + np.arange(180, dtype=float)
REG_LON = -179.5 + np.arange(360, dtype=float)


def _find_dim(ds, candidates, std_name=None):
    for c in candidates:
        if c in ds.dims:
            return c
    if std_name is not None:
        for c in ds.dims:
            if ds[c].attrs.get("standard_name") == std_name:
                return c
    raise KeyError(f"None of {candidates} found in dimensions {list(ds.dims)}")


def load_woa(prefix, subdir):
    fn = os.path.join(WOA_ROOT, subdir, f"woa23_all_{prefix}00_01.nc")
    ds = xr.open_dataset(fn, decode_times=False)
    vname = f"{prefix}_an"
    if vname not in ds.data_vars:
        cand = [v for v in ds.data_vars
                if v.lower().startswith(prefix.lower()) and v.lower().endswith("_an")]
        vname = cand[0] if cand else list(ds.data_vars)[0]
    da = ds[vname]
    for tname in ["time", "t"]:
        if tname in da.dims:
            da = da.isel({tname: 0})
    da = da.squeeze(drop=True)
    dname = _find_dim(ds, ["depth", "deptht", "z", "lev"], "depth")
    latname = _find_dim(ds, ["lat", "latitude"], "latitude")
    lonname = _find_dim(ds, ["lon", "longitude"], "longitude")
    da = da.transpose(dname, latname, lonname)
    out = da.compute()
    ds.close()
    return out


def load_basin_mask():
    path = os.path.join(WOA_MASKS, "basinmask_01.msk")
    rows = []
    with open(path, "r", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s[0] in "#!%":
                continue
            try:
                arr = np.fromstring(s, sep=" ")
            except Exception:
                continue
            if arr.size >= 104 and np.all(np.isfinite(arr)):
                rows.append(arr)
    if not rows:
        raise RuntimeError("Could not parse basinmask_01.msk")
    lens = Counter(len(a) for a in rows)
    L = lens.most_common(1)[0][0]
    rows = [a for a in rows if len(a) == L]
    arr = np.vstack(rows)

    # Detect lat/lon columns among first few columns.
    ncol = arr.shape[1]
    best = None
    for i in range(min(4, ncol)):
        for j in range(min(4, ncol)):
            if i == j:
                continue
            lat = arr[:, i]
            lon = arr[:, j]
            if (np.nanmin(lat) >= -90 and np.nanmax(lat) <= 90 and
                    np.nanmax(np.abs(lon)) <= 400):
                score = np.unique(np.round(lat)).size + np.unique(np.round(lon)).size
                if best is None or score > best[0]:
                    best = (score, i, j)
    if best is not None:
        _, latcol, loncol = best
    else:
        latcol, loncol = 0, 1

    lat = arr[:, latcol].astype(float)
    lon = arr[:, loncol].astype(float)
    vals = arr[:, -102:].astype(float)

    basin = np.zeros((102, 180, 360), dtype=np.int16)
    lon180 = ((lon + 180.0) % 360.0) - 180.0
    ilat = np.clip(np.rint(lat - (-89.5)).astype(int), 0, 179)
    ilon = np.rint(lon180 - (-179.5)).astype(int) % 360

    for k in range(arr.shape[0]):
        v = vals[k]
        good = np.isfinite(v) & (v >= 1) & (v <= 5)
        if np.any(good):
            basin[good, ilat[k], ilon[k]] = v[good].astype(np.int16)
    return basin


def woa_layer_volume(depth_all, d_idx, z0=200.0, z1=1000.0):
    d = np.asarray(depth_all, dtype=float)
    edges = np.empty(len(d) + 1, dtype=float)
    edges[1:-1] = 0.5 * (d[:-1] + d[1:])
    edges[0] = max(0.0, d[0] - (d[1] - d[0]) / 2.0)
    edges[-1] = d[-1] + (d[-1] - d[-2]) / 2.0
    thick = np.minimum(edges[1:], z1) - np.maximum(edges[:-1], z0)
    thick = np.maximum(thick, 0.0)
    return thick[d_idx]


def load_cmip_field(model, var):
    info = MODELS[model]
    files = sorted(glob.glob(os.path.join(CMIP_ROOT, model, var, "*.nc")))
    if not files:
        raise FileNotFoundError(f"No files for {model}/{var}")
    ds = xr.open_mfdataset(files, combine="by_coords", use_cftime=True,
                           decode_times=True, chunks={"time": 12})
    if "time" in ds.dims:
        yrs = ds.time.dt.year.values
        idx = np.where((yrs >= 1971) & (yrs <= 2000))[0]
        if len(idx) == 0:
            raise ValueError(f"No 1971-2000 times for {model}/{var}")
        ds = ds.isel(time=idx)

    dname = info["depth"]
    lev = ds[dname].values.astype(float)
    if info["depth_units"] == "centimeters":
        lev_m = lev * info["depth_scale"]
    else:
        lev_m = lev
    lev_idx = np.where((lev_m >= 200.0) & (lev_m <= 1000.0))[0]
    if len(lev_idx) == 0:
        raise ValueError(f"No 200-1000 m levels for {model}/{var}")

    da = ds[var].isel({dname: lev_idx})
    if "time" in da.dims:
        da = da.mean("time")
    da = da.squeeze(drop=True)
    da = da.transpose(dname, info["y"], info["x"])
    units = ds[var].attrs.get("units", "")

    latv = ds[info["lat"]].values
    lonv = ds[info["lon"]].values
    if latv.ndim == 1 and lonv.ndim == 1:
        lon2d, lat2d = np.meshgrid(lonv, latv)
    else:
        lat2d, lon2d = np.broadcast_arrays(latv, lonv)

    da = da.compute()
    ds.close()
    return da, lat2d, lon2d, lev_m[lev_idx], units


def load_cmip_volcello(model):
    info = MODELS[model]
    files = sorted(glob.glob(os.path.join(CMIP_ROOT, model, "volcello", "*.nc")))
    if not files:
        raise FileNotFoundError(f"No volcello files for {model}")
    ds = xr.open_mfdataset(files, combine="by_coords", use_cftime=True,
                           decode_times=True)

    dname = info["depth"]
    lev = ds[dname].values.astype(float)
    if info["depth_units"] == "centimeters":
        lev_m = lev * info["depth_scale"]
    else:
        lev_m = lev
    lev_idx = np.where((lev_m >= 200.0) & (lev_m <= 1000.0))[0]

    vol = ds["volcello"]
    if "time" in vol.dims:
        vol = vol.isel(time=0)
    vol = vol.isel({dname: lev_idx}).squeeze(drop=True)
    vol = vol.transpose(dname, info["y"], info["x"])
    vol = vol.compute()
    ds.close()
    return vol


def to_umol_per_kg(da, units):
    u = (units or "").lower().replace(" ", "")
    arr = da.values if hasattr(da, "values") else np.asarray(da)
    if "micromol" in u or "umol" in u:
        return arr
    if "mmol" in u:
        return arr * 1e3 / RHO
    # CMIP6 Omon no3/po4 are mol m-3.
    return arr * 1e6 / RHO


def regrid_to_regular(lat2d, lon2d, data2d):
    lat2d, lon2d = np.broadcast_arrays(np.asarray(lat2d, dtype=float),
                                       np.asarray(lon2d, dtype=float))
    data2d = np.asarray(data2d, dtype=float)
    valid = np.isfinite(data2d) & np.isfinite(lat2d) & np.isfinite(lon2d)
    if not np.any(valid):
        return np.full((180, 360), np.nan)

    lat_v = lat2d[valid]
    lon_v = lon2d[valid]
    dat_v = data2d[valid]

    lon_v = ((lon_v + 180.0) % 360.0) - 180.0
    ilat = np.clip(np.rint(lat_v - REG_LAT[0]).astype(int), 0, 179)
    ilon = np.rint(lon_v - REG_LON[0]).astype(int) % 360
    idx = ilat * 360 + ilon

    sums = np.bincount(idx, weights=dat_v, minlength=180 * 360)
    cnts = np.bincount(idx, minlength=180 * 360)
    grid = np.full(180 * 360, np.nan)
    m = cnts > 0
    grid[m] = sums[m] / cnts[m]
    return grid.reshape(180, 360)


def model_basin_array(lat2d, lon2d, lev_m, basin3d, woa_depth):
    lat2d, lon2d = np.broadcast_arrays(np.asarray(lat2d, dtype=float),
                                       np.asarray(lon2d, dtype=float))
    lon180 = ((lon2d + 180.0) % 360.0) - 180.0
    ilat = np.clip(np.rint(lat2d - (-89.5)).astype(int), 0, 179)
    ilon = np.rint(lon180 - (-179.5)).astype(int) % 360

    out = np.zeros((len(lev_m),) + lat2d.shape, dtype=np.int16)
    woa_depth = np.asarray(woa_depth, dtype=float)
    for k, z in enumerate(lev_m):
        di = int(np.argmin(np.abs(woa_depth - z)))
        out[k] = basin3d[di][ilat, ilon]
    return out


def analyze_source(name, no3, po4, vol, basin, lat2d, lon2d):
    no3 = np.asarray(no3, dtype=float)
    po4 = np.asarray(po4, dtype=float)
    vol = np.asarray(vol, dtype=float)
    basin = np.asarray(basin)

    valid = (np.isfinite(no3) & np.isfinite(po4) & (po4 > 0) &
             np.isfinite(vol) & (vol > 0))
    ratio = np.full(no3.shape, np.nan, dtype=float)
    ratio[valid] = no3[valid] / po4[valid]

    x = po4[valid]
    y = no3[valid]
    r = ratio[valid]
    w = vol[valid]
    b = basin[valid]

    mean_ratio_w = np.sum(r * w) / np.sum(w)
    mean_ratio_u = np.mean(r)
    ratio_means_w = np.sum(y * w) / np.sum(x * w)
    ratio_means_u = np.mean(y) / np.mean(x)

    a, bint = np.polyfit(x, y, 1)
    yhat = a * x + bint
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    pct_14_18 = 100.0 * np.sum(w[(r > 14) & (r < 18)]) / np.sum(w)
    pct_below12_global = 100.0 * np.sum(w[r < 12]) / np.sum(w)

    basin_means = {}
    basin_pct_below12 = {}
    for code, bname in [(1, "Atlantic"), (2, "Pacific"), (3, "Indian")]:
        m = (b == code)
        if np.any(m):
            basin_means[bname] = np.sum(r[m] * w[m]) / np.sum(w[m])
            basin_pct_below12[bname] = (
                100.0 * np.sum(w[m & (r < 12)]) / np.sum(w[m])
            )
        else:
            basin_means[bname] = np.nan
            basin_pct_below12[bname] = np.nan

    num = np.nansum(np.where(valid, ratio * vol, 0.0), axis=0)
    den = np.nansum(np.where(valid, vol, 0.0), axis=0)
    col = np.where(den > 0, num / den - 16.0, np.nan)
    map_np = regrid_to_regular(lat2d, lon2d, col)

    rng = np.random.default_rng(0)
    n = len(x)
    if n > 20000:
        ii = rng.choice(n, size=20000, replace=False)
    else:
        ii = np.arange(n)

    return {
        "mean_np": float(mean_ratio_w),
        "mean_np_unweighted": float(mean_ratio_u),
        "ratio_of_mean_concentrations": float(ratio_means_w),
        "ratio_of_mean_concentrations_unweighted": float(ratio_means_u),
        "a": float(a),
        "b": float(bint),
        "r2": float(r2),
        "pct_14_18": float(pct_14_18),
        "pct_below12_global": float(pct_below12_global),
        "basin_mean_np": {k: float(v) for k, v in basin_means.items()},
        "basin_pct_below12": {k: float(v) for k, v in basin_pct_below12.items()},
        "pct_below12_pacific": float(basin_pct_below12["Pacific"]),
        "pct_below12_indian": float(basin_pct_below12["Indian"]),
        "scatter_x": x[ii].copy(),
        "scatter_y": y[ii].copy(),
        "map_np": map_np,
        "map_lat": REG_LAT.copy(),
        "map_lon": REG_LON.copy(),
        "n_cells": int(n),
        "volume_total": float(np.sum(w)),
    }


basin3d = load_basin_mask()

# WOA23 observed annual objectively analysed nutrients.
woa_no3_da = load_woa("n", "nitrate/all_time_periods")
woa_po4_da = load_woa("p", "phosphate/all_time_periods")
woa_no3_da, woa_po4_da = xr.align(woa_no3_da, woa_po4_da, join="inner")

woa_depth = woa_no3_da["depth"].values.astype(float)
woa_d_idx = np.where((woa_depth >= 200.0) & (woa_depth <= 1000.0))[0]
woa_lat = woa_no3_da["lat"].values.astype(float)
woa_lon = woa_no3_da["lon"].values.astype(float)

no3_woa = woa_no3_da.isel(depth=woa_d_idx).values
po4_woa = woa_po4_da.isel(depth=woa_d_idx).values
basin_woa = basin3d[woa_d_idx, :, :]

thick = woa_layer_volume(woa_depth, woa_d_idx)
vol_woa = np.broadcast_to(
    thick[:, None, None] * np.cos(np.deg2rad(woa_lat))[None, :, None],
    no3_woa.shape
).copy()
lat2d_woa, lon2d_woa = np.meshgrid(woa_lon, woa_lat)

results = {}
results["WOA23"] = analyze_source(
    "WOA23", no3_woa, po4_woa, vol_woa, basin_woa, lat2d_woa, lon2d_woa
)

# CMIP6 historical model nutrients, 1971-2000, each model's own volume.
for model in ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]:
    no3_da, lat2d, lon2d, lev_m, units_no3 = load_cmip_field(model, "no3")
    po4_da, _, _, _, units_po4 = load_cmip_field(model, "po4")
    vol_da = load_cmip_volcello(model)

    dname = MODELS[model]["depth"]
    yname = MODELS[model]["y"]
    xname = MODELS[model]["x"]

    no3_da, po4_da, vol_da = xr.align(no3_da, po4_da, vol_da, join="inner")
    no3_da = no3_da.transpose(dname, yname, xname)
    po4_da = po4_da.transpose(dname, yname, xname)
    vol_da = vol_da.transpose(dname, yname, xname)

    no3 = to_umol_per_kg(no3_da, units_no3)
    po4 = to_umol_per_kg(po4_da, units_po4)
    vol = vol_da.values

    basin_model = model_basin_array(lat2d, lon2d, lev_m, basin3d, woa_depth)
    results[model] = analyze_source(
        model, no3, po4, vol, basin_model, lat2d, lon2d
    )

summary_rows = []
for s in SOURCES:
    d = results[s]
    summary_rows.append({
        "source": s,
        "mean_NP_cell_volume_weighted": d["mean_np"],
        "mean_NP_cell_unweighted": d["mean_np_unweighted"],
        "NO3_mean_over_PO4_mean_volume_weighted": d["ratio_of_mean_concentrations"],
        "NO3_mean_over_PO4_mean_unweighted": d["ratio_of_mean_concentrations_unweighted"],
        "a": d["a"],
        "b": d["b"],
        "R2": d["r2"],
        "pct_14_18": d["pct_14_18"],
        "pct_below12_global": d["pct_below12_global"],
        "pct_below12_Pacific": d["pct_below12_pacific"],
        "pct_below12_Indian": d["pct_below12_indian"],
        "NP_Atlantic": d["basin_mean_np"]["Atlantic"],
        "NP_Pacific": d["basin_mean_np"]["Pacific"],
        "NP_Indian": d["basin_mean_np"]["Indian"],
    })

summary = pd.DataFrame(summary_rows).set_index("source")

redfield_position = {
    s: ("above" if results[s]["mean_np"] > 16.0 else
        "below" if results[s]["mean_np"] < 16.0 else "at")
    for s in SOURCES
}

low_np_reproduction = {}
for s in SOURCES[1:]:
    low_np_reproduction[s] = {
        "global_pct_below12_bias_vs_WOA": (
            results[s]["pct_below12_global"] -
            results["WOA23"]["pct_below12_global"]
        ),
        "Pacific_pct_below12_bias_vs_WOA": (
            results[s]["pct_below12_pacific"] -
            results["WOA23"]["pct_below12_pacific"]
        ),
        "Indian_pct_below12_bias_vs_WOA": (
            results[s]["pct_below12_indian"] -
            results["WOA23"]["pct_below12_indian"]
        ),
    }

print(summary)
print("\nRedfield position:")
print(redfield_position)
print("\nLow-N:P volume biases vs WOA23:")
print(low_np_reproduction)
