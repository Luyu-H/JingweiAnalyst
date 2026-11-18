import os
import glob
import re
import warnings
import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import griddata

warnings.filterwarnings("ignore")

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
CMIP6_HIST_ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"
ERSST_FILE = "/data1/lorraine-20419/ocean_data/ERSSTv5/sst.mnmean.nc"
CPC_ROOT = "/data1/lorraine-20419/ocean_data/CPC"

TSTART = "1950-01-01"
TEND = "2014-12-31"
YEAR_MIN, YEAR_MAX = 1950, 2014

# Target 1-degree grid over the PDO domain 20N-70N, 110E-100W.
target_lat = np.arange(20.0, 71.0, 1.0)
target_lon = np.arange(110.0, 261.0, 1.0)

# Area weights (cosine latitude).  Used for EOF weighting and weighted RMS.
coslat = np.cos(np.deg2rad(target_lat))
w2d = np.broadcast_to(coslat[:, None], (len(target_lat), len(target_lon))).copy()

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def select_1950_2014(da):
    """Select 1950-01 through 2014-12 using year/month, robust to cftime."""
    yrs = da["time"].dt.year.values
    mask = (yrs >= YEAR_MIN) & (yrs <= YEAR_MAX)
    return da.isel(time=mask)


def monthly_anomalies(da):
    """Subtract the 1950-2014 monthly climatology."""
    da = select_1950_2014(da)
    da = da.sortby("time")
    clim = da.groupby("time.month").mean("time")
    anom = da.groupby("time.month") - clim
    return anom


def regrid_regular(da, var, lat_name="lat", lon_name="lon"):
    """Bilinear regrid a regular lat/lon field onto the target 1-degree grid."""
    if isinstance(da, xr.Dataset):
        da = da[var]

    da = select_1950_2014(da)
    da = da.sortby(lat_name)

    lon_vals = da[lon_name].values
    if np.nanmin(lon_vals) < 0 or (np.nanmax(lon_vals) <= 180 and target_lon.max() > 180):
        da = da.assign_coords({lon_name: ((da[lon_name] + 360.0) % 360.0)}).sortby(lon_name)

    da = da.compute()

    out = da.interp({lat_name: target_lat, lon_name: target_lon}, method="linear")

    # Interpolate a 0/1 ocean mask so that land cells do not get filled by the
    # bilinear SST interpolation.
    valid = (da.notnull().mean("time") > 0.5).astype(float)
    mask_out = valid.interp({lat_name: target_lat, lon_name: target_lon}, method="linear")
    out = out.where(mask_out > 0.5)

    if lat_name != "lat" or lon_name != "lon":
        out = out.rename({lat_name: "lat", lon_name: "lon"})

    return out


def regrid_curvilinear(da, var, lat2d, lon2d):
    """Bilinear regrid a curvilinear field onto the target 1-degree grid."""
    if isinstance(da, xr.Dataset):
        da = da[var]

    da = select_1950_2014(da)
    times = da["time"].values
    ntime = len(times)

    tlon, tlat = np.meshgrid(target_lon, target_lat)

    lon_use = lon2d.copy()
    if np.nanmin(lon_use) < 0:
        lon_use = lon_use % 360.0

    # Ocean mask fraction over time, using all grid points (ocean and land).
    mask_frac = da.notnull().mean("time").values.astype(float)

    good_pts = np.isfinite(lat2d) & np.isfinite(lon_use)
    # Work on the PDO domain plus a buffer; this avoids global griddata calls.
    dom = (
        good_pts
        & (lat2d >= 10.0)
        & (lat2d <= 80.0)
        & (lon_use >= 100.0)
        & (lon_use <= 270.0)
    )

    if dom.sum() > 4:
        mask_out = griddata(
            (lon_use[dom], lat2d[dom]),
            mask_frac[dom],
            (tlon, tlat),
            method="linear",
        )
    else:
        mask_out = np.zeros_like(tlat)
    mask_out = np.nan_to_num(mask_out, nan=0.0)

    out = np.full((ntime, len(target_lat), len(target_lon)), np.nan, dtype=float)

    for i in range(ntime):
        vals = da.isel(time=i).values.astype(float)
        valid = np.isfinite(vals) & dom
        if valid.sum() < 4:
            continue
        out[i] = griddata(
            (lon_use[valid], lat2d[valid]),
            vals[valid],
            (tlon, tlat),
            method="linear",
        )

    out[:, mask_out < 0.5] = np.nan

    return xr.DataArray(
        out,
        dims=("time", "lat", "lon"),
        coords={"time": times, "lat": target_lat, "lon": target_lon},
    )


def pdo_eof(anom_da, label=""):
    """Leading EOF/PDO analysis on monthly SST anomalies."""
    da = anom_da.transpose("time", "lat", "lon").compute()
    X = da.values.astype(float)
    ntime, nlat, nlon = X.shape

    # Cells with a complete finite time series.
    valid = np.isfinite(X).all(axis=0)
    Xv = X[:, valid]
    w = w2d[valid]

    # Area-weighted EOF: multiply by sqrt(cos lat).
    Xw = Xv * np.sqrt(w)[None, :]
    Xw = Xw - Xw.mean(axis=0, keepdims=True)

    U, s, Vt = np.linalg.svd(Xw, full_matrices=False)

    var_explained = float((s[0] ** 2 / np.sum(s ** 2)) * 100.0)

    pc = U[:, 0] * s[0]
    index = pc / pc.std(ddof=1)

    # Regression pattern in degrees C per standard deviation of the index.
    pat = np.full((nlat, nlon), np.nan, dtype=float)
    pat[valid] = (Xv * index[:, None]).mean(axis=0)

    # Weighted RMS of the regression pattern = mode amplitude in degC per SD.
    amplitude = float(np.sqrt(np.average(pat[valid] ** 2, weights=w)))

    # Lag-1-year (12-month) autocorrelation of the standardized index.
    autocorr = float(np.corrcoef(index[:-12], index[12:])[0, 1])

    yrs = da["time"].dt.year.values
    mos = da["time"].dt.month.values
    time_frac = yrs + (mos - 0.5) / 12.0

    return {
        "label": label,
        "pattern": pat,
        "index": index,
        "time": da["time"].values,
        "time_frac": time_frac,
        "var_explained": var_explained,
        "amplitude": amplitude,
        "autocorr": autocorr,
        "pc_std": float(pc.std(ddof=1)),
        "valid": valid,
    }


def spatial_corr(a, b, weights=None):
    """Pearson correlation over common finite cells, optionally area-weighted."""
    valid = np.isfinite(a) & np.isfinite(b)
    if valid.sum() < 3:
        return np.nan

    aa = a[valid].astype(float)
    bb = b[valid].astype(float)

    if weights is None:
        return float(np.corrcoef(aa, bb)[0, 1])

    ww = weights[valid].astype(float)
    ma = np.average(aa, weights=ww)
    mb = np.average(bb, weights=ww)
    cov = np.average((aa - ma) * (bb - mb), weights=ww)
    va = np.average((aa - ma) ** 2, weights=ww)
    vb = np.average((bb - mb) ** 2, weights=ww)
    return float(cov / np.sqrt(va * vb))


def find_transitions(time, index, window=108, merge_years=3):
    """
    Find sign changes in the centred 9-year (108-month) running mean.
    Crossings less than 3 years after the last counted transition are merged
    into it, keeping the earlier year.
    """
    s = pd.Series(index, index=pd.to_datetime(time)).sort_index()
    smooth = s.rolling(window=window, center=True, min_periods=window).mean()

    vals = smooth.values
    dates = smooth.index
    transitions = []
    last_counted = None

    for i in range(1, len(vals)):
        if not np.isfinite(vals[i - 1]) or not np.isfinite(vals[i]):
            continue
        if vals[i - 1] * vals[i] < 0.0:
            yr = int(dates[i].year)
            if last_counted is None or (yr - last_counted) >= merge_years:
                transitions.append(yr)
                last_counted = yr
            # else: merge into last counted transition, keep earlier year

    # Remove accidental duplicates while preserving order.
    transitions = list(dict.fromkeys(transitions))
    return smooth, transitions


def load_pdo_index():
    """Read the published NCEI/CPC PDO monthly index file."""
    path = os.path.join(CPC_ROOT, "ersst.v6.pdo.dat")
    rows = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("%") or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 13:
                continue
            if parts[0].lower() == "year":
                continue
            try:
                yr = int(parts[0])
            except ValueError:
                continue
            vals = []
            for p in parts[1:13]:
                try:
                    vals.append(float(p))
                except ValueError:
                    vals.append(np.nan)
            rows.append([yr] + vals)

    df = pd.DataFrame(rows, columns=["year"] + list(range(1, 13)))
    df = df.replace(99.99, np.nan)
    long = df.melt(id_vars="year", var_name="month", value_name="pdo")
    long["time"] = pd.to_datetime(
        dict(year=long["year"].astype(int), month=long["month"].astype(int), day=1)
    )
    return long.set_index("time")["pdo"].sort_index()


def open_model_tos(model):
    """Open the CMIP6 historical tos files for one model."""
    var = "tos"
    exp = "historical"
    files = glob.glob(os.path.join(CMIP6_HIST_ROOT, model, var, "*.nc"))
    files = [
        f
        for f in files
        if os.path.basename(f).startswith(var + "_")
        and ("_" + exp + "_") in os.path.basename(f)
    ]
    if not files:
        raise FileNotFoundError(f"No {var} files found for {model} in {CMIP6_HIST_ROOT}")

    by_variant = {}
    for f in files:
        m = re.search(r"_(r\d+i\d+p\d+f\d+)_", os.path.basename(f))
        if m:
            by_variant.setdefault(m.group(1), []).append(f)

    if not by_variant:
        raise FileNotFoundError(f"Could not parse variant labels for {model} {var}")

    # Prefer r1i1p1f1 if present, otherwise the variant with the most files.
    variant = (
        "r1i1p1f1"
        if "r1i1p1f1" in by_variant
        else max(by_variant, key=lambda k: len(by_variant[k]))
    )
    files = sorted(by_variant[variant])

    ds = xr.open_mfdataset(
        files,
        combine="by_coords",
        use_cftime=True,
        decode_times=True,
        chunks={"time": 120},
    )
    da = ds[var]
    if "lev" in da.dims:
        da = da.isel(lev=0)

    units = str(da.attrs.get("units", "")).lower()
    if units in ("k", "kelvin", "degrees_k", "degree_kelvin"):
        da = da - 273.15

    return da, ds, variant


# -----------------------------------------------------------------------------
# ERSSTv5 observed SST
# -----------------------------------------------------------------------------
ersst_ds = xr.open_dataset(ERSST_FILE, decode_times=True)
ersst_reg = regrid_regular(ersst_ds["sst"], "sst", lat_name="lat", lon_name="lon")
ersst_anom = monthly_anomalies(ersst_reg)
ersst_res = pdo_eof(ersst_anom, "ERSSTv5")

# Published PDO index, 1950-2014.
pub_pdo = load_pdo_index()
pub_1950_2014 = pub_pdo.loc["1950-01-01":"2014-12-01"].copy()

# Orient the ERSST-derived EOF to the published PDO series.
obs_series = pd.Series(ersst_res["index"], index=pd.to_datetime(ersst_res["time"]))
pub_aligned = pub_1950_2014.reindex(obs_series.index)
corr_obs_pub = float(obs_series.corr(pub_aligned))
if corr_obs_pub < 0:
    ersst_res["index"] *= -1.0
    ersst_res["pattern"] *= -1.0
    corr_obs_pub = -corr_obs_pub
ersst_res["corr_published"] = corr_obs_pub

# Standardized published PDO for plotting.
pub_std = (pub_1950_2014 - pub_1950_2014.mean()) / pub_1950_2014.std()
pub_time_frac = pub_std.index.year + (pub_std.index.month - 0.5) / 12.0

# Observed PDO phase transitions.
smooth_ersst, transition_years = find_transitions(
    ersst_res["time"], ersst_res["index"], window=108, merge_years=3
)
smooth_time_frac = smooth_ersst.index.year + (smooth_ersst.index.month - 0.5) / 12.0

# Documented Pacific climate shifts around the PDO transition years.
documented_shifts = [1976, 1977, 1988, 1989, 1998, 1999]
documented_transition_matches = [
    yr for yr in transition_years
    if min(abs(yr - d) for d in documented_shifts) <= 1
]

# -----------------------------------------------------------------------------
# CMIP6 historical tos for the three models
# -----------------------------------------------------------------------------
model_results = {}
model_variants = {}

for model in ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]:
    da, ds, variant = open_model_tos(model)
    model_variants[model] = variant

    if model == "CESM2":
        lat2d = ds["lat"].values
        lon2d = ds["lon"].values
        reg = regrid_curvilinear(da, "tos", lat2d, lon2d)
    elif model == "MPI-ESM1-2-LR":
        lat2d = ds["latitude"].values
        lon2d = ds["longitude"].values
        reg = regrid_curvilinear(da, "tos", lat2d, lon2d)
    else:  # GFDL-ESM4, regular grid
        reg = regrid_regular(da, "tos", lat_name="lat", lon_name="lon")

    anom = monthly_anomalies(reg)
    res = pdo_eof(anom, model)

    # Orient each model's EOF to the observed ERSSTv5 PDO pattern.
    pc_corr = spatial_corr(ersst_res["pattern"], res["pattern"])
    if pc_corr < 0:
        res["index"] *= -1.0
        res["pattern"] *= -1.0
        pc_corr = -pc_corr

    res["pattern_corr_vs_ERSSTv5"] = float(pc_corr)
    res["pattern_corr_area_weighted_vs_ERSSTv5"] = spatial_corr(
        ersst_res["pattern"], res["pattern"], weights=w2d
    )
    res["variant"] = variant
    model_results[model] = res

    ds.close()

ersst_ds.close()

# -----------------------------------------------------------------------------
# Assemble results
# -----------------------------------------------------------------------------
all_results = {"ERSSTv5": ersst_res}
all_results.update(model_results)

summary_rows = []
for name, res in all_results.items():
    row = {
        "source": name,
        "var_explained_pct": round(float(res["var_explained"]), 1),
        "amplitude_C_per_sd": round(float(res["amplitude"]), 3),
        "lag12_autocorr": round(float(res["autocorr"]), 3),
    }
    if name == "ERSSTv5":
        row["pattern_corr_vs_ERSSTv5"] = 1.0
    else:
        row["pattern_corr_vs_ERSSTv5"] = round(
            float(res["pattern_corr_vs_ERSSTv5"]), 3
        )
    summary_rows.append(row)

summary = pd.DataFrame(summary_rows)

obs_pub_corr = round(float(corr_obs_pub), 3)

model_names = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
best_model = max(
    model_names, key=lambda m: model_results[m]["pattern_corr_vs_ERSSTv5"]
)
best_pattern_corr = round(float(model_results[best_model]["pattern_corr_vs_ERSSTv5"]), 3)

amplitude_diff = {
    m: round(float(model_results[m]["amplitude"] - ersst_res["amplitude"]), 3)
    for m in model_names
}
autocorr_diff = {
    m: round(float(model_results[m]["autocorr"] - ersst_res["autocorr"]), 3)
    for m in model_names
}

model_consistency = {
    "best_model": best_model,
    "best_pattern_corr": best_pattern_corr,
    "amplitude_difference_C_per_sd": amplitude_diff,
    "lag12_autocorr_difference": autocorr_diff,
    "model_variants": model_variants,
}

# -----------------------------------------------------------------------------
# Values kept for the answer / plotting
# -----------------------------------------------------------------------------
var_explained = {
    k: round(float(v["var_explained"]), 1) for k, v in all_results.items()
}
amplitudes = {
    k: round(float(v["amplitude"]), 3) for k, v in all_results.items()
}
lag12_autocorr = {
    k: round(float(v["autocorr"]), 3) for k, v in all_results.items()
}
pattern_correlations = {
    k: (
        1.0
        if k == "ERSSTv5"
        else round(float(v["pattern_corr_vs_ERSSTv5"]), 3)
    )
    for k, v in all_results.items()
}

print("Summary")
print(summary.to_string(index=False))
print()
print("ERSSTv5-derived vs published PDO correlation:", obs_pub_corr)
print("PDO transition years:", transition_years)
print("Matches to documented Pacific shifts:", documented_transition_matches)
print("Best model by pattern correlation:", best_model, best_pattern_corr)
print("Amplitude differences (model - ERSSTv5):", amplitude_diff)
print("Lag-12 autocorrelation differences (model - ERSSTv5):", autocorr_diff)
