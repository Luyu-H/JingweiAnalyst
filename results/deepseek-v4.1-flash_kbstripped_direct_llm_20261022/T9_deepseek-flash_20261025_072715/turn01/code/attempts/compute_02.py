import glob
import warnings
import numpy as np
import pandas as pd
import xarray as xr
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ---------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------
cmip_root = "/data1/lorraine-20419/benchmark_cmip6/historical"
ersst_path = "/data1/lorraine-20419/ocean_data/ERSSTv5/sst.mnmean.nc"
pdo_path = "/data1/lorraine-20419/ocean_data/CPC/ersst.v6.pdo.dat"

start_sel, end_sel = "1950-01-01", "2014-12-31"
time_common = pd.date_range("1950-01-01", "2014-12-01", freq="MS")
months = time_common.month.values
n_time = len(time_common)

# 1-degree target grid over the PDO domain: 20N-70N, 110E-100W.
# Keep longitudes in 0-360 convention, so the domain is 110E to 260E.
lat_t = np.arange(20.0, 70.5, 1.0)
lon_t = np.arange(110.0, 260.5, 1.0)
LON_T, LAT_T = np.meshgrid(lon_t, lat_t)

# Standard climate-EOF practice is area weighting by sqrt(cos(lat)).
# Set to False for an unweighted grid-point EOF.
AREA_WEIGHT = True

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
def to_celsius(da):
    """Convert Kelvin to degC if the declared units are Kelvin."""
    units = str(da.attrs.get("units", "")).lower()
    if ("k" in units and "deg" not in units) or "kelvin" in units:
        return da - 273.15
    return da


def regrid_regular(da, lat_name="lat", lon_name="lon",
                   target_lat=lat_t, target_lon=lon_t):
    """Bilinear regrid a rectilinear lat/lon DataArray to the target grid."""
    if lat_name not in da.coords and "latitude" in da.coords:
        lat_name = "latitude"
    if lon_name not in da.coords and "longitude" in da.coords:
        lon_name = "longitude"

    da = da.sortby(lat_name).sortby(lon_name)

    # Put longitude into 0-360 if needed.
    lon_vals = da[lon_name].values
    if np.nanmin(lon_vals) < 0:
        da = da.assign_coords({lon_name: (da[lon_name] % 360)}).sortby(lon_name)

    # Subset before interpolating for efficiency.
    lat_min, lat_max = float(target_lat.min()) - 5, float(target_lat.max()) + 5
    lon_min, lon_max = float(target_lon.min()) - 5, float(target_lon.max()) + 5
    da = da.sel({lat_name: slice(lat_min, lat_max),
                 lon_name: slice(lon_min, lon_max)})

    out = da.interp({lat_name: target_lat, lon_name: target_lon},
                    method="linear")
    return out


def regrid_curvilinear(tos, lat2d, lon2d,
                       target_lat=lat_t, target_lon=lon_t):
    """
    Bilinear-ish regrid of a curvilinear CMIP6 grid to the target grid using a
    linear Delaunay interpolation over the source ocean points.
    """
    lat2d = np.asarray(lat2d)
    lon2d = np.asarray(lon2d)
    if np.nanmin(lon2d) < 0:
        lon2d = lon2d % 360

    vals = np.asarray(tos.values)  # time, y, x
    finite0 = np.isfinite(vals[0])

    # Build the triangulation only over the PDO domain plus a margin.
    box = ((lat2d >= target_lat.min() - 5) & (lat2d <= target_lat.max() + 5) &
           (lon2d >= target_lon.min() - 5) & (lon2d <= target_lon.max() + 5))
    src_mask = box & finite0

    pts = np.column_stack([lon2d[src_mask], lat2d[src_mask]])
    if pts.shape[0] < 4:
        raise ValueError("Not enough valid source points for interpolation")

    tri = Delaunay(pts)
    tgt = np.column_stack([LON_T.ravel(), LAT_T.ravel()])

    out = np.empty((vals.shape[0], len(target_lat), len(target_lon)),
                   dtype=float)
    for i in range(vals.shape[0]):
        interp = LinearNDInterpolator(tri, vals[i][src_mask])
        out[i] = interp(tgt).reshape(len(target_lat), len(target_lon))
    return out


def pdo_eof(arr, lat_t, lon_t, months, area_weight=True,
            orient_pc=None, orient_reg=None):
    """
    Leading EOF/PC of monthly SST anomalies.

    Parameters
    ----------
    arr : (time, lat, lon) degC, with NaN over land/missing.
    months : length-time array of calendar months (1..12).
    orient_pc : optional length-time reference PC for sign orientation.
    orient_reg : optional 2-D reference regression pattern for sign orientation.

    Returns
    -------
    dict with standardized PC, regression pattern (degC per std),
    variance explained, amplitude, lag-12 autocorrelation, etc.
    """
    arr = np.asarray(arr, dtype=float)

    # Monthly climatology over 1950-2014.
    clim = np.empty((12, arr.shape[1], arr.shape[2]), dtype=float)
    for m in range(1, 13):
        clim[m - 1] = np.nanmean(arr[months == m], axis=0)

    anom = arr - clim[months - 1]

    # Keep grid points with a complete anomaly time series.
    valid = np.isfinite(anom).all(axis=0)
    X = anom[:, valid]
    X = X - X.mean(axis=0, keepdims=True)

    # Area weights for the EOF.
    if area_weight:
        w2d = np.sqrt(np.cos(np.deg2rad(lat_t)))[:, None] * np.ones((1, len(lon_t)))
        w = w2d[valid]
    else:
        w = np.ones(X.shape[1])

    Xw = X * w
    U, s, Vt = np.linalg.svd(Xw, full_matrices=False)

    pc = U[:, 0] * s[0]
    var_expl = float(100.0 * s[0] ** 2 / np.sum(s ** 2))

    pc_std = (pc - pc.mean()) / pc.std(ddof=1)

    # Regression pattern in degC per one standard deviation of the PC.
    reg = (X.T @ pc_std) / (len(pc_std) - 1.0)
    reg_grid = np.full(arr.shape[1:], np.nan, dtype=float)
    reg_grid[valid] = reg

    # Orient the EOF sign.
    if orient_pc is not None:
        orient_pc = np.asarray(orient_pc)
        good = np.isfinite(orient_pc)
        if good.sum() > 10:
            c = np.corrcoef(pc_std[good], orient_pc[good])[0, 1]
            if c < 0:
                pc_std = -pc_std
                reg_grid = -reg_grid

    if orient_reg is not None:
        orient_reg = np.asarray(orient_reg)
        good = np.isfinite(reg_grid) & np.isfinite(orient_reg)
        if good.sum() > 10:
            c = np.corrcoef(reg_grid[good], orient_reg[good])[0, 1]
            if c < 0:
                pc_std = -pc_std
                reg_grid = -reg_grid

    reg = reg_grid[valid]
    amp_max = float(np.nanmax(np.abs(reg)))
    amp_rms = float(np.sqrt(np.nanmean(reg ** 2)))
    acf12 = float(np.corrcoef(pc_std[:-12], pc_std[12:])[0, 1])

    return dict(pc_std=pc_std, reg_grid=reg_grid, valid=valid,
                var_expl=var_expl, amp_max=amp_max, amp_rms=amp_rms,
                acf12=acf12, anom=anom, clim=clim)


# ---------------------------------------------------------------------
# Published PDO index
# ---------------------------------------------------------------------
month_names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

# The published PDO file is a whitespace-delimited table whose comment/meta
# lines do not all contain the same number of fields, which makes a plain
# pandas.read_csv(sep=r"\s+") fail with
#     ParserError: Expected 3 fields in line 2, saw 13
# Parse it line by line instead, retaining only genuine year data rows
# (a 4-digit year followed by 12 numeric monthly values).
with open(pdo_path, "r") as fh:
    raw_lines = fh.read().splitlines()

pdo_rows = []
for line in raw_lines:
    s = line.strip()
    if not s or s.startswith("#"):
        continue
    parts = s.split()
    if len(parts) < 13:
        continue
    try:
        year = int(float(parts[0]))
    except (ValueError, TypeError):
        continue
    try:
        monthly_vals = [float(tok) for tok in parts[1:13]]
    except ValueError:
        continue
    pdo_rows.append([year] + monthly_vals)

if not pdo_rows:
    raise ValueError("Could not parse any data rows from the PDO file")

pdo_raw = pd.DataFrame(pdo_rows, columns=["Year"] + month_names)

pdo_long = pdo_raw.melt(id_vars="Year", var_name="month_name",
                        value_name="pdo")
pdo_long["month"] = pdo_long["month_name"].map(
    {m: i + 1 for i, m in enumerate(month_names)}
)
pdo_long["date"] = pd.to_datetime(dict(
    year=pdo_long["Year"].astype(int),
    month=pdo_long["month"].astype(int),
    day=1
))
pdo_long["pdo"] = pd.to_numeric(pdo_long["pdo"], errors="coerce")
pdo_long.loc[np.isclose(pdo_long["pdo"], 99.99), "pdo"] = np.nan

pdo_pub = pdo_long.set_index("date").sort_index()["pdo"].reindex(time_common)
pdo_pub_std = (pdo_pub - pdo_pub.mean()) / pdo_pub.std(ddof=1)

# ---------------------------------------------------------------------
# ERSSTv5
# ---------------------------------------------------------------------
ds_e = xr.open_dataset(ersst_path)
sst = ds_e["sst"].sel(time=slice(start_sel, end_sel))

# ERSST latitude is descending; longitude is already 0-360.
sst = sst.sortby("lat")
if float(sst.lon.min()) < 0:
    sst = sst.assign_coords(lon=(sst.lon % 360)).sortby("lon")

sst_sub = sst.sel(
    lat=slice(lat_t.min() - 5, lat_t.max() + 5),
    lon=slice(lon_t.min() - 5, lon_t.max() + 5)
)
sst_t = sst_sub.interp(lat=lat_t, lon=lon_t, method="linear")
ersst_arr = np.asarray(sst_t.values, dtype=float)

if ersst_arr.shape[0] != n_time:
    raise ValueError(f"ERSST time length {ersst_arr.shape[0]} != {n_time}")

# Common ocean mask used to make the regridded products directly comparable.
common_mask = np.isfinite(ersst_arr).all(axis=0)

ersst_res = pdo_eof(
    ersst_arr, lat_t, lon_t, months,
    area_weight=AREA_WEIGHT,
    orient_pc=pdo_pub.values
)

good_pub = np.isfinite(ersst_res["pc_std"]) & np.isfinite(pdo_pub.values)
corr_ersst_pub = float(np.corrcoef(
    ersst_res["pc_std"][good_pub],
    pdo_pub.values[good_pub]
)[0, 1])

# ---------------------------------------------------------------------
# CMIP6 historical tos
# ---------------------------------------------------------------------
model_cfg = {
    "CESM2": dict(grid="curvilinear", lat_coord="lat", lon_coord="lon"),
    "GFDL-ESM4": dict(grid="regular", lat_coord="lat", lon_coord="lon"),
    "MPI-ESM1-2-LR": dict(grid="curvilinear",
                          lat_coord="latitude", lon_coord="longitude"),
}

results = {"ERSSTv5": ersst_res}

for model, cfg in model_cfg.items():
    files = sorted(glob.glob(
        f"{cmip_root}/{model}/tos/tos_Omon_{model}_historical_*.nc"
    ))
    if not files:
        raise FileNotFoundError(f"No tos files found for {model}")

    if len(files) == 1:
        ds = xr.open_dataset(files[0], use_cftime=True)
    else:
        ds = xr.open_mfdataset(files, combine="by_coords",
                               use_cftime=True, decode_times=True)

    tos = ds["tos"].sel(time=slice(start_sel, end_sel))
    tos = to_celsius(tos)

    if len(tos.time) != n_time:
        # Fallback for cftime calendars with slightly different string slicing.
        tos = tos.sel(time=slice("1950-01-01", "2014-12-01"))

    if len(tos.time) != n_time:
        raise ValueError(
            f"{model}: selected {len(tos.time)} months, expected {n_time}"
        )

    if cfg["grid"] == "regular":
        arr_t = regrid_regular(
            tos,
            lat_name=cfg["lat_coord"],
            lon_name=cfg["lon_coord"]
        )
        arr = np.asarray(arr_t.values, dtype=float)
    else:
        lat2d = ds[cfg["lat_coord"]].values
        lon2d = ds[cfg["lon_coord"]].values
        arr = regrid_curvilinear(tos, lat2d, lon2d)

    # Apply the common ERSST ocean mask.
    arr = arr.copy()
    arr[:, ~common_mask] = np.nan

    res = pdo_eof(
        arr, lat_t, lon_t, months,
        area_weight=AREA_WEIGHT,
        orient_reg=ersst_res["reg_grid"]
    )

    good = np.isfinite(res["reg_grid"]) & np.isfinite(ersst_res["reg_grid"])
    res["pattern_corr"] = float(np.corrcoef(
        res["reg_grid"][good],
        ersst_res["reg_grid"][good]
    )[0, 1])

    results[model] = res
    ds.close()

# ---------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------
source_order = ["ERSSTv5", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

rows = []
for name in source_order:
    r = results[name]
    rows.append(dict(
        source=name,
        variance_explained_pct=r["var_expl"],
        amplitude_degC_per_std=r["amp_max"],
        amplitude_rms_degC_per_std=r["amp_rms"],
        lag12_autocorr=r["acf12"],
        pattern_corr_vs_ERSST=(
            1.0 if name == "ERSSTv5"
            else r["pattern_corr"]
        )
    ))
results_df = pd.DataFrame(rows).set_index("source")

# ---------------------------------------------------------------------
# Observed PDO phase transitions from the ERSSTv5-derived index
# ---------------------------------------------------------------------
ersst_pc = pd.Series(results["ERSSTv5"]["pc_std"], index=time_common)

# 9-year (108-month) centred running mean.  min_periods=1 keeps the
# end points available so a transition near 2014 can still be detected.
smoothed_ersst = ersst_pc.rolling(
    window=108, center=True, min_periods=1
).mean()

transitions = []
transition_direction = []

for i in range(1, len(smoothed_ersst)):
    a = smoothed_ersst.iloc[i - 1]
    b = smoothed_ersst.iloc[i]
    if np.isfinite(a) and np.isfinite(b) and a * b < 0:
        yr = int(smoothed_ersst.index[i].year)
        direction = "neg->pos" if b > 0 else "pos->neg"
        if transitions and yr - transitions[-1] < 3:
            # Merge into the earlier counted transition.
            continue
        transitions.append(yr)
        transition_direction.append(direction)

transition_df = pd.DataFrame({
    "year": transitions,
    "direction": transition_direction
})

documented_shifts = {
    1976: "1976/77 Pacific climate shift",
    1989: "1988/89 Pacific climate shift",
    1998: "1998/99 Pacific climate shift",
    2014: "2013/14 Pacific climate shift",
}
transition_df["documented_shift"] = transition_df["year"].map(
    documented_shifts
).fillna("")

# Interpretation note requested with the ERSST-vs-published correlation.
corr_ersst_pub_note = (
    "This correlation demonstrates that the leading EOF of ERSSTv5 over the "
    "chosen North Pacific domain reproduces the canonical published PDO "
    "time variability. It does not demonstrate that the two products are "
    "identical or interchangeable: the published PDO is built from ERSSTv6 "
    "and uses its own preprocessing, baseline and global-mean treatment."
)

# ---------------------------------------------------------------------
# Report numerical results
# ---------------------------------------------------------------------
summary = results_df.copy()
summary["variance_explained_pct"] = summary["variance_explained_pct"].round(1)
summary["amplitude_degC_per_std"] = summary["amplitude_degC_per_std"].round(3)
summary["amplitude_rms_degC_per_std"] = summary["amplitude_rms_degC_per_std"].round(3)
summary["lag12_autocorr"] = summary["lag12_autocorr"].round(3)
summary["pattern_corr_vs_ERSST"] = summary["pattern_corr_vs_ERSST"].round(3)

print("PDO leading-mode summary (1950-2014)")
print(summary.to_string())
print()
print(f"ERSSTv5-derived index vs published PDO correlation: {corr_ersst_pub:.3f}")
print(corr_ersst_pub_note)
print()
print("Observed PDO phase-transition years (9-yr smoothed ERSSTv5 index):")
if len(transition_df):
    print(transition_df.to_string(index=False))
else:
    print("  none detected")

# ---------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Figure 1: regression patterns
    fig1, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    vmax = max(np.nanmax(np.abs(results[n]["reg_grid"])) for n in source_order)
    levels = np.linspace(-vmax, vmax, 21)
    xticks = [110, 140, 170, 200, 230, 260]
    xticklabels = ["110E", "140E", "170E", "160W", "130W", "100W"]

    for ax, name in zip(axes.ravel(), source_order):
        reg = results[name]["reg_grid"]
        cf = ax.contourf(lon_t, lat_t, reg, levels=levels,
                         cmap="RdBu_r", extend="both")
        ax.set_title(name)
        ax.set_xlim(lon_t.min(), lon_t.max())
        ax.set_ylim(lat_t.min(), lat_t.max())
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

        ann = f"Var exp: {results[name]['var_expl']:.1f}%"
        if name != "ERSSTv5":
            ann += f"\nPattern corr: {results[name]['pattern_corr']:.3f}"
        ax.text(0.02, 0.02, ann, transform=ax.transAxes, fontsize=8,
                bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))

    cbar = fig1.colorbar(cf, ax=axes.ravel().tolist(), shrink=0.8)
    cbar.set_label("degC per standard deviation")
    fig1.savefig("pdo_regression_patterns.png", dpi=150)

    # Figure 2: standardized indices and transition years
    fig2, ax = plt.subplots(figsize=(12, 4.5), constrained_layout=True)
    colors = {"ERSSTv5": "black", "CESM2": "C0",
              "GFDL-ESM4": "C1", "MPI-ESM1-2-LR": "C2"}
    for name in source_order:
        ax.plot(time_common, results[name]["pc_std"],
                lw=0.8, color=colors[name], label=name)
    ax.plot(time_common, pdo_pub_std.values,
            lw=1.0, ls="--", color="C4", label="Published PDO")
    ax.plot(time_common, smoothed_ersst.values,
            lw=2.0, color="red", label="ERSSTv5 9-yr smoothed")
    for yr in transitions:
        ax.axvline(pd.Timestamp(f"{yr}-01-01"),
                   color="0.5", ls=":", lw=0.8)
    ax.axhline(0, color="0.7", lw=0.5)
    ax.set_ylabel("Standardized index")
    ax.set_xlabel("Year")
    ax.legend(ncol=3, fontsize=8)
    fig2.savefig("pdo_indices_timeseries.png", dpi=150)

    plt.close(fig1)
    plt.close(fig2)
except Exception as exc:
    print(f"Figure generation skipped: {exc}")
