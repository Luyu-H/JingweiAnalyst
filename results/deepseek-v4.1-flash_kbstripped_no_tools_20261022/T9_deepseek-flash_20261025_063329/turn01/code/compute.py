# ============================================================
# 1. Common grid, time axis, and helpers
# ============================================================
import os
import glob
from scipy.interpolate import LinearNDInterpolator

ersst_path = "/data1/lorraine-20419/ocean_data/ERSSTv5/sst.mnmean.nc"
cmip6_root = "/data1/lorraine-20419/benchmark_cmip6/historical"
cpc_root = "/data1/lorraine-20419/ocean_data/CPC"

common_time = pd.date_range("1950-01-01", "2014-12-01", freq="MS")
target_lat = np.arange(20.0, 70.0 + 1.0, 1.0)
target_lon = np.arange(110.0, 260.0 + 1.0, 1.0)

def select_years(da, y0=1950, y1=2014):
    da = da.sortby("time")
    m = (da.time.dt.year >= y0) & (da.time.dt.year <= y1)
    idx = np.where(np.asarray(m.values))[0]
    da = da.isel(time=idx)
    n_expected = 12 * (y1 - y0 + 1)
    if da.sizes["time"] != n_expected:
        raise ValueError(f"Expected {n_expected} months, got {da.sizes['time']}")
    da = da.assign_coords(time=common_time)
    return da

def sort_lon_0360(da, lonname):
    if float(da[lonname].max()) <= 180.0:
        da = da.assign_coords({lonname: (da[lonname] % 360.0)}).sortby(lonname)
    else:
        da = da.sortby(lonname)
    return da

# ============================================================
# 2. ERSSTv5: open, select period, remove monthly climatology
#    and global-mean SST anomaly, then regrid to 1 deg.
# ============================================================
ds_ersst = xr.open_dataset(ersst_path)
sst = ds_ersst["sst"]
if "lev" in sst.dims:
    sst = sst.isel(lev=0, drop=True)

sst = select_years(sst)
sst = sst.sortby("lat")
sst = sort_lon_0360(sst, "lon")

clim_e = sst.groupby("time.month").mean("time")
anom_e = sst.groupby("time.month") - clim_e

lat_e = anom_e.lat.values
lon_e = anom_e.lon.values
w_e = np.cos(np.deg2rad(lat_e))[:, None] * np.ones((len(lat_e), len(lon_e)))
w_e = xr.DataArray(w_e, dims=("lat", "lon"), coords={"lat": anom_e.lat, "lon": anom_e.lon})
w_e = w_e.where(np.isfinite(anom_e))
gm_e = (anom_e * w_e).sum(dim=("lat", "lon")) / w_e.sum(dim=("lat", "lon"))
anom_e_gm = (anom_e - gm_e).compute()

reg_ersst = anom_e_gm.interp(
    lat=target_lat,
    lon=target_lon,
    method="linear"
).compute()

# ============================================================
# 3. CMIP6: open tos and areacello, remove monthly climatology,
#    remove global-mean SST anomaly, and regrid to 1 deg.
# ============================================================
model_info = {
    "CESM2": {
        "y_dim": "nlat", "x_dim": "nlon",
        "lat_coord": "lat", "lon_coord": "lon",
    },
    "GFDL-ESM4": {
        "y_dim": "lat", "x_dim": "lon",
        "lat_coord": "lat", "lon_coord": "lon",
    },
    "MPI-ESM1-2-LR": {
        "y_dim": "j", "x_dim": "i",
        "lat_coord": "latitude", "lon_coord": "longitude",
    },
}

regridded = {"ERSSTv5": reg_ersst}

for model, info in model_info.items():
    tos_files = sorted(glob.glob(
        f"{cmip6_root}/{model}/tos/tos_Omon_{model}_historical_*.nc"
    ))
    if not tos_files:
        raise FileNotFoundError(f"No tos files found for {model}")

    # Parse the actual variant label from the filename; never assume r1i1p1f1.
    variant = os.path.basename(tos_files[0]).split("_")[4]

    ds_tos = xr.open_mfdataset(tos_files, combine="by_coords", use_cftime=True)
    tos = ds_tos["tos"]
    for d in ("lev", "depth", "olevel"):
        if d in tos.dims:
            tos = tos.isel({d: 0}, drop=True)

    tos = select_years(tos)
    clim = tos.groupby("time.month").mean("time")
    anom = tos.groupby("time.month") - clim

    area_files = sorted(glob.glob(
        f"{cmip6_root}/{model}/areacello/areacello_Ofx_{model}_historical_{variant}_*.nc"
    ))
    if not area_files:
        area_files = sorted(glob.glob(
            f"{cmip6_root}/{model}/areacello/areacello_Ofx_{model}_historical_*.nc"
        ))
    if not area_files:
        raise FileNotFoundError(f"No areacello files found for {model}")

    area = xr.open_dataset(area_files[0])["areacello"]
    if "time" in area.dims:
        area = area.isel(time=0, drop=True)

    if model == "CESM2":
        sftof_files = sorted(glob.glob(
            f"{cmip6_root}/CESM2/sftof/sftof_Ofx_CESM2_historical_{variant}_*.nc"
        ))
        if not sftof_files:
            sftof_files = sorted(glob.glob(
                f"{cmip6_root}/CESM2/sftof/sftof_Ofx_CESM2_historical_*.nc"
            ))
        if not sftof_files:
            raise FileNotFoundError("No CESM2 sftof file found")
        sftof = xr.open_dataset(sftof_files[0])["sftof"]
        if "time" in sftof.dims:
            sftof = sftof.isel(time=0, drop=True)
        area = area * (sftof > 0)
    else:
        area = area.where(area > 0)

    w = area.expand_dims(time=anom.time).where(np.isfinite(anom))
    gm = (anom * w).sum(dim=(info["y_dim"], info["x_dim"])) / w.sum(
        dim=(info["y_dim"], info["x_dim"])
    )
    anom_gm = anom - gm

    if model == "GFDL-ESM4":
        # Already gr: regular lat/lon remapping, still interpolate to common 1 deg.
        anom_gm = anom_gm.sortby(info["lat_coord"]).sortby(info["lon_coord"])
        anom_gm = sort_lon_0360(anom_gm, info["lon_coord"])

        latv = anom_gm[info["lat_coord"]].values
        lonv = anom_gm[info["lon_coord"]].values
        sub = anom_gm.sel({
            info["lat_coord"]: slice(
                max(19.0, float(latv.min())),
                min(71.0, float(latv.max()))
            ),
            info["lon_coord"]: slice(
                max(109.0, float(lonv.min())),
                min(261.0, float(lonv.max()))
            ),
        })

        reg = sub.interp({
            info["lat_coord"]: target_lat,
            info["lon_coord"]: target_lon,
        }, method="linear").compute()

        if info["lat_coord"] != "lat":
            reg = reg.rename({info["lat_coord"]: "lat"})
        if info["lon_coord"] != "lon":
            reg = reg.rename({info["lon_coord"]: "lon"})

    else:
        # CESM2 and MPI-ESM1-2-LR are curvilinear.  Interpolate from the
        # native 2-D lat/lon coordinates to the common 1 deg grid.
        anom_gm = anom_gm.transpose("time", info["y_dim"], info["x_dim"]).compute()

        lat2 = anom_gm[info["lat_coord"]].values
        lon2 = anom_gm[info["lon_coord"]].values % 360.0
        data = anom_gm.values

        sel = (
            (lat2 >= 15.0) & (lat2 <= 75.0) &
            (lon2 >= 105.0) & (lon2 <= 265.0)
        )

        pts = np.column_stack([lon2[sel], lat2[sel]])
        vals = data[:, sel].T  # (npoints, ntime)

        good = np.isfinite(vals).all(axis=1)
        pts = pts[good]
        vals = vals[good]

        interp = LinearNDInterpolator(pts, vals)

        glon, glat = np.meshgrid(target_lon, target_lat)
        targ = np.column_stack([glon.ravel(), glat.ravel()])
        out = interp(targ)  # (ntarget, ntime)
        out = out.T.reshape(
            len(common_time), len(target_lat), len(target_lon)
        )

        reg = xr.DataArray(
            out,
            dims=("time", "lat", "lon"),
            coords={
                "time": common_time,
                "lat": target_lat,
                "lon": target_lon,
            },
        )

    regridded[model] = reg

# ============================================================
# 4. Build one shared valid-data mask and apply it to all four
#    regridded anomaly fields.
# ============================================================
source_names = ["ERSSTv5", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

arrays = {}
for name in source_names:
    da = regridded[name].transpose("time", "lat", "lon")
    arrays[name] = da.values

valid_masks = [np.isfinite(arrays[name]).all(axis=0) for name in source_names]
shared_mask = np.logical_and.reduce(valid_masks)

# ============================================================
# 5. Leading EOF, regression pattern, amplitude, variance,
#    and lag-1-year autocorrelation for each source.
# ============================================================
lat_t = target_lat
lon_t = target_lon
lon2d, lat2d = np.meshgrid(lon_t, lat_t)
area_w = np.cos(np.deg2rad(lat_t))[:, None] * np.ones(
    (len(lat_t), len(lon_t)), dtype=float
)

metrics = {}
patterns = {}
indices = {}

for name in source_names:
    arr = arrays[name].copy()
    arr[:, ~shared_mask] = np.nan

    X = arr[:, shared_mask]              # (time, ncells)
    w_cells = area_w[shared_mask]
    sw = np.sqrt(w_cells)
    Xw = X * sw[None, :]                 # sqrt(cos lat) weighting before EOF

    U, S, Vh = np.linalg.svd(Xw, full_matrices=False)

    variance_explained = 100.0 * S[0] ** 2 / np.sum(S ** 2)
    pc = U[:, 0] * S[0]

    # Standardize the leading PC over 1950-2014.
    index = (pc - pc.mean()) / pc.std(ddof=0)

    # Regression pattern: anomalies regressed on the standardized index.
    pattern_cells = (X * index[:, None]).mean(axis=0)

    # Fix sign: positive index = negative SST anomalies in the central
    # and western North Pacific.
    sign_region = (
        (lon2d >= 140.0) & (lon2d <= 220.0) &
        (lat2d >= 30.0) & (lat2d <= 50.0) &
        shared_mask
    )
    sign_cells = sign_region[shared_mask]
    if np.any(sign_cells) and np.nanmean(pattern_cells[sign_cells]) > 0.0:
        index = -index
        pattern_cells = -pattern_cells

    # Area-weighted RMS of the regression pattern over the common domain.
    amplitude = np.sqrt(
        np.sum(w_cells * pattern_cells ** 2) / np.sum(w_cells)
    )

    # Lag-1-year (12-month) autocorrelation.
    lag1yr = np.corrcoef(index[12:], index[:-12])[0, 1]

    pattern = np.full(shared_mask.shape, np.nan)
    pattern[shared_mask] = pattern_cells

    metrics[name] = {
        "variance_explained_pct": float(variance_explained),
        "amplitude_C_per_SD": float(amplitude),
        "lag1yr_autocorr": float(lag1yr),
    }
    patterns[name] = pattern
    indices[name] = index

# ============================================================
# 6. Spatial pattern correlation of each model with ERSSTv5.
#    Plain unweighted Pearson correlation over shared valid cells.
# ============================================================
pattern_corr = {}
for model in ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]:
    x = patterns[model][shared_mask]
    y = patterns["ERSSTv5"][shared_mask]
    pattern_corr[model] = float(np.corrcoef(x, y)[0, 1])

# ============================================================
# 7. Read the published PDO index.
# ============================================================
pdo_path = os.path.join(cpc_root, "ersst.v6.pdo.dat")
with open(pdo_path, "r") as f:
    pdo_lines = f.readlines()

pdo_base_period = None
for line in pdo_lines:
    if "base" in line.lower():
        pdo_base_period = line.strip()
        break

header_idx = None
for i, line in enumerate(pdo_lines):
    if ("Year" in line) and ("Jan" in line):
        header_idx = i
        break

data_start = header_idx + 1 if header_idx is not None else 0

pdo_records = []
for line in pdo_lines[data_start:]:
    s = line.strip()
    if not s or s.startswith("#"):
        continue
    parts = s.split()
    if len(parts) < 13:
        continue
    try:
        yr = int(parts[0])
        vals = [float(x) for x in parts[1:13]]
    except ValueError:
        continue
    pdo_records.append((yr, vals))

pdo_rows = []
for yr, vals in pdo_records:
    for mo, val in enumerate(vals, start=1):
        pdo_rows.append((yr, mo, val))

pdo_df = pd.DataFrame(pdo_rows, columns=["year", "month", "pdo"])
pdo_df["time"] = pd.to_datetime(dict(
    year=pdo_df["year"], month=pdo_df["month"], day=1
))
pdo_series = pdo_df.set_index("time")["pdo"]
pdo_series = pdo_series.replace(99.99, np.nan)
pdo_series = pdo_series.mask(pdo_series.abs() > 90.0)
pdo_series = pdo_series.reindex(common_time)

ersst_idx = indices["ERSSTv5"]
pdo_vals = pdo_series.values
good_pdo = np.isfinite(ersst_idx) & np.isfinite(pdo_vals)
pdo_corr_ersst = float(
    np.corrcoef(ersst_idx[good_pdo], pdo_vals[good_pdo])[0, 1]
)

# ============================================================
# 8. 9-year (108-month) centred running mean of the ERSSTv5
#    standardized index, and PDO phase-transition years.
# ============================================================
ersst_idx_series = pd.Series(indices["ERSSTv5"], index=common_time)
ersst_smoothed = ersst_idx_series.rolling(
    window=108, center=True, min_periods=108
).mean()

sm = ersst_smoothed.dropna()
transitions_raw = []
prev_sign = None
for date, val in sm.items():
    s = 1 if val > 0 else (-1 if val < 0 else 0)
    if s == 0:
        continue
    if prev_sign is None:
        prev_sign = s
    elif s != prev_sign:
        transitions_raw.append(date.year)
        prev_sign = s

# Merge any crossing less than 3 years after the last counted transition
# into that earlier transition, keeping the earlier year.
pdo_transition_years = []
for y in transitions_raw:
    if not pdo_transition_years or (y - pdo_transition_years[-1]) >= 3:
        pdo_transition_years.append(y)

# Compare with commonly documented Pacific climate shift years.
documented_shift_years = [1957, 1976, 1988, 1998]
documented_shift_comparison = {}
for y in pdo_transition_years:
    documented_shift_comparison[y] = [
        d for d in documented_shift_years if abs(y - d) <= 1
    ]

# ============================================================
# 9. Summary table and plotting helper variables.
# ============================================================
best_model = max(pattern_corr, key=pattern_corr.get)

summary_df = pd.DataFrame([
    {
        "source": name,
        "variance_explained_pct": metrics[name]["variance_explained_pct"],
        "amplitude_C_per_SD": metrics[name]["amplitude_C_per_SD"],
        "lag1yr_autocorr": metrics[name]["lag1yr_autocorr"],
        "pattern_corr_with_ERSSTv5": (
            1.0 if name == "ERSSTv5" else pattern_corr[name]
        ),
    }
    for name in source_names
])

pdo_std_series = (pdo_series - pdo_series.mean()) / pdo_series.std()

print("PDO metrics, 1950-2014, common 1-degree North Pacific grid")
print(summary_df.to_string(index=False))
print()
print(f"ERSSTv5-derived index vs published PDO correlation: {pdo_corr_ersst:.3f}")
print(f"Published PDO base_period field: {pdo_base_period}")
print(f"Best model by spatial pattern correlation: {best_model}")
print(f"PDO phase-transition years: {pdo_transition_years}")
print("Documented shift comparison:")
for y, matches in documented_shift_comparison.items():
    print(f"  {y}: {matches}")
