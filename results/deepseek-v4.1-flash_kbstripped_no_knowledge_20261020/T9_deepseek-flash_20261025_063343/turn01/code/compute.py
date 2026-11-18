# =============================================================================
# PDO assessment: CMIP6 (CESM2, GFDL-ESM4, MPI-ESM1-2-LR) vs ERSSTv5, 1950-2014
# =============================================================================
tr = ["1950-01", "2014-12"]
bp = ["1950-01", "2014-12"]
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
source_names = ["ERSSTv5"] + models

# -----------------------------------------------------------------------------
# 1. Load the data (global fields: compute_pdo_index needs the whole planet,
#    because it removes the GLOBAL mean anomaly before taking the EOF)
# -----------------------------------------------------------------------------
ersst_raw = load_regular_data(dataset="ersst", variable="sst", time_range=tr)
cmip_raw = {m: load_regular_data(dataset="cmip6_historical", variable="tos",
                                 model=m, time_range=tr) for m in models}
pdo_pub = load_regular_data(dataset="cpc_indices", product="pdo", time_range=tr)

# -----------------------------------------------------------------------------
# 2. Regrid every source to the regular 1-degree grid (bilinear) and mask land
# -----------------------------------------------------------------------------
def to_degc(da):
    u = str(da.attrs.get("units", "degC")).strip().lower()
    if u in ("k", "kelvin", "degrees_k", "degk", "degreesk", "degrees_kelvin"):
        return transform_unit(da, from_unit="K", to_unit="degC")
    return da

fields = {"ERSSTv5": regrid_horizontal(to_degc(ersst_raw), target="1deg", method="bilinear")}
for m in models:
    fields[m] = regrid_horizontal(to_degc(cmip_raw[m]), target="1deg", method="bilinear")

ref = fields["ERSSTv5"]
ocean = (load_landsea_mask() > 1)
ocean = ocean.reindex(lat=ref["lat"], lon=ref["lon"], method="nearest").fillna(False).astype(bool)
fields = {k: v.where(ocean) for k, v in fields.items()}

# -----------------------------------------------------------------------------
# 3. Monthly anomalies relative to each source's own 1950-2014 monthly climatology
# -----------------------------------------------------------------------------
anoms = {k: compute_anomaly_vs_own_climatology(v, freq="month", base_period=bp)
         for k, v in fields.items()}

# -----------------------------------------------------------------------------
# 4. PDO = leading EOF of the 20N-70N, 110E-100W SST anomalies (tool does the
#    box selection and the global-mean removal internally)
# -----------------------------------------------------------------------------
def orient_mode(idx, reg, eof):
    """Fix the arbitrary EOF sign: positive loadings in the eastern North Pacific."""
    pat = reg.sortby("lat").sel(lat=slice(20, 45))
    pat = pat.assign_coords(lon=pat["lon"] % 360).sortby("lon").sel(lon=slice(200, 232))
    v = float(pat.mean(skipna=True).compute())
    if np.isfinite(v) and v < 0:
        return -idx, -reg, -eof
    return idx, reg, eof

idx_s, reg_s, ev_s = {}, {}, {}
for k in source_names:
    idx, eof, ve, reg = compute_pdo_index(anoms[k], standardize=True)
    idx, reg, eof = orient_mode(idx, reg, eof)
    idx_s[k] = idx
    reg_s[k] = reg
    ev_s[k] = float(ve.compute()) if hasattr(ve, "compute") else float(ve)

# mode amplitude = spatial standard deviation of the regression pattern
amp_s = {k: float(compute_basic_stat(reg_s[k], op="std", dim=["lat", "lon"]).compute())
         for k in source_names}

# lag-1-year (12-month) autocorrelation of the standardized monthly index
ac_s = {k: float(analyze_correlation(idx_s[k], idx_s[k], dim="time", lag=12)["correlation"].compute())
        for k in source_names}

# ERSSTv5-derived index vs the published PDO series
corr_pub = float(analyze_correlation(idx_s["ERSSTv5"], pdo_pub, dim="time")["correlation"].compute())

# spatial pattern correlation of each model's regression pattern with ERSSTv5
pcorr = {m: float(analyze_correlation(reg_s[m], reg_s["ERSSTv5"], dim=["lat", "lon"],
                                      method="pearson", weighting="area")["correlation"].compute())
         for m in models}

# -----------------------------------------------------------------------------
# 5. 9-year (108-month) centred running mean of the ERSSTv5 index and the
#    zero-crossing transition years, merging crossings < 3 years apart
# -----------------------------------------------------------------------------
sm_da = temporal_aggregate(idx_s["ERSSTv5"], freq="rolling", window=108, center=True).compute()
tvals = pd.to_datetime(sm_da["time"].values)
frac_year = tvals.year + (tvals.month - 0.5) / 12.0
v = np.asarray(sm_da.values, dtype=float)

sgn = np.sign(v)
for i in range(1, len(sgn)):
    if np.isfinite(sgn[i]) and sgn[i] == 0:
        sgn[i] = sgn[i - 1]

cross_years = []
for i in range(len(v) - 1):
    if np.isfinite(v[i]) and np.isfinite(v[i + 1]) and sgn[i] * sgn[i + 1] < 0:
        cross_years.append(int(np.floor(frac_year[i + 1])))

transitions = []
for y in cross_years:
    if not transitions or (y - transitions[-1]) >= 3:
        transitions.append(y)
transitions = sorted(set(int(y) for y in transitions))

documented_shifts = [1976, 1988, 1998, 2014]
matched_shifts = {y: min(documented_shifts, key=lambda d: abs(d - y)) for y in transitions}
shift_matches = {d: [y for y in transitions if min(documented_shifts, key=lambda dd: abs(dd - y)) == d]
                 for d in documented_shifts}

# concrete arrays for plotting
idx_plot = {k: idx_s[k].compute() for k in source_names}
pdo_pub_c = pdo_pub.compute()

best_model = max(pcorr, key=pcorr.get)

summary = pd.DataFrame({
    "variance_explained_pct": pd.Series(ev_s).round(1),
    "amplitude_degC_per_SD":  pd.Series(amp_s).round(3),
    "lag12_autocorr":         pd.Series(ac_s).round(3),
    "pattern_corr_vs_ERSST":  pd.Series(pcorr).reindex(source_names).round(3),
})

print(summary)
print("\nERSSTv5-derived index vs published PDO, Pearson r = %.3f" % corr_pub)
print("Best pattern reproduction:", best_model, "(r = %.3f)" % pcorr[best_model])
print("\nSmoothed ERSSTv5 PDO transitions (merging < 3 yr):", transitions)
print("Nearest documented shift for each transition:", matched_shifts)
