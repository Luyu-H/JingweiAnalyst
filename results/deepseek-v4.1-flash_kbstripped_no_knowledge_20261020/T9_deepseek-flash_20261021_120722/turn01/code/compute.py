# ---------------------------------------------------------------------------
# PDO assessment: three CMIP6 models vs ERSSTv5, 1950-01 to 2014-12
# ---------------------------------------------------------------------------
time_range  = ["1950-01", "2014-12"]
base_period = ["1950-01", "2014-12"]
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
order  = ["ERSSTv5"] + models

# ---- 1. load the GLOBAL SST fields (compute_pdo_index needs global coverage
#         because it removes the global-mean anomaly itself) -----------------
raw = {}
raw["ERSSTv5"] = load_regular_data(dataset="ersst", variable="sst",
                                   time_range=time_range)
for m in models:
    raw[m] = load_regular_data(dataset="cmip6_historical", variable="tos",
                               model=m, time_range=time_range)

# ---- 2. bilinear regridding onto the common 1-degree grid -----------------
sst_1deg = {k: regrid_horizontal(v, target="1deg", method="bilinear")
            for k, v in raw.items()}

# ---- 3. monthly anomalies vs each source's own 1950-2014 climatology ------
sst_anom = {k: compute_anomaly_vs_own_climatology(v, freq="month",
                                                  base_period=base_period)
            for k, v in sst_1deg.items()}

# ---- 4. leading North Pacific mode (tool cuts 20N-70N, 110E-100W itself) --
pdo_index, eof_pattern, explained_variance, regression_pattern = {}, {}, {}, {}
for k in order:
    idx, eof, ev, reg = compute_pdo_index(sst_anom[k])
    pdo_index[k]          = idx.compute()
    eof_pattern[k]        = eof.compute()
    explained_variance[k] = float(np.asarray(ev).ravel()[0])
    regression_pattern[k] = reg.compute()

# ---- 5. amplitude of the mode (deg C per standard deviation of the index) --
amplitude     = {}
amplitude_rms = {}
for k in order:
    v = np.asarray(regression_pattern[k].values, dtype=float)
    amplitude[k]     = float(np.nanmax(np.abs(v)))          # peak amplitude
    amplitude_rms[k] = float(np.sqrt(np.nanmean(v ** 2)))   # domain RMS

# ---- 6. lag-1-year (12-month) autocorrelation of each standardized index --
lag12_autocorr = {}
lag12_pvalue   = {}
for k in order:
    r = analyze_correlation(pdo_index[k], pdo_index[k], dim="time", lag=12)
    lag12_autocorr[k] = float(np.asarray(r["correlation"]).ravel()[0])
    lag12_pvalue[k]   = float(np.asarray(r["pvalue"]).ravel()[0])

# ---- 7. published NOAA/NCEI PDO index and its agreement with ERSSTv5 ------
# pdo is a basin-wide pattern index defined on its own domain -> no region arg
pdo_published = load_regular_data(dataset="cpc_indices", product="pdo",
                                  time_range=time_range)
pdo_published = pdo_published.squeeze().compute()

_r = analyze_correlation(pdo_index["ERSSTv5"], pdo_published, dim="time")
corr_published        = float(np.asarray(_r["correlation"]).ravel()[0])
corr_published_pvalue = float(np.asarray(_r["pvalue"]).ravel()[0])
corr_published_n      = int(np.asarray(_r["n_samples"]).ravel()[0])

# ---- 8. spatial (area weighted) pattern correlation vs ERSSTv5 ------------
pattern_corr, pattern_corr_p = {}, {}
for m in models:
    r = analyze_correlation(regression_pattern["ERSSTv5"], regression_pattern[m],
                            dim=["lat", "lon"], method="pearson",
                            weighting="area")
    pattern_corr[m]   = float(np.asarray(r["correlation"]).ravel()[0])
    pattern_corr_p[m] = float(np.asarray(r["pvalue"]).ravel()[0])

# ---- 9. 9-year (108-month) centred running mean of the ERSSTv5 index ------
smoothed_ersst = temporal_aggregate(pdo_index["ERSSTv5"], freq="rolling",
                                    op="mean", window=108, center=True).compute()
sv = np.asarray(smoothed_ersst.values, dtype=float)
st = pd.to_datetime(np.asarray(smoothed_ersst["time"].values))
sg = np.sign(sv)

crossings = []
for i in range(1, len(sv)):
    if not np.isfinite(sv[i]) or not np.isfinite(sv[i - 1]):
        continue
    if sg[i] == 0 or sg[i - 1] == 0:
        continue
    if sg[i] != sg[i - 1]:
        crossings.append(int(st[i].year))

# merge crossings less than 3 years after the last counted transition
transitions = []
for y in crossings:
    if (not transitions) or (y - transitions[-1] >= 3):
        transitions.append(y)

# documented Pacific climate shifts over 1950-2014
documented_shifts = {"1976/77": (1976, 1977),
                     "1988/89": (1988, 1989),
                     "1998/99": (1998, 1999)}
shift_match = {}
for name, yrs in documented_shifts.items():
    shift_match[name] = [y for y in transitions
                         if min(abs(y - yrs[0]), abs(y - yrs[1])) <= 1]

# ---- 10. summary tables ---------------------------------------------------
summary = pd.DataFrame({
    "explained_variance_pct": {k: round(explained_variance[k], 1) for k in order},
    "amplitude_degC_per_SD":  {k: round(amplitude[k], 3)         for k in order},
    "amplitude_RMS_degC_SD":  {k: round(amplitude_rms[k], 3)     for k in order},
    "lag12_autocorr":         {k: round(lag12_autocorr[k], 3)    for k in order},
})
summary.loc[[m for m in models], "pattern_corr_vs_ERSSTv5"] = \
    [round(pattern_corr[m], 3) for m in models]

best_model = max(models, key=lambda m: pattern_corr[m])

print(summary.to_string())
print()
print(f"ERSSTv5 index vs published NOAA PDO : r = {corr_published:.3f} "
      f"(p = {corr_published_pvalue:.3g}, n = {corr_published_n})")
print(f"PDO phase-transition years (108-month running mean) : {transitions}")
print(f"Matched documented shifts : {shift_match}")
print(f"Best pattern reproduction : {best_model} "
      f"(r = {pattern_corr[best_model]:.3f})")
