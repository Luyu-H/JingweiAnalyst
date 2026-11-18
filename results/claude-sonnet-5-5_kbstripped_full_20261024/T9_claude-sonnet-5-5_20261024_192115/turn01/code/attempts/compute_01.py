# ---------------- PDO assessment: ERSSTv5 vs three CMIP6 models, 1950-2014 ----------------
# Baseline: each field's own 1950-2014 monthly climatology. All fields are bilinearly regridded
# to the 1-degree grid first. compute_pdo_index receives the GLOBAL field and removes the global-mean
# anomaly itself, then takes the leading EOF over 20N-70N, 110E-100W.

t0, t1 = "1950-01", "2014-12"
base = ["1950-01", "2014-12"]
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
names = ["ERSSTv5"] + models


def month_start(da):
    """Put the time axis on month-start timestamps so indices from different sources align."""
    tt = pd.to_datetime(da["time"].values).to_period("M").to_timestamp()
    return da.assign_coords(time=tt)


def in_pdo_domain(da):
    lon360 = da["lon"] % 360
    return (da["lat"] >= 20) & (da["lat"] <= 70) & (lon360 >= 110) & (lon360 <= 260)


raw = {}
raw["ERSSTv5"] = load_regular_data(dataset="ersst", variable="sst", time_range=[t0, t1])
for m in models:
    raw[m] = load_regular_data(dataset="cmip6_historical", variable="tos", model=m,
                               time_range=[t0, t1])

pdo_index = {}
eof_pattern = {}
expl_var = {}
reg_pattern = {}
sst_1deg_anom = {}
global_mean_anom = {}
amplitude = {}
ac12 = {}

for name in names:
    f = raw[name]
    # mask unmasked fill values (~1e20)
    f = f.where(np.abs(f) < 1e10)
    units = str(f.attrs.get("units", "")).lower()
    if units in ("k", "kelvin"):
        f = transform_unit(f, from_unit="K", to_unit="degC")
    # global field onto the common 1-degree grid (bilinear, as prescribed)
    f1 = regrid_horizontal(f, target="1deg", method="bilinear").compute()
    # anomaly relative to own 1950-2014 monthly climatology
    an = compute_anomaly_vs_own_climatology(f1, freq="month", base_period=base)
    sst_1deg_anom[name] = an
    # area-weighted global-mean SST anomaly on the 1-degree field (reported; the tool removes it internally)
    global_mean_anom[name] = compute_basic_stat(an, op="area_weighted_mean", dim=["lat", "lon"]).compute()

    # PDO: leading EOF on the global field (box cut and global-mean removal inside the tool)
    idx, eof, ev, reg = compute_pdo_index(an, standardize=True)
    idx = month_start(idx.compute())
    pdo_index[name] = idx
    eof_pattern[name] = eof.compute()
    expl_var[name] = float(ev)
    reg = reg.compute()
    reg_pattern[name] = reg.where(in_pdo_domain(reg))

    # amplitude: area-weighted RMS of the regression pattern (degC per standard deviation)
    msq = compute_basic_stat(reg_pattern[name] ** 2, op="area_weighted_mean", dim=["lat", "lon"])
    amplitude[name] = float(np.sqrt(float(msq.compute())))

    # lag-12-month autocorrelation of the standardised monthly index
    ac = analyze_correlation(idx, idx, dim="time", lag=12)
    ac12[name] = float(ac["correlation"])

# ---------------- ERSSTv5 index vs published PDO (as is, no detrending) ----------------
pub = load_regular_data(dataset="cpc_indices", product="pdo", time_range=[t0, t1]).compute()
pub = month_start(pub)
pub = pub.where(pub > -90)  # drop any missing-value flags
obs_idx, pub_a = xr.align(pdo_index["ERSSTv5"], pub, join="inner")
cc = analyze_correlation(obs_idx, pub_a, dim="time", detrend=False)
corr_pub = float(cc["correlation"])
p_pub_corrected = float(cc["pvalue"])
p_pub_raw = float(cc["pvalue_unadjusted"])
n_pub = int(cc["n_samples"])
n_eff_pub = float(cc["n_effective"])
pub_std = (pub_a - pub_a.mean("time")) / pub_a.std("time")

# ---------------- spatial pattern correlation with ERSSTv5 (area weighted, no shared mask) ----------------
pattern_corr = {}
pattern_ncells = {}
for m in models:
    sc = analyze_correlation(reg_pattern[m], reg_pattern["ERSSTv5"], dim=["lat", "lon"],
                             weighting="area")
    pattern_corr[m] = float(sc["correlation"])
    pattern_ncells[m] = int(sc["n_samples"])
best_model = max(pattern_corr, key=pattern_corr.get)

# ---------------- PDO phase transitions of the ERSSTv5 index ----------------
smooth = temporal_aggregate(pdo_index["ERSSTv5"], freq="rolling", window=108, center=True).compute()
up = detect_threshold_events(smooth, 0.0, direction="above", min_duration=1)
dn = detect_threshold_events(smooth, 0.0, direction="below", min_duration=1)
runs = []
for _, r in up.iterrows():
    runs.append((pd.Timestamp(r["start"]), pd.Timestamp(r["end"]), "positive"))
for _, r in dn.iterrows():
    runs.append((pd.Timestamp(r["start"]), pd.Timestamp(r["end"]), "negative"))
runs.sort(key=lambda x: x[0])
# first run starts at the edge of the valid record, so it is not a crossing
crossings = [(r[0], r[2]) for r in runs[1:]]

counted = []
for t, sign in crossings:
    if counted and (t - counted[-1][0]).days / 365.25 < 3.0:
        continue  # merged into the earlier transition
    counted.append((t, sign))
transition_times = [c[0] for c in counted]
transition_years = [t.year for t in transition_times]
transition_signs = [c[1] for c in counted]

documented = []
for y in transition_years:
    if abs(y - 1977) <= 1:
        documented.append("1976/77 shift")
    elif abs(y - 1998) <= 2:
        documented.append("late-1990s (1998/99) shift")
    else:
        documented.append("no documented match")
transition_table = pd.DataFrame({"year": transition_years, "new_phase": transition_signs,
                                 "documented_shift": documented})

# ---------------- model vs observation consistency ----------------
consistency = {}
for m in models:
    consistency[m] = {
        "amplitude": amplitude[m],
        "amplitude_ratio_to_obs": amplitude[m] / amplitude["ERSSTv5"],
        "ac12": ac12[m],
        "ac12_diff_from_obs": ac12[m] - ac12["ERSSTv5"],
        "amplitude_consistent(+-25%)": abs(amplitude[m] / amplitude["ERSSTv5"] - 1) <= 0.25,
        "persistence_consistent(+-0.1)": abs(ac12[m] - ac12["ERSSTv5"]) <= 0.1,
    }

summary = pd.DataFrame({
    "explained_variance_%": {n: round(expl_var[n], 1) for n in names},
    "amplitude_degC_per_sd": {n: round(amplitude[n], 3) for n in names},
    "lag12_autocorr": {n: round(ac12[n], 3) for n in names},
    "pattern_corr_vs_ERSSTv5": {n: (round(pattern_corr[n], 3) if n in pattern_corr else np.nan) for n in names},
    "n_cells": {n: (pattern_ncells[n] if n in pattern_ncells else np.nan) for n in names},
})
print(summary)
print("ERSSTv5 index vs published PDO: r = %.3f, corrected p = %.3g (raw p = %.3g), n = %d, n_eff = %.1f"
      % (corr_pub, p_pub_corrected, p_pub_raw, n_pub, n_eff_pub))
print(transition_table)
print("Best pattern:", best_model)
