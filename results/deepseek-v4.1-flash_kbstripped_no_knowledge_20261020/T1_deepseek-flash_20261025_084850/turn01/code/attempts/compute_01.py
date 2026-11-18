# ---------------------------------------------------------------------------
# 1. Load the four SST sources over 1950-01 .. 2014-12, trimmed to a box
#    around the Nino 3.4 region (5S-5N, 170W-120W) on each native grid.
# ---------------------------------------------------------------------------
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
sources = ["ERSSTv5"] + models

time_range = ["1950-01", "2014-12"]
lat_box = [-10.0, 10.0]      # a little wider than 5S-5N
lon_box = [-180.0, -110.0]   # a little wider than 170W-120W

fields = {}
for m in models:
    fields[m] = load_regular_data(
        dataset="cmip6_historical", variable="tos", model=m,
        time_range=time_range, lat_range=lat_box, lon_range=lon_box,
    )
fields["ERSSTv5"] = load_regular_data(
    dataset="ersst", variable="sst",
    time_range=time_range, lat_range=lat_box, lon_range=lon_box,
)

# ---------------------------------------------------------------------------
# 2. Monthly Nino 3.4 index for each source, anomaly w.r.t. that source's own
#    1950-2014 monthly climatology, area weighted on its own grid.
# ---------------------------------------------------------------------------
nino = {}
for s in sources:
    nino[s] = compute_enso_nino_index(
        data=fields[s],
        nino_region="nino3.4",
        anomaly_freq="month",
        climatology="self",
        base_period=time_range,
    ).compute()

# ---------------------------------------------------------------------------
# 3. Full-period population standard deviation (ddof = 0, unweighted in time)
# ---------------------------------------------------------------------------
std_full = {}
for s in sources:
    std_full[s] = float(compute_basic_stat(data=nino[s], op="std", dim="time").compute())

# ---------------------------------------------------------------------------
# 4. Amplitude bias / ratio against ERSSTv5 and ranking
# ---------------------------------------------------------------------------
obs_std = std_full["ERSSTv5"]
ampl_bias = {m: std_full[m] - obs_std for m in models}
ampl_ratio = {m: std_full[m] / obs_std for m in models}
ampl_ranking = sorted(models, key=lambda m: abs(ampl_bias[m]))

# ---------------------------------------------------------------------------
# 5. Calendar-month standard deviations (12 groups of 65 values per source)
# ---------------------------------------------------------------------------
std_month = {}
for s in sources:
    vals = []
    for mm in range(1, 13):
        sub = nino[s].where(nino[s].time.dt.month == mm, drop=True)
        vals.append(float(compute_basic_stat(data=sub, op="std", dim="time").compute()))
    std_month[s] = np.array(vals)

# ---------------------------------------------------------------------------
# 6. Peak / trough month of the seasonal amplitude cycle
# ---------------------------------------------------------------------------
peak_month = {s: int(np.argmax(std_month[s])) + 1 for s in sources}
peak_value = {s: float(np.max(std_month[s])) for s in sources}
trough_month = {s: int(np.argmin(std_month[s])) + 1 for s in sources}
trough_value = {s: float(np.min(std_month[s])) for s in sources}
peak_trough_ratio = {s: peak_value[s] / trough_value[s] for s in sources}

# ---------------------------------------------------------------------------
# 7. Fraction / number of months with |anomaly| strictly > 1.0 degC
# ---------------------------------------------------------------------------
exc_count = {}
exc_percent = {}
for s in sources:
    r = compute_exceedance_fraction(
        data=nino[s], threshold=1.0, comparison="abs_greater",
        dim=None, weighting="none",
    )
    exc_count[s] = int(np.asarray(r["n_exceeding"]).item())
    exc_percent[s] = float(np.asarray(r["fraction_percent"]).item())
n_months = int(nino["ERSSTv5"].sizes["time"])

# ---------------------------------------------------------------------------
# 8. Report
# ---------------------------------------------------------------------------
month_names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

print(f"Months per series: {n_months}")
print("\n--- Full-period Nino 3.4 anomaly standard deviation (degC) ---")
for s in sources:
    print(f"  {s:14s} {std_full[s]:.2f}")

print("\n--- Amplitude bias and ratio vs ERSSTv5 ---")
for m in models:
    print(f"  {m:14s} bias = {ampl_bias[m]:+.2f} degC   ratio = {ampl_ratio[m]:.2f}")
print("  Ranking (smallest -> largest |bias|): " + " < ".join(ampl_ranking))

print("\n--- Calendar-month standard deviation (degC) ---")
hdr = "Month  " + "".join(f"{s:>14s}" for s in sources)
print(hdr)
for i, mn in enumerate(month_names):
    print(f"{mn:5s}  " + "".join(f"{std_month[s][i]:14.2f}" for s in sources))

print("\n--- Seasonal amplitude extremes ---")
for s in sources:
    print(f"  {s:14s} peak {month_names[peak_month[s]-1]} ({peak_value[s]:.2f} degC), "
          f"trough {month_names[trough_month[s]-1]} ({trough_value[s]:.2f} degC), "
          f"ratio {peak_trough_ratio[s]:.2f}")

print("\n--- |anomaly| > 1.0 degC ---")
for s in sources:
    print(f"  {s:14s} {exc_count[s]:4d} months  ({exc_percent[s]:.1f} %)")

# ---------------------------------------------------------------------------
# 9. Bundles for plotting
# ---------------------------------------------------------------------------
series_time = {s: np.asarray(nino[s].time.values) for s in sources}
series_values = {s: np.asarray(nino[s].values, dtype=float) for s in sources}
months_axis = np.arange(1, 13)
