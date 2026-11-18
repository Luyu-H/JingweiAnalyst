# Load each source on its native grid, subset to the Nino 3.4 region, and compute
# the monthly Nino 3.4 anomaly relative to each source's own 1950-2014 monthly climatology.
sources = ["ERSSTv5", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

time_range = ["1950-01", "2014-12"]
lat_range = [-5.0, 5.0]
lon_range = [-170.0, -120.0]

index_series = {}

ersst = load_regular_data(
    dataset="ersst",
    variable="sst",
    time_range=time_range,
    lat_range=lat_range,
    lon_range=lon_range,
)
index_series["ERSSTv5"] = compute_enso_nino_index(
    data=ersst,
    nino_region="nino3.4",
    anomaly_freq="month",
    climatology="self",
    base_period=time_range,
).compute()

for model in models:
    da = load_regular_data(
        dataset="cmip6_historical",
        variable="tos",
        model=model,
        time_range=time_range,
        lat_range=lat_range,
        lon_range=lon_range,
    )
    index_series[model] = compute_enso_nino_index(
        data=da,
        nino_region="nino3.4",
        anomaly_freq="month",
        climatology="self",
        base_period=time_range,
    ).compute()

# Number of months per source.
n_months = {src: int(index_series[src].sizes["time"]) for src in sources}

# Overall population standard deviation (ddof=0) over all months.
overall_std = {}
for src in sources:
    overall_std[src] = float(
        compute_basic_stat(
            data=index_series[src],
            op="std",
            dim="time",
        ).compute()
    )

# Population standard deviation for each calendar month.
monthly_std = {src: np.empty(12, dtype=float) for src in sources}
for src in sources:
    da = index_series[src]
    for m in range(1, 13):
        sub = da.where(da["time"].dt.month == m, drop=True)
        monthly_std[src][m - 1] = float(
            compute_basic_stat(
                data=sub,
                op="std",
                dim="time",
            ).compute()
        )

# Peak/minimum calendar months, values, and max-to-min ratio.
peak_month = {}
peak_value = {}
min_month = {}
min_value = {}
max_min_ratio = {}

for src in sources:
    vals = monthly_std[src]
    ip = int(np.nanargmax(vals))
    im = int(np.nanargmin(vals))
    peak_month[src] = ip + 1
    peak_value[src] = float(vals[ip])
    min_month[src] = im + 1
    min_value[src] = float(vals[im])
    max_min_ratio[src] = float(vals[ip] / vals[im])

# Count individual months with |Nino 3.4 anomaly| strictly greater than 1.0 degC.
exceed_count = {}
exceed_pct = {}
for src in sources:
    cnt = int((abs(index_series[src]) > 1.0).sum().compute())
    exceed_count[src] = cnt
    exceed_pct[src] = 100.0 * cnt / n_months[src]

# Model amplitude bias and ratio relative to ERSSTv5, using unrounded overall std values.
amp_bias = {m: float(overall_std[m] - overall_std["ERSSTv5"]) for m in models}
amp_ratio = {m: float(overall_std[m] / overall_std["ERSSTv5"]) for m in models}
rank_abs_bias = sorted(models, key=lambda m: abs(amp_bias[m]))

# Peak-month and phase-locking comparison with ERSSTv5.
obs_peak_month = peak_month["ERSSTv5"]
obs_max_min_ratio = max_min_ratio["ERSSTv5"]
model_peak_match = {m: (peak_month[m] == obs_peak_month) for m in models}
model_ratio_difference = {m: float(max_min_ratio[m] - obs_max_min_ratio) for m in models}

# Print a compact summary.
print("Nino 3.4 index lengths:")
for src in sources:
    print(f"  {src:12s}: {n_months[src]} months")

print("\nOverall Nino 3.4 anomaly standard deviation (ddof=0):")
for src in sources:
    print(f"  {src:12s}: {overall_std[src]:.2f} degC")

print("\nCalendar-month standard deviations (ddof=0), January through December:")
print("Source       " + " ".join(f"{m:>5d}" for m in range(1, 13)))
for src in sources:
    print(f"{src:12s} " + " ".join(f"{v:5.2f}" for v in monthly_std[src]))

print("\nPeak/minimum calendar-month standard deviations:")
for src in sources:
    print(
        f"  {src:12s}: peak month={peak_month[src]:2d} "
        f"({peak_value[src]:.2f} degC), min month={min_month[src]:2d} "
        f"({min_value[src]:.2f} degC), ratio={max_min_ratio[src]:.2f}"
    )

print("\nAmplitude bias and ratio relative to ERSSTv5:")
for m in models:
    print(
        f"  {m:12s}: bias={amp_bias[m]:+.2f} degC, "
        f"ratio={amp_ratio[m]:.2f}"
    )
print("Ranking by absolute amplitude bias (smallest to largest):")
print("  " + " < ".join(rank_abs_bias))

print("\nExceedance of |Nino 3.4 anomaly| > 1.0 degC:")
for src in sources:
    print(
        f"  {src:12s}: {exceed_count[src]:3d} months "
        f"({exceed_pct[src]:.1f}% of {n_months[src]})"
    )

print("\nPeak-month match with ERSSTv5:")
for m in models:
    print(
        f"  {m:12s}: {'yes' if model_peak_match[m] else 'no'} "
        f"(peak month={peak_month[m]}, ERSSTv5 peak month={obs_peak_month}; "
        f"max/min ratio difference={model_ratio_difference[m]:+.2f})"
    )
