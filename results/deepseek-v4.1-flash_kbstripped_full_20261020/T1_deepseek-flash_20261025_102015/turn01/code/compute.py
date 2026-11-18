# Load ERSSTv5 SST over the Nino 3.4 region
ersst_sst = load_regular_data(
    dataset="ersst",
    variable="sst",
    time_range=["1950-01", "2014-12"],
    lat_range=[-5.0, 5.0],
    lon_range=[-170.0, -120.0],
)
ersst_sst = ersst_sst.where(np.abs(ersst_sst) < 1e10)

# Load CMIP6 historical tos for the three models over the same region/window
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
model_sst = {}
for model in models:
    da = load_regular_data(
        dataset="cmip6_historical",
        variable="tos",
        model=model,
        time_range=["1950-01", "2014-12"],
        lat_range=[-5.0, 5.0],
        lon_range=[-170.0, -120.0],
    )
    da = da.where(np.abs(da) < 1e10)
    model_sst[model] = da

# Compute monthly Nino 3.4 indices against each source's own 1950-2014 monthly climatology
base_period = ["1950-01", "2014-12"]
indices = {}

indices["ERSSTv5"] = compute_enso_nino_index(
    data=ersst_sst,
    nino_region="nino3.4",
    anomaly_freq="month",
    climatology="self",
    base_period=base_period,
).compute()

for model in models:
    indices[model] = compute_enso_nino_index(
        data=model_sst[model],
        nino_region="nino3.4",
        anomaly_freq="month",
        climatology="self",
        base_period=base_period,
    ).compute()

# Overall population standard deviation (ddof=0) over all months
std_all = {}
for name, idx in indices.items():
    std_all[name] = float(
        compute_basic_stat(data=idx, op="std", dim="time").compute()
    )

# Standard deviation separately for each calendar month
month_labels = [
    "Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
]
monthly_std = {}
for name, idx in indices.items():
    vals = []
    for m in range(1, 13):
        sub = idx.where(idx.time.dt.month == m, drop=True)
        vals.append(
            float(compute_basic_stat(data=sub, op="std", dim="time").compute())
        )
    monthly_std[name] = np.array(vals)

# Peak and minimum months of the seasonal amplitude cycle
peak_month = {name: int(np.argmax(monthly_std[name])) + 1 for name in indices}
peak_value = {name: float(np.max(monthly_std[name])) for name in indices}
min_month = {name: int(np.argmin(monthly_std[name])) + 1 for name in indices}
min_value = {name: float(np.min(monthly_std[name])) for name in indices}
max_min_ratio = {
    name: float(peak_value[name] / min_value[name]) for name in indices
}

# Amplitude bias and ratio against ERSSTv5
ersst_std = std_all["ERSSTv5"]
bias = {m: float(std_all[m] - ersst_std) for m in models}
amp_ratio = {m: float(std_all[m] / ersst_std) for m in models}
ranking = sorted(models, key=lambda m: abs(bias[m]))

# Strict |anomaly| > 1.0 degC exceedances
exceed_count = {}
exceed_pct = {}
for name, idx in indices.items():
    vals = idx.values
    cnt = int(np.sum(np.abs(vals) > 1.0))
    n_months = int(idx.sizes["time"])
    exceed_count[name] = cnt
    exceed_pct[name] = float(100.0 * cnt / n_months)

# Seasonal phase-locking comparison with ERSSTv5
peak_match = {m: (peak_month[m] == peak_month["ERSSTv5"]) for m in models}
peak_strength_ratio_to_ersst = {
    m: float(peak_value[m] / peak_value["ERSSTv5"]) for m in models
}

# Print numerical results
print("Overall Nino 3.4 anomaly standard deviation (ddof=0):")
for name in ["ERSSTv5"] + models:
    print(f"  {name}: {std_all[name]:.2f} degC")

print("\nAmplitude bias and ratio (model - ERSSTv5, model / ERSSTv5):")
for m in models:
    print(f"  {m}: bias = {bias[m]:+.2f} degC, ratio = {amp_ratio[m]:.2f}")
print("Ranking by smallest absolute amplitude bias:")
for rank, m in enumerate(ranking, start=1):
    print(f"  {rank}. {m}: |bias| = {abs(bias[m]):.2f} degC")

print("\nMonthly Nino 3.4 anomaly standard deviation (degC):")
monthly_df = pd.DataFrame(monthly_std, index=month_labels)
print(monthly_df.round(2).to_string())

print("\nSeasonal amplitude peak/minimum and ratio:")
for name in ["ERSSTv5"] + models:
    print(
        f"  {name}: max {month_labels[peak_month[name] - 1]} "
        f"{peak_value[name]:.2f} degC; "
        f"min {month_labels[min_month[name] - 1]} "
        f"{min_value[name]:.2f} degC; "
        f"max/min ratio {max_min_ratio[name]:.2f}"
    )

print("\nPeak-month and peak-strength comparison with ERSSTv5:")
for m in models:
    print(
        f"  {m}: peak month matches ERSSTv5 = {peak_match[m]}, "
        f"peak strength relative to ERSSTv5 = {peak_strength_ratio_to_ersst[m]:.2f}"
    )

print("\nExceedances (strict |Nino 3.4 anomaly| > 1.0 degC):")
for name in ["ERSSTv5"] + models:
    print(
        f"  {name}: {exceed_count[name]} months "
        f"({exceed_pct[name]:.1f}% of 780 months)"
    )
