# Computation code only
# ENSO amplitude assessment over Niño 3.4 (5S–5N, 170W–120W) for 1950-01 to 2014-12

time_range = ["1950-01", "2014-12"]
lat_range = [-5, 5]
lon_range = [-170, -120]  # 170W to 120W in -180..180 convention
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
sources = ["ERSSTv5"] + models

# 1) Load monthly SST on each source's own native grid over the Niño 3.4 box
cmip_sst = {}
for model in models:
    cmip_sst[model] = load_regular_data(
        dataset="cmip6_historical",
        variable="tos",
        model=model,
        time_range=time_range,
        lat_range=lat_range,
        lon_range=lon_range,
        chunk_time=12,
    )

ersst_sst = load_regular_data(
    dataset="ersst",
    variable="sst",
    time_range=time_range,
    lat_range=lat_range,
    lon_range=lon_range,
    chunk_time=12,
)

# 2) Load model area fields for the same regional subset
cmip_areas = {}
for model in models:
    cmip_areas[model] = load_areas(
        dataset="cmip6_historical",
        model=model,
        variable="tos",
        lat_range=lat_range,
        lon_range=lon_range,
    )

# 3) Compute monthly Niño 3.4 regional mean SST on each source's own grid
nino34_total = {}

nino34_total["ERSSTv5"] = compute_basic_stat(
    ersst_sst,
    op="area_weighted_mean",
    dim=["lat", "lon"],
).compute()

for model in models:
    nino34_total[model] = compute_basic_stat(
        cmip_sst[model],
        op="area_weighted_mean",
        dim=["lat", "lon"],
    ).compute()

# 4) Convert to Niño 3.4 anomalies relative to each source's own 1950–2014 monthly climatology
nino34_anom = {}
nino34_clim_monthly = {}

for src in sources:
    nino34_clim_monthly[src] = temporal_aggregate(
        nino34_total[src],
        freq="climatology",
        op="mean",
    ).compute()

    nino34_anom[src] = compute_anomaly_vs_own_climatology(
        nino34_total[src],
        freq="month",
        base_period=time_range,
        method="subtract",
    ).compute()

# Materialize exact month counts
n_months = {}
for src in sources:
    n_months[src] = int(nino34_anom[src].sizes["time"])

# 5) Full-period population standard deviation, amplitude bias/ratio, ranking
std_full = {}
for src in sources:
    std_full[src] = float(
        compute_basic_stat(nino34_anom[src], op="std", dim="time").compute().item()
    )

amplitude_bias = {}
amplitude_ratio = {}
for model in models:
    amplitude_bias[model] = std_full[model] - std_full["ERSSTv5"]
    amplitude_ratio[model] = std_full[model] / std_full["ERSSTv5"]

ranking_by_abs_bias = sorted(models, key=lambda m: abs(amplitude_bias[m]))

# 6) Calendar-month standard deviation across years for each source
month_std = {}
month_std_values = {}
peak_month = {}
peak_value = {}
min_month = {}
min_value = {}
phase_lock_ratio = {}

month_names = {
    1: "Jan", 2: "Feb", 3: "Mar", 4: "Apr", 5: "May", 6: "Jun",
    7: "Jul", 8: "Aug", 9: "Sep", 10: "Oct", 11: "Nov", 12: "Dec"
}

for src in sources:
    monthly_groups = temporal_aggregate(
        nino34_anom[src],
        freq="climatology",
        op="std",
    ).compute()

    month_std[src] = monthly_groups
    month_std_values[src] = [float(monthly_groups.sel(month=m).item()) for m in range(1, 13)]

    peak_idx = int(np.nanargmax(month_std_values[src])) + 1
    min_idx = int(np.nanargmin(month_std_values[src])) + 1

    peak_month[src] = peak_idx
    peak_value[src] = float(monthly_groups.sel(month=peak_idx).item())
    min_month[src] = min_idx
    min_value[src] = float(monthly_groups.sel(month=min_idx).item())
    phase_lock_ratio[src] = peak_value[src] / min_value[src]

observed_peak_month = peak_month["ERSSTv5"]
observed_phase_lock_ratio = phase_lock_ratio["ERSSTv5"]

reproduces_peak_month = {model: (peak_month[model] == observed_peak_month) for model in models}
phase_lock_ratio_bias = {model: phase_lock_ratio[model] - observed_phase_lock_ratio for model in models}
phase_lock_ratio_abs_error = {model: abs(phase_lock_ratio_bias[model]) for model in models}

# 7) Count and percentage of months with |anomaly| > 1.0 °C
exceedance = {}
exceed_count = {}
exceed_percent = {}

for src in sources:
    exceedance[src] = compute_exceedance_fraction(
        nino34_anom[src],
        threshold=1.0,
        comparison="abs_greater",
        dim="time",
        weighting="none",
        reference=0.0,
    )
    exceed_count[src] = int(exceedance[src]["n_exceeding"])
    exceed_percent[src] = float(exceedance[src]["fraction_percent"])

# Tabular summaries
full_period_summary = pd.DataFrame({
    "source": sources,
    "std_degC": [std_full[s] for s in sources],
    "n_months": [n_months[s] for s in sources],
    "months_abs_gt_1C": [exceed_count[s] for s in sources],
    "percent_abs_gt_1C": [exceed_percent[s] for s in sources],
}).set_index("source")

model_bias_summary = pd.DataFrame({
    "model": models,
    "std_degC": [std_full[m] for m in models],
    "amplitude_bias_degC": [amplitude_bias[m] for m in models],
    "abs_amplitude_bias_degC": [abs(amplitude_bias[m]) for m in models],
    "amplitude_ratio": [amplitude_ratio[m] for m in models],
    "peak_month": [month_names[peak_month[m]] for m in models],
    "reproduces_observed_peak_month": [reproduces_peak_month[m] for m in models],
    "phase_lock_ratio": [phase_lock_ratio[m] for m in models],
    "phase_lock_ratio_bias": [phase_lock_ratio_bias[m] for m in models],
    "phase_lock_ratio_abs_error": [phase_lock_ratio_abs_error[m] for m in models],
}).set_index("model").loc[ranking_by_abs_bias]

monthly_std_table = pd.DataFrame(
    {src: month_std_values[src] for src in sources},
    index=pd.Index(range(1, 13), name="month")
)
monthly_std_table.index = [month_names[m] for m in range(1, 13)]

seasonal_phase_summary = pd.DataFrame({
    "source": sources,
    "peak_month_num": [peak_month[s] for s in sources],
    "peak_month": [month_names[peak_month[s]] for s in sources],
    "peak_std_degC": [peak_value[s] for s in sources],
    "min_month_num": [min_month[s] for s in sources],
    "min_month": [month_names[min_month[s]] for s in sources],
    "min_std_degC": [min_value[s] for s in sources],
    "max_to_min_ratio": [phase_lock_ratio[s] for s in sources],
}).set_index("source")

# Rounded versions for reporting
full_period_summary_rounded = full_period_summary.copy()
full_period_summary_rounded["std_degC"] = full_period_summary_rounded["std_degC"].round(2)
full_period_summary_rounded["percent_abs_gt_1C"] = full_period_summary_rounded["percent_abs_gt_1C"].round(1)

model_bias_summary_rounded = model_bias_summary.copy()
for col in ["std_degC", "amplitude_bias_degC", "abs_amplitude_bias_degC", "amplitude_ratio",
            "phase_lock_ratio", "phase_lock_ratio_bias", "phase_lock_ratio_abs_error"]:
    model_bias_summary_rounded[col] = model_bias_summary_rounded[col].round(2)

monthly_std_table_rounded = monthly_std_table.round(2)

seasonal_phase_summary_rounded = seasonal_phase_summary.copy()
for col in ["peak_std_degC", "min_std_degC", "max_to_min_ratio"]:
    seasonal_phase_summary_rounded[col] = seasonal_phase_summary_rounded[col].round(2)

# Explicit result variables
result_full_period_std = {src: round(std_full[src], 2) for src in sources}
result_model_amplitude_bias = {m: round(amplitude_bias[m], 2) for m in models}
result_model_amplitude_ratio = {m: round(amplitude_ratio[m], 2) for m in models}
result_ranking_by_abs_bias = ranking_by_abs_bias
result_monthly_std = {src: [round(v, 2) for v in month_std_values[src]] for src in sources}
result_peak_month = {src: month_names[peak_month[src]] for src in sources}
result_peak_std = {src: round(peak_value[src], 2) for src in sources}
result_min_month = {src: month_names[min_month[src]] for src in sources}
result_phase_lock_ratio = {src: round(phase_lock_ratio[src], 2) for src in sources}
result_reproduces_peak_month = reproduces_peak_month
result_exceed_count = exceed_count
result_exceed_percent = {src: round(exceed_percent[src], 1) for src in sources}
