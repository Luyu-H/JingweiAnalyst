# Load ERSSTv5 SST over the Niño 3.4 box: 5S–5N, 170W–120W for 1950-01 to 2023-12
ersst_box = load_regular_data(
    dataset="ersst",
    variable="sst",
    time_range=["1950-01", "2023-12"],
    lat_range=[-5, 5],
    lon_range=[-170, -120],
    chunk_time=12,
)

# Area-weighted monthly regional mean SST series on the native ERSST grid
nino34_sst = compute_basic_stat(
    data=ersst_box,
    op="area_weighted_mean",
    dim=["lat", "lon"],
)

# Reconstruct Niño 3.4 anomalies from the same SST series with two baselines
RECON_9120 = compute_anomaly_vs_own_climatology(
    data=nino34_sst,
    freq="month",
    base_period=["1991-01", "2020-12"],
    method="subtract",
).assign_attrs(units="degC")

RECON_8110 = compute_anomaly_vs_own_climatology(
    data=nino34_sst,
    freq="month",
    base_period=["1981-01", "2010-12"],
    method="subtract",
).assign_attrs(units="degC")

# Load published CPC monthly Niño 3.4 anomaly series and official ONI
CPC_9120 = load_regular_data(
    dataset="cpc_indices",
    variable="NINO3.4",
    product="ersst5_91-20",
    field="anomaly",
    time_range=["1950-01", "2023-12"],
    chunk_time=12,
).assign_attrs(units="degC")

CPC_8110 = load_regular_data(
    dataset="cpc_indices",
    variable="NINO3.4",
    product="ersst5_81-10",
    field="anomaly",
    time_range=["1950-01", "2023-12"],
    chunk_time=12,
).assign_attrs(units="degC")

CPC_ONI = load_regular_data(
    dataset="cpc_indices",
    variable="NINO3.4",
    product="oni",
    field="anomaly",
    time_range=["1950-01", "2023-12"],
    chunk_time=12,
).assign_attrs(units="degC")

# Materialize the series
RECON_9120 = RECON_9120.compute()
RECON_8110 = RECON_8110.compute()
CPC_9120 = CPC_9120.compute()
CPC_8110 = CPC_8110.compute()
CPC_ONI = CPC_ONI.compute()

# Agreement metrics between reconstructed 1991-2020 baseline index and CPC published series
corr_out = analyze_correlation(
    var1=RECON_9120,
    var2=CPC_9120,
    dim="time",
    method="pearson",
    lag=0,
    detrend=False,
    min_periods=10,
)

corr_value = float(corr_out["correlation"].compute().item())
corr_pvalue = float(corr_out["pvalue"].compute().item())
corr_pvalue_unadjusted = float(corr_out["pvalue_unadjusted"].compute().item())
corr_n = int(corr_out["n_samples"].compute().item())
corr_neff = float(corr_out["n_effective"].compute().item())

diff_metrics = compute_difference_metrics(
    data=RECON_9120,
    reference=CPC_9120,
    dim=None,
    weighting="none",
)

rmse_value = float(diff_metrics["rmse"])
max_abs_dev_value = float(diff_metrics["max_abs_error"])
max_abs_dev_at = diff_metrics["max_abs_error_at"]

# Robust conversion of returned location metadata to a Timestamp
max_abs_dev_time = None
if max_abs_dev_at is not None:
    if isinstance(max_abs_dev_at, dict):
        time_like = max_abs_dev_at.get("time", max_abs_dev_at)
    else:
        time_like = max_abs_dev_at
    max_abs_dev_time = pd.Timestamp(time_like)

RECON_minus_CPC = diff_metrics["difference"].assign_attrs(units="degC").compute()

# ONI-style running mean from the reconstructed monthly series
RECON_9120_3mo = temporal_aggregate(
    data=RECON_9120,
    freq="rolling",
    op="mean",
    window=3,
    center=True,
).assign_attrs(units="degC").compute()

# Detect El Niño and La Niña events from reconstructed running mean
recon_warm_events = detect_threshold_events(
    data=RECON_9120_3mo,
    threshold=0.5,
    direction="above",
    min_duration=5,
)

recon_cool_events = detect_threshold_events(
    data=RECON_9120_3mo,
    threshold=-0.5,
    direction="below",
    min_duration=5,
)

recon_el_nino_count = int(len(recon_warm_events))
recon_top3_warm = recon_warm_events.reindex(
    recon_warm_events["peak_value"].sort_values(ascending=False).index
).head(3).copy()

# Apply identical event detection to published CPC ONI
oni_warm_events = detect_threshold_events(
    data=CPC_ONI,
    threshold=0.5,
    direction="above",
    min_duration=5,
)

oni_el_nino_count = int(len(oni_warm_events))
oni_top3_warm = oni_warm_events.reindex(
    oni_warm_events["peak_value"].sort_values(ascending=False).index
).head(3).copy()

# Baseline-period offsets over 1950-01 to 2020-12
offset_slice = slice("1950-01", "2020-12")

recon_offset = (RECON_8110 - RECON_9120).assign_attrs(units="degC").sel(time=offset_slice).compute()
published_offset = (CPC_8110 - CPC_9120).assign_attrs(units="degC").sel(time=offset_slice).compute()

recon_offset_mean = float(
    compute_basic_stat(recon_offset, op="mean", dim="time").compute().item()
)
published_offset_mean = float(
    compute_basic_stat(published_offset, op="mean", dim="time").compute().item()
)
offset_agreement = published_offset_mean - recon_offset_mean

# Prepare concise event summaries
recon_top3_summary = [
    {
        "start": pd.Timestamp(row["start"]),
        "end": pd.Timestamp(row["end"]),
        "duration": int(row["duration"]),
        "peak_time": pd.Timestamp(row["peak_time"]),
        "peak_value": float(row["peak_value"]),
        "mean_value": float(row["mean_value"]),
    }
    for _, row in recon_top3_warm.iterrows()
]

oni_top3_summary = [
    {
        "start": pd.Timestamp(row["start"]),
        "end": pd.Timestamp(row["end"]),
        "duration": int(row["duration"]),
        "peak_time": pd.Timestamp(row["peak_time"]),
        "peak_value": float(row["peak_value"]),
        "mean_value": float(row["mean_value"]),
    }
    for _, row in oni_top3_warm.iterrows()
]

# Print requested numerical outputs
print(f"RECON_9120 vs CPC ersst5_91-20, 1950-01 to 2023-12")
print(f"Pearson correlation = {corr_value:.6f}")
print(f"RMSE = {rmse_value:.4f} degC")
if max_abs_dev_time is not None:
    print(
        f"Largest absolute deviation = {max_abs_dev_value:.4f} degC at "
        f"{max_abs_dev_time.strftime('%Y-%m')}"
    )
else:
    print(f"Largest absolute deviation = {max_abs_dev_value:.4f} degC")

print("\nRECON_9120 ONI-style El Niño events, 1950-2023")
print(f"Count = {recon_el_nino_count}")
for i, ev in enumerate(recon_top3_summary, start=1):
    print(
        f"{i}. peak {ev['peak_time'].strftime('%Y-%m')} : "
        f"{ev['peak_value']:.2f} degC "
        f"(start {ev['start'].strftime('%Y-%m')}, end {ev['end'].strftime('%Y-%m')}, "
        f"duration {ev['duration']})"
    )

print("\nCPC ONI El Niño events, 1950-2023")
print(f"Count = {oni_el_nino_count}")
for i, ev in enumerate(oni_top3_summary, start=1):
    print(
        f"{i}. peak {ev['peak_time'].strftime('%Y-%m')} : "
        f"{ev['peak_value']:.2f} degC "
        f"(start {ev['start'].strftime('%Y-%m')}, end {ev['end'].strftime('%Y-%m')}, "
        f"duration {ev['duration']})"
    )

print("\nBaseline offsets over 1950-01 to 2020-12")
print(f"Mean reconstructed offset (RECON_8110 - RECON_9120) = {recon_offset_mean:.3f} degC")
print(f"Mean published offset (CPC_8110 - CPC_9120) = {published_offset_mean:.3f} degC")
print(f"Published - reconstructed mean offset difference = {offset_agreement:.3f} degC")

oni_count_explanation = (
    "Any count difference between the reconstructed ONI-style events and the published CPC ONI "
    "comes from definition/production details: the reconstruction applies a centered 3-month running "
    "mean to monthly Niño 3.4 anomalies rebuilt from ERSSTv5 area-averaged SST using a fixed "
    "1991-2020 monthly climatology, whereas the official ONI is the published operational series "
    "and may differ slightly because of CPC production conventions, exact averaging/rounding, and "
    "the fact that event classification is being recomputed here from the monthly ONI values with "
    "the same threshold rule rather than taken from CPC's separate historical episode listings."
)
print("\nExplanation:")
print(oni_count_explanation)
