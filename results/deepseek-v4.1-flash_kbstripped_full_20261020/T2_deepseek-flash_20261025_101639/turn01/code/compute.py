# Load ERSSTv5 SST
sst = load_regular_data(
    dataset="ersst",
    variable="sst",
    field="total",
    time_range=["1950-01", "2023-12"],
    lat_range=[-5, 5],
    lon_range=[-170, -120],
)

# Reconstruct Nino 3.4 anomalies against two baselines
recon_9120_field = compute_anomaly_vs_own_climatology(
    sst, freq="month", base_period=["1991-01", "2020-12"], method="subtract"
)
recon_9120 = compute_basic_stat(
    recon_9120_field, op="area_weighted_mean", dim=["lat", "lon"]
).compute()

recon_8110_field = compute_anomaly_vs_own_climatology(
    sst, freq="month", base_period=["1981-01", "2010-12"], method="subtract"
)
recon_8110 = compute_basic_stat(
    recon_8110_field, op="area_weighted_mean", dim=["lat", "lon"]
).compute()

# Load CPC published series
cpc_9120 = load_regular_data(
    dataset="cpc_indices",
    variable="NINO3.4",
    product="ersst5_91-20",
    field="anomaly",
    time_range=["1950-01", "2023-12"],
).compute()

cpc_8110 = load_regular_data(
    dataset="cpc_indices",
    variable="NINO3.4",
    product="ersst5_81-10",
    field="anomaly",
    time_range=["1950-01", "2023-12"],
).compute()

oni = load_regular_data(
    dataset="cpc_indices",
    variable="NINO3.4",
    product="oni",
    field="anomaly",
    time_range=["1950-01", "2023-12"],
).compute()

# Put compared series on a common time axis
recon_9120, cpc_9120 = xr.align(recon_9120, cpc_9120, join="inner")
recon_8110, cpc_8110 = xr.align(recon_8110, cpc_8110, join="inner")

# 3-month centred running means
recon_9120_sm = temporal_aggregate(
    recon_9120, freq="rolling", window=3, center=True, op="mean"
).compute()
cpc_9120_sm = temporal_aggregate(
    cpc_9120, freq="rolling", window=3, center=True, op="mean"
).compute()

# Event detection: ONI criterion, >=5 consecutive 3-month mean values beyond threshold
recon_el_nino_events = detect_threshold_events(
    recon_9120_sm, threshold=0.5, direction="above", min_duration=5
)
recon_la_nina_events = detect_threshold_events(
    recon_9120_sm, threshold=-0.5, direction="below", min_duration=5
)
oni_el_nino_events = detect_threshold_events(
    oni, threshold=0.5, direction="above", min_duration=5
)
oni_la_nina_events = detect_threshold_events(
    oni, threshold=-0.5, direction="below", min_duration=5
)

n_recon_el_nino = len(recon_el_nino_events)
n_recon_la_nina = len(recon_la_nina_events)
n_oni_el_nino = len(oni_el_nino_events)
n_oni_la_nina = len(oni_la_nina_events)

def top_events(df, n=3):
    if df.empty:
        return pd.DataFrame(columns=["peak_year", "peak_month", "peak_value", "duration"])
    d = df.sort_values("peak_value", ascending=False).head(n).copy()
    times = pd.to_datetime(d["peak_time"])
    d["peak_year"] = times.dt.year.values
    d["peak_month"] = times.dt.month.values
    return d[["peak_year", "peak_month", "peak_value", "duration"]].reset_index(drop=True)

recon_top3 = top_events(recon_el_nino_events)
oni_top3 = top_events(oni_el_nino_events)

# Correlation
corr_dict = analyze_correlation(
    recon_9120, cpc_9120, dim="time", method="pearson", detrend=False
)
corr_value = float(corr_dict["correlation"].compute()) if hasattr(corr_dict["correlation"], "compute") else float(corr_dict["correlation"])
corr_pvalue = float(corr_dict["pvalue"].compute()) if hasattr(corr_dict["pvalue"], "compute") else float(corr_dict["pvalue"])
corr_n_eff = float(corr_dict["n_effective"].compute()) if hasattr(corr_dict["n_effective"], "compute") else float(corr_dict["n_effective"])

# Difference metrics
metrics = compute_difference_metrics(
    recon_9120, cpc_9120, dim=None, weighting="auto"
)

def _scalar(x):
    return float(x.compute()) if hasattr(x, "compute") else float(x)

rmse_value = _scalar(metrics["rmse"])
max_abs_error_metric = _scalar(metrics["max_abs_error"])

# Direct difference for the location of the largest deviation
diff_recon_cpc = (recon_9120 - cpc_9120).compute()
abs_diff = np.abs(diff_recon_cpc)
imax = int(abs_diff.argmax())
max_abs_error = float(abs_diff.isel(time=imax))
max_abs_time = pd.Timestamp(diff_recon_cpc.time.isel(time=imax).values)
max_abs_year = max_abs_time.year
max_abs_month = max_abs_time.month

# Baseline offsets over 1950-01 to 2020-12, unweighted in time
recon_offset_series = (recon_8110 - recon_9120).sel(time=slice("1950-01", "2020-12"))
recon_offset = compute_basic_stat(
    recon_offset_series, op="mean", dim="time"
).compute()
recon_offset_value = float(recon_offset)

cpc_offset_series = (cpc_8110 - cpc_9120).sel(time=slice("1950-01", "2020-12"))
cpc_offset = compute_basic_stat(
    cpc_offset_series, op="mean", dim="time"
).compute()
cpc_offset_value = float(cpc_offset)

offset_agreement = recon_offset_value - cpc_offset_value

print(f"RECON_9120 vs CPC ersst5_91-20 Pearson correlation: {corr_value:.6f}")
print(f"Effective sample size: {corr_n_eff:.1f}; p-value: {corr_pvalue:.3g}")
print(f"RMSE: {rmse_value:.4f} degC")
print(f"Largest absolute deviation: {max_abs_error:.4f} degC at {max_abs_year}-{max_abs_month:02d}")
print(f"RECON_9120 El Nino count (3-month mean >= 0.5 degC for >=5 steps): {n_recon_el_nino}")
print("Top 3 RECON_9120 El Nino events:")
print(recon_top3.to_string(index=False))
print(f"CPC ONI El Nino count (3-month mean >= 0.5 degC for >=5 steps): {n_oni_el_nino}")
print("Top 3 CPC ONI El Nino events:")
print(oni_top3.to_string(index=False))
print(f"Reconstructed baseline offset (1981-2010 minus 1991-2020): {recon_offset_value:.3f} degC")
print(f"Published baseline offset (1981-2010 minus 1991-2020): {cpc_offset_value:.3f} degC")
print(f"Reconstructed minus published offset: {offset_agreement:.3f} degC")
print("Counts can differ because ONI is CPC's operational 3-month running index on its own baseline, while RECON_9120 is an area-weighted ERSSTv5 anomaly on a fixed 1991-2020 climatology; near-threshold events may gain or lose a fifth consecutive season.")
