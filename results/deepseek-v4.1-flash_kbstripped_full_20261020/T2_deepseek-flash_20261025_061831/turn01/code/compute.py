# Load ERSSTv5 monthly SST over the Nino 3.4 region
sst = load_regular_data(
    dataset="ersst",
    variable="sst",
    time_range=["1950-01", "2023-12"],
    lat_range=[-5, 5],
    lon_range=[-170, -120],
)

# Load the three CPC Nino 3.4 series
cpc_9120 = load_regular_data(
    dataset="cpc_indices",
    variable="NINO3.4",
    product="ersst5_91-20",
    time_range=["1950-01", "2023-12"],
    field="anomaly",
)
cpc_8110 = load_regular_data(
    dataset="cpc_indices",
    variable="NINO3.4",
    product="ersst5_81-10",
    time_range=["1950-01", "2023-12"],
    field="anomaly",
)
oni = load_regular_data(
    dataset="cpc_indices",
    variable="NINO3.4",
    product="oni",
    time_range=["1950-01", "2023-12"],
    field="anomaly",
)

# Build reconstructed anomalies on the two baselines
anom_9120 = compute_anomaly_vs_own_climatology(
    data=sst,
    freq="month",
    base_period=["1991-01", "2020-12"],
    method="subtract",
)
anom_8110 = compute_anomaly_vs_own_climatology(
    data=sst,
    freq="month",
    base_period=["1981-01", "2010-12"],
    method="subtract",
)

# Area-weighted spatial means over the Nino 3.4 box
recon_9120 = compute_basic_stat(
    data=anom_9120,
    op="area_weighted_mean",
    dim=["lat", "lon"],
).compute()
recon_8110 = compute_basic_stat(
    data=anom_8110,
    op="area_weighted_mean",
    dim=["lat", "lon"],
).compute()

# Materialise CPC series
cpc_9120 = cpc_9120.compute()
cpc_8110 = cpc_8110.compute()
oni = oni.compute()

# Align the main comparison pair only, preserving the full common record
recon_9120, cpc_9120 = xr.align(recon_9120, cpc_9120, join="inner")

# Correlation and difference metrics
corr = analyze_correlation(
    var1=recon_9120,
    var2=cpc_9120,
    dim="time",
    method="pearson",
    detrend=False,
)
r = float(corr["correlation"].compute())
pval = float(corr["pvalue"].compute())
pval_unadj = float(corr["pvalue_unadjusted"].compute())
n_eff = float(corr["n_effective"].compute())
n_samples = int(corr["n_samples"].compute())

metrics = compute_difference_metrics(
    data=recon_9120,
    reference=cpc_9120,
)
bias = float(metrics["bias"])
rmse = float(metrics["rmse"])
max_abs = float(metrics["max_abs_error"])

diff = metrics["difference"].compute()
abs_diff = np.abs(diff.values)
imax = int(np.nanargmax(abs_diff))
max_time = pd.Timestamp(diff.time.values[imax])
max_year = int(max_time.year)
max_month = int(max_time.month)

# 3-month running mean and event detection on the reconstruction
recon_run3 = temporal_aggregate(
    data=recon_9120,
    freq="rolling",
    window=3,
    center=True,
    op="mean",
).compute()

el_events = detect_threshold_events(
    data=recon_run3,
    threshold=0.5,
    direction="above",
    min_duration=5,
)
la_events = detect_threshold_events(
    data=recon_run3,
    threshold=-0.5,
    direction="below",
    min_duration=5,
)

el_nino_count = len(el_events)
la_nina_count = len(la_events)

el_top3 = el_events.sort_values("peak_value", ascending=False).head(3).copy()
el_top3["peak_year"] = pd.to_datetime(el_top3["peak_time"]).dt.year
el_top3["peak_month"] = pd.to_datetime(el_top3["peak_time"]).dt.month
el_top3_report = el_top3[["peak_year", "peak_month", "peak_value"]].copy()
el_top3_report["peak_value"] = el_top3_report["peak_value"].round(2)

# Event detection on CPC ONI
oni_events = detect_threshold_events(
    data=oni,
    threshold=0.5,
    direction="above",
    min_duration=5,
)
oni_el_nino_count = len(oni_events)
oni_top3 = oni_events.sort_values("peak_value", ascending=False).head(3).copy()
oni_top3["peak_year"] = pd.to_datetime(oni_top3["peak_time"]).dt.year
oni_top3["peak_month"] = pd.to_datetime(oni_top3["peak_time"]).dt.month
oni_top3_report = oni_top3[["peak_year", "peak_month", "peak_value"]].copy()
oni_top3_report["peak_value"] = oni_top3_report["peak_value"].round(2)

# Baseline offsets over 1950-01 to 2020-12 (unweighted in time)
recon_8110_off, recon_9120_off, cpc_8110_off, cpc_9120_off = xr.align(
    recon_8110, recon_9120, cpc_8110, cpc_9120, join="inner"
)

recon_offset_series = (recon_8110_off - recon_9120_off).sel(
    time=slice("1950-01", "2020-12")
)
recon_offset_mean = float(
    compute_basic_stat(
        data=recon_offset_series,
        op="mean",
        dim="time",
    ).compute()
)

cpc_offset_series = (cpc_8110_off - cpc_9120_off).sel(
    time=slice("1950-01", "2020-12")
)
cpc_offset_mean = float(
    compute_basic_stat(
        data=cpc_offset_series,
        op="mean",
        dim="time",
    ).compute()
)

offset_difference = recon_offset_mean - cpc_offset_mean

# Print summary
print(f"Pearson r = {r:.6f}")
print(f"p-value (effective N) = {pval:.6g}")
print(f"p-value (unadjusted) = {pval_unadj:.6g}")
print(f"Effective N = {n_eff:.2f}, samples = {n_samples}")
print(f"RMSE = {rmse:.4f} degC")
print(f"Max |deviation| = {max_abs:.4f} degC at {max_year}-{max_month:02d}")
print(f"Bias = {bias:.4f} degC")
print(f"RECON El Nino count = {el_nino_count}")
print("RECON top 3 El Nino peaks:")
print(el_top3_report.to_string(index=False))
print(f"RECON La Nina count = {la_nina_count}")
print(f"CPC ONI El Nino count = {oni_el_nino_count}")
print("CPC ONI top 3 El Nino peaks:")
print(oni_top3_report.to_string(index=False))
print(f"RECON baseline offset (8110-9120) = {recon_offset_mean:.3f} degC")
print(f"CPC baseline offset (8110-9120) = {cpc_offset_mean:.3f} degC")
print(f"Offset difference (RECON - CPC) = {offset_difference:.3f} degC")
