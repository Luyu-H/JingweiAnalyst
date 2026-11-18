# Load raw ERSSTv5 SST for the Nino 3.4 region
sst_nino34 = load_regular_data(
    dataset="ersst",
    variable="sst",
    time_range=["1950-01", "2023-12"],
    lat_range=[-5.0, 5.0],
    lon_range=[-170.0, -120.0],
)

# Reconstruct the Nino 3.4 anomaly relative to 1991-2020
anom_9120 = compute_anomaly_vs_own_climatology(
    data=sst_nino34,
    freq="month",
    base_period=["1991-01", "2020-12"],
    method="subtract",
)
recon_9120 = compute_basic_stat(
    data=anom_9120,
    op="area_weighted_mean",
    dim=["lat", "lon"],
).compute()
recon_9120.name = "RECON_9120"

# Same reconstruction relative to 1981-2010
anom_8110 = compute_anomaly_vs_own_climatology(
    data=sst_nino34,
    freq="month",
    base_period=["1981-01", "2010-12"],
    method="subtract",
)
recon_8110 = compute_basic_stat(
    data=anom_8110,
    op="area_weighted_mean",
    dim=["lat", "lon"],
).compute()
recon_8110.name = "RECON_8110"

# CPC published Nino 3.4 anomaly series
cpc_9120 = load_regular_data(
    dataset="cpc_indices",
    variable="NINO3.4",
    product="ersst5_91-20",
    time_range=["1950-01", "2023-12"],
).compute()
cpc_9120.name = "CPC_9120"

cpc_8110 = load_regular_data(
    dataset="cpc_indices",
    variable="NINO3.4",
    product="ersst5_81-10",
    time_range=["1950-01", "2023-12"],
).compute()
cpc_8110.name = "CPC_8110"

# CPC ONI (already a 3-month running index)
oni = load_regular_data(
    dataset="cpc_indices",
    variable="NINO3.4",
    product="oni",
    time_range=["1950-01", "2023-12"],
).compute()
oni.name = "ONI"

# Correlation and error metrics against CPC 1991-2020
corr_result = analyze_correlation(
    var1=recon_9120,
    var2=cpc_9120,
    dim="time",
    method="pearson",
)
pearson_r = float(corr_result["correlation"].compute())
pearson_r_rounded = round(pearson_r, 6)

diff_metrics = compute_difference_metrics(
    data=recon_9120,
    reference=cpc_9120,
    dim="time",
    weighting="none",
)
rmse = float(diff_metrics["rmse"].compute())
rmse_rounded = round(rmse, 4)
max_abs_error = float(diff_metrics["max_abs_error"].compute())
max_abs_error_rounded = round(max_abs_error, 4)

diff_series = diff_metrics["difference"].compute()
abs_diff_series = np.abs(diff_series)
max_idx = int(abs_diff_series.argmax(dim="time").compute())
max_time = pd.Timestamp(diff_series["time"].isel(time=max_idx).values)
max_year = int(max_time.year)
max_month = int(max_time.month)

# 3-month centred running mean and event detection
run3_recon = temporal_aggregate(
    data=recon_9120,
    freq="rolling",
    op="mean",
    window=3,
    center=True,
).compute()
run3_recon.name = "RECON_9120_3mo"

elnino_recon = detect_threshold_events(
    data=run3_recon,
    threshold=0.5,
    direction="above",
    min_duration=5,
)
lanina_recon = detect_threshold_events(
    data=run3_recon,
    threshold=-0.5,
    direction="below",
    min_duration=5,
)

n_elnino_recon = int(len(elnino_recon))
n_lanina_recon = int(len(lanina_recon))

elnino_recon_top = elnino_recon.sort_values("peak_value", ascending=False).head(3).copy()
elnino_recon_top["peak_time"] = pd.to_datetime(elnino_recon_top["peak_time"])
elnino_recon_top["peak_year"] = elnino_recon_top["peak_time"].dt.year
elnino_recon_top["peak_month"] = elnino_recon_top["peak_time"].dt.month
recon_top3 = elnino_recon_top[["peak_year", "peak_month", "peak_value"]].copy()
recon_top3["peak_value"] = recon_top3["peak_value"].round(2)
recon_top3 = recon_top3.rename(
    columns={"peak_year": "year", "peak_month": "month", "peak_value": "peak_degC"}
).reset_index(drop=True)

# ONI event detection under the identical criterion
elnino_oni = detect_threshold_events(
    data=oni,
    threshold=0.5,
    direction="above",
    min_duration=5,
)
n_elnino_oni = int(len(elnino_oni))

elnino_oni_top = elnino_oni.sort_values("peak_value", ascending=False).head(3).copy()
elnino_oni_top["peak_time"] = pd.to_datetime(elnino_oni_top["peak_time"])
elnino_oni_top["peak_year"] = elnino_oni_top["peak_time"].dt.year
elnino_oni_top["peak_month"] = elnino_oni_top["peak_time"].dt.month
oni_top3 = elnino_oni_top[["peak_year", "peak_month", "peak_value"]].copy()
oni_top3["peak_value"] = oni_top3["peak_value"].round(2)
oni_top3 = oni_top3.rename(
    columns={"peak_year": "year", "peak_month": "month", "peak_value": "peak_degC"}
).reset_index(drop=True)

# Baseline offsets over 1950-01 to 2020-12
recon_8110_2020 = recon_8110.sel(time=slice("1950-01", "2020-12"))
recon_9120_2020 = recon_9120.sel(time=slice("1950-01", "2020-12"))
offset_recon_metrics = compute_difference_metrics(
    data=recon_8110_2020,
    reference=recon_9120_2020,
    dim="time",
    weighting="none",
)
recon_baseline_offset = float(offset_recon_metrics["bias"].compute())
recon_baseline_offset_rounded = round(recon_baseline_offset, 3)

cpc_8110_2020 = cpc_8110.sel(time=slice("1950-01", "2020-12"))
cpc_9120_2020 = cpc_9120.sel(time=slice("1950-01", "2020-12"))
offset_cpc_metrics = compute_difference_metrics(
    data=cpc_8110_2020,
    reference=cpc_9120_2020,
    dim="time",
    weighting="none",
)
cpc_baseline_offset = float(offset_cpc_metrics["bias"].compute())
cpc_baseline_offset_rounded = round(cpc_baseline_offset, 3)

offset_difference = recon_baseline_offset - cpc_baseline_offset
offset_difference_rounded = round(offset_difference, 3)
abs_offset_difference_rounded = round(abs(offset_difference), 3)

explanation = (
    "RECON and CPC/ONI El Nino counts can differ because CPC's published index uses "
    "its own product version, quality control, baseline, and spatial averaging; small "
    "differences near the +0.5 degC threshold or near the 5-season duration can move "
    "a run in or out of the event set."
)

print(f"Pearson correlation (RECON_9120 vs CPC 1991-2020): {pearson_r_rounded:.6f}")
print(f"RMSE: {rmse_rounded:.4f} degC")
print(
    f"Largest absolute deviation: {max_abs_error_rounded:.4f} degC "
    f"at {max_year}-{max_month:02d}"
)
print(f"RECON_9120 El Nino events (>= +0.5 degC, >=5 seasons): {n_elnino_recon}")
print(recon_top3.to_string(index=False))
print(f"RECON_9120 La Nina events (<= -0.5 degC, >=5 seasons): {n_lanina_recon}")
print(f"CPC ONI El Nino events (>= +0.5 degC, >=5 seasons): {n_elnino_oni}")
print(oni_top3.to_string(index=False))
print(explanation)
print(
    f"Reconstructed baseline offset (RECON_8110 - RECON_9120, 1950-2020): "
    f"{recon_baseline_offset_rounded:.3f} degC"
)
print(
    f"Published baseline offset (CPC 1981-2010 - CPC 1991-2020, 1950-2020): "
    f"{cpc_baseline_offset_rounded:.3f} degC"
)
print(
    f"Offset difference (reconstructed - published): "
    f"{offset_difference_rounded:.3f} degC (abs {abs_offset_difference_rounded:.3f} degC)"
)
