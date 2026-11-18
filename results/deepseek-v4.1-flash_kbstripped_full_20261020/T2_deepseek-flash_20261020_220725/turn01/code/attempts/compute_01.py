import numpy as np
import pandas as pd
import xarray as xr

# ---------------------------------------------------------------------------
# Step 1 - ERSSTv5 monthly SST over the Nino 3.4 box (5S-5N, 170W-120W)
# ---------------------------------------------------------------------------
sst_n34 = load_regular_data(
    dataset="ersst",
    variable="sst",
    time_range=["1950-01", "2023-12"],
    lat_range=[-5.0, 5.0],
    lon_range=[-170.0, -120.0],
)

# ---------------------------------------------------------------------------
# Step 2 - area-weighted box mean -> a single monthly regional SST series
# ---------------------------------------------------------------------------
n34_sst_series = compute_basic_stat(
    data=sst_n34, op="area_weighted_mean", dim=["lat", "lon"]
).compute()

# ---------------------------------------------------------------------------
# Step 3 - two reconstructions, each on its own monthly climatology
# ---------------------------------------------------------------------------
recon_9120 = compute_anomaly_vs_own_climatology(
    data=n34_sst_series, freq="month", base_period=["1991-01", "2020-12"]
).compute()

recon_8110 = compute_anomaly_vs_own_climatology(
    data=n34_sst_series, freq="month", base_period=["1981-01", "2010-12"]
).compute()

# ---------------------------------------------------------------------------
# Step 4 - the published CPC series
# ---------------------------------------------------------------------------
cpc_9120 = load_regular_data(
    dataset="cpc_indices", variable="NINO3.4", product="ersst5_91-20",
    time_range=["1950-01", "2023-12"],
).compute()

cpc_8110 = load_regular_data(
    dataset="cpc_indices", variable="NINO3.4", product="ersst5_81-10",
    time_range=["1950-01", "2020-12"],     # this control series ends in 2020
).compute()

oni = load_regular_data(
    dataset="cpc_indices", variable="NINO3.4", product="oni",
    time_range=["1950-01", "2023-12"],
).compute()

# ---------------------------------------------------------------------------
# Step 5 - correlation of RECON_9120 with the CPC 1991-2020 index
# ---------------------------------------------------------------------------
corr_out = analyze_correlation(
    var1=recon_9120, var2=cpc_9120, dim="time", method="pearson"
)
correlation = float(corr_out["correlation"].compute())
n_pairs = int(corr_out["n_samples"].compute())

# ---------------------------------------------------------------------------
# Step 6 - difference metrics
# ---------------------------------------------------------------------------
dm = compute_difference_metrics(data=recon_9120, reference=cpc_9120, dim="time")
rmse = float(dm["rmse"])
bias = float(dm["bias"])
mae = float(dm["mae"])
centered_rmse = float(dm["centered_rmse"])
max_abs_error = float(dm["max_abs_error"])

diff_series = dm["difference"]
if hasattr(diff_series, "compute"):
    diff_series = diff_series.compute()
abs_diff = abs(diff_series)
peak_abs_time = pd.Timestamp(abs_diff.idxmax(dim="time").values)
max_abs_error_year = int(peak_abs_time.year)
max_abs_error_month = int(peak_abs_time.month)

# ---------------------------------------------------------------------------
# Step 7 - 3-month centred running means
# (the CPC 'oni' product is already the official 3-month running mean, so it
#  is used as published rather than smoothed a second time)
# ---------------------------------------------------------------------------
rm_recon_9120 = temporal_aggregate(
    data=recon_9120, freq="rolling", window=3, center=True
).compute()
rm_recon_8110 = temporal_aggregate(
    data=recon_8110, freq="rolling", window=3, center=True
).compute()
rm_oni = oni

# ---------------------------------------------------------------------------
# Step 8 - operational ONI criterion (>= +0.5 / <= -0.5 for >= 5 seasons)
# ---------------------------------------------------------------------------
recon_elnino = detect_threshold_events(
    data=rm_recon_9120, threshold=0.5, direction="above", min_duration=5
)
recon_lanina = detect_threshold_events(
    data=rm_recon_9120, threshold=-0.5, direction="below", min_duration=5
)
oni_elnino = detect_threshold_events(
    data=rm_oni, threshold=0.5, direction="above", min_duration=5
)
oni_lanina = detect_threshold_events(
    data=rm_oni, threshold=-0.5, direction="below", min_duration=5
)

recon_elnino_count = int(len(recon_elnino))
oni_elnino_count = int(len(oni_elnino))
recon_lanina_count = int(len(recon_lanina))
oni_lanina_count = int(len(oni_lanina))


def top3_elnino(df):
    """Three strongest El Nino episodes: peak month/year and peak value."""
    out = []
    for _, row in df.sort_values("peak_value", ascending=False).head(3).iterrows():
        ts = pd.Timestamp(row["peak_time"])
        out.append(
            dict(year=int(ts.year), month=int(ts.month),
                 peak_value=round(float(row["peak_value"]), 2),
                 start=str(pd.Timestamp(row["start"]).date()),
                 end=str(pd.Timestamp(row["end"]).date()),
                 duration=int(row["duration"]))
        )
    return out


recon_top3 = top3_elnino(recon_elnino)
oni_top3 = top3_elnino(oni_elnino)

# ---------------------------------------------------------------------------
# Steps 9-10 - baseline offsets over 1950-01 .. 2020-12
# ---------------------------------------------------------------------------
recon_offset_series = (recon_8110 - recon_9120).sel(time=slice("1950-01", "2020-12"))
offset_recon = float(
    compute_basic_stat(data=recon_offset_series, op="mean", dim="time").compute()
)

cpc_offset_series = (cpc_8110 - cpc_9120).sel(time=slice("1950-01", "2020-12"))
offset_cpc = float(
    compute_basic_stat(data=cpc_offset_series, op="mean", dim="time").compute()
)

offset_difference = offset_recon - offset_cpc          # reconstructed minus published

# ---------------------------------------------------------------------------
# Assemble the reported numbers
# ---------------------------------------------------------------------------
summary = {
    "correlation_recon9120_vs_cpc9120": round(correlation, 6),
    "n_months_correlated": n_pairs,
    "rmse_degC": round(rmse, 4),
    "bias_degC": round(bias, 4),
    "mae_degC": round(mae, 4),
    "centered_rmse_degC": round(centered_rmse, 4),
    "max_abs_error_degC": round(max_abs_error, 4),
    "max_abs_error_year": max_abs_error_year,
    "max_abs_error_month": max_abs_error_month,
    "recon_elnino_count_1950_2023": recon_elnino_count,
    "recon_elnino_top3": recon_top3,
    "oni_elnino_count_1950_2023": oni_elnino_count,
    "oni_elnino_top3": oni_top3,
    "recon_lanina_count": recon_lanina_count,
    "oni_lanina_count": oni_lanina_count,
    "recon_baseline_offset_degC": round(offset_recon, 3),
    "cpc_baseline_offset_degC": round(offset_cpc, 3),
    "offset_difference_degC": round(offset_difference, 3),
}

for k, v in summary.items():
    print(f"{k}: {v}")
