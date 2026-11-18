# Rebuild the Nino 3.4 index from raw ERSSTv5 SST and compare with CPC products.

def _as_float(x):
    if isinstance(x, xr.DataArray):
        return float(x.compute().values)
    return float(x)

# ------------------------------------------------------------
# 1. Raw ERSSTv5 total SST over the Nino 3.4 box
# ------------------------------------------------------------
sst_full = load_regular_data(
    dataset="ersst",
    variable="sst",
    field="total",
    time_range=["1950-01", "2023-12"],
    lat_range=[-5.0, 5.0],
    lon_range=[-170.0, -120.0],
)

# Area-weighted regional mean SST, one value per month
sst_series = compute_basic_stat(
    data=sst_full,
    op="area_weighted_mean",
    dim=["lat", "lon"],
).compute()
sst_series.name = "nino34_sst"
sst_series = sst_series.assign_attrs(
    units="degC",
    long_name="Nino 3.4 area-weighted SST",
)

# ------------------------------------------------------------
# 2. Monthly climatologies from the two CPC baseline periods
# ------------------------------------------------------------
clim_9120 = temporal_aggregate(
    data=sst_series.sel(time=slice("1991-01", "2020-12")),
    freq="climatology",
    op="mean",
).compute()

clim_8110 = temporal_aggregate(
    data=sst_series.sel(time=slice("1981-01", "2010-12")),
    freq="climatology",
    op="mean",
).compute()

# ------------------------------------------------------------
# 3. Reconstructed Nino 3.4 anomalies on each baseline
# ------------------------------------------------------------
recon_9120 = (sst_series.groupby("time.month") - clim_9120).compute()
recon_9120.name = "RECON_9120"
recon_9120 = recon_9120.assign_attrs(
    units="degC",
    long_name="Nino 3.4 SST anomaly, 1991-2020 baseline",
)

recon_8110 = (sst_series.groupby("time.month") - clim_8110).compute()
recon_8110.name = "RECON_8110"
recon_8110 = recon_8110.assign_attrs(
    units="degC",
    long_name="Nino 3.4 SST anomaly, 1981-2010 baseline",
)

# ------------------------------------------------------------
# 4. CPC published Nino 3.4 indices
# ------------------------------------------------------------
cpc_9120 = load_regular_data(
    dataset="cpc_indices",
    variable="NINO3.4",
    product="ersst5_91-20",
    time_range=["1950-01", "2023-12"],
).squeeze(drop=True).compute()
cpc_9120.name = "CPC_9120"
cpc_9120 = cpc_9120.assign_attrs(units="degC")

cpc_8110 = load_regular_data(
    dataset="cpc_indices",
    variable="NINO3.4",
    product="ersst5_81-10",
    time_range=["1950-01", "2023-12"],
).squeeze(drop=True).compute()
cpc_8110.name = "CPC_8110"
cpc_8110 = cpc_8110.assign_attrs(units="degC")

cpc_oni = load_regular_data(
    dataset="cpc_indices",
    variable="NINO3.4",
    product="oni",
    time_range=["1950-01", "2023-12"],
).squeeze(drop=True).compute()
cpc_oni.name = "CPC_ONI"
cpc_oni = cpc_oni.assign_attrs(units="degC")

# ------------------------------------------------------------
# 5. Align reconstruction and CPC 1991-2020 series
# ------------------------------------------------------------
recon_9120_al, cpc_9120_al = xr.align(recon_9120, cpc_9120, join="inner")

# Pearson correlation
corr_res = analyze_correlation(
    var1=recon_9120_al,
    var2=cpc_9120_al,
    dim="time",
    method="pearson",
)
pearson_r = _as_float(corr_res["correlation"])

# Difference metrics, including RMSE and largest absolute deviation
diff_metrics = compute_difference_metrics(
    data=recon_9120_al,
    reference=cpc_9120_al,
    dim="time",
    weighting="none",
)
rmse = _as_float(diff_metrics["rmse"])
bias = _as_float(diff_metrics["bias"])

diff_field = diff_metrics["difference"].compute()
diff_field = diff_field.assign_attrs(units="degC")

max_abs_dev = _as_float(np.abs(diff_field).max())
max_abs_time = pd.Timestamp(diff_field.abs().idxmax(dim="time").values)
max_abs_year = int(max_abs_time.year)
max_abs_month = int(max_abs_time.month)

# ------------------------------------------------------------
# 6. 3-month centred running mean and event detection
# ------------------------------------------------------------
recon_run3 = temporal_aggregate(
    data=recon_9120,
    freq="rolling",
    window=3,
    center=True,
    op="mean",
).compute()
recon_run3.name = "RECON_9120_3mo_running_mean"

el_nino_recon = detect_threshold_events(
    data=recon_run3,
    threshold=0.5,
    direction="above",
    min_duration=5,
)
n_el_nino_recon = int(len(el_nino_recon))
top3_el_nino_recon = el_nino_recon.sort_values("peak_value", ascending=False).head(3).copy()

la_nina_recon = detect_threshold_events(
    data=recon_run3,
    threshold=-0.5,
    direction="below",
    min_duration=5,
)
n_la_nina_recon = int(len(la_nina_recon))

# CPC ONI is already the official 3-month running index
el_nino_oni = detect_threshold_events(
    data=cpc_oni,
    threshold=0.5,
    direction="above",
    min_duration=5,
)
n_el_nino_oni = int(len(el_nino_oni))
top3_el_nino_oni = el_nino_oni.sort_values("peak_value", ascending=False).head(3).copy()

# Format event tables for reporting
def _event_table(df):
    out = df.copy()
    out["peak_time"] = pd.to_datetime(out["peak_time"])
    out["peak_year"] = out["peak_time"].dt.year
    out["peak_month"] = out["peak_time"].dt.month
    out["peak_value"] = out["peak_value"].astype(float).round(2)
    return out

top3_el_nino_recon_tbl = _event_table(top3_el_nino_recon)
top3_el_nino_oni_tbl = _event_table(top3_el_nino_oni)

# ------------------------------------------------------------
# 7. Baseline offsets, 1950-01 to 2020-12, unweighted in time
# ------------------------------------------------------------
offset_recon_da = compute_basic_stat(
    data=(recon_8110 - recon_9120).sel(time=slice("1950-01", "2020-12")),
    op="mean",
    dim="time",
)
offset_recon = _as_float(offset_recon_da)

cpc_8110_al, cpc_9120_al_offset = xr.align(cpc_8110, cpc_9120, join="inner")
offset_cpc_da = compute_basic_stat(
    data=(cpc_8110_al - cpc_9120_al_offset).sel(time=slice("1950-01", "2020-12")),
    op="mean",
    dim="time",
)
offset_cpc = _as_float(offset_cpc_da)

offset_gap = offset_recon - offset_cpc
offset_abs_gap = abs(offset_gap)

offset_recon_3 = round(offset_recon, 3)
offset_cpc_3 = round(offset_cpc, 3)
offset_gap_3 = round(offset_gap, 3)

count_explanation = (
    "The counts differ because CPC ONI is the official 3-month running mean of CPC's own "
    "ERSSTv5 Nino 3.4 anomaly product, computed with CPC's base-period and processing choices. "
    "RECON_9120 is rebuilt here from raw ERSSTv5 SST using an area-weighted average over "
    "5S-5N, 170W-120W and a fixed 1991-2020 monthly climatology. Even under the same "
    "+0.5 degC and 5-season rule, small differences in base period, spatial averaging, and "
    "monthly anomaly construction can move a marginal event across the threshold or change "
    "its duration."
)

# ------------------------------------------------------------
# 8. Compact result summary
# ------------------------------------------------------------
results_summary = {
    "pearson_r": pearson_r,
    "rmse": rmse,
    "bias": bias,
    "max_abs_dev": max_abs_dev,
    "max_abs_year": max_abs_year,
    "max_abs_month": max_abs_month,
    "n_el_nino_recon": n_el_nino_recon,
    "top3_el_nino_recon": top3_el_nino_recon_tbl,
    "n_la_nina_recon": n_la_nina_recon,
    "n_el_nino_oni": n_el_nino_oni,
    "top3_el_nino_oni": top3_el_nino_oni_tbl,
    "offset_recon": offset_recon,
    "offset_cpc": offset_cpc,
    "offset_gap": offset_gap,
    "count_explanation": count_explanation,
}
