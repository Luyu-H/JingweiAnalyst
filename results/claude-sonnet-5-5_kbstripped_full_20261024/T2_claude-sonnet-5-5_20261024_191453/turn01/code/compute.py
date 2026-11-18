# ---------- 1. Load ERSSTv5 SST over a box around Nino 3.4 (native 2-degree grid) ----------
sst = load_regular_data(
    dataset="ersst", variable="sst",
    time_range=["1950-01", "2023-12"],
    lat_range=[-10, 10], lon_range=[-175, -115],
)

# ---------- 2. Rebuild the index on ERSSTv5's own climatology ----------
recon_9120 = compute_enso_nino_index(
    sst, nino_region="nino3.4", anomaly_freq="month",
    climatology="self", base_period=["1991-01", "2020-12"],
).compute()
recon_8110 = compute_enso_nino_index(
    sst, nino_region="nino3.4", anomaly_freq="month",
    climatology="self", base_period=["1981-01", "2010-12"],
).compute()

def norm_time(da):
    t = pd.DatetimeIndex(da["time"].values).to_period("M").to_timestamp()
    return da.assign_coords(time=t)

recon_9120 = norm_time(recon_9120)
recon_8110 = norm_time(recon_8110)

# ---------- 3. Published CPC series ----------
tr = ["1950-01", "2023-12"]
cpc_9120 = norm_time(load_regular_data(dataset="cpc_indices", variable="NINO3.4",
                                       product="ersst5_91-20", field="anomaly",
                                       time_range=tr).compute())
cpc_8110 = norm_time(load_regular_data(dataset="cpc_indices", variable="NINO3.4",
                                       product="ersst5_81-10", field="anomaly",
                                       time_range=tr).compute())
cpc_oni = norm_time(load_regular_data(dataset="cpc_indices", variable="NINO3.4",
                                      product="oni", field="anomaly",
                                      time_range=tr).compute())

# ---------- 4. RECON_9120 vs CPC 1991-2020 on a common mask ----------
r_a, c_a = xr.align(recon_9120, cpc_9120, join="inner")
valid = r_a.notnull() & c_a.notnull()
r_cmp = r_a.where(valid).dropna("time")
c_cmp = c_a.where(valid).dropna("time")

corr_res = analyze_correlation(r_cmp, c_cmp, dim="time", detrend=False)
pearson_r = float(corr_res["correlation"])
n_samples = int(corr_res["n_samples"])
n_effective = float(corr_res["n_effective"])
pvalue_adjusted = float(corr_res["pvalue"])

metrics = compute_difference_metrics(r_cmp, c_cmp, dim="time")
rmse = float(metrics["rmse"])
bias = float(metrics["bias"])
diff_series = metrics["difference"].compute()
imax = int(np.nanargmax(np.abs(diff_series.values)))
max_abs_dev = float(np.abs(diff_series.values[imax]))
max_dev_signed = float(diff_series.values[imax])
max_dev_time = pd.Timestamp(diff_series["time"].values[imax])
max_dev_ym = max_dev_time.strftime("%Y-%m")

print(f"Common months: {n_samples}, window {str(r_cmp.time.values[0])[:7]} to {str(r_cmp.time.values[-1])[:7]}")
print(f"Pearson r = {pearson_r:.6f} (n_eff = {n_effective:.1f}, autocorrelation-adjusted p = {pvalue_adjusted:.3g})")
print(f"RMSE = {rmse:.4f} degC, bias = {bias:.4f} degC")
print(f"Max |deviation| = {max_abs_dev:.4f} degC (signed {max_dev_signed:.4f}) at {max_dev_ym}")

# ---------- 5. El Nino / La Nina events ----------
recon_rm = temporal_aggregate(recon_9120, freq="rolling", window=3, center=True).compute()
oni = cpc_oni

def top3(df):
    d = df.sort_values("peak_value", ascending=False).head(3)
    return [(pd.Timestamp(r.peak_time).strftime("%Y-%m"), round(float(r.peak_value), 2))
            for r in d.itertuples()]

ev_recon = detect_threshold_events(recon_rm, 0.5, direction="above", min_duration=5)
ev_oni = detect_threshold_events(oni, 0.5, direction="above", min_duration=5)
ln_recon = detect_threshold_events(recon_rm, -0.5, direction="below", min_duration=5)

n_el_recon = len(ev_recon)
n_el_oni = len(ev_oni)
top3_recon = top3(ev_recon)
top3_oni = top3(ev_oni)
print(f"El Nino events RECON_9120 (3-mo running mean): {n_el_recon}; top 3 peaks: {top3_recon}")
print(f"El Nino events CPC ONI: {n_el_oni}; top 3 peaks: {top3_oni}")
print(f"La Nina episodes RECON_9120: {len(ln_recon)}")
print("Counts can differ because CPC ONI uses centred 30-year base periods updated every 5 years, "
      "while RECON_9120 uses a single fixed 1991-2020 baseline; early years are therefore "
      "relatively warmer/cooler in ONI, shifting values across the +/-0.5 threshold.")

# ---------- 6. Baseline offsets, 1950-01 to 2020-12 ----------
sl = slice("1950-01", "2020-12")
rec_off_series = (recon_8110.sel(time=sl) - recon_9120.sel(time=sl)).assign_attrs(units="degC")
c81, c91 = xr.align(cpc_8110.sel(time=sl), cpc_9120.sel(time=sl), join="inner")
pub_off_series = (c81 - c91).assign_attrs(units="degC")

rec_off_series, pub_off_series = xr.align(rec_off_series, pub_off_series, join="inner")
rec_offset = float(compute_basic_stat(rec_off_series, op="mean", dim="time").compute())
pub_offset = float(compute_basic_stat(pub_off_series, op="mean", dim="time").compute())

off_metrics = compute_difference_metrics(rec_off_series, pub_off_series, dim="time")
off_rmse = float(off_metrics["rmse"])
off_maxdiff = float(off_metrics["max_abs_error"])
offset_gap = rec_offset - pub_offset
print(f"Reconstructed offset (8110-9120) = {rec_offset:.3f} degC")
print(f"Published offset (CPC 81-10 minus 91-20) = {pub_offset:.3f} degC")
print(f"Difference of means = {offset_gap:.3f} degC; monthly-series RMSE = {off_rmse:.4f}, max |diff| = {off_maxdiff:.4f}")

# Arrays for plotting
t_all = r_cmp["time"].values
r_vals = r_cmp.values
c_vals = c_cmp.values
diff_vals = diff_series.values
rm_t = recon_rm["time"].values
rm_vals = recon_rm.values
