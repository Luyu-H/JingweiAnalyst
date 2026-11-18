# ---------------------------------------------------------------------------
# 1. Load the raw fields / published index series
# ---------------------------------------------------------------------------
sst_box = load_regular_data(
    dataset="ersst", variable="sst",
    time_range=["1950-01", "2023-12"],
    lat_range=[-6, 6], lon_range=[-172, -118],
)

cpc_9120_da = load_regular_data(dataset="cpc_indices", variable="NINO3.4",
                                product="ersst5_91-20",
                                time_range=["1950-01", "2023-12"])
cpc_8110_da = load_regular_data(dataset="cpc_indices", variable="NINO3.4",
                                product="ersst5_81-10",
                                time_range=["1950-01", "2023-12"])
oni_da = load_regular_data(dataset="cpc_indices", variable="NINO3.4",
                           product="oni",
                           time_range=["1950-01", "2023-12"])

# ---------------------------------------------------------------------------
# 2. Rebuild the Nino 3.4 index from ERSSTv5 with two different base periods
# ---------------------------------------------------------------------------
recon_9120_da = compute_enso_nino_index(
    data=sst_box, nino_region="nino3.4", anomaly_freq="month",
    climatology="self", base_period=["1991-01", "2020-12"],
)
recon_8110_da = compute_enso_nino_index(
    data=sst_box, nino_region="nino3.4", anomaly_freq="month",
    climatology="self", base_period=["1981-01", "2010-12"],
)

recon_9120_da = recon_9120_da.compute()
recon_8110_da = recon_8110_da.compute()
cpc_9120_da = cpc_9120_da.compute()
cpc_8110_da = cpc_8110_da.compute()
oni_da = oni_da.compute()


# --- helper: force every series onto a clean month-start axis and join -------
def _monthly_series(da, name):
    d = da.squeeze(drop=True)
    vals = np.asarray(d.values, dtype=float)
    idx = pd.DatetimeIndex(np.asarray(d["time"].values))
    s = pd.Series(vals, index=idx.to_period("M").to_timestamp(), name=name)
    return s


s_recon = _monthly_series(recon_9120_da, "recon_9120")
s_cpc = _monthly_series(cpc_9120_da, "cpc_9120")

df_cmp = pd.concat([s_recon, s_cpc], axis=1, join="inner").dropna()
df_cmp = df_cmp.loc["1950-01":"2023-12"]

recon_m = xr.DataArray(df_cmp["recon_9120"].values, dims="time",
                       coords={"time": df_cmp.index.values},
                       name="nino34_recon_9120")
cpc_m = xr.DataArray(df_cmp["cpc_9120"].values, dims="time",
                     coords={"time": df_cmp.index.values},
                     name="nino34_cpc_9120")

# ---------------------------------------------------------------------------
# 3. Agreement between the reconstruction and the published Nino 3.4 anomaly
# ---------------------------------------------------------------------------
corr_res = analyze_correlation(var1=recon_m, var2=cpc_m, dim="time",
                               method="pearson")
corr_r = float(corr_res["correlation"].compute())

diff_res = compute_difference_metrics(data=recon_m, reference=cpc_m, dim="time")
rmse = float(diff_res["rmse"])
max_abs_dev = float(diff_res["max_abs_error"])

# --- robust extraction of the location of the largest departure ------------
_at = diff_res["max_abs_error_at"]

if isinstance(_at, dict):
    # compute_difference_metrics returns the coordinates as a dict of {dim: coord}
    if "time" in _at:
        _at_val = _at["time"]
    else:
        _at_val = next(iter(_at.values()))
else:
    _at_val = _at
    # unwrap a scalar DataArray / numpy array down to a single value
    if hasattr(_at_val, "values") and not callable(getattr(_at_val, "values")):
        _at_val = np.asarray(_at_val.values).ravel()[0]

_ts = pd.Timestamp(_at_val)
max_dev_ym = f"{_ts.year}-{_ts.month:02d}"

# ---------------------------------------------------------------------------
# 4. 3-month centred running mean of the reconstruction + ONI-style events
# ---------------------------------------------------------------------------
recon_run = temporal_aggregate(data=recon_m, freq="rolling", window=3,
                               center=True, op="mean")

elnino_recon = detect_threshold_events(data=recon_run, threshold=0.5,
                                       direction="above", min_duration=5)
lanina_recon = detect_threshold_events(data=recon_run, threshold=-0.5,
                                       direction="below", min_duration=5)

n_elnino_recon = int(len(elnino_recon))
n_lanina_recon = int(len(lanina_recon))


def _top_events(ev, k=3):
    d = ev.copy()
    d["peak_time"] = pd.to_datetime(d["peak_time"])
    d = d.sort_values("peak_value", ascending=False).head(k)
    return pd.DataFrame({
        "peak_year_month": d["peak_time"].dt.strftime("%Y-%m").values,
        "peak_value_degC": np.round(d["peak_value"].values, 2),
        "duration_months": d["duration"].values,
    }).reset_index(drop=True)


recon_top3 = _top_events(elnino_recon)

# ---------------------------------------------------------------------------
# 5. Same criterion applied directly to the operational CPC ONI
# ---------------------------------------------------------------------------
oni = oni_da.squeeze(drop=True).dropna(dim="time")
elnino_oni = detect_threshold_events(data=oni, threshold=0.5,
                                     direction="above", min_duration=5)
n_elnino_oni = int(len(elnino_oni))
oni_top3 = _top_events(elnino_oni)

# ---------------------------------------------------------------------------
# 6. Baseline offsets: 1981-2010 minus 1991-2020
# ---------------------------------------------------------------------------
s_recon81 = _monthly_series(recon_8110_da, "recon_8110")
df_r = pd.concat([s_recon81, s_recon], axis=1, join="inner")
d_recon = (df_r["recon_8110"] - df_r["recon_9120"]).loc["1950-01":"2020-12"]
d_recon_da = xr.DataArray(d_recon.values, dims="time",
                          coords={"time": d_recon.index.values})
offset_recon = float(compute_basic_stat(data=d_recon_da, op="mean",
                                        dim="time").compute())

s_cpc81 = _monthly_series(cpc_8110_da, "cpc_8110")
df_c = pd.concat([s_cpc81, s_cpc], axis=1, join="inner")
d_cpc = (df_c["cpc_8110"] - df_c["cpc_9120"]).loc["1950-01":"2020-12"]
d_cpc_da = xr.DataArray(d_cpc.values, dims="time",
                        coords={"time": d_cpc.index.values})
offset_cpc = float(compute_basic_stat(data=d_cpc_da, op="mean",
                                      dim="time").compute())

offset_gap = offset_cpc - offset_recon

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
print("=== RECON_9120 vs CPC ersst5_91-20 (1950-01 to 2023-12) ===")
print(f"months compared          : {len(df_cmp)}")
print(f"Pearson correlation      : {corr_r:.6f}")
print(f"RMSE (degC)              : {rmse:.4f}")
print(f"Largest abs deviation    : {max_abs_dev:.4f} degC at {max_dev_ym}")
print()
print("=== El Nino events, RECON_9120 (3-month running mean, >= +0.5, >=5 seasons) ===")
print(f"count : {n_elnino_recon}")
print(recon_top3.to_string(index=False))
print(f"(La Nina events in RECON_9120 for shading: {n_lanina_recon})")
print()
print("=== El Nino events, CPC ONI (>= +0.5, >=5 seasons) ===")
print(f"count : {n_elnino_oni}")
print(oni_top3.to_string(index=False))
print()
print("=== Baseline offsets, mean 1950-01 to 2020-12 ===")
print(f"RECON 1981-2010 minus 1991-2020 : {offset_recon:.3f} degC")
print(f"CPC   1981-2010 minus 1991-2020 : {offset_cpc:.3f} degC")
print(f"difference (CPC - RECON)        : {offset_gap:.3f} degC")
