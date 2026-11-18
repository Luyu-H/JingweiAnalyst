# Load SST data for 1950-01 to 2014-12
time_range = ["1950-01", "2014-12"]
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

ersst_sst = load_regular_data(
    dataset="ersst",
    variable="sst",
    time_range=time_range,
    chunk_time=12,
)

model_sst_native = {
    model: load_regular_data(
        dataset="cmip6_historical",
        variable="tos",
        model=model,
        time_range=time_range,
        chunk_time=12,
    )
    for model in models
}

published_pdo = load_regular_data(
    dataset="cpc_indices",
    product="pdo",
    time_range=time_range,
    chunk_time=12,
)

# Regrid all SST fields to a common 1-degree grid with bilinear interpolation
ersst_rg = regrid_horizontal(ersst_sst, target="1deg", method="bilinear")
model_sst_rg = {
    model: regrid_horizontal(model_sst_native[model], target="1deg", method="bilinear")
    for model in models
}

# Build monthly anomalies relative to each source's own 1950-2014 monthly climatology
base_period = ["1950-01", "2014-12"]
ersst_anom = compute_anomaly_vs_own_climatology(
    ersst_rg, freq="month", base_period=base_period, method="subtract"
)
model_anom = {
    model: compute_anomaly_vs_own_climatology(
        model_sst_rg[model], freq="month", base_period=base_period, method="subtract"
    )
    for model in models
}

# Compute monthly global-mean SST anomalies from the regridded full fields
ersst_gm = compute_basic_stat(ersst_anom, op="area_weighted_mean", dim=["lat", "lon"])
model_gm = {
    model: compute_basic_stat(model_anom[model], op="area_weighted_mean", dim=["lat", "lon"])
    for model in models
}

# Remove global-mean anomaly before PDO calculation
ersst_anom_demeaned = (ersst_anom - ersst_gm).assign_attrs(units=ersst_anom.attrs.get("units", "degC"))
model_anom_demeaned = {
    model: (model_anom[model] - model_gm[model]).assign_attrs(units=model_anom[model].attrs.get("units", "degC"))
    for model in models
}

# Compute PDO index, EOF pattern, explained variance, regression pattern
ersst_pdo_idx, ersst_eof, ersst_expl_var, ersst_reg = compute_pdo_index(
    ersst_anom_demeaned, standardize=True
)

model_pdo = {}
for model in models:
    model_pdo[model] = compute_pdo_index(
        model_anom_demeaned[model], standardize=True
    )

# Materialize main PDO outputs
ersst_pdo_idx = ersst_pdo_idx.compute()
ersst_expl_var_value = float(ersst_expl_var)
ersst_reg = ersst_reg.compute()

model_pdo_idx = {}
model_expl_var_value = {}
model_reg = {}
for model in models:
    idx, eofp, expl, regp = model_pdo[model]
    model_pdo_idx[model] = idx.compute()
    model_expl_var_value[model] = float(expl)
    model_reg[model] = regp.compute()

published_pdo = published_pdo.compute()

# Restrict regression patterns to the PDO domain for reporting/plotting
pdo_lat_range = [20, 70]
pdo_lon_range = [110, 260]  # 110E to 100W on 0-360 plotting axis; data remain on native opened convention

ersst_reg_np = ersst_reg.where(
    (ersst_reg["lat"] >= pdo_lat_range[0]) & (ersst_reg["lat"] <= pdo_lat_range[1]),
    drop=True
)
ersst_reg_np = ersst_reg_np.where(
    ((ersst_reg_np["lon"] % 360) >= pdo_lon_range[0]) & ((ersst_reg_np["lon"] % 360) <= pdo_lon_range[1]),
    drop=True
)

model_reg_np = {}
for model in models:
    da = model_reg[model]
    da = da.where((da["lat"] >= pdo_lat_range[0]) & (da["lat"] <= pdo_lat_range[1]), drop=True)
    da = da.where(((da["lon"] % 360) >= pdo_lon_range[0]) & ((da["lon"] % 360) <= pdo_lon_range[1]), drop=True)
    model_reg_np[model] = da

# Amplitude of the mode from the regression pattern, in degC per std
ersst_amp = float(
    np.sqrt(
        compute_basic_stat((ersst_reg_np ** 2), op="area_weighted_mean", dim=["lat", "lon"]).compute().item()
    )
)

model_amp = {}
for model in models:
    model_amp[model] = float(
        np.sqrt(
            compute_basic_stat((model_reg_np[model] ** 2), op="area_weighted_mean", dim=["lat", "lon"]).compute().item()
        )
    )

# Build one common comparison mask across observation and all model regression patterns
common_mask = build_comparison_mask(
    fields=[ersst_reg_np] + [model_reg_np[m] for m in models],
    ocean_only=True,
    reduce_other_dims="all",
).compute()

ersst_reg_masked = ersst_reg_np.where(common_mask)
model_reg_masked = {model: model_reg_np[model].where(common_mask) for model in models}

# Spatial Pearson correlation of each model regression pattern with ERSSTv5 regression pattern
pattern_corr = {}
for model in models:
    corr_dict = analyze_correlation(
        var1=model_reg_masked[model],
        var2=ersst_reg_masked,
        dim=["lat", "lon"],
        method="pearson",
        weighting="area",
        min_periods=10,
    )
    pattern_corr[model] = float(corr_dict["correlation"].compute().item())

# 12-month autocorrelation of each standardized monthly PDO index
def lag12_autocorr(series):
    out = analyze_correlation(
        var1=series,
        var2=series,
        dim="time",
        method="pearson",
        lag=12,
        detrend=False,
        min_periods=10,
    )
    return {
        "correlation": float(out["correlation"].compute().item()),
        "pvalue": float(out["pvalue"].compute().item()),
        "pvalue_unadjusted": float(out["pvalue_unadjusted"].compute().item()),
        "n_samples": int(out["n_samples"].compute().item()),
        "n_effective": float(out["n_effective"].compute().item()),
    }

ersst_acf12 = lag12_autocorr(ersst_pdo_idx)
model_acf12 = {model: lag12_autocorr(model_pdo_idx[model]) for model in models}

# Correlation between ERSSTv5-derived PDO and published PDO over 1950-2014
ersst_vs_published = analyze_correlation(
    var1=ersst_pdo_idx,
    var2=published_pdo,
    dim="time",
    method="pearson",
    detrend=False,
    min_periods=10,
)
ersst_vs_published_corr = float(ersst_vs_published["correlation"].compute().item())
ersst_vs_published_p = float(ersst_vs_published["pvalue"].compute().item())
ersst_vs_published_p_naive = float(ersst_vs_published["pvalue_unadjusted"].compute().item())
ersst_vs_published_n = int(ersst_vs_published["n_samples"].compute().item())
ersst_vs_published_neff = float(ersst_vs_published["n_effective"].compute().item())

# 9-year (108-month) centred running mean of ERSSTv5 PDO index
ersst_pdo_9yr = temporal_aggregate(
    ersst_pdo_idx, freq="rolling", op="mean", window=108, center=True
).compute()

# Identify sign-change years; merge any crossing within 3 years of prior counted transition
sm = ersst_pdo_9yr
sm_vals = sm.values
sm_time = pd.to_datetime(sm["time"].values)

valid = np.isfinite(sm_vals)
cross_years_raw = []
for i in range(1, len(sm_vals)):
    if not (valid[i - 1] and valid[i]):
        continue
    y1 = sm_vals[i - 1]
    y2 = sm_vals[i]
    if y1 == 0:
        cross_years_raw.append(int(sm_time[i - 1].year))
    elif y1 * y2 < 0:
        cross_time = sm_time[i]
        cross_years_raw.append(int(cross_time.year))

transition_years = []
for yr in cross_years_raw:
    if not transition_years:
        transition_years.append(yr)
    elif yr - transition_years[-1] < 3:
        continue
    else:
        transition_years.append(yr)

documented_shift_years = []
for yr in transition_years:
    if min(abs(yr - 1947), abs(yr - 1976), abs(yr - 1988), abs(yr - 1998)) <= 1:
        documented_shift_years.append(yr)

# Determine best-matching model by highest pattern correlation
best_model = max(pattern_corr, key=pattern_corr.get)

# Summaries
summary_rows = [{
    "source": "ERSSTv5",
    "explained_variance_pct": round(ersst_expl_var_value, 1),
    "amplitude_degC_per_std": round(ersst_amp, 3),
    "acf_12mo": round(ersst_acf12["correlation"], 3),
    "pattern_corr_vs_ersst": np.nan,
}]
for model in models:
    summary_rows.append({
        "source": model,
        "explained_variance_pct": round(model_expl_var_value[model], 1),
        "amplitude_degC_per_std": round(model_amp[model], 3),
        "acf_12mo": round(model_acf12[model]["correlation"], 3),
        "pattern_corr_vs_ersst": round(pattern_corr[model], 3),
    })
summary_table = pd.DataFrame(summary_rows)

ersst_vs_published_summary = {
    "correlation": round(ersst_vs_published_corr, 3),
    "pvalue_effective": ersst_vs_published_p,
    "pvalue_unadjusted": ersst_vs_published_p_naive,
    "n_samples": ersst_vs_published_n,
    "n_effective": ersst_vs_published_neff,
}

transition_summary = {
    "transition_years": transition_years,
    "documented_shift_matches": documented_shift_years,
}

# Qualitative consistency checks against observations
obs_amp = ersst_amp
obs_persist = ersst_acf12["correlation"]
model_assessment = {}
for model in models:
    amp_ratio = model_amp[model] / obs_amp if obs_amp != 0 else np.nan
    pers_diff = model_acf12[model]["correlation"] - obs_persist
    model_assessment[model] = {
        "pattern_corr_vs_obs": pattern_corr[model],
        "amp_ratio_to_obs": amp_ratio,
        "persistence_minus_obs": pers_diff,
    }

# Variables intended for direct downstream use / inspection:
# summary_table
# ersst_vs_published_summary
# transition_summary
# best_model
# model_assessment
# ersst_reg_np, model_reg_np, ersst_pdo_idx, model_pdo_idx, published_pdo, ersst_pdo_9yr
