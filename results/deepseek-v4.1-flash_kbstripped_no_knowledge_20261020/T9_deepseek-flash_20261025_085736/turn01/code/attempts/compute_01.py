time_range = ["1950-01", "2014-12"]
target_grid = "1deg"
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
all_names = ["ERSSTv5"] + models

def _extract_scalar(x):
    return float(x.compute()) if hasattr(x, "compute") else float(x)

def _pdo_diagnostics(field):
    pdo_idx, eof, var_expl, reg = compute_pdo_index(
        data=field,
        depth_level=None,
        standardize=True,
    )
    pdo_idx_c = pdo_idx.compute()
    reg_c = reg.compute()
    var_expl_v = _extract_scalar(var_expl)

    # Mode amplitude: area-weighted RMS of the regression pattern, in degC per SD.
    reg2 = reg_c.copy(data=reg_c.values ** 2)
    rms2 = compute_basic_stat(data=reg2, op="area_weighted_mean", dim=["lat", "lon"])
    amp_rms = float(np.sqrt(_extract_scalar(rms2)))
    amp_peak = float(np.abs(reg_c).max().compute())

    return {
        "index": pdo_idx_c,
        "reg": reg_c,
        "var_expl": var_expl_v,
        "amplitude_rms": amp_rms,
        "amplitude_peak": amp_peak,
    }

# ERSSTv5: load global, regrid globally, then let compute_pdo_index cut the PDO box.
ersst = load_regular_data(dataset="ersst", variable="sst", time_range=time_range)
ersst_reg = regrid_horizontal(ersst, target=target_grid, method="bilinear")
results = {"ERSSTv5": _pdo_diagnostics(ersst_reg)}

# CMIP6 historical tos: load global, regrid globally, then compute PDO.
for model in models:
    da = load_regular_data(
        dataset="cmip6_historical",
        variable="tos",
        model=model,
        time_range=time_range,
    )
    da_reg = regrid_horizontal(da, target=target_grid, method="bilinear")
    results[model] = _pdo_diagnostics(da_reg)

# Pattern correlations against ERSSTv5 regression pattern on the common ocean mask.
regs = {name: results[name]["reg"] for name in all_names}
mask = build_comparison_mask(
    fields=[regs[name] for name in all_names],
    ocean_only=True,
)

ersst_masked = regs["ERSSTv5"].where(mask)
pattern_corr = {}
for model in models:
    model_masked = regs[model].where(mask)
    corr_res = analyze_correlation(
        var1=model_masked,
        var2=ersst_masked,
        dim=["lat", "lon"],
        method="pearson",
        weighting="area",
    )
    pattern_corr[model] = _extract_scalar(corr_res["correlation"])

# Lag-12-month autocorrelation of each standardized monthly index.
autocorr_lag12 = {}
for name in all_names:
    idx = results[name]["index"]
    ac_res = analyze_correlation(
        var1=idx,
        var2=idx,
        dim="time",
        method="pearson",
        lag=12,
    )
    autocorr_lag12[name] = _extract_scalar(ac_res["correlation"])

# Published CPC/NCEI PDO series and its correlation with the ERSSTv5-derived index.
pdo_pub = load_regular_data(
    dataset="cpc_indices",
    product="pdo",
    time_range=time_range,
)
pdo_corr_res = analyze_correlation(
    var1=results["ERSSTv5"]["index"],
    var2=pdo_pub,
    dim="time",
    method="pearson",
)
pdo_pub_corr = _extract_scalar(pdo_corr_res["correlation"])

# Standardize published PDO for plotting on the same axis.
pdo_pub_std = (
    (pdo_pub - pdo_pub.mean(dim="time")) / pdo_pub.std(dim="time")
).compute()

# 9-year (108-month) centred running mean of the ERSSTv5 standardized index.
smoothed = temporal_aggregate(
    data=results["ERSSTv5"]["index"],
    freq="rolling",
    window=108,
    center=True,
    op="mean",
)
smoothed_c = smoothed.compute()

# Zero crossings of the smoothed ERSSTv5 index, then merge crossings less than 3 years apart.
times = pd.to_datetime(smoothed_c.time.values)
vals = smoothed_c.values.astype(float)
crossing_years = []
for i in range(1, len(vals)):
    if np.isfinite(vals[i]) and np.isfinite(vals[i - 1]) and vals[i - 1] * vals[i] < 0:
        crossing_years.append(times[i].year)

transition_years = []
for yr in crossing_years:
    if not transition_years or yr - transition_years[-1] >= 3:
        transition_years.append(yr)

documented_shift_map = {
    1976: "1976/77", 1977: "1976/77",
    1988: "1988/89", 1989: "1988/89",
    1998: "1998/99", 1999: "1998/99",
}
documented_matches = {
    yr: documented_shift_map[yr]
    for yr in transition_years
    if yr in documented_shift_map
}

best_model = max(models, key=lambda m: pattern_corr[m])

summary = pd.DataFrame([
    {
        "source": name,
        "variance_explained_pct": results[name]["var_expl"],
        "amplitude_rms_degC_per_SD": results[name]["amplitude_rms"],
        "amplitude_peak_degC_per_SD": results[name]["amplitude_peak"],
        "lag12_autocorr": autocorr_lag12[name],
        "pattern_corr_with_ERSSTv5": 1.0 if name == "ERSSTv5" else pattern_corr[name],
    }
    for name in all_names
])

print(summary.round(3).to_string(index=False))
print("Correlation ERSSTv5-derived index vs published PDO:", round(pdo_pub_corr, 3))
print("Transition years:", transition_years)
print("Documented shift matches:", documented_matches)
print("Best pattern match:", best_model, round(pattern_corr[best_model], 3))
