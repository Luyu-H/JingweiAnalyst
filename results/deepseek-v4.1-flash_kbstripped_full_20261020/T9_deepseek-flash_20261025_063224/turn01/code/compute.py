# Computation
model_names = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
time_range = ["1950-01", "2014-12"]


def _month_start(da):
    t = pd.to_datetime(da.time.values).to_period("M").to_timestamp()
    return da.assign_coords(time=t)


def _scalar_float(x):
    if hasattr(x, "compute"):
        x = x.compute()
    return float(np.asarray(x).squeeze())


def _pdo_from_global(field, source_name):
    idx, eof, ev, reg = compute_pdo_index(
        field, depth_level=None, standardize=True
    )
    idx = _month_start(idx.squeeze(drop=True).compute())
    eof = eof.compute()
    reg = reg.compute()

    ev_val = _scalar_float(ev)

    # Area-weighted RMS amplitude of the regression pattern.
    reg_sq = (reg**2).assign_attrs(reg.attrs).assign_attrs(units="degC2")
    rms = compute_basic_stat(
        reg_sq, op="area_weighted_mean", dim=["lat", "lon"]
    )
    amp_val = float(np.sqrt(_scalar_float(rms)))

    # Lag-12-month autocorrelation of the standardized monthly index.
    lagged = idx.shift(time=12)
    a = np.asarray(idx.values, dtype=float).ravel()
    b = np.asarray(lagged.values, dtype=float).ravel()
    mask = np.isfinite(a) & np.isfinite(b)
    lag12 = float(np.corrcoef(a[mask], b[mask])[0, 1])

    return {
        "source": source_name,
        "index": idx,
        "eof": eof,
        "explained_variance": ev_val,
        "reg_pattern": reg,
        "amplitude": amp_val,
        "lag12": lag12,
    }


# ERSSTv5: load global SST, regrid to 1 degree, compute PDO.
ersst = load_regular_data(
    dataset="ersst", variable="sst", time_range=time_range
)
ersst_1deg = regrid_horizontal(ersst, target="1deg", method="bilinear")

results = {}
results["ERSSTv5"] = _pdo_from_global(ersst_1deg, "ERSSTv5")

# CMIP6 historical tos: load global field, regrid to 1 degree, compute PDO.
for model in model_names:
    tos = load_regular_data(
        dataset="cmip6_historical",
        variable="tos",
        model=model,
        time_range=time_range,
    )
    tos_1deg = regrid_horizontal(tos, target="1deg", method="bilinear")
    results[model] = _pdo_from_global(tos_1deg, model)

# Published PDO index.
pub_pdo = load_regular_data(
    dataset="cpc_indices", product="pdo", time_range=time_range
)
pub_pdo = _month_start(pub_pdo.squeeze(drop=True).compute())

# Correlation of ERSSTv5-derived index with published PDO, as-is.
pub_corr_dict = analyze_correlation(
    var1=results["ERSSTv5"]["index"],
    var2=pub_pdo,
    dim="time",
    detrend=False,
)
pub_corr = _scalar_float(pub_corr_dict["correlation"])
pub_corr_pvalue = _scalar_float(pub_corr_dict["pvalue"])
pub_corr_pvalue_unadjusted = _scalar_float(
    pub_corr_dict["pvalue_unadjusted"]
)
pub_corr_n_eff = _scalar_float(pub_corr_dict["n_effective"])
pub_corr_n = _scalar_float(pub_corr_dict["n_samples"])

# Spatial pattern correlations against ERSSTv5.
pattern_corr = {}
for model in model_names:
    cdict = analyze_correlation(
        var1=results[model]["reg_pattern"],
        var2=results["ERSSTv5"]["reg_pattern"],
        dim=["lat", "lon"],
        weighting="area",
    )
    pattern_corr[model] = _scalar_float(cdict["correlation"])
pattern_corr["ERSSTv5"] = 1.0

# 9-year (108-month) centred running mean of the ERSSTv5 index.
ersst_smoothed = temporal_aggregate(
    results["ERSSTv5"]["index"],
    freq="rolling",
    window=108,
    center=True,
    op="mean",
).compute()


def _sign_crossings(smoothed, merge_gap_years=3):
    s = smoothed.squeeze(drop=True).sortby("time")
    times = pd.to_datetime(s.time.values)
    vals = np.asarray(s.values, dtype=float).ravel()

    crossings = []
    prev_sign = None
    for t, v in zip(times, vals):
        if not np.isfinite(v) or v == 0:
            continue
        sign = 1 if v > 0 else -1
        if prev_sign is not None and sign != prev_sign:
            crossings.append(int(t.year))
        prev_sign = sign

    merged = []
    for y in crossings:
        if not merged or (y - merged[-1]) >= merge_gap_years:
            merged.append(y)
    return merged, crossings


transition_years, all_crossing_years = _sign_crossings(
    ersst_smoothed, merge_gap_years=3
)

shift_labels = {}
for y in transition_years:
    if 1975 <= y <= 1977:
        shift_labels[y] = "mid-1970s shift"
    elif 1987 <= y <= 1989:
        shift_labels[y] = "late-1980s shift"
    elif 1997 <= y <= 1999:
        shift_labels[y] = "late-1990s shift"
    else:
        shift_labels[y] = "other / not one of the listed shifts"

summary_rows = []
for name in ["ERSSTv5"] + model_names:
    summary_rows.append(
        {
            "Source": name,
            "Explained variance (%)": round(
                results[name]["explained_variance"], 1
            ),
            "Amplitude (degC per std)": round(
                results[name]["amplitude"], 3
            ),
            "Lag-12 autocorrelation": round(results[name]["lag12"], 3),
            "Pattern corr. with ERSSTv5": round(pattern_corr[name], 3),
        }
    )
summary_table = pd.DataFrame(summary_rows)

best_model = max(model_names, key=lambda m: pattern_corr[m])
amplitude_ratio = {
    m: results[m]["amplitude"] / results["ERSSTv5"]["amplitude"]
    for m in model_names
}
lag12_difference = {
    m: results[m]["lag12"] - results["ERSSTv5"]["lag12"]
    for m in model_names
}
