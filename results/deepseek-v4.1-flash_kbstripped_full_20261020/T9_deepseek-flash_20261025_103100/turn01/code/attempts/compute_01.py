import numpy as np
import pandas as pd
import xarray as xr

time_range = ["1950-01", "2014-12"]
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
plot_order = ["ERSSTv5"] + models


def as_float(da):
    """Materialise a scalar DataArray to a Python float."""
    return float(np.asarray(da.compute()).item())


def area_weighted_rms(pattern):
    """Area-weighted RMS of a regression pattern (degC per standard deviation)."""
    sq = pattern ** 2
    sq.attrs = dict(pattern.attrs)
    sq.attrs["units"] = "degC2"
    wm = compute_basic_stat(data=sq, op="area_weighted_mean", dim=["lat", "lon"])
    return float(np.sqrt(as_float(wm)))


def lag12_autocorr(index):
    """Lag-12 (1-year) autocorrelation of a monthly standardized index."""
    out = analyze_correlation(var1=index, var2=index, dim="time", lag=12)
    return as_float(out["correlation"])


# ------------------------------------------------------------------
# 1. ERSSTv5: load global, regrid bilinearly to 1 degree, anomalies,
#    leading North Pacific EOF / PDO mode.
# ------------------------------------------------------------------
ersst_raw = load_regular_data(dataset="ersst", variable="sst", time_range=time_range)
ersst_1deg = regrid_horizontal(ersst_raw, target="1deg", method="bilinear")
ersst_anom = compute_anomaly_vs_own_climatology(
    ersst_1deg, freq="month", base_period=time_range, method="subtract"
)
ersst_index, _, ersst_var, ersst_reg = compute_pdo_index(
    ersst_anom, depth_level=None, standardize=True
)
ersst_index = ersst_index.compute()
ersst_reg = ersst_reg.compute()
ersst_var = as_float(ersst_var)

results = {}
results["ERSSTv5"] = {
    "explained_variance": ersst_var,
    "amplitude": area_weighted_rms(ersst_reg),
    "autocorr12": lag12_autocorr(ersst_index),
    "index": ersst_index,
    "reg_pattern": ersst_reg,
}

# ------------------------------------------------------------------
# 2. CMIP6 `tos`: load global, regrid bilinearly to 1 degree,
#    anomalies, leading North Pacific EOF / PDO mode.
# ------------------------------------------------------------------
for model in models:
    raw = load_regular_data(
        dataset="cmip6_historical",
        variable="tos",
        model=model,
        time_range=time_range,
    )
    r1 = regrid_horizontal(raw, target="1deg", method="bilinear")
    anom = compute_anomaly_vs_own_climatology(
        r1, freq="month", base_period=time_range, method="subtract"
    )
    idx, _, var, reg = compute_pdo_index(anom, depth_level=None, standardize=True)
    idx = idx.compute()
    reg = reg.compute()
    var = as_float(var)
    results[model] = {
        "explained_variance": var,
        "amplitude": area_weighted_rms(reg),
        "autocorr12": lag12_autocorr(idx),
        "index": idx,
        "reg_pattern": reg,
    }

# ------------------------------------------------------------------
# 3. Published PDO index.
# ------------------------------------------------------------------
pdo_pub = load_regular_data(
    dataset="cpc_indices", product="pdo", time_range=time_range
).squeeze()
pdo_pub = pdo_pub.compute()
pdo_pub_std = ((pdo_pub - pdo_pub.mean(dim="time")) / pdo_pub.std(dim="time")).compute()

# Agreement between ERSSTv5-derived PDO index and published PDO.
ersst_idx_al, pdo_al = xr.align(results["ERSSTv5"]["index"], pdo_pub, join="inner")
corr_ersst_pub = analyze_correlation(
    var1=ersst_idx_al, var2=pdo_al, dim="time", method="pearson"
)
results["ERSSTv5"]["corr_with_published_pdo"] = as_float(corr_ersst_pub["correlation"])

# ------------------------------------------------------------------
# 4. Spatial pattern correlations: each model vs ERSSTv5.
#    Use only the cells the two regression patterns share.
# ------------------------------------------------------------------
for model in models:
    mask = build_comparison_mask(
        fields=[results["ERSSTv5"]["reg_pattern"], results[model]["reg_pattern"]],
        ocean_only=True,
    ).compute()
    n_shared = int(mask.sum().compute().item())
    a = results["ERSSTv5"]["reg_pattern"].where(mask)
    b = results[model]["reg_pattern"].where(mask)
    pc = analyze_correlation(
        var1=a,
        var2=b,
        dim=["lat", "lon"],
        method="pearson",
        weighting="none",
    )
    results[model]["pattern_corr_vs_ersst"] = as_float(pc["correlation"])
    results[model]["shared_cells"] = n_shared

# ------------------------------------------------------------------
# 5. Nine-year centred running mean of the ERSSTv5 index and
#    phase-transition years (crossings merged if < 3 years apart).
# ------------------------------------------------------------------
ersst_smoothed = temporal_aggregate(
    data=results["ERSSTv5"]["index"],
    freq="rolling",
    op="mean",
    window=108,
    center=True,
).compute()

sm = ersst_smoothed.dropna(dim="time")
sm_times = pd.to_datetime(sm["time"].values)
sm_vals = np.asarray(sm.values)
sm_signs = np.sign(sm_vals)

transition_years = []
last_counted_year = None
for i in range(1, len(sm)):
    if sm_signs[i] == 0 or sm_signs[i - 1] == 0:
        continue
    if sm_signs[i] != sm_signs[i - 1]:
        yr = int(sm_times[i].year)
        if last_counted_year is None or (yr - last_counted_year) >= 3:
            transition_years.append(yr)
            last_counted_year = yr
transition_years = sorted(set(transition_years))

# ------------------------------------------------------------------
# 6. Summary.
# ------------------------------------------------------------------
print("PDO assessment, 1950-2014, bilinear regridding to 1 degree")
print("Source        Var.expl.(%)  Amplitude(degC/s.d.)  AC12   Pattern r vs ERSSTv5  Shared cells")
for name in plot_order:
    r = results[name]
    pc_str = f"{r['pattern_corr_vs_ersst']:.3f}" if "pattern_corr_vs_ersst" in r else "-"
    n_str = str(r["shared_cells"]) if "shared_cells" in r else "-"
    print(
        f"{name:12s} {r['explained_variance']:11.1f} "
        f"{r['amplitude']:20.3f} {r['autocorr12']:8.3f} "
        f"{pc_str:>20s} {n_str:>12s}"
    )
print(
    "ERSSTv5-derived index vs published PDO correlation: "
    f"{results['ERSSTv5']['corr_with_published_pdo']:.3f}"
)
print(f"Detected PDO transition years: {transition_years}")
best_model = max(models, key=lambda m: results[m]["pattern_corr_vs_ersst"])
print(
    f"Best pattern match to ERSSTv5: {best_model} "
    f"(r = {results[best_model]['pattern_corr_vs_ersst']:.3f})"
)
