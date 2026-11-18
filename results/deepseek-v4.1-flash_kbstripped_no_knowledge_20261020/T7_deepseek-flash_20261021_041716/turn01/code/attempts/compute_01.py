# ------------------------------------------------------------------
# Nitrate : phosphate stoichiometry of intermediate water (200-1000 m)
# CMIP6 historical (3 models, 1971-2000 climatology) vs WOA23 annual
# ------------------------------------------------------------------
import numpy as np
import pandas as pd
import xarray as xr

models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
order  = ["WOA23"] + models
DEPTH  = [200.0, 1000.0]
TR     = ["1971-01", "2000-12"]
BASINS = ["Pacific", "Atlantic", "Indian"]
RHO    = 1.025          # kg/L  == 1025 kg m-3

# ---------------- 1. load & prepare each source ----------------
src = {}
for m in models:
    no3_c = temporal_aggregate(
        load_regular_data(dataset="cmip6_historical", variable="no3", model=m,
                          time_range=TR, depth_range=DEPTH),
        freq="mean", op="mean")                       # day-weighted 1971-2000 climatology
    po4_c = temporal_aggregate(
        load_regular_data(dataset="cmip6_historical", variable="po4", model=m,
                          time_range=TR, depth_range=DEPTH),
        freq="mean", op="mean")
    src[m] = dict(
        no3=transform_unit(no3_c, from_unit=no3_c.attrs.get("units", "mol m-3"),
                           to_unit="umol/kg", density=RHO).assign_attrs(units="umol kg-1"),
        po4=transform_unit(po4_c, from_unit=po4_c.attrs.get("units", "mol m-3"),
                           to_unit="umol/kg", density=RHO).assign_attrs(units="umol kg-1"),
    )

src["WOA23"] = dict(
    no3=load_climatology_data(dataset="woa23", variable="nitrate",
                              time_slice="annual", depth_range=DEPTH),
    po4=load_climatology_data(dataset="woa23", variable="phosphate",
                              time_slice="annual", depth_range=DEPTH),
)

# ---------------- helpers ----------------
def _fit_terms(res):
    """Return (slope, intercept) from fit_multivariate_relation output."""
    c = res["coefficients"]
    if isinstance(c, dict):
        inter = slope = None
        for k, v in c.items():
            val = float(np.asarray(v, dtype=float).ravel()[0])
            kl = str(k).lower()
            if "intercept" in kl or kl in ("const", "b", "constant"):
                inter = val
            else:
                slope = val
        if slope is None or inter is None:
            raise ValueError(f"unparsed coefficients: {c}")
        return slope, inter
    arr = np.asarray(c, dtype=float).ravel()
    if arr.size != 2:
        raise ValueError(f"unparsed coefficients: {c}")
    # an N:P-vs-PO4 slope is O(10); the intercept is O(1) -> disambiguate order
    if abs(arr[0]) > 5.0 and abs(arr[1]) < 5.0:
        return float(arr[0]), float(arr[1])
    return float(arr[1]), float(arr[0])

def _subsample(x, y, n=20000, seed=0):
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if x.size > n:
        idx = np.random.default_rng(seed).choice(x.size, n, replace=False)
        x, y = x[idx], y[idx]
    return x, y

# ---------------- 2. per-source analysis ----------------
results, stats, maps = {}, {}, {}

for name in order:
    no3 = src[name]["no3"]
    po4 = src[name]["po4"]

    valid = no3.notnull() & po4.notnull() & (po4 > 0)   # ratio needs PO4 > 0
    no3v  = no3.where(valid)
    po4v  = po4.where(valid)

    npf = compute_ratio(no3v, po4v)
    npf = npf.assign_attrs(**{**no3v.attrs, "units": "1"})   # keep grid provenance

    # ---- headline: volume-weighted mean of the cell-by-cell ratio
    headline = float(compute_basic_stat(npf, op="volume_weighted_mean",
                                        dim=["depth", "lat", "lon"]).compute())

    # ---- ratio of the volume-weighted means
    mean_no3 = float(compute_basic_stat(no3v, op="volume_weighted_mean",
                                        dim=["depth", "lat", "lon"]).compute())
    mean_po4 = float(compute_basic_stat(po4v, op="volume_weighted_mean",
                                        dim=["depth", "lat", "lon"]).compute())
    ratio_of_means = mean_no3 / mean_po4

    # ---- unweighted OLS NO3 = a*PO4 + b
    fit = fit_multivariate_relation(target=no3v, predictors=[po4v], method="linear")
    a, b = _fit_terms(fit)
    r2   = float(fit["metrics"]["r2"])

    # ---- volume fractions
    # |N:P - 16| < 2  ==  14 < N:P < 18 (strict)
    band = compute_exceedance_fraction(npf, threshold=2.0, comparison="abs_less",
                                       reference=16.0, weighting="volume")
    low  = compute_exceedance_fraction(npf, threshold=12.0, comparison="less",
                                       weighting="volume")
    pct_14_18 = float(band["fraction_percent"])
    pct_below12_global = float(low["fraction_percent"])

    # ---- basins (marginal seas folded in; Southern/Arctic are separate basins)
    basin_mean = {}
    for bname in BASINS:
        sub = select_basin(npf, basin=bname, include_marginal_seas=True)
        basin_mean[bname] = float(compute_basic_stat(
            sub, op="volume_weighted_mean", dim=["depth", "lat", "lon"]).compute())

    pct_below12 = {}
    for bname in ["Pacific", "Indian"]:
        sub = select_basin(npf, basin=bname, include_marginal_seas=True)
        pct_below12[bname] = float(compute_exceedance_fraction(
            sub, threshold=12.0, comparison="less", weighting="volume")["fraction_percent"])

    results[name] = dict(
        mean_np_ratio=headline,
        ratio_of_means=ratio_of_means,
        mean_no3=mean_no3, mean_po4=mean_po4,
        slope=a, intercept=b, r2=r2,
        pct_14_18=pct_14_18,
        pct_below12_global=pct_below12_global,
        basin_mean=basin_mean,
        pct_below12_basin=pct_below12,
    )

    # ---- scatter data (subsample) and map field (200-1000 m column mean)
    xv, yv = _subsample(po4v.values.ravel(), no3v.values.ravel())
    colmean = compute_basic_stat(npf, op="volume_weighted_mean", dim="depth")
    colmean = colmean.assign_attrs(**{**npf.attrs, "units": "1"})
    anom = (colmean - 16.0).assign_attrs(**{**npf.attrs, "units": "1"})
    stats[name] = dict(po4=xv, no3=yv, np_field=npf)

    if name != "WOA23":
        anom = regrid_horizontal(anom, target="1deg", method="bilinear")  # display only
    maps[name] = anom.compute()

# ---------------- 3. report ----------------
rows = []
for name in order:
    r = results[name]
    rows.append(dict(
        source=name,
        NP_cellwise=round(r["mean_np_ratio"], 2),
        NP_ratio_of_means=round(r["ratio_of_means"], 2),
        mean_NO3=round(r["mean_no3"], 3),
        mean_PO4=round(r["mean_po4"], 3),
        slope_a=round(r["slope"], 2),
        intercept_b=round(r["intercept"], 2),
        R2=round(r["r2"], 3),
        pct_14_18=round(r["pct_14_18"], 1),
        pct_below12=round(r["pct_below12_global"], 1),
        Pacific=round(r["basin_mean"]["Pacific"], 2),
        Atlantic=round(r["basin_mean"]["Atlantic"], 2),
        Indian=round(r["basin_mean"]["Indian"], 2),
        below12_Pacific=round(r["pct_below12_basin"]["Pacific"], 1),
        below12_Indian=round(r["pct_below12_basin"]["Indian"], 1),
    ))
summary = pd.DataFrame(rows).set_index("source")
print("Nitrate : phosphate stoichiometry, 200-1000 m (umol kg-1 basis)")
print(summary.to_string())
print()
for name in order:
    r = results[name]
    tag = "below" if r["mean_np_ratio"] < 16 else ("above" if r["mean_np_ratio"] > 16 else "at")
    print(f"{name:14s} headline N:P = {r['mean_np_ratio']:.2f} -> {tag} Redfield 16")
