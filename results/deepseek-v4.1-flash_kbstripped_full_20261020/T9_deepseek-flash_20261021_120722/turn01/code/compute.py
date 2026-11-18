import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt


def _scalar(x):
    """Materialise a scalar-like xarray/dask/numpy object."""
    if hasattr(x, "compute"):
        x = x.compute()
    size = int(getattr(x, "size", 1))
    if size == 1 and hasattr(x, "item"):
        return float(x.item())
    if size == 1:
        return float(np.asarray(x).reshape(-1)[0])
    return float(x)


def regrid_tos(model):
    da = load_regular_data(
        dataset="cmip6_historical",
        variable="tos",
        model=model,
        time_range=["1950-01", "2014-12"],
    )
    # CMIP6 tos is usually in K; convert to degC when needed.
    units = str(da.attrs.get("units", "")).lower()
    if units in ("k", "kelvin", "degrees_kelvin", "degree_kelvin"):
        da = transform_unit(data=da, from_unit="K", to_unit="degC")
    da = regrid_horizontal(data=da, target="1deg", method="bilinear")
    if "depth" in da.dims:
        da = da.isel(depth=0, drop=True)
    return da


# Load and regrid all SST sources onto the common 1-degree grid.
ersst_reg = regrid_horizontal(
    data=load_regular_data(
        dataset="ersst",
        variable="sst",
        time_range=["1950-01", "2014-12"],
    ),
    target="1deg",
    method="bilinear",
)

fields = {"ERSSTv5": ersst_reg}
for model in ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]:
    fields[model] = regrid_tos(model)

# Compute PDO modes for each source.
sources = {}
for name, field in fields.items():
    anom = compute_anomaly_vs_own_climatology(
        data=field,
        freq="month",
        base_period=["1950-01", "2014-12"],
        method="subtract",
    )
    idx, eof, ev, pat = compute_pdo_index(
        data=anom,
        depth_level=None,
        standardize=True,
    )
    pat = pat.squeeze(drop=True)
    idx = idx.sel(time=slice("1950-01", "2014-12"))
    idx_std = (idx - idx.mean("time")) / idx.std("time")

    sources[name] = {
        "index_raw": idx,
        "index": idx_std,
        "eof": eof.squeeze(drop=True) if hasattr(eof, "squeeze") else eof,
        "var_exp": ev,
        "regpat": pat,
    }

# Explained variance, area-weighted RMS amplitude, and lag-12 autocorrelation.
for name in sources:
    pat = sources[name]["regpat"]
    sources[name]["var_exp_pct"] = _scalar(sources[name]["var_exp"])

    pat_sq = (pat ** 2).assign_attrs(pat.attrs)
    rms2 = compute_basic_stat(
        data=pat_sq,
        op="area_weighted_mean",
        dim=["lat", "lon"],
    )
    sources[name]["amplitude"] = float(np.sqrt(_scalar(rms2)))

    ac = analyze_correlation(
        var1=sources[name]["index"],
        var2=sources[name]["index"],
        dim="time",
        method="pearson",
        lag=12,
    )
    sources[name]["lag12"] = _scalar(ac["correlation"])
    sources[name]["lag12_pvalue"] = _scalar(ac["pvalue"])

# Spatial pattern correlations with ERSSTv5.
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
for name in models:
    cc = analyze_correlation(
        var1=sources[name]["regpat"],
        var2=sources["ERSSTv5"]["regpat"],
        dim=["lat", "lon"],
        method="pearson",
        weighting="area",
    )
    sources[name]["pattern_r"] = _scalar(cc["correlation"])
    sources[name]["pattern_pvalue"] = _scalar(cc["pvalue"])
    sources[name]["pattern_n"] = _scalar(cc["n_samples"])

# Published PDO index. The 'pdo' product is a basin-wide pattern index on its
# own domain, so it takes no region/box argument: only product is given.
pdo_pub = load_regular_data(
    dataset="cpc_indices",
    product="pdo",
    time_range=["1950-01", "2014-12"],
)
pdo_pub = pdo_pub.squeeze(drop=True).sel(time=slice("1950-01", "2014-12"))
pdo_pub_base = pdo_pub.attrs.get(
    "climatology_period",
    pdo_pub.attrs.get(
        "base_period",
        pdo_pub.attrs.get("baseline", "see product metadata"),
    ),
)
pdo_pub_plot = (pdo_pub - pdo_pub.mean("time")) / pdo_pub.std("time")

# Correlation of the ERSSTv5-derived standardized index with the published PDO.
ersst_idx_aligned, pdo_pub_aligned = xr.align(
    sources["ERSSTv5"]["index"],
    pdo_pub,
    join="inner",
)
pub_corr = analyze_correlation(
    var1=ersst_idx_aligned,
    var2=pdo_pub_aligned,
    dim="time",
    method="pearson",
)
pub_r = _scalar(pub_corr["correlation"])
pub_p = _scalar(pub_corr["pvalue"])
pub_neff = _scalar(pub_corr["n_effective"])
pub_n = _scalar(pub_corr["n_samples"])

# ERSSTv5 9-year (108-month) centred running mean and sign-change years.
ersst_smooth = temporal_aggregate(
    data=sources["ERSSTv5"]["index"],
    freq="rolling",
    op="mean",
    window=108,
    center=True,
).sel(time=slice("1950-01", "2014-12"))
ersst_smooth_c = ersst_smooth.compute()

sv = np.asarray(ersst_smooth_c.values, dtype=float)
tv = pd.to_datetime(ersst_smooth_c.time.values)
sg = np.sign(sv)

raw_crossing_years = []
for i in range(1, len(sv)):
    if (
        np.isfinite(sv[i])
        and np.isfinite(sv[i - 1])
        and sg[i] != 0
        and sg[i - 1] != 0
        and sg[i] != sg[i - 1]
    ):
        raw_crossing_years.append(int(tv[i].year))

transition_years = []
for year in raw_crossing_years:
    if not transition_years or (year - transition_years[-1]) >= 3:
        transition_years.append(year)

# Compute final lazy objects for plotting / reuse.
for name in sources:
    sources[name]["index"] = sources[name]["index"].compute()
    sources[name]["regpat"] = sources[name]["regpat"].compute()
pdo_pub_plot = pdo_pub_plot.compute()

# Report results.
print("PDO assessment, 1950-2014; bilinear regridding to the common 1-degree grid")
print("(Bilinear is the prescribed comparison method; conservative regridding would be preferred when coarsening before an area-weighted step.)")
print()
print(f"{'source':12s} {'var_expl_%':>10s} {'amplitude_C_per_SD':>20s} {'lag12_autocorr':>16s}")
for name in ["ERSSTv5", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]:
    print(
        f"{name:12s} "
        f"{sources[name]['var_exp_pct']:10.1f} "
        f"{sources[name]['amplitude']:20.3f} "
        f"{sources[name]['lag12']:16.3f}"
    )

print("\nSpatial pattern correlations with ERSSTv5 (area-weighted Pearson over shared cells):")
for name in models:
    print(
        f"{name:12s} r={sources[name]['pattern_r']:.3f}  "
        f"p={sources[name]['pattern_pvalue']:.3g}  "
        f"n_shared={sources[name]['pattern_n']:.0f}"
    )

print("\nERSSTv5-derived standardized index vs published PDO:")
print(
    f"  Pearson r = {pub_r:.3f}; p = {pub_p:.3g}; "
    f"n = {pub_n:.0f}; n_effective = {pub_neff:.1f}"
)
print("  This demonstrates temporal agreement between two PDO estimates; it does not demonstrate equal variance or amplitude.")
print(f"  Published PDO base/climatology metadata: {pdo_pub_base}")

print("\nERSSTv5 9-year (108-month) centred running-mean sign changes:")
print(f"  raw crossing years: {raw_crossing_years}")
print(f"  merged transition years (<3 yr after last kept merged into last): {transition_years}")

documented_shifts = {
    1976: "1976/77 PDO regime shift",
    1977: "1976/77 PDO regime shift",
    1989: "late-1980s Pacific shift",
    1998: "1998/99 Pacific shift",
    1999: "1998/99 Pacific shift",
    2014: "2013/14-2014 North Pacific shift (near record end; may be unresolved)",
}
print("  Documented-shift comparison:")
for year in transition_years:
    if year in documented_shifts:
        print(f"    {year}: {documented_shifts[year]}")
    else:
        near = [d for d in documented_shifts if abs(year - d) <= 1]
        if near:
            d = near[0]
            print(f"    {year}: near {d} ({documented_shifts[d]})")
        else:
            print(f"    {year}: not one of the commonly cited shifts listed here")

best_model = max(models, key=lambda m: sources[m]["pattern_r"])
print(
    f"\nBest model by spatial pattern correlation: {best_model} "
    f"(r={sources[best_model]['pattern_r']:.3f})"
)

print("\nAmplitude and persistence comparison (no threshold imposed):")
print(
    f"  ERSSTv5: amplitude {sources['ERSSTv5']['amplitude']:.3f} °C/SD, "
    f"lag-12 {sources['ERSSTv5']['lag12']:.3f}"
)
for name in models:
    print(
        f"  {name}: amplitude {sources[name]['amplitude']:.3f} °C/SD, "
        f"lag-12 {sources[name]['lag12']:.3f}"
    )
