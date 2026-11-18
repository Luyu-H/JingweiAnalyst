# ------------------------------------------------------------------
# 1. WOA23 reference: annual objective analysis, 200-1000 m
# ------------------------------------------------------------------
woa_no3 = load_climatology_data(dataset="woa23", variable="nitrate",
                                time_slice="annual", depth_range=[200, 1000])
woa_po4 = load_climatology_data(dataset="woa23", variable="phosphate",
                                time_slice="annual", depth_range=[200, 1000])
woa_lev, woa_bnd = load_depth_levels(dataset="woa23", variable="nitrate",
                                     time_slice="annual", depth_range=[200, 1000])
woa_area = load_areas(dataset="woa23").compute()
basin_mask = load_basin_mask()          # WOA23 1-deg basin codes
woa_thick = np.abs(np.diff(woa_bnd, axis=1)).ravel()
woa_vol = (woa_area * xr.DataArray(woa_thick, dims="depth",
                                   coords={"depth": np.asarray(woa_lev)})).transpose("depth", "lat", "lon")

# N:P ratio field (mol mol-1). Re-attach grid provenance so the weighting tools
# can find the WOA23 cell areas / layer thicknesses.
woa_np = compute_ratio(woa_no3, woa_po4).assign_attrs(**{**woa_no3.attrs, "units": "1"})

# ------------------------------------------------------------------
# 2-5. CMIP6 historical no3 / po4, monthly 1971-2000 -> day-weighted
#      climatology -> umol kg-1 (rho = 1025) -> N:P
# ------------------------------------------------------------------
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
src = {}
src["WOA23"] = dict(no3=woa_no3, po4=woa_po4, np=woa_np, vol=woa_vol)

for m in models:
    n3 = load_regular_data(dataset="cmip6_historical", variable="no3", model=m,
                           time_range=["1971-01", "2000-12"], depth_range=[200, 1000])
    p4 = load_regular_data(dataset="cmip6_historical", variable="po4", model=m,
                           time_range=["1971-01", "2000-12"], depth_range=[200, 1000])
    lev, bnd = load_depth_levels(dataset="cmip6_historical", model=m,
                                 variable="no3", depth_range=[200, 1000])

    n3 = temporal_aggregate(n3, freq="mean")          # day-weighted 1971-2000 mean
    p4 = temporal_aggregate(p4, freq="mean")
    n3 = transform_unit(n3, from_unit="mol m-3", to_unit="umol/kg", density=1025.0)
    p4 = transform_unit(p4, from_unit="mol m-3", to_unit="umol/kg", density=1025.0)
    n3 = n3.assign_attrs(units="umol kg-1")
    p4 = p4.assign_attrs(units="umol kg-1")

    npf = compute_ratio(n3, p4).assign_attrs(**{**n3.attrs, "units": "1"})

    vol = load_volumes_for_cmip6(dataset="cmip6_historical", model=m,
                                 variable="no3", depth_range=[200, 1000])
    if "time" in vol.dims:
        vol = vol.mean("time")

    src[m] = dict(no3=n3, po4=p4, np=npf, vol=vol)

src_names = ["WOA23", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]


def _num(v):
    if isinstance(v, xr.DataArray):
        v = v.compute()
    return float(v)


# ------------------------------------------------------------------
# 7. Native-grid volume-weighted layer statistics, 200-1000 m
# ------------------------------------------------------------------
layer = {}
for s in src_names:
    d = src[s]
    mn  = compute_basic_stat(d["np"],  op="volume_weighted_mean").compute()
    m3  = compute_basic_stat(d["no3"], op="volume_weighted_mean").compute()
    mp  = compute_basic_stat(d["po4"], op="volume_weighted_mean").compute()
    layer[s] = dict(mean_np=_num(mn), mean_no3=_num(m3), mean_po4=_num(mp),
                    ratio_of_means=_num(m3) / _num(mp))

# ------------------------------------------------------------------
# 9. Volume percentage with 14 < N:P < 18 (native grid)
# ------------------------------------------------------------------
band = {}
for s in src_names:
    ex = compute_exceedance_fraction(src[s]["np"], threshold=2.0,
                                     comparison="abs_less", reference=16.0,
                                     weighting="volume")
    band[s] = dict(pct=_num(ex["fraction_percent"]),
                   m3=_num(ex["exceeding_amount"]))

# ------------------------------------------------------------------
# 10 & 11. Basin means of N:P and volume percentage below 12
# ------------------------------------------------------------------
basins3 = ["Pacific", "Atlantic", "Indian"]
basin_mean, below12, below12_basin = {}, {}, {}
for s in src_names:
    for b in basins3:
        nb = select_basin(src[s]["np"], basin=b, include_marginal_seas=True,
                          depth_dependent=False)
        basin_mean[(s, b)] = _num(compute_basic_stat(nb, op="volume_weighted_mean").compute())
    ex = compute_exceedance_fraction(src[s]["np"], threshold=12.0,
                                     comparison="less", weighting="volume")
    below12[s] = dict(pct=_num(ex["fraction_percent"]), m3=_num(ex["exceeding_amount"]))
    for b in ["Pacific", "Indian"]:
        nb = select_basin(src[s]["np"], basin=b, include_marginal_seas=True,
                          depth_dependent=False)
        exb = compute_exceedance_fraction(nb, threshold=12.0, comparison="less",
                                          weighting="volume")
        below12_basin[(s, b)] = dict(pct=_num(exb["fraction_percent"]),
                                     m3=_num(exb["exceeding_amount"]))

# ------------------------------------------------------------------
# 6 & 8. Common comparison sample (WOA23 grid) and per-source OLS
# ------------------------------------------------------------------
reg_no3 = {"WOA23": woa_no3}
reg_po4 = {"WOA23": woa_po4}
for m in models:
    reg_no3[m] = regrid_horizontal(src[m]["no3"], target="woa23", method="bilinear")
    reg_po4[m] = regrid_horizontal(src[m]["po4"], target="woa23", method="bilinear")

mask_fields = [woa_no3, woa_po4] + [reg_no3[m] for m in models] + [reg_po4[m] for m in models]
common = build_comparison_mask(mask_fields, ocean_only=True, reduce_other_dims="any")

ols = {}
scatter_data = {}
for s in src_names:
    n = reg_no3[s].where(common)
    p = reg_po4[s].where(common)
    ok = (n.notnull() & p.notnull()).compute()
    nv = n.where(ok).values.ravel()
    pv = p.where(ok).values.ravel()
    good = np.isfinite(nv) & np.isfinite(pv)
    nv, pv = nv[good], pv[good]

    # equal-weight OLS: NO3 = a * PO4 + b
    A = np.vstack([pv, np.ones_like(pv)]).T
    a_ols, b_ols = np.linalg.lstsq(A, nv, rcond=None)[0]
    fit_vals = a_ols * pv + b_ols
    ss_res = np.sum((nv - fit_vals) ** 2)
    ss_tot = np.sum((nv - nv.mean()) ** 2)
    r2_ols = 1.0 - ss_res / ss_tot

    # the same fit through the regression tool
    fit = fit_multivariate_relation(
        target=xr.DataArray(nv, dims="sample"),
        predictors=[xr.DataArray(pv, dims="sample")],
        method="linear", validation="random")

    ols[s] = dict(a=float(a_ols), b=float(b_ols), r2=float(r2_ols),
                  n=int(nv.size),
                  tool_r2=float(fit["metrics"]["r2"]))
    scatter_data[s] = (pv, nv)

# ------------------------------------------------------------------
# 13. 200-1000 m vertical volume-weighted mean N:P minus 16 (maps)
# ------------------------------------------------------------------
np_map = {}
for s in src_names:
    h = compute_basic_stat(src[s]["np"], op="volume_weighted_mean", dim="depth").compute()
    if s != "WOA23":
        h = regrid_horizontal(h, target="woa23", method="bilinear").compute()
    np_map[s] = (h - 16.0).assign_attrs(units="1")

# ------------------------------------------------------------------
# Report
# ------------------------------------------------------------------
print("=" * 108)
print("Intermediate-water (200-1000 m) N:P stoichiometry, 1971-2000 climatologies")
print("CMIP6 models converted mol m-3 -> umol kg-1 at constant rho = 1025 kg m-3")
print("=" * 108)
print(f"{'source':14s} {'meanNP':>7s} {'NO3/PO4':>8s} {'meanNO3':>8s} {'meanPO4':>8s} "
      f"{'a':>7s} {'b':>8s} {'R2':>7s} {'n':>7s} {'14-18%':>7s}")
for s in src_names:
    L, O, B = layer[s], ols[s], band[s]
    print(f"{s:14s} {L['mean_np']:7.2f} {L['ratio_of_means']:8.2f} {L['mean_no3']:8.2f} "
          f"{L['mean_po4']:8.2f} {O['a']:7.2f} {O['b']:8.2f} {O['r2']:7.3f} "
          f"{O['n']:7d} {B['pct']:7.1f}")

print("\nVolume with 14 < N:P < 18 (m3 behind the percentage):")
for s in src_names:
    print(f"  {s:14s} {band[s]['pct']:5.1f} %   ({band[s]['m3']:.3e} m3)")

print("\nHeadline mean N:P by basin (WOA23 basin mask, marginal seas folded in):")
print(f"{'source':14s} {'Pacific':>8s} {'Atlantic':>9s} {'Indian':>8s}")
for s in src_names:
    print(f"{s:14s} {basin_mean[(s,'Pacific')]:8.2f} {basin_mean[(s,'Atlantic')]:9.2f} "
          f"{basin_mean[(s,'Indian')]:8.2f}")

print("\nVolume with N:P < 12, global and by basin:")
for s in src_names:
    print(f"  {s:14s} global {below12[s]['pct']:5.1f} % ({below12[s]['m3']:.3e} m3)   "
          f"Pacific {below12_basin[(s,'Pacific')]['pct']:5.1f} % "
          f"({below12_basin[(s,'Pacific')]['m3']:.3e} m3)   "
          f"Indian {below12_basin[(s,'Indian')]['pct']:5.1f} % "
          f"({below12_basin[(s,'Indian')]['m3']:.3e} m3)")

print("\nCommon-mask cells used in the scatter/regression: " +
      ", ".join(f"{s}={ols[s]['n']}" for s in src_names))
