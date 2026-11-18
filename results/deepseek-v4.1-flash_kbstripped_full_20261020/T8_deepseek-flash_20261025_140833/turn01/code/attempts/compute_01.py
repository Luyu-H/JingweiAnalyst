# =============================================================================
# Ocean carbon system of three CMIP6 models vs GLODAPv2.2016b
# DIC (dissic / TCO2), total alkalinity (talk / TAlk), pH (ph / pHtsinsitutp)
# =============================================================================
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
dsid = "cmip6_historical"
VAR2002 = ["dissic", "talk", "ph", "thetao", "so"]
VARHIST = ["dissic", "thetao", "so"]


def scalar(da):
    """Materialise a (possibly dask-backed) reduction into a plain float."""
    return float(np.asarray(da.compute()).squeeze())


# -----------------------------------------------------------------------------
# 1. Load the sources (global fields only: the conservative regridding needs
#    the full source grid)
# -----------------------------------------------------------------------------
raw = {m: {} for m in models}
for m in models:
    for v in VAR2002:
        raw[m][v] = load_regular_data(dataset=dsid, variable=v, model=m,
                                      time_range=["2002-01", "2002-12"])
    for v in VARHIST:
        raw[m]["hist_" + v] = load_regular_data(dataset=dsid, variable=v, model=m,
                                                time_range=["1971-01", "2000-12"])

# GLODAP mapped climatology (no time axis)
g_tco2 = load_climatology_data(dataset="glodap", variable="TCO2")
g_talk = load_climatology_data(dataset="glodap", variable="TAlk")
g_ph = load_climatology_data(dataset="glodap", variable="pHtsinsitutp")

# -----------------------------------------------------------------------------
# 2. Time aggregation (2002 annual mean; 1971-2000 day-weighted long-term mean)
#    and unit conversion of the volume-basis CMIP6 carbon fields to umol/kg.
# -----------------------------------------------------------------------------
ann = {m: {} for m in models}
for m in models:
    for v in VAR2002:
        ann[m][v] = temporal_aggregate(raw[m][v], freq="mean", op="mean")
    for v in VARHIST:
        ann[m]["hist_" + v] = temporal_aggregate(raw[m]["hist_" + v], freq="mean", op="mean")

for m in models:
    # in-situ density anomaly from the matching 2002 thetao/so
    rho02 = compute_potential_density_anomaly(ann[m]["thetao"], ann[m]["so"], ref_pressure=None)
    ann[m]["rho2002"] = rho02
    ann[m]["dissic_umol"] = transform_unit(
        ann[m]["dissic"], from_unit=ann[m]["dissic"].attrs.get("units", "mol m-3"),
        to_unit="umol/kg", density=rho02)
    ann[m]["talk_umol"] = transform_unit(
        ann[m]["talk"], from_unit=ann[m]["talk"].attrs.get("units", "mol m-3"),
        to_unit="umol/kg", density=rho02)
    ann[m]["ph"] = ann[m]["ph"]                      # no unit conversion (total scale, in situ)

    # same for the 1971-2000 mean state (used only for the period change)
    rhoh = compute_potential_density_anomaly(ann[m]["hist_thetao"], ann[m]["hist_so"], ref_pressure=None)
    ann[m]["hist_rho"] = rhoh
    ann[m]["hist_dissic_umol"] = transform_unit(
        ann[m]["hist_dissic"], from_unit=ann[m]["hist_dissic"].attrs.get("units", "mol m-3"),
        to_unit="umol/kg", density=rhoh)
    ann[m]["hist_dissic_umol"] = ann[m]["hist_dissic_umol"].assign_attrs(
        density_assumption=ann[m]["dissic_umol"].attrs.get("density_assumption", "in-situ density field"))

# -----------------------------------------------------------------------------
# 3. Global 0-50 m means, each source over its own valid cells
# -----------------------------------------------------------------------------
DIC_0_50, TALK_0_50, PH_0_50 = {}, {}, {}
for m in models:
    DIC_0_50[m] = scalar(compute_basic_stat(ann[m]["dissic_umol"].sel(depth=slice(0, 50)),
                                            op="volume_weighted_mean"))
    TALK_0_50[m] = scalar(compute_basic_stat(ann[m]["talk_umol"].sel(depth=slice(0, 50)),
                                             op="volume_weighted_mean"))
    PH_0_50[m] = scalar(compute_basic_stat(ann[m]["ph"].sel(depth=slice(0, 50)),
                                           op="volume_weighted_mean"))

DIC_0_50["GLODAP"] = scalar(compute_basic_stat(g_tco2.sel(depth=slice(0, 50)), op="volume_weighted_mean"))
TALK_0_50["GLODAP"] = scalar(compute_basic_stat(g_talk.sel(depth=slice(0, 50)), op="volume_weighted_mean"))
PH_0_50["GLODAP"] = scalar(compute_basic_stat(g_ph.sel(depth=slice(0, 50)), op="volume_weighted_mean"))

DIC_0_50_bias = {m: DIC_0_50[m] - DIC_0_50["GLODAP"] for m in models}
TALK_0_50_bias = {m: TALK_0_50[m] - TALK_0_50["GLODAP"] for m in models}
PH_0_50_bias = {m: PH_0_50[m] - PH_0_50["GLODAP"] for m in models}

# -----------------------------------------------------------------------------
# 4. Map comparison: per-cell 0-50 m mean DIC, regridded conservatively onto the
#    GLODAP 1-degree grid, one shared mask, then bias / RMSE / spatial pattern
# -----------------------------------------------------------------------------
map_nat = {}
for m in models:
    map_nat[m] = compute_basic_stat(ann[m]["dissic_umol"].sel(depth=slice(0, 50)),
                                    op="volume_weighted_mean", dim="depth")
map_nat["GLODAP"] = compute_basic_stat(g_tco2.sel(depth=slice(0, 50)),
                                       op="volume_weighted_mean", dim="depth")

map_reg = {}
for m in models:
    map_reg[m] = regrid_horizontal(map_nat[m], target="glodap", method="conservative")
map_reg["GLODAP"] = map_nat["GLODAP"]

common_mask = build_comparison_mask([map_reg[m] for m in models] + [map_reg["GLODAP"]],
                                    ocean_only=False)
mask_n = int(common_mask.sum().compute())
mask_total = int(common_mask.size)

map_metrics = {}
bias_maps = {}
for m in models:
    a = map_reg[m].where(common_mask)
    b = map_reg["GLODAP"].where(common_mask)
    dm = compute_difference_metrics(data=a, reference=b, dim=["lat", "lon"])
    cor = analyze_correlation(a, b, dim=["lat", "lon"], weighting="area")
    map_metrics[m] = {
        "bias": float(dm["bias"]),
        "rmse": float(dm["rmse"]),
        "centered_rmse": float(dm["centered_rmse"]),
        "mae": float(dm["mae"]),
        "std_ratio": float(dm["std_ratio"]),
        "n_valid": int(dm["n_valid"]),
        "correlation": float(cor["correlation"].compute()),
        "pvalue": float(cor["pvalue"].compute()),
    }
    bias_maps[m] = (map_reg[m].where(common_mask) - map_reg["GLODAP"].where(common_mask)).compute()

# -----------------------------------------------------------------------------
# 5. Vertical DIC bias profile on the GLODAP 1-degree grid
# -----------------------------------------------------------------------------
profile_depths = [0, 100, 500, 1000, 2000]
prof = {}
for m in models:
    p = regrid_vertical(ann[m]["dissic_umol"], target=profile_depths, method="linear",
                        fill_surface=True)
    prof[m] = regrid_horizontal(p, target="glodap", method="conservative")
prof["GLODAP"] = regrid_vertical(g_tco2, target=profile_depths, method="linear", fill_surface=True)

vert_bias = {m: {} for m in models}
vert_ncells = {m: {} for m in models}
for m in models:
    for d in profile_depths:
        a = prof[m].sel(depth=d, method="nearest")
        b = prof["GLODAP"].sel(depth=d, method="nearest")
        valid = a.notnull() & b.notnull()
        vert_ncells[m][d] = int(valid.sum().compute())
        dm = compute_difference_metrics(data=a, reference=b, dim=["lat", "lon"])
        vert_bias[m][d] = float(dm["bias"])

# -----------------------------------------------------------------------------
# 6. DIC period change (native grid, own cell measures) and its ratio to the
#    2002 0-50 m DIC bias
# -----------------------------------------------------------------------------
DIC_2002_glob = {}
DIC_hist_glob = {}
for m in models:
    DIC_2002_glob[m] = scalar(compute_basic_stat(ann[m]["dissic_umol"].sel(depth=slice(0, 50)),
                                                 op="volume_weighted_mean"))
    DIC_hist_glob[m] = scalar(compute_basic_stat(ann[m]["hist_dissic_umol"].sel(depth=slice(0, 50)),
                                                 op="volume_weighted_mean"))

DIC_change = {m: DIC_2002_glob[m] - DIC_hist_glob[m] for m in models}
DIC_change_ratio = {m: DIC_change[m] / DIC_0_50_bias[m] for m in models}
DIC_change_ratio_via_tool = {
    m: scalar(compute_ratio(xr.DataArray(DIC_change[m], attrs={"units": "umol/kg"}),
                            xr.DataArray(DIC_0_50_bias[m], attrs={"units": "umol/kg"})))
    for m in models
}

# -----------------------------------------------------------------------------
# 7. Synthesis: rank models by overall carbon-system error, check sign agreement
# -----------------------------------------------------------------------------
def _zrank(values, higher_is_worse=True):
    order = np.argsort(values)
    return order

dic_bias_abs = {m: abs(DIC_0_50_bias[m]) for m in models}
talk_bias_abs = {m: abs(TALK_0_50_bias[m]) for m in models}
ph_bias_abs = {m: abs(PH_0_50_bias[m]) for m in models}
rmse_rank = {m: map_metrics[m]["rmse"] for m in models}
corr_rank = {m: map_metrics[m]["correlation"] for m in models}

# integer rank score (1 = best) summed over the three independent error measures
def _rank_score(d):
    order = sorted(d, key=lambda k: d[k])
    return {k: i + 1 for i, k in enumerate(order)}

score = {m: 0 for m in models}
for d in (dic_bias_abs, talk_bias_abs, ph_bias_abs, rmse_rank, {m: -corr_rank[m] for m in models}):
    for k, v in _rank_score(d).items():
        score[k] += v
best_model = min(score, key=score.get)

sign_agreement = {
    "DIC":  (len({np.sign(round(DIC_0_50_bias[m], 6)) for m in models}) == 1,
             {m: float(np.sign(DIC_0_50_bias[m])) for m in models}),
    "TAlk": (len({np.sign(round(TALK_0_50_bias[m], 6)) for m in models}) == 1,
             {m: float(np.sign(TALK_0_50_bias[m])) for m in models}),
    "pH":   (len({np.sign(round(PH_0_50_bias[m], 6)) for m in models}) == 1,
             {m: float(np.sign(PH_0_50_bias[m])) for m in models}),
}

# -----------------------------------------------------------------------------
# Reporting summary (unrounded intermediates kept above)
# -----------------------------------------------------------------------------
summary = {}
summary["DIC_0_50"] = {k: round(v, 1) for k, v in DIC_0_50.items()}
summary["TAlk_0_50"] = {k: round(v, 1) for k, v in TALK_0_50.items()}
summary["pH_0_50"] = {k: round(v, 3) for k, v in PH_0_50.items()}
summary["DIC_0_50_bias"] = {k: round(v, 1) for k, v in DIC_0_50_bias.items()}
summary["TAlk_0_50_bias"] = {k: round(v, 1) for k, v in TALK_0_50_bias.items()}
summary["pH_0_50_bias"] = {k: round(v, 3) for k, v in PH_0_50_bias.items()}
summary["map_bias"] = {k: round(v["bias"], 1) for k, v in map_metrics.items()}
summary["map_rmse"] = {k: round(v["rmse"], 1) for k, v in map_metrics.items()}
summary["map_corr"] = {k: round(v["correlation"], 3) for k, v in map_metrics.items()}
summary["mask_cells"] = {"kept": mask_n, "total": mask_total}
summary["vertical_bias"] = {m: {d: round(vert_bias[m][d], 1) for d in profile_depths} for m in models}
summary["DIC_2002"] = {m: round(DIC_2002_glob[m], 1) for m in models}
summary["DIC_1971_2000"] = {m: round(DIC_hist_glob[m], 1) for m in models}
summary["DIC_change"] = {m: round(DIC_change[m], 1) for m in models}
summary["DIC_change_ratio"] = {m: round(DIC_change_ratio[m], 2) for m in models}
summary["best_model"] = best_model
summary["rank_score"] = score
summary["sign_agreement"] = {k: v[0] for k, v in sign_agreement.items()}

print("=== 0-50 m global means (volume-weighted, each source's own cells) ===")
print("DIC  (umol/kg):", {k: round(v, 1) for k, v in DIC_0_50.items()})
print("TAlk (umol/kg):", {k: round(v, 1) for k, v in TALK_0_50.items()})
print("pH   (total)  :", {k: round(v, 3) for k, v in PH_0_50.items()})
print("DIC  bias:", {k: round(v, 1) for k, v in DIC_0_50_bias.items()})
print("TAlk bias:", {k: round(v, 1) for k, v in TALK_0_50_bias.items()})
print("pH   bias:", {k: round(v, 3) for k, v in PH_0_50_bias.items()})
print("=== 0-50 m DIC map comparison on GLODAP 1-deg grid ===")
print("shared mask cells kept:", mask_n, "of", mask_total)
for m in models:
    print(m, "bias %.1f  rmse %.1f  r %.3f  n=%d" % (map_metrics[m]["bias"],
          map_metrics[m]["rmse"], map_metrics[m]["correlation"], map_metrics[m]["n_valid"]))
print("=== vertical DIC bias profile (umol/kg) ===")
for m in models:
    print(m, {d: round(vert_bias[m][d], 1) for d in profile_depths})
print("=== DIC period change 2002 - 1971/2000 ===")
for m in models:
    print(m, "2002 %.1f  hist %.1f  change %.1f  bias %.1f  ratio %.2f" % (
        DIC_2002_glob[m], DIC_hist_glob[m], DIC_change[m], DIC_0_50_bias[m], DIC_change_ratio[m]))
print("=== synthesis ===")
print("rank scores:", score, "-> best:", best_model)
print("sign agreement:", {k: v[0] for k, v in sign_agreement.items()}, sign_agreement)
