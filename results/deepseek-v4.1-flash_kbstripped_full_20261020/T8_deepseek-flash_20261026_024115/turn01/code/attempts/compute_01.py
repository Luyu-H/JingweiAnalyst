# Models and periods
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
period_2002 = ["2002-01", "2002-12"]
period_base = ["1971-01", "2000-12"]
profile_depths = [0, 100, 500, 1000, 2000]

# ---------------------------------------------------------------------------
# GLODAPv2.2016b mapped climatology
# ---------------------------------------------------------------------------
glodap_tco2 = load_climatology_data(dataset="glodap", variable="TCO2")
glodap_talk = load_climatology_data(dataset="glodap", variable="TAlk",
                                   depth_range=[0, 50])
glodap_ph = load_climatology_data(dataset="glodap", variable="pHtsinsitutp",
                                  depth_range=[0, 50])

glodap_tco2_0_50 = glodap_tco2.sel(depth=slice(0, 50))

# Volume-weighted 0-50 m means on the GLODAP grid
glodap_dic_mean = compute_basic_stat(glodap_tco2_0_50,
                                     op="volume_weighted_mean",
                                     dim=None).compute()
glodap_talk_mean = compute_basic_stat(glodap_talk,
                                      op="volume_weighted_mean",
                                      dim=None).compute()
glodap_ph_mean = compute_basic_stat(glodap_ph,
                                    op="volume_weighted_mean",
                                    dim=None).compute()

# 0-50 m DIC map (each cell's mean over its own 0-50 m levels)
glodap_dic_map = compute_basic_stat(glodap_tco2_0_50,
                                    op="mean",
                                    dim="depth").compute()

# GLODAP profile interpolated to the five comparison depths
glodap_tco2_prof = regrid_vertical(glodap_tco2,
                                   target=profile_depths,
                                   method="linear")

# ---------------------------------------------------------------------------
# Each CMIP6 model
# ---------------------------------------------------------------------------
results = {}
dic_diff_fields = {}

for m in models:
    # --- 2002 monthly fields -------------------------------------------------
    dissic_2002 = load_regular_data(dataset="cmip6_historical",
                                    variable="dissic",
                                    model=m,
                                    time_range=period_2002,
                                    field="total")
    thetao_2002 = load_regular_data(dataset="cmip6_historical",
                                    variable="thetao",
                                    model=m,
                                    time_range=period_2002,
                                    field="total")
    so_2002 = load_regular_data(dataset="cmip6_historical",
                                variable="so",
                                model=m,
                                time_range=period_2002,
                                field="total")
    talk_2002 = load_regular_data(dataset="cmip6_historical",
                                  variable="talk",
                                  model=m,
                                  time_range=period_2002,
                                  depth_range=[0, 50],
                                  field="total")
    ph_2002 = load_regular_data(dataset="cmip6_historical",
                                variable="ph",
                                model=m,
                                time_range=period_2002,
                                depth_range=[0, 50],
                                field="total")

    # In-situ density (sigma) and conversion mol m-3 -> umol kg-1
    sigma_2002 = compute_potential_density_anomaly(thetao_2002,
                                                   so_2002,
                                                   ref_pressure=None)
    dissic_2002_umol = transform_unit(dissic_2002,
                                      from_unit="mol m-3",
                                      to_unit="umol kg-1",
                                      density=sigma_2002)

    sigma_2002_0_50 = sigma_2002.sel(depth=slice(0, 50))
    talk_2002_umol = transform_unit(talk_2002,
                                    from_unit="mol m-3",
                                    to_unit="umol kg-1",
                                    density=sigma_2002_0_50)

    # Annual means for 2002 (day-weighted in the model calendar)
    dissic_2002_annual = temporal_aggregate(dissic_2002_umol,
                                            freq="mean",
                                            op="mean")
    talk_2002_annual = temporal_aggregate(talk_2002_umol,
                                          freq="mean",
                                          op="mean")
    ph_2002_annual = temporal_aggregate(ph_2002,
                                        freq="mean",
                                        op="mean")

    # --- Global 0-50 m volume-weighted means ---------------------------------
    dic_0_50 = dissic_2002_annual.sel(depth=slice(0, 50))
    dic_mean = float(compute_basic_stat(dic_0_50,
                                        op="volume_weighted_mean",
                                        dim=None).compute())
    talk_mean = float(compute_basic_stat(talk_2002_annual,
                                         op="volume_weighted_mean",
                                         dim=None).compute())
    ph_mean = float(compute_basic_stat(ph_2002_annual,
                                       op="volume_weighted_mean",
                                       dim=None).compute())

    # --- 0-50 m DIC map and comparison on the GLODAP 1-degree grid -----------
    dic_map_native = compute_basic_stat(dic_0_50,
                                        op="mean",
                                        dim="depth").compute()
    dic_map_regrid = regrid_horizontal(dic_map_native,
                                       target="glodap",
                                       method="conservative")

    mask_map = build_comparison_mask([dic_map_regrid, glodap_dic_map],
                                     ocean_only=False)
    mod_map_m = dic_map_regrid.where(mask_map)
    glo_map_m = glodap_dic_map.where(mask_map)

    met_map = compute_difference_metrics(mod_map_m,
                                         glo_map_m,
                                         dim=["lat", "lon"],
                                         weighting="area")
    corr_map = analyze_correlation(mod_map_m,
                                   glo_map_m,
                                   dim=["lat", "lon"],
                                   method="pearson",
                                   weighting="area")

    dic_bias = float(met_map["bias"].compute())
    dic_rmse = float(met_map["rmse"].compute())
    dic_r = float(corr_map["correlation"].compute())
    dic_diff_field = met_map["difference"].compute()

    # --- Vertical DIC bias profile ------------------------------------------
    dic_prof_native = regrid_vertical(dissic_2002_annual,
                                      target=profile_depths,
                                      method="linear")
    dic_prof_regrid = regrid_horizontal(dic_prof_native,
                                        target="glodap",
                                        method="conservative")

    prof_bias = []
    prof_n = []
    for d in profile_depths:
        mod_d = dic_prof_regrid.sel(depth=d)
        glo_d = glodap_tco2_prof.sel(depth=d)

        mask_d = build_comparison_mask([mod_d, glo_d],
                                       ocean_only=False)
        mod_d_m = mod_d.where(mask_d)
        glo_d_m = glo_d.where(mask_d)

        met_d = compute_difference_metrics(mod_d_m,
                                           glo_d_m,
                                           dim=["lat", "lon"],
                                           weighting="area")
        prof_bias.append(float(met_d["bias"].compute()))
        prof_n.append(int(mask_d.sum().compute()))

    # --- 2002 minus 1971-2000 0-50 m DIC change -----------------------------
    dissic_base = load_regular_data(dataset="cmip6_historical",
                                    variable="dissic",
                                    model=m,
                                    time_range=period_base,
                                    depth_range=[0, 50],
                                    field="total")
    thetao_base = load_regular_data(dataset="cmip6_historical",
                                    variable="thetao",
                                    model=m,
                                    time_range=period_base,
                                    depth_range=[0, 50],
                                    field="total")
    so_base = load_regular_data(dataset="cmip6_historical",
                                variable="so",
                                model=m,
                                time_range=period_base,
                                depth_range=[0, 50],
                                field="total")

    sigma_base = compute_potential_density_anomaly(thetao_base,
                                                   so_base,
                                                   ref_pressure=None)
    dissic_base_umol = transform_unit(dissic_base,
                                      from_unit="mol m-3",
                                      to_unit="umol kg-1",
                                      density=sigma_base)

    dissic_base_mean_field = temporal_aggregate(dissic_base_umol,
                                                freq="mean",
                                                op="mean")
    dic_base_mean = float(compute_basic_stat(dissic_base_mean_field,
                                             op="volume_weighted_mean",
                                             dim=None).compute())

    dic_change = dic_mean - dic_base_mean
    dic_diff_glodap = dic_mean - float(glodap_dic_mean)
    ratio = dic_change / dic_diff_glodap if dic_diff_glodap != 0 else np.nan

    results[m] = {
        "dic_mean": dic_mean,
        "talk_mean": talk_mean,
        "ph_mean": ph_mean,
        "dic_diff_glodap": dic_diff_glodap,
        "talk_diff_glodap": talk_mean - float(glodap_talk_mean),
        "ph_diff_glodap": ph_mean - float(glodap_ph_mean),
        "dic_bias": dic_bias,
        "dic_rmse": dic_rmse,
        "dic_r": dic_r,
        "prof_bias": prof_bias,
        "prof_n": prof_n,
        "dic_change": dic_change,
        "ratio": ratio,
        "prof_mae": float(np.mean(np.abs(prof_bias))),
    }
    dic_diff_fields[m] = dic_diff_field

# ---------------------------------------------------------------------------
# Summary output
# ---------------------------------------------------------------------------
print("GLODAPv2.2016b 0-50 m volume-weighted means:")
print(f"  DIC  = {float(glodap_dic_mean):.0f} umol kg-1")
print(f"  TAlk = {float(glodap_talk_mean):.0f} umol kg-1")
print(f"  pH   = {float(glodap_ph_mean):.3f}")

for m in models:
    r = results[m]
    print(f"\n{m}:")
    print(f"  0-50 m DIC  = {r['dic_mean']:.0f} umol kg-1 "
          f"(model - GLODAP = {r['dic_diff_glodap']:+.0f})")
    print(f"  0-50 m TAlk = {r['talk_mean']:.0f} umol kg-1 "
          f"(model - GLODAP = {r['talk_diff_glodap']:+.0f})")
    print(f"  0-50 m pH   = {r['ph_mean']:.3f} "
          f"(model - GLODAP = {r['ph_diff_glodap']:+.3f})")
    print(f"  DIC map vs GLODAP: bias = {r['dic_bias']:+.0f} umol kg-1, "
          f"RMSE = {r['dic_rmse']:.0f} umol kg-1, r = {r['dic_r']:.3f}")
    print("  DIC bias profile (0,100,500,1000,2000 m): "
          + ", ".join(f"{x:+.0f}" for x in r["prof_bias"]) + " umol kg-1")
    print("  valid comparison cells per depth: "
          + ", ".join(str(n) for n in r["prof_n"]))
    print(f"  2002 minus 1971-2000 DIC change = {r['dic_change']:+.0f} umol kg-1")
    print(f"  change / (2002 DIC - GLODAP) = {r['ratio']:.2f}")

# Ranking
ranking = sorted(
    models,
    key=lambda mm: (
        results[mm]["dic_rmse"],
        results[mm]["prof_mae"],
        abs(results[mm]["dic_bias"]),
        -results[mm]["dic_r"],
    ),
)

print("\nRanking by DIC map RMSE, then depth-profile MAE, then |bias|, then correlation:")
for i, mm in enumerate(ranking, 1):
    r = results[mm]
    print(f"  {i}. {mm}: RMSE={r['dic_rmse']:.0f}, "
          f"profile MAE={r['prof_mae']:.0f}, "
          f"bias={r['dic_bias']:+.0f}, r={r['dic_r']:.3f}")

print("\nSign consistency of model - GLODAP errors:")
for m in models:
    r = results[m]
    signs = np.sign([r["dic_diff_glodap"],
                     r["talk_diff_glodap"],
                     r["ph_diff_glodap"]])
    print(f"  {m}: DIC {signs[0]:+.0f}, TAlk {signs[1]:+.0f}, pH {signs[2]:+.0f} "
          + ("-> consistent" if len(set(signs)) == 1 else "-> not all the same sign"))

print("\nNote: the 2002-minus-1971-2000 change uses the model's own fields, "
      "so its mean-state bias cancels; it should not be read as a model error.")
