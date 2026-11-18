# ============================================================================
# OMZ assessment: CMIP6 historical / ssp245 against WOA23
# ============================================================================

hist_models = ["GFDL-ESM4", "MPI-ESM1-2-LR"]
ssp_models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
hist_time = ["1971-01", "2000-12"]
ssp_time = ["2071-01", "2100-12"]


def load_period_o2(dataset, model, time_range):
    """Load monthly o2/thetao/so, convert o2 to umol kg-1 with in-situ density,
    and form day-weighted period means on both mass and volume bases."""
    o2 = load_regular_data(dataset=dataset, variable="o2", model=model,
                           time_range=time_range)
    thetao = load_regular_data(dataset=dataset, variable="thetao", model=model,
                               time_range=time_range)
    so = load_regular_data(dataset=dataset, variable="so", model=model,
                           time_range=time_range)

    rho = compute_potential_density_anomaly(
        temp_data=thetao, sal_data=so, ref_pressure=None,
        temperature_type="potential"
    )
    o2_mass = transform_unit(
        data=o2, from_unit="mol m-3", to_unit="umol kg-1", density=rho
    )
    o2_mass_mean = temporal_aggregate(o2_mass, freq="mean", op="mean").compute()
    o2_raw_mean = temporal_aggregate(o2, freq="mean", op="mean").compute()
    dens_assump = o2_mass.attrs.get("density_assumption",
                                    "in-situ density from TEOS-10")
    return o2_mass_mean, o2_raw_mean, dens_assump


cmip = {}

for dataset, models, time_range, tag in [
    ("cmip6_historical", hist_models, hist_time, "hist"),
    ("cmip6_ssp245", ssp_models, ssp_time, "ssp245"),
]:
    for model in models:
        o2_mean, o2_raw_mean, dens_assump = load_period_o2(
            dataset, model, time_range
        )

        thick20, mask20 = identify_oxygen_minimum_zone(o2_mean, threshold=20.0)
        thick60, mask60 = identify_oxygen_minimum_zone(o2_mean, threshold=60.0)
        thick20 = thick20.compute()
        thick60 = thick60.compute()
        mask20 = mask20.compute()
        mask60 = mask60.compute()

        res20 = compute_exceedance_fraction(
            o2_mean, threshold=20.0, comparison="less",
            dim=None, weighting="volume"
        )
        res60 = compute_exceedance_fraction(
            o2_mean, threshold=60.0, comparison="less",
            dim=None, weighting="volume"
        )
        vol20 = float(res20["exceeding_amount"]) / 1e15
        vol60 = float(res60["exceeding_amount"]) / 1e15

        content = convert_concentration_to_content(concentration=o2_raw_mean)
        inventory = float(content.sum(skipna=True).compute())

        ext60 = locate_extreme(thick60, op="max")
        ext_coords = ext60["coordinates"]
        max_lat = ext_coords.get("latitude", ext_coords.get("lat"))
        max_lon = ext_coords.get("longitude", ext_coords.get("lon"))

        cmip[(tag, model)] = {
            "o2_mean": o2_mean,
            "o2_raw_mean": o2_raw_mean,
            "density_assumption": dens_assump,
            "thick20": thick20,
            "thick60": thick60,
            "mask20": mask20,
            "mask60": mask60,
            "vol20": vol20,
            "vol60": vol60,
            "inventory": inventory,
            "max_thick60": ext60["value"],
            "max_thick60_lat": max_lat,
            "max_thick60_lon": max_lon,
        }


# ----------------------------------------------------------------------------
# WOA23 observational reference
# ----------------------------------------------------------------------------
woa_o2 = load_climatology_data(dataset="woa23", variable="oxygen",
                               time_slice="annual")
woa_aou = load_climatology_data(dataset="woa23",
                                variable="apparent_oxygen_utilization",
                                time_slice="annual")
woa_sat = load_climatology_data(dataset="woa23",
                                variable="percent_oxygen_saturation",
                                time_slice="annual")

woa_thick20, woa_mask20 = identify_oxygen_minimum_zone(woa_o2, threshold=20.0)
woa_thick60, woa_mask60 = identify_oxygen_minimum_zone(woa_o2, threshold=60.0)
woa_thick20 = woa_thick20.compute()
woa_thick60 = woa_thick60.compute()
woa_mask20 = woa_mask20.compute()
woa_mask60 = woa_mask60.compute()

woa_res20 = compute_exceedance_fraction(
    woa_o2, threshold=20.0, comparison="less", dim=None, weighting="volume"
)
woa_res60 = compute_exceedance_fraction(
    woa_o2, threshold=60.0, comparison="less", dim=None, weighting="volume"
)
woa_vol20 = float(woa_res20["exceeding_amount"]) / 1e15
woa_vol60 = float(woa_res60["exceeding_amount"]) / 1e15

woa_ext60 = locate_extreme(woa_thick60, op="max")
woa_ext_coords = woa_ext60["coordinates"]
woa_max_lat = woa_ext_coords.get("latitude", woa_ext_coords.get("lat"))
woa_max_lon = woa_ext_coords.get("longitude", woa_ext_coords.get("lon"))

# AOU and saturation inside the WOA23 <60 umol kg-1 OMZ
woa_aou_omz = woa_aou.where(woa_mask60)
woa_sat_omz = woa_sat.where(woa_mask60)
woa_aou_mean = float(
    compute_basic_stat(
        woa_aou_omz, op="volume_weighted_mean", dim=["depth", "lat", "lon"]
    ).compute()
)
woa_sat_mean = float(
    compute_basic_stat(
        woa_sat_omz, op="volume_weighted_mean", dim=["depth", "lat", "lon"]
    ).compute()
)


# ----------------------------------------------------------------------------
# Basin volumes at the 60 umol kg-1 threshold
# ----------------------------------------------------------------------------
basin_names = ["Pacific", "Indian", "Atlantic"]


def basin_vol60(o2_field):
    out = {}
    for b in basin_names:
        o2_b = select_basin(
            o2_field, basin=b, include_marginal_seas=True,
            depth_dependent=False
        )
        r = compute_exceedance_fraction(
            o2_b, threshold=60.0, comparison="less",
            dim=None, weighting="volume"
        )
        out[b] = float(r["exceeding_amount"]) / 1e15
    return out


woa_basin_vol60 = basin_vol60(woa_o2)
for model in hist_models:
    cmip[("hist", model)]["basin_vol60"] = basin_vol60(
        cmip[("hist", model)]["o2_mean"]
    )


# ----------------------------------------------------------------------------
# Historical -> ssp245 changes for models with both experiments
# ----------------------------------------------------------------------------
pct_change = {}
for model in hist_models:
    h = cmip[("hist", model)]
    s = cmip[("ssp245", model)]
    pct_change[model] = {
        "pct20": 100.0 * (s["vol20"] - h["vol20"]) / h["vol20"],
        "pct60": 100.0 * (s["vol60"] - h["vol60"]) / h["vol60"],
    }


# ----------------------------------------------------------------------------
# Comparison of historical OMZ volume with WOA23
# ----------------------------------------------------------------------------
woa_comparison = {}
for model in hist_models:
    h = cmip[("hist", model)]
    err20 = 100.0 * abs(h["vol20"] - woa_vol20) / woa_vol20
    err60 = 100.0 * abs(h["vol60"] - woa_vol60) / woa_vol60
    woa_comparison[model] = {
        "err20_pct": err20,
        "err60_pct": err60,
        "mean_abs_err_pct": 0.5 * (err20 + err60),
    }

best_model = min(woa_comparison,
                 key=lambda m: woa_comparison[m]["mean_abs_err_pct"])


# ----------------------------------------------------------------------------
# Regrid historical 60-threshold thickness fields to WOA23 and compare
# ----------------------------------------------------------------------------
gfdl_hist_reg = regrid_horizontal(
    cmip[("hist", "GFDL-ESM4")]["thick60"],
    target="woa23", method="conservative", skipna=True
)
mpi_hist_reg = regrid_horizontal(
    cmip[("hist", "MPI-ESM1-2-LR")]["thick60"],
    target="woa23", method="conservative", skipna=True
)

comparison_mask = build_comparison_mask(
    fields=[gfdl_hist_reg, mpi_hist_reg, woa_thick60],
    ocean_only=True
)
gfdl_hist_m = gfdl_hist_reg.where(comparison_mask)
mpi_hist_m = mpi_hist_reg.where(comparison_mask)
woa_thick60_m = woa_thick60.where(comparison_mask)

metrics_gfdl = compute_difference_metrics(
    data=gfdl_hist_m, reference=woa_thick60_m,
    dim=["lat", "lon"], weighting="area"
)
metrics_mpi = compute_difference_metrics(
    data=mpi_hist_m, reference=woa_thick60_m,
    dim=["lat", "lon"], weighting="area"
)


# ----------------------------------------------------------------------------
# Fields for the scenario-change figure
# ----------------------------------------------------------------------------
gfdl_ssp_reg = regrid_horizontal(
    cmip[("ssp245", "GFDL-ESM4")]["thick60"],
    target="woa23", method="conservative", skipna=True
)
mpi_ssp_reg = regrid_horizontal(
    cmip[("ssp245", "MPI-ESM1-2-LR")]["thick60"],
    target="woa23", method="conservative", skipna=True
)
cesm_ssp_reg = regrid_horizontal(
    cmip[("ssp245", "CESM2")]["thick60"],
    target="woa23", method="conservative", skipna=True
)

change_gfdl = (gfdl_ssp_reg - gfdl_hist_reg).compute()
change_mpi = (mpi_ssp_reg - mpi_hist_reg).compute()


# ----------------------------------------------------------------------------
# Printed summary
# ----------------------------------------------------------------------------
print("=== OMZ volumes and maximum thickness at O2 < 60 umol kg-1 ===")
print(
    f"WOA23: vol<20={woa_vol20:.3g}e15 m3, vol<60={woa_vol60:.3g}e15 m3, "
    f"max thickness={woa_ext60['value']:.1f} m "
    f"at ({woa_max_lat:.2f}, {woa_max_lon:.2f})"
)

for tag, models in [("hist", hist_models), ("ssp245", ssp_models)]:
    for model in models:
        d = cmip[(tag, model)]
        print(
            f"{model} {tag}: vol<20={d['vol20']:.3g}e15 m3, "
            f"vol<60={d['vol60']:.3g}e15 m3, "
            f"max thickness={d['max_thick60']:.1f} m "
            f"at ({d['max_thick60_lat']:.2f}, {d['max_thick60_lon']:.2f})"
        )

print("\n=== Dissolved oxygen inventory (mol) ===")
for tag, models in [("hist", hist_models), ("ssp245", ssp_models)]:
    for model in models:
        print(f"{model} {tag}: {cmip[(tag, model)]['inventory']:.3g} mol")

print("\n=== Percentage change historical -> ssp245 ===")
for model in hist_models:
    print(
        f"{model}: <20 {pct_change[model]['pct20']:.2f}%, "
        f"<60 {pct_change[model]['pct60']:.2f}%"
    )

print("\n=== Basin OMZ volumes at O2 < 60 umol kg-1 ===")
print("WOA23:")
for b in basin_names:
    print(f"  {b}: {woa_basin_vol60[b]:.3g}e15 m3")
for model in hist_models:
    print(f"{model} historical:")
    for b in basin_names:
        print(f"  {b}: {cmip[('hist', model)]['basin_vol60'][b]:.3g}e15 m3")

print("\n=== WOA23 AOU and saturation inside the <60 umol kg-1 OMZ ===")
print(f"AOU = {woa_aou_mean:.1f} umol kg-1")
print(f"Percent oxygen saturation = {woa_sat_mean:.1f}%")

print("\n=== Historical model vs WOA23 OMZ volume ===")
for model in hist_models:
    c = woa_comparison[model]
    print(
        f"{model}: err<20={c['err20_pct']:.2f}%, "
        f"err<60={c['err60_pct']:.2f}%, "
        f"mean abs err={c['mean_abs_err_pct']:.2f}%"
    )
print(f"Better historical OMZ volume by mean absolute error: {best_model}")

print("\n=== Historical 60-threshold OMZ thickness metrics vs WOA23 ===")
for name, met in [("GFDL-ESM4", metrics_gfdl),
                  ("MPI-ESM1-2-LR", metrics_mpi)]:
    print(
        f"{name}: n={int(met['n_valid'])}, "
        f"bias={float(met['bias']):.2f} m, "
        f"rmse={float(met['rmse']):.2f} m, "
        f"max_abs_error={float(met['max_abs_error']):.2f} m, "
        f"mean_data={float(met['mean_data']):.2f} m, "
        f"mean_reference={float(met['mean_reference']):.2f} m"
    )
    print(f"  max_abs_error_at = {met['max_abs_error_at']}")

print("\nNote: WOA23 oxygen is used as the available all-period composite; "
      "no separate 1971-2000 oxygen normal is assumed.")
print("Note: CMIP6 o2 was converted monthly to umol kg-1 using in-situ "
      "TEOS-10 density, then period-averaged with month-length weights.")
print("Note: CESM2 historical o2 is not published; no CESM2 historical-to-"
      "ssp245 change is formed (ssp245 is r10, historical r10 unavailable).")
