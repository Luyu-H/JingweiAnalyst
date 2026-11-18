# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def as_float(x):
    if hasattr(x, "compute"):
        x = x.compute()
    return float(np.asarray(x))


def coord_lat_lon(loc):
    coords = loc["coordinates"]
    lat = coords.get("latitude", coords.get("lat"))
    lon = coords.get("longitude", coords.get("lon"))
    return lat, lon


def period_mean_volume(dataset, model, time_range):
    # Not every model publishes a time-varying volcello: MPI-ESM1-2-LR (and
    # CESM2) have a fixed field, and the loader refuses a time_range for those.
    # Only GFDL-ESM4 publishes volcello monthly, so ask for a window only where
    # the loader accepts one and fall back to the static field otherwise.
    try:
        volc = load_volumes_for_cmip6(
            dataset=dataset, model=model, time_range=time_range
        )
    except ValueError:
        volc = load_volumes_for_cmip6(dataset=dataset, model=model)
    if "time" in volc.dims:
        volc = temporal_aggregate(volc, freq="mean")
    return volc


def regrid_for_display(field, source_field, method="bilinear"):
    # Preserve the grid provenance from the source field, but keep thickness units.
    out = field.copy()
    attrs = dict(source_field.attrs)
    attrs["units"] = field.attrs.get("units", "m")
    out.attrs.update(attrs)
    return regrid_horizontal(out, target="woa23", method=method)


# ---------------------------------------------------------------------------
# WOA23 observed climatology
# ---------------------------------------------------------------------------
woa_o2 = load_climatology_data(
    dataset="woa23", variable="oxygen", period="1971-2000", time_slice="annual"
)
woa_aou = load_climatology_data(
    dataset="woa23", variable="apparent_oxygen_utilization",
    period="1971-2000", time_slice="annual"
)
woa_o2sat = load_climatology_data(
    dataset="woa23", variable="percent_oxygen_saturation",
    period="1971-2000", time_slice="annual"
)

landsea = load_landsea_mask()
basinmask = load_basin_mask()
woa_levels, woa_bounds = load_depth_levels(
    dataset="woa23", variable="oxygen", time_slice="annual", period="1971-2000"
)
woa_areas = load_areas(dataset="woa23", variable="oxygen")

woa_ocean = landsea > 1
woa_o2 = woa_o2.where(woa_ocean)
woa_aou = woa_aou.where(woa_ocean)
woa_o2sat = woa_o2sat.where(woa_ocean)

woa_thick20, woa_mask20 = identify_oxygen_minimum_zone(woa_o2, threshold=20.0)
woa_thick60, woa_mask60 = identify_oxygen_minimum_zone(woa_o2, threshold=60.0)

woa_ex20 = compute_exceedance_fraction(
    data=woa_o2, threshold=20.0, comparison="less", dim=None, weighting="volume"
)
woa_ex60 = compute_exceedance_fraction(
    data=woa_o2, threshold=60.0, comparison="less", dim=None, weighting="volume"
)
woa_vol20 = as_float(woa_ex20["exceeding_amount"])
woa_vol60 = as_float(woa_ex60["exceeding_amount"])

woa_loc60 = locate_extreme(woa_thick60, op="max", dim=None)
woa_max60 = float(woa_loc60["value"])
woa_max60_lat, woa_max60_lon = coord_lat_lon(woa_loc60)

woa_aou_omz = compute_basic_stat(
    woa_aou.where(woa_mask60), op="volume_weighted_mean", dim=None
).compute()
woa_o2sat_omz = compute_basic_stat(
    woa_o2sat.where(woa_mask60), op="volume_weighted_mean", dim=None
).compute()

woa_results = {
    "vol20": woa_vol20,
    "vol60": woa_vol60,
    "thick20": woa_thick20,
    "thick60": woa_thick60,
    "mask20": woa_mask20,
    "mask60": woa_mask60,
    "loc60": woa_loc60,
    "max60": woa_max60,
    "max60_lat": woa_max60_lat,
    "max60_lon": woa_max60_lon,
    "aou_omz60": as_float(woa_aou_omz),
    "o2sat_omz60": as_float(woa_o2sat_omz),
}


# ---------------------------------------------------------------------------
# CMIP6 processing
# ---------------------------------------------------------------------------
def process_cmip6(dataset, model, time_range):
    o2 = load_regular_data(
        dataset=dataset, variable="o2", model=model, time_range=time_range
    )
    thetao = load_regular_data(
        dataset=dataset, variable="thetao", model=model, time_range=time_range
    )
    so = load_regular_data(
        dataset=dataset, variable="so", model=model, time_range=time_range
    )

    rho = compute_potential_density_anomaly(
        temp_data=thetao, sal_data=so, ref_pressure=None, temperature_type="auto"
    )
    o2_umolkg = transform_unit(
        o2, from_unit="mol m-3", to_unit="umol/kg", density=rho
    )

    o2_mean = temporal_aggregate(o2_umolkg, freq="mean")
    thetao_mean = temporal_aggregate(thetao, freq="mean")
    so_mean = temporal_aggregate(so, freq="mean")
    volc_mean = period_mean_volume(dataset, model, time_range)

    # Dissolved oxygen inventory on the volume basis: sum(o2 * volcello), no density.
    o2_content = convert_concentration_to_content(concentration=o2, density=1025.0)
    o2_content_mean = temporal_aggregate(o2_content, freq="mean")
    inventory = compute_basic_stat(o2_content_mean, op="sum", dim=None).compute()

    thick20, mask20 = identify_oxygen_minimum_zone(o2_mean, threshold=20.0)
    thick60, mask60 = identify_oxygen_minimum_zone(o2_mean, threshold=60.0)

    ex20 = compute_exceedance_fraction(
        data=o2_mean, threshold=20.0, comparison="less", dim=None, weighting="volume"
    )
    ex60 = compute_exceedance_fraction(
        data=o2_mean, threshold=60.0, comparison="less", dim=None, weighting="volume"
    )
    vol20 = as_float(ex20["exceeding_amount"])
    vol60 = as_float(ex60["exceeding_amount"])

    loc60 = locate_extreme(thick60, op="max", dim=None)
    max60 = float(loc60["value"])
    max60_lat, max60_lon = coord_lat_lon(loc60)

    return {
        "o2": o2,
        "thetao": thetao,
        "so": so,
        "rho": rho,
        "o2_umolkg": o2_umolkg,
        "o2_mean": o2_mean,
        "thetao_mean": thetao_mean,
        "so_mean": so_mean,
        "volc_mean": volc_mean,
        "inventory": as_float(inventory),
        "thick20": thick20,
        "thick60": thick60,
        "mask20": mask20,
        "mask60": mask60,
        "vol20": vol20,
        "vol60": vol60,
        "loc60": loc60,
        "max60": max60,
        "max60_lat": max60_lat,
        "max60_lon": max60_lon,
    }


hist_period = ["1971-01", "2000-12"]
ssp_period = ["2071-01", "2100-12"]
hist_models = ["GFDL-ESM4", "MPI-ESM1-2-LR"]
ssp_models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

hist_results = {}
for model in hist_models:
    hist_results[model] = process_cmip6("cmip6_historical", model, hist_period)

ssp_results = {}
for model in ssp_models:
    ssp_results[model] = process_cmip6("cmip6_ssp245", model, ssp_period)


# ---------------------------------------------------------------------------
# Basin OMZ volumes at the 60 umol/kg threshold
# ---------------------------------------------------------------------------
basins = ["Pacific", "Indian", "Atlantic"]


def basin_omz_volume(o2_field, basin):
    f = select_basin(
        data=o2_field,
        basin=basin,
        include_marginal_seas=True,
        depth_dependent=False,
    )
    ex = compute_exceedance_fraction(
        data=f, threshold=60.0, comparison="less", dim=None, weighting="volume"
    )
    return as_float(ex["exceeding_amount"])


woa_basin_vol60 = {}
for b in basins:
    woa_basin_vol60[b] = basin_omz_volume(woa_o2, b)

hist_basin_vol60 = {}
for model in hist_models:
    hist_basin_vol60[model] = {}
    for b in basins:
        hist_basin_vol60[model][b] = basin_omz_volume(hist_results[model]["o2_mean"], b)


# ---------------------------------------------------------------------------
# Scenario changes at 20 and 60 umol/kg
# ---------------------------------------------------------------------------
scenario_change_pct = {}
for model in hist_models:
    scenario_change_pct[model] = {}
    for thr, key in [(20, "vol20"), (60, "vol60")]:
        h = hist_results[model][key]
        s = ssp_results[model][key]
        scenario_change_pct[model][thr] = (s - h) / h * 100.0


# ---------------------------------------------------------------------------
# Historical model vs WOA23 volume comparison
# ---------------------------------------------------------------------------
model_vs_woa = {}
for model in hist_models:
    model_vs_woa[model] = {}
    for thr, key in [(20, "vol20"), (60, "vol60")]:
        mv = hist_results[model][key]
        ov = woa_results[key]
        rel = (mv - ov) / ov * 100.0
        model_vs_woa[model][thr] = {
            "model_volume": mv,
            "woa_volume": ov,
            "relative_difference_pct": rel,
            "absolute_relative_difference_pct": abs(rel),
        }

# Mean absolute relative difference across the two thresholds; lower is closer.
model_closeness = {
    model: float(
        np.mean([model_vs_woa[model][20]["absolute_relative_difference_pct"],
                model_vs_woa[model][60]["absolute_relative_difference_pct"]])
    )
    for model in hist_models
}
closer_model = min(model_closeness, key=model_closeness.get)

# ---------------------------------------------------------------------------
# Display fields (figures only; all reported values remain native-grid)
# ---------------------------------------------------------------------------
fig1_fields = {
    "GFDL-ESM4 historical": regrid_for_display(
        hist_results["GFDL-ESM4"]["thick60"], hist_results["GFDL-ESM4"]["o2_mean"]
    ),
    "MPI-ESM1-2-LR historical": regrid_for_display(
        hist_results["MPI-ESM1-2-LR"]["thick60"], hist_results["MPI-ESM1-2-LR"]["o2_mean"]
    ),
    "WOA23 1971-2000": woa_thick60,
}

fig2_change_fields = {}
for model in hist_models:
    ssp_display = regrid_for_display(
        ssp_results[model]["thick60"], ssp_results[model]["o2_mean"]
    )
    hist_display = regrid_for_display(
        hist_results[model]["thick60"], hist_results[model]["o2_mean"]
    )
    dm = compute_difference_metrics(
        data=ssp_display, reference=hist_display, dim=None, weighting="none"
    )
    fig2_change_fields[model] = dm["difference"]

cesm2_ssp_display = regrid_for_display(
    ssp_results["CESM2"]["thick60"], ssp_results["CESM2"]["o2_mean"]
)
