"""
Category 7: safeguards, by fault injection, and regression tests.
"""
from __future__ import annotations

import contextlib
import functools
import io

import numpy as np
import xarray as xr

from src.eval.tool_verification.common import case, coord_like, monthly_axis, refill
from src.tools.basic_calculator import (analyze_correlation, compute_basic_stat, compute_difference_metrics,
                                        compute_exceedance_fraction, detect_outliers, detect_threshold_events,
                                        fit_1D_trend, fit_multivariate_relation, fit_trend_map, locate_extreme,
                                        temporal_aggregate, transform_unit)
from src.tools.data_loader import (load_basin_mask, load_climatology_data, load_data_with_geolocator,
                                   load_depth_levels, load_regular_data, load_volumes_for_cmip6)
from src.tools.preprocessing import build_comparison_mask, regrid_horizontal, regrid_vertical, select_basin
from src.tools.scientific_calculator import (compute_anomaly_vs_own_climatology, compute_anomaly_vs_woa23,
                                             compute_enso_nino_index, compute_mixed_layer_depth,
                                             compute_nutricline_depth, compute_pdo_index,
                                             compute_potential_density_anomaly, compute_ratio,
                                             compute_vertical_gradient, convert_concentration_to_content,
                                             identify_oxygen_minimum_zone)

CAT = "cat7_safeguards"


@functools.lru_cache(None)
def woa_sst():
    return load_climatology_data("woa23", "temperature", depth_range=[0, 0]).isel(depth=0).load()


@functools.lru_cache(None)
def woa_profile(variable="temperature"):
    return load_climatology_data("woa23", variable, lat_range=[10, 12], lon_range=[-30, -28]).load()


@functools.lru_cache(None)
def ersst_short():
    return load_regular_data("ersst", "sst", time_range=["2000-01", "2001-12"]).load()


@functools.lru_cache(None)
def cesm_tos():
    return load_regular_data("cmip6_historical", "tos", model="CESM2", time_range=["2000-01", "2000-02"]).load()


def series(n=60, values=None):
    values = np.arange(n, dtype=float) if values is None else values
    return xr.DataArray(values, dims="time", coords={"time": monthly_axis("2000-01", len(values))})


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

@case("load_regular_data", "refusals")
def guard_load_regular(rec):
    rec.raises("index dataset with a lat_range", lambda: load_regular_data("cpc_indices", "NINO3.4", lat_range=[-5, 5]),
               match="ungridded")
    rec.raises("reversed lat_range", lambda: load_regular_data("ersst", "sst", lat_range=[10, -10]),
               match="low to high")
    rec.raises("time_range outside the record", lambda: load_regular_data("ersst", "sst", time_range=["1700-01", "1701-01"]),
               match="selects no steps")
    rec.raises("a climatology through the time-series loader", lambda: load_regular_data("woa23", "temperature"),
               match="load_climatology_data")
    rec.raises("a multi-model dataset without a model", lambda: load_regular_data("cmip6_historical", "tos"),
               match="needs a model")


@case("load_regular_data", "regression: a dateline window comes back on a continuous, increasing axis")
def regression_dateline(rec):
    da = load_regular_data("ersst", "sst", time_range=["2000-01", "2000-01"], lat_range=[-5, 5], lon_range=[160, -150])
    lon = da["lon"].values
    rec.metric("longitudes strictly increasing", bool(np.all(np.diff(lon) > 0)), op="true", tier="regression",
               note=f"{lon.min():g}..{lon.max():g}")
    rec.metric("window is contiguous (2-degree steps only)", bool(np.allclose(np.diff(lon), 2.0)), op="true",
               tier="regression")
    mean = compute_basic_stat(da, "area_weighted_mean", dim=["lat", "lon"])
    rec.metric("and can still be area-weighted", bool(np.isfinite(mean).all()), op="true", tier="regression")


@case("load_climatology_data", "refusals")
def guard_load_clim(rec):
    rec.raises("months on an annual climatology", lambda: load_climatology_data("woa23", "temperature", months=[1]),
               match="month")
    rec.raises("a time series through the climatology loader", lambda: load_climatology_data("ersst", "sst"),
               match="load_regular_data")


@case("load_data_with_geolocator", "refusals and the western-hemisphere regression")
def guard_geolocator(rec):
    rec.raises("no loc_name", lambda: load_data_with_geolocator("woa23", "temperature"), match="loc_name")
    rec.raises("an index dataset", lambda: load_data_with_geolocator("cpc_indices", "NINO3.4", loc_name="Red Sea"),
               match="no grid to mask")
    gulf = load_data_with_geolocator("woa23", "temperature", loc_name="Gulf of Mexico", depth_range=[0, 0])
    rec.metric("regression: a western-hemisphere sea is not empty (0-360 vs -180..180 grid)",
               int(gulf.notnull().sum()), threshold=50, op=">=", tier="regression")


@case("load_volumes_for_cmip6", "refusals")
def guard_volumes(rec):
    rec.raises("time_range for a fixed volcello", lambda: load_volumes_for_cmip6("cmip6_historical", "CESM2",
                                                                                time_range=["2000-01", "2000-12"]),
               match="fixed field")
    rec.raises("a non-CMIP6 dataset", lambda: load_volumes_for_cmip6("woa23", "CESM2"), match="not a CMIP6")


@case("load_depth_levels", "refusals")
def guard_depth_levels(rec):
    rec.raises("a depth_range between levels", lambda: load_depth_levels(depth_range=[1.0, 2.0]), match="selects no levels")
    rec.raises("a surface-only dataset", lambda: load_depth_levels(dataset="ersst"), match="no vertical")


@case("load_basin_mask", "refusals")
def guard_basin_mask(rec):
    rec.raises("a depth that is not a standard level", lambda: load_basin_mask("7m"), match="not found")


# ---------------------------------------------------------------------------
# Preprocessing
# ---------------------------------------------------------------------------

@case("regrid_horizontal", "refusals")
def guard_regrid_h(rec):
    uo = load_regular_data("cmip6_historical", "uo", model="CESM2", time_range=["2000-01", "2000-01"],
                           depth_range=[0, 10])
    rec.raises("CESM2 uo (a grid-relative vector component)", lambda: regrid_horizontal(uo.isel(time=0)),
               match="component")
    box = load_regular_data("cmip6_historical", "tos", model="CESM2", time_range=["2000-01", "2000-01"],
                            lat_range=[0, 20], lon_range=[200, 230])
    rec.raises("a curvilinear subset with method='conservative'", lambda: regrid_horizontal(box.isel(time=0)),
               match="global field")
    ok = regrid_horizontal(box.isel(time=0), method="bilinear")
    rec.metric("no false alarm: the same subset with method='bilinear' is accepted", int(ok.notnull().sum()),
               threshold=1, op=">=", tier="safeguard")
    rec.raises("an unknown method", lambda: regrid_horizontal(woa_sst(), method="spline"), match="method")
    rec.raises("na_thres outside 0..1", lambda: regrid_horizontal(woa_sst(), na_thres=2.0), match="na_thres")
    rec.raises("a curvilinear template", lambda: regrid_horizontal(woa_sst(), target=cesm_tos().isel(time=0)),
               match="rectilinear")


@case("regrid_vertical", "refusals and the on-level regression")
def guard_regrid_v(rec):
    rec.raises("negative target depths", lambda: regrid_vertical(woa_profile(), target=[-10.0, 5.0]), match="positive")
    rec.raises("a field without depth", lambda: regrid_vertical(woa_sst()), match="depth")
    rec.raises("an unknown method", lambda: regrid_vertical(woa_profile(), method="cubic"), match="method")
    p = woa_profile()
    z = coord_like(p, "depth")
    field = refill(p, np.where(z <= 200, 5.0 + 0.0 * z, np.nan))
    out = regrid_vertical(field, target=[200.0])
    rec.metric("regression: a target exactly on the deepest valid level keeps that level's value",
               bool(np.all(out.values == 5.0)), op="true", tier="regression",
               note="linear interpolation at a level would otherwise multiply the NaN below by zero")


@case("select_basin", "refusals")
def guard_select_basin(rec):
    rec.raises("an unknown basin", lambda: select_basin(woa_sst(), "Atlantis"), match="not a basin")
    rec.raises("depth_dependent without a depth dimension", lambda: select_basin(woa_sst(), "Atlantic", depth_dependent=True),
               match="depth")
    pacific_box = load_climatology_data("woa23", "temperature", depth_range=[0, 0], lat_range=[-10, 10],
                                        lon_range=[-170, -120]).isel(depth=0)
    rec.raises("a basin absent from a regional field", lambda: select_basin(pacific_box, "Atlantic"),
               match="covers no cell")
    rec.raises("an unknown code", lambda: select_basin(woa_sst(), [99]), match="Unknown")


@case("build_comparison_mask", "refusals and the third-product coastline warning")
def guard_comparison_mask(rec):
    rec.raises("a single field instead of a list", lambda: build_comparison_mask(woa_sst()), match="LIST")
    rec.raises("fields on different grids", lambda: build_comparison_mask([woa_sst(), ersst_short().isel(time=0)]),
               match="different")
    rec.raises("a curvilinear field", lambda: build_comparison_mask([cesm_tos().isel(time=0), woa_sst()]),
               match="curvilinear")
    glodap = load_climatology_data("glodap", "temperature", depth_range=[0, 0]).isel(depth=0)
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        mask = build_comparison_mask([glodap, glodap.copy()])
    rec.metric("GLODAP-only comparison: the WOA23 coastline use is flagged", "WARNING" in mask.attrs["ocean_mask"]
               and "Warning" in buffer.getvalue(), op="true", tier="safeguard", note=mask.attrs["ocean_mask"][:160])
    with contextlib.redirect_stdout(io.StringIO()):
        mask = build_comparison_mask([woa_sst(), woa_sst().copy()])
    rec.metric("no false alarm when WOA23 is in the comparison", "WARNING" not in mask.attrs["ocean_mask"], op="true",
               tier="safeguard")


# ---------------------------------------------------------------------------
# Basic calculators
# ---------------------------------------------------------------------------

@case("compute_basic_stat", "refusals")
def guard_basic_stat(rec):
    sst = woa_sst()
    rec.raises("area_weighted_mean over longitude only", lambda: compute_basic_stat(sst, "area_weighted_mean", dim="lon"),
               match="has to reduce")
    shifted = sst.assign_coords(lat=sst["lat"] + 0.25)
    rec.raises("a field moved off its own grid", lambda: compute_basic_stat(shifted, "area_weighted_mean", dim=["lat", "lon"]),
               match="original grid")
    bare = xr.DataArray(sst.values, dims=sst.dims, coords=sst.coords)
    rec.raises("a field that lost its provenance", lambda: compute_basic_stat(bare, "area_weighted_mean", dim=["lat", "lon"]),
               match="source_dataset")
    rec.raises("argmax over two dimensions", lambda: compute_basic_stat(sst, "argmax", dim=["lat", "lon"]),
               match="ONE axis")
    rec.raises("volume_weighted_mean without depth", lambda: compute_basic_stat(sst, "volume_weighted_mean"),
               match="depth")


@case("compute_basic_stat", "regression: CESM2 depth subsets (float32 levels) keep their deepest level")
def regression_cesm_depth(rec):
    thetao = load_regular_data("cmip6_historical", "thetao", model="CESM2", time_range=["2000-01", "2000-01"],
                               depth_range=[0, 1000], lat_range=[0, 3], lon_range=[200, 203]).isel(time=0).load()
    for label, call in (("volume_weighted_mean", lambda: compute_basic_stat(thetao, "volume_weighted_mean")),
                        ("compute_vertical_gradient", lambda: compute_vertical_gradient(thetao)),
                        ("compute_nutricline_depth", lambda: compute_nutricline_depth(refill(thetao, -thetao.values)))):
        try:
            out = call()
            rec.metric(f"{label} on a 0-1000 m CESM2 subset runs", bool(np.isfinite(out).any()), op="true",
                       tier="regression", note=f"deepest level {float(thetao['depth'].max()):.6f} m")
        except Exception as exc:
            rec.metric(f"{label} on a 0-1000 m CESM2 subset runs", False, op="true", tier="regression",
                       note=f"{type(exc).__name__}: {exc}"[:300])


@case("temporal_aggregate", "refusals")
def guard_temporal(rec):
    clim = load_climatology_data("woa23", "temperature", time_slice="monthly", depth_range=[0, 0])
    rec.raises("a climatology (month dimension)", lambda: temporal_aggregate(clim), match="climatology")
    rec.raises("rolling without a window", lambda: temporal_aggregate(series(), freq="rolling"), match="window")
    rec.raises("a window longer than the series", lambda: temporal_aggregate(series(10), freq="rolling", window=12),
               match="does not fit")
    rec.raises("an unknown op", lambda: temporal_aggregate(series(), op="sum"), match="not supported")
    partial = temporal_aggregate(series(18), freq="year")
    rec.metric("an incomplete year is flagged", "incomplete_periods" in partial.attrs, op="true", tier="safeguard",
               note=partial.attrs.get("incomplete_periods", "")[:120])


@case("fit_1D_trend", "refusals")
def guard_fit_1d(rec):
    rec.raises("fewer than three finite values", lambda: fit_1D_trend(series(10, np.r_[1.0, 2.0, [np.nan] * 8])),
               match="too few")
    rec.raises("a 2-D input", lambda: fit_1D_trend(woa_sst()), match="1D")
    rec.raises("exponential on a series that changes sign", lambda: fit_1D_trend(series(20, np.linspace(-1, 1, 20)),
                                                                                method="exponential"),
               match="exponential")


@case("fit_trend_map", "refusals")
def guard_trend_map(rec):
    rec.raises("a single series", lambda: fit_trend_map(series()), match="fit_1D_trend")
    rec.raises("a missing dimension", lambda: fit_trend_map(ersst_short(), dim="year"), match="not found")


@case("transform_unit", "refusals and the documented double-1000 limitation")
def guard_transform(rec):
    x = xr.DataArray([1.0, 2.0], dims="n")
    rec.raises("density given as 1025 (kg m-3, not kg/L)", lambda: transform_unit(x, "umol/kg", "mol m-3", density=1025.0),
               match="plausible")
    rec.raises("a cross-category conversion", lambda: transform_unit(x, "degC", "m"), match="no conversion")
    rec.raises("an unknown unit", lambda: transform_unit(x, "furlong", "m"), match="Unsupported")
    sigma = xr.DataArray([25.0, 26.0], dims="n", attrs={"long_name": "Potential density anomaly (sigma_0)"})
    right = transform_unit(x, "umol/kg", "mmol m-3", density=sigma)
    pre_added = sigma + 1000.0
    pre_added.attrs = dict(sigma.attrs)
    wrong = transform_unit(x, "umol/kg", "mmol m-3", density=pre_added)
    rec.metric("documented limitation: sigma + 1000 still labelled 'anomaly' -> converted value ratio",
               float(wrong[0] / right[0]), op="info", tier="safeguard",
               note="the docstring warns of this; ~2 means 1000 was added twice")


@case("fit_multivariate_relation", "refusals")
def guard_multivariate(rec):
    a = xr.DataArray(np.random.default_rng(0).normal(size=(5, 6)), dims=("lat", "lon"),
                     coords={"lat": np.arange(5.0), "lon": np.arange(6.0)})
    b = a.assign_coords(lat=a["lat"] + 1)
    rec.raises("predictor on other coordinates", lambda: fit_multivariate_relation(a, [b]), match="coordinates")
    rec.raises("no predictors", lambda: fit_multivariate_relation(a, []), match="predictor")


@case("analyze_correlation", "refusals")
def guard_correlation(rec):
    sst = woa_sst()
    rec.raises("area weighting with spearman", lambda: analyze_correlation(sst, sst, dim=["lat", "lon"], method="spearman",
                                                                          weighting="area"), match="Pearson")
    rec.raises("a lag as long as the series", lambda: analyze_correlation(series(12), series(12), lag=12), match="exceeds")
    rec.raises("a lag on a spatial correlation", lambda: analyze_correlation(sst, sst, dim=["lat", "lon"], lag=1),
               match="time")


@case("detect_outliers", "refusals")
def guard_outliers(rec):
    rec.raises("percentile threshold 100", lambda: detect_outliers(series(), method="percentile", threshold=100),
               match="0 < threshold < 100")
    rec.raises("a negative z-score threshold", lambda: detect_outliers(series(), threshold=-1), match="positive")
    flat = detect_outliers(series(50, np.ones(50)))
    rec.metric("no false alarm on a constant series", int(flat.sum()), threshold=0, op="==", tier="safeguard")


@case("compute_difference_metrics", "refusals")
def guard_difference(rec):
    a = series(12)
    b = xr.DataArray(np.arange(12.0), dims="time", coords={"time": monthly_axis("2010-01", 12)})
    rec.raises("no overlap", lambda: compute_difference_metrics(a, b), match="nothing in common")
    bare = xr.DataArray(woa_sst().values, dims=("lat", "lon"), coords=woa_sst().coords)
    rec.raises("area weighting forced on a field without provenance",
               lambda: compute_difference_metrics(bare, bare, weighting="area"), match="grid")


@case("detect_threshold_events", "refusals")
def guard_events(rec):
    rec.raises("a 2-D field", lambda: detect_threshold_events(woa_sst(), 20.0), match="1-D series")
    rec.raises("min_duration 0", lambda: detect_threshold_events(series(), 1.0, min_duration=0), match="at least 1")


@case("compute_exceedance_fraction", "refusals")
def guard_exceedance(rec):
    rec.raises("a negative threshold", lambda: compute_exceedance_fraction(woa_sst(), threshold=-1.0), match="negative")
    rec.raises("an unknown comparison", lambda: compute_exceedance_fraction(woa_sst(), 1.0, comparison="between"),
               match="comparison")


@case("locate_extreme", "refusals")
def guard_extreme(rec):
    rec.raises("a search leaving a dimension behind", lambda: locate_extreme(ersst_short(), dim=["lat", "lon"]),
               match="leave")
    rec.raises("an all-missing field", lambda: locate_extreme(woa_sst() * np.nan), match="missing")


# ---------------------------------------------------------------------------
# Scientific calculators
# ---------------------------------------------------------------------------

@case("compute_anomaly_vs_woa23", "refusals")
def guard_anomaly_woa(rec):
    o2 = load_regular_data("cmip6_historical", "o2", model="GFDL-ESM4", time_range=["2000-01", "2000-02"],
                           depth_range=[0, 10])
    rec.raises("CMIP6 o2 in mol m-3 against WOA23 umol/kg", lambda: compute_anomaly_vs_woa23(o2, var="o2"),
               match="Unit mismatch")
    rec.raises("ERSSTv5 on its 2-degree grid", lambda: compute_anomaly_vs_woa23(ersst_short(), var="temp"),
               match="1-degree grid")
    rec.raises("a curvilinear model", lambda: compute_anomaly_vs_woa23(cesm_tos(), var="temp"), match="curvilinear")
    sst = load_regular_data("cmip6_historical", "tos", model="GFDL-ESM4", time_range=["2000-01", "2000-02"])
    rec.raises("a percentage anomaly of temperature", lambda: compute_anomaly_vs_woa23(sst, var="temp", method="ratio"),
               match="true zero")


@case("compute_anomaly_vs_own_climatology", "refusals")
def guard_anomaly_own(rec):
    rec.raises("a base period shorter than a year", lambda: compute_anomaly_vs_own_climatology(
        series(48), base_period=["2000-01", "2000-06"]), match="Widen")
    rec.raises("a base period outside the series", lambda: compute_anomaly_vs_own_climatology(
        series(48), base_period=["1950-01", "1960-12"]), match="contains none")
    rec.raises("no time axis", lambda: compute_anomaly_vs_own_climatology(woa_sst()), match="time axis")
    degc = series(48)
    degc.attrs["units"] = "degC"
    rec.raises("a percentage anomaly of degC", lambda: compute_anomaly_vs_own_climatology(degc, method="ratio"),
               match="true zero")


@case("convert_concentration_to_content", "refusals")
def guard_content(rec):
    rec.raises("a temperature field", lambda: convert_concentration_to_content(woa_profile("temperature")), match="units")
    o2 = woa_profile("oxygen")
    moved = regrid_vertical(o2, target=[12.0, 55.0, 230.0])
    rec.raises("a field interpolated off its levels", lambda: convert_concentration_to_content(moved), match="depth")


@case("compute_vertical_gradient", "refusals")
def guard_gradient(rec):
    moved = regrid_vertical(woa_profile(), target=[12.0, 55.0, 230.0, 700.0])
    rec.raises("levels interpolated off the dataset's grid", lambda: compute_vertical_gradient(moved), match="do not match")
    rec.raises("no depth dimension", lambda: compute_vertical_gradient(woa_sst()), match="depth")


@case("compute_ratio", "refusals")
def guard_ratio(rec):
    a = xr.DataArray([1.0], dims="n", attrs={"units": "umol kg-1"})
    b = xr.DataArray([1.0], dims="n", attrs={"units": "mol m-3"})
    rec.raises("a per-mass over a per-volume concentration", lambda: compute_ratio(a, b), match="molar ratio")


@case("identify_oxygen_minimum_zone", "refusals")
def guard_omz(rec):
    o2 = load_regular_data("cmip6_historical", "o2", model="GFDL-ESM4", time_range=["2000-01", "2000-01"],
                           lat_range=[0, 3], lon_range=[200, 203])
    rec.raises("oxygen in mol m-3 with an umol/kg threshold", lambda: identify_oxygen_minimum_zone(o2), match="per-volume")


@case("compute_potential_density_anomaly", "refusals")
def guard_density(rec):
    t, s = woa_profile("temperature"), woa_profile("salinity")
    kelvin = t + 273.15
    kelvin.attrs = dict(t.attrs, units="K")
    rec.raises("temperature in kelvin", lambda: compute_potential_density_anomaly(kelvin, s), match="Celsius")
    absolute = s.copy()
    absolute.attrs = dict(s.attrs, units="g/kg")
    rec.raises("absolute salinity in g/kg", lambda: compute_potential_density_anomaly(t, absolute), match="practical")
    rec.raises("mismatched dimensions", lambda: compute_potential_density_anomaly(t, s.isel(lat=0)), match="same dimensions")


@case("compute_mixed_layer_depth", "refusals")
def guard_mld(rec):
    t, s = woa_profile("temperature"), woa_profile("salinity")
    rec.raises("the method passed where ref_depth goes", lambda: compute_mixed_layer_depth(t, s, "density"),
               exc=TypeError, match="ref_depth")
    rec.raises("a reference depth below the profile", lambda: compute_mixed_layer_depth(t, s, ref_depth=9000.0),
               match="outside")
    rec.raises("an unknown method", lambda: compute_mixed_layer_depth(t, s, method="potential"), match="Method")


@case("compute_nutricline_depth", "refusals")
def guard_nutricline(rec):
    rec.raises("no depth dimension", lambda: compute_nutricline_depth(woa_sst()), match="depth")


@case("compute_enso_nino_index", "refusals and the model climatology rule")
def guard_nino(rec):
    sst = ersst_short()
    rec.raises("WOA23 climatology with a base_period", lambda: compute_enso_nino_index(
        sst, climatology="woa23", base_period=["2000-01", "2000-12"]), match="base_period")
    rec.raises("WOA23 climatology on ERSST's 2-degree grid", lambda: compute_enso_nino_index(sst, climatology="woa23"),
               match="1-degree")
    rec.raises("an unknown region", lambda: compute_enso_nino_index(sst, nino_region="nino5"), match="nino_region")
    gfdl = load_regular_data("cmip6_historical", "tos", model="GFDL-ESM4", time_range=["2000-01", "2001-12"],
                             lat_range=[-10, 10], lon_range=[180, 260])
    index = compute_enso_nino_index(gfdl)
    rec.metric("GFDL-ESM4 on the WOA23 grid still uses its own climatology under 'auto'",
               "own" in index.attrs["climatology_baseline"], op="true", tier="safeguard",
               note=index.attrs["climatology_baseline"][:100])


@case("compute_pdo_index", "refusals")
def guard_pdo(rec):
    north_pacific = load_regular_data("ersst", "sst", time_range=["1990-01", "1999-12"], lat_range=[15, 75],
                                      lon_range=[100, 270])
    rec.raises("a North Pacific subset instead of the global field", lambda: compute_pdo_index(north_pacific),
               match="GLOBAL")
    rec.raises("a curvilinear field", lambda: compute_pdo_index(cesm_tos()), match="curvilinear")
    rec.raises("no time axis", lambda: compute_pdo_index(woa_sst()), match="time")
