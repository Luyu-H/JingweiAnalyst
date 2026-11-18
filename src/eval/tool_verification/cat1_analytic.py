"""
Category 1: analytic solutions on synthetic fields.
"""
from __future__ import annotations

import functools

import cftime
import numpy as np
import pandas as pd
import xarray as xr

from src.eval.tool_verification import refs
from src.eval.tool_verification.common import (ar1, case, coord_like, lat_lon_like, max_abs, monthly_axis,
                                               refill, rel_err, zone_area)
from src.tools.basic_calculator import (compute_basic_stat, compute_difference_metrics,
                                        compute_exceedance_fraction, detect_outliers,
                                        detect_threshold_events, fit_1D_trend, fit_trend_map,
                                        locate_extreme, temporal_aggregate, transform_unit)
from src.tools.data_loader import load_areas, load_climatology_data, load_regular_data
from src.tools.data_openers import woa23_files
from src.tools.preprocessing import regrid_horizontal, regrid_vertical
from src.tools.scientific_calculator import (compute_anomaly_vs_own_climatology, compute_anomaly_vs_woa23,
                                             compute_enso_nino_index, compute_mixed_layer_depth,
                                             compute_nutricline_depth, compute_pdo_index, compute_ratio,
                                             compute_vertical_gradient, convert_concentration_to_content,
                                             identify_oxygen_minimum_zone)

CAT = "cat1_analytic"
MODELS = ("CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR")


# ---------------------------------------------------------------------------
# Templates (loaded once per process)
# ---------------------------------------------------------------------------

@functools.lru_cache(None)
def woa_surface():
    return load_climatology_data("woa23", "temperature", depth_range=[0, 0]).isel(depth=0).load()


@functools.lru_cache(None)
def ersst_surface():
    return load_regular_data("ersst", "sst", time_range=["2000-01", "2000-01"]).isel(time=0).load()


@functools.lru_cache(None)
def woa_column(variable="temperature", time_slice="annual"):
    """A 2x2 block of whole WOA23 columns, every level present."""
    return load_climatology_data("woa23", variable, time_slice=time_slice,
                                 lat_range=[10, 12], lon_range=[-30, -28]).load()


@functools.lru_cache(None)
def model_column(model, variable="thetao"):
    """One month of whole model columns over a small open-ocean box."""
    return load_regular_data("cmip6_historical", variable, model=model, time_range=["2000-01", "2000-01"],
                             lat_range=[0, 3], lon_range=[200, 203]).isel(time=0).load()


def raw_woa_bounds(variable, time_slice="annual"):
    path = woa23_files(variable, time_slice)["files"][0][1]
    bounds, _ = refs.raw_nc(path, "depth_bnds")
    depth, _ = refs.raw_nc(path, "depth")
    return depth, bounds


# ---------------------------------------------------------------------------
# load_areas: closed-form quadrangle areas on the WGS84 ellipsoid
# ---------------------------------------------------------------------------

@case("load_areas", "closed-form WGS84 band areas, 1-degree and ERSST 2-degree grids")
def areas_closed_form(rec):
    for label, areas, half in (("woa23 1deg", load_areas(dataset="woa23"), 0.5),
                               ("default 1deg", load_areas(), 0.5),
                               ("ersst 2deg", load_areas(dataset="ersst"), 1.0)):
        lat = areas["lat"].values
        lo = np.clip(lat - half, -89.99, 89.99)
        hi = np.clip(lat + half, -89.99, 89.99)
        dlon = float(np.diff(areas["lon"].values[:2])[0])
        expected = zone_area(lo, hi, dlon)
        band = areas.isel(lon=0).values
        # The tool's cells are polygons with geodesic edges, which bulge poleward of
        # the parallels by an amount that grows as dlon^2; 1e-4 per square degree
        # bounds that difference with a factor of ~2 to spare.
        rec.metric(f"{label}: max relative error of a band's cell area", rel_err(band, expected),
                   threshold=1e-4 * dlon ** 2, tier="T1", method="authalic-latitude closed form (areaquad)",
                   note="geodesic vs parallel cell edges; the docstring's 'equivalent to areaquad' holds to this level")
        rec.metric(f"{label}: cells in a band identical", float(np.ptp(areas.values, axis=1).max()),
                   threshold=1e-6, tier="T0")
        total = float(areas.sum())
        covered = float(zone_area(lo.min(), hi.max(), 360.0))
        rec.metric(f"{label}: global sum vs closed form over the covered latitudes",
                   abs(total / covered - 1), threshold=1e-7, tier="T1",
                   note="shared geodesic edges cancel in the sum; only the two outermost edges remain "
                        "(~3e-8 for ERSST's edges at +-89, ~1e-12 at +-89.99)")


def geodesic_bulge(lat1, lat2, dlon_deg):
    """
    Relative area by which a band of cells with geodesic top and bottom edges
    differs from the parallel-bounded band, to leading order on a sphere. A
    geodesic between two points dlon apart on the parallel lat bulges poleward
    by dlon^2/8 sin(lat) cos(lat), which over a full circle encloses
    2 pi R^2 cos^2(lat) |sin(lat)| dlon^2 / 12.
    """
    d = np.radians(dlon_deg)
    edges = [lat for lat in (lat1, lat2) if abs(lat) < 89.999]
    excess = sum(np.cos(np.radians(lat)) ** 2 * abs(np.sin(np.radians(lat))) * d ** 2 / 12 for lat in edges)
    return float(excess / abs(np.sin(np.radians(lat2)) - np.sin(np.radians(lat1))))


# ---------------------------------------------------------------------------
# compute_basic_stat: area- and volume-weighted means of analytic fields
# ---------------------------------------------------------------------------

@case("compute_basic_stat", "area-weighted means of analytic fields on full-cover grids")
def basic_stat_area_analytic(rec):
    # Band edges are placed on cell edges (integers on the 1-degree grid, odd
    # integers on ERSST's 2-degree grid) so the indicator has no partial cells.
    for label, template, bands in (("woa23 1deg", woa_surface(), ((-30, 30), (60, 90), (-45, -10))),
                                   ("ersst 2deg", ersst_surface(), ((-29, 29), (59, 89), (-45, -9)))):
        lat, lon = lat_lon_like(template)
        lat_lo, lat_hi = float(lat.min()) - (0.5 if "woa" in label else 1.0), float(lat.max()) + (0.5 if "woa" in label else 1.0)
        lat_lo, lat_hi = max(lat_lo, -89.99), min(lat_hi, 89.99)
        whole = float(zone_area(lat_lo, lat_hi))
        mean = lambda values: float(compute_basic_stat(refill(template, values), "area_weighted_mean",
                                                       dim=["lat", "lon"]))
        rec.metric(f"{label}: mean of f=1", abs(mean(np.ones_like(lat)) - 1), threshold=1e-12, tier="T0")
        rec.metric(f"{label}: mean of f=sin(lat) (0 by symmetry)", abs(mean(np.sin(np.radians(lat)))),
                   threshold=1e-12, tier="T0")
        dlon = float(np.diff(template["lon"].values[:2])[0])
        for lo, hi in bands:
            lo_c, hi_c = max(lo, lat_lo), min(hi, lat_hi)
            expected = float(zone_area(lo_c, hi_c)) / whole
            got = mean(((lat > lo) & (lat < hi)).astype(float))
            bound = 1.5 * geodesic_bulge(lo_c, hi_c, dlon) + 1e-12
            rec.metric(f"{label}: area fraction of the band ({lo:g}, {hi:g})", abs(got / expected - 1),
                       threshold=bound, tier="T1", expected=expected,
                       note=f"tool {got:.10f}; threshold = 1.5 x the geodesic-edge bulge at the band's edges")

    # A window crossing the dateline comes back on a continuous 160..200 axis; its
    # two halves have equal areas, so the eastern half's share is exactly 1/2.
    box = load_climatology_data("woa23", "temperature", depth_range=[0, 0], lat_range=[-10, 10],
                                lon_range=[160, -160]).isel(depth=0)
    _, lon = lat_lon_like(box)
    got = float(compute_basic_stat(refill(box, (lon > 180).astype(float)), "area_weighted_mean",
                                   dim=["lat", "lon"]))
    rec.metric("dateline window: share of the half east of 180", abs(got - 0.5), threshold=1e-12, tier="T0",
               note=f"lon axis {float(lon.min()):g}..{float(lon.max()):g}")


@case("compute_basic_stat", "volume-weighted mean of a depth indicator on WOA23 columns")
def basic_stat_volume_analytic(rec):
    template = load_climatology_data("woa23", "temperature", depth_range=[0, 1000],
                                     lat_range=[10, 14], lon_range=[-30, -26]).load()
    depth_raw, bounds = raw_woa_bounds("temperature")
    keep = depth_raw <= 1000
    dz = (bounds[:, 1] - bounds[:, 0])[keep]
    depth = coord_like(template, "depth")
    rec.metric("mean of f=1", abs(float(compute_basic_stat(refill(template, 1.0), "volume_weighted_mean")) - 1),
               threshold=1e-12, tier="T0")
    for cut in (100.0, 500.0):
        got = float(compute_basic_stat(refill(template, (depth <= cut).astype(float)), "volume_weighted_mean"))
        expected = dz[depth_raw[keep] <= cut].sum() / dz.sum()
        rec.metric(f"volume fraction above {cut:g} m", abs(got - expected), threshold=1e-12, tier="T0",
                   expected=expected, method="raw WOA23 depth_bnds read with netCDF4",
                   note="areas cancel on a full-cover block, leaving the layer thicknesses")


@case("compute_basic_stat", "argmax/argmin return the coordinate label of the extreme")
def basic_stat_argmax(rec):
    time = monthly_axis("2001-01", 36)
    values = np.sin(np.arange(36) / 3.0)
    series = xr.DataArray(values, coords={"time": time}, dims="time")
    got = compute_basic_stat(series, "argmax", dim="time")
    rec.metric("argmax label equals the numpy argmax label",
               pd.Timestamp(got.values) == time[int(np.argmax(values))], op="true", tier="T0")
    got = compute_basic_stat(series, "argmin", dim="time")
    rec.metric("argmin label equals the numpy argmin label",
               pd.Timestamp(got.values) == time[int(np.argmin(values))], op="true", tier="T0")


# ---------------------------------------------------------------------------
# compute_exceedance_fraction
# ---------------------------------------------------------------------------

@case("compute_exceedance_fraction", "area and volume fractions of analytic fields")
def exceedance_analytic(rec):
    template = woa_surface()
    lat, _ = lat_lon_like(template)
    field = refill(template, lat, units="degrees_north")
    whole = float(zone_area(-89.99, 89.99))
    result = compute_exceedance_fraction(field, threshold=30.0, comparison="abs_less")
    expected = float(zone_area(-30, 30)) / whole
    rec.metric("abs_less 30 on f=lat: area fraction of |lat|<30", abs(result["fraction"] / expected - 1),
               threshold=1e-4, tier="T1", expected=expected)
    rec.metric("auto weighting chose area", result["weighting"].startswith("area"), op="true", tier="T0",
               note=result["weighting"][:120])
    result = compute_exceedance_fraction(field, threshold=0.0, comparison="greater")
    rec.metric("greater 0 on f=lat: northern half", abs(result["fraction"] - 0.5), threshold=1e-12, tier="T0")
    result = compute_exceedance_fraction(field, threshold=20.0, comparison="abs_greater", reference=10.0)
    expected = 1 - float(zone_area(-10, 30)) / whole
    rec.metric("abs_greater 20 about reference 10: outside (-10, 30)", abs(result["fraction"] / expected - 1),
               threshold=1e-4, tier="T1")

    constant = refill(template, 2.0)
    rec.metric("strict inequality: x=2 is not > 2",
               compute_exceedance_fraction(constant, threshold=2.0, comparison="greater")["fraction"],
               threshold=0.0, op="==", tier="T0")
    rec.metric("strict inequality: |x-2| < 0 holds nowhere",
               compute_exceedance_fraction(constant, threshold=0.0, comparison="abs_less", reference=2.0)["fraction"],
               threshold=0.0, op="==", tier="T0")

    column = load_climatology_data("woa23", "temperature", depth_range=[0, 1000],
                                   lat_range=[10, 14], lon_range=[-30, -26]).load()
    depth_raw, bounds = raw_woa_bounds("temperature")
    keep = depth_raw <= 1000
    dz = (bounds[:, 1] - bounds[:, 0])[keep]
    result = compute_exceedance_fraction(refill(column, coord_like(column, "depth"), units="m"),
                                         threshold=100.0, comparison="less")
    expected = dz[depth_raw[keep] < 100].sum() / dz.sum()
    rec.metric("volume fraction with depth < 100 m (strict, the 100 m layer excluded)",
               abs(result["fraction"] - expected), threshold=1e-12, tier="T0", expected=expected)
    rec.metric("auto weighting chose volume when depth is reduced", result["weighting"].startswith("volume"),
               op="true", tier="T0")


# ---------------------------------------------------------------------------
# temporal_aggregate: day weighting under each calendar
# ---------------------------------------------------------------------------

CALENDARS = ("noleap", "365_day", "all_leap", "366_day", "360_day", "standard", "proleptic_gregorian",
             "julian")


def _cftime_series(calendar, values_of_month, years=(1999, 2000, 2001, 2002)):
    stamps = [cftime.datetime(y, m, 15, calendar=calendar) for y in years for m in range(1, 13)]
    values = np.array([values_of_month(y, m) for y in years for m in range(1, 13)], dtype="float64")
    time = xr.DataArray(xr.CFTimeIndex(stamps), dims="time")
    return xr.DataArray(values, coords={"time": time}, dims="time", name="f")


@case("temporal_aggregate", "day-weighted means under eight CF calendars (cftime axis)")
def temporal_calendars(rec):
    years = (1999, 2000, 2001, 2002)
    for calendar in CALENDARS:
        for label, func in (("f=1[Feb]", lambda y, m: float(m == 2)), ("f=month", lambda y, m: float(m))):
            series = _cftime_series(calendar, func, years)
            try:
                annual = temporal_aggregate(series, freq="year")
            except Exception as exc:
                rec.metric(f"{calendar} {label}: freq='year' runs", False, op="true", tier="T0",
                           note=f"{type(exc).__name__}: {exc}"[:300])
                continue
            expected = []
            for y in years:
                d = np.array([refs.days_in_month(y, m, calendar) for m in range(1, 13)], dtype=float)
                v = np.array([func(y, m) for m in range(1, 13)])
                expected.append(np.sum(v * d) / d.sum())
            rec.metric(f"{calendar} {label}: annual means vs calendar rules", max_abs(annual.values, expected),
                       threshold=1e-12, tier="T0", method="month lengths written out from the CF calendar rules",
                       note=f"year 2000: tool {annual.values[1]:.10f}, expected {expected[1]:.10f}")
        series = _cftime_series(calendar, lambda y, m: float(m), years)
        whole = temporal_aggregate(series, freq="mean")
        d = np.array([refs.days_in_month(y, m, calendar) for y in years for m in range(1, 13)], dtype=float)
        v = np.array([m for y in years for m in range(1, 13)], dtype=float)
        rec.metric(f"{calendar} f=month: freq='mean' vs calendar rules", abs(float(whole) - np.sum(v * d) / d.sum()),
                   threshold=1e-12, tier="T0")
        seasons = temporal_aggregate(series, freq="season")
        djf = seasons.where((seasons["season"] == "DJF") & (seasons["year"] == 2000), drop=True)
        dd = [refs.days_in_month(1999, 12, calendar), refs.days_in_month(2000, 1, calendar),
              refs.days_in_month(2000, 2, calendar)]
        expected = (12 * dd[0] + 1 * dd[1] + 2 * dd[2]) / sum(dd)
        rec.metric(f"{calendar}: DJF 2000 (Dec 1999-Feb 2000) day-weighted, labelled 2000",
                   abs(float(djf.values[0]) - expected) if djf.size == 1 else f"{djf.size} DJF-2000 entries",
                   threshold=1e-12, tier="T0")


@case("temporal_aggregate", "day weighting end to end: model calendar after the loader")
def temporal_calendar_end_to_end(rec):
    for model, calendar in (("CESM2", "noleap"), ("GFDL-ESM4", "noleap"), ("MPI-ESM1-2-LR", "proleptic_gregorian")):
        tos = load_regular_data("cmip6_historical", "tos", model=model, time_range=["2000-01", "2000-12"],
                                lat_range=[40, 50], lon_range=[150, 200])
        series = compute_basic_stat(tos, "area_weighted_mean", dim=["lat", "lon"]).load()
        annual = float(temporal_aggregate(series, freq="year").sel(year=2000))
        own = np.array([refs.days_in_month(2000, m, calendar) for m in range(1, 13)], dtype=float)
        greg = np.array([refs.days_in_month(2000, m, "standard") for m in range(1, 13)], dtype=float)
        x = series.values
        expected = np.sum(x * own) / own.sum()
        rec.metric(f"{model} ({calendar}): 2000 annual mean vs its own calendar's month lengths",
                   abs(annual - expected), threshold=1e-9, tier="T0",
                   note=f"tool {annual:.6f}; own-calendar {expected:.6f}; Gregorian weights would give "
                        f"{np.sum(x * greg) / greg.sum():.6f}; time dtype after loading: {series['time'].dtype}")


# ---------------------------------------------------------------------------
# Trends
# ---------------------------------------------------------------------------

@case("fit_1D_trend", "known slopes on year, date and depth axes")
def trend_known_slope(rec):
    years = np.arange(1950, 2020)
    series = xr.DataArray(3.0 + 0.02 * (years - 1950), coords={"year": years}, dims="year", attrs={"units": "degC"})
    result = fit_1D_trend(series, method="linear")
    rec.metric("integer year axis: slope 0.02 per year", rel_err(result["trend"]["slope"], 0.02),
               threshold=1e-10, tier="T0", note=result["trend"]["units"])
    result = fit_1D_trend(xr.DataArray(1.0 + 0.5 * (years - 1950) + 0.003 * (years - 1950) ** 2,
                                       coords={"year": years}, dims="year"), method="quadratic")
    rec.metric("quadratic: coefficients (0.003, 0.5, 1)", rel_err([result["params"][k] for k in "abc"],
                                                                 [0.003, 0.5, 1.0]), threshold=1e-8, tier="T0")
    result = fit_1D_trend(xr.DataArray(2.0 * np.exp(0.01 * (years - 1950)), coords={"year": years}, dims="year"),
                          method="exponential")
    rec.metric("exponential: (a, b) = (2, 0.01)", rel_err([result["params"]["a"], result["params"]["b"]],
                                                          [2.0, 0.01]), threshold=1e-6, tier="T0")

    time = monthly_axis("1980-01", 480)
    series = xr.DataArray(0.03 * np.arange(480) / 12.0, coords={"time": time}, dims="time", attrs={"units": "degC"})
    result = fit_1D_trend(series, method="linear")
    rec.metric("monthly date axis: 0.03 per year (month = 1/12 year)", rel_err(result["trend"]["slope"], 0.03),
               threshold=1e-3, tier="T1", note=f"x_unit={result['x_unit']}; exact months vs 365.2425-day years")
    rec.metric("monthly date axis: slope reported per year", result["x_unit"] == "year", op="true", tier="T0")

    depth = np.array([0., 10, 20, 50, 100, 200, 500, 1000])
    series = xr.DataArray(5 + 0.01 * depth, coords={"depth": ("depth", depth, {"units": "m"})}, dims="depth")
    result = fit_1D_trend(series, method="linear")
    rec.metric("depth axis: 0.01 per metre", rel_err(result["trend"]["slope"], 0.01), threshold=1e-10, tier="T0",
               note=result["trend"]["units"])

    rng = np.random.default_rng(1)
    noisy = xr.DataArray(0.02 * np.arange(480) / 12.0 + ar1(rng, 480, 0.5) * 0.2,
                         coords={"time": time}, dims="time")
    result = fit_1D_trend(noisy, method="linear")
    t = result["trend"]
    rec.metric("AR(1) noise: |slope - 0.02| in standard errors", abs(t["slope"] - 0.02) / t["slope_se"],
               threshold=4.0, tier="T1", note=f"slope {t['slope']:.5f} [{t['ci95_low']:.5f}, {t['ci95_high']:.5f}]")


@case("fit_trend_map", "known slope field on year and date axes")
def trend_map_known(rec):
    lat = np.linspace(-20, 20, 5)
    lon = np.linspace(0, 50, 6)
    a = 0.01 + 0.002 * lat[:, None] + 0.0005 * lon[None, :]
    years = np.arange(1960, 2020)
    field = xr.DataArray(1.0 + a[None] * (years - 1960)[:, None, None], dims=("year", "lat", "lon"),
                         coords={"year": years, "lat": lat, "lon": lon}, attrs={"units": "degC"})
    result = fit_trend_map(field, dim="year")
    rec.metric("year axis, per_decade=True: slope = 10 a", rel_err(result["slope"].values, 10 * a),
               threshold=1e-10, tier="T0", note=result["units"])
    result = fit_trend_map(field, dim="year", per_decade=False)
    rec.metric("year axis, per_decade=False: slope = a", rel_err(result["slope"].values, a), threshold=1e-10,
               tier="T0")
    time = monthly_axis("1980-01", 480)
    field = xr.DataArray(a[None] * (np.arange(480) / 12.0)[:, None, None], dims=("time", "lat", "lon"),
                         coords={"time": time, "lat": lat, "lon": lon})
    result = fit_trend_map(field)
    rec.metric("monthly date axis: slope per decade = 10 a", rel_err(result["slope"].values, 10 * a),
               threshold=1e-3, tier="T1")
    gappy = field.copy()
    gappy[::7, 1, 1] = np.nan
    result = fit_trend_map(gappy)
    rec.metric("a cell with gaps still recovers its slope", rel_err(float(result["slope"][1, 1]), 10 * a[1, 1]),
               threshold=1e-3, tier="T1")


# ---------------------------------------------------------------------------
# Vertical operators
# ---------------------------------------------------------------------------

def _column_templates():
    yield "woa23 temperature", woa_column("temperature")
    yield "glodap", load_climatology_data("glodap", "temperature", lat_range=[10, 12],
                                          lon_range=[-30, -28]).load()
    for model in MODELS:
        yield model, model_column(model)


@case("compute_vertical_gradient", "linear profiles are exact at every level, end layers included")
def vertical_gradient_linear(rec):
    for label, template in _column_templates():
        z = coord_like(template, "depth")
        grad = compute_vertical_gradient(refill(template, 3.0 + 0.01 * z))
        rec.metric(f"{label}: c = 3 + 0.01 z gives 0.01 everywhere", rel_err(grad.values, 0.01), threshold=1e-9,
                   tier="T0", note=f"{template.sizes['depth']} levels; regression test for the end layers")


@case("compute_vertical_gradient", "uniform 5 m levels reduce to the central difference")
def vertical_gradient_uniform(rec):
    template = woa_column("temperature")
    z = template["depth"].values
    values = np.sin(z / 37.0) + 1e-4 * z ** 2
    grad = compute_vertical_gradient(refill(template, values.reshape(-1, 1, 1))).isel(lat=0, lon=0).values
    interior = slice(1, 20)                       # 5 m spacing, 0-100 m
    central = (values[2:21] - values[0:19]) / (z[2:21] - z[0:19])
    rec.metric("levels 5-95 m equal (c[k+1]-c[k-1]) / 10", max_abs(grad[interior], central), threshold=1e-12,
               tier="T0", method="np.gradient's interior stencil")
    exact = np.cos(z / 37.0) / 37.0 + 2e-4 * z
    deep = z > 500
    rec.metric("smooth profile, levels below 500 m: max relative error vs dc/dz", rel_err(grad[deep], exact[deep]),
               op="info", tier="T1", note="report only: coarse 50-100 m levels")


@case("regrid_vertical", "linear profiles exact, quadratic profiles match the interpolation error formula")
def regrid_vertical_analytic(rec):
    template = woa_column("temperature")
    z = template["depth"].values
    targets = [3.3, 47.1, 333.3, 1234.5, 4999.0]
    result = regrid_vertical(refill(template, (2 + 0.5 * z).reshape(-1, 1, 1)), target=targets)
    rec.metric("linear profile reproduced at arbitrary depths",
               max_abs(result.isel(lat=0, lon=0).values, 2 + 0.5 * np.array(targets)), threshold=1e-9, tier="T0")
    result = regrid_vertical(refill(template, (z ** 2).reshape(-1, 1, 1)), target=targets).isel(lat=0, lon=0).values
    k = np.searchsorted(z, targets) - 1
    zt = np.array(targets)
    expected = zt ** 2 + (zt - z[k]) * (z[k + 1] - zt)
    rec.metric("quadratic: equals z^2 + (z - z_k)(z_k+1 - z)", rel_err(result, expected), threshold=1e-12,
               tier="T0")
    gfdl = model_column("GFDL-ESM4")
    zg = gfdl["depth"].values
    result = regrid_vertical(refill(gfdl, coord_like(gfdl, "depth")), target=[1.0, zg[0], 30.0], fill_surface=True)
    column = result.isel({gfdl.attrs.get("y_dim", "lat"): 0, gfdl.attrs.get("x_dim", "lon"): 0}).values
    rec.metric("fill_surface: 1 m above GFDL's 2.5 m level takes that level's value", abs(column[0] - zg[0]),
               threshold=1e-12, tier="T0")
    shallow = refill(template, (2 + 0.5 * z).reshape(-1, 1, 1))
    shallow = shallow.where(shallow["depth"] <= 200)
    result = regrid_vertical(shallow, target=[150.0, 210.0]).isel(lat=0, lon=0).values
    rec.metric("no extrapolation below the deepest valid level", bool(np.isfinite(result[0]) and np.isnan(result[1])),
               op="true", tier="T0")


@case("compute_nutricline_depth", "tanh profiles: maximum gradient at the inflection")
def nutricline_tanh(rec):
    template = woa_column("nitrate")
    z = coord_like(template, "depth")
    got = compute_nutricline_depth(refill(template, 15 + 15 * np.tanh((z - 150.0) / 40.0)))
    rec.metric("inflection on a level (150 m) is returned exactly", max_abs(got.values, 150.0), threshold=0.0,
               op="<=", tier="T0")
    got = compute_nutricline_depth(refill(template, 15 + 15 * np.tanh((z - 160.0) / 40.0)))
    rec.metric("inflection between levels (160 m): within half the 25 m spacing", max_abs(got.values, 160.0),
               threshold=12.5, tier="T1", note=f"returned {np.unique(got.values)} (level centres, not interpolated)")
    got = compute_nutricline_depth(refill(template, 15 + 15 * np.tanh((z - 150.0) / 40.0)), min_depth=200.0)
    rec.metric("min_depth=200 keeps the answer at or below 200 m", float(np.nanmin(got.values)), threshold=200.0,
               op=">=", tier="T0")
    got = compute_nutricline_depth(refill(template, 15 - 15 * np.tanh((z - 150.0) / 40.0)))
    rec.metric("decreasing profile (no positive gradient) is missing", bool(np.all(np.isnan(got.values))),
               op="true", tier="T0")


def _mld_truth(lat, lon, h, ref, kind, criterion, threshold, t0=20.0, gamma=0.01, s0=35.0):
    """The crossing depth found on a 1 cm grid with the same TEOS-10 route the tool takes."""
    import gsw
    zf = np.arange(0.0, h + 400.0, 0.01)
    tf = np.where(zf <= h, t0, t0 - gamma * (zf - h))
    if criterion == "temperature":
        delta = np.abs(tf - t0)
    else:
        p = gsw.p_from_z(-zf, lat)
        sa = gsw.SA_from_SP(np.full_like(zf, s0), p, lon, lat)
        ct = gsw.CT_from_pt(sa, tf) if kind == "potential" else gsw.CT_from_t(sa, tf, p)
        sigma = gsw.rho(sa, ct, 0.0)
        p_ref = gsw.p_from_z(-ref, lat)
        sa_ref = gsw.SA_from_SP(s0, p_ref, lon, lat)
        ct_ref = gsw.CT_from_pt(sa_ref, t0) if kind == "potential" else gsw.CT_from_t(sa_ref, t0, p_ref)
        delta = sigma - gsw.rho(sa_ref, ct_ref, 0.0)
    k = int(np.argmax(delta >= threshold))
    return zf[k - 1] + (threshold - delta[k - 1]) / (delta[k] - delta[k - 1]) * (zf[k] - zf[k - 1])


@case("compute_mixed_layer_depth", "idealised mixed layer over a linear thermocline on four vertical grids")
def mld_analytic(rec):
    pairs = [("woa23", woa_column("temperature"), woa_column("salinity"), "in_situ")]
    for model in MODELS:
        pairs.append((model, model_column(model, "thetao"), model_column(model, "so"), "potential"))
    for label, temp, sal, kind in pairs:
        z = temp["depth"].values
        h = float(z[z >= 40][0])
        zz = coord_like(temp, "depth")
        t = np.where(zz <= h, 20.0, 20.0 - 0.01 * (zz - h))
        temp_f, sal_f = refill(temp, t), refill(sal, 35.0)
        lat, lon = lat_lon_like(temp)
        mld_t = compute_mixed_layer_depth(temp_f, sal_f, method="temperature")
        ref_used = float(str(mld_t.attrs["reference_depth"]).split()[0])
        rec.metric(f"{label}: temperature criterion gives h + 0.2/0.01 = {h + 20:g} m", max_abs(mld_t.values, h + 20),
                   threshold=1e-9, tier="T0", note=f"reference level used: {ref_used:g} m")
        mld_d = compute_mixed_layer_depth(temp_f, sal_f, method="density")
        flat_lat, flat_lon = lat.reshape(temp.sizes["depth"], -1)[0], lon.reshape(temp.sizes["depth"], -1)[0]
        truth = np.array([_mld_truth(a, b, h, ref_used, kind, "density", 0.03) for a, b in zip(flat_lat, flat_lon)])
        rec.metric(f"{label}: density criterion vs 1 cm TEOS-10 truth", max_abs(mld_d.values.ravel(), truth),
                   threshold=0.05, tier="T1", note=f"truth {truth.min():.3f}-{truth.max():.3f} m")
        both = compute_mixed_layer_depth(temp_f, sal_f, method="both")
        rec.metric(f"{label}: 'both' is the shallower of the two",
                   max_abs(both.values, np.minimum(mld_t.values, mld_d.values)), threshold=1e-12, tier="T0")
        # A column mixed to the bottom has uniform SA and CT. Uniform in-situ
        # temperature is NOT that (its potential temperature falls with depth), so
        # the column is built backwards from constant SA and CT with TEOS-10.
        import gsw
        p = gsw.p_from_z(-zz, lat)
        sp = gsw.SP_from_SA(35.16504, p, lon, lat)
        t_mixed = gsw.t_from_CT(35.16504, 20.0, p) if kind == "in_situ" else gsw.pt_from_CT(35.16504, 20.0)
        never = compute_mixed_layer_depth(refill(temp, t_mixed), refill(sal, sp), method="density")
        rec.metric(f"{label}: a column mixed to the bottom (uniform SA, CT) is missing",
                   bool(np.all(np.isnan(never.values))), op="true", tier="T0")


@case("identify_oxygen_minimum_zone", "thickness of a prescribed low-oxygen band on WOA23 levels")
def omz_analytic(rec):
    template = woa_column("oxygen")
    depth_raw, bounds = raw_woa_bounds("oxygen")
    z = coord_like(template, "depth")
    band = (z >= 200) & (z <= 600)
    thickness, mask = identify_oxygen_minimum_zone(refill(template, np.where(band, 10.0, 150.0)), threshold=20.0)
    in_band = (depth_raw >= 200) & (depth_raw <= 600)
    expected = float(np.sum(bounds[in_band, 1] - bounds[in_band, 0]))
    rec.metric("thickness = sum of the band's layer thicknesses", max_abs(thickness.values, expected),
               threshold=1e-9, tier="T0", expected=expected, method="raw WOA23 depth_bnds")
    values = np.where(band, 10.0, 150.0)
    values = np.where(z == 300.0, 20.0, values)
    thickness, _ = identify_oxygen_minimum_zone(refill(template, values), threshold=20.0)
    k = int(np.where(depth_raw == 300.0)[0][0])
    rec.metric("a layer exactly at the threshold is not counted (strict <)",
               max_abs(thickness.values, expected - (bounds[k, 1] - bounds[k, 0])), threshold=1e-9, tier="T0")
    land = refill(template, np.where(band, 10.0, 150.0))
    land[:, 0, 0] = np.nan
    thickness, _ = identify_oxygen_minimum_zone(land, threshold=20.0)
    rec.metric("an all-missing column is NaN, not 0 m", bool(np.isnan(thickness.values[0, 0])), op="true", tier="T0")


@case("convert_concentration_to_content", "constant concentration over a WOA23 box")
def content_analytic(rec):
    template = load_climatology_data("woa23", "oxygen", lat_range=[10, 20], lon_range=[-30, -20],
                                     depth_range=[0, 2000]).load()
    depth_raw, bounds = raw_woa_bounds("oxygen")
    keep = depth_raw <= 2000
    dz = float(np.sum(bounds[keep, 1] - bounds[keep, 0]))
    area = float(zone_area(10.0, 20.0, 10.0))
    content = convert_concentration_to_content(refill(template, 200.0))
    expected = 200.0e-6 * 1025.0 * area * dz
    rec.metric("total moles = c * rho * area * depth", abs(float(content.sum()) / expected - 1), threshold=1e-4,
               tier="T1", expected=expected, method="closed-form zone area x raw depth_bnds",
               note="1e-4 absorbs the geodesic-edge difference in the tool's areas")
    rec.metric("mass-basis input recorded", content.attrs.get("molar_basis_of_input") == "mass", op="true", tier="T0")


# ---------------------------------------------------------------------------
# Anomalies and indices
# ---------------------------------------------------------------------------

@case("compute_anomaly_vs_own_climatology", "seasonal cycle removed exactly, base-period means zero")
def own_anomaly_analytic(rec):
    time = monthly_axis("1981-01", 480)
    m = time.month.values
    t = np.arange(480) / 12.0
    rng = np.random.default_rng(2)
    seasonal = 3 * np.sin(2 * np.pi * (m - 1) / 12)[:, None] * np.ones((1, 4))
    field = xr.DataArray(seasonal, coords={"time": time}, dims=("time", "x"))
    anomaly = compute_anomaly_vs_own_climatology(field)
    rec.metric("pure seasonal cycle -> anomaly 0", float(np.abs(anomaly).max()), threshold=1e-12, tier="T0")
    noisy = field + 0.02 * t[:, None] + rng.standard_normal((480, 4))
    anomaly = compute_anomaly_vs_own_climatology(noisy, base_period=["1991-01", "2020-12"])
    base = anomaly.sel(time=slice("1991-01", "2020-12"))
    rec.metric("base-period monthly means of the anomaly are 0",
               float(np.abs(base.groupby("time.month").mean()).max()), threshold=1e-12, tier="T0")
    clim = noisy.sel(time=slice("1991-01", "2020-12")).groupby("time.month").mean()
    expected = noisy.values - clim.sel(month=xr.DataArray(m, dims="time")).values
    rec.metric("steps outside the base period also have the base climatology removed",
               max_abs(anomaly.values, expected), threshold=1e-12, tier="T0")
    anomaly = compute_anomaly_vs_own_climatology(noisy, freq="year")
    rec.metric("freq='year': anomaly has zero mean", float(np.abs(anomaly.mean("time")).max()), threshold=1e-12,
               tier="T0")


@case("compute_anomaly_vs_woa23", "a tiled WOA23 climatology has zero anomaly, plus c gives c")
def woa_anomaly_analytic(rec):
    for freq, slice_name, group in (("month", "monthly", "month"), ("season", "seasonal", "season"),
                                    ("year", "annual", None)):
        clim = load_climatology_data("woa23", "temperature", time_slice=slice_name, depth_range=[0, 20],
                                     lat_range=[-5, 5], lon_range=[170, -170]).load()
        time = monthly_axis("2000-01", 24)
        if group == "month":
            values = clim.sel(month=xr.DataArray(time.month.values, dims="time")).values
        elif group == "season":
            values = clim.isel(season=xr.DataArray((time.month.values - 1) // 3, dims="time")).values
        else:
            values = np.broadcast_to(clim.values, (24,) + clim.shape)
        dims = ("time",) + tuple(d for d in clim.dims if d not in ("month", "season"))
        coords = {d: clim[d] for d in dims if d != "time"}
        coords["time"] = time
        data = xr.DataArray(values, dims=dims, coords=coords, attrs=dict(clim.attrs), name="temperature")
        anomaly = compute_anomaly_vs_woa23(data, var="temp", freq=freq)
        rec.metric(f"freq={freq}: tiled climatology -> 0", float(np.nanmax(np.abs(anomaly.values))),
                   threshold=1e-6, tier="T0", note="WOA23 is float32")
        anomaly = compute_anomaly_vs_woa23(data.copy(data=data.values + 1.5), var="temp", freq=freq)
        rec.metric(f"freq={freq}: +1.5 -> 1.5", float(np.nanmax(np.abs(anomaly.values - 1.5))), threshold=1e-5,
                   tier="T0")
        rec.metric(f"freq={freq}: dateline window keeps its continuous lon labels",
                   bool(np.array_equal(anomaly["lon"].values, data["lon"].values)), op="true", tier="T0",
                   note=f"{data['lon'].values.min():g}..{data['lon'].values.max():g}")


NINO_BOXES = {"nino1+2": (-10, 0, -90, -80), "nino3": (-5, 5, -150, -90), "nino3.4": (-5, 5, -170, -120),
              "nino4": (-5, 5, 160, -150)}


@case("compute_enso_nino_index", "a uniform signal planted in each Nino box is returned exactly")
def nino_synthetic(rec):
    template = load_regular_data("ersst", "sst", time_range=["1990-01", "2019-12"]).load()
    lat, lon = lat_lon_like(template)
    time = template["time"].to_index()
    tt = np.arange(time.size)
    seasonal = 26 + 2 * np.sin(2 * np.pi * (time.month.values - 1) / 12)[:, None, None] * np.cos(np.radians(lat))
    signal = 1.5 * np.sin(2 * np.pi * tt / 43.0) + 0.5 * np.cos(2 * np.pi * tt / 17.0)
    for region, (la0, la1, lo0, lo1) in NINO_BOXES.items():
        in_lat = (lat >= la0) & (lat <= la1)
        in_lon = ((lon >= lo0) | (lon <= lo1)) if lo1 < lo0 else ((lon >= lo0) & (lon <= lo1))
        inside = in_lat & in_lon
        values = seasonal + np.where(inside, signal[:, None, None], 50.0)   # outside: a value that must not leak
        index = compute_enso_nino_index(refill(template, values), nino_region=region, climatology="self")
        clim = pd.Series(signal).groupby(time.month.values).transform("mean").values
        rec.metric(f"{region}: index = signal - its monthly climatology", max_abs(index.values, signal - clim),
                   threshold=1e-10, tier="T0", note=f"{int(inside[0].sum())} ERSST cells in the box")


@case("compute_pdo_index", "a planted North Pacific mode with a uniform global trend")
def pdo_synthetic(rec):
    template = load_regular_data("ersst", "sst", time_range=["1980-01", "2019-12"]).load()
    lat, lon = lat_lon_like(template.isel(time=0))
    land = np.isnan(template.isel(time=0).values)
    lon360 = lon % 360
    west = np.exp(-(((lat - 40) / 10) ** 2 + ((lon360 - 180) / 25) ** 2))
    coast = np.exp(-(((lat - 45) / 8) ** 2 + ((lon360 - 225) / 8) ** 2))
    box = (lat >= 20) & (lat <= 70) & (lon360 >= 110) & (lon360 <= 260)
    pattern = np.where(box, -west + 0.8 * coast, 0.0)
    areas = load_areas(dataset="ersst").values
    ocean = ~land
    # A compensating constant outside the box gives the pattern a zero global mean, so
    # removing the global-mean anomaly cannot alter it.
    outside = ocean & ~box
    pattern = np.where(outside, -np.sum((areas * pattern)[ocean]) / np.sum(areas[outside]), pattern)
    pattern = np.where(land, np.nan, pattern)

    time = template["time"].to_index()
    n = time.size
    rng = np.random.default_rng(3)
    pc = ar1(rng, n, 0.9)
    pc = (pc - pc.mean()) / pc.std()
    seasonal = 15 + 10 * np.cos(np.radians(lat)) + 3 * np.sin(2 * np.pi * (time.month.values - 1) / 12)[:, None, None] * np.sin(np.radians(lat))
    noise = 0.05 * rng.standard_normal((n,) + lat.shape)
    base = seasonal + pattern[None] * pc[:, None, None] + noise
    trend = 0.02 * (np.arange(n) / 12.0)[:, None, None]
    index, eof, explained, regression = compute_pdo_index(refill(template, np.where(land, np.nan, base + trend)))
    index0, eof0, _, _ = compute_pdo_index(refill(template, np.where(land, np.nan, base)))

    rec.metric("index correlates with the planted PC", float(np.corrcoef(index.values, pc)[0, 1]), threshold=0.99,
               op=">=", tier="T1", note="positive: the planted pattern is cool in the west, as the convention asks")
    npac = refill(template.isel(time=0), pattern).sel(lat=eof["lat"], lon=eof["lon"]).values
    valid = np.isfinite(eof.values) & np.isfinite(npac)
    rec.metric("EOF pattern correlates with the planted pattern",
               float(np.corrcoef(eof.values[valid], npac[valid])[0, 1]), threshold=0.99, op=">=", tier="T1")
    reg_valid = np.isfinite(regression.values) & np.isfinite(npac)
    rms = np.sqrt(np.mean((regression.values[reg_valid] - npac[reg_valid]) ** 2) / np.mean(npac[reg_valid] ** 2))
    rec.metric("regression pattern recovers the planted amplitude (relative RMS)", float(rms), threshold=0.05,
               tier="T1")
    rec.metric("a uniform global trend does not change the index", max_abs(index.values, index0.values),
               threshold=1e-6, tier="T0")
    rec.metric("explained variance", explained, op="info", tier="T3")


# ---------------------------------------------------------------------------
# Horizontal regridding
# ---------------------------------------------------------------------------

def _target_cell_areas(out):
    """Unit-sphere areas of a regular target grid's cells, with great-circle edges as ESMF has them."""
    lat_e = np.clip(refs._edges(out["lat"].values), -90, 90)
    lon_e = refs._edges(out["lon"].values)
    return refs.spherical_quad_area(*refs.regular_cell_vertices(lat_e, lon_e))


@case("regrid_horizontal", "nested ERSST 2deg -> 1deg: constants, meridian steps, conservation, leakage")
def regrid_nested(rec):
    template = ersst_surface()
    lat, lon = lat_lon_like(template)
    out = regrid_horizontal(refill(template, 1.0), target="woa23")
    rec.metric("f=1 -> 1 on every valid target cell", float(np.nanmax(np.abs(out.values - 1))), threshold=1e-12,
               tier="T0")
    rec.metric("full-cover source leaves no target cell inside its latitudes empty",
               int(np.isnan(out.values[np.abs(out["lat"].values) < 89]).sum()), threshold=0, op="==", tier="T0")

    # Meridians are great circles, so a step placed on a source cell edge (odd
    # longitudes on ERSST) must come through with no leakage at all.
    step = regrid_horizontal(refill(template, (lon >= 2).astype(float)), target="woa23")
    tlon = step["lon"].values
    # Source cells centred at 2..178 cover longitudes 1..179; the one centred at
    # -180 spans the dateline (179..181) and is 0.
    expected = ((tlon > 1) & (tlon < 179)).astype(float)[None, :] * np.ones((step.sizes["lat"], 1))
    rec.metric("step on a meridian (a great circle) is reproduced exactly",
               float(np.nanmax(np.abs(step.values - expected))), threshold=1e-12, tier="T0")

    src_lat_e = refs._edges(template["lat"].values)
    src_lon_e = refs._edges(template["lon"].values)
    src_area = refs.spherical_quad_area(*refs.regular_cell_vertices(src_lat_e, src_lon_e))
    f = np.sin(np.radians(lat)) * np.cos(np.radians(3 * lon)) + np.cos(np.radians(lat)) ** 3
    out = regrid_horizontal(refill(template, f), target="woa23")
    tgt_area = _target_cell_areas(out)
    total_src = float(np.sum(src_area * f))
    total_tgt = float(np.nansum(tgt_area * out.values))
    rec.metric("area integral conserved (unit-sphere areas, great-circle edges)", abs(total_tgt / total_src - 1),
               threshold=1e-6, tier="T0", method="spherical excess of each cell (Van Oosterom & Strackee)")

    row = np.round((lat - template["lat"].values.min()) / 2).astype(int)
    stripes = regrid_horizontal(refill(template, (row % 2).astype(float)), target="woa23")
    tl = stripes["lat"].values
    src_row = np.floor((tl - (template["lat"].values.min() - 1)) / 2).astype(int)
    leak = np.abs(stripes.values - (src_row % 2).astype(float)[:, None])
    worst = np.nanmax(leak, axis=1)
    rec.metric("alternating rows: largest fraction leaked from a neighbouring source row", float(np.nanmax(worst)),
               threshold=0.01, tier="T1",
               note=f"ESMF bounds cells by great circles, not parallels; leak at 30N {worst[np.argmin(abs(tl - 30.5))]:.1e}, "
                    f"60N {worst[np.argmin(abs(tl - 60.5))]:.1e}, 80N {worst[np.argmin(abs(tl - 80.5))]:.1e}")


@case("regrid_horizontal", "curvilinear CESM2 and MPI-ESM1-2-LR -> 1deg: constants and conservation")
def regrid_curvilinear(rec):
    from src.tools.data_openers import open_cmip6_grid_file
    for model in ("CESM2", "MPI-ESM1-2-LR"):
        tos = load_regular_data("cmip6_historical", "tos", model=model, time_range=["2000-01", "2000-01"]).isel(time=0).load()
        lat, lon = lat_lon_like(tos)
        out = regrid_horizontal(refill(tos, 1.0), target="woa23")
        rec.metric(f"{model}: f=1 -> 1", float(np.nanmax(np.abs(out.values - 1))), threshold=1e-9, tier="T0",
                   note=f"{int(np.isfinite(out.values).sum())} of {out.size} target cells filled")

        grid, info = open_cmip6_grid_file("historical", model, "tos")
        grid.close()
        names = [("vertices_latitude", "vertices_longitude"), ("lat_bnds", "lon_bnds")]
        for vlat_name, vlon_name in names:
            try:
                vlat, _ = refs.raw_nc(info["source_file"], vlat_name)
                vlon, _ = refs.raw_nc(info["source_file"], vlon_name)
                if vlat.ndim == 3:
                    break
            except (KeyError, IndexError):
                continue
        # The file is read raw, so it still holds every column it was written with;
        # the loader may have dropped cyclic halo columns. The loader's grid is
        # located inside the raw one by its coordinates, the field is built on the
        # loader's grid, and the source integral is taken over the same raw columns.
        cfg = info["model_cfg"]
        raw_lat, _ = refs.raw_nc(info["source_file"], cfg["lat_coord"])
        raw_lon, _ = refs.raw_nc(info["source_file"], cfg["lon_coord"])
        ncols = tos.sizes[cfg["x_dim"]]
        offsets = [o for o in range(raw_lat.shape[1] - ncols + 1)
                   if np.allclose(raw_lat[:, o:o + ncols], lat)
                   and np.allclose(((raw_lon[:, o:o + ncols] - lon + 180) % 360) - 180, 0)]
        rec.metric(f"{model}: loader grid found inside the file's {raw_lat.shape[1]} columns",
                   bool(offsets), op="true", tier="T0",
                   note=f"loader columns {offsets[0]}..{offsets[0] + ncols - 1}" if offsets else "not found")
        off = offsets[0] if offsets else 0
        cols = slice(off, off + ncols)
        src_area = refs.spherical_quad_area(vlat[:, cols], vlon[:, cols])
        # Neither ocean grid reaches the South Pole (CESM2 stops near 79S), so the
        # field is tapered to zero poleward of 60-70 degrees, where coverage ends.
        func = lambda la, lo: (np.cos(np.radians(la)) * np.cos(np.radians(2 * lo)) + 0.3 * np.sin(np.radians(la))
                               + 1.0) * np.clip((70 - np.abs(la)) / 10, 0, 1)
        out = regrid_horizontal(refill(tos, func(lat, lon)), target="woa23")
        tgt_area = _target_cell_areas(out)
        rec.metric(f"{model}: source cells' share of the sphere (sum of areas / 4 pi - 1)",
                   src_area.sum() / (4 * np.pi) - 1, op="info", tier="T3", note="the grid stops short of the South Pole")
        total_src = np.sum(src_area * func(lat, lon))
        rec.metric(f"{model}: area integral conserved",
                   float(abs(np.nansum(tgt_area * out.values) / total_src - 1)), threshold=1e-6, tier="T1",
                   method="unit-sphere areas from the raw vertices of the loader's columns (spherical excess)")
        blob = lambda la, lo: np.exp(-(la / 5) ** 2 - ((((lo - 200) + 180) % 360 - 180) / 5) ** 2)
        out_blob = regrid_horizontal(refill(tos, blob(lat, lon)), target="woa23")
        rec.metric(f"{model}: localized blob (0N, 160W) conserved exactly",
                   float(abs(np.nansum(tgt_area * out_blob.values) / np.sum(src_area * blob(lat, lon)) - 1)),
                   threshold=1e-9, tier="T0")
        # The strip where MPI's halo column carried wrong vertices (the last kept
        # column, next to the cyclic seam): a regression check for that fix.
        seam_lat, seam_lon = float(lat[lat.shape[0] // 2, -1]), float(lon[lat.shape[0] // 2, -1])
        seam = lambda la, lo: np.exp(-((la - seam_lat) / 5) ** 2 - ((((lo - seam_lon) + 180) % 360 - 180) / 5) ** 2)
        out_seam = regrid_horizontal(refill(tos, seam(lat, lon)), target="woa23")
        rec.metric(f"{model}: blob on the cyclic seam ({seam_lat:.0f}N, {seam_lon:.0f}E) conserved exactly",
                   float(abs(np.nansum(tgt_area * out_seam.values) / np.sum(src_area * seam(lat, lon)) - 1)),
                   threshold=1e-5, tier="regression",
                   note="bounded by the file's own seam: MPI's raw vertices on either side of it disagree by "
                        "~7.5e-5 deg; before the halo fix this blob was mis-conserved by several percent")

        mask = regrid_horizontal(tos, target="woa23")
        wet = np.isfinite(mask.values)
        tlat, tlon = np.meshgrid(out["lat"].values, out["lon"].values, indexing="ij")
        truth = (np.cos(np.radians(tlat)) * np.cos(np.radians(2 * tlon)) + 0.3 * np.sin(np.radians(tlat)) + 1.0) \
            * np.clip((70 - np.abs(tlat)) / 10, 0, 1)
        err = np.abs(out.values - truth)[wet]
        bil = regrid_horizontal(refill(tos, func(lat, lon)), target="woa23", method="bilinear")
        err_b = np.abs(bil.values - truth)[wet]
        rec.metric(f"{model}: ocean cells, error vs centre values (conservative / bilinear p99)",
                   f"{np.nanpercentile(err, 99):.2e} / {np.nanpercentile(err_b, 99):.2e}", op="info", tier="T3",
                   note="first-order conservative onto a finer grid is O(grad f * dx_source), not O(dx^2); "
                        f"bilinear max over all cells {np.nanmax(np.abs(bil.values - truth)):.2e}")


# ---------------------------------------------------------------------------
# Known answers for the element-wise and event tools
# ---------------------------------------------------------------------------

@case("detect_threshold_events", "hand-built series with known runs")
def events_known(rec):
    v = np.array([1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, np.nan, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1], dtype=float)
    time = monthly_axis("2000-01", v.size)
    series = xr.DataArray(v * np.linspace(1, 2, v.size), coords={"time": time}, dims="time")
    table = detect_threshold_events(series, threshold=0.5, direction="above", min_duration=5)
    expected = [(5, 9), (12, 17), (19, 23)]
    got = [(time.get_loc(pd.Timestamp(s)), time.get_loc(pd.Timestamp(e))) for s, e in zip(table["start"], table["end"])]
    rec.metric("runs of >= 5 found, NaN breaks a run, run at the record end kept", got == expected, op="true",
               tier="T0", note=f"got {got}")
    rec.metric("durations", list(table["duration"]) == [5, 6, 5], op="true", tier="T0")
    table = detect_threshold_events(series, threshold=0.5, min_duration=6)
    rec.metric("min_duration=6 keeps only the 6-step run", len(table) == 1, op="true", tier="T0")
    table = detect_threshold_events(-series, threshold=-0.5, direction="below", min_duration=5)
    rec.metric("direction='below' mirrors 'above'", len(table) == 3 and np.allclose(table["peak_value"],
               -detect_threshold_events(series, 0.5, min_duration=5)["peak_value"]), op="true", tier="T0")
    edge = xr.DataArray([0.5] * 5 + [0.0], coords={"time": monthly_axis("2000-01", 6)}, dims="time")
    rec.metric("value equal to the threshold counts (>=)", len(detect_threshold_events(edge, 0.5, min_duration=5)) == 1,
               op="true", tier="T0")


@case("compute_difference_metrics", "offset and scaled fields give known metrics")
def difference_known(rec):
    rng = np.random.default_rng(4)
    ref = xr.DataArray(rng.standard_normal((50, 40)), dims=("y", "x"))
    m = compute_difference_metrics(ref + 0.7, ref, weighting="none")
    rec.metric("offset 0.7: bias", abs(m["bias"] - 0.7), threshold=1e-12, tier="T0")
    rec.metric("offset 0.7: centred RMSE 0", m["centered_rmse"], threshold=1e-12, tier="T0")
    rec.metric("offset 0.7: RMSE 0.7", abs(m["rmse"] - 0.7), threshold=1e-12, tier="T0")
    rec.metric("offset 0.7: std ratio 1", abs(m["std_ratio"] - 1), threshold=1e-12, tier="T0")
    mean = float(ref.mean())
    m = compute_difference_metrics(mean + 2 * (ref - mean), ref, weighting="none")
    rec.metric("scaled by 2 about the mean: std ratio 2", abs(m["std_ratio"] - 2), threshold=1e-12, tier="T0")
    rec.metric("scaled by 2 about the mean: bias 0", abs(m["bias"]), threshold=1e-12, tier="T0")
    template = woa_surface()
    base = refill(template, 1.0)
    spike = base.copy()
    spike.loc[{"lat": 30.5, "lon": -40.5}] = 5.0
    m = compute_difference_metrics(spike, base)
    rec.metric("planted spike located", m["max_abs_error_at"] == {"lat": 30.5, "lon": -40.5}, op="true", tier="T0",
               note=str(m["max_abs_error_at"]))
    area = load_areas(dataset="woa23").sel(lat=30.5, lon=-40.5)
    expected_bias = 4.0 * float(area) / float(load_areas(dataset="woa23").sum())
    rec.metric("area-weighted bias of the spike = 4 * its area share", rel_err(m["bias"], expected_bias),
               threshold=1e-10, tier="T0")


@case("compute_ratio", "known ratios, prefixes and zero denominators")
def ratio_known(rec):
    x = xr.DataArray(np.array([1.0, 2.0, 0.5, 3.0]), dims="n", attrs={"units": "umol kg-1"})
    r = compute_ratio(16 * x, x)
    rec.metric("16x / x = 16", max_abs(r.values, 16.0), threshold=1e-12, tier="T0")
    n = xr.DataArray(np.array([1600.0]), dims="n", attrs={"units": "umol kg-1"})
    p = xr.DataArray(np.array([0.1]), dims="n", attrs={"units": "mmol kg-1"})
    rec.metric("umol/kg over mmol/kg rescaled: 1600 / 0.1 mmol = 16", abs(float(compute_ratio(n, p)) - 16),
               threshold=1e-12, tier="T0")
    zero = xr.DataArray(np.array([1.0, 0.0]), dims="n", attrs={"units": "umol kg-1"})
    rec.metric("zero denominator -> NaN", bool(np.isnan(compute_ratio(zero, zero).values[1])), op="true", tier="T0")


@case("transform_unit", "conversion factors, round trips and density fields")
def transform_known(rec):
    x = xr.DataArray(np.array([1.0, 250.0]), dims="n")
    rec.metric("1 umol/kg -> 0.001025 mol m-3 at 1.025 kg/L",
               rel_err(transform_unit(x, "umol/kg", "mol m-3").values, x.values * 1.025e-3), threshold=1e-14, tier="T0")
    rec.metric("mmol/m3 equals umol/L", rel_err(transform_unit(x, "mmol/m3", "umol/L").values, x.values),
               threshold=1e-14, tier="T0")
    rec.metric("K -> degC", max_abs(transform_unit(x + 273.15, "K", "degC").values, x.values), threshold=1e-12,
               tier="T0")
    rec.metric("cm/s -> m/s", rel_err(transform_unit(x, "cm/s", "m/s").values, x.values * 0.01), threshold=1e-14,
               tier="T0")
    there = transform_unit(transform_unit(x, "umol/kg", "mol m-3"), "mol m-3", "umol/kg")
    rec.metric("round trip umol/kg -> mol m-3 -> umol/kg", rel_err(there.values, x.values), threshold=1e-14, tier="T0")
    rho_abs = xr.DataArray(np.array([1020.0, 1030.0]), dims="n", attrs={"long_name": "density"})
    rho_sig = xr.DataArray(np.array([20.0, 30.0]), dims="n",
                           attrs={"long_name": "Potential density anomaly (sigma_0)"})
    a = transform_unit(x, "umol/kg", "mmol m-3", density=rho_abs).values
    b = transform_unit(x, "umol/kg", "mmol m-3", density=rho_sig).values
    rec.metric("density field: each cell at its own density (c * rho / 1000)", rel_err(a, x.values * rho_abs.values / 1000),
               threshold=1e-14, tier="T0")
    rec.metric("sigma field and absolute density give the same result", rel_err(b, a), threshold=1e-14, tier="T0")


@case("locate_extreme", "planted extremes, signed abs_max and ties")
def extreme_known(rec):
    rng = np.random.default_rng(5)
    time = monthly_axis("2000-01", 6)
    field = xr.DataArray(rng.uniform(0, 1, (6, 10, 12)), dims=("time", "lat", "lon"),
                         coords={"time": time, "lat": np.arange(10.0), "lon": np.arange(12.0) * 10})
    field[3, 4, 7] = 5.0
    field[1, 2, 2] = -7.0
    result = locate_extreme(field)
    rec.metric("max located", result["coordinates"] == {"time": "2000-04", "lat": 4.0, "lon": 70.0}, op="true",
               tier="T0", note=str(result["coordinates"]))
    result = locate_extreme(field, op="abs_max")
    rec.metric("abs_max keeps the sign", result["value"] == -7.0, op="true", tier="T0")
    field[5, 9, 11] = 5.0
    rec.metric("two cells at the maximum -> ties = 2", locate_extreme(field)["ties"], threshold=2, op="==", tier="T0")


@case("detect_outliers", "planted outliers are flagged by every method")
def outliers_planted(rec):
    rng = np.random.default_rng(6)
    values = rng.standard_normal(5000)
    planted = rng.choice(5000, 10, replace=False)
    values[planted] = np.where(np.arange(10) % 2 == 0, 8.0, -8.0)
    series = xr.DataArray(values, dims="time")
    for method in ("zscore", "iqr", "mad", "percentile"):
        mask = detect_outliers(series, method=method, dim="time").values
        rec.metric(f"{method}: all 10 planted +-8 sigma values flagged", int(mask[planted].sum()), threshold=10,
                   op="==", tier="T0", note=f"{int(mask.sum())} flagged in total")
