"""
Category 4: physical magnitudes against published benchmarks.
"""
from __future__ import annotations

import numpy as np

from src.eval.tool_verification import refs
from src.eval.tool_verification.common import case, zone_area
from src.tools.basic_calculator import compute_exceedance_fraction
from src.tools.data_loader import (cell_areas_from_vertices, load_areas, load_climatology_data, load_depth_levels,
                                   load_landsea_mask, load_regular_data, load_volumes_for_cmip6)
from src.tools.data_openers import open_cmip6_grid_file
from src.tools.scientific_calculator import identify_oxygen_minimum_zone

CAT = "cat4_magnitude"
MODELS = ("CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR")
OCEAN_AREA, OCEAN_VOLUME, MEAN_DEPTH = 3.61e14, 1.332e18, 3682.0


@case("load_areas", "global ocean area: WOA23 1-degree grid and the three models' areacello")
def ocean_area(rec):
    areas = load_areas()
    ocean = (load_landsea_mask() > 1).values
    total = float(areas.values[ocean].sum())
    rec.metric("WOA23: area of cells the land-sea mask calls ocean / 3.61e14", total / OCEAN_AREA, op="info",
               tier="T3", note=f"{total:.4e} m2; a binary 1-degree mask has no partial coastal cells, so it is "
                               f"not held to the +-1% meant for a model's fractional areacello")
    rec.metric("1-degree grid total / WGS84 surface area", float(areas.sum()) / float(zone_area(-90, 90)),
               op="info", tier="T3", note="short only by the two 0.01-degree polar caps the area solver clips")
    for model in MODELS:
        a = load_areas(dataset="cmip6_historical", model=model)
        total = float(a.sum())
        rec.metric(f"{model}: sum of areacello / 3.61e14", total / OCEAN_AREA, threshold=(0.99, 1.01), op="in",
                   tier="T3", note=f"{total:.4e} m2; {a.attrs.get('area_source', '')[:80]}")
        if model == "MPI-ESM1-2-LR":
            lat, lon = a[a.attrs["lat_coord"]].values, a[a.attrs["lon_coord"]].values
            dup = refs.duplicate_columns(lat, lon)
            keep = np.ones(a.shape[1], dtype=bool)
            keep[dup] = False
            rec.metric(f"{model}: the same with its repeated columns {dup} counted once",
                       float(np.nansum(a.values[:, keep])) / OCEAN_AREA, op="info", tier="T3")


@case("load_volumes_for_cmip6", "global ocean volume and mean depth of each model")
def ocean_volume(rec):
    for model in MODELS:
        kw = {"time_range": ["2000-01", "2000-01"]} if model == "GFDL-ESM4" else {}
        v = load_volumes_for_cmip6("cmip6_historical", model, **kw)
        if "time" in v.dims:
            v = v.isel(time=0)
        volume = float(v.sum())
        area = float(load_areas(dataset="cmip6_historical", model=model).sum())
        rec.metric(f"{model}: global volume / 1.332e18", volume / OCEAN_VOLUME, threshold=(0.97, 1.03), op="in",
                   tier="T3", note=f"{volume:.4e} m3; {v.attrs.get('volume_source', '')[:80]}")
        rec.metric(f"{model}: mean depth (volume / area) / 3682 m", volume / area / MEAN_DEPTH,
                   threshold=(0.95, 1.05), op="in", tier="T3", note=f"{volume / area:.0f} m")


@case("load_volumes_for_cmip6", "published volcello vs the rebuild from areacello, layer bounds and deptho")
def volume_rebuild(rec):
    for model in MODELS:
        kw = {"time_range": ["2000-01", "2000-01"]} if model == "GFDL-ESM4" else {}
        published = load_volumes_for_cmip6("cmip6_historical", model, **kw)
        if "time" in published.dims:
            published = published.isel(time=0)
        areas = load_areas(dataset="cmip6_historical", model=model)
        _, bounds = load_depth_levels(dataset="cmip6_historical", model=model)
        deptho = load_regular_data("cmip6_historical", "deptho", model=model)
        if "time" in deptho.dims:
            deptho = deptho.isel(time=0)
        top = bounds[:, 0].values[:, None, None]
        bottom = bounds[:, 1].values[:, None, None]
        floor = deptho.values[None]
        thickness = np.clip(np.minimum(bottom, floor) - top, 0, None)
        rebuilt = float(np.nansum(areas.values[None] * thickness))
        rec.metric(f"{model}: rebuilt / published global volume", rebuilt / float(published.sum()), op="info",
                   tier="T3", note="docstring: GFDL-ESM4 1.53e18 vs 1.33e18; MPI partial bottom cells; CESM2 agree")


@case("load_areas", "the fallback that computes areas from cell vertices vs each model's published areacello")
def area_fallback(rec):
    for model in ("CESM2", "MPI-ESM1-2-LR"):
        published = load_areas(dataset="cmip6_historical", model=model)
        grid, info = open_cmip6_grid_file("historical", model, "tos")
        names = ("vertices_latitude", "vertices_longitude") if "vertices_latitude" in grid else ("lat_bnds", "lon_bnds")
        computed = cell_areas_from_vertices(grid[names[0]].values, grid[names[1]].values)
        grid.close()
        ratio = computed / published.values
        ok = np.isfinite(ratio)
        rec.metric(f"{model}: computed / published cell area, median", float(np.median(ratio[ok])),
                   threshold=(0.99, 1.01), op="in", tier="T3",
                   note=f"p1 {np.percentile(ratio[ok], 1):.4f}, p99 {np.percentile(ratio[ok], 99):.4f}; "
                        f"docstring: about 0.5% (WGS84 vs the model's sphere)")


@case("identify_oxygen_minimum_zone", "volume of low-oxygen water in WOA23")
def omz_volume(rec):
    o2 = load_climatology_data("woa23", "oxygen")
    for threshold in (20.0, 60.0):
        result = compute_exceedance_fraction(o2, threshold=threshold, comparison="less")
        rec.metric(f"share of the WOA23 ocean volume with O2 < {threshold:g} umol/kg", result["fraction_percent"],
                   op="info", tier="T3", note=f"{result['exceeding_amount']:.3e} m3; literature puts the < 20 "
                                              f"umol/kg core at the order of 1% of the ocean")
    thickness, _ = identify_oxygen_minimum_zone(o2, threshold=20.0)
    for label, box in (("eastern tropical North Pacific", dict(lat=slice(5, 25), lon=slice(-120, -90))),
                       ("Arabian Sea", dict(lat=slice(10, 25), lon=slice(55, 72)))):
        rec.metric(f"{label}: maximum OMZ (<20) thickness (m)", float(thickness.sel(**box).max()), op="info",
                   tier="T3", note="these are the two thickest OMZs of the world ocean")
