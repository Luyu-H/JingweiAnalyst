"""
Category 5: identities and internal consistency.
"""
from __future__ import annotations

import numpy as np

from src.eval.tool_verification import refs
from src.eval.tool_verification.common import case, lat_lon_like, max_abs, refill
from src.tools.basic_calculator import compute_basic_stat, compute_difference_metrics, compute_exceedance_fraction
from src.tools.data_loader import (areas_like, load_areas, load_basin_mask, load_climatology_data,
                                   load_data_with_geolocator, load_depth_levels, load_landsea_mask, load_regular_data)
from src.tools.preprocessing import BASIN_GROUPS, regrid_horizontal, select_basin
from src.tools.scientific_calculator import compute_enso_nino_index, identify_oxygen_minimum_zone

CAT = "cat5_consistency"
MODELS = ("CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR")


@case("compute_difference_metrics", "rmse^2 = bias^2 + centred_rmse^2 for full and partial reductions")
def difference_identity(rec):
    a = load_regular_data("cmip6_historical", "tos", model="GFDL-ESM4", time_range=["2000-01", "2000-12"]).load()
    b = load_regular_data("cmip6_historical", "tos", model="GFDL-ESM4", time_range=["2001-01", "2001-12"]).load()
    b = b.assign_coords(time=a["time"])
    for label, dim in (("whole field, area-weighted", None), ("per time step over lat/lon", ["lat", "lon"]),
                       ("per cell over time", "time")):
        m = compute_difference_metrics(a, b, dim=dim)
        rmse2 = np.asarray(m["rmse"], dtype="float64") ** 2
        residual = np.abs(rmse2 - np.asarray(m["bias"]) ** 2 - np.asarray(m["centered_rmse"]) ** 2)
        rec.metric(f"{label}: max |rmse^2 - bias^2 - crmse^2| / rmse^2", float(np.nanmax(residual / rmse2)),
                   threshold=1e-5, tier="T0", note=f"float32 input; {m['weighting'][:80]}")


@case("select_basin", "the ocean groupings partition WOA23's 60 codes and the ocean area")
def basin_partition(rec):
    codes = [c for g in BASIN_GROUPS.values() for c in g["main"] + g["marginal"]] + [53]
    rec.metric("every code 1-60 assigned exactly once (Caspian on its own)", sorted(codes) == list(range(1, 61)),
               op="true", tier="T0", note=f"{len(codes)} assignments, {len(set(codes))} distinct")
    template = load_climatology_data("woa23", "temperature", depth_range=[0, 0]).isel(depth=0).load()
    surface_codes = load_basin_mask().isel(depth=0).values
    ocean = np.isfinite(surface_codes)
    full = refill(template, 1.0)
    areas = load_areas(dataset="woa23").values
    covered = np.zeros(full.shape, dtype=int)
    for basin in BASIN_GROUPS:
        covered += np.isfinite(select_basin(full, basin).values)
    covered += np.isfinite(select_basin(full, [53]).values)
    rec.metric("WOA23 grid: no cell in two basins", int(covered.max()), threshold=1, op="<=", tier="T0")
    rec.metric("WOA23 grid: area in some basin / area with a basin code", float(areas[covered > 0].sum() / areas[ocean].sum()),
               threshold=1 - 1e-12, op=">=", tier="T0")
    atlantic = select_basin(full, "Atlantic")
    direct = np.isin(surface_codes, BASIN_GROUPS["Atlantic"]["main"] + BASIN_GROUPS["Atlantic"]["marginal"])
    rec.metric("WOA23 grid: Atlantic selection equals the code mask", bool(np.array_equal(np.isfinite(atlantic.values), direct)),
               op="true", tier="T0")


@case("select_basin", "basin areas on the models' own grids vs the WOA23 grid")
def basin_on_models(rec):
    woa = refill(load_climatology_data("woa23", "temperature", depth_range=[0, 0]).isel(depth=0).load(), 1.0)
    woa_area = load_areas(dataset="woa23").values
    reference = {b: float(woa_area[np.isfinite(select_basin(woa, b).values)].sum()) for b in BASIN_GROUPS}
    for model in MODELS:
        tos = load_regular_data("cmip6_historical", "tos", model=model, time_range=["2000-01", "2000-01"]).isel(time=0).load()
        areas = load_areas(dataset="cmip6_historical", model=model).values
        covered = np.zeros(tos.shape, dtype=int)
        shares = {}
        for basin in BASIN_GROUPS:
            sel = np.isfinite(select_basin(tos, basin).values)
            covered += sel
            shares[basin] = float(np.nansum(areas[sel])) / reference[basin]
        rec.metric(f"{model}: no model cell in two basins", int(covered.max()), threshold=1, op="<=", tier="T0")
        ocean = np.isfinite(tos.values)
        rec.metric(f"{model}: ocean area left out of every basin (%)",
                   100 * float(np.nansum(areas[ocean & (covered == 0)]) / np.nansum(areas[ocean])), op="info",
                   tier="T3", note="cells whose nearest 1-degree cell is land in the WOA23 mask")
        rec.metric(f"{model}: basin area / WOA23 basin area",
                   {k: round(v, 4) for k, v in shares.items()}, op="info", tier="T3")


@case("load_basin_mask", "surface basin codes cover exactly the land-sea mask's ocean; ocean shrinks with depth")
def basin_vs_landsea(rec):
    basin = load_basin_mask()
    landsea = load_landsea_mask()
    surface = np.isfinite(basin.isel(depth=0).values)
    rec.metric("surface cells with a code == cells with landsea > 1", bool(np.array_equal(surface, (landsea > 1).values)),
               op="true", tier="T0", note=f"{int(surface.sum())} vs {int((landsea > 1).sum())}")
    counts = np.isfinite(basin.values).sum(axis=(1, 2))
    rec.metric("ocean cell count never increases with depth", bool(np.all(np.diff(counts) <= 0)), op="true", tier="T0")
    values = basin.values[np.isfinite(basin.values)]
    rec.metric("codes lie in 1..60", bool(values.min() >= 1 and values.max() <= 60), op="true", tier="T0")


@case("load_landsea_mask", "(mask - 1) equals the number of WOA23 annual levels holding data")
def landsea_vs_woa(rec):
    landsea = load_landsea_mask().values.astype(int)
    temp = load_climatology_data("woa23", "temperature").values
    levels = np.isfinite(temp).sum(axis=0)
    expected = np.clip(landsea - 1, 0, temp.shape[0])
    agree = float(np.mean(levels == expected))
    rec.metric("share of cells where the two agree", agree, threshold=0.999, op=">=", tier="T0",
               note=f"{int((levels != expected).sum())} cells differ")


@case("load_regular_data", "model grids hold each cell once (no repeated cyclic columns)")
def unique_cells(rec):
    for model in ("CESM2", "MPI-ESM1-2-LR"):
        tos = load_regular_data("cmip6_historical", "tos", model=model, time_range=["2000-01", "2000-01"]).isel(time=0)
        lat, lon = lat_lon_like(tos)
        dup = refs.duplicate_columns(lat, lon)
        areas = load_areas(dataset="cmip6_historical", model=model).values
        share = float(np.nansum(areas[:, dup]) / np.nansum(areas)) if dup else 0.0
        rec.metric(f"{model}: columns whose every cell repeats another column", len(dup), threshold=0, op="==",
                   tier="T0", note=f"columns {dup}; their share of the ocean area {100 * share:.2f}%" if dup else "")


@case("load_depth_levels", "layer bounds are contiguous, start at the surface and contain their centres")
def depth_invariants(rec):
    grids = [("woa23", dict(dataset="woa23", variable=v, time_slice=t))
             for v in ("temperature", "salinity", "oxygen", "nitrate") for t in ("annual", "monthly")]
    grids += [("glodap", dict(dataset="glodap"))]
    grids += [(m, dict(dataset="cmip6_historical", model=m)) for m in MODELS]
    for label, kw in grids:
        levels, bounds = load_depth_levels(**kw)
        z, b = levels.values, bounds.values
        name = f"{label} {kw.get('variable', '')} {kw.get('time_slice', '')}".strip()
        ok = (np.allclose(b[1:, 0], b[:-1, 1]) and abs(b[0, 0]) < 1e-9 and np.all(b[:, 0] <= z) and np.all(z <= b[:, 1])
              and np.all(np.diff(b, axis=1) > 0))
        rec.metric(f"{name}: contiguous, from 0 m, centres inside, thickness > 0", bool(ok), op="true", tier="T0",
                   note=f"{z.size} levels to {b[-1, 1]:g} m")


@case("compute_enso_nino_index", "a model's Nino3.4 on its native grid vs after conservative regridding to 1 degree")
def nino_native_vs_regridded(rec):
    for model in ("CESM2", "MPI-ESM1-2-LR"):
        tos = load_regular_data("cmip6_historical", "tos", model=model, time_range=["1980-01", "2009-12"],
                                lat_range=[-15, 15], lon_range=[180, 260]).load()
        native = compute_enso_nino_index(tos, climatology="self")
        full = load_regular_data("cmip6_historical", "tos", model=model, time_range=["1980-01", "2009-12"]).load()
        regridded = compute_enso_nino_index(regrid_horizontal(full, target="woa23"), climatology="self")
        r = float(np.corrcoef(native.values, regridded.values)[0, 1])
        rec.metric(f"{model}: correlation native / regridded", r, threshold=0.99, op=">=", tier="T3",
                   note=f"RMS difference {np.sqrt(np.mean((native.values - regridded.values) ** 2)):.4f} degC")


@case("compute_basic_stat", "a model's global-mean SST on its native grid vs after conservative regridding")
def global_mean_native_vs_regridded(rec):
    for model in MODELS:
        tos = load_regular_data("cmip6_historical", "tos", model=model, time_range=["2000-01", "2000-12"]).load()
        native = compute_basic_stat(tos, "area_weighted_mean", dim=["lat", "lon"])
        regridded = compute_basic_stat(regrid_horizontal(tos, target="woa23"), "area_weighted_mean", dim=["lat", "lon"])
        if model != "GFDL-ESM4":
            rec.metric(f"{model}: max |native - regridded| global-mean SST (degC)",
                       max_abs(native.values, regridded.values), threshold=0.05, tier="T3",
                       note="coastal cells and the na_thres=0.5 rule move the coastline slightly")
        else:
            # GFDL-ESM4's areacello holds only the ocean part of its coastal cells, and
            # a regridded field is weighted with the target's whole-cell areas, so its
            # global mean moves by ~0.2 degC (an accepted limitation, not tested).
            # GFDL-ESM4 is archived on the WOA23 1-degree grid itself, so the regrid is
            # an identity on the values; what changes is the weight of coastal cells.
            first = tos.isel(time=0)
            same = regrid_horizontal(first, target="woa23")
            both = np.isfinite(first.values) & np.isfinite(same.values)
            rec.metric("GFDL-ESM4: regridding onto its own grid leaves every value unchanged",
                       max_abs(first.values[both], same.values[both]), threshold=0.0, tier="T0")
            full = load_areas(dataset="woa23").values
            own = load_areas(dataset="cmip6_historical", model=model).values.astype("float64")
            ok = np.isfinite(first.values)
            v = first.values.astype("float64")[ok]
            partial = int((own[ok] / full[ok] < 0.99).sum())
            rec.metric("GFDL-ESM4: the regridded mean equals the native field weighted with whole-cell areas",
                       abs(float(np.sum(v * full[ok]) / np.sum(full[ok])) - float(regridded.isel(time=0))),
                       threshold=1e-5, tier="T0",
                       note=f"{partial} coastal cells hold a partial ocean area in areacello; after regridding they "
                            f"carry the target's whole-cell area")


@case("identify_oxygen_minimum_zone", "sum of OMZ thickness x area equals the volume compute_exceedance_fraction gives")
def omz_thickness_vs_volume(rec):
    o2 = load_climatology_data("woa23", "oxygen", lat_range=[-30, 30], lon_range=[40, 120]).load()
    thickness, _ = identify_oxygen_minimum_zone(o2, threshold=20.0)
    areas = load_areas(dataset="woa23", lat_range=[-30, 30], lon_range=[40, 120]).values
    volume_from_thickness = float(np.nansum(thickness.values * areas))
    volume = compute_exceedance_fraction(o2, threshold=20.0, comparison="less")["exceeding_amount"]
    rec.metric("relative difference", abs(volume_from_thickness / volume - 1), threshold=1e-10, tier="T0",
               note=f"{volume:.4e} m3")


@case("load_data_with_geolocator", "one named sea on four grids: area of the masked ocean cells")
def region_across_grids(rec):
    name = "Mediterranean Sea"
    rows = {}
    rows["woa23"] = load_data_with_geolocator("woa23", "temperature", loc_name=name, depth_range=[0, 0]).isel(depth=0)
    rows["ersst"] = load_data_with_geolocator("ersst", "sst", loc_name=name, time_range=["2000-01", "2000-01"]).isel(time=0)
    for model in MODELS:
        rows[model] = load_data_with_geolocator("cmip6_historical", "tos", model=model, loc_name=name,
                                                time_range=["2000-01", "2000-01"]).isel(time=0)
    result = {label: float(areas_like(da).where(da.notnull()).sum()) / 1e12 for label, da in rows.items()}
    rec.metric(f"{name}: wet area on each grid (10^6 km2)", {k: round(v, 3) for k, v in result.items()}, op="info",
               tier="T3", note="Natural Earth polygon 2.5e6 km2 approx.; coarse grids gain coastal cells")
