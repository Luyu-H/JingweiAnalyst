# tools: load_regular_data, load_climatology_data, load_data_with_geolocator, load_areas, load_volumes_for_cmip6, load_depth_levels, load_basin_mask, load_landsea_mask

# Two entry points, by the kind of data:
#   load_regular_data         anything with a real time axis: the gridded records
#                             (CMIP6 historical and ssp245, ERSSTv5) and the
#                             ungridded published index series (CPC/NCEI)
#   load_climatology_data     climatologies (WOA23, GLODAP, SEANOE MLD)
#   load_data_with_geolocator either gridded kind, selected by named region

import os

import numpy as np
import pandas as pd
import xarray as xr
import yaml
from pyproj import Geod

from src.tools.registry import register_tools
from src.tools.data_openers import (
    cmip6_experiments,
    depth_scale_from_units,
    get_dataset_cfg,
    load_dataset_config,
    open_cmip6,
    open_cmip6_grid_file,
    open_cpc_index,
    open_ersst,
    open_glodap,
    open_seanoe_mld,
    open_woa23,
    snap_coords,
    sort_ascending,
    to_lon180,
    woa23_files,
)
from src.utils.geolocator import Geolocator


PATHS_CONFIG = os.path.normpath(
    os.path.join(os.path.dirname(__file__), "..", "..", "configs", "default_paths.yaml")
)

# Geolocator grid.
GEOLOCATOR_LATS = np.arange(-89.5, 90.5, 1.0)
GEOLOCATOR_LONS = np.arange(0.5, 360.5, 1.0)


# ---------------------------------------------------------------------------
# Subsetting primitives, called only by the loaders below
# ---------------------------------------------------------------------------

# How closely a stored depth coordinate must match the dataset's own level definitions to count as the same level. 
DEPTH_MATCH_TOL = 1e-4


def normalize_lon_range(lon_range):
    """
    Accept a longitude window in any convention and return it as
    (low, high, wraps) in the -180~180 convention, or None for a full circle.
    """
    low, high = float(lon_range[0]), float(lon_range[1])
    if abs(high - low) >= 360:
        return None
    low = ((low + 180) % 360) - 180
    high = ((high + 180) % 360) - 180
    return low, high, low > high


def box_mask_2d(lat2d, lon2d, lat_range=None, lon_range=None):
    """Boolean mask over a curvilinear grid's 2-D lat/lon coordinates."""
    mask = xr.ones_like(lat2d, dtype=bool)
    if lat_range is not None:
        mask = mask & (lat2d >= lat_range[0]) & (lat2d <= lat_range[1])
    if lon_range is not None:
        window = normalize_lon_range(lon_range)
        if window is not None:
            low, high, wraps = window
            if wraps:
                mask = mask & ((lon2d >= low) | (lon2d <= high))
            else:
                mask = mask & (lon2d >= low) & (lon2d <= high)
    return mask


def trim_to_mask(da, mask, y_dim, x_dim):
    """Shrink the index dimensions of a curvilinear array to the mask's
    bounding box, so the caller does not carry a global array of mostly NaN."""
    y_hits = np.where(mask.any(dim=x_dim).values)[0]
    x_hits = np.where(mask.any(dim=y_dim).values)[0]
    if y_hits.size == 0 or x_hits.size == 0:
        raise ValueError("The requested region selects no grid cells on this model's native grid.")
    return da.isel({
        y_dim: slice(int(y_hits.min()), int(y_hits.max()) + 1),
        x_dim: slice(int(x_hits.min()), int(x_hits.max()) + 1),
    })


def mask_1deg_to_grid(mask_1deg, da, grid_type, lat_coord, lon_coord):
    """
    Put a 1-degree (lat, lon) boolean mask onto whatever grid `da` lives on.
    """
    if grid_type == "regular":
        return mask_1deg.reindex(
            {"lat": da[lat_coord].values, "lon": da[lon_coord].values}, method="nearest",
        )

    lat2d, lon2d = da[lat_coord], da[lon_coord]
    row = np.clip(np.floor(lat2d.values + 90).astype(int), 0, mask_1deg.sizes["lat"] - 1)
    col = np.clip(np.floor(lon2d.values + 180).astype(int), 0, mask_1deg.sizes["lon"] - 1)
    return xr.DataArray(
        mask_1deg.values[row, col], dims=lat2d.dims, coords=lat2d.coords,
    )


def subset_regular(da, lat_range=None, lon_range=None):
    """
    Slice a regular grid by latitude and longitude.
    """
    check_ascending("lat", lat_range)
    if lat_range is not None:
        da = da.sel(lat=slice(lat_range[0], lat_range[1]))
    if lon_range is not None:
        window = normalize_lon_range(lon_range)
        if window is not None:
            low, high, wraps = window
            if wraps:
                east = da.sel(lon=slice(low, float(da.lon.max())))
                west = da.sel(lon=slice(float(da.lon.min()), high))
                west = west.assign_coords(lon=west["lon"] + 360)
                da = xr.concat([east, west], dim="lon")
                da.attrs["lon_convention"] = (
                    "continuous across the dateline: longitudes run past 180 (e.g. 160 to 210) so that the "
                    "window stays contiguous and the axis stays monotonic"
                )
            else:
                da = da.sel(lon=slice(low, high))
    return da


def subset_curvilinear(da, lat_range=None, lon_range=None):
    """Mask and trim a curvilinear grid by a lat/lon box."""
    lat_coord, lon_coord = da.attrs["lat_coord"], da.attrs["lon_coord"]
    y_dim, x_dim = da.attrs["y_dim"], da.attrs["x_dim"]
    mask = box_mask_2d(da[lat_coord], da[lon_coord], lat_range, lon_range)
    return trim_to_mask(da.where(mask), mask, y_dim, x_dim)


def check_ascending(name, value_range):
    """
    Reject a reversed latitude or depth window.
    """
    if value_range is None:
        return
    low, high = value_range[0], value_range[1]
    if low > high:
        raise ValueError(
            f"{name}_range must run from low to high, got [{low}, {high}]. "
            f"(Only lon_range may read high-to-low, to cross the dateline.)"
        )


def subset_spatial(da, lat_range=None, lon_range=None):
    """Dispatch spatial subsetting on the grid type recorded by the opener."""
    check_ascending("lat", lat_range)
    if lat_range is None and lon_range is None:
        return da
    if da.attrs.get("grid_type") == "curvilinear":
        return subset_curvilinear(da, lat_range, lon_range)
    return subset_regular(da, lat_range, lon_range)


def subset_depth(da, depth_range=None):
    check_ascending("depth", depth_range)
    if depth_range is None or "depth" not in da.coords:
        return da
    subset = da.sel(depth=slice(depth_range[0], depth_range[1]))
    if subset.sizes.get("depth", 0) == 0:
        available = da["depth"].values
        raise ValueError(
            f"depth_range {depth_range} selects no levels; this field has {available.size} levels from "
            f"{available.min():g} to {available.max():g} m."
        )
    return subset


# ---------------------------------------------------------------------------
# Resolving a dataset and routing it to its opener
# ---------------------------------------------------------------------------

def resolve_dataset(dataset: str) -> str:
    """Accept a registry id or a filesystem path under a registered root."""
    registry = load_dataset_config()["datasets"]
    if dataset in registry:
        return dataset
    candidate = os.path.normpath(dataset)
    matches = [
        (name, cfg["root"]) for name, cfg in registry.items()
        if candidate == os.path.normpath(cfg["root"])
        or candidate.startswith(os.path.normpath(cfg["root"]) + os.sep)
    ]
    if matches:
        # Prefer the deepest matching root (historical vs ssp245 share a parent).
        return max(matches, key=lambda item: len(item[1]))[0]
    raise KeyError(f"Unknown dataset '{dataset}'. Available: {sorted(registry)}")


def open_dataset_array(dataset: str, variable: str, model: str = None, time_slice: str = "annual",
                       period: str = None, field: str = None, version: str = None,
                       product: str = None, chunk_time: int = 12) -> xr.DataArray:
    """Route to the opener that owns this dataset."""
    dataset = resolve_dataset(dataset)

    if dataset.startswith("cmip6_"):
        if model is None:
            raise ValueError(
                f"'{dataset}' needs a model: one of {sorted(get_dataset_cfg(dataset)['models'])}"
            )
        return open_cmip6(dataset.replace("cmip6_", ""), model, variable, chunk_time=chunk_time)
    if dataset == "ersst":
        return open_ersst(variable or "sst")
    if dataset == "woa23":
        return open_woa23(variable, time_slice=time_slice, period=period, field=field or "an")
    if dataset == "glodap":
        return open_glodap(variable, field=field or "main")
    if dataset == "seanoe_mld":
        return open_seanoe_mld(version=version, variable=variable or "mld_dr003")
    if dataset == "cpc_indices":
        # `variable` names the region the series is published for. It stays None when
        # unstated: the products defined on a single region or on a whole basin reject
        # a region outright, and cannot tell a default apart from a real request.
        return open_cpc_index(product=product, region=variable, field=field or "anomaly")
    raise KeyError(f"No opener registered for dataset '{dataset}'")


def require_kind(dataset: str, expected, tool: str, alternative: str):
    dataset = resolve_dataset(dataset)
    kind = get_dataset_cfg(dataset)["kind"]
    accepted = (expected,) if isinstance(expected, str) else tuple(expected)
    if kind not in accepted:
        raise ValueError(
            f"'{dataset}' is a '{kind}' dataset; {tool} handles "
            f"{' or '.join(repr(k) for k in accepted)} datasets. Use {alternative}."
        )
    return dataset


# ---------------------------------------------------------------------------
# 1. Time series
# ---------------------------------------------------------------------------

@register_tools(
    name="load_regular_data",
    description="Load a rectangular subset of a time series dataset.\n"
        "This is the entry point for any question about change over time — trends, variability, scenario "
        "projections — and for the published index series (cpc_indices), which have a time axis but no "
        "grid. Give the window you need: reading a global field and slicing it afterwards is slower "
        "and, on a curvilinear model grid, wrong, because index dimensions cannot be sliced by latitude "
        "and longitude. A window crossing the dateline is handled, and comes back on a continuous "
        "longitude axis.\n"
        "The result stays lazy, so a large request is cheap until it is reduced.",
    args={
        "dataset": "str, required. The dataset id, from the registry of available datasets.",
        "variable": "str, required. The variable's short name in that dataset. For an index dataset it "
            "names the region the series is published for, e.g. 'NINO3.4'; leave it out for an index "
            "that is defined on its own domain rather than on a region box.",
        "model": "str, required for multi-model datasets. Which model to read.",
        "product": "str, required for a dataset that publishes several series of the same quantity. "
            "Which one to read; the registry lists each product with its baseline period. Default is "
            "the dataset's operational series.",
        "field": "str, optional, default 'anomaly' where it applies. Which published field to read: "
            "'anomaly' for the departure from the baseline climatology, or 'total' for the absolute "
            "value.",
        "time_range": "list[str, str], optional. Continuous window as ['start', 'end'] in ISO 8601, e.g. "
            "['1990-01', '2014-12']. Default is the whole record.",
        "depth_range": "list[float, float], optional. [min, max] in metres, positive downward. Levels are "
            "selected by their centre. Default is every level.",
        "lat_range": "list[float, float], optional. [min, max] in degrees North.",
        "lon_range": "list[float, float], optional. [min, max] in degrees East, in either the 0-360 or the "
            "-180-180 convention; a window crossing the dateline is handled.",
        "chunk_time": "int, optional, default 12. Chunk size along time for lazy reading."
    },
    returns="xr.DataArray on the dataset's own grid, with units and the provenance attributes that let the "
        "weighting and derived-quantity tools find the right cell areas, volumes and layer thicknesses. On "
        "a regular grid the dimensions are (time, depth, lat, lon); on a curvilinear model grid they are "
        "the model's index dimensions with 2-D latitude and longitude coordinates, trimmed to the "
        "requested box. An index dataset comes back with 'time' alone and no spatial dimensions."
)
def load_regular_data(dataset: str, variable: str = None, model: str = None, time_range: list = None,
                      depth_range: list = None, lat_range: list = None, lon_range: list = None,
                      product: str = None, field: str = None, chunk_time: int = 12) -> xr.DataArray:
    dataset = require_kind(dataset, ("timeseries", "index"), "load_regular_data",
                           "load_climatology_data")
    if get_dataset_cfg(dataset)["kind"] == "index":
        given = [name for name, value in (("depth_range", depth_range), ("lat_range", lat_range),
                                          ("lon_range", lon_range)) if value is not None]
        if given:
            raise ValueError(
                f"'{dataset}' is an ungridded index series with no spatial dimensions, so "
                f"{', '.join(given)} cannot be applied. Select the region with 'variable' "
                f"(e.g. variable='NINO3.4'), or read a gridded dataset to cut out a box."
            )
    da = open_dataset_array(dataset, variable, model=model, field=field, product=product,
                            chunk_time=chunk_time)

    if time_range is not None:
        if "time" not in da.coords:
            raise ValueError(f"'{variable}' is a fixed field with no time axis; drop time_range.")
        da = da.sel(time=slice(time_range[0], time_range[1]))
        if da.sizes.get("time", 0) == 0:
            raise ValueError(f"time_range {time_range} selects no steps from {dataset}.")

    da = subset_depth(da, depth_range)
    return subset_spatial(da, lat_range, lon_range)


# ---------------------------------------------------------------------------
# 2. Climatologies
# ---------------------------------------------------------------------------

@register_tools(
    name="load_climatology_data",
    description="Load a subset of an observational climatology.\n"
        "This is the entry point for a mean state, a seasonal cycle or a vertical structure, and for the "
        "reference field a model is evaluated against. A climatology is already an average over many "
        "years, so it has no time axis to take a trend along.",
    args={
        "dataset": "str, required. The climatology's dataset id.",
        "variable": "str, required for most datasets. The variable name that dataset uses; common aliases "
            "are accepted.",
        "time_slice": "str, optional, default 'annual'. 'annual' for a single all-year field, which "
            "reaches the deepest levels; 'monthly' for twelve fields, adding a 'month' dimension; "
            "'seasonal' for four, adding a 'season' dimension. The monthly and seasonal fields are "
            "shallower than the annual one.",
        "months": "list[int], optional. Select particular months from a monthly climatology, e.g. [1, 2, "
            "12].",
        "period": "str, optional. Which climate normal to read, where the variable offers more than one.",
        "field": "str, optional. Which statistical field to read, where the product publishes several "
            "(the objectively analysed value, the raw mean, the standard deviation, the number of "
            "observations). Default is the analysed value.",
        "version": "str, optional. Product version, where more than one exists.",
        "depth_range": "list[float, float], optional. [min, max] in metres, positive downward.",
        "lat_range": "list[float, float], optional. [min, max] in degrees North.",
        "lon_range": "list[float, float], optional. [min, max] in degrees East, either convention; a "
            "window crossing the dateline is handled."
    },
    returns="xr.DataArray with dimensions (depth, lat, lon), plus 'month' or 'season' for the "
        "corresponding slices, carrying units and the provenance attributes the weighting tools need."
)
def load_climatology_data(dataset: str, variable: str = None, time_slice: str = "annual",
                          months: list = None, period: str = None, field: str = None,
                          version: str = None, depth_range: list = None,
                          lat_range: list = None, lon_range: list = None) -> xr.DataArray:
    dataset = require_kind(dataset, "climatology", "load_climatology_data", "load_regular_data")
    da = open_dataset_array(
        dataset, variable, time_slice=time_slice, period=period, field=field, version=version
    )

    if months is not None:
        if "month" not in da.coords:
            raise ValueError(
                "This field has no 'month' dimension; for WOA23 pass time_slice='monthly' to use months."
            )
        da = da.sel(month=months)

    da = subset_depth(da, depth_range)
    return subset_spatial(da, lat_range, lon_range)


# ---------------------------------------------------------------------------
# 3. Named regions
# ---------------------------------------------------------------------------

SUPPORTED_WATER_BODIES = ['Arctic Ocean', 'SOUTHERN OCEAN', 'North Atlantic Ocean', 'North Pacific Ocean', 'South Pacific Ocean', 'INDIAN OCEAN', 'South Atlantic Ocean', 'Black Sea', 'Philippine Sea', 'Great Barrier Reef', 'Tasman Sea', 'Bay of Bengal', 'South China Sea', 'Arabian Sea', 'Beaufort Sea', 'Caribbean Sea', 'Gulf of Mexico', 'Labrador Sea', 'Hudson Bay', 'Caspian Sea', 'Baffin Bay', 'Gulf of Alaska', 'Red Sea', 'Sea of Okhotsk', 'Weddell Sea', 'Persian Gulf', 'Celebes Sea', 'Sulu Sea', 'Norwegian Sea', 'Greenland Sea', 'Banda Sea', 'Luzon Strait', 'Bay of Biscay', 'Mozambique Channel', 'Gulf of Guinea', 'Scotia Sea', 'Baltic Sea', 'Barents Sea', 'North Sea', 'Bristol Channel', 'Inner Seas', 'Irish Sea', 'Java Sea', 'Andaman Sea', 'Yellow Sea', 'East China Sea', 'Chukchi Sea', 'Bahía de Campeche', 'Arafura Sea', 'Timor Sea', 'Gulf of Thailand', 'Laccadive Sea', 'Bellingshausen Sea', 'Amundsen Sea', 'Davis Strait', 'Kara Sea', 'Laptev Sea', 'Drake Passage', 'Great Australian Bight', 'Gulf of Aden', 'Gulf of Oman', 'Tyrrhenian Sea', 'White Sea', 'Sea of Azov', 'James Bay', 'Gulf of Carpentaria', 'Golfo de California', 'Hudson Strait', 'Gulf of Finland', 'English Channel', 'Gulf of Honduras', 'Golfe du Lion', 'Gulf of Bothnia', 'Adriatic Sea', 'Straits of Florida', 'Bay of Plenty', 'The North Western Passages', 'Queen Charlotte Sound', 'Ionian Sea', 'Gulf of Saint Lawrence', 'Molucca Sea', 'Bismarck Sea', 'Solomon Sea', 'Gulf of Tonkin', 'Bay of Fundy', 'Río de la Plata', 'Strait of Malacca', 'Strait of Singapore', 'Makassar Strait', 'Ceram Sea', 'Inner Sea', 'Taiwan Strait', 'Cook Inlet', 'Bristol Bay', 'Melville Bay', 'Gulf of Maine', 'Chesapeake Bay', 'Strait of Gibraltar', 'Golfo San Jorge', 'Balearic Sea', 'Aegean Sea', 'Gulf of Mannar', 'Shelikhova Gulf', 'Bo Hai', 'Gulf of Kutch', 'Golfo de Panamá', 'Amundsen Gulf', 'Viscount Melville Sound', 'Ungava Bay', 'Yucatan Channel', 'Bering Sea', 'East Siberian Sea', 'Lincoln Sea', 'Bass Strait', 'Franklin Bay', 'Dixon Entrance', "M'Clure Strait", 'Boca Grande', 'Baía de Marajó', 'Golfo San Matías', 'Golfo Corcovado', 'Wrigley Gulf', 'Sulzberger Bay', 'Ronne Entrance', 'Bransfield Strait', 'McMurdo Sound', 'Baía de São Marcos', 'Cumberland Sound', 'Frobisher Bay', 'Disko Bay', 'Hamilton Inlet', 'Storfjorden', 'Vestfjorden', 'Skagerrak', 'Karskiye Strait', "Vil'kitskogo Strait", 'Sognefjorden', 'Trondheimsfjorden', 'Kattegat', 'Gulf of Khambhät', 'Palk Strait', 'Kronotskiy Gulf', 'Uda Bay', 'Uchiura Bay', 'Tsugaru Strait', 'Tatar Strait', 'Shark Bay', 'Gulf of Gabès', 'Mackenzie Bay', 'Norton Sound', 'Kotzebue Sound', 'Gulf of Boothia', 'Kane Basin', 'Foxe Basin', 'Golfo de Tehuantepec', 'Golfo de Guayaquil', 'Gulf of Masira', 'Gulf of Yana', 'Dmitriy Laptev Strait', 'La Pérouse Strait', 'East Korea Bay', 'Qiongzhou Strait', 'Cook Strait', 'Torres Strait', 'Geographe Bay', 'Gulf of Papua', 'Gulf of Olen‘k', 'Gulf of Riga', 'Waddenzee', 'Gulf of Sidra', 'Sea of Crete', 'Ligurian Sea', 'Alboran Sea', 'Bight of Benin', 'Bight of Biafra', 'Strait of Belle Isle', 'Prince William Sound', 'Estrecho de Magellanes', 'Bahía Blanca', 'Golfo de Penas', 'Bahía Grande', 'Marguerite Bay', 'Prydz Bay', 'Vincennes Bay', 'Porpoise Bay', 'Davis Sea', 'Lützow-Holm Bay', 'Antongila Bay', 'Hangzhou Bay', 'Karaginskiy Gulf', 'Gulf of Kamchatka', 'Bosporus', 'Bab el Mandeb', 'Gulf St. Vincent', 'Golfo de Urabá', 'Baia de Maputo', 'St. Helena Bay', 'Boknafjorden', 'Joseph Bonaparte Gulf', 'Gulf of Martaban', 'Gulf of Sakhalin', 'Kangertittivaq', 'Saint Lawrence River', 'Chaun Bay', 'Khatanga Gulf', 'Garabogaz Bay', 'Gulf of Suez', 'Gulf of Ob', 'Yenisey Gulf', 'Dardanelles', 'Sea of Marmara', 'Gulf of Aqaba', 'Lago de Maracaibo', 'Denmark Strait', 'Massachusetts Bay', 'Delaware Bay', 'Long Island Sound', 'Albemarle Sound', 'Pamlico Sound', 'Smith Sound', 'Queen Charlotte Strait', 'Baird Inlet', 'Husky Lakes', 'Lake Pontchartrain', 'Stettiner Haff', 'Darnley Bay', 'Prince of Wales Strait', 'Minto Inlet', 'Richard Collinson Inlet', 'Prince ALbert Sound', 'Liddon Gulf', 'Wynniatt Bay', 'Hadley Bay', 'Bathurst Inlet', 'Goldsmith Channel', 'Sherman Basin', 'Wager Bay', 'Fury and Hecla Strait', 'Jones Sound', 'Eclipse Sound', 'Hall Basin', 'Uummannaq Fjord', 'Peacock Sound', 'George VI Sound', 'Seno de Skyring', 'Seno Otway', 'Bay Inútil', 'Lagoa dos Patos', 'Amazon River', "Bras d'Or Lake", 'Kaliningrad', 'Matochkin Shar Strait', 'Ozero Mogotoyevo', 'Guba Gusinaya', 'Bali Sea', 'Davao Gulf', 'Halmahera Sea', 'Selat Bali', 'Gulf of Tomini', 'Flores Sea', 'Sibuyan Sea', 'Selat Dampier', 'Gulf of Buli', 'Gulf of Kau', 'Bohol Sea', 'Surigao Strait', 'Ragay Gulf', 'Samar Sea', 'Tayabas Bay', 'Leyte Gulf', 'Visayan Sea', 'Savu Sea', 'Hecate Strait', 'Cordova Bay', 'Murchison Sound', 'Robeson Channel', 'Kennedy Channel', 'Canal Perigoso', 'Canal do Sul', 'Canal do Norte', 'Øresund', 'Mecklenburger Bucht', 'Sargasso Sea', 'San Francisco Bay', 'Monterey Bay', 'Yangtze River', 'Columbia River', 'Salish Sea', "Gulf of Anadyr'", 'Mediterranean Sea', 'Ross Sea', 'Coral Sea', 'Sea of Japan', 'Korea Strait']

_GEOLOCATOR = None


def get_geolocator() -> Geolocator:
    """Build the Geolocator once per process; its caches take ~10 s to read."""
    global _GEOLOCATOR
    if _GEOLOCATOR is not None:
        return _GEOLOCATOR

    if not os.path.exists(PATHS_CONFIG):
        raise FileNotFoundError(f"Config file not found: {PATHS_CONFIG}")
    with open(PATHS_CONFIG, "r") as f:
        cfg = yaml.safe_load(f)
    cache_dir, natural_earth_path = cfg.get("cache_dir"), cfg.get("natural_earth_path")
    if not cache_dir or not natural_earth_path:
        raise ValueError("Both 'cache_dir' and 'natural_earth_path' must be set in configs/default_paths.yaml")
    for path, label in ((cache_dir, "cache_dir"), (natural_earth_path, "natural_earth_path")):
        if not os.path.isdir(path):
            raise FileNotFoundError(f"{label} does not exist: {path}")

    _GEOLOCATOR = Geolocator(GEOLOCATOR_LATS, GEOLOCATOR_LONS, cache_dir, natural_earth_path)
    return _GEOLOCATOR


def build_region_mask(geofeature: str, loc_names) -> xr.DataArray:
    """
    Union the Natural Earth masks for one or more named regions into a
    1-degree (lat, lon) boolean field in the -180~180 convention.
    """
    geolocator = get_geolocator()
    names = [loc_names] if isinstance(loc_names, str) else list(loc_names)

    combined = None
    for name in names:
        try:
            mask = geolocator.get_mask_by_geofeature(geofeature, name)
        except Exception as exc:
            raise RuntimeError(f"Failed to get mask for {geofeature}='{name}': {exc}")
        combined = mask if combined is None else (combined | mask)

    # The cache is stored as (lon, lat) on the 0~360 grid the geometries use.
    mask_da = xr.DataArray(
        combined.T,
        coords={"lat": GEOLOCATOR_LATS, "lon": ((GEOLOCATOR_LONS + 180) % 360) - 180},
        dims=["lat", "lon"],
    )
    return mask_da.sortby("lon")


@register_tools(
    name="load_data_with_geolocator",
    description="Load data for a NAMED geographic region instead of a coordinate box, using Natural Earth "
        "boundaries.\n"
        "Use it when the question names a sea, a bay or an ocean: the outline is more faithful than any "
        "box, and no coordinates have to be guessed. For a region that is not in the supported list, fall "
        "back to a loader with explicit ranges. To split a field you already have by ocean basin, use "
        "select_basin instead — this tool masks while reading, so it only works on a stored variable.",
    args={
        "dataset": "str, required. The dataset id to read from.",
        "variable": "str, required for most datasets. The variable name that dataset uses.",
        "model": "str, required for multi-model datasets.",
        "geofeature": "str, optional, default 'water_body'. Which kind of boundary to use: 'water_body', "
            "'country', 'state' or 'continent'. Ocean analysis almost always wants 'water_body'.",
        "loc_name": "str or list[str], required. Name(s) of the region, case-sensitive, from the supported "
            "list below. Pass a list to merge several regions into one field for a single calculation; "
            "make separate calls when they are to be compared.\n"
            f"Supported water bodies: {SUPPORTED_WATER_BODIES}",
        "time_range": "list[str, str], optional. Time series datasets only: ['start', 'end'] in ISO 8601.",
        "depth_range": "list[float, float], optional. [min, max] in metres, positive downward.",
        "time_slice": "str, optional, default 'annual'. Climatologies only: 'annual', 'monthly' or "
            "'seasonal'.",
        "months": "list[int], optional. Select particular months from a monthly climatology.",
        "period": "str, optional. Which climate normal to read.",
        "field": "str, optional. Which statistical field to read.",
        "version": "str, optional. Product version, where more than one exists."
    },
    returns="xr.DataArray with the dimensions the dataset would otherwise return and cells outside the "
        "region set to NaN. A regular grid keeps its full coordinates; a curvilinear grid keeps its index "
        "dimensions and is trimmed to the region's bounding block. Average it with an area-weighted mean, "
        "which picks the right cell areas for the grid on its own."
)
def load_data_with_geolocator(dataset: str, variable: str = None, model: str = None,
                              geofeature: str = "water_body", loc_name=None,
                              time_range: list = None, depth_range: list = None,
                              time_slice: str = "annual", months: list = None,
                              period: str = None, field: str = None,
                              version: str = None) -> xr.DataArray:
    if loc_name is None:
        raise ValueError("loc_name is required; pass a region name from the supported list.")

    dataset = resolve_dataset(dataset)
    kind = get_dataset_cfg(dataset)["kind"]
    if kind not in ("timeseries", "climatology"):
        raise ValueError(
            f"'{dataset}' is a '{kind}' dataset and has no grid to mask. Read it with "
            f"load_regular_data, which takes its region as an argument."
        )

    da = open_dataset_array(
        dataset, variable, model=model, time_slice=time_slice, period=period,
        field=field, version=version,
    )

    if time_range is not None:
        if "time" not in da.coords:
            raise ValueError(f"'{dataset}' has no time axis; use time_slice/months instead of time_range.")
        da = da.sel(time=slice(time_range[0], time_range[1]))
    if months is not None:
        if "month" not in da.coords:
            raise ValueError("This field has no 'month' dimension; for WOA23 pass time_slice='monthly'.")
        da = da.sel(month=months)

    da = subset_depth(da, depth_range)

    mask_1deg = build_region_mask(geofeature, loc_name)
    grid_type = da.attrs.get("grid_type", "regular")
    lat_coord = da.attrs.get("lat_coord", "lat")
    lon_coord = da.attrs.get("lon_coord", "lon")
    y_dim, x_dim = da.attrs.get("y_dim"), da.attrs.get("x_dim")
    attrs = dict(da.attrs)

    mask = mask_1deg_to_grid(mask_1deg, da, grid_type, lat_coord, lon_coord)
    da = da.where(mask)
    if grid_type == "curvilinear":
        da = trim_to_mask(da, mask, y_dim, x_dim)

    attrs["region"] = loc_name if isinstance(loc_name, str) else ", ".join(loc_name)
    attrs["region_geofeature"] = geofeature
    da.attrs = attrs
    return da


# ---------------------------------------------------------------------------
# 4. Grid quantities
# ---------------------------------------------------------------------------

WGS84 = Geod(ellps="WGS84")

# The variable each observational product's grid is read from
GRID_TEMPLATE_VARIABLE = {
    "woa23": "temperature",
    "glodap": "temperature",
    "seanoe_mld": "mld_dr003",
    "ersst": "sst",
}


def bounds_from_centres(centres):
    """
    Infer 1-D cell edges from cell centres, for the products that publish no
    bounds of their own (GLODAP, SEANOE, ERSSTv5). Interior edges are midpoints
    and the two outer edges extend by the adjacent spacing.
    """
    centres = np.asarray(centres, dtype="float64")
    if centres.size < 2:
        raise ValueError("Need at least two cell centres to infer bounds.")
    midpoints = (centres[:-1] + centres[1:]) / 2
    edges = np.empty((centres.size, 2))
    edges[1:, 0] = midpoints
    edges[:-1, 1] = midpoints
    edges[0, 0] = centres[0] - (centres[1] - centres[0]) / 2
    edges[-1, 1] = centres[-1] + (centres[-1] - centres[-2]) / 2
    return edges


def cell_areas_from_bounds(lat_bnds, lon_bnds):
    """
    Areas of a regular lat/lon grid on the WGS84 ellipsoid, equivalent to
    MATLAB's areaquad.
    """
    lat_bnds = np.asarray(lat_bnds, dtype="float64")
    lon_bnds = np.asarray(lon_bnds, dtype="float64")
    widths = lon_bnds[:, 1] - lon_bnds[:, 0]
    uniform = np.allclose(widths, widths[0])

    areas = np.empty((lat_bnds.shape[0], lon_bnds.shape[0]))
    for i, (lower, upper) in enumerate(lat_bnds):
        # A pole-touching cell is degenerate for the polygon solver.
        lower = np.clip(lower, -89.99, 89.99)
        upper = np.clip(upper, -89.99, 89.99)
        if uniform:
            lon1, lon2 = lon_bnds[0]
            area, _ = WGS84.polygon_area_perimeter(
                [lon1, lon2, lon2, lon1], [lower, lower, upper, upper]
            )
            areas[i, :] = abs(area)
        else:
            for j, (lon1, lon2) in enumerate(lon_bnds):
                area, _ = WGS84.polygon_area_perimeter(
                    [lon1, lon2, lon2, lon1], [lower, lower, upper, upper]
                )
                areas[i, j] = abs(area)
    return areas


def cell_areas_from_vertices(lat_vertices, lon_vertices):
    """
    Areas of a curvilinear grid from its per-cell corner vertices, shaped
    (ny, nx, 4).
    """
    lat_vertices = np.asarray(lat_vertices, dtype="float64")
    lon_vertices = np.asarray(lon_vertices, dtype="float64")
    ny, nx = lat_vertices.shape[:2]
    areas = np.empty((ny, nx))
    for j in range(ny):
        for i in range(nx):
            area, _ = WGS84.polygon_area_perimeter(lon_vertices[j, i], lat_vertices[j, i])
            areas[j, i] = abs(area)
    return areas


def area_dataarray(areas, coords, dims, source, **attrs):
    """Wrap a computed area field with the attributes every caller expects."""
    return xr.DataArray(
        areas, dims=dims, coords=coords, name="areacello",
        attrs={
            "units": "m2",
            "long_name": "grid cell area",
            "description": "Area of each grid cell",
            "area_source": source,
            "method": "WGS84 ellipsoid (equivalent to MATLAB areaquad)",
            "ellipsoid": "WGS84 (a=6378137.0 m, f=1/298.257223563)",
            **attrs,
        },
    )


def default_1deg_areas():
    """The 1-degree WGS84 grid shared by WOA23, GLODAP and SEANOE."""
    lat = np.arange(-89.5, 90.0, 1.0)
    lon = np.arange(-179.5, 180.0, 1.0)
    areas = cell_areas_from_bounds(bounds_from_centres(lat), bounds_from_centres(lon))
    return area_dataarray(
        areas,
        {"lat": ("lat", lat, {"units": "degrees_north", "long_name": "latitude"}),
         "lon": ("lon", lon, {"units": "degrees_east", "long_name": "longitude"})},
        ["lat", "lon"], "1-degree WGS84 grid", grid_resolution="1 degree",
    )


def model_areas(dataset, model, variable, lat_range, lon_range):
    """
    Areas for a CMIP6 model, preferring the model's own published cell areas.
    """
    experiment = dataset.replace("cmip6_", "")
    others = [e for e in cmip6_experiments() if e != experiment]
    attempts = []

    for candidate in [experiment] + others:
        try:
            areas = open_cmip6(candidate, model, "areacello")
        except (FileNotFoundError, KeyError, ValueError, OSError) as exc:
            attempts.append(f"cmip6_{candidate}: {exc.__class__.__name__}")
            continue
        borrowed = areas.attrs.get("experiment", candidate) != experiment
        areas.attrs["area_source"] = (
            f"areacello published by cmip6_{candidate}"
            + (" (borrowed: this experiment publishes none)" if borrowed else "")
        )
        return subset_spatial(areas, lat_range, lon_range)

    # Nothing published, so derive the areas from the file's own grid geometry.
    grid_file, info = open_cmip6_grid_file(experiment, model, variable)
    model_cfg = info["model_cfg"]
    template = open_cmip6(experiment, model, variable)
    source = (
        f"computed from cell bounds in {os.path.basename(info['source_file'])} "
        f"(no areacello published; tried {attempts}). These are whole-cell areas on the WGS84 ellipsoid, so "
        f"they run about 0.5% from a model's own areas, which are usually spherical, and they cannot "
        f"reproduce the partial coastal cells of a regridded product such as GFDL-ESM4, whose published "
        f"areacello holds only the ocean part of each cell"
    )

    try:
        if model_cfg["grid_type"] == "regular":
            if "lat_bnds" not in grid_file or "lon_bnds" not in grid_file:
                raise ValueError(
                    f"{model} publishes no areacello for any experiment and "
                    f"{os.path.basename(info['source_file'])} carries no lat_bnds/lon_bnds, so cell areas "
                    f"cannot be computed. Supply an areacello file for this model."
                )
            areas = cell_areas_from_bounds(grid_file["lat_bnds"].values, grid_file["lon_bnds"].values)
            areas_da = area_dataarray(
                areas,
                {"lat": ("lat", grid_file["lat"].values), "lon": ("lon", grid_file["lon"].values)},
                ["lat", "lon"], source,
            )
            # The bounds are in the file's native longitude convention.
            areas_da = sort_ascending(to_lon180(areas_da, "lon"), ["lat"])
            areas_da = snap_coords(areas_da)
        else:
            corners = {
                "lat": next((n for n in ("vertices_latitude", "lat_bnds") if n in grid_file), None),
                "lon": next((n for n in ("vertices_longitude", "lon_bnds") if n in grid_file), None),
            }
            if not corners["lat"] or not corners["lon"] or grid_file[corners["lat"]].ndim != 3:
                raise ValueError(
                    f"{model} publishes no areacello for any experiment and "
                    f"{os.path.basename(info['source_file'])} carries no per-cell corner vertices, so cell "
                    f"areas cannot be computed on this curvilinear grid. Supply an areacello file."
                )
            areas = cell_areas_from_vertices(
                grid_file[corners["lat"]].values, grid_file[corners["lon"]].values
            )
            y_dim, x_dim = model_cfg["y_dim"], model_cfg["x_dim"]
            lat_coord, lon_coord = model_cfg["lat_coord"], model_cfg["lon_coord"]
            areas_da = area_dataarray(
                areas,
                {lat_coord: template[lat_coord], lon_coord: template[lon_coord]},
                [y_dim, x_dim], source,
            )
    finally:
        grid_file.close()

    for key in ("grid_type", "y_dim", "x_dim", "lat_coord", "lon_coord", "model", "experiment"):
        if key in template.attrs:
            areas_da.attrs[key] = template.attrs[key]

    # Cell bounds cover the whole grid, land included, so the derived areas have
    # to be masked to ocean to mean the same thing as a published areacello.
    try:
        sea_floor = open_cmip6(experiment, model, "deptho")
        areas_da = areas_da.where(sea_floor.notnull())
        areas_da.attrs["area_source"] += ", masked to ocean with deptho"
    except (FileNotFoundError, KeyError, ValueError, OSError):
        areas_da.attrs["area_source"] += (
            ", NOT masked to ocean because this model publishes no deptho, so land cells are included"
        )
    return subset_spatial(areas_da, lat_range, lon_range)


def observational_areas(dataset, lat_range, lon_range):
    """Areas on an observational product's own grid, from its published cell
    bounds where it has them and from its cell spacing where it does not."""
    if dataset == "woa23":
        path = woa23_files(GRID_TEMPLATE_VARIABLE[dataset], "annual")["files"][0][1]
        with xr.open_dataset(path, decode_times=False) as ds:
            lat = ds["lat"].values.astype("float64").round(6)
            lon = ds["lon"].values.astype("float64").round(6)
            lat_bnds = ds["lat_bnds"].values.astype("float64")
            lon_bnds = ds["lon_bnds"].values.astype("float64")
        source = f"WOA23 lat_bnds/lon_bnds ({os.path.basename(path)})"
    else:
        template = open_dataset_array(dataset, GRID_TEMPLATE_VARIABLE[dataset])
        lat = template["lat"].values
        lon = template["lon"].values
        lat_bnds = bounds_from_centres(lat)
        lon_bnds = bounds_from_centres(lon)
        source = f"computed from {dataset} cell centres (the product publishes no cell bounds)"

    areas = cell_areas_from_bounds(lat_bnds, lon_bnds)
    areas_da = area_dataarray(
        areas,
        {"lat": ("lat", lat, {"units": "degrees_north", "long_name": "latitude"}),
         "lon": ("lon", lon, {"units": "degrees_east", "long_name": "longitude"})},
        ["lat", "lon"], source, source_dataset=dataset,
    )
    return subset_regular(areas_da, lat_range, lon_range)


@register_tools(
    name="load_areas",
    description="Load the area of each grid cell in square metres.\n"
        "You usually do NOT need this to take an average: an area-weighted mean assembles the right "
        "areas for whatever grid the data is on and masks them to the valid cells. Reach for it when "
        "you want the area field itself — to report the size of a region, to turn a per-square-metre "
        "quantity into a total, or to weight something by hand.\n"
        "Say which grid you mean, because the areas must match the data being weighted cell for cell. "
        "A model's own published cell areas are used where they exist; areas that have to be computed "
        "are computed on the WGS84 ellipsoid.",
    args={
        "lat_range": "list[float, float], optional. [min, max] in degrees North. Default is global.",
        "lon_range": "list[float, float], optional. [min, max] in degrees East, either convention; a "
            "window crossing the dateline is handled. Default is global.",
        "dataset": "str, optional. Which dataset's grid to describe. If None, the shared 1-degree "
            "observational grid is returned, which is not correct for products on any other grid.",
        "model": "str, required when the dataset has several models.",
        "variable": "str, optional, default 'thetao'. Only consulted when a model publishes no cell areas "
            "of its own: it says which variable's file to take the cell corner bounds from."
    },
    returns="xr.DataArray of cell areas in m2, already masked to ocean, on the requested grid and "
        "matching what the loaders return for it so the two multiply directly. attrs['area_source'] "
        "records where the areas came from."
)
def load_areas(lat_range: list = None, lon_range: list = None, dataset: str = None,
               model: str = None, variable: str = "thetao") -> xr.DataArray:
    if dataset is None:
        return subset_regular(default_1deg_areas(), lat_range, lon_range)

    dataset = resolve_dataset(dataset)
    if dataset.startswith("cmip6_"):
        if model is None:
            raise ValueError(
                f"'{dataset}' needs a model: one of {sorted(get_dataset_cfg(dataset)['models'])}"
            )
        return model_areas(dataset, model, variable, lat_range, lon_range)

    if dataset not in GRID_TEMPLATE_VARIABLE:
        raise ValueError(f"'{dataset}' has no spatial grid, so it has no cell areas.")
    return observational_areas(dataset, lat_range, lon_range)


def model_volumes(dataset, model, variable, time_range, depth_range, lat_range, lon_range,
                  strict_time: bool = True):
    """
    Volumes for a CMIP6 model, preferring the model's own published cell volumes.
    """
    experiment = dataset.replace("cmip6_", "")
    others = [e for e in cmip6_experiments() if e != experiment]
    attempts = []

    volumes = None
    for candidate in [experiment] + others:
        try:
            volumes = open_cmip6(candidate, model, "volcello")
        except (FileNotFoundError, KeyError, ValueError, OSError) as exc:
            attempts.append(f"cmip6_{candidate}: {exc.__class__.__name__}")
            continue
        borrowed = volumes.attrs.get("experiment", candidate) != experiment
        volumes.attrs["volume_source"] = (
            f"volcello published by cmip6_{candidate}"
            + (" (borrowed: this experiment publishes none)" if borrowed else "")
        )
        break

    if volumes is None:
        areas = load_areas(dataset=dataset, model=model, variable=variable)
        _, depth_bnds = load_depth_levels(dataset=dataset, model=model, variable=variable)
        thickness = xr.DataArray(
            (depth_bnds[:, 1] - depth_bnds[:, 0]).values,
            coords={"depth": depth_bnds["depth"].values}, dims=["depth"],
        )
        note = f"computed as cell area x layer thickness (no volcello published; tried {attempts})"

        try:
            sea_floor = open_cmip6(experiment, model, "deptho")
            upper = xr.DataArray(
                depth_bnds[:, 1].values, coords={"depth": depth_bnds["depth"].values}, dims=["depth"]
            )
            lower = xr.DataArray(
                depth_bnds[:, 0].values, coords={"depth": depth_bnds["depth"].values}, dims=["depth"]
            )
            thickness = np.minimum(upper, sea_floor) - lower
            thickness = thickness.where(thickness > 0, 0.0)
            note += ", with the bottom layer clipped at deptho"
        except (FileNotFoundError, KeyError, ValueError, OSError):
            note += ", WITHOUT clipping the bottom layer (this model publishes no deptho), so cells at the sea floor are overestimated"

        volumes = areas * thickness
        volumes = volumes.transpose("depth", *[d for d in volumes.dims if d != "depth"])
        volumes.name = "volcello"
        volumes.attrs = dict(areas.attrs)
        volumes.attrs.update({"units": "m3", "long_name": "ocean grid-cell volume",
                              "volume_source": note})
        volumes.attrs.pop("area_source", None)

    if time_range is not None:
        if "time" not in volumes.coords:
            # A fixed field has no window to take.
            if strict_time:
                raise ValueError(
                    f"volcello for {model} is a fixed field with no time axis; drop time_range. "
                    f"(Only GFDL-ESM4 publishes a time-varying volcello.)"
                )
        else:
            volumes = volumes.sel(time=slice(time_range[0], time_range[1]))

    volumes = subset_depth(volumes, depth_range)
    return subset_spatial(volumes, lat_range, lon_range)


@register_tools(
    name="load_volumes_for_cmip6",
    description="Load the volume of each grid cell of a model in cubic metres.\n"
        "You usually do NOT need this to take a volume-weighted average, which assembles the right "
        "volumes itself. Reach for it to report the volume of a region, to turn a per-cubic-metre "
        "quantity into a total, or to weight something by hand.\n"
        "A model's own published cell volumes are used where they exist; otherwise volumes are rebuilt "
        "from cell areas and layer thicknesses.",
    args={
        "dataset": "str, required. The experiment's dataset id.",
        "model": "str, required. Which model's grid to describe.",
        "time_range": "list[str, str], optional. ['start', 'end'] in ISO 8601. Only accepted for a model "
            "that publishes a time-varying cell volume; for a fixed field it is an error.",
        "depth_range": "list[float, float], optional. [min, max] in metres, positive downward.",
        "lat_range": "list[float, float], optional. [min, max] in degrees North.",
        "lon_range": "list[float, float], optional. [min, max] in degrees East, either convention.",
        "variable": "str, optional, default 'thetao'. Only consulted when the volumes have to be rebuilt: "
            "it says which variable's vertical grid to use."
    },
    returns="xr.DataArray of cell volumes in m3, masked to ocean, on the model's own grid and matching "
        "what the loaders return for it. Where a model publishes its cell volume monthly the result keeps "
        "a time axis; reduce over time, or pass a time_range, before using it as a static weight. "
        "attrs['volume_source'] records whether the volumes are the model's own or were rebuilt."
)
def load_volumes_for_cmip6(dataset: str, model: str, time_range: list = None, depth_range: list = None,
                           lat_range: list = None, lon_range: list = None,
                           variable: str = "thetao") -> xr.DataArray:
    dataset = resolve_dataset(dataset)
    if not dataset.startswith("cmip6_"):
        raise ValueError(
            f"'{dataset}' is not a CMIP6 experiment and publishes no cell volumes. "
            f"Build volumes for it from load_areas and load_depth_levels instead."
        )
    if model is None:
        raise ValueError(f"'{dataset}' needs a model: one of {sorted(get_dataset_cfg(dataset)['models'])}")

    return model_volumes(dataset, model, variable, time_range, depth_range, lat_range, lon_range)


@register_tools(
    name="load_depth_levels",
    description="Load a dataset's depth levels and the boundaries of each layer, in metres.\n"
        "You usually do NOT need this to average or integrate over depth, which applies the correct layer "
        "thicknesses itself. Reach for it when you want the level definitions themselves — to report which "
        "depths a dataset resolves, to check whether two products share a vertical grid, or to weight "
        "something by hand.\n"
        "Vertical grids are not shared between products, and within a product they can differ by variable "
        "and by time slice, so say which you mean.",
    args={
        "depth_range": "list[float, float], optional. [min, max] in metres, positive downward. Levels are "
            "selected by their CENTRE and their layers returned whole, so the column covered reaches past "
            "the range asked for; check the returned bounds and report those.",
        "dataset": "str, optional, default the shared observational climatology grid. Which dataset's "
            "vertical grid to read.",
        "model": "str, required when the dataset has several models.",
        "variable": "str, optional. Which variable's vertical grid to read, since levels can differ by "
            "variable within a product.",
        "time_slice": "str, optional, default 'annual'. Which time slice's vertical axis to read, for "
            "climatologies whose monthly and seasonal fields are shallower than the annual one.",
        "period": "str, optional. Which climate normal to read the axis from."
    },
    returns="tuple (levels, bounds): the level centres in metres, and the upper and lower boundary of each "
        "layer as an (n, 2) array. Layer thickness is the difference between the two."
)
def load_depth_levels(depth_range: list = None, dataset: str = None, model: str = None,
                      variable: str = None, time_slice: str = "annual",
                      period: str = None) -> tuple:
    dataset = resolve_dataset(dataset or "woa23")

    if dataset.startswith("cmip6_"):
        experiment = dataset.replace("cmip6_", "")
        if model is None:
            raise ValueError(
                f"'{dataset}' needs a model: one of {sorted(get_dataset_cfg(dataset)['models'])}"
            )
        grid_file, info = open_cmip6_grid_file(experiment, model, variable or "thetao")
        model_cfg = info["model_cfg"]
        depth_coord = model_cfg["depth_coord"]
        try:
            if depth_coord not in grid_file.coords:
                raise ValueError(
                    f"'{variable or 'thetao'}' in {model}/{experiment} has no vertical dimension "
                    f"(it is a surface or depth-integrated field), so it has no depth levels."
                )
            bounds_name = grid_file[depth_coord].attrs.get("bounds", f"{depth_coord}_bnds")
            if bounds_name not in grid_file:
                raise ValueError(
                    f"{os.path.basename(info['source_file'])} carries no '{bounds_name}', so layer "
                    f"boundaries cannot be determined for {model}. Only the level centres are available."
                )
            level_scale = depth_scale_from_units(
                grid_file[depth_coord].attrs.get("units"), model_cfg["depth_scale"]
            )
            bounds_scale = depth_scale_from_units(
                grid_file[bounds_name].attrs.get("units"), level_scale
            )
            depth_values = grid_file[depth_coord].values.astype("float64") * level_scale
            bounds_values = grid_file[bounds_name].values.astype("float64") * bounds_scale
            level_units = grid_file[depth_coord].attrs.get("units")
            bounds_units = grid_file[bounds_name].attrs.get("units") or level_units
        finally:
            grid_file.close()
        source = (
            f"{depth_coord} ({level_units}) / {bounds_name} ({bounds_units}) in "
            f"{os.path.basename(info['source_file'])}, converted to metres"
        )

    elif dataset == "woa23":
        resolved = woa23_files(variable or GRID_TEMPLATE_VARIABLE["woa23"], time_slice, period)
        path = resolved["files"][0][1]
        with xr.open_dataset(path, decode_times=False) as ds:
            depth_values = ds["depth"].values.astype("float64")
            bounds_values = ds["depth_bnds"].values.astype("float64")
        source = f"WOA23 depth_bnds ({os.path.basename(path)})"

    elif dataset == "glodap":
        depth_values = np.asarray(get_dataset_cfg("glodap")["depth"]["levels"], dtype="float64")
        bounds_values = bounds_from_centres(depth_values)
        bounds_values[0, 0] = 0.0          # the surface level starts at the surface
        source = "derived as midpoints between levels; GLODAP publishes no depth bounds"

    else:
        raise ValueError(f"'{dataset}' has no vertical dimension, so it has no depth levels.")

    if depth_range is not None:
        keep = (depth_values >= depth_range[0]) & (depth_values <= depth_range[1])
        depth_values = depth_values[keep]
        bounds_values = bounds_values[keep, :]
        if depth_values.size == 0:
            raise ValueError(f"depth_range {depth_range} selects no levels from '{dataset}'.")

    depth_attrs = {
        "standard_name": "depth", "bounds": "depth_bnds", "positive": "down",
        "units": "meters", "axis": "Z", "long_name": "depth", "source": source,
    }
    depth_da = xr.DataArray(
        depth_values,
        coords={"depth": ("depth", depth_values, depth_attrs)},
        dims=["depth"], attrs=depth_attrs,
    )
    depth_bnds_da = xr.DataArray(
        bounds_values,
        coords={"depth": ("depth", depth_values), "nbounds": ("nbounds", [0, 1])},
        dims=["depth", "nbounds"],
        attrs={
            "comment": "depth bounds", "units": "meters", "long_name": "depth bounds",
            "description": "Upper and lower boundaries of each depth level",
            "bounds_source": source,
        },
    )
    return depth_da, depth_bnds_da



# ---------------------------------------------------------------------------
# 5. Grid quantities taken from the data's own provenance
# ---------------------------------------------------------------------------

def dataset_of(da: xr.DataArray) -> tuple:
    """The (dataset, model) an array came from, as recorded by its opener."""
    dataset = da.attrs.get("source_dataset")
    if dataset is None:
        raise ValueError(
            "This array carries no 'source_dataset' attribute, so the grid it belongs to cannot be "
            "determined. Load it with load_regular_data, load_climatology_data or "
            "load_data_with_geolocator, and avoid rebuilding it in a way that drops attrs "
            "(xr.DataArray(...) without attrs=, or arithmetic with keep_attrs off)."
        )
    return dataset, da.attrs.get("model")


def unmatched_labels(source_values, target_values, tolerance):
    """
    Which of `target_values` have no neighbour in `source_values` within tolerance.
    """
    source = np.unique(np.asarray(source_values, dtype="float64"))
    target = np.asarray(target_values, dtype="float64")
    if source.size == 0:
        return target
    if source.size == 1:
        return target[np.abs(target - source[0]) > tolerance]
    position = np.clip(np.searchsorted(source, target), 1, source.size - 1)
    below, above = source[position - 1], source[position]
    nearest = np.where(np.abs(target - below) <= np.abs(target - above), below, above)
    return target[np.abs(target - nearest) > tolerance]


def _align_horizontal(field, da):
    """
    Put a grid field (areas, volumes) onto exactly `da`'s horizontal coordinates.
    """
    if da.attrs.get("grid_type") == "curvilinear":
        y_dim, x_dim = da.attrs["y_dim"], da.attrs["x_dim"]
        lat_coord, lon_coord = da.attrs["lat_coord"], da.attrs["lon_coord"]
        if field.sizes.get(y_dim) != da.sizes.get(y_dim) or field.sizes.get(x_dim) != da.sizes.get(x_dim):
            # A lat/lon box and load_data_with_geolocator both trim the index
            # dimensions, so put the field through the same trim by position.
            y_slice, x_slice = _locate_block(
                field[lat_coord].values, field[lon_coord].values,
                da[lat_coord].values, da[lon_coord].values,
            )
            field = field.isel({y_dim: y_slice, x_dim: x_slice})

        horizontal = {y_dim, x_dim}
        return field.assign_coords({
            c: da[c] for c in field.coords
            if c in da.coords and set(da[c].dims) <= horizontal and set(field[c].dims) <= horizontal
        })

    if "lat" not in da.coords or "lon" not in da.coords:
        raise ValueError(
            f"This array needs 'lat' and 'lon' coordinates to be weighted by its grid, but it has "
            f"{list(da.coords)}. Grid weighting cannot be applied after the horizontal dimensions have "
            f"been aggregated away — weight first, then reduce."
        )

    targets = {
        "lat": np.asarray(da["lat"].values, dtype=float),
        "lon": ((np.asarray(da["lon"].values, dtype=float) + 180.0) % 360.0) - 180.0,
    }
    for name, target in targets.items():
        missing = unmatched_labels(field[name].values, target, 1e-4)
        if missing.size:
            raise ValueError(
                f"{missing.size} of this array's {name} values have no cell on the dataset's own grid "
                f"(first few: {np.round(missing[:5], 4).tolist()}). The array has been regridded, "
                f"interpolated or shifted since it was loaded, and grid weighting needs the dataset's "
                f"original grid. Weight the field before any such step, not after."
            )
    aligned = field.reindex(lat=targets["lat"], lon=targets["lon"], method="nearest", tolerance=1e-4)
    return aligned.assign_coords(lat=da["lat"], lon=da["lon"])


def areas_like(da: xr.DataArray) -> xr.DataArray:
    """
    Cell areas on the grid `da` actually lives on, aligned to its coordinates.
    """
    dataset, model = dataset_of(da)
    areas = _align_horizontal(load_areas(dataset=dataset, model=model), da)
    if da.attrs.get("grid_point") == "velocity":
        areas.attrs["area_source"] = (
            f"{areas.attrs.get('area_source', 'cell areas')}, tracer-cell areas applied to the "
            f"velocity points of the same index (offset by half a cell)"
        )
    return areas


def volumes_like(da: xr.DataArray) -> xr.DataArray:
    """
    Cell volumes in m3 for `da`, on its own grid and its own depth levels.
    """
    dataset, model = dataset_of(da)
    if not dataset.startswith("cmip6_"):
        areas = areas_like(da)
        thickness = layer_thickness_like(da)
        volumes = areas * thickness
        # Binary arithmetic drops attrs, and the caller reports the provenance.
        volumes.name = "volcello"
        volumes.attrs = {
            "units": "m3", "long_name": "grid cell volume",
            "volume_source": f"{areas.attrs.get('area_source', 'cell areas')} x "
                             f"{thickness.attrs.get('bounds_source', 'layer thickness')}",
        }
        return volumes

    time_range = None
    if "time" in da.coords and da["time"].size:
        stamps = np.atleast_1d(da["time"].values)
        time_range = [str(stamps.min())[:10], str(stamps.max())[:10]]
    volumes = model_volumes(dataset, model, "thetao", time_range, None, None, None,
                            strict_time=False)

    if "time" in volumes.dims:
        if "time" in da.dims:
            volumes = volumes.sel(time=da["time"], method="nearest")
            volumes = volumes.assign_coords(time=da["time"])
        else:
            volumes = volumes.isel(time=0, drop=True)
    volumes = _align_horizontal(volumes, da)
    if "depth" in da.coords:
        missing = unmatched_labels(volumes["depth"].values, np.atleast_1d(da["depth"].values), 1e-4)
        if missing.size:
            raise ValueError(
                f"{missing.size} of this array's depth levels are not on {model}'s vertical grid "
                f"(first few: {np.round(missing[:5], 3).tolist()} m), so its cell volumes cannot be matched. "
                f"The array must keep the model's original depth coordinate: a depth subset is fine, "
                f"interpolating or averaging over depth is not."
            )
        volumes = volumes.reindex(depth=np.atleast_1d(da["depth"].values), method="nearest", tolerance=1e-4)
        volumes = volumes.assign_coords(depth=da["depth"])
    volumes.attrs.setdefault("units", "m3")
    return volumes


def _locate_block(full_lat, full_lon, sub_lat, sub_lon):
    """
    Where a trimmed curvilinear block sits inside the model's full grid.
    """
    ny, nx = sub_lat.shape
    corners = np.argwhere(
        np.isclose(full_lat, sub_lat[0, 0]) & np.isclose(full_lon, sub_lon[0, 0])
    )
    for y0, x0 in corners:
        if y0 + ny > full_lat.shape[0] or x0 + nx > full_lat.shape[1]:
            continue
        window = (slice(y0, y0 + ny), slice(x0, x0 + nx))
        if np.allclose(full_lat[window], sub_lat) and np.allclose(full_lon[window], sub_lon):
            return window
    raise ValueError(
        "Could not locate this subset inside the model's full grid, so its cell areas cannot be matched. "
        "Load the areas with load_areas(dataset=..., model=...) and subset both of them the same way."
    )


def layer_bounds_like(da: xr.DataArray) -> xr.DataArray:
    """
    The upper and lower boundary of each of `da`'s layers, in metres.
    """
    if "depth" not in da.coords:
        raise ValueError("This array has no 'depth' coordinate, so it has no layer boundaries.")
    dataset, model = dataset_of(da)
    depths = np.atleast_1d(da["depth"].values).astype("float64")

    pad = DEPTH_MATCH_TOL
    kwargs = {"dataset": dataset, "model": model,
              "depth_range": [float(depths.min()) - pad, float(depths.max()) + pad]}
    if dataset == "woa23":
        variable = da.attrs.get("woa_variable")
        time_slice = da.attrs.get("woa_time_slice")
        if variable is None or time_slice is None:
            raise ValueError(
                f"Cannot tell which WOA23 vertical axis this array belongs to: it is missing the "
                f"'woa_variable' and/or 'woa_time_slice' attributes the loader records. The axes of "
                f"different variables and time slices share their level centres and differ only in how the "
                f"deepest layer closes, so none can stand in for another and guessing would quietly "
                f"misweight the bottom layer. Reload the field with load_climatology_data and keep its "
                f"attributes through the calculation."
            )
        kwargs["variable"] = variable
        kwargs["time_slice"] = time_slice
        kwargs["period"] = da.attrs.get("woa_period")
    levels, bounds = load_depth_levels(**kwargs)

    if levels.size != depths.size or not np.allclose(levels.values, depths, atol=DEPTH_MATCH_TOL):
        raise ValueError(
            f"The depth levels of this array do not match the '{dataset}' vertical grid it claims to come "
            f"from ({depths.size} levels against {levels.size}). Layer boundaries can only be taken from the "
            f"dataset's own level definitions, so the array must keep its original depth coordinate: "
            f"subsetting a depth range is fine, interpolating or averaging over depth is not."
        )
    return bounds.assign_coords(depth=da["depth"])


def layer_thickness_like(da: xr.DataArray) -> xr.DataArray:
    """Thickness in metres of each of `da`'s layers, on `da`'s depth coordinate."""
    bounds = layer_bounds_like(da)
    thickness = (bounds[:, 1] - bounds[:, 0]).astype("float64")
    thickness.name = "layer_thickness"
    thickness.attrs = {
        "units": "m", "long_name": "layer thickness",
        "bounds_source": bounds.attrs.get("bounds_source", ""),
    }
    return thickness


def cell_thickness_like(da: xr.DataArray) -> xr.DataArray:
    """
    The thickness of water in each of `da`'s cells, in metres.

    For a CMIP6 model this is volcello / areacello. Observational 
    products have no partial cells, so their layer bounds are used. 
    """
    dataset, _ = dataset_of(da)
    if not dataset.startswith("cmip6_"):
        return layer_thickness_like(da)

    volumes = volumes_like(da)
    areas = areas_like(da)
    thickness = (volumes / areas.where(areas > 0)).astype("float64")
    thickness.name = "cell_thickness"
    thickness.attrs = {
        "units": "m", "long_name": "cell thickness (volume per unit area)",
        "bounds_source": f"volcello / areacello ({volumes.attrs.get('volume_source', 'cell volumes')}), "
                         f"which accounts for partial cells at the sea floor",
    }
    return thickness


BASIN_NAMES = {
    1: "Atlantic Ocean",
    2: "Pacific Ocean",
    3: "Indian Ocean",
    4: "Mediterranean Sea",
    5: "Baltic Sea",
    6: "Black Sea",
    7: "Red Sea",
    8: "Persian Gulf",
    9: "Hudson Bay",
    10: "Southern Ocean",
    11: "Arctic Ocean",
    12: "Sea of Japan",
    13: "Kara Sea",
    14: "Sulu Sea",
    15: "Baffin Bay",
    16: "East Mediterranean",
    17: "West Mediterranean",
    18: "Sea of Okhotsk",
    19: "Banda Sea",
    20: "Caribbean Sea",
    21: "Andaman Basin",
    22: "North Caribbean",
    23: "Gulf of Mexico",
    24: "Beaufort Sea",
    25: "South China Sea",
    26: "Barents Sea",
    27: "Celebes Sea",
    28: "Aleutian Basin",
    29: "Fiji Basin",
    30: "North American Basin",
    31: "West European Basin",
    32: "Southeast Indian Basin",
    33: "Coral Sea",
    34: "East Indian Basin",
    35: "Central Indian Basin",
    36: "Southwest Atlantic Basin",
    37: "Southeast Atlantic Basin",
    38: "Southeast Pacific Basin",
    39: "Guatemala Basin",
    40: "East Caroline Basin",
    41: "Marianas Basin",
    42: "Philippine Sea",
    43: "Arabian Sea",
    44: "Chile Basin",
    45: "Somali Basin",
    46: "Mascarene Basin",
    47: "Crozet Basin",
    48: "Guinea Basin",
    49: "Brazil Basin",
    50: "Argentine Basin",
    51: "Tasman Sea",
    52: "Atlantic Indian Basin",
    53: "Caspian Sea",
    54: "Sulu Sea II",
    55: "Venezuela Basin",
    56: "Bay of Bengal",
    57: "Java Sea",
    58: "East Indian Atlantic Basin",
    59: "Chiloe",
    60: "Bransfield Strait",
}
BASIN_CODES = {v: k for k, v in BASIN_NAMES.items()}


def _read_woa_mask_table(filename: str) -> pd.DataFrame:
    """Read one of the WOA23 ASCII masks. The first line is a comment and the
    row order is not a raster scan, so the grid is rebuilt from the
    Latitude/Longitude columns rather than from the row order."""
    ds_cfg = get_dataset_cfg("woa23")
    path = os.path.join(ds_cfg["masks_root"], filename)
    if not os.path.exists(path):
        raise FileNotFoundError(f"WOA23 mask file not found: {path}")
    return pd.read_csv(path, comment="#")


@register_tools(
    name="load_basin_mask",
    description="Load the ocean basin each 1-degree cell belongs to, as numeric codes.\n"
        "Codes are given separately at each standard depth level, because a marginal sea's label stops at "
        "its sill and below that the water belongs to the open ocean. Use it to separate basins in a "
        "calculation by hand; to mask a field to a named basin, select_basin does this and carries the "
        "mask onto the field's own grid.",
    args={
        "depth_levels": "str or list[str], optional. Depth level(s) to select, written as the depth with "
            "an 'm' suffix, e.g. '0m', '100m'. Must be exact standard levels. If None, the full 3-D mask "
            "is returned."
    },
    returns="xr.DataArray of basin codes on the 1-degree grid, 2-D (lat, lon) for a single level and 3-D "
        "(depth, lat, lon) otherwise. Land and below-sea-floor cells are missing. Use "
        "attrs['basin_names'] to decode the numbers; the codes distinguish each marginal sea from the "
        "ocean it opens onto, so an ocean taken as its main code alone excludes them."
)
def load_basin_mask(depth_levels=None) -> xr.DataArray:
    df = _read_woa_mask_table("basinmask_01.msk")

    basin_columns = [c for c in df.columns if c.startswith("Basin_")]
    depths = np.array([int(c.replace("Basin_", "").replace("m", "")) for c in basin_columns])
    lats = np.arange(-89.5, 90.0, 1.0)
    lons = np.arange(-179.5, 180.0, 1.0)

    # The table lists only the cells that have ocean, and not in raster order,
    # so the grid is filled by index rather than reshaped.
    lat_idx = np.rint(df["Latitude"].to_numpy() + 89.5).astype(int)
    lon_idx = np.rint(df["Longitude"].to_numpy() + 179.5).astype(int)
    basin_data = np.full((len(depths), len(lats), len(lons)), np.nan)
    basin_data[:, lat_idx, lon_idx] = df[basin_columns].to_numpy(dtype=float).T

    full_mask = xr.DataArray(
        basin_data,
        dims=["depth", "lat", "lon"],
        coords={"depth": depths, "lat": lats, "lon": lons},
        name="basin_mask",
        attrs={
            "long_name": "WOA23 Basin Mask",
            "description": "Ocean basin identification codes",
            "basin_names": BASIN_NAMES,
            "depth_levels": [f"{d}m" for d in depths],
            "valid_range": [1, 60],
            "source": "WOA23 masks/basinmask_01.msk (labels inherited from WOA13)",
        },
    )

    if depth_levels is None:
        return full_mask

    if isinstance(depth_levels, str):
        depth_levels = [depth_levels]
    selected = []
    for level in depth_levels:
        value = int(str(level).replace("m", ""))
        if value not in depths:
            raise ValueError(f"Depth {level} not found. Available: {[f'{d}m' for d in depths]}")
        selected.append(value)

    return full_mask.sel(depth=selected[0] if len(selected) == 1 else selected)


@register_tools(
    name="load_landsea_mask",
    description="Load the land/sea mask on the 1-degree grid, which doubles as a bathymetry index.\n"
        "The value counts how deep the ocean reaches at each cell, so it can exclude land before "
        "averaging or restrict a calculation to cells deep enough for the levels being analysed.",
    args={},
    returns="xr.DataArray of integer bottom level numbers with dims (lat, lon) on the 1-degree grid, "
        "aligned with the observational climatologies.\n"
        "This is NOT a 0/1 field: a value of 1 means land, and ocean is (mask > 1). The value is offset by "
        "one from the depth axis, so (mask - 1) is the number of levels of ocean, i.e. the 1-based index "
        "of the deepest level holding data. Values beyond the climatology's own axis index deeper standard "
        "levels in the trenches; clip them before using the mask as an index."
)
def load_landsea_mask() -> xr.DataArray:
    df = _read_woa_mask_table("landsea_01.msk")

    lats = np.arange(-89.5, 90.0, 1.0)
    lons = np.arange(-179.5, 180.0, 1.0)
    lat_idx = np.rint(df["Latitude"].to_numpy() + 89.5).astype(int)
    lon_idx = np.rint(df["Longitude"].to_numpy() + 179.5).astype(int)

    column = [c for c in df.columns if c.lower().startswith("bottom_standard")][0]
    values = np.zeros((len(lats), len(lons)), dtype=np.int16)
    values[lat_idx, lon_idx] = df[column].to_numpy(dtype=int)

    return xr.DataArray(
        values,
        dims=["lat", "lon"],
        coords={"lat": lats, "lon": lons},
        name="bottom_standard_level",
        attrs={
            "long_name": "WOA23 bottom standard level",
            "description": "1 = land; otherwise (value - 1) is the number of WOA23 standard levels of ocean, "
                "i.e. the 1-based index of the deepest level holding data (verified against the WOA23 "
                "fields). Values above 103 lie below the 5500 m climatology floor.",
            "ocean_test": "value > 1",
            "level_index": "value - 1, clipped to 102",
            "source": "WOA23 masks/landsea_01.msk (labels inherited from WOA13)",
        },
    )
