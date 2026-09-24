# Internal readers: open_cmip6, open_woa23, open_glodap, open_seanoe_mld, open_ersst,
# open_cpc_index. None of them is registered in src/tools/registry, so none is in the
# agent's namespace: every one is reachable through a load_* tool in data_loader.py, which
# is the only entry point the agent has. Registering them too would offer the planning step
# a second, shorter-looking route that skips the subsetting and the dataset-kind checks the
# loaders do.
#
# One opener per data source, each doing a single job: read the raw files
# correctly and normalise them.
#   - resolve the file layout (chunked files, variable prefixes, time-slice
#     numbering, member ids)
#   - correct each product's own metadata defects (non-CF time units, units
#     written under the attribute name `unit`, mislabelled units)
#   - converge on one coordinate contract: longitude -180~180, latitude
#     ascending, depth in metres positive down and ascending, monthly time axes
#     stamped on the first of the month as datetime64
# Subsetting (time, space, named regions) does not happen here; that is
# src/tools/data_loader.py.
#
# Every path, variable prefix, level count and grid convention is read from
# configs/datasets.yaml; nothing is hard-coded here.

import os
import re
import glob
import functools

import numpy as np
import pandas as pd
import xarray as xr
import yaml



CONFIG_PATH = os.path.normpath(
    os.path.join(os.path.dirname(__file__), "..", "..", "configs", "datasets.yaml")
)

SEASON_ORDER = ["winter", "spring", "summer", "autumn"]

# ONI's three-month overlapping seasons map to their centre month: DJF 1950
# is centred on January 1950.
ONI_SEASON_TO_CENTER_MONTH = {
    "DJF": 1, "JFM": 2, "FMA": 3, "MAM": 4, "AMJ": 5, "MJJ": 6,
    "JJA": 7, "JAS": 8, "ASO": 9, "SON": 10, "OND": 11, "NDJ": 12,
}


# ---------------------------------------------------------------------------
# Unit normalisation
#
# The same physical unit is spelled differently in every product: degrees
# Celsius is 'degC' in CMIP6 and ERSSTv5, 'degrees_celsius' in WOA23 and the
# misspelled 'degrees celcius' in GLODAP; practical salinity is '1' in WOA23,
# '0.001' in CMIP6 and an empty string in GLODAP. Units are therefore compared
# after normalisation, never literally.
#
# What normalisation cannot remove is the molar basis: WOA23 and GLODAP report
# per unit mass (umol kg-1), CMIP6 per unit volume (mol m-3). The two convert
# into each other only through in-situ density, so molar_basis leaves that
# judgement to the caller.
# ---------------------------------------------------------------------------

UNIT_EQUIVALENTS = {
    # temperature
    "degc": "degC", "degrees_celsius": "degC", "degrees celcius": "degC",
    "degrees celsius": "degC", "celsius": "degC", "c": "degC", "deg_c": "degC",
    "kelvin": "K", "k": "K", "degk": "K", "degrees_kelvin": "K",
    # concentration, mass basis
    "micromoles_per_kilogram": "umol kg-1", "micro-mol kg-1": "umol kg-1",
    "umol kg-1": "umol kg-1", "umol/kg": "umol kg-1", "μmol/kg": "umol kg-1",
    "µmol/kg": "umol kg-1", "umolkg^-1": "umol kg-1", "micromol/kg": "umol kg-1",
    "mmol kg-1": "mmol kg-1", "mmol/kg": "mmol kg-1",
    # concentration, volume basis
    "mol m-3": "mol m-3", "mol/m3": "mol m-3", "mol m^-3": "mol m-3",
    "molm^-3": "mol m-3", "mol/m^3": "mol m-3",
    "mmol m-3": "mmol m-3", "mmol/m3": "mmol m-3", "mmol/m^3": "mmol m-3",
    "mmolm^-3": "mmol m-3",
    "umol l-1": "mmol m-3", "umol/l": "mmol m-3", "μmol/l": "mmol m-3",
    "µmol/l": "mmol m-3", "micromol/l": "mmol m-3", "umoll^-1": "mmol m-3",
    # dimensionless (practical salinity, pH, saturation state)
    "1": "1", "": "1", "psu": "1", "0.001": "1", "pss-78": "1", "none": "1",
    # length and speed
    "m": "m", "meter": "m", "meters": "m", "metre": "m", "metres": "m",
    "cm": "cm", "centimeter": "cm", "centimeters": "cm",
    "m s-1": "m s-1", "m/s": "m s-1", "ms^-1": "m s-1",
    "cm s-1": "cm s-1", "cm/s": "cm s-1", "cms^-1": "cm s-1",
    "mm s-1": "mm s-1", "mm/s": "mm s-1", "mms^-1": "mm s-1",
}

# Molar basis. Converting between mass and volume basis multiplies or divides
# by seawater density; it is not a decimal conversion.
MASS_BASIS = {"umol kg-1", "mmol kg-1"}
VOLUME_BASIS = {"mol m-3", "mmol m-3"}


def normalize_units(units: str) -> str:
    """Map a units string onto a canonical spelling for comparison."""
    if units is None:
        return None
    key = str(units).strip().lower()
    return UNIT_EQUIVALENTS.get(key, key)


def molar_basis(units: str) -> str:
    """'mass', 'volume' or None — which basis a concentration unit is on."""
    canonical = normalize_units(units)
    if canonical in MASS_BASIS:
        return "mass"
    if canonical in VOLUME_BASIS:
        return "volume"
    return None


# ---------------------------------------------------------------------------
# Reading the registry
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=1)
def load_dataset_config() -> dict:
    """Read configs/datasets.yaml once and cache it."""
    if not os.path.exists(CONFIG_PATH):
        raise FileNotFoundError(f"Dataset registry not found: {CONFIG_PATH}")
    with open(CONFIG_PATH, "r") as f:
        return yaml.safe_load(f)


def get_dataset_cfg(dataset: str) -> dict:
    """Return the registry entry for one dataset id."""
    cfg = load_dataset_config()["datasets"]
    if dataset not in cfg:
        raise KeyError(f"Unknown dataset '{dataset}'. Available: {sorted(cfg)}")
    return cfg[dataset]


def list_datasets() -> list:
    """
    Describe the registered datasets for the agent's planning step, as
    [{"name", "path", "description"}].

    The name is the dataset id the loader tools take, not a file path, since a
    single dataset holds many variables (and, for CMIP6, several models). Built
    from configs/datasets.yaml so the agent's catalogue cannot drift from what
    the openers actually read.
    """
    entries = []
    for name, cfg in load_dataset_config()["datasets"].items():
        parts = [f"[{cfg['kind']}]", " ".join(cfg.get("description", "").split())]

        if "models" in cfg:
            parts.append(f"Models: {', '.join(sorted(cfg['models']))}.")
            span = cfg.get("time", {}).get("range")
            if span:
                parts.append(f"Monthly, {span[0]} to {span[1]}.")
            parts.append("Variables differ by model; a variable a model does not publish "
                         "is reported as missing when it is loaded.")
        elif "variables" in cfg:
            parts.append(f"Variables: {', '.join(sorted(available_variables(cfg)))}.")
        elif "products" in cfg:
            parts.append(f"Products: {', '.join(sorted(cfg['products']))}.")

        if cfg.get("climatology_years"):
            parts.append(f"Climatology period {cfg['climatology_years']}.")

        entries.append({
            "name": name,
            "path": cfg["root"],
            "description": " ".join(p for p in parts if p),
        })
    return entries


def dataset_notes() -> dict:
    """
    The registry's cross-dataset conventions (longitude, depth, grid, time, unit
    basis, climatology periods), as {topic: note}.

    These are the rules a caller has to respect when combining two products, so
    the agent injects them into its code-generation prompt rather than restating
    them in a docstring that could drift from configs/datasets.yaml.
    """
    return load_dataset_config().get("cross_dataset_notes", {}) or {}


def available_variables(ds_cfg: dict) -> list:
    """The registered variables a dataset can actually deliver, leaving out any
    marked unavailable (WOA23's apparent oxygen utilization was never downloaded)."""
    return [
        name for name, cfg in ds_cfg.get("variables", {}).items()
        if not isinstance(cfg, dict) or cfg.get("available", True)
    ]


def resolve_variable(ds_cfg: dict, variable: str) -> str:
    """Map a user-supplied variable name onto a registered one via the alias table."""
    if variable is None:
        raise ValueError(
            f"A variable is required for this dataset. Choose one of "
            f"{sorted(available_variables(ds_cfg))}."
        )
    if variable in ds_cfg.get("variables", {}):
        return variable
    aliases = ds_cfg.get("aliases", {}) or {}
    key = variable.lower()
    if key in aliases:
        return aliases[key]
    for name in ds_cfg.get("variables", {}):
        if name.lower() == key:
            return name
    raise KeyError(
        f"Variable '{variable}' not available. "
        f"Known variables: {sorted(ds_cfg.get('variables', {}))}; "
        f"aliases: {sorted(aliases)}"
    )


# ---------------------------------------------------------------------------
# Normalisation primitives, called only by the openers below
# ---------------------------------------------------------------------------

def to_lon180(obj, lon_name: str = "lon", sort: bool = True):
    """
    Convert a longitude coordinate to the -180~180 convention.

    Handles every convention present in the registry (0~360, 20.5~379.5) with the
    same expression. 1-D coordinates are re-sorted ascending; 2-D curvilinear
    coordinates keep their native layout (they cannot be sorted).
    """
    lon = obj[lon_name]
    converted = ((lon + 180) % 360) - 180
    converted.attrs = dict(lon.attrs)
    converted.attrs["units"] = "degrees_east"
    obj = obj.assign_coords({lon_name: converted})
    if sort and lon.ndim == 1:
        obj = obj.sortby(lon_name)
    return obj


def snap_coords(obj, names=("lat", "lon", "depth"), decimals=6):
    """
    Round 1-D coordinates to a fixed precision, as float64.

    Grids that are nominally identical do not compare equal cell by cell:
    WOA23 stores its coordinates as float32, GFDL-ESM4's float64 latitudes
    carry ~1e-16 of representation error, and the modulo in to_lon180 adds
    ~5e-14 of its own. Any of that is enough for xarray's alignment to drop
    the offending cells silently when two products are combined, so every
    opener snaps its coordinates onto the same values.
    """
    updates = {}
    for name in names:
        if name in obj.coords and obj[name].ndim == 1 and obj[name].dtype.kind == "f":
            snapped = obj[name].astype("float64").round(decimals)
            snapped.attrs = dict(obj[name].attrs)
            updates[name] = snapped
    return obj.assign_coords(updates) if updates else obj


def sort_ascending(obj, names):
    """Sort 1-D coordinates ascending, skipping ones that are absent or 2-D."""
    for name in names:
        if name in obj.coords and obj[name].ndim == 1 and obj[name].size > 1:
            values = obj[name].values
            if values[0] > values[-1]:
                obj = obj.sortby(name)
    return obj


def normalize_depth(obj, depth_coord: str = "depth", scale: float = 1.0):
    """
    Rename the vertical coordinate to 'depth', convert it to metres and sort it
    downward. `scale` is 0.01 for CESM2 (native centimetres) and 1.0 elsewhere.
    """
    if depth_coord not in obj.coords and depth_coord not in obj.dims:
        return obj
    if depth_coord != "depth":
        obj = obj.rename({depth_coord: "depth"})
    if scale != 1.0:
        obj = obj.assign_coords(depth=obj["depth"] * scale)
    obj["depth"].attrs.update(
        {"units": "m", "positive": "down", "axis": "Z", "standard_name": "depth"}
    )
    return sort_ascending(obj, ["depth"])


def normalize_monthly_time(obj, time_name: str = "time"):
    """
    Put a monthly time axis on a common footing: plain datetime64 stamped on the
    first day of each month, whatever calendar (noleap, 365_day,
    proleptic_gregorian) or day-of-month convention the file used.

    Monthly means carry no information in the day-of-month, so this is lossless,
    and it makes CMIP6, ERSSTv5 and the CPC indices align exactly.
    """
    if time_name not in obj.coords:
        return obj
    values = obj[time_name].values
    if np.issubdtype(obj[time_name].dtype, np.datetime64):
        index = pd.DatetimeIndex(values)
        years, months = index.year, index.month
    else:  # cftime objects
        years = [t.year for t in values]
        months = [t.month for t in values]
    stamps = pd.to_datetime([f"{y:04d}-{m:02d}-01" for y, m in zip(years, months)])
    obj = obj.assign_coords({time_name: stamps})
    obj[time_name].attrs["stamp_convention"] = "month_start"
    return obj


def fix_unit_attribute(obj):
    """
    SEANOE writes the non-CF attribute `unit` instead of `units` on every
    coordinate and variable. Copy it across so CF-aware code sees a unit.
    """
    for name in list(obj.coords) + list(getattr(obj, "data_vars", [])):
        attrs = obj[name].attrs
        if "unit" in attrs and "units" not in attrs:
            attrs["units"] = attrs["unit"]
    return obj


# The attributes that identify which grid a field lives on. A derived quantity
# stays on the same grid as the field it came from, so it has to keep these or
# the grid-aware tools downstream (area and volume weighting, layer thickness)
# can no longer tell what it is — an OMZ thickness map that cannot be
# area-averaged is of little use.
PROVENANCE_ATTRS = (
    "source_dataset", "model", "experiment", "variant_label", "grid_type",
    "y_dim", "x_dim", "lat_coord", "lon_coord",
    "cmip6_variable", "grid_point", "vector_component",
    "woa_variable", "woa_time_slice", "woa_period",
    "region", "region_geofeature",
)

# The horizontal vector components the CMIP6 archive holds, and which grid axis
# each runs along. Their CMOR standard names are sea_water_x/y_velocity and
# ocean_mass_x/y_transport: x-ward and y-ward along the model's own grid lines,
# which on the displaced-pole (CESM2) and bipolar (MPI-ESM1-2-LR) grids are not
# east and north. They also sit on the velocity points of a staggered grid, not
# on the tracer points areacello and the tracer masks describe; checked on the
# local files, CESM2's u/v points are offset by a median 0.22 deg in latitude and
# 0.56 deg in longitude, and MPI's C-grid puts uo/umo on x-faces and vo/vmo on
# y-faces. wo and wmo, the vertical components, sit on the tracer points.
HORIZONTAL_VECTOR_COMPONENTS = {"uo": "x", "umo": "x", "vo": "y", "vmo": "y"}


def carry_provenance(result, source):
    """Copy the grid-identity attributes from `source` onto a derived field."""
    for key in PROVENANCE_ATTRS:
        if key in source.attrs and key not in result.attrs:
            result.attrs[key] = source.attrs[key]
    return result


def tag_provenance(da, dataset: str, ds_cfg: dict, **extra):
    """Record where a DataArray came from, so downstream code and the agent can
    see the climatology period, the molar basis and any product caveats."""
    da.attrs["source_dataset"] = dataset
    if ds_cfg.get("climatology_years"):
        da.attrs["climatology_period"] = ds_cfg["climatology_years"]
    for key, value in extra.items():
        if value is not None:
            da.attrs[key] = value
    notes = ds_cfg.get("notes")
    if notes:
        da.attrs["dataset_notes"] = " | ".join(notes)
    return da


# ---------------------------------------------------------------------------
# 1. WOA23
# ---------------------------------------------------------------------------

def _woa_period(var_cfg: dict, ds_cfg: dict, period, time_slice: str) -> str:
    """
    Pick the climatology period: the one requested, or else the first period in
    the registry's priority order that provides the requested time slice.

    That order prefers the normal with all three time slices over the newest
    one, so that a variable's annual, monthly and seasonal fields share a
    baseline and their anomalies stay comparable.
    """
    periods = var_cfg["periods"]
    if period is not None:
        if period not in periods:
            raise ValueError(
                f"Period '{period}' not available for this variable. "
                f"Available: {sorted(periods)}"
            )
        if time_slice not in periods[period]["depth_levels"]:
            raise ValueError(
                f"Period '{period}' has no '{time_slice}' files. "
                f"Available slices: {sorted(periods[period]['depth_levels'])}"
            )
        return period
    for candidate in ds_cfg["period_priority"]:
        if candidate in periods and time_slice in periods[candidate]["depth_levels"]:
            return candidate
    raise ValueError(
        f"No period provides time_slice '{time_slice}'. "
        f"Available: { {p: sorted(c['depth_levels']) for p, c in periods.items()} }"
    )


def woa23_files(variable: str, time_slice: str = "annual", period: str = None) -> dict:
    """
    Resolve which WOA23 files hold one variable's time slice.

    Returns the file list as [(nn, path)] along with the resolved variable name,
    the period actually chosen and its config. Shared by open_woa23, which reads
    the fields, and load_depth_levels, which needs the depth bounds those files
    carry (the bounds differ between time slices, so the file matters).
    """
    ds_cfg = get_dataset_cfg("woa23")
    variable = resolve_variable(ds_cfg, variable)
    var_cfg = ds_cfg["variables"][variable]

    if var_cfg.get("available", True) is False:
        raise FileNotFoundError(f"WOA23 variable '{variable}' is not available: {var_cfg.get('note', '')}")
    if time_slice not in ds_cfg["time"]["slices"]:
        raise ValueError(f"time_slice must be one of {sorted(ds_cfg['time']['slices'])}, got '{time_slice}'")

    period = _woa_period(var_cfg, ds_cfg, period, time_slice)
    period_cfg = var_cfg["periods"][period]
    prefix = var_cfg["prefix"]

    files = []
    for nn in ds_cfg["time"]["slices"][time_slice]:
        path = os.path.join(
            ds_cfg["root"], period_cfg["subdir"],
            ds_cfg["filename_template"].format(tag=period_cfg["tag"], prefix=prefix, nn=nn),
        )
        if not os.path.exists(path):
            raise FileNotFoundError(f"WOA23 file not found: {path}")
        files.append((nn, path))

    return {"files": files, "variable": variable, "period": period,
            "period_cfg": period_cfg, "prefix": prefix, "ds_cfg": ds_cfg}


def open_woa23(variable: str, time_slice: str = "annual", period: str = None,
               field: str = "an") -> xr.DataArray:
    """Open a World Ocean Atlas objectively analysed climatology field.

    A low-level reader returning the whole global field; use load_climatology_data to read a
    subset, which is the normal entry point.

    Args:
        variable: str, required. One of the physical or chemical variables the product
            publishes; common aliases are accepted.
        time_slice: str, optional, default 'annual'. 'annual' for a single all-year field
            reaching the deepest levels, 'monthly' for twelve fields with a 'month' dimension,
            'seasonal' for four with a 'season' dimension. The monthly and seasonal fields stop well
            short of the annual one's deepest level.
        period: str, optional. Which climate normal to read, where the variable offers more than
            one; otherwise a single all-years field is returned.
        field: str, optional, default 'an'. Which statistical field: 'an' the objectively
            analysed value, 'mn' the raw mean of observations, 'sd' their standard deviation, 'dd'
            the number of observations, and others the product defines.

    Returns:
        xr.DataArray with dims (depth, lat, lon), plus 'month' or 'season' for those slices, on
        a 1-degree grid with longitude in the -180 to 180 convention and depth in metres
        positive down. attrs carry the units, the climatology period and the provenance the
        weighting tools need.
    """
    resolved = woa23_files(variable, time_slice, period)
    ds_cfg, variable = resolved["ds_cfg"], resolved["variable"]
    period, period_cfg, prefix = resolved["period"], resolved["period_cfg"], resolved["prefix"]

    if field not in ds_cfg["fields"]:
        raise ValueError(f"field must be one of {sorted(ds_cfg['fields'])}, got '{field}'")
    var_name = f"{prefix}_{field}"

    arrays = []
    for nn, path in resolved["files"]:
        ds = xr.open_dataset(path, decode_times=False)
        if var_name not in ds:
            available = sorted(v for v in ds.data_vars if v.startswith(f"{prefix}_"))
            raise KeyError(
                f"Field '{var_name}' not in {os.path.basename(path)}. Available: {available}. "
                f"(Note: '_ma' exists only in monthly/seasonal files, '_sdo' only for temperature/salinity "
                f"and annual oxygen.)"
            )
        arrays.append(ds[var_name].squeeze("time", drop=True))

    if time_slice == "annual":
        da = arrays[0]
    elif time_slice == "monthly":
        da = xr.concat(arrays, dim=pd.Index(range(1, 13), name="month"))
        da["month"].attrs.update({"long_name": "month of year", "units": "1"})
    else:
        da = xr.concat(arrays, dim=pd.Index(SEASON_ORDER, name="season"))
        da["season"].attrs["long_name"] = "climatological season"

    da = normalize_depth(da, "depth")
    da = sort_ascending(da, ["lat", "lon"])
    da = snap_coords(da)
    da.name = variable
    return tag_provenance(
        da, "woa23", ds_cfg,
        climatology_period=period_cfg.get("climatology_years", period),
        woa_field=field, woa_period=period, woa_time_slice=time_slice,
        # The vertical axis depends on the variable, and not only through the
        # number of levels: the deepest layer of a monthly or seasonal axis is
        # closed at the level itself, so nitrate's 800 m layer is 25 m thick
        # where temperature's is 50 m. Recording the variable here keeps that
        # recoverable after a calculation has renamed the field.
        woa_variable=variable,
        original_units=da.attrs.get("units"),
    )


# ---------------------------------------------------------------------------
# 2. GLODAPv2.2016b
# ---------------------------------------------------------------------------

def open_glodap(variable: str, field: str = "main") -> xr.DataArray:
    """Open a GLODAP mapped climatology field of ocean carbon chemistry and hydrography.

    A low-level reader returning the whole global field; use load_climatology_data to read a
    subset.

    Args:
        variable: str, required. One of the physical, nutrient or carbon variables the product
            publishes, including dissolved inorganic carbon, alkalinity, pH and the saturation
            states.
        field: str, optional, default 'main'. Which field to read: the mapped value, or one of
            the accompanying error and input-count fields.

    Returns:
        xr.DataArray with dims (depth, lat, lon) and no time dimension, on a 1-degree grid with
        longitude rolled to the -180 to 180 convention. Concentrations are on a mass basis, so
        they must be converted before being compared with a volume-basis product.
    """
    ds_cfg = get_dataset_cfg("glodap")
    variable = resolve_variable(ds_cfg, variable)
    var_cfg = ds_cfg["variables"][variable]

    if field not in ds_cfg["fields"]:
        raise ValueError(f"field must be one of {sorted(ds_cfg['fields'])}, got '{field}'")

    path = os.path.join(ds_cfg["root"], ds_cfg["filename_template"].format(varfile=var_cfg["varfile"]))
    if not os.path.exists(path):
        raise FileNotFoundError(f"GLODAP file not found: {path}")

    ds = xr.open_dataset(path)
    depth_dim = ds_cfg["depth"]["dim"]
    ds = ds.assign_coords({depth_dim: ds[ds_cfg["depth"]["source_variable"]].values})
    ds = ds.rename({depth_dim: "depth"})
    ds = to_lon180(ds, "lon")
    ds = sort_ascending(ds, ["lat"])
    ds = normalize_depth(ds, "depth")
    ds = snap_coords(ds)

    var_name = ds_cfg["fields"][field].format(var=variable)
    if var_name not in ds:
        raise KeyError(f"Field '{var_name}' not in {os.path.basename(path)}. Available: {list(ds.data_vars)}")
    da = ds[var_name]

    # Input_mean / Input_std carry a hard-coded units attribute that is wrong for
    # every variable; the main variable's units are the authoritative ones.
    if field in ("input_mean", "input_std"):
        da.attrs["units"] = ds[variable].attrs.get("units", "")
        da.attrs["units_note"] = "corrected from the main variable; the file's own units attribute is wrong"

    da.name = f"{variable}_{field}" if field != "main" else variable
    return tag_provenance(
        da, "glodap", ds_cfg,
        glodap_field=field, original_units=da.attrs.get("units"),
    )


# ---------------------------------------------------------------------------
# 3. SEANOE mixed layer depth
# ---------------------------------------------------------------------------

def open_seanoe_mld(version: str = None, variable: str = "mld_dr003") -> xr.DataArray:
    """Open a global monthly climatology of ocean mixed layer depth.

    A low-level reader returning the whole global field; use load_climatology_data to read a
    subset. This is the observed MLD, built profile by profile, and is the reference to compare
    a computed MLD against.

    Args:
        version: str, optional. Which version to read; versions differ in the reference depth
            their criterion uses.
        variable: str, optional, default the density-threshold MLD. Which field to read: the
            mixed layer depth itself, or the accompanying mask.

    Returns:
        xr.DataArray with dims (month, lat, lon), month running 1-12, on a 1-degree grid with
        depth in metres. attrs record the criterion the product used.
    """
    ds_cfg = get_dataset_cfg("seanoe_mld")
    version = version or ds_cfg["default_version"]
    if version not in ds_cfg["versions"]:
        raise ValueError(f"version must be one of {sorted(ds_cfg['versions'])}, got '{version}'")
    version_cfg = ds_cfg["versions"][version]

    path = os.path.join(ds_cfg["root"], version_cfg["file"])
    if not os.path.exists(path):
        raise FileNotFoundError(f"SEANOE file not found: {path}")

    # The time axis says 'days since 1-1-1' under the non-CF attribute name
    # `unit`, so CF decoding neither fires nor would be meaningful.
    ds = xr.open_dataset(path, decode_times=False)
    ds = fix_unit_attribute(ds)
    ds = ds.assign_coords(month=("time", np.arange(1, 13))).swap_dims({"time": "month"})
    ds = ds.drop_vars("time")
    ds["month"].attrs.update({"long_name": "month of year", "units": "1"})
    ds = to_lon180(ds, "lon")
    ds = sort_ascending(ds, ["lat"])
    ds = snap_coords(ds)

    if variable not in ds:
        raise KeyError(f"Variable '{variable}' not in {os.path.basename(path)}. Available: {list(ds.data_vars)}")
    da = ds[variable]

    return tag_provenance(
        da, "seanoe_mld", ds_cfg,
        mld_version=version,
        reference_depth_m=version_cfg["reference_depth_m"],
        doi=version_cfg["doi"],
        original_units=da.attrs.get("units"),
    )


# ---------------------------------------------------------------------------
# 4. ERSSTv5
# ---------------------------------------------------------------------------

def open_ersst(variable: str = "sst") -> xr.DataArray:
    """Open the extended reconstructed monthly sea surface temperature record.

    A low-level reader returning the whole global field; use load_regular_data to read a subset.
    This is the long observed SST record, and the usual observational reference for SST
    variability and for climate indices.

    Args:
        variable: str, optional, default 'sst'. The dataset provides only sea surface
            temperature.

    Returns:
        xr.DataArray with dims (time, lat, lon) in degrees Celsius, monthly, stamped on the
        first of each month. It is on a coarser grid than the 1-degree climatologies, so it must
        be regridded before any cell-by-cell comparison with them.
    """
    ds_cfg = get_dataset_cfg("ersst")
    path = os.path.join(ds_cfg["root"], ds_cfg["file"])
    if not os.path.exists(path):
        raise FileNotFoundError(f"ERSSTv5 file not found: {path}")

    ds = xr.open_dataset(path)
    ds = to_lon180(ds, "lon")
    ds = sort_ascending(ds, ["lat"])          # native latitude runs 88 -> -88
    ds = snap_coords(ds)
    ds = normalize_monthly_time(ds)

    if variable not in ds:
        raise KeyError(f"Variable '{variable}' not in {os.path.basename(path)}. Available: {list(ds.data_vars)}")
    da = ds[variable]

    return tag_provenance(da, "ersst", ds_cfg, original_units=da.attrs.get("units"))


# ---------------------------------------------------------------------------
# 5. CPC ENSO indices
# ---------------------------------------------------------------------------

def _read_cpc_table(path: str, layout: str) -> pd.DataFrame:
    """Parse one whitespace-delimited CPC table, giving the four repeated ANOM
    columns names that say which Nino region they belong to."""
    if layout == "year_rows":
        # One row per year and one column per month, behind a title line. The
        # PDO index comes this way; 99.99 is its missing value, and the trailing
        # months of the current year carry it.
        table = pd.read_csv(path, sep=r"\s+", skiprows=1)
        months = {name: number for number, name in enumerate(
            ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
             "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"], start=1)}
        long = table.melt(id_vars="Year", value_vars=list(months),
                          var_name="MONTH_NAME", value_name="VALUE")
        long["MON"] = long["MONTH_NAME"].map(months)
        long = long.rename(columns={"Year": "YR"})
        long = long[long["VALUE"] < 99.0]        # drop the 99.99 sentinels
        return long.sort_values(["YR", "MON"]).reset_index(drop=True)

    df = pd.read_csv(path, sep=r"\s+")
    if layout == "seasonal":
        return df
    # Monthly files: YR MON NINO1+2 ANOM NINO3 ANOM NINO4 ANOM NINO3.4 ANOM,
    # which pandas de-duplicates into ANOM, ANOM.1, ANOM.2, ANOM.3.
    renamed, last_region = {}, None
    for column in df.columns:
        if column.startswith("NINO"):
            last_region = column
        elif column.startswith("ANOM") and last_region is not None:
            renamed[column] = f"{last_region}_ANOM"
    return df.rename(columns=renamed)


def open_cpc_index(product: str = None, region: str = None, field: str = "anomaly") -> xr.DataArray:
    """Open an official ENSO monitoring index as a monthly series.

    These are published area-averaged SST values and anomalies over the Nino regions, with no
    spatial dimensions: use them as the observed reference signal to correlate a field or a
    reconstruction against. Products differ in baseline period and are NOT interchangeable, and
    an index built on another baseline is offset from them by a constant.

    Args:
        product: str, optional. Which published series to read; they differ in the baseline
            period their anomalies are taken against, and one of them is a running-mean index rather
            than a monthly one. Default is the current operational series.
        region: str, optional. Which Nino region: 'NINO1+2', 'NINO3', 'NINO3.4' or 'NINO4';
            default 'NINO3.4' for the products that publish all four. Leave it unset for a
            product defined on a single region or on its own basin-wide domain, which reject
            a region they cannot provide rather than quietly returning another one.
        field: str, optional, default 'anomaly'. 'anomaly' for the departure from the baseline
            climatology, or 'total' for the absolute area-averaged temperature.

    Returns:
        xr.DataArray with a single 'time' dimension, monthly and stamped on the first of each
        month, in degrees Celsius. A running-mean product has its overlapping seasons placed at
        their centre month so it aligns with the monthly series. attrs record the baseline
        period and the region's box.
    """
    ds_cfg = get_dataset_cfg("cpc_indices")
    product = product or ds_cfg["default_product"]
    if product not in ds_cfg["products"]:
        raise ValueError(f"product must be one of {sorted(ds_cfg['products'])}, got '{product}'")
    if field not in ("anomaly", "total"):
        raise ValueError(f"field must be 'anomaly' or 'total', got '{field}'")
    product_cfg = ds_cfg["products"][product]

    path = os.path.join(ds_cfg["root"], product_cfg["file"])
    if not os.path.exists(path):
        raise FileNotFoundError(f"CPC file not found: {path}")

    df = _read_cpc_table(path, product_cfg["layout"])

    if product_cfg["layout"] == "year_rows":
        # A basin-wide EOF index, not a regional box average: there is no
        # absolute value to report and no region to choose.
        if field == "total":
            raise ValueError(
                f"product='{product}' publishes an index, not an area-averaged temperature, so it "
                f"has no 'total' field. Use field='anomaly'."
            )
        if region is not None:
            raise ValueError(
                f"product='{product}' is a basin-wide pattern index defined over its own domain, "
                f"so it cannot provide the '{region}' box. Drop the region argument, or use one of "
                f"the Nino products for a regional index."
            )
        time = pd.to_datetime(
            [f"{y:04d}-{m:02d}-01" for y, m in zip(df["YR"].astype(int), df["MON"].astype(int))]
        )
        values = df["VALUE"]
        region_used = product_cfg.get("domain", product)

    elif product_cfg["layout"] == "seasonal":
        months = df["SEAS"].map(ONI_SEASON_TO_CENTER_MONTH)
        if months.isna().any():
            raise ValueError(f"Unrecognised season labels in {path}: {sorted(set(df['SEAS'][months.isna()]))}")
        time = pd.to_datetime(
            [f"{y:04d}-{m:02d}-01" for y, m in zip(df["YR"].astype(int), months.astype(int))]
        )
        if region not in (None, "NINO3.4"):
            raise ValueError(
                f"product='{product}' is the ONI, which is defined on Nino3.4 alone; it cannot provide "
                f"'{region}'. Use one of the monthly products "
                f"({[p for p, c in ds_cfg['products'].items() if c['layout'] == 'monthly']}) for other regions."
            )
        values = df["ANOM"] if field == "anomaly" else df["TOTAL"]
        region_used = "NINO3.4"
    else:
        region = region or "NINO3.4"
        time = pd.to_datetime(
            [f"{y:04d}-{m:02d}-01" for y, m in zip(df["YR"].astype(int), df["MON"].astype(int))]
        )
        column = f"{region}_ANOM" if field == "anomaly" else region
        if column not in df.columns:
            raise KeyError(f"Region '{region}' not in {os.path.basename(path)}. Columns: {list(df.columns)}")
        values = df[column]
        region_used = region

    da = xr.DataArray(
        values.to_numpy(dtype="float64"),
        coords={"time": ("time", time, {"stamp_convention": "month_start"})},
        dims=["time"],
        name=f"{region_used}_{'anomaly' if field == 'anomaly' else 'sst'}",
        # The Nino products are temperatures; the PDO is a dimensionless pattern
        # index, so the unit comes from the product entry where it states one.
        attrs={"units": product_cfg.get("units", "degC"),
               "long_name": f"{region_used} {field} from {product_cfg.get('source', 'CPC')} {product}"},
    )
    return tag_provenance(
        da, "cpc_indices", ds_cfg,
        cpc_product=product, region=region_used, field=field,
        baseline=product_cfg.get("baseline"),
        region_box=str(ds_cfg["regions"].get(region_used)),
        product_note=product_cfg.get("note"),
    )


# ---------------------------------------------------------------------------
# 6. CMIP6 benchmark (historical / ssp245)
# ---------------------------------------------------------------------------

def available_cmip6_variables(experiment: str, model: str) -> list:
    """List the variable directories actually present for one model."""
    ds_cfg = get_dataset_cfg(f"cmip6_{experiment}")
    model_root = os.path.join(ds_cfg["root"], model)
    if not os.path.isdir(model_root):
        return []
    return sorted(d for d in os.listdir(model_root) if os.path.isdir(os.path.join(model_root, d)))


def _cmip6_files(ds_cfg: dict, experiment: str, model: str, variable: str):
    """
    Locate the files for one variable, returning (files, variant, experiment_used).

    Variant labels are read off the filenames rather than assumed: they are not
    uniform across this archive (historical/CESM2/phycos is r6i1p1f1,
    historical/MPI-ESM1-2-LR/ph is r10i1p1f1, most of ssp245/CESM2 is r10i1p1f1
    while its fixed fields are r1i1p1f1).
    """
    var_dir = os.path.join(ds_cfg["root"], model, variable)
    used_experiment = experiment

    if not os.path.isdir(var_dir) and variable in ds_cfg.get("fixed_fields", []):
        # ssp245/GFDL-ESM4 ships no Ofx fields; the grid is identical to
        # historical, so borrow them from there.
        fallback = ds_cfg.get("fixed_field_fallback")
        if fallback:
            fallback_cfg = get_dataset_cfg(fallback)
            var_dir = os.path.join(fallback_cfg["root"], model, variable)
            used_experiment = fallback_cfg["experiment"]

    if not os.path.isdir(var_dir):
        raise FileNotFoundError(
            f"Variable '{variable}' not available for {model} in {experiment}. "
            f"Available: {available_cmip6_variables(experiment, model)}"
        )

    files = sorted(glob.glob(os.path.join(var_dir, "*.nc")))
    if not files:
        raise FileNotFoundError(f"No NetCDF files in {var_dir}")

    variants = set()
    for path in files:
        match = re.match(rf"{re.escape(variable)}_\w+_{re.escape(model)}_\w+_(r\d+i\d+p\d+f\d+)_", os.path.basename(path))
        if match:
            variants.add(match.group(1))
    if len(variants) > 1:
        raise ValueError(f"Mixed ensemble members in {var_dir}: {sorted(variants)}")

    return files, (variants.pop() if variants else "unknown"), used_experiment


# CESM2 labels `lev` as centimetres but its `lev_bnds` as metres, in the same
# file, so level centres and their bounds must each be converted according to
# the units they actually declare. Taking the bounds to share the level's units
# makes every CESM2 layer 100 times too thin.
DEPTH_UNIT_TO_METRES = {
    "m": 1.0, "meter": 1.0, "meters": 1.0, "metre": 1.0, "metres": 1.0,
    "cm": 0.01, "centimeter": 0.01, "centimeters": 0.01, "centimetre": 0.01, "centimetres": 0.01,
}


def depth_scale_from_units(units, default: float) -> float:
    """Metres per unit for a declared depth unit, falling back to `default`
    when the variable declares none."""
    if not units:
        return default
    key = str(units).strip().lower()
    if key not in DEPTH_UNIT_TO_METRES:
        raise ValueError(f"Unrecognised depth unit '{units}'; cannot convert to metres.")
    return DEPTH_UNIT_TO_METRES[key]


def cmip6_level_bounds_m(ds, depth_coord: str, default_scale: float):
    """Layer bounds in metres from an open CMIP6 file, or None if it has none."""
    if depth_coord not in ds.coords:
        return None
    bounds_name = ds[depth_coord].attrs.get("bounds", f"{depth_coord}_bnds")
    if bounds_name not in ds:
        return None
    level_scale = depth_scale_from_units(ds[depth_coord].attrs.get("units"), default_scale)
    bounds_scale = depth_scale_from_units(ds[bounds_name].attrs.get("units"), level_scale)
    return ds[bounds_name].values.astype("float64") * bounds_scale


def cmip6_ocean_mask(ds_cfg: dict, experiment: str, model: str, model_cfg: dict,
                     level_bounds_m=None):
    """
    The model's ocean mask, for the models whose fixed fields ship unmasked.

    CESM2 publishes areacello, volcello and deptho over its whole global grid,
    land included, while GFDL-ESM4 and MPI-ESM1-2-LR mask theirs to ocean. This
    project is ocean-only, so the registry names a mask source for such models
    and the fixed fields are masked on the way out, making all three behave the
    same.

    Pass `level_bounds_m` for a 3-D field: the land/sea split alone is not
    enough there, because CESM2's volcello also fills every level below the sea
    floor — summed raw it comes to the ocean area times the full 5500 m column,
    about 148% of the real ocean volume. Cells below deptho are dropped as well,
    which reproduces the model's own 3-D ocean extent exactly (it agrees with
    thetao's valid cells to the last cell) and brings the volume in line with
    the other two models.

    Returns (mask_values, source_name), or (None, None) when no masking applies.
    """
    source = model_cfg.get("ocean_mask_from")
    if not source:
        return None, None

    files, _, _ = _cmip6_files(ds_cfg, experiment, model, source)
    with xr.open_dataset(files[0]) as ds:
        # CESM2's sftof is strictly 0 or 100, so this is a clean land/sea split
        # rather than a fractional one.
        land_sea = ds[source].values > 0
    if level_bounds_m is None:
        return land_sea, source

    floor_files, _, _ = _cmip6_files(ds_cfg, experiment, model, "deptho")
    with xr.open_dataset(floor_files[0]) as ds:
        sea_floor = ds["deptho"].values.astype("float64")
    sea_floor = np.where(np.isnan(sea_floor), -1.0, sea_floor)

    above_floor = level_bounds_m[:, 0][:, None, None] < sea_floor[None, :, :]
    return (above_floor & land_sea[None, :, :]), f"{source} and deptho"


def open_cmip6_grid_file(experiment: str, model: str, variable: str = "thetao"):
    """
    Open one chunk file of a CMIP6 variable for its grid metadata alone.

    The cell corner bounds and level bounds are data variables, so selecting the
    DataArray in open_cmip6 drops them; load_areas and load_depth_levels need
    them and open the file through here instead. Returns (dataset, info) with
    the raw, un-normalised coordinates — the caller decides what to convert.
    """
    ds_cfg = get_dataset_cfg(f"cmip6_{experiment}")
    if model not in ds_cfg["models"]:
        raise ValueError(f"model must be one of {sorted(ds_cfg['models'])}, got '{model}'")

    files, variant, used_experiment = _cmip6_files(ds_cfg, experiment, model, variable)
    info = {
        "model": model,
        "experiment": used_experiment,
        "variant_label": variant,
        "model_cfg": ds_cfg["models"][model],
        "source_file": files[0],
    }
    return xr.open_dataset(files[0], use_cftime=True), info


def cmip6_experiments() -> list:
    """The registered CMIP6 experiment names, for cross-experiment fallbacks."""
    return [
        name.replace("cmip6_", "")
        for name in load_dataset_config()["datasets"]
        if name.startswith("cmip6_")
    ]


def open_cmip6(experiment: str, model: str, variable: str, chunk_time: int = 12) -> xr.DataArray:
    """Open a CMIP6 ocean variable from the model archive as a monthly time series.

    A low-level reader returning the whole global field; use load_regular_data or
    load_data_with_geolocator to read a subset, which is the normal entry point and is the only
    way to subset a curvilinear model grid by latitude and longitude.
    Multi-file variables are concatenated along time, the models' different calendars are
    unified, and depth is normalised to metres. Variable coverage differs by model and
    experiment; the error message lists what is actually present.

    Args:
        experiment: str, required. Which experiment to read.
        model: str, required. Which model to read.
        variable: str, required. The CMIP6 short name — a 3-D tracer, a 2-D surface or
            depth-integrated field, or one of the fixed grid fields (cell area, cell volume, sea
            floor depth, sea area fraction).
        chunk_time: int, optional, default 12. Chunk size along time. The array stays lazy, so
            opening even a very large variable is cheap; computation happens when it is reduced.

    Returns:
        xr.DataArray, monthly and stamped on the first of each month, with depth in metres
        positive down and longitude in the -180 to 180 convention. On a regular model grid the
        dimensions are (time, depth, lat, lon); on a curvilinear one they are the model's index
        dimensions with 2-D latitude and longitude coordinates, and those index dimensions are
        NOT geographic and must not be sliced by value. Concentrations are on a volume basis, so
        they must be converted before being compared with a mass-basis product. attrs carry the
        grid identity that lets the weighting tools work without any manual use of cell areas.
    """
    dataset = f"cmip6_{experiment}"
    ds_cfg = get_dataset_cfg(dataset)
    if model not in ds_cfg["models"]:
        raise ValueError(f"model must be one of {sorted(ds_cfg['models'])}, got '{model}'")
    if variable is None:
        raise ValueError(
            f"A variable is required. {model} in {experiment} provides: "
            f"{available_cmip6_variables(experiment, model)}"
        )
    model_cfg = ds_cfg["models"][model]

    files, variant, used_experiment = _cmip6_files(ds_cfg, experiment, model, variable)

    open_kwargs = {"use_cftime": True}
    try:
        if len(files) == 1:
            ds = xr.open_dataset(files[0], chunks={"time": chunk_time}, **open_kwargs)
        else:
            ds = xr.open_mfdataset(
                files, combine="by_coords", chunks={"time": chunk_time},
                parallel=False, data_vars="minimal", coords="minimal", compat="override",
                **open_kwargs,
            )
    except OSError as exc:
        # A truncated or corrupt chunk file surfaces as a bare HDF error that
        # names no file, so say which one of the chunks cannot be read.
        unreadable = []
        for path in files:
            try:
                xr.open_dataset(path, **open_kwargs).close()
            except OSError:
                unreadable.append(f"{os.path.basename(path)} ({os.path.getsize(path) / 1e6:.1f} MB)")
        if unreadable:
            raise OSError(
                f"Cannot read {variable} for {model}/{experiment}: "
                f"{len(unreadable)} of {len(files)} chunk files are corrupt and need re-downloading: "
                f"{unreadable}"
            ) from exc
        raise

    if variable not in ds:
        raise KeyError(f"Variable '{variable}' not in {os.path.basename(files[0])}. Available: {list(ds.data_vars)}")
    da = ds[variable]

    # Fixed grid fields are masked to ocean for the models that ship them
    # unmasked, so that areas, volumes and depths mean the same thing across
    # models. Data variables are already NaN over land and are left alone.
    depth_coord = model_cfg["depth_coord"]
    ocean_mask, mask_source = (None, None)
    if variable in ds_cfg.get("fixed_fields", []):
        # The three models disagree on how a fixed field marks land: CESM2
        # leaves areacello and volcello unmasked and writes 0 in deptho,
        # GFDL-ESM4 writes 0 in areacello, MPI-ESM1-2-LR uses NaN throughout.
        # None of these fields can legitimately be zero (the smallest real
        # values are 6e4 m2, 3e2 m3 and 9.5 m), so a zero means no ocean.
        da = da.where(da != 0)
        # The mask must be looked up in the experiment the field came from,
        # which is not the requested one when a fixed field was borrowed.
        mask_cfg = get_dataset_cfg(f"cmip6_{used_experiment}")
        ocean_mask, mask_source = cmip6_ocean_mask(
            mask_cfg, used_experiment, model, model_cfg,
            cmip6_level_bounds_m(ds, depth_coord, model_cfg["depth_scale"]),
        )
    if ocean_mask is not None:
        dims = (model_cfg["y_dim"], model_cfg["x_dim"])
        if ocean_mask.ndim == 3:
            dims = (depth_coord,) + dims
        da = da.where(xr.DataArray(ocean_mask, dims=dims))

    if depth_coord in da.coords:
        da = normalize_depth(da, depth_coord, scale=model_cfg["depth_scale"])

    # Longitude goes to -180~180 for every model. On the curvilinear grids the
    # coordinate is 2-D, so it is converted in place and not re-sorted.
    lon_coord = model_cfg["lon_coord"]
    if lon_coord in da.coords:
        da = to_lon180(da, lon_coord, sort=(model_cfg["grid_type"] == "regular"))
    if model_cfg["grid_type"] == "regular":
        da = sort_ascending(da, [model_cfg["lat_coord"]])

    da = snap_coords(da)
    da = normalize_monthly_time(da)

    return tag_provenance(
        da, dataset, ds_cfg,
        model=model, experiment=used_experiment, variant_label=variant,
        grid_type=model_cfg["grid_type"],
        y_dim=model_cfg["y_dim"], x_dim=model_cfg["x_dim"],
        lat_coord=model_cfg["lat_coord"], lon_coord=model_cfg["lon_coord"],
        n_source_files=len(files),
        borrowed_from=used_experiment if used_experiment != experiment else None,
        original_units=da.attrs.get("units"),
        ocean_masked_with=mask_source,
        # Which file describes this field's points, carried onto everything
        # derived from it: a velocity-point field regridded against the tracer
        # corner mesh lands half a cell out of place.
        cmip6_variable=variable,
        grid_point="velocity" if variable in HORIZONTAL_VECTOR_COMPONENTS else "tracer",
        vector_component=HORIZONTAL_VECTOR_COMPONENTS.get(variable),
    )
