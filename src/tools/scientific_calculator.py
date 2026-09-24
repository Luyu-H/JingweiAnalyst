# tools: compute_anomaly_vs_woa23, compute_anomaly_vs_own_climatology,
# convert_concentration_to_content, compute_vertical_gradient,
# compute_ratio, identify_oxygen_minimum_zone, compute_potential_density_anomaly, compute_mixed_layer_depth,
# compute_nutricline_depth, compute_enso_nino_index

import xarray as xr
import numpy as np
import gsw
from sklearn.decomposition import PCA

from src.tools.registry import register_tools
from src.tools.data_openers import (SEASON_ORDER, carry_provenance, molar_basis,
                                    normalize_units, open_woa23)
from src.tools.data_loader import (areas_like, cell_thickness_like, layer_bounds_like,
                                   layer_thickness_like, volumes_like)
from src.tools.basic_calculator import compute_basic_stat


def _broadcastable(coord: xr.DataArray, data: xr.DataArray) -> np.ndarray:
    """
    Reshape a coordinate so it broadcasts against `data`'s array.

    gsw works on plain arrays, so depth, latitude and longitude have to be laid
    out on the axes they belong to before they are handed over. Passing the
    bare 1-D coordinates instead only works when the field is a single profile:
    on a (depth, lat, lon) field, depth of length 102 and latitude of length 180
    cannot broadcast against each other at all.
    """
    shape = [1] * data.ndim
    for dim in coord.dims:
        shape[data.dims.index(dim)] = coord.sizes[dim]
    return np.asarray(coord.values, dtype="float64").reshape(shape)


def temperature_kind(temp_data: xr.DataArray, declared: str = "auto") -> tuple:
    """
    Whether a temperature field is in-situ or potential temperature, and why.

    TEOS-10 needs to know: in-situ temperature is converted to Conservative
    Temperature with CT_from_t, potential temperature with CT_from_pt, and using
    the wrong one biases density. CMIP6's 'thetao' is potential temperature by
    definition, while WOA23 and GLODAP report in-situ temperature, so mixing
    them up is easy and the error grows with depth: about 0.01 kg/m3 in sigma0
    at 2000 m and 0.026 kg/m3 at 4000 m, the latter comparable to the 0.03
    kg/m3 threshold used to define a mixed layer.
    """
    if declared in ("in_situ", "potential"):
        return declared, "specified by the caller"
    if declared != "auto":
        raise ValueError(
            f"temperature_type must be 'auto', 'in_situ' or 'potential', got '{declared}'"
        )
    standard_name = str(temp_data.attrs.get("standard_name", "")).lower()
    if "potential_temperature" in standard_name:
        return "potential", f"standard_name '{standard_name}'"
    if str(temp_data.name).lower() in ("thetao", "theta", "pt", "potential_temperature"):
        return "potential", f"variable name '{temp_data.name}'"
    if "sea_water_temperature" in standard_name:
        return "in_situ", f"standard_name '{standard_name}'"
    return "in_situ", "assumed: nothing in the metadata identifies it as potential temperature"


def surface_layer(data: xr.DataArray, depth_level=None) -> xr.DataArray:
    """
    The sea-surface value of a field that has a depth dimension.

    By default the shallowest level alone, because that is how CMIP6 defines
    its surface fields: tos is "the prognostic temperature in the upper-most
    ocean grid cell" (Griffies et al. 2016, H13), and sos likewise. The top cell
    is 0-10 m in CESM2, 0-2 m natively in GFDL-ESM4 (2.5 m level as archived)
    and 0-12 m in MPI-ESM1-2-LR. A number instead averages every level
    shallower than it, thickness-weighted - which is a different quantity, and
    covers two levels in GFDL-ESM4 but one in the other models at 10 m.
    """
    if depth_level is None:
        ordered = data.sortby("depth")
        top = float(ordered["depth"].values[0])
        result = ordered.isel(depth=0, drop=True)
        result.attrs = dict(data.attrs)
        result.attrs["depth_averaging"] = f"shallowest level only ({top:g} m), the top cell as tos is defined"
        return result
    return depth_average(data, depth_level)


def depth_average(data: xr.DataArray, max_depth: float) -> xr.DataArray:
    """
    Average the layers above `max_depth`, weighted by layer thickness.

    Layers are not equally thick anywhere in this archive, so a plain mean over
    depth weights a 5 m layer the same as a 100 m one. A CMIP6 field is weighted
    by the water each cell actually holds (volcello / areacello), so a partial
    cell at the sea floor counts for what it contains rather than a full layer.
    Where the thicknesses cannot be established — after an anomaly has
    intersected two products' vertical grids, for instance — the average falls
    back to unweighted and says so in the attributes rather than failing.
    """
    upper = data.sel(depth=slice(0, max_depth))
    if upper.sizes.get("depth", 0) == 0:
        raise ValueError(
            f"No levels at or above {max_depth} m; this field spans "
            f"{float(data.depth.min()):g} to {float(data.depth.max()):g} m."
        )
    if upper.sizes["depth"] == 1:
        result = upper.isel(depth=0, drop=True)
        result.attrs = dict(data.attrs)
        result.attrs["depth_averaging"] = f"single level at {float(upper.depth[0]):g} m"
        return result

    try:
        thickness = cell_thickness_like(upper)
        weights = thickness.where(upper.notnull())
        result = (upper * weights).sum("depth", skipna=True) / weights.sum("depth")
        note = f"thickness-weighted ({thickness.attrs.get('bounds_source', 'layer thickness')})"
    except Exception as exc:
        result = upper.mean("depth", skipna=True)
        note = f"unweighted (layer thicknesses unavailable: {exc.__class__.__name__})"

    result.attrs = dict(data.attrs)
    result.attrs["depth_averaging"] = f"{note} over 0-{max_depth} m"
    return result


@register_tools(
    name="compute_anomaly_vs_woa23",
    description="Compute the anomaly of a gridded time series against the WOA23 observed climatology.\n"
        "Each time step is matched to the climatology of its own calendar month, season or to the annual "
        "mean, so the seasonal cycle is removed along with the mean state. Because the reference is "
        "observed, the result carries the field's bias against observations inside it: this is the tool "
        "for 'how far is this field from the observed ocean'.\n"
        "Its sibling compute_anomaly_vs_own_climatology removes the field's OWN climatology instead, "
        "leaving variability with the bias taken out — that is the one for an index, a trend or a "
        "comparison between two products, and the only one that works off this climatology's grid or on a "
        "chosen baseline period.\n"
        "The input must sit on the climatology's grid and share its units; a spatial subset is fine, but "
        "aggregation, interpolation or regridding beforehand is not.",
    args={
        "data": "xr.DataArray, required. A field on the climatology's grid, with a 'time' coordinate.",
        "var": "str, required. Which climatology to use: 'temp', 'so', 'si', 'o2', 'no3' or 'po4'. Full "
            "names such as 'temperature' or 'nitrate' are also accepted.",
        "freq": "str, optional, default 'month'. Which climatology each step is matched against: 'month' "
            "(12 fields), 'season' (4 fields) or 'year' (one annual field). The seasons here are the "
            "calendar quarters the seasonal files are built on, NOT the meteorological DJF/MAM/JJA/SON; "
            "the two are offset by a month.",
        "method": "str, optional, default 'subtract'. 'subtract' gives data minus climatology in the "
            "data's units; 'ratio' gives the difference as a percentage of the climatology, which is only "
            "accepted for concentrations, since a percentage of temperature, salinity or pH has no "
            "physical meaning.",
        "period": "str, optional. Which climate normal to use where the variable offers more than one. If "
            "None, a default period is chosen so that anomalies computed with different freq values share "
            "a baseline."
    },
    returns="xr.DataArray of anomalies, in the input's units for 'subtract' and in percent for 'ratio'. "
        "The monthly and seasonal climatologies are shallower than the annual one, so the result is "
        "trimmed to the depth levels the two share. attrs record which climatology period was used, which "
        "generally differs from the input's own period."
)
def compute_anomaly_vs_woa23(data: xr.DataArray, var: str, freq: str = "month",
                             method: str = "subtract", period: str = None) -> xr.DataArray:
    if freq not in ["month", "season", "year"]:
        raise ValueError(f"freq must be 'month', 'season', or 'year', got '{freq}'")
    if method not in ["subtract", "ratio"]:
        raise ValueError(f"method must be 'subtract' or 'ratio', got '{method}'")
    if data.attrs.get("grid_type") == "curvilinear":
        raise ValueError(
            f"{data.attrs.get('model', 'This dataset')} is on a curvilinear grid, whose index dimensions do "
            f"not align with the WOA23 1-degree grid. Put it on the WOA23 grid first with "
            f"regrid_horizontal, use GFDL-ESM4 or an observational product, or remove the field's own "
            f"climatology with compute_anomaly_vs_own_climatology."
        )
    if "time" not in data.coords:
        raise ValueError("data must have a 'time' coordinate to be matched against the climatology.")

    slice_for_freq = {"year": "annual", "month": "monthly", "season": "seasonal"}
    climatology = open_woa23(var, time_slice=slice_for_freq[freq], period=period)

    data_units = data.attrs.get("units")
    clim_units = climatology.attrs.get("units")
    if normalize_units(data_units) != normalize_units(clim_units):
        raise ValueError(
            f"Unit mismatch: the data is in '{data_units}' but the WOA23 climatology is in '{clim_units}'. "
            f"CMIP6 reports concentrations per unit volume (mol m-3) while WOA23 reports them per unit mass "
            f"(micromoles per kilogram); convert with in-situ density before differencing them "
            f"(mol m-3 = micromol kg-1 * rho / 1000)."
        )

    # A single level selected with isel leaves depth behind as a scalar
    # coordinate, which cannot be aligned against the climatology's depth axis.
    # Take the matching climatology level instead and drop the coordinate.
    if "depth" in data.coords and data["depth"].ndim == 0:
        if "depth" in climatology.dims:
            climatology = climatology.sel(depth=data["depth"], method="nearest")
        data = data.drop_vars("depth")
        climatology = climatology.drop_vars("depth", errors="ignore")

    # A surface field (CMIP6 'tos', ERSSTv5 'sst') has no depth dimension at all.
    # Left alone, it would broadcast against the climatology's 57 levels and the
    # "anomaly" would be the surface value minus the temperature at every depth —
    # a plausible-looking array that is not an anomaly. A surface field is
    # comparable to the shallowest climatology level, so take that one.
    surface_note = None
    if "depth" not in data.dims and "depth" in climatology.dims:
        surface_level = float(climatology["depth"].min())
        climatology = climatology.sel(depth=surface_level).drop_vars("depth", errors="ignore")
        surface_note = (f"input has no depth dimension, so it was compared against the shallowest "
                        f"WOA23 level ({surface_level:g} m)")

    # A dateline-crossing subset arrives on a continuous axis running past 180
    # (the loaders shift the western piece by +360 to keep the window contiguous
    # and the axis monotonic). Those labels are absent from WOA23's -180~180
    # grid, so match on the wrapped values and put the array's own labels back on
    # the result. The cells are the same cells either way.
    original_lon = None
    if "lon" in data.coords and float(np.max(np.atleast_1d(data["lon"].values))) > 180.0:
        original_lon = data["lon"]
        data = data.assign_coords(
            lon=((np.asarray(original_lon.values, dtype=float) + 180.0) % 360.0) - 180.0
        )

    # The climatology is on WOA23's 1-degree grid. A field on any other regular
    # grid — ERSSTv5's 2 degrees, or anything interpolated or regridded — shares
    # no cell centres with it, and the inner join below would then return an
    # empty array rather than raising.
    for axis in ("lat", "lon"):
        if axis not in data.coords or axis not in climatology.coords:
            continue
        data_labels = np.round(np.atleast_1d(data[axis].values).astype(float), 4)
        clim_labels = np.round(np.asarray(climatology[axis].values, dtype=float), 4)
        missing = data_labels[~np.isin(data_labels, clim_labels)]
        if missing.size:
            raise ValueError(
                f"{missing.size} of this array's {len(data_labels)} '{axis}' values have no cell on the "
                f"WOA23 1-degree grid (first few: {missing[:5].tolist()}). This tool matches the data "
                f"cell by cell against WOA23, so the field must be on that grid: ERSSTv5 is on a 2-degree "
                f"grid, and interpolating or regridding takes a field off it too. A spatial subset of the "
                f"WOA23 grid is fine. Either put the field on the WOA23 grid with regrid_horizontal, or "
                f"remove its own climatology with compute_anomaly_vs_own_climatology, which works on any "
                f"grid."
            )

    # Trim both to the levels and cells they share before any arithmetic. The
    # monthly and seasonal climatologies are shallower than the annual one, and
    # groupby arithmetic reuses the data's chunking for the climatology, which
    # fails outright if the two disagree on a dimension's length.
    data, climatology = xr.align(data, climatology, join="inner", exclude=("time", "month", "season"))
    if data.sizes.get("depth", 1) == 0:
        raise ValueError(
            f"The data and the WOA23 '{slice_for_freq[freq]}' climatology share no depth levels. "
            f"The monthly and seasonal climatologies only reach 800 m (nutrients) or 1500 m "
            f"(temperature/salinity/oxygen), and levels must match exactly."
        )

    if method == "ratio":
        # A percentage anomaly only means something on a ratio scale, i.e. one
        # with a true zero. Degrees Celsius is an interval scale: a percentage
        # of it diverges as the climatology approaches 0 C and changes sign
        # below it, so "10 percent warmer" is arithmetic without physical
        # content. Concentrations are fine.
        if molar_basis(clim_units) is None:
            raise ValueError(
                f"method='ratio' gives a percentage anomaly, which requires a quantity with a true zero. "
                f"'{var}' is in '{clim_units}', an interval scale (degrees Celsius, practical salinity, pH), "
                f"where a percentage of the climatology diverges near zero and flips sign below it. Use "
                f"method='subtract' for this variable; 'ratio' is meant for concentrations."
            )
        # Zero climatology values would blow up the ratio; NaN propagates instead.
        climatology = climatology.where(climatology != 0)

    if freq == "year":
        difference = data - climatology
        anomaly = difference if method == "subtract" else difference / climatology * 100
    else:
        if freq == "month":
            group_coord = "month"
            labels = data["time"].dt.month
        else:
            # WOA23 seasons are Jan-Mar, Apr-Jun, Jul-Sep, Oct-Dec, so the season
            # index is simply the quarter of the year, not the usual DJF grouping.
            group_coord = "season"
            labels = xr.DataArray(
                np.asarray(SEASON_ORDER)[(data["time"].dt.month.values - 1) // 3],
                dims="time", coords={"time": data["time"]},
            )
        data = data.assign_coords({group_coord: labels})
        difference = data.groupby(group_coord) - climatology
        if method == "subtract":
            anomaly = difference
        else:
            anomaly = difference.groupby(group_coord) / climatology * 100
        anomaly = anomaly.drop_vars(group_coord)

    anomaly.attrs = data.attrs.copy()
    anomaly.attrs["anomaly_method"] = method
    anomaly.attrs["anomaly_freq"] = freq
    if freq == "season":
        # WOA23's seasonal files are calendar quarters, so this tool is bound to
        # that convention. temporal_aggregate defaults to DJF/MAM/JJA/SON, which
        # is a different three months; pass season_convention='woa' there to
        # match this, and never composite one against the other.
        anomaly.attrs["season_convention"] = (
            "woa (JFM/AMJ/JAS/OND calendar quarters, as WOA23 defines them) - NOT the meteorological "
            "DJF/MAM/JJA/SON; use temporal_aggregate(season_convention='woa') to stay consistent"
        )
    anomaly.attrs["climatology_reference"] = (
        f"WOA23 {climatology.attrs.get('woa_period', '')} {climatology.attrs.get('woa_time_slice', '')}".strip()
    )
    if original_lon is not None:
        # Back to the caller's longitude convention, in the caller's order: the
        # inner join above may have reordered the axis while intersecting it.
        wrapped = ((np.asarray(original_lon.values, dtype=float) + 180.0) % 360.0) - 180.0
        anomaly = anomaly.sel(lon=wrapped).assign_coords(lon=original_lon)

    anomaly.attrs["climatology_period"] = climatology.attrs.get("climatology_period", "")
    if surface_note:
        anomaly.attrs["climatology_depth"] = surface_note

    long_name = data.attrs.get("long_name", var)
    if method == "subtract":
        anomaly.attrs["long_name"] = f"{long_name} anomaly"
    else:
        anomaly.attrs["long_name"] = f"{long_name} percentage anomaly"
        anomaly.attrs["units"] = "percent"

    return anomaly


@register_tools(
    name="compute_anomaly_vs_own_climatology",
    description="Remove a field's own climatology, leaving its variability about its own mean state.\n"
        "Each time step is matched to the mean of its own calendar month, season or of the whole record, "
        "computed from the same field over a baseline period you choose. Whatever bias the field carries "
        "cancels, which is what an index, a trend, a composite or a comparison between two products "
        "needs. Its sibling compute_anomaly_vs_woa23 measures the departure from observations instead, "
        "keeping the bias in.\n"
        "It works on any grid and any product, including a series that is already an area average with no "
        "spatial dimensions left.\n"
        "THE BASELINE IS PART OF THE ANSWER: an anomaly is only defined relative to the period its "
        "climatology came from, and two series built on different baselines differ by a constant offset "
        "even when their variability is identical. Match the baseline of any published series being "
        "compared against.",
    args={
        "data": "xr.DataArray, required. A field or series with a 'time' dimension.",
        "freq": "str, optional, default 'month'. What each step is matched against: 'month' (the "
            "climatology of its own calendar month, removing the seasonal cycle), 'season' (the "
            "meteorological DJF/MAM/JJA/SON), or 'year' (a single long-term mean, leaving the seasonal "
            "cycle in place).",
        "base_period": "list[str, str], optional. The baseline window as ['start', 'end']. The climatology "
            "is computed over this window only and subtracted from every step of the series, including "
            "those outside it. Default None: the whole record is the baseline.",
        "method": "str, optional, default 'subtract'. 'subtract' gives the anomaly in the data's units; "
            "'ratio' gives it as a percentage of the climatology, which is refused for quantities without "
            "a true zero."
    },
    returns="xr.DataArray of anomalies on the input's own grid and time axis. attrs record the exact "
        "baseline window and how many years stand behind each climatological value; state the baseline "
        "whenever the series is compared with a published index."
)
def compute_anomaly_vs_own_climatology(data: xr.DataArray, freq: str = "month",
                                       base_period: list = None,
                                       method: str = "subtract") -> xr.DataArray:
    if freq not in ("month", "season", "year"):
        raise ValueError(f"freq must be 'month', 'season' or 'year', got '{freq}'")
    if method not in ("subtract", "ratio"):
        raise ValueError(f"method must be 'subtract' or 'ratio', got '{method}'")
    if "time" not in data.dims:
        raise ValueError(
            f"An anomaly against a climatology needs a time axis to build the climatology from; this field "
            f"has dimensions {list(data.dims)}. A climatology product (WOA23, GLODAP) is already a "
            f"long-term mean and has no anomaly of its own."
        )

    baseline = data
    if base_period is not None:
        if len(base_period) != 2:
            raise ValueError(f"base_period must be ['start', 'end'], got {base_period}")
        start, end = str(base_period[0]), str(base_period[1])
        baseline = data.sel(time=slice(start, end))
        if baseline.sizes["time"] == 0:
            stamps = data["time"].values
            raise ValueError(
                f"The baseline window {start} to {end} contains none of this series, which runs "
                f"{str(stamps.min())[:7]} to {str(stamps.max())[:7]}. Load a period that covers the "
                f"baseline, or leave base_period out to use the record itself."
            )

    stamps = baseline["time"].values
    span = f"{str(stamps.min())[:7]} to {str(stamps.max())[:7]}"

    if freq == "year":
        climatology = baseline.mean("time", skipna=True)
        coverage = f"{baseline.sizes['time']} time steps"
    else:
        group = "time.month" if freq == "month" else "time.season"
        grouped = baseline.groupby(group)
        climatology = grouped.mean("time", skipna=True)
        sizes = {str(key): len(index) for key, index in grouped.groups.items()}
        expected = 12 if freq == "month" else 4
        if len(sizes) < expected:
            missing = expected - len(sizes)
            raise ValueError(
                f"The baseline {span} covers only {len(sizes)} of the {expected} {freq}s, so {missing} of "
                f"them would have no climatology and their anomalies would be NaN. Widen base_period to at "
                f"least one full year."
            )
        coverage = f"{min(sizes.values())} to {max(sizes.values())} samples per {freq}"

    if method == "ratio":
        # A percentage anomaly needs a ratio scale: degrees Celsius is not one,
        # and a percentage of a climatology near 0 C diverges and changes sign.
        units = data.attrs.get("units")
        if molar_basis(normalize_units(units)) is None and molar_basis(units) is None:
            raise ValueError(
                f"method='ratio' gives a percentage anomaly, which needs a quantity with a true zero. This "
                f"field is in '{units}', an interval scale (degrees Celsius, practical salinity, pH), where "
                f"a percentage of the climatology diverges near zero and flips sign below it. Use "
                f"method='subtract'."
            )
        climatology = climatology.where(climatology != 0)

    if freq == "year":
        difference = data - climatology
        anomaly = difference if method == "subtract" else difference / climatology * 100
    else:
        group = "time.month" if freq == "month" else "time.season"
        difference = data.groupby(group) - climatology
        if method == "subtract":
            anomaly = difference
        else:
            anomaly = difference.groupby(group) / climatology * 100
        anomaly = anomaly.drop_vars(group.split(".")[1], errors="ignore")

    # groupby arithmetic drops attrs, and the grid identity recorded there is
    # what lets a later area-weighted mean find the right cell areas.
    anomaly.attrs = dict(data.attrs)
    anomaly.attrs["anomaly_method"] = method
    anomaly.attrs["anomaly_freq"] = freq
    anomaly.attrs["climatology_reference"] = f"the field's own {freq} climatology over {span}"
    anomaly.attrs["climatology_baseline"] = span
    anomaly.attrs["baseline_coverage"] = coverage
    if freq == "season":
        anomaly.attrs["season_convention"] = (
            "meteorological (DJF/MAM/JJA/SON) - NOT the WOA23 calendar quarters used by "
            "compute_anomaly_vs_woa23"
        )
    long_name = data.attrs.get("long_name", data.name or "field")
    if method == "subtract":
        anomaly.attrs["long_name"] = f"{long_name} anomaly"
    else:
        anomaly.attrs["long_name"] = f"{long_name} percentage anomaly"
        anomaly.attrs["units"] = "percent"
    anomaly.name = f"{data.name}_anomaly" if data.name else "anomaly"
    return anomaly


@register_tools(
    name="convert_concentration_to_content",
    description="Integrate a concentration over grid-cell volume to get content in moles, for inventories "
        "such as the total oxygen or nitrate held in a region.\n"
        "Cell volumes come from the grid the data came from, preferring a model's own published cell "
        "volumes where they exist, since those alone account for partial cells at the sea floor.\n"
        "Both molar bases are handled: a mass-basis concentration needs a density to reach a volume, a "
        "volume-basis one needs none, and the units attribute decides which path is taken.\n"
        "Prefer the volume basis. A model's native mol m-3 field integrates exactly, with no density in "
        "the arithmetic at all. Handing over a umol/kg field without an explicit density applies this "
        "function's CONSTANT default instead, which silently undoes any in-situ density field used to "
        "make that field and biases the inventory by about one per cent.\n"
        "The field must still be on its original grid and levels: subsetting is fine, interpolation, "
        "regridding or averaging over depth is not.",
    args={
        "concentration": "xr.DataArray, required. A concentration with 'depth' and horizontal dimensions, "
            "optionally 'time', in either molar basis.",
        "density": "float, optional, default 1025. Seawater density in kg/m3, used only for a mass-basis "
            "input. It is ignored for a volume-basis field, which is already per unit volume."
    },
    returns="xr.DataArray of content in moles on the input's coordinates. Sum it over depth and the "
        "horizontal to get a region's inventory. attrs record where the volumes came from and, for a "
        "mass-basis input, the density assumed."
)
def convert_concentration_to_content(concentration: xr.DataArray, density: float = 1025.0) -> xr.DataArray:
    units = concentration.attrs.get("units")
    basis = molar_basis(units)
    if basis is None:
        raise ValueError(
            f"Cannot integrate a field whose units are '{units}': a content in moles needs a concentration "
            f"in umol kg-1 (WOA23, GLODAP) or mol m-3 (CMIP6). If the array lost its units attribute, set it "
            f"back before integrating."
        )

    volume = volumes_like(concentration)

    canonical = normalize_units(units)
    if basis == "mass":
        # umol kg-1 -> mol m-3 needs a density; mmol kg-1 is a thousand times that.
        per_volume = concentration * density * (1e-6 if canonical == "umol kg-1" else 1e-3)
        density_note = f"constant seawater density {density} kg/m3 (in-situ density not used)"
    else:
        # mol m-3 is already per unit volume; mmol m-3 is a thousandth of it.
        per_volume = concentration * (1.0 if canonical == "mol m-3" else 1e-3)
        density_note = "not needed: the input is already per unit volume"

    content = per_volume * volume
    content.attrs = dict(concentration.attrs)
    content.attrs.update({
        "units": "mol",
        "long_name": f"{concentration.attrs.get('long_name', concentration.name or 'tracer')} content",
        "description": "Concentration integrated over grid-cell volume",
        "original_units": units,
        "molar_basis_of_input": basis,
        "density_assumption": density_note,
        "volume_source": volume.attrs.get(
            "volume_source", volume.attrs.get("area_source", "cell area x layer thickness")
        ),
    })
    if concentration.name:
        content.name = f"{concentration.name}_content"
    return content


@register_tools(
    name="compute_vertical_gradient",
    description="Compute the vertical gradient of a field with respect to depth, layer by layer.\n"
        "The value is interpolated onto each layer's interfaces and differenced across them, so an uneven "
        "vertical grid is handled correctly: differencing neighbouring levels directly would attribute a "
        "near-surface difference and a deep-ocean one to the same nominal spacing. The topmost and "
        "bottommost gradients are one-sided.\n"
        "Layer boundaries come from the grid the data came from, so the field must keep its original "
        "levels.",
    args={
        "data": "xr.DataArray, required. A field with a 'depth' dimension, plus any other dimensions."
    },
    returns="xr.DataArray of d(value)/d(depth) at the same levels and dimensions as the input, in the "
        "field's units per metre. Depth is positive downward, so a positive gradient means the quantity "
        "increases with depth."
)
def compute_vertical_gradient(data: xr.DataArray) -> xr.DataArray:
    """
    Compute the vertical gradient layer by layer, across each layer's own bounds.

    The value is linearly interpolated from the two neighbouring level centres
    onto each interface, and the gradient of a layer is the difference across
    its two interfaces divided by its thickness. This respects an uneven
    vertical grid, which every one of these products has: differencing
    neighbouring levels directly would attribute a 5 m difference near the
    surface and a 500 m one in the deep ocean to the same nominal spacing.
    The topmost and bottommost interfaces take the adjacent level's value, so
    the gradient there is one-sided.
    """
    if 'depth' not in data.dims:
        raise ValueError("Input DataArray must have a 'depth' dimension.")

    # The layer boundaries come from the vertical grid the data itself came from.
    # They are not shared between products: GLODAP's levels are far thicker than
    # WOA23's at the same depths, and each CMIP6 model has its own, so taking a
    # default set would rescale every gradient.
    depth_vals = data['depth'].values
    depth_bnds = layer_bounds_like(data).values

    interfaces = np.empty(len(depth_bnds) + 1, dtype=np.float64)
    interfaces[0] = depth_bnds[0, 0]
    for i in range(len(depth_bnds)):
        interfaces[i + 1] = depth_bnds[i, 1]

    # The span each difference is taken over. For an interior layer that is the
    # layer thickness, but the top and bottom differences run only from the
    # outermost level CENTRE to the first interior interface, which is part of a
    # layer. Dividing those by the full thickness understates the gradient at
    # both ends: on a perfectly linear profile the deepest layer comes out six
    # times too small.
    span = (depth_bnds[:, 1] - depth_bnds[:, 0]).astype("float64")
    if len(span) > 1:
        span[0] = interfaces[1] - depth_vals[0]
        span[-1] = depth_vals[-1] - interfaces[-2]
    dz_safe = np.where(span <= 0, np.nan, span)

    transposed = data.transpose(..., 'depth')
    v_np = transposed.values.astype("float64")  # shape: (..., depth)
    n_depth = v_np.shape[-1]

    v_int = np.empty(v_np.shape[:-1] + (n_depth + 1,), dtype=np.float64)
    v_int[..., 0] = v_np[..., 0]
    v_int[..., -1] = v_np[..., -1]
    for k in range(1, n_depth):
        z_if = interfaces[k]
        z_up, z_dn = depth_vals[k - 1], depth_vals[k]
        denom = (z_dn - z_up)
        if denom == 0:
            v_int[..., k] = 0.5 * (v_np[..., k - 1] + v_np[..., k])
        else:
            w_up = (z_dn - z_if) / denom
            w_dn = (z_if - z_up) / denom
            v_int[..., k] = w_up * v_np[..., k - 1] + w_dn * v_np[..., k]

    dv = v_int[..., 1:] - v_int[..., :-1]  # (..., depth)
    dvdz = dv / dz_safe.reshape((1,) * (v_np.ndim - 1) + (-1,))

    grad = xr.DataArray(
        dvdz,
        dims=transposed.dims,
        coords={name: coord for name, coord in transposed.coords.items()
                if set(coord.dims).issubset(transposed.dims)},
        attrs={
            "long_name": f"Vertical gradient of {data.name or 'variable'}",
            "units": f"{data.attrs.get('units', '')} / m",
            "method": "interface",
            "positive": "down"
        }
    ).transpose(*data.dims)

    return carry_provenance(grad, data)


@register_tools(
    name="compute_ratio",
    description="Compute the element-wise ratio of two fields, for stoichiometric quantities such as N:P.\n"
        "Two concentrations on different molar bases are refused rather than divided, because their ratio "
        "would not be a molar ratio; convert one of them first. Two concentrations on the same basis but "
        "with different prefixes are rescaled so the ratio is dimensionless.",
    args={
        "data1": "xr.DataArray, required. The numerator.",
        "data2": "xr.DataArray, required. The denominator. Cells where it is zero come back missing rather "
            "than infinite."
    },
    returns="xr.DataArray of the ratio at every cell. Note that averaging this field is not the same as "
        "dividing two averages: the mean of a ratio weights every cell's composition equally, while the "
        "ratio of means weights each cell by its content and is pulled towards the richest water. Both are "
        "defensible, so state which was computed; this tool gives the first, and the second is a weighted "
        "mean of each field before the ratio. attrs record the units of both inputs and whether the result "
        "is a true molar ratio."
)
def compute_ratio(data1: xr.DataArray, data2: xr.DataArray) -> xr.DataArray:
    units1, units2 = data1.attrs.get("units"), data2.attrs.get("units")
    basis1, basis2 = molar_basis(units1), molar_basis(units2)
    canonical1, canonical2 = normalize_units(units1), normalize_units(units2)

    scale = 1.0
    if basis1 and basis2:
        if basis1 != basis2:
            raise ValueError(
                f"'{units1}' is a per-{basis1} concentration and '{units2}' is a per-{basis2} one, so their "
                f"ratio is not a molar ratio: WOA23 and GLODAP report per unit mass (umol kg-1) while CMIP6 "
                f"reports per unit volume (mol m-3). Convert one of them with transform_unit first — left as "
                f"they are, a Redfield N:P of about 16 comes out near 0.008."
            )
        # Same basis, possibly different prefixes (umol kg-1 against mmol kg-1,
        # or mol m-3 against mmol m-3); rescale so the ratio is dimensionless.
        prefix = {"umol kg-1": 1.0, "mmol kg-1": 1000.0, "mmol m-3": 1.0, "mol m-3": 1000.0}
        scale = prefix.get(canonical1, 1.0) / prefix.get(canonical2, 1.0)
        ratio_units, ratio_basis = "1", f"both per unit {basis1}"
    elif canonical1 is not None and canonical1 == canonical2:
        ratio_units, ratio_basis = "1", f"both in {units1}"
    else:
        ratio_units = f"({units1 or '?'})/({units2 or '?'})"
        ratio_basis = "not a molar ratio; the two fields carry different kinds of unit"

    ratio = xr.where(data2 != 0, (data1 * scale) / data2, np.nan)
    ratio.attrs.update({
        'long_name': f'Ratio of {data1.name} to {data2.name}',
        'description': 'Computed as (data1 / data2)',
        'units': ratio_units,
        'ratio_basis': ratio_basis,
        'numerator_units': units1,
        'denominator_units': units2,
    })
    ratio.name = f"{data1.name}_to_{data2.name}_ratio" if data1.name and data2.name else "ratio"
    return carry_provenance(ratio, data1)


@register_tools(
    name="identify_oxygen_minimum_zone",
    description="Identify oxygen minimum zones from a dissolved oxygen field and measure their total "
        "vertical thickness.\n"
        "Every layer below the threshold counts, and the thickness is the summed thickness of those "
        "layers, giving a map of how much of the water column is oxygen-depleted. Thicknesses come from "
        "the grid the data came from, so the field must keep its original levels. For a CMIP6 model each "
        "cell counts the water it actually holds (volcello / areacello), so a partial cell at the sea "
        "floor is not counted as a full layer; this matters most for GFDL-ESM4, whose archived levels "
        "carry values in bottom cells that are only partly water.\n"
        "The threshold must be in the same units as the field. A mass-basis threshold applied to a "
        "volume-basis field marks the entire ocean as depleted; convert one of them first.",
    args={
        "o2_data": "xr.DataArray, required. Dissolved oxygen with a 'depth' dimension, optionally with "
            "time and horizontal dimensions.",
        "threshold": "float, optional, default 20.0. The concentration below which a layer counts, in the "
            "data's own units. Common alternatives are 60 for hypoxia and 5 for suboxia in umol kg-1."
    },
    returns="tuple (omz_thickness, omz_mask):\n"
        "- omz_thickness: total thickness in metres of the layers meeting the criterion at each horizontal "
        "point, missing where there is no ocean.\n"
        "- omz_mask: boolean array with the input's dimensions, True inside the zone.\n"
        "For the volume rather than the thickness, use compute_exceedance_fraction with volume weighting, "
        "which returns the cubic metres directly. attrs record the threshold and where the layer "
        "thicknesses came from."
)
def identify_oxygen_minimum_zone(o2_data: xr.DataArray, threshold: float = 20.0) -> tuple[xr.DataArray, xr.DataArray]:
    units = o2_data.attrs.get("units")
    if molar_basis(units) == "volume":
        raise ValueError(
            f"This oxygen field is in '{units}', a per-volume unit, but the threshold {threshold} is read in "
            f"the data's own units and the conventional OMZ thresholds (20, 60, 5) are per mass (umol kg-1). "
            f"Convert the field with transform_unit(data, '{units}', 'umol/kg') first, or pass a threshold "
            f"already in {units} ({threshold} umol kg-1 is about {threshold * 1.025e-3:.4g} mol m-3). "
            f"Left as it is, every ocean cell would be counted as an oxygen minimum zone."
        )

    thickness = cell_thickness_like(o2_data)
    omz_mask = o2_data < threshold

    omz_thickness = (thickness * omz_mask).sum(dim="depth", skipna=True)
    # Cells with no ocean at all must stay NaN rather than reading as a zero-metre
    # OMZ, so the result is masked wherever the column holds no data.
    omz_thickness = omz_thickness.where(o2_data.notnull().any(dim="depth"))

    omz_thickness.attrs = {
        "long_name": "Oxygen minimum zone thickness",
        "units": "m",
        "description": f"Total thickness of the layers where oxygen < {threshold} {units}",
        "threshold": threshold,
        "threshold_units": units,
        "thickness_source": thickness.attrs.get("bounds_source", ""),
    }
    omz_mask.attrs = {
        "long_name": "Oxygen minimum zone mask",
        "description": f"True where oxygen < {threshold} {units}",
        "threshold": threshold,
        "threshold_units": units,
    }
    return carry_provenance(omz_thickness, o2_data), carry_provenance(omz_mask, o2_data)


@register_tools(
    name="compute_potential_density_anomaly",
    description="Compute density anomaly (sigma) from temperature and practical salinity using TEOS-10.\n"
        "Sigma referenced to a fixed pressure is what labels a water mass, because it removes the "
        "compression that makes any water denser with depth. Referenced to each cell's own pressure it is "
        "the in-situ density, which is what mass-per-volume arithmetic needs — converting a concentration "
        "between molar bases, or a volume into a mass.\n"
        "The result is an ANOMALY (roughly 25-35 kg/m3) and carries the metadata that every tool taking a "
        "density argument reads. Hand it over unchanged: those tools add the 1000 themselves when they need "
        "an absolute density, so adding it here makes them do it twice.\n"
        "The two temperature conventions take different routes to Conservative Temperature, and using the "
        "wrong one biases density increasingly with depth, so the kind of temperature supplied is read "
        "from its metadata and recorded in the result.",
    args={
        "temp_data": "xr.DataArray, required. Temperature in degrees Celsius, in-situ or potential, with a "
            "'depth' coordinate in metres.",
        "sal_data": "xr.DataArray, required. Practical salinity (PSS-78), on the same dimensions as "
            "temp_data.",
        "ref_pressure": "float or None, optional, default 0. Reference pressure in dbar: 0 gives sigma0, "
            "the usual water-mass label at the surface, and 2000 or 4000 are used for the deep ocean. Pass "
            "None for the in-situ density anomaly, evaluated at each cell's own pressure — use that for "
            "unit conversion, and a fixed reference whenever water masses are compared.",
        "temperature_type": "str, optional, default 'auto'. Whether temp_data is 'in_situ' or 'potential' "
            "temperature; 'auto' reads it from the metadata."
    },
    returns="xr.DataArray of sigma in kg/m3, i.e. density minus 1000, on the input's coordinates. attrs "
        "record the reference pressure and which temperature convention was detected."
)
def compute_potential_density_anomaly(temp_data: xr.DataArray, sal_data: xr.DataArray, ref_pressure=0.0,
                                      temperature_type: str = "auto") -> xr.DataArray:
    if temp_data.dims != sal_data.dims:
        raise ValueError(f"Temperature and salinity must have same dimensions. "
                        f"Got temp: {temp_data.dims}, sal: {sal_data.dims}")

    if 'depth' not in temp_data.coords:
        raise ValueError("Temperature and salinity data must have 'depth' coordinate")

    # Depth, latitude and longitude are laid out on their own axes first, so
    # that a full (depth, lat, lon) field broadcasts. The 2-D lat/lon of the
    # curvilinear models are handled by the same reshape.
    depth_arr = _broadcastable(temp_data["depth"], temp_data)
    lat_coord = temp_data.attrs.get("lat_coord", "lat")
    lon_coord = temp_data.attrs.get("lon_coord", "lon")
    lat_arr = _broadcastable(temp_data[lat_coord], temp_data) if lat_coord in temp_data.coords else 0.0
    lon_arr = _broadcastable(temp_data[lon_coord], temp_data) if lon_coord in temp_data.coords else 0.0

    salinity_units = normalize_units(sal_data.attrs.get("units"))
    if salinity_units not in (None, "1"):
        raise ValueError(
            f"Salinity must be practical salinity (PSS-78, dimensionless), but sal_data is in "
            f"'{sal_data.attrs.get('units')}'. WOA23 reports it as '1', CMIP6 as '0.001'; both are PSS-78."
        )
    temperature_units = normalize_units(temp_data.attrs.get("units"))
    if temperature_units not in (None, "degC"):
        raise ValueError(
            f"Temperature must be in degrees Celsius, but temp_data is in "
            f"'{temp_data.attrs.get('units')}'. Convert it with transform_unit first."
        )

    kind, kind_reason = temperature_kind(temp_data, temperature_type)

    pressure = gsw.p_from_z(-depth_arr, lat_arr)
    absolute_salinity = gsw.SA_from_SP(sal_data.values, pressure, lon_arr, lat_arr)
    # CMIP6's thetao is potential temperature and WOA23's is in-situ; TEOS-10
    # takes a different route from each to Conservative Temperature, and using
    # the wrong one biases density increasingly with depth.
    if kind == "potential":
        conservative_temperature = gsw.CT_from_pt(absolute_salinity, temp_data.values)
    else:
        conservative_temperature = gsw.CT_from_t(absolute_salinity, temp_data.values, pressure)

    # rho(SA, CT, p_ref) - 1000 is the potential density anomaly referenced to
    # p_ref, and reduces to sigma0 at p_ref = 0. Evaluated at each cell's own
    # pressure instead, it is the in-situ density anomaly: not a water-mass
    # label, but the density the water actually has, which is what converts
    # between a concentration per unit volume and one per unit mass.
    evaluation_pressure = pressure if ref_pressure is None else ref_pressure
    sigma_theta = gsw.rho(absolute_salinity, conservative_temperature, evaluation_pressure) - 1000.0

    rho_data = xr.DataArray(
        sigma_theta,
        coords=temp_data.coords,
        dims=temp_data.dims,
        attrs={
            'long_name': ('In-situ density anomaly' if ref_pressure is None
                          else f'Potential density anomaly (sigma_{int(ref_pressure)})'),
            'units': 'kg/m³',
            'reference_pressure': ("each cell's own in-situ pressure" if ref_pressure is None
                                   else f'{ref_pressure} dbar'),
            'description': (
                'In-situ density anomaly (rho - 1000 kg/m³), evaluated at each cell\'s own pressure'
                if ref_pressure is None else
                f'Potential density anomaly (rho - 1000 kg/m³) referenced to {ref_pressure} dbar'),
            'method': 'TEOS-10 (GSW library)',
            'source_temperature': f'{"potential" if kind == "potential" else "in-situ"} temperature '
                                  f'({kind_reason})',
            'source_salinity': 'practical salinity (PSS-78)'
        }
    )

    return carry_provenance(rho_data, temp_data)


def _first_crossing_depth(data: xr.DataArray, depth: np.ndarray, threshold: float,
                          ref_depth: float) -> xr.DataArray:
    """
    Depth at which `data` first reaches `threshold`, interpolated between levels.

    Returning the discrete level instead quantises the answer onto the vertical
    grid, and these grids are coarse where mixed layers actually sit: WOA23
    jumps 25 m at a time below 100 m, so a deep winter mixed layer could only
    ever come back as one of a handful of values. Linear interpolation between
    the bracketing levels is what the community definition asks for
    (de Boyer Montegut et al. 2004).

    Vectorised over every other dimension: the profile-by-profile Python loop
    this replaces took one pass per horizontal cell, which is millions of passes
    on a global monthly field.
    """
    transposed = data.transpose(..., "depth")
    values = transposed.values
    z = np.asarray(depth, dtype="float64")

    # Only levels at or below the reference depth can define the crossing.
    eligible = z >= ref_depth
    exceeds = (values >= threshold) & eligible
    found = exceeds.any(axis=-1)

    first = np.argmax(exceeds, axis=-1)
    # Interpolate between the level that crosses and the one above it.
    upper_index = np.clip(first - 1, 0, len(z) - 1)
    value_at = np.take_along_axis(values, first[..., None], axis=-1)[..., 0]
    value_above = np.take_along_axis(values, upper_index[..., None], axis=-1)[..., 0]
    z_at, z_above = z[first], z[upper_index]

    span = value_at - value_above
    with np.errstate(invalid="ignore", divide="ignore"):
        fraction = np.where(np.abs(span) > 0, (threshold - value_above) / span, 0.0)
    fraction = np.clip(np.nan_to_num(fraction, nan=0.0), 0.0, 1.0)
    crossing = z_above + fraction * (z_at - z_above)

    # The very first eligible level already exceeding the threshold leaves
    # nothing to interpolate from, and a mixed layer cannot be shallower than
    # the reference depth it is measured against.
    crossing = np.where(first == 0, z_at, crossing)
    crossing = np.maximum(crossing, ref_depth)
    crossing = np.where(found, crossing, np.nan)

    remaining = [d for d in transposed.dims if d != "depth"]
    return xr.DataArray(
        crossing,
        coords={k: v for k, v in data.coords.items() if k != "depth" and set(v.dims) <= set(remaining)},
        dims=remaining,
    )


@register_tools(
    name="compute_mixed_layer_depth",
    description="Compute mixed layer depth from a density or temperature threshold.\n"
        "The MLD is the depth at which potential density, or temperature, first departs from its value at "
        "a near-surface reference depth by more than the threshold. The defaults are the de Boyer "
        "Montegut criterion, which is what observational MLD climatologies use, so the two are directly "
        "comparable. The crossing depth is interpolated between the bracketing levels rather than snapped "
        "to one, because these vertical grids are coarse where mixed layers sit.\n"
        "COMPARING WITH A MODEL'S OWN mlotst: pass ref_depth equal to the model's shallowest level, since "
        "CMIP6 mlotst is referenced to the ocean surface, not to 10 m. That removes the largest definitional "
        "difference but NOT all of it, so the result is close to mlotst, never identical, and the residual "
        "must not be read as model error. The CMIP6 protocol (Griffies et al. 2016, H24) computes mlotst "
        "ONLINE at every model time step, as a buoyancy difference between a surface parcel moved down to "
        "each level's own pressure and the local water (0.0003 m s-2, about 0.03 kg/m3), and archives the "
        "monthly MEAN of those instantaneous depths. This tool applies a sigma0 difference to whatever "
        "profile it is given, and a monthly field is itself a time-averaged profile, so the operator is "
        "biased shallow against mlotst for the same reason given below. GFDL-ESM4 also computes mlotst on "
        "its native 75 hybrid layers, while its archived temperature and salinity were remapped to 35 z "
        "levels, which adds a resolution difference only that model has.\n"
        "COLUMNS THAT NEVER REACH THE THRESHOLD come back missing. The models' own mlotst does not agree "
        "on these: GFDL-ESM4 and MPI-ESM1-2-LR report roughly the full column depth there, CESM2 reports 0 "
        "(checked on March 1990). Exclude them from any comparison with mlotst and report how much area "
        "was excluded.\n"
        "ORDER OF OPERATIONS: the MLD of an averaged profile is not the average of the MLDs and is "
        "systematically shallower, because averaging smooths the pycnocline away. Compute it step by step "
        "and average the result; never average temperature and salinity first.",
    args={
        "temp_data": "xr.DataArray, required. Temperature in degrees Celsius with a 'depth' coordinate in "
            "metres.",
        "sal_data": "xr.DataArray, required. Practical salinity (PSS-78) on identical dimensions.",
        "ref_depth": "float, optional, default 10.0. The near-surface reference depth in metres. Use the "
            "field's shallowest level when matching a model's own surface-referenced output.",
        "method": "str, optional, default 'density'. 'density' uses the density threshold, 'temperature' "
            "the temperature threshold, 'both' takes whichever is shallower.",
        "density_threshold": "float, optional, default 0.03. Density difference in kg/m3 defining the base "
            "of the mixed layer.",
        "temp_threshold": "float, optional, default 0.2. Temperature difference in degrees Celsius, for "
            "the temperature criterion.",
        "temperature_type": "str, optional, default 'auto'. Whether temp_data is in-situ or potential "
            "temperature; 'auto' reads it from the metadata."
    },
    returns="xr.DataArray of mixed layer depth in metres on all non-depth dimensions of the input, "
        "interpolated between levels. Missing where the threshold is never reached — a column mixed all "
        "the way down, or one with no data."
)
def compute_mixed_layer_depth(temp_data: xr.DataArray, sal_data: xr.DataArray,
                              ref_depth: float = 10.0, method: str = 'density',
                              density_threshold: float = 0.03, temp_threshold: float = 0.2,
                              temperature_type: str = "auto") -> xr.DataArray:
    if temp_data.dims != sal_data.dims:
        raise ValueError(f"Temperature and salinity must have same dimensions. "
                        f"Got temp: {temp_data.dims}, sal: {sal_data.dims}")
    
    if 'depth' not in temp_data.coords:
        raise ValueError("Data must have 'depth' coordinate")
    
    if method not in ['density', 'temperature', 'both']:
        raise ValueError(f"Method must be 'density', 'temperature', or 'both', got '{method}'")
    # The third positional parameter is ref_depth, not method, so a call written
    # as compute_mixed_layer_depth(t, s, 'density') lands the method name here
    # and would otherwise fail deep inside numpy with an unreadable dtype error.
    if not isinstance(ref_depth, (int, float)) or isinstance(ref_depth, bool):
        raise TypeError(
            f"ref_depth must be a depth in metres, got {ref_depth!r}. Note the parameter order: "
            f"compute_mixed_layer_depth(temp_data, sal_data, ref_depth, method, ...) — pass the criterion "
            f"by keyword, as method='{ref_depth}', if that is what you meant."
        )

    depth = temp_data['depth'].values
    if ref_depth < depth.min() or ref_depth > depth.max():
        raise ValueError(f"Reference depth {ref_depth}m is outside data depth range "
                        f"[{depth.min()}, {depth.max()}]m")
    
    sigma_theta = compute_potential_density_anomaly(
        temp_data, sal_data, ref_pressure=0.0, temperature_type=temperature_type
    )

    ref_idx = np.argmin(np.abs(depth - ref_depth))
    actual_ref_depth = depth[ref_idx]

    sigma_ref = sigma_theta.isel(depth=ref_idx)
    temp_ref = temp_data.isel(depth=ref_idx)

    if method in ['density', 'both']:
        delta_sigma = sigma_theta - sigma_ref
        mld_density = _first_crossing_depth(
            delta_sigma, depth, threshold=density_threshold, ref_depth=actual_ref_depth
        )
    if method in ['temperature', 'both']:
        delta_temp = np.abs(temp_data - temp_ref)
        mld_temperature = _first_crossing_depth(
            delta_temp, depth, threshold=temp_threshold, ref_depth=actual_ref_depth
        )
    
    if method == 'density':
        mld = mld_density
        method_description = f'Density threshold: Δσ_θ = {density_threshold} kg/m³'
    elif method == 'temperature':
        mld = mld_temperature
        method_description = f'Temperature threshold: ΔT = {temp_threshold} °C'
    else:  # both
        mld = xr.where(
            np.isnan(mld_density), 
            mld_temperature,
            xr.where(
                np.isnan(mld_temperature),
                mld_density,
                np.minimum(mld_density, mld_temperature)
            )
        )
        method_description = f'Minimum of density (Δσ_θ = {density_threshold} kg/m³) and temperature (ΔT = {temp_threshold} °C)'
    
    mld_data = xr.DataArray(
        mld,
        coords={k: v for k, v in temp_data.coords.items() if k != 'depth'},
        dims=[d for d in temp_data.dims if d != 'depth'],
        attrs={
            'long_name': 'Mixed layer depth',
            'standard_name': 'ocean_mixed_layer_thickness',
            'units': 'm',
            'reference_depth': f'{actual_ref_depth} m',
            'method': method_description
        }
    )

    return carry_provenance(mld_data, temp_data)


@register_tools(
    name="compute_nutricline_depth",
    description="Find the nutricline: the depth at which nutrient concentration increases most steeply "
        "with depth.\n"
        "It is the depth of the maximum vertical gradient, computed across each layer's own bounds so an "
        "uneven vertical grid is handled correctly. Where the profile is smooth — in the subtropical gyres "
        "especially — the position of the maximum can jump between levels, so read a nutricline map "
        "together with the gradient itself rather than alone.\n"
        "Restrict the search with min_depth and max_depth when only a particular part of the water column "
        "is of interest.",
    args={
        "data": "xr.DataArray, required. Nutrient concentration with a 'depth' dimension, plus any other "
            "dimensions.",
        "min_depth": "float, optional, default 0. Ignore levels shallower than this when looking for the "
            "maximum gradient.",
        "max_depth": "float, optional, default None (every level). Ignore levels deeper than this."
    },
    returns="xr.DataArray of nutricline depth in metres on the input's dimensions with 'depth' removed. "
        "Missing where the gradient is never positive or the column has no data."
)
def compute_nutricline_depth(data: xr.DataArray, min_depth: float = 0.0, max_depth: float = None) -> xr.DataArray:
    if 'depth' not in data.dims:
        raise ValueError("Input DataArray must have a 'depth' dimension.")
    
    gradient = compute_vertical_gradient(data)
    depth_vals = gradient['depth'].values

    depth_mask = np.ones(len(depth_vals), dtype=bool)
    if min_depth is not None:
        depth_mask &= (depth_vals >= min_depth)
    if max_depth is not None:
        depth_mask &= (depth_vals <= max_depth)

    gradient_masked = gradient.copy()
    if not depth_mask.all():
        mask_full = xr.DataArray(
            depth_mask,
            coords={'depth': gradient['depth']},
            dims=['depth']
        )
        gradient_masked = gradient_masked.where(mask_full, np.nan)

    # argmax raises on an all-NaN profile, and over any ocean field most columns
    # are land, so the land columns are filled before the search and masked out
    # again afterwards.
    all_nan = gradient_masked.isnull().all(dim='depth')
    max_idx = gradient_masked.fillna(-np.inf).argmax(dim='depth')
    nutricline_depth = gradient['depth'].isel(depth=max_idx)

    max_gradient = gradient_masked.max(dim='depth', skipna=True)
    nutricline_depth = nutricline_depth.where(
        (max_gradient > 0) & (~all_nan),
        np.nan
    )
    
    nutricline_depth.attrs = {
        'long_name': f'Nutricline depth ({data.name or "nutrient"})',
        # No CF standard name exists for the nutricline, so none is claimed;
        # 'ocean_mixed_layer_thickness' on the MLD below is a real one.
        'units': 'm',
        'description': 'Depth of maximum vertical gradient in nutrient concentration',
        'method': 'maximum_gradient',
        'positive': 'down',
        'source_variable': data.name or 'unknown'
    }
    nutricline_depth.name = f'{data.name}_nutricline_depth' if data.name else 'nutricline_depth'

    # Selecting from the depth coordinate builds the answer out of that
    # coordinate, so it arrives carrying 'depth' (now meaningless: the values ARE
    # the depths) and without the horizontal coordinates it was indexed by. The
    # 2-D lat/lon of a curvilinear model are lost that way, and with them any
    # chance of regridding or plotting the result.
    nutricline_depth = nutricline_depth.drop_vars('depth', errors='ignore')
    nutricline_depth = nutricline_depth.assign_coords({
        name: coord for name, coord in data.coords.items()
        if name != 'depth' and set(coord.dims).issubset(set(nutricline_depth.dims))
    })

    return carry_provenance(nutricline_depth, data)


def _on_woa_grid(data: xr.DataArray) -> bool:
    """Whether a field sits on the WOA23 1-degree cell centres (x.5 spacing)."""
    if data.attrs.get("grid_type") == "curvilinear":
        return False
    if "lat" not in data.coords or "lon" not in data.coords:
        return False
    lat, lon = np.atleast_1d(data["lat"].values), np.atleast_1d(data["lon"].values)
    return bool(
        np.allclose(np.abs(lat % 1), 0.5, atol=1e-3) and np.allclose(np.abs(lon % 1), 0.5, atol=1e-3)
    )


def _select_region(data: xr.DataArray, lat_min, lat_max, lon_min, lon_max) -> xr.DataArray:
    """
    Cut a lat/lon box out of a field, including one that crosses the dateline.

    Masking and dropping the empty rows, rather than slicing and concatenating,
    keeps the longitude axis in the -180~180 convention and strictly ascending.
    A concatenated Nino4 box would otherwise run 160..180, -180..-150, which is
    non-monotonic, so nothing could be selected from it afterwards and the cell
    areas could not be matched to it.

    The curvilinear models are handled by the same mask. Their lat/lon are 2-D
    fields over index dimensions, so the box picks out a ragged set of cells
    rather than a rectangle; dropping the index rows and columns it never
    touches trims the array to the box's bounding block and leaves the cells
    outside the box as NaN, which is what an area-weighted mean over it needs —
    the weights are masked by the field itself. This keeps a regional index on
    the model's OWN grid and its own areacello, with no interpolation in the way.
    """
    curvilinear = data.attrs.get("grid_type") == "curvilinear"
    lat_name = data.attrs.get("lat_coord", "lat") if curvilinear else "lat"
    lon_name = data.attrs.get("lon_coord", "lon") if curvilinear else "lon"
    for name in (lat_name, lon_name):
        if name not in data.coords:
            raise ValueError(
                f"This field has no '{name}' coordinate, so a latitude/longitude box cannot be cut out of "
                f"it. Load it with load_regular_data or load_climatology_data."
            )

    in_lat = (data[lat_name] >= lat_min) & (data[lat_name] <= lat_max)
    if lon_max < lon_min:                      # crosses the dateline, e.g. Nino4
        in_lon = (data[lon_name] >= lon_min) | (data[lon_name] <= lon_max)
    else:
        in_lon = (data[lon_name] >= lon_min) & (data[lon_name] <= lon_max)
    # The models' 2-D lat/lon arrive chunked, and xarray cannot drop rows with a
    # lazy boolean mask because the resulting shape is not known until it is
    # computed. The mask is one horizontal slice, so computing it here is cheap.
    mask = in_lat & in_lon
    if hasattr(mask.data, "compute"):
        mask = mask.compute()
    region = data.where(mask, drop=True)

    if curvilinear:
        y_dim, x_dim = data.attrs.get("y_dim"), data.attrs.get("x_dim")
        if region.sizes.get(y_dim, 0) == 0 or region.sizes.get(x_dim, 0) == 0:
            raise ValueError(
                f"The box lat [{lat_min}, {lat_max}], lon [{lon_min}, {lon_max}] selects no cells from "
                f"{data.attrs.get('model', 'this model')}'s grid."
            )
        return region

    if region.sizes.get("lat", 0) == 0 or region.sizes.get("lon", 0) == 0:
        raise ValueError(
            f"The box lat [{lat_min}, {lat_max}], lon [{lon_min}, {lon_max}] selects no cells from this "
            f"field, which covers lat [{float(data.lat.min())}, {float(data.lat.max())}] and lon "
            f"[{float(data.lon.min())}, {float(data.lon.max())}]."
        )
    return region


@register_tools(
    name="compute_enso_nino_index",
    description="Compute an ENSO Nino index: the area-weighted sea surface temperature anomaly over one of "
        "the Nino regions, as a monthly series.\n"
        "Use it to derive the index from a gridded field, which is what a model needs, or to reproduce the "
        "calculation from an observational product. Where an official observed index is already available "
        "as a series, there is nothing to compute.\n"
        "Any grid is accepted, curvilinear included: the region is cut out with a mask and averaged with "
        "the field's own cell areas, so do not regrid first.\n"
        "The baseline matters. An index built on one baseline period differs from one built on another by "
        "a constant offset, even when the variability matches month for month; set base_period to the "
        "baseline of whatever series it will be compared against.",
    args={
        "data": "xr.DataArray, required. Temperature with 'time' and horizontal dimensions, optionally "
            "'depth', in degrees Celsius, on its original grid. Pass the global field, or one already "
            "subset to a box containing the region.",
        "nino_region": "str, optional, default 'nino3.4'. 'nino1+2', 'nino3', 'nino3.4' or 'nino4'. "
            "Nino3.4 is the standard choice for ENSO monitoring.",
        "depth_level": "float or None, optional, default None. Only used when the field has a depth "
            "dimension. None takes the shallowest level, which is how CMIP6 defines SST (tos is the "
            "top-cell temperature) and so is the like-for-like substitute where tos is missing. A number "
            "averages every level shallower than it instead, thickness-weighted; at 10 m that is two "
            "levels in GFDL-ESM4 but one in CESM2 and MPI-ESM1-2-LR, so it is not the same quantity "
            "across models.",
        "anomaly_freq": "str, optional, default 'month'. Which climatology each step is matched against: "
            "'month' (the usual choice, removing the seasonal cycle), 'season' or 'year'.",
        "climatology": "str, optional, default 'auto'. 'auto' removes the field's own climatology for "
            "CMIP6 model output, always, and uses the observed WOA23 one only for an observational field "
            "on its grid; 'woa23' forces the observed one and fails if the grid does not match; 'self' "
            "always removes the field's own. Forcing 'woa23' on a model puts the model's mean-state bias "
            "into the index, which is then not an ENSO index.",
        "base_period": "list[str, str], optional. The baseline window for the field's own climatology, as "
            "['start', 'end']. Passing it forces the self-climatology, since observed climate normals are "
            "fixed periods. Default None: the whole record loaded is the baseline."
    },
    returns="xr.DataArray, a monthly series of the index in degrees Celsius with a single 'time' "
        "dimension. attrs record the region box, the climatology used and its baseline, which must be "
        "stated when the values are compared with a published series."
)
def compute_enso_nino_index(data: xr.DataArray, nino_region: str = "nino3.4", depth_level: float = None,
                            anomaly_freq: str = "month", climatology: str = "auto",
                            base_period: list = None) -> xr.DataArray:
    # (lat_min, lat_max, lon_min, lon_max), matching the boxes in configs/datasets.yaml
    nino_regions = {
        "nino1+2": (-10, 0, -90, -80),
        "nino3": (-5, 5, -150, -90),
        "nino3.4": (-5, 5, -170, -120),
        "nino4": (-5, 5, 160, -150),   # crosses the dateline
    }
    if nino_region not in nino_regions:
        raise ValueError(f"nino_region must be one of {list(nino_regions.keys())}, got '{nino_region}'")
    if climatology not in ("auto", "woa23", "self"):
        raise ValueError(f"climatology must be 'auto', 'woa23' or 'self', got '{climatology}'")
    if "time" not in data.dims:
        raise ValueError("An ENSO index needs a time series; this field has no 'time' dimension.")

    if base_period is not None and climatology == "woa23":
        raise ValueError(
            "base_period sets the window for the field's OWN climatology; WOA23's climate normals are "
            "fixed periods chosen with the 'period' of the WOA23 tools. Use climatology='self' (or leave "
            "it at 'auto') to build the index on your own baseline."
        )

    lat_min, lat_max, lon_min, lon_max = nino_regions[nino_region]
    data_region = _select_region(data, lat_min, lat_max, lon_min, lon_max)

    # A model's index is its anomaly from its OWN climatology. GFDL-ESM4 happens
    # to sit on the WOA23 1-degree grid, and subtracting the observed climatology
    # from it would put the model's monthly mean-state bias into every value of
    # the index (and drop its 2.5 m level, which WOA23 does not have).
    is_model = str(data.attrs.get("source_dataset", "")).startswith("cmip6_")
    use_woa = (climatology == "woa23") or (
        climatology == "auto" and base_period is None and not is_model and _on_woa_grid(data)
    )
    if use_woa:
        anomaly = compute_anomaly_vs_woa23(data=data_region, var="temp", freq=anomaly_freq,
                                          method="subtract")
        baseline = anomaly.attrs.get("climatology_reference", "WOA23")
    else:
        if climatology == "woa23":
            raise ValueError(
                "climatology='woa23' needs the field to be on the WOA23 1-degree grid, and this one is not "
                "(ERSSTv5 is 2-degree). Use climatology='self' to remove the field's own climatology."
            )
        anomaly = compute_anomaly_vs_own_climatology(
            data_region, freq=anomaly_freq, base_period=base_period, method="subtract"
        )
        baseline = anomaly.attrs.get("climatology_reference", "the field's own climatology")

    surface_note = "surface"
    if "depth" in anomaly.dims:
        # The top level by default, matching CMIP6's tos; a depth_level averages
        # the levels above it, weighted by the water each cell holds.
        anomaly = surface_layer(anomaly, depth_level)
        surface_note = anomaly.attrs.get("depth_averaging", surface_note)

    nino_index = compute_basic_stat(data=anomaly, op="area_weighted_mean", dim=["lat", "lon"])

    nino_index.attrs = {
        "long_name": f"{nino_region.upper()} Index",
        "units": "degC",
        "description": f"Area-weighted SST anomaly over the {nino_region.upper()} region",
        "region_bounds": f"lat: [{lat_min}, {lat_max}], lon: [{lon_min}, {lon_max}]",
        "depth_range": surface_note,
        "anomaly_frequency": anomaly_freq,
        "anomaly_method": "subtract",
        "weighting": "area weighted on the dataset's own grid",
        "climatology_baseline": baseline,
    }
    nino_index.name = f"{nino_region}_index"
    return nino_index


@register_tools(
    name="compute_pdo_index",
    description="Compute the Pacific Decadal Oscillation index: the leading EOF of North Pacific sea "
        "surface temperature anomalies.\n"
        "The global mean anomaly is removed first, so the index measures the North Pacific pattern rather "
        "than global warming, and cells are weighted by the square root of their area so the EOF is taken "
        "on a field whose variance is area-proportional.\n"
        "The sign of an EOF is arbitrary, so it is fixed by the convention that the positive phase is cool "
        "in the central and western North Pacific. Both the pattern and the index follow it, which makes a "
        "negative correlation between two patterns a real disagreement rather than an arbitrary flip.\n"
        "The EOF is taken over a rectangular latitude-longitude block, so a field on a curvilinear grid "
        "must be regridded first.",
    args={
        "data": "xr.DataArray, required. Temperature with 'time' and regular 'lat'/'lon' dimensions, "
            "optionally 'depth', in degrees Celsius.",
        "depth_level": "float or None, optional, default None. Only used when the field has a depth "
            "dimension. None takes the shallowest level, matching CMIP6's tos; a number averages every "
            "level shallower than it, thickness-weighted, which is a different quantity.",
        "standardize": "bool, optional, default True. Scale the index to zero mean and unit variance, "
            "which is how it is normally published."
    },
    returns="tuple (pdo_index, eof_pattern, explained_variance, regression_pattern):\n"
        "- pdo_index: the leading principal component as a monthly series.\n"
        "- eof_pattern: the unit-norm eigenvector. It says WHERE the mode acts, but it is "
        "dimensionless and its magnitude means nothing on its own; do not read an amplitude off it.\n"
        "- explained_variance: the percentage of variance the mode accounts for. It depends on the "
        "grid and on the length of the record, so compare it only between fields on the same grid "
        "over the same period.\n"
        "- regression_pattern: the anomaly field regressed on the STANDARDISED index, in the field's "
        "own units per standard deviation. This is the map with a physical amplitude, and the form a "
        "PDO pattern is normally drawn in; its area-weighted RMS is the mode's strength in degrees.\n"
        "Both patterns carry the grid they were computed on, so two of them can be compared with an "
        "area-weighted correlation — worth doing, since the domain reaches high latitudes where cells "
        "are much smaller. Two patterns on DIFFERENT grids cannot be correlated at all; put every "
        "field on one grid before comparing."
)
def compute_pdo_index(data: xr.DataArray, depth_level: float = None,
                      standardize: bool = True) -> tuple[xr.DataArray, xr.DataArray, float, xr.DataArray]:
    if "time" not in data.dims:
        raise ValueError("A PDO index needs a time series; this field has no 'time' dimension.")
    # The EOF is taken on a rectangular lat/lon block and the pattern comes back
    # as a map on those axes, neither of which a curvilinear model's index
    # dimensions provide. An area-weighted regional mean tolerates a ragged set
    # of cells; an EOF does not.
    if data.attrs.get("grid_type") == "curvilinear":
        raise ValueError(
            f"{data.attrs.get('model', 'This model')} is on a curvilinear grid, and the PDO is defined as "
            f"an EOF over a North Pacific lat/lon block, which its index dimensions cannot describe. Put "
            f"the field on a regular grid first with regrid_horizontal, then compute the index."
        )

    # A surface field (ERSSTv5 'sst', CMIP6 'tos') has no depth dimension at all,
    # so it is used as it stands rather than being averaged over one.
    if 'depth' in data.dims:
        sst_data = surface_layer(data, depth_level)
    else:
        sst_data = data

    # The PDO is defined on monthly ANOMALIES. Removing only the time mean, as an
    # EOF of raw SST, leaves the seasonal cycle as by far the largest mode: it
    # then comes out as the leading EOF with over 90 percent of the variance and
    # the result is a seasonal cycle, not the PDO.
    sst_anomaly = compute_anomaly_vs_own_climatology(sst_data, freq="month", method="subtract")

    # North Pacific, 20N-70N, 110E-100W: a box that crosses the dateline.
    north_pacific = _select_region(sst_anomaly, 20, 70, 110, -100)

    # The PDO is defined on SST anomalies with the global mean removed, so that
    # it measures the North Pacific pattern rather than the global warming signal.
    global_mean = compute_basic_stat(sst_anomaly, op='area_weighted_mean', dim=['lat', 'lon'])
    north_pacific = north_pacific - global_mean
    north_pacific.attrs = dict(sst_anomaly.attrs)

    # Cells shrink towards the pole, and this box reaches 70N, so each cell is
    # weighted by the square root of its area: the EOF is then computed on a
    # field whose variance is area-proportional.
    areas = areas_like(north_pacific)
    weights = np.sqrt(areas)
    weighted_anomalies = north_pacific * weights
    weighted_anomalies = weighted_anomalies - weighted_anomalies.mean(dim='time')

    n_times = weighted_anomalies.sizes['time']
    n_lat = weighted_anomalies.sizes['lat']
    n_lon = weighted_anomalies.sizes['lon']
    data_matrix = weighted_anomalies.transpose('time', 'lat', 'lon').values.reshape(n_times, n_lat * n_lon)

    # Land columns are dropped rather than zero-filled, so they cannot dilute the
    # leading mode; any remaining gap in an ocean column is filled with zero,
    # which is the mean of an anomaly field.
    valid_points = ~np.isnan(data_matrix).any(axis=0)
    if valid_points.sum() < 2:
        raise ValueError(
            "Fewer than two North Pacific cells have a complete time series, so no EOF can be computed. "
            "Check that the field covers 20N-70N, 110E-100W over the whole period."
        )
    data_matrix_clean = np.nan_to_num(data_matrix[:, valid_points], nan=0.0)

    pca = PCA(n_components=1)
    pdo_index_values = pca.fit_transform(data_matrix_clean).flatten()

    eof_pattern_1d = np.full(n_lat * n_lon, np.nan)
    eof_pattern_1d[valid_points] = pca.components_[0, :]
    eof_pattern_2d = eof_pattern_1d.reshape(n_lat, n_lon)

    # Undo the area weighting so the pattern is in the units of the field itself.
    weight_values = weights.transpose('lat', 'lon').values
    with np.errstate(invalid='ignore', divide='ignore'):
        eof_pattern_2d = eof_pattern_2d / weight_values

    explained_var = float(pca.explained_variance_ratio_[0] * 100)

    # The sign of an EOF is arbitrary, so it is fixed by the convention that the
    # PDO's positive phase is cool in the central/western North Pacific: the
    # pattern is flipped if the region west of the dateline is not negative.
    western = eof_pattern_2d[:, (north_pacific.lon.values >= 150) | (north_pacific.lon.values <= -170)]
    if np.nanmean(western) > 0:
        eof_pattern_2d = -eof_pattern_2d
        pdo_index_values = -pdo_index_values

    if standardize:
        pdo_index_values = (pdo_index_values - np.nanmean(pdo_index_values)) / np.nanstd(pdo_index_values)

    # The EOF vector above has unit norm, so it says where the mode acts but not
    # how strongly in degrees. Regressing the unweighted anomaly field on the
    # STANDARDISED index gives the amplitude the mode actually has, in the
    # field's own units per standard deviation of the index, which is how a PDO
    # pattern is normally drawn and the only form with a physical magnitude.
    standardised_index = (pdo_index_values - np.nanmean(pdo_index_values)) / np.nanstd(pdo_index_values)
    unweighted = north_pacific.transpose('time', 'lat', 'lon').values.reshape(
        n_times, n_lat * n_lon)
    centred = unweighted - np.nanmean(unweighted, axis=0, keepdims=True)
    with np.errstate(invalid='ignore'):
        regression_1d = np.nansum(centred * standardised_index[:, None], axis=0) / n_times
    regression_1d[~valid_points] = np.nan
    regression_2d = regression_1d.reshape(n_lat, n_lon)

    pdo_index = xr.DataArray(
        pdo_index_values,
        coords={'time': north_pacific.time},
        dims=['time'],
        attrs={
            'long_name': 'Pacific Decadal Oscillation Index',
            'description': 'Leading principal component of North Pacific SST anomalies (20N-70N)',
            'units': 'standardized' if standardize else 'dimensionless',
            'method': 'EOF of area-weighted monthly SST anomalies with the global mean anomaly removed',
            'region': 'North Pacific (20N-70N, 110E-100W)',
            'standardized': standardize,
            'sign_convention': 'positive phase is cool in the central and western North Pacific',
        }
    )

    eof_pattern = xr.DataArray(
        eof_pattern_2d,
        coords={'lat': north_pacific.lat, 'lon': north_pacific.lon},
        dims=['lat', 'lon'],
        attrs={
            'long_name': 'First EOF pattern of North Pacific SST',
            'description': 'Spatial pattern of the leading mode of SST variability',
            'units': 'dimensionless',
            'explained_variance_percent': explained_var,
            'sign_convention': 'positive phase is cool in the central and western North Pacific',
        }
    )

    regression_pattern = xr.DataArray(
        regression_2d,
        coords={'lat': north_pacific.lat, 'lon': north_pacific.lon},
        dims=['lat', 'lon'],
        attrs={
            'long_name': 'SST anomaly regressed on the standardized PDO index',
            'description': 'Amplitude of the leading mode, in the field units per standard '
                           'deviation of the index. This is the map with a physical magnitude; '
                           'eof_pattern is the unit-norm eigenvector and carries none.',
            'units': sst_data.attrs.get('units', 'degC'),
            'explained_variance_percent': explained_var,
            'sign_convention': 'positive phase is cool in the central and western North Pacific',
        }
    )

    # The patterns are maps on the grid the SST came from, so they keep that
    # grid's identity: without it the cell areas cannot be found, and comparing
    # two models' patterns is exactly what needs an AREA-WEIGHTED correlation —
    # the domain reaches 70N, where a cell covers a third of a tropical one.
    eof_pattern = carry_provenance(eof_pattern, north_pacific)
    regression_pattern = carry_provenance(regression_pattern, north_pacific)

    return pdo_index, eof_pattern, explained_var, regression_pattern
