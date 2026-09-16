# tools: regrid_horizontal, regrid_vertical, select_basin
#
# Preprocessing steps that put a field on a different grid before it can be
# compared with anything else. Nothing here changes what a field means: the
# values are the same field, carried onto another set of cells.

import numpy as np
import xarray as xr

from src.tools.registry import register_tools
from src.tools.data_openers import PROVENANCE_ATTRS, open_cmip6_grid_file
from src.tools.data_loader import bounds_from_centres, load_areas


# Building the ESMF weights for a curvilinear model takes seconds, and a
# multi-model comparison regrids the same pair of grids over and over. The
# weights depend only on the two grids and the method, so they are kept.
_REGRIDDER_CACHE = {}

# The observational grids a field can be put onto by name. 'woa23' is the
# 1-degree grid WOA23, GLODAP and SEANOE share, and is the usual meeting point
# for a model-observation comparison.
_NAMED_TARGETS = {
    "woa23": "woa23",
    "1deg": "woa23",
    "glodap": "glodap",
    "seanoe_mld": "seanoe_mld",
    "ersst": "ersst",
}

_METHODS = ("conservative", "bilinear", "nearest_s2d", "patch")


def _edges_from_centres(centres):
    """Cell edges as a 1-D array of length n+1, from n cell centres."""
    pairs = bounds_from_centres(np.asarray(centres, dtype="float64"))
    return np.concatenate([pairs[:, 0], pairs[-1:, 1]])


def _corners_from_vertices(vertex_lat, vertex_lon):
    """
    Assemble a curvilinear grid's (ny+1, nx+1) corner mesh from CMIP6's per-cell
    vertices, which are stored as (ny, nx, 4).

    A single shared corner mesh is only possible because neighbouring cells
    share the corners between them, but the models do not agree on the order
    their four vertices are written in: CESM2 starts at the south-west corner
    and winds anticlockwise, while MPI-ESM1-2-LR winds the other way on a grid
    whose j index runs southward. The winding is therefore worked out from the
    vertices themselves — which rotation makes every cell's corners coincide
    with its neighbours' — rather than assumed, because the wrong one produces a
    plausible-looking mesh describing cells that are not the model's.

    Returns (corner_lat, corner_lon, note). A handful of cells may still
    disagree where the grid closes on itself (the periodic seam and the fold
    over the Arctic); how many is reported in the note so it can be stated.
    """
    vertex_lat = np.asarray(vertex_lat, dtype="float64")
    vertex_lon = np.asarray(vertex_lon, dtype="float64")
    if vertex_lat.ndim != 3 or vertex_lat.shape[-1] != 4:
        raise ValueError(
            f"Expected per-cell corner vertices shaped (ny, nx, 4), got {vertex_lat.shape}."
        )

    def mismatch(lat_a, lon_a, lat_b, lon_b):
        """Distance between two vertices in degrees, insensitive to the 0/360 seam."""
        return np.maximum(np.abs(lat_a - lat_b), np.abs((lon_a - lon_b + 180.0) % 360.0 - 180.0))

    # The four vertices in index order: the corner shared with the previous cell
    # in both directions, then along i, then the far corner, then along j.
    best = None
    for first in range(4):
        for winding in (1, -1):
            order = tuple((first + step * winding) % 4 for step in range(4))
            low_low, low_high, high_high, high_low = order
            along_i = mismatch(vertex_lat[:, :-1, low_high], vertex_lon[:, :-1, low_high],
                               vertex_lat[:, 1:, low_low], vertex_lon[:, 1:, low_low])
            along_j = mismatch(vertex_lat[:-1, :, high_low], vertex_lon[:-1, :, high_low],
                               vertex_lat[1:, :, low_low], vertex_lon[1:, :, low_low])
            broken = int(np.sum(along_i > 1e-3) + np.sum(along_j > 1e-3))
            if best is None or broken < best[0]:
                best = (broken, order, along_i.size + along_j.size)

    broken, (low_low, low_high, high_high, high_low), total = best
    if broken > 0.01 * total:
        raise ValueError(
            f"This model's cell vertices do not form a single shared corner mesh: {broken} of {total} "
            f"neighbouring cell pairs disagree about the corner between them under every vertex ordering. "
            f"Conservative regridding intersects cells and needs that mesh; use method='bilinear' for this "
            f"field."
        )

    ny, nx = vertex_lat.shape[:2]
    corner_lat = np.empty((ny + 1, nx + 1), dtype="float64")
    corner_lon = np.empty((ny + 1, nx + 1), dtype="float64")
    for corners, vertices in ((corner_lat, vertex_lat), (corner_lon, vertex_lon)):
        corners[:-1, :-1] = vertices[:, :, low_low]
        corners[:-1, -1] = vertices[:, -1, low_high]
        corners[-1, :-1] = vertices[-1, :, high_low]
        corners[-1, -1] = vertices[-1, -1, high_high]

    note = (f"corner mesh assembled from per-cell vertices in index order {(low_low, low_high, high_high, high_low)}"
            + (f"; {broken} of {total} neighbouring cell pairs disagree about their shared corner, at the "
               f"grid's periodic seam and polar fold" if broken else ""))
    return corner_lat, corner_lon, note


def _open_grid_file(experiment, model, preferred_variable):
    """
    One of a model's files, opened for its grid geometry alone.

    The field being regridded is often a DERIVED one — a mixed layer depth, a
    nutricline, a trend map — whose name is no CMIP6 variable, so asking for a
    file of that name fails. Any variable of the same model is on the same
    horizontal grid, so the field's own name is tried first and anything the
    model publishes stands in for it.
    """
    candidates = [preferred_variable, "thetao", "so", "tos", "areacello", "deptho"]
    attempts = []
    seen = set()
    for variable in candidates:
        if not variable or variable in seen:
            continue
        seen.add(variable)
        try:
            return open_cmip6_grid_file(experiment, model, variable)
        except (FileNotFoundError, KeyError, ValueError, OSError) as exc:
            attempts.append(f"{variable}: {exc.__class__.__name__}")
    raise FileNotFoundError(
        f"No file of {model} in {experiment} could be opened for its grid geometry (tried {attempts}), "
        f"so the model's corner mesh is unavailable. Use method='bilinear', which needs only the field's "
        f"own 2-D coordinates."
    )


def _horizontal_dims(data):
    """The (y, x) dimension names of a field, curvilinear grids included."""
    y_dim = data.attrs.get("y_dim", "lat")
    x_dim = data.attrs.get("x_dim", "lon")
    if y_dim in data.dims and x_dim in data.dims:
        return y_dim, x_dim
    if "lat" in data.dims and "lon" in data.dims:
        return "lat", "lon"
    raise ValueError(
        f"This field has no horizontal dimensions to regrid. Dimensions: {list(data.dims)}."
    )


def _regular_grid(lat, lon, with_bounds):
    """An xESMF grid description for a rectilinear lat/lon grid."""
    grid = xr.Dataset({"lat": ("lat", np.asarray(lat, dtype="float64")),
                       "lon": ("lon", np.asarray(lon, dtype="float64"))})
    if with_bounds:
        grid["lat_b"] = ("lat_b", np.clip(_edges_from_centres(lat), -90.0, 90.0))
        grid["lon_b"] = ("lon_b", _edges_from_centres(lon))
    return grid


def _source_grid(data, method):
    """
    The grid `data` sits on, described the way xESMF wants it.

    A rectilinear field is described from its own coordinates, so a spatial
    subset regrids just as well as a global field. A curvilinear model is taken
    from its own grid file instead: conservative regridding needs the corner
    mesh, which only the file carries, and the centres must come from the same
    place and in the same longitude convention as the corners or the cells they
    describe are not the model's cells.
    """
    y_dim, x_dim = _horizontal_dims(data)
    needs_bounds = method == "conservative"

    if data.attrs.get("grid_type") != "curvilinear":
        grid = _regular_grid(data[y_dim].values, data[x_dim].values, needs_bounds)
        return grid, (y_dim, x_dim), "the field's own rectilinear cell centres"

    lat_coord = data.attrs.get("lat_coord", "lat")
    lon_coord = data.attrs.get("lon_coord", "lon")
    shape = (data.sizes[y_dim], data.sizes[x_dim])

    if not needs_bounds and lat_coord in data.coords and lon_coord in data.coords:
        grid = xr.Dataset({
            "lat": (("y", "x"), np.asarray(data[lat_coord].values, dtype="float64")),
            "lon": (("y", "x"), np.asarray(data[lon_coord].values, dtype="float64")),
        })
        return grid, (y_dim, x_dim), "the model's own 2-D cell centres"

    model = data.attrs.get("model")
    experiment = data.attrs.get("experiment")
    if not model or not experiment:
        raise ValueError(
            f"This curvilinear field carries neither its 2-D '{lat_coord}'/'{lon_coord}' coordinates nor "
            f"the 'model'/'experiment' attributes that would let them be read from the model's grid file, "
            f"so there is nothing to locate its cells with. It has {list(data.coords)}. Work from a field "
            f"that kept them, or regrid the field this one was derived from."
        )

    grid_file, info = _open_grid_file(experiment, model, data.name)
    try:
        model_cfg = info["model_cfg"]
        names = [("vertices_latitude", "vertices_longitude"), ("lat_bnds", "lon_bnds")]
        pair = next((p for p in names if p[0] in grid_file and grid_file[p[0]].ndim == 3), None)
        if pair is None and needs_bounds:
            raise ValueError(
                f"{model}'s grid file publishes no per-cell corner vertices, so its cells cannot be "
                f"intersected with the target grid. Use method='bilinear'."
            )
        centre_lat = np.asarray(grid_file[model_cfg["lat_coord"]].values, dtype="float64")
        centre_lon = np.asarray(grid_file[model_cfg["lon_coord"]].values, dtype="float64")
        if needs_bounds:
            corner_lat, corner_lon, note = _corners_from_vertices(grid_file[pair[0]].values,
                                                                  grid_file[pair[1]].values)
        else:
            corner_lat = corner_lon = None
            note = f"{model}'s cell centres, read from its grid file"
    finally:
        grid_file.close()

    if centre_lat.shape != shape:
        raise ValueError(
            f"This field covers {shape} cells but {model}'s grid is {centre_lat.shape}: the grid geometry "
            f"is read from the model's own file, which only describes the whole grid, and a field that has "
            f"lost its 2-D coordinates cannot say which part of it this is. Regrid the global field and "
            f"subset afterwards, or keep the field's lat/lon coordinates and use method='bilinear'."
        )

    variables = {
        "lat": (("y", "x"), centre_lat),
        "lon": (("y", "x"), centre_lon),
    }
    if needs_bounds:
        variables["lat_b"] = (("y_b", "x_b"), np.clip(corner_lat, -90.0, 90.0))
        variables["lon_b"] = (("y_b", "x_b"), corner_lon)
    return xr.Dataset(variables), (y_dim, x_dim), note


def _target_grid(target, method):
    """The grid to regrid onto, plus the attrs that identify it downstream."""
    needs_bounds = method == "conservative"

    if isinstance(target, (xr.DataArray, xr.Dataset)):
        template = target if isinstance(target, xr.DataArray) else target[list(target.data_vars)[0]]
        if "lat" not in template.dims or "lon" not in template.dims:
            raise ValueError(
                "A template must be a field on a rectilinear grid, with 'lat' and 'lon' dimensions. "
                "Regridding ONTO a curvilinear model grid is not supported; put both fields on a regular "
                "grid instead."
            )
        grid = _regular_grid(template["lat"].values, template["lon"].values, needs_bounds)
        identity = {k: v for k, v in template.attrs.items() if k in PROVENANCE_ATTRS}
        return grid, identity, f"the grid of the supplied template ({template.name or 'unnamed'})"

    if not isinstance(target, str) or target.lower() not in _NAMED_TARGETS:
        raise ValueError(
            f"target must be one of {sorted(_NAMED_TARGETS)} or a DataArray to use as a template, "
            f"got {target!r}."
        )
    dataset = _NAMED_TARGETS[target.lower()]
    areas = load_areas(dataset=dataset)
    grid = _regular_grid(areas["lat"].values, areas["lon"].values, needs_bounds)
    resolution = "1-degree" if dataset != "ersst" else "2-degree"
    return grid, {"source_dataset": dataset}, f"the {dataset} {resolution} grid"


def _get_regridder(source_grid, target_grid, method, periodic):
    """A cached xESMF regridder for this pair of grids."""
    def signature(grid):
        return tuple(
            (name, grid[name].shape, float(np.nansum(grid[name].values)))
            for name in sorted(grid.data_vars) + sorted(grid.coords)
        )

    key = (signature(source_grid), signature(target_grid), method, periodic)
    if key not in _REGRIDDER_CACHE:
        try:
            import xesmf
        except ImportError as exc:
            raise ImportError(
                "Regridding needs xESMF (with ESMPy). Install it with "
                "'conda install -c conda-forge xesmf'."
            ) from exc
        # A model grid closes on itself over the Arctic, and the cells along that
        # fold collapse to zero area. ESMF refuses to build conservative weights
        # when it meets one, which would fail the whole regridding over a handful
        # of cells; they are skipped instead and simply contribute nothing.
        _REGRIDDER_CACHE[key] = xesmf.Regridder(
            source_grid, target_grid, method, periodic=periodic, ignore_degenerate=True,
            unmapped_to_nan=True,
        )
    return _REGRIDDER_CACHE[key]


@register_tools(
    name="regrid_horizontal",
    description="Put a field on another horizontal grid so two products can be compared cell by cell.\n"
        "A difference map, a pattern correlation or an RMSE needs both fields on the same cells. Model "
        "grids are often curvilinear, with index dimensions that are not geographic, and observational "
        "products come on grids of their own; moving one onto the other is what makes them comparable.\n"
        "Regrid only for the cell-by-cell comparison. An average, a trend or an integral belongs on the "
        "field's own grid with its own cell areas, because regridding replaces those with the target's.\n"
        "The result carries the target grid's identity, so later weighting uses the target's cell areas.",
    args={
        "data": "xr.DataArray, required. The field to move, on a regular or a curvilinear grid. Other "
            "dimensions (time, depth) pass through untouched. It must still carry the attrs its loader set, "
            "which is how its own grid is identified.",
        "target": "str or xr.DataArray, optional, default 'woa23'. The grid to move onto: a registered "
            "observational grid by name ('woa23', '1deg', 'glodap', 'seanoe_mld', 'ersst'), or a DataArray "
            "on a regular grid to use as a template. The target must be rectilinear.",
        "method": "str, optional, default 'conservative'. 'conservative' area-averages the source cells "
            "that overlap each target cell, conserving the field's area integral; it is the right choice "
            "going to a coarser grid and before any area-weighted step. 'bilinear' interpolates between "
            "cell centres, and is the fallback when a model's corner mesh is unusable. 'nearest_s2d' "
            "copies the nearest source cell, for masks and other categorical fields. 'patch' is a "
            "higher-order smooth interpolation.",
        "skipna": "bool, optional, default True. Ignore missing cells instead of letting them spread. Land "
            "is missing in every ocean field, so with skipna=False every coastal target cell comes back "
            "empty.",
        "na_thres": "float, optional, default 0.5. How much of a target cell may be missing before the "
            "cell itself is set to NaN, from 0 (any missing source data empties it) to 1 (keep it if any "
            "valid data contributes)."
    },
    returns="xr.DataArray on the target grid, with ('lat', 'lon') in place of the source's horizontal "
        "dimensions and everything else unchanged. Missing data stays missing. attrs record the method and "
        "the source grid, and carry the TARGET grid's identity so later weighting uses the target's areas. "
        "State the method whenever a regridded RMSE or pattern correlation is quoted: the number depends "
        "on it."
)
def regrid_horizontal(data: xr.DataArray, target="woa23", method: str = "conservative",
                      skipna: bool = True, na_thres: float = 0.5) -> xr.DataArray:
    if not isinstance(data, xr.DataArray):
        raise TypeError(f"data must be an xr.DataArray, got {type(data).__name__}.")
    method = str(method).lower()
    if method not in _METHODS:
        raise ValueError(f"method must be one of {list(_METHODS)}, got '{method}'.")
    if not 0.0 <= na_thres <= 1.0:
        raise ValueError(f"na_thres must lie between 0 and 1, got {na_thres}.")

    source_grid, (y_dim, x_dim), source_note = _source_grid(data, method)
    target_grid, identity, target_name = _target_grid(target, method)

    # ESMF wraps the source grid round the globe only when told to, and a
    # longitude subset must not be wrapped: its two edges are not neighbours.
    # Conservative regridding intersects cells and needs no wrapping at all.
    span = float(np.nanmax(source_grid["lon"].values) - np.nanmin(source_grid["lon"].values))
    periodic = method != "conservative" and span > 350.0

    regridder = _get_regridder(source_grid, target_grid, method, periodic)

    # xESMF matches the field to the grid by the shape of its last two
    # dimensions, so the horizontal pair has to be the trailing one.
    ordered = data.transpose(*[d for d in data.dims if d not in (y_dim, x_dim)], y_dim, x_dim)
    regridded = regridder(ordered, keep_attrs=True, skipna=skipna, na_thres=na_thres)

    # The field now lives on the target's cells, so it must answer to the
    # target's grid identity: areas_like and the weighted means read these attrs
    # to decide which cell areas to use, and leaving the source's behind would
    # weight a 1-degree field with a curvilinear model's areacello.
    source_identity = ", ".join(
        f"{key}={data.attrs[key]}" for key in ("source_dataset", "model", "experiment", "grid_type")
        if key in data.attrs
    )
    for key in PROVENANCE_ATTRS:
        regridded.attrs.pop(key, None)
    regridded.attrs.update(identity)
    regridded.attrs["regridded_from"] = source_identity or "an array carrying no grid provenance"
    regridded.attrs["regrid_method"] = method
    regridded.attrs["regrid_target"] = target_name
    regridded.attrs["regrid_source_grid"] = source_note
    regridded.attrs["regrid_na_thres"] = na_thres if skipna else "skipna=False (any missing source cell gives NaN)"
    regridded.name = data.name
    return regridded


# Which WOA23 basin codes make up each ocean. The three great oceans are codes
# 1-3, but the mask gives every marginal sea and deep basin a code of its own,
# so an "Atlantic" written as code 1 alone quietly drops the Gulf of Mexico, the
# Caribbean, the Mediterranean and the whole north-west Atlantic shelf. The
# groupings below fold each marginal sea into the ocean it opens onto; they are
# reported in the result's attrs so a basin mean can be reproduced exactly.
# The Southern and Arctic Oceans stay separate: they are oceans in their own
# right, and folding them into their neighbours would double-count them.
BASIN_GROUPS = {
    "Atlantic": {
        "main": [1],
        "marginal": [4, 5, 6, 9, 15, 16, 17, 20, 22, 23, 30, 31, 36, 37, 48, 49, 50, 55],
    },
    "Pacific": {
        "main": [2],
        "marginal": [12, 14, 18, 19, 25, 27, 28, 29, 33, 38, 39, 40, 41, 42, 44, 51, 54, 57, 59],
    },
    "Indian": {
        "main": [3],
        "marginal": [7, 8, 21, 32, 34, 35, 43, 45, 46, 47, 56],
    },
    "Southern": {
        "main": [10],
        "marginal": [52, 58, 60],
    },
    "Arctic": {
        "main": [11],
        "marginal": [13, 24, 26],
    },
}

# The Caspian Sea (53) belongs to no ocean and is never folded into one.


def _basin_codes(basin, include_marginal_seas):
    """The WOA23 codes a basin request stands for, and a description of them."""
    from src.tools.data_loader import BASIN_CODES, BASIN_NAMES

    if isinstance(basin, (int, np.integer)):
        basin = [int(basin)]
    if isinstance(basin, (list, tuple, set)) and all(
        isinstance(code, (int, np.integer)) for code in basin
    ):
        codes = sorted(int(code) for code in basin)
        unknown = [code for code in codes if code not in BASIN_NAMES]
        if unknown:
            raise ValueError(f"Unknown WOA23 basin code(s) {unknown}. Valid codes are 1-60.")
        return codes, "explicit codes: " + ", ".join(f"{c} {BASIN_NAMES[c]}" for c in codes)

    if not isinstance(basin, str):
        raise TypeError(f"basin must be a name, a WOA23 code or a list of codes, got {basin!r}.")

    key = basin.strip().replace(" Ocean", "").title()
    if key in BASIN_GROUPS:
        group = BASIN_GROUPS[key]
        codes = sorted(group["main"] + (group["marginal"] if include_marginal_seas else []))
        folded = ", ".join(BASIN_NAMES[c] for c in group["marginal"])
        description = (
            f"{key} Ocean, WOA23 codes {codes}"
            + (f", with its marginal seas and sub-basins folded in ({folded})"
               if include_marginal_seas else ", main basin only; marginal seas excluded")
        )
        return codes, description

    # Any other name from the mask's own table, e.g. 'Mediterranean Sea'.
    match = next((name for name in BASIN_CODES if name.lower() == basin.strip().lower()), None)
    if match is None:
        raise ValueError(
            f"'{basin}' is not a basin this mask knows. Use one of {sorted(BASIN_GROUPS)} for a whole "
            f"ocean, or one of the WOA23 basin names: {sorted(BASIN_CODES)}."
        )
    return [BASIN_CODES[match]], f"{match}, WOA23 code {BASIN_CODES[match]}"


@register_tools(
    name="select_basin",
    description="Keep the part of a field that lies in one ocean basin and blank the rest.\n"
        "This is what a per-basin number needs — an Atlantic mean, a Pacific bias, an Indian Ocean RMSE. "
        "It masks a field you already have, so a difference, a trend map or any other computed quantity "
        "can be split by basin, which a loader that cuts a region while reading cannot do.\n"
        "Marginal seas each carry their own code in the mask, so an ocean taken as its main code alone "
        "silently excludes a large part of what the name usually means. They are folded into the ocean "
        "they open onto by default; either choice is defensible as long as it is stated, and the codes "
        "used are recorded in the result.\n"
        "The mask is defined on a 1-degree grid but is carried onto whatever grid the field is on by "
        "nearest cell, so a curvilinear model can be split on its own grid and weighted with its own cell "
        "areas, with no regridding.",
    args={
        "data": "xr.DataArray, required. Any field with horizontal dimensions — loaded, derived or "
            "differenced. It must still carry the attrs its loader set, which is how its grid is "
            "identified.",
        "basin": "str or int or list[int], required. 'Atlantic', 'Pacific', 'Indian', 'Southern' or "
            "'Arctic' for a whole ocean; the name of a single sea in the mask, such as 'Mediterranean Sea', "
            "for that basin alone; or raw mask codes for full control.",
        "include_marginal_seas": "bool, optional, default True. Fold each ocean's marginal seas and "
            "sub-basins into it. Ignored when a single named sea or explicit codes are requested. The "
            "Southern and Arctic Oceans are always separate basins, never folded into their neighbours.",
        "depth_dependent": "bool, optional, default False. The mask labels every standard level "
            "separately, because a marginal sea stops at its sill and below that depth the water belongs "
            "to the open ocean. False applies the surface labels through the whole water column, which is "
            "what a named ocean usually means. True applies each level's own labels, which is what a "
            "water-mass analysis wants; the field must then have a 'depth' dimension."
    },
    returns="xr.DataArray with the same dimensions, coordinates and grid as the input, and everything "
        "outside the basin set to NaN. Because the grid is untouched, a later area-weighted mean uses the "
        "same cell areas as the global field and masks them to the basin automatically. attrs record "
        "'basin', 'basin_codes' and 'basin_definition' — state the definition with any basin number, "
        "because it is part of the answer."
)
def select_basin(data: xr.DataArray, basin, include_marginal_seas: bool = True,
                 depth_dependent: bool = False) -> xr.DataArray:
    from src.tools.data_loader import load_basin_mask, mask_1deg_to_grid

    if not isinstance(data, xr.DataArray):
        raise TypeError(f"data must be an xr.DataArray, got {type(data).__name__}.")
    codes, definition = _basin_codes(basin, include_marginal_seas)

    grid_type = data.attrs.get("grid_type", "regular")
    lat_coord = data.attrs.get("lat_coord", "lat")
    lon_coord = data.attrs.get("lon_coord", "lon")
    for name in (lat_coord, lon_coord):
        if name not in data.coords:
            raise ValueError(
                f"This field has no '{name}' coordinate, so it cannot be located against the basin mask. "
                f"Load it with the data tools, or work from a field that kept its coordinates."
            )

    basin_mask = load_basin_mask()
    if depth_dependent:
        if "depth" not in data.dims:
            raise ValueError(
                "depth_dependent=True applies each level's own basin labels, so the field must have a "
                "'depth' dimension. Use the default for a surface or vertically collapsed field."
            )
        # Each data level takes the labels of the nearest WOA standard level.
        codes_on_grid = basin_mask.sel(depth=data["depth"], method="nearest")
        codes_on_grid = codes_on_grid.assign_coords(depth=data["depth"])
        level_note = "each level's own basin labels (nearest WOA23 standard level)"
    else:
        codes_on_grid = basin_mask.isel(depth=0, drop=True)
        level_note = "the surface basin labels, applied through the whole water column"

    in_basin = codes_on_grid.isin(codes)

    if "depth" in in_basin.dims:
        # mask_1deg_to_grid works one horizontal field at a time, so the levels
        # are carried across one by one and stacked back together.
        levels = [
            mask_1deg_to_grid(in_basin.isel(depth=i), data, grid_type, lat_coord, lon_coord)
            for i in range(in_basin.sizes["depth"])
        ]
        mask = xr.concat(levels, dim=in_basin["depth"])
    else:
        mask = mask_1deg_to_grid(in_basin, data, grid_type, lat_coord, lon_coord)

    if not bool(mask.any()):
        raise ValueError(
            f"The basin '{basin}' ({definition}) covers no cell of this field. Check that the field is "
            f"global, or that it has not already been masked to a region elsewhere."
        )

    selected = data.where(mask)
    selected.attrs = dict(data.attrs)
    selected.attrs["basin"] = basin if isinstance(basin, str) else str(codes)
    selected.attrs["basin_codes"] = str(codes)
    selected.attrs["basin_definition"] = f"WOA23 basinmask_01: {definition}; {level_note}"
    selected.name = data.name
    return selected


# The vertical grids a field can be put onto by name.
_NAMED_DEPTH_TARGETS = ("glodap", "woa23")


@register_tools(
    name="regrid_vertical",
    description="Interpolate a field onto another set of depth levels so profiles can be compared level by "
        "level.\n"
        "Vertical grids differ between products, and rarely line up. A bias asked for 'at 100, 500 and "
        "1000 m' cannot be answered by selecting levels, because each product would answer at a different "
        "depth, and those differ most where the gradients are steepest.\n"
        "This is for comparison, not for budgets: an integral or a volume-weighted mean belongs on the "
        "field's own levels with its own layer thicknesses, because interpolation does not conserve the "
        "vertical integral.",
    args={
        "data": "xr.DataArray, required. A field with a 'depth' dimension in metres, positive downward. "
            "Other dimensions pass through untouched.",
        "target": "str, list[float] or xr.DataArray, optional, default 'glodap'. Where to interpolate to: "
            "a registered vertical grid by name ('glodap', 'woa23'), an explicit list of depths in metres, "
            "or another field whose 'depth' coordinate is to be matched.",
        "method": "str, optional, default 'linear'. 'linear' interpolates between the two bracketing "
            "levels, which suits a smooth profile. 'nearest' takes the closer level without averaging, for "
            "a mask or a categorical field.",
        "fill_surface": "bool, optional, default True. Give target depths above the field's shallowest "
            "level that level's value, since a model's first level is its surface. Nothing is ever "
            "extrapolated downwards: below a column's deepest valid level the result is missing."
    },
    returns="xr.DataArray on the target levels, with every other dimension, coordinate and the grid "
        "identity unchanged. Missing below each column's deepest valid level, and wherever a gap in the "
        "profile straddles a target depth. attrs record the method and the axis the field was put onto."
)
def regrid_vertical(data: xr.DataArray, target="glodap", method: str = "linear",
                    fill_surface: bool = True) -> xr.DataArray:
    if not isinstance(data, xr.DataArray):
        raise TypeError(f"data must be an xr.DataArray, got {type(data).__name__}.")
    if "depth" not in data.dims:
        raise ValueError(
            f"This field has no 'depth' dimension to interpolate along; it has {list(data.dims)}. A "
            f"surface field or one already averaged over depth has no profile left to put on another axis."
        )
    if method not in ("linear", "nearest"):
        raise ValueError(f"method must be 'linear' or 'nearest', got '{method}'.")

    if isinstance(target, str):
        if target.lower() not in _NAMED_DEPTH_TARGETS:
            raise ValueError(
                f"target must be one of {list(_NAMED_DEPTH_TARGETS)}, a list of depths in metres, or a "
                f"DataArray to take the levels from, got '{target}'."
            )
        from src.tools.data_loader import load_depth_levels
        levels = np.asarray(load_depth_levels(dataset=target.lower())[0].values, dtype="float64")
        target_name = f"the {target.lower()} standard levels"
    elif isinstance(target, (xr.DataArray, xr.Dataset)):
        if "depth" not in target.coords:
            raise ValueError("The template has no 'depth' coordinate to take levels from.")
        levels = np.asarray(target["depth"].values, dtype="float64")
        target_name = "the levels of the supplied template"
    else:
        levels = np.atleast_1d(np.asarray(target, dtype="float64"))
        target_name = f"{levels.size} requested level(s)"

    if levels.size == 0:
        raise ValueError("No target depths were given.")
    if np.any(levels < 0):
        raise ValueError(f"Depths are positive downward; {levels[levels < 0].tolist()} are negative.")

    # interp needs the whole profile in one chunk, and a lazily-read model field
    # is chunked along time, not depth, so only this axis is gathered.
    if data.chunks is not None and len(data.chunksizes.get("depth", ())) > 1:
        data = data.chunk({"depth": -1})

    interpolated = data.interp(depth=levels, method=method, kwargs={"fill_value": np.nan})

    shallowest = float(np.asarray(data["depth"].values).min())
    surface_note = "no extrapolation"
    if fill_surface and np.any(levels < shallowest):
        first_level = data.sel(depth=shallowest, drop=True)
        interpolated = xr.where(interpolated["depth"] < shallowest, first_level, interpolated)
        filled = levels[levels < shallowest]
        surface_note = (f"target depth(s) {filled.tolist()} above the field's shallowest level "
                        f"({shallowest:g} m) took that level's value")

    interpolated.attrs = dict(data.attrs)
    interpolated.attrs["vertical_regrid"] = (
        f"{method} interpolation onto {target_name} from {data.sizes['depth']} native levels "
        f"({shallowest:g} to {float(np.asarray(data['depth'].values).max()):g} m); {surface_note}. "
        f"Interpolated for comparison: it does not conserve a vertical integral."
    )
    interpolated["depth"].attrs.update({"units": "m", "positive": "down"})
    interpolated.name = data.name
    return interpolated
