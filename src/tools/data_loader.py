# tools: load_regular_data, load_data_with_geolocator

import re
import os
import yaml
import numpy as np
import xarray as xr

from src.tools.registry import register_tools
from src.utils.geolocator import Geolocator


@register_tools(
    name="load_regular_data",
    description="Load regional ocean data from a selected zarr file with dimensions (time, depth, lat, lon). " \
    "The variable name is parsed automatically from the file name. Supports subsetting by 1 – 4 dimensions.",
    args={
        "zarr_path": "str, required. Path to the zarr dataset. File name should follow the format 'recon_{var}_sjtu_jingwei_{ver}.zarr'.",
        "time_range": "tuple[str, str], optional. Time range (start, end) in unified date format (e.g., '2000-01-01').",
        "depth_range": "tuple[float, float], optional. Depth range (min_depth, max_depth) in meters.",
        "lat_range": "tuple[float, float], optional. Latitude range (min_lat, max_lat) in degrees.",
        "lon_range": "tuple[float, float], optional. Longitude range (min_lon, max_lon) in degrees."
    },
    returns="xr.DataArray. Regular subset of the requested variable data."
)
def load_regular_data(zarr_path: str, time_range: tuple[str, str] = None, depth_range: tuple[float, float] = None, 
                      lat_range: tuple[float, float] = None, lon_range: tuple[float, float] = None) -> xr.DataArray:
    fname = os.path.basename(zarr_path)
    match = re.match(r"recon_([a-zA-Z0-9]+)_", fname)
    if not match:
        raise ValueError(f"Cannot parse variable name from file: {fname}")
    var_name = match.group(1)

    ds = xr.open_zarr(zarr_path)
    if var_name not in ds:
        raise KeyError(f"Variable '{var_name}' not found in {fname}. "
                       f"Available: {list(ds.data_vars)}")
    data = ds[var_name]

    if time_range is not None:
        data = data.sel(time=slice(time_range[0], time_range[1]))
    if depth_range is not None:
        data = data.sel(depth=slice(depth_range[0], depth_range[1]))
    if lat_range is not None:
        data = data.sel(lat=slice(lat_range[0], lat_range[1]))
    if lon_range is not None:
        data = data.sel(lon=slice(lon_range[0], lon_range[1]))

    return data


@register_tools(
    name="load_data_with_geolocator",
    description="Load regional ocean data from a Zarr file using geographic area selection powered by Geolocator. ",
    args={
        "zarr_path": "str, required. Path to the zarr dataset. File name should follow the format 'recon_{var}_sjtu_jingwei_{ver}.zarr'.",
        "time_range": "tuple[str, str], optional. Time range (start, end) in unified date format (e.g., '2000-01-01').",
        "depth_range": "tuple[float, float], optional. Depth range (min_depth, max_depth) in meters.",
        "geofeature": "str, required. Must be one of 'country', 'state', 'continent', 'water_body'.",
        "loc_name": "str, required. Name of the geofeature, e.g. 'Arctic Ocean'."
    },
    returns="xr.DataArray. Regular subset of the requested variable data."
)
def load_data_with_geolocator(zarr_path: str, time_range: tuple[str, str] = None, depth_range: tuple[float, float] = None, 
                                geofeature: str = "water_body", loc_name: str = None) -> xr.DataArray:
    fname = os.path.basename(zarr_path)
    match = re.match(r"recon_([a-zA-Z0-9]+)_", fname)
    if not match:
        raise ValueError(f"Cannot parse variable name from file: {fname}")
    var_name = match.group(1)

    ds = xr.open_zarr(zarr_path)
    if var_name not in ds:
        raise KeyError(f"Variable '{var_name}' not found in {fname}. "
                       f"Available: {list(ds.data_vars)}")
    data = ds[var_name]

    config_path = "/home/lorraine/jingwei_QA/configs/default_paths.yaml"
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)
    
    cache_dir = cfg.get("cache_dir")
    natural_earth_path = cfg.get("natural_earth_path")
    if not cache_dir or not natural_earth_path:
        raise ValueError("Both 'cache_dir' and 'natural_earth_path' must be defined in configs/default_paths.yaml")

    latitudes = ds.lat.data
    longitudes = ds.lon.data
    geolocator = Geolocator(latitudes, longitudes, cache_dir, natural_earth_path)

    if time_range is not None:
        data = data.sel(time=slice(time_range[0], time_range[1]))
    if depth_range is not None:
        data = data.sel(depth=slice(depth_range[0], depth_range[1]))

    if geofeature is not None and loc_name is not None:
        try:
            mask = geolocator.get_mask_by_geofeature(geofeature, loc_name)
        except Exception as e:
            raise RuntimeError(f"Failed to get mask for {geofeature}='{loc_name}': {e}")
        
        mask_da = xr.DataArray(mask, coords={"lat": latitudes, "lon": longitudes}, dims=["lat", "lon"])
        data = data.where(mask_da)
    
    return data

