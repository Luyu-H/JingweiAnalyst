# tools: load_regular_data, load_data_with_geolocator, load_areas, load_depth_levels

import re
import os
import yaml
import numpy as np
import xarray as xr
import pandas as pd
import xesmf as xe
from pyproj import Geod
from typing import Optional, Union, List

from src.tools.registry import register_tools
from src.utils.geolocator import Geolocator


@register_tools(
    name="load_regular_data",
    description="Load a regular subset of ocean data from a Zarr dataset with dimensions (time, depth, lat, lon). "
        "This tool supports only continuous range-based subsetting along 1–4 dimensions, "
        "rather than irregular or point-based selections. "
        "The target variable name is automatically inferred from the file name.",
    args={
        "zarr_path": "str, required. Path to the zarr dataset. File name should follow the format 'recon_{var}_sjtu_jingwei_{ver}.zarr'.",
        "time_range": "list[str, str], optional. Time range (start, end) in unified date format (e.g., '2000-01-01').",
        "depth_range": "list[float, float], optional. Depth range (min_depth, max_depth) in meters.",
        "lat_range": "list[float, float], optional. Latitude range (min_lat, max_lat) in degrees.",
        "lon_range": "list[float, float], optional. Longitude range (min_lon, max_lon) in degrees."
    },
    returns="xr.DataArray. Regular subset of the requested variable data."
)
def load_regular_data(zarr_path: str, time_range: list[str, str] = None, depth_range: list[float, float] = None, 
                      lat_range: list[float, float] = None, lon_range: list[float, float] = None) -> xr.DataArray:
    fname = os.path.basename(zarr_path)
    match = re.match(r"recon_([a-zA-Z0-9]+)_", fname)
    if not match:
        raise ValueError(f"Cannot parse variable name from file: {fname}")
    var_name = match.group(1)

    ds = xr.open_zarr(zarr_path)
    if var_name not in ds:
        if "n" in ds and var_name == "no3":
            var_name = "n"
        elif "p" in ds and var_name == "po4":
            var_name = "p"
        elif "salinity" in ds and var_name == "so":
            var_name = "salinity"
        elif "temperature" in ds and var_name == "temp":
            var_name = "temperature"
        else:
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
        lon_min, lon_max = lon_range[0], lon_range[1]
        data_lon_min = float(data.lon.min())
        data_lon_max = float(data.lon.max())
        if lon_min > 180 or lon_max > 180:
            lon_min_converted = lon_min if lon_min <= 180 else lon_min - 360
            lon_max_converted = lon_max if lon_max <= 180 else lon_max - 360
            
            if lon_min_converted <= lon_max_converted:
                data = data.sel(lon=slice(lon_min_converted, lon_max_converted))
            else:
                part1 = data.sel(lon=slice(lon_min_converted, data_lon_max))
                part2 = data.sel(lon=slice(data_lon_min, lon_max_converted))
                
                data = xr.concat([part1, part2], dim='lon')
                
        else:
            data = data.sel(lon=slice(lon_min, lon_max))

    return data


@register_tools(
    name="load_data_with_geolocator",
    description="Load regional ocean data from a Zarr file using geographic area selection powered by Geolocator. ",
    args={
        "zarr_path": "str, required. Path to the zarr dataset. File name should follow the format 'recon_{var}_sjtu_jingwei_{ver}.zarr'.",
        "time_range": "tuple[str, str], optional. Time range (start, end) in unified date format (e.g., '2000-01-01').",
        "depth_range": "tuple[float, float], optional. Depth range (min_depth, max_depth) in meters.",
        "geofeature": "str, required. Must be one of 'country', 'state', 'continent', 'water_body'.",
        "loc_name": "str or list[str], required. Name(s) of the geofeature, e.g. 'Arctic Ocean' or ['North Atlantic Ocean', 'South Atlantic Ocean']."
    },
    returns="xr.DataArray. Regular subset of the requested variable data."
)
def load_data_with_geolocator(zarr_path: str, time_range: tuple[str, str] = None, depth_range: tuple[float, float] = None, 
                                geofeature: str = "water_body", loc_name: str | list[str] = None) -> xr.DataArray:
    fname = os.path.basename(zarr_path)
    match = re.match(r"recon_([a-zA-Z0-9]+)_", fname)
    if not match:
        raise ValueError(f"Cannot parse variable name from file: {fname}")
    var_name = match.group(1)

    ds = xr.open_zarr(zarr_path)
    if var_name not in ds:
        if "n" in ds and var_name == "no3":
            var_name = "n"
        elif "p" in ds and var_name == "po4":
            var_name = "p"
        elif "salinity" in ds and var_name == "so":
            var_name = "salinity"
        elif "temperature" in ds and var_name == "temp":
            var_name = "temperature"
        else:
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

    geolocator_lats = np.arange(-89.5, 90.5, 1.0)
    geolocator_lons = np.arange(-179.5, 180.5, 1.0)
    geolocator = Geolocator(geolocator_lats, geolocator_lons, cache_dir, natural_earth_path)

    if time_range is not None:
        data = data.sel(time=slice(time_range[0], time_range[1]))
    if depth_range is not None:
        data = data.sel(depth=slice(depth_range[0], depth_range[1]))

    if geofeature is not None and loc_name is not None:
        loc_names = [loc_name] if isinstance(loc_name, str) else loc_name

        combined_mask = None
        for name in loc_names:
            try:
                mask = geolocator.get_mask_by_geofeature(geofeature, name)
            except Exception as e:
                raise RuntimeError(f"Failed to get mask for {geofeature}='{name}': {e}")
            
            if combined_mask is None:
                combined_mask = mask
            else:
                combined_mask = combined_mask | mask
        
        combined_mask = combined_mask.T
        geolocator_lons_converted = np.where(geolocator_lons >= 180, geolocator_lons - 360, geolocator_lons)
        mask_da = xr.DataArray(combined_mask.astype(float), coords={"lat": geolocator_lats, "lon": geolocator_lons_converted}, dims=["lat", "lon"])
        mask_da = mask_da.interp(lat=data.lat, lon=data.lon, method="nearest", kwargs={"fill_value": 0})
        mask_da = mask_da > 0.5
        data = data.where(mask_da)
    
    return data


@register_tools(
    name="load_areas",
    description="Load grid cell area data for latitude-longitude grids. Calculates precise grid cell areas on the WGS84 ellipsoid, equivalent to MATLAB's areaquad function.",
    args={
        "lat_range": "list[float, float], optional. Latitude range [min_lat, max_lat] in degrees. Default is global range [-90, 90].",
        "lon_range": "list[float, float], optional. Longitude range [min_lon, max_lon] in degrees. Default is global range [-180, 180]."
    },
    returns="xr.DataArray. Grid cell areas in square meters (m²) with coordinates 'lat' and 'lon' (both float32, 1-degree resolution). Calculated using pyproj.Geod with WGS84 ellipsoid (a=6378137.0 m, f=1/298.257223563) for high-precision geodetic area computation."
)
def load_areas(lat_range: list[float, float] = None, lon_range: list[float, float] = None) -> xr.DataArray:
    if lat_range is None: lat_range = [-90, 90]
    if lon_range is None: lon_range = [-180, 180]
    full_lat = np.arange(-89.5, 90, 1.0, dtype=np.float32)
    full_lon = np.arange(-179.5, 180, 1.0, dtype=np.float32)
    lat_mask = (full_lat >= lat_range[0]) & (full_lat <= lat_range[1])
    lon_mask = (full_lon >= lon_range[0]) & (full_lon <= lon_range[1])
    lat = full_lat[lat_mask]
    lon = full_lon[lon_mask]

    lat_b = np.arange(-90, 90.5, 1.0, dtype=np.float64)
    lon_b = np.arange(-180, 180.5, 1.0, dtype=np.float64)
    lat_b_mask = (lat_b >= lat_range[0] - 0.5) & (lat_b <= lat_range[1] + 0.5)
    lon_b_mask = (lon_b >= lon_range[0] - 0.5) & (lon_b <= lon_range[1] + 0.5)
    lat_b_subset = lat_b[lat_b_mask]
    lon_b_subset = lon_b[lon_b_mask]

    geod = Geod(ellps='WGS84')
    n_lat = len(lat)
    n_lon = len(lon)
    areas = np.zeros((n_lat, n_lon))
    
    for i in range(n_lat):
        lat1 = lat_b_subset[i]
        lat2 = lat_b_subset[i + 1]
        
        lat1 = np.clip(lat1, -89.99, 89.99)
        lat2 = np.clip(lat2, -89.99, 89.99)
        
        for j in range(n_lon):
            lon1 = lon_b_subset[j]
            lon2 = lon_b_subset[j + 1]
        
            lons = np.array([lon1, lon2, lon2, lon1])
            lats = np.array([lat1, lat1, lat2, lat2])
            
            area, _ = geod.polygon_area_perimeter(lons, lats)
            areas[i, j] = abs(area)

    areas_da = xr.DataArray(
        areas,
        coords={
            'lat': ('lat', lat, {'units': 'degrees_north', 'long_name': 'latitude'}),
            'lon': ('lon', lon, {'units': 'degrees_east', 'long_name': 'longitude'})
        },
        dims=['lat', 'lon'],
        attrs={
            'units': 'm²',
            'long_name': 'grid cell area',
            'description': 'Area of each grid cell on WGS84 ellipsoid',
            'grid_resolution': '1 degree',
            'method': 'WGS84 ellipsoid (equivalent to MATLAB areaquad)',
            'ellipsoid': 'WGS84 (a=6378137.0 m, f=1/298.257223563)'
        }
    )
    
    return areas_da


@register_tools(
    name="load_depth_levels",
    description="Load depth level and boundary information for reconstruction data.",
    args={"depth_range": "list[float, float], optional. Depth range [min_depth, max_depth] in meters. Default is None (all depth levels)."},
    returns="tuple[xr.DataArray, xr.DataArray]. Returns (depth, depth_bnds) where depth is 1D array of depth levels (m) and depth_bnds is 2D array of depth boundaries (m) with shape (depth, 2)."
)
def load_depth_levels(depth_range: list[float, float] = None) -> tuple[xr.DataArray, xr.DataArray]:
    # Depth levels and boundaries of reconstruction data are the same as WOA23 climatology data
    source_file = "/data1/lorraine-20419/jingwei_QA/woa23_climatology/nitrate_all_time_periods/woa23_all_n00_01.nc"
    
    with xr.open_dataset(source_file, decode_times=False) as ds:
        depth_full = ds['depth'].astype(np.float32)
        depth_bnds_full = ds['depth_bnds'].astype(np.float32)
        
        if depth_range is not None:
            depth_mask = (depth_full.values >= depth_range[0]) & (depth_full.values <= depth_range[1])
            depth_values = depth_full.values[depth_mask]
            depth_bnds_values = depth_bnds_full.values[depth_mask, :]
        else:
            depth_values = depth_full.values
            depth_bnds_values = depth_bnds_full.values
        
        depth_da = xr.DataArray(
            depth_values,
            coords={
                'depth': ('depth', depth_values, {
                    'standard_name': 'depth',
                    'bounds': 'depth_bnds',
                    'positive': 'down',
                    'units': 'meters',
                    'axis': 'Z',
                    'long_name': 'depth'
                })
            },
            dims=['depth'],
            attrs={
                'standard_name': 'depth',
                'bounds': 'depth_bnds',
                'positive': 'down',
                'units': 'meters',
                'axis': 'Z',
                'long_name': 'depth',
                'source': 'WOA23 climatology'
            }
        )
        
        depth_bnds_da = xr.DataArray(
            depth_bnds_values,
            coords={
                'depth': ('depth', depth_values),
                'nbounds': ('nbounds', [0, 1])
            },
            dims=['depth', 'nbounds'],
            attrs={
                'comment': 'depth bounds',
                'units': 'meters',
                'long_name': 'depth bounds',
                'description': 'Upper and lower boundaries of each depth level',
                'source': 'WOA23 climatology'
            }
        )
    
    return depth_da, depth_bnds_da


# https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DOCUMENTATION/WOA23_Product_Documentation.pdf#page=15&zoom=100,93,86
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


def load_basin_mask(depth_levels: Optional[Union[str, List[str]]] = None) -> xr.DataArray:
    filepath = "/data1/lorraine-20419/jingwei_QA/others/basin_01.csv"
    df = pd.read_csv(filepath, comment='#')
    
    lats = df['Latitude'].unique()
    lons = df['Longitude'].unique()
    lats = np.sort(lats)
    lons = np.sort(lons)
    
    depth_cols = [col for col in df.columns if col.startswith('Basin_')]
    depths = []
    for col in depth_cols:
        depth_str = col.replace('Basin_', '').replace('m', '')
        depths.append(int(depth_str))
    depths = np.array(depths)
    
    basin_data = np.full((len(depths), len(lats), len(lons)), np.nan)
    
    # Fill data
    for _, row in df.iterrows():
        lat_idx = np.searchsorted(lats, row['Latitude'])
        lon_idx = np.searchsorted(lons, row['Longitude'])
        
        for d_idx, col in enumerate(depth_cols):
            value = row[col]
            if pd.notna(value) and value != '':
                try:
                    basin_data[d_idx, lat_idx, lon_idx] = float(value)
                except (ValueError, TypeError):
                    pass
    
    # Build full 3D mask first
    full_mask = xr.DataArray(
        data=basin_data,
        dims=['depth', 'lat', 'lon'],
        coords={
            'depth': depths,
            'lat': lats,
            'lon': lons
        },
        name='basin_mask',
        attrs={
            'long_name': 'WOA23 Basin Mask',
            'description': 'Ocean basin identification codes',
            'basin_names': BASIN_NAMES,
            'depth_levels': [f"{d}m" for d in depths],
            'valid_range': [1, 60]
        }
    )
    
    if depth_levels is not None:
        if isinstance(depth_levels, str):
            depth_levels = [depth_levels]
        
        selected_depths = []
        for dl in depth_levels:
            depth_val = int(dl.replace('m', ''))
            if depth_val not in depths:
                available = [f"{d}m" for d in depths]
                raise ValueError(f"Depth {dl} not found. Available: {available}")
            selected_depths.append(depth_val)
        
        if len(selected_depths) == 1:
            # Single depth, return 2D
            return full_mask.sel(depth=selected_depths[0])
        else:
            # Multiple depths, return 3D subset
            return full_mask.sel(depth=selected_depths)
    
    return full_mask