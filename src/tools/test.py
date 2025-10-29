import re
import os
import yaml
import numpy as np
import xarray as xr
import xesmf as xe
from pyproj import Geod

def load_regular_data(zarr_path: str, time_range: list[str, str] = None, depth_range: list[float, float] = None, 
                      lat_range: list[float, float] = None, lon_range: list[float, float] = None) -> xr.DataArray:
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

def convert_concentration_to_content(concentration: xr.DataArray, density: float = 1025.0) -> xr.DataArray:
    lat_values = concentration.lat.values
    lat_range = [float(lat_values.min()), float(lat_values.max())]
    lon_values = concentration.lon.values
    lon_range = [float(lon_values.min()), float(lon_values.max())]
    depth_values = concentration.depth.values
    depth_range = [float(depth_values.min()), float(depth_values.max())]

    areas = load_areas(lat_range=lat_range, lon_range=lon_range)
    depth, depth_bnds = load_depth_levels(depth_range=depth_range)
    thickness = depth_bnds[:, 1] - depth_bnds[:, 0]
    thickness_da = xr.DataArray(
        thickness.values,
        coords={'depth': depth.values},
        dims=['depth']
    )

    # in case data dimensions are not aligned, though it should not happen
    thickness_aligned = thickness_da.reindex(depth=concentration.depth, method='nearest', tolerance=0.1)
    areas_aligned = areas.reindex(lat=concentration.lat, lon=concentration.lon, method='nearest', tolerance=0.01)
    volume = thickness_aligned * areas_aligned

    # concentration (μmol/kg), density (kg/m³), volumn (m³), content (mol)
    content = concentration * density * volume * 1e-6
    content.attrs.update({
        'units': 'mol',
        'description': 'Content calculated from concentration * density * volume',
        'original_units_of_concentration': 'μmol/kg',
        'density_used': f'{density} kg/m³',
        'conversion_factor': '1e-6 (μmol to mol)'
    })

    if hasattr(concentration, 'name') and concentration.name:
        content.name = f'{concentration.name}_content'
    
    return content

path = "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_o2_sjtu_jingwei_v1_1.zarr"
da = load_regular_data(path, time_range=("2000-01-01", "2010-12-31"), depth_range=(0, 100), lat_range=(20, 30), lon_range=(120, 130))
# print(da)

content = convert_concentration_to_content(da, density=1025.0)
print(content)
