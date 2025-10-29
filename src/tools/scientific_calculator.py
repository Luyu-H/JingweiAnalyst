# tools: compute_anomaly, convert_concentration_to_content

import os
import yaml
from pathlib import Path
import xarray as xr
import numpy as np

from src.tools.registry import register_tools
from src.tools.data_loader import load_areas, load_depth_levels


@register_tools(
    name="compute_anomaly",
    description="Compute anomaly of reconstruction data based on WOA23 climatology data.",
    args={
        "data": "xr.DataArray, required. Input reconstruction data.",
        "var": "str, required. Need to select from 'thetao', 'so', 'si', 'o2', 'no3', 'po4'.",
        "freq": "str, optional. Frequency of climatology data. Choose from 'month', 'season', 'year'. Default is 'month'.",
        "method": "str, optional. Method to compute anomaly. Choose from 'subtract' (data - climatology) "
                  "or 'ratio' ((data - climatology) / climatology * 100). Default is 'subtract'."
    },
    returns="xr.DataArray. Anomaly data with same dimensions as input."
)
def compute_anomaly(data: xr.DataArray, var: str, freq: str = "month", method: str = "subtract") -> xr.DataArray:
    if freq not in ["month", "season", "year"]:
        raise ValueError(f"freq must be 'month', 'season', or 'year', got '{freq}'")
    if method not in ["subtract", "ratio"]:
        raise ValueError(f"method must be 'subtract' or 'ratio', got '{method}'")
    
    config_path = "/home/lorraine/jingwei_QA/configs/default_paths.yaml"
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)
    clim_base_path = cfg.get("climatology_data_path")
    
    var_path_map = {
        "thetao": os.path.join(clim_base_path, "temperature_1971-2000"),
        "so": os.path.join(clim_base_path, "salinity_1971-2000"),
        "si": os.path.join(clim_base_path, "silicate_all_time_periods"),
        "o2": os.path.join(clim_base_path, "oxygen_all_time_periods"),
        "no3": os.path.join(clim_base_path, "nitrate_all_time_periods"),
        "po4": os.path.join(clim_base_path, "phosphate_all_time_periods")
    }
    if var not in var_path_map:
        raise ValueError(f"Variable '{var}' not supported. Available: {list(var_path_map.keys())}")
    
    # Map frequency to file suffix codes
    freq_to_codes = {
        "year": ["00"], # annual
        "month": [f"{i:02d}" for i in range(1, 13)], # 01-12
        "season": ["13", "14", "15", "16"] # winter, spring, summer, autumn
    }
    
    clim_dir = var_path_map[var]
    clim_files = list(Path(clim_dir).glob("*.nc"))
    if not clim_files:
        raise FileNotFoundError(f"No climatology files found in {clim_dir}")
    
    sample_file = clim_files[0].stem # woa23_all_o00_01
    # Pattern: woa23_*_01.nc -> extract prefix before last underscore and number
    parts = sample_file.rsplit('_', 2) # ['woa23_all', 'o00', '01']

    if len(parts) >= 3:
        base_name = '_'.join(parts[:-2]) # 'woa23_all'
        var_code = parts[-2] # 'o00'
        clim_pattern = f"{base_name}_{var_code[:-2]}" # 'woa23_all_o'
    else:
        # Fallback: try to find pattern
        clim_pattern = sample_file.rsplit('_', 1)[0] + '_'
    
    # Load climatology data
    climatology_dict = {}
    for code in freq_to_codes[freq]:
        clim_file = os.path.join(clim_dir, f"{clim_pattern}{code}_01.nc")
        if os.path.exists(clim_file):
            ds = xr.open_dataset(clim_file, decode_times=False)
            # Find the main variable (usually ends with _an for analyzed mean)
            var_name = [var for var in ds.data_vars if var.endswith("_an")]
            if var_name:
                climatology_dict[code] = ds[var_name[0]]
            ds.close()
    
    if not climatology_dict:
        raise FileNotFoundError(f"No valid climatology files found for freq='{freq}'")
    
    # Create time-based mapping for anomaly calculation
    def get_time_code(time_val):
        """Map datetime to climatology code based on frequency"""
        if freq == "year":
            return "00"
        elif freq == "month":
            return f"{time_val.month:02d}"
        elif freq == "season":
            month = time_val.month
            if month in [12, 1, 2]:
                return "13" # winter
            elif month in [3, 4, 5]:
                return "14" # spring
            elif month in [6, 7, 8]:
                return "15" # summer
            else:
                return "16" # autumn
    
    anomaly_list = []
    
    for time_idx in range(len(data.time)):
        time_val = data.time.values[time_idx]
        time_dt = xr.DataArray(time_val).dt
        code = get_time_code(time_dt)
        
        if code not in climatology_dict:
            raise ValueError(f"Climatology for code '{code}' not found")
        
        clim = climatology_dict[code]
        data_slice = data.isel(time=time_idx)
        
        # Calculate anomaly
        if method == "subtract":
            anom = data_slice - clim
        elif method == "ratio":
            anom = xr.where(
                clim != 0,
                (data_slice - clim) / clim * 100,
                np.nan
            )
        
        anomaly_list.append(anom)
    
    # Concatenate along time dimension
    anomaly_data = xr.concat(anomaly_list, dim=data.time)
    
    anomaly_data.attrs = data.attrs.copy()
    anomaly_data.attrs['anomaly_method'] = method
    anomaly_data.attrs['anomaly_freq'] = freq
    anomaly_data.attrs['climatology_reference'] = clim_dir
    
    if method == "subtract":
        anomaly_data.attrs['long_name'] = f"{data.attrs.get('long_name', var)} anomaly"
    else:
        anomaly_data.attrs['long_name'] = f"{data.attrs.get('long_name', var)} percentage anomaly"
        anomaly_data.attrs['units'] = 'percent'
    
    return anomaly_data


@register_tools(
    name="convert_concentration_to_content",
    description="Convert concentration data (μmol/kg) to content data (mol) by multiplying with grid cell volume. ",
    args={
        "concentration": "xr.DataArray, required. Concentration data with dimensions (time, depth, lat, lon) in μmol/kg.",
        "density": "float, optional. Seawater density in kg/m³. Default is 1025 kg/m³."
    },
    returns="xr.DataArray. Content data in moles (mol) with same dimensions as input."
)
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