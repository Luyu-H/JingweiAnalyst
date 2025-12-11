# tools: compute_anomaly, convert_concentration_to_content, compute_vertical_gradient, 
# compute_ratio, identify_oxygen_minimum_zone, compute_potential_density_anomaly, compute_mixed_layer_depth,
# compute_nutricline_depth, compute_enso_nino_index

import os
import yaml
from pathlib import Path
import xarray as xr
import numpy as np
import gsw
from sklearn.decomposition import PCA

from src.tools.registry import register_tools
from src.tools.data_loader import load_areas, load_depth_levels
from src.tools.basic_calculator import compute_basic_stat


@register_tools(
    name="compute_anomaly",
    description="Compute anomaly of reconstruction data based on WOA23 climatology data.",
    args={
        "data": "xr.DataArray, required. Input reconstruction data.",
        "var": "str, required. Need to select from 'temp', 'so', 'si', 'o2', 'no3', 'po4'.",
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
        "temp": os.path.join(clim_base_path, "temperature_1971-2000"),
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


@register_tools(
    name="compute_vertical_gradient",
    description="Compute vertical gradient of a variable in an xarray.DataArray along the depth dimension. "
                "Gradient is computed as Δvar / Δz, consistent with interface method.",
    args={
        "data": "xr.DataArray, required. Input data with 'depth' dimension. Should have dimensions like (time, depth, lat, lon) or subset thereof."
    },
    returns="xr.DataArray. Vertical gradient (dVar/dDepth) with same dimensions as input. Units are [original_units]/meters."
)
def compute_vertical_gradient(data: xr.DataArray) -> xr.DataArray:
    """
    Compute vertical gradient using centered finite differences.
    - Interior points: gradient[i] = (data[i+1] - data[i-1]) / (depth[i+1] - depth[i-1])
    - Top boundary: gradient[0] = (data[1] - data[0]) / (depth[1] - depth[0])
    - Bottom boundary: gradient[-1] = (data[-1] - data[-2]) / (depth[-1] - depth[-2])
    """
    if 'depth' not in data.dims:
        raise ValueError("Input DataArray must have a 'depth' dimension.")

    # Load corresponding depth levels and boundaries
    depth_vals = data['depth'].values
    depth_range = [float(depth_vals.min()), float(depth_vals.max())]
    depth_da, depth_bnds_da = load_depth_levels(depth_range=depth_range)

    depth_levels = depth_da.values
    depth_bnds = depth_bnds_da.values

    if not np.allclose(depth_levels, depth_vals, atol=1e-4):
        lower_interp = np.interp(depth_vals, depth_levels, depth_bnds[:, 0])
        upper_interp = np.interp(depth_vals, depth_levels, depth_bnds[:, 1])
        depth_bnds = np.stack([lower_interp, upper_interp], axis=-1)

    # Compute interface values
    interfaces = np.empty(len(depth_bnds) + 1, dtype=np.float32)
    interfaces[0] = depth_bnds[0, 0]
    for i in range(len(depth_bnds)):
        interfaces[i + 1] = depth_bnds[i, 1]

    dz = depth_bnds[:, 1] - depth_bnds[:, 0]
    dz_safe = np.where(dz == 0, np.nan, dz)

    transposed = data.transpose(..., 'depth')
    v_np = transposed.values  # shape: (..., depth)
    n_depth = v_np.shape[-1]

    # Compute interface values for variable using linear interpolation between centers
    v_int = np.empty(v_np.shape[:-1] + (n_depth + 1,), dtype=np.float32)
    # Top/bottom extrapolation
    v_int[..., 0] = v_np[..., 0]
    v_int[..., -1] = v_np[..., -1]
    # Internal interfaces interpolation
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

    # Compute vertical gradients for each layer
    dv = v_int[..., 1:] - v_int[..., :-1]  # (..., depth)
    dvdz = dv / dz_safe.reshape((1,) * (v_np.ndim - 1) + (-1,))

    grad = xr.DataArray(
        dvdz,
        dims=transposed.dims,
        coords={dim: transposed.coords[dim] for dim in transposed.dims},
        attrs={
            "long_name": f"Vertical gradient of {data.name or 'variable'}",
            "units": f"{data.attrs.get('units', '')} / m",
            "method": "interface",
            "positive": "down"
        }
    ).transpose(*data.dims)

    return grad


@register_tools(
    name="compute_ratio",
    description="Compute ratio between two DataArrays: (data1 / data2). ",
    args={
        "data1": "xr.DataArray, required. Numerator DataArray.",
        "data2": "xr.DataArray, required. Denominator DataArray."
    },
    returns="xr.DataArray. Ratio DataArray with same dimensions as input."
)
def compute_ratio(data1: xr.DataArray, data2: xr.DataArray) -> xr.DataArray:
    """
    Compute ratio between two DataArrays: (data1 / data2)
    Handles division by zero by assigning NaN where data2 is zero.
    """
    ratio = xr.where(
        data2 != 0,
        (data1 / data2),
        np.nan
    )
    ratio.attrs.update({
        'long_name': f'Ratio of {data1.name} to {data2.name}',
        'description': 'Computed as (data1 / data2)'
    })
    return ratio


@register_tools(
    name="identify_oxygen_minimum_zone",
    description="Identify Oxygen Minimum Zones (OMZ) with dissolved oxygen concentration data and calculate their thickness. "
                "An OMZ is defined as areas where oxygen concentration is below the threshold.",
    args={
        "data": "xr.DataArray, required. Oxygen concentration data.",
        "threshold": "float, optional. Oxygen concentration threshold for OMZ identification in μmol/kg. Default is 20.0 μmol/kg.",
    },
    returns="tuple[xr.DataArray, xr.DataArray]. Returns (omz_thickness, omz_mask) where:\n"
            "  - omz_thickness: Thickness of OMZ in meters with dimensions (time, lat, lon) or subset\n"
            "  - omz_mask: Boolean array indicating OMZ locations with same dimensions as input data"
)
def identify_oxygen_minimum_zone(o2_data: xr.DataArray, threshold: float = 20.0) -> tuple[xr.DataArray, xr.DataArray]:
    depth_values = o2_data['depth'].values
    depth_range = [float(depth_values.min()), float(depth_values.max())]
    depth_da, depth_bnds_da = load_depth_levels(depth_range=depth_range)
    
    omz_mask = xr.where(o2_data.notnull(), o2_data < threshold, np.nan)

    layer_thickness = depth_bnds_da.values[:, 1] - depth_bnds_da.values[:, 0]
    thickness_da = xr.DataArray(
        layer_thickness,
        coords={'depth': depth_da},
        dims=['depth']
    )
    
    omz_thickness = (thickness_da * omz_mask).sum(dim='depth', skipna=True)
    omz_thickness = omz_thickness.where(o2_data.isel(depth=0).notnull())
    omz_thickness.attrs = {
        'long_name': 'Oxygen minimum zone thickness',
        'units': 'meters',
        'description': f'Total vertical thickness of layers where oxygen < {threshold} μmol/kg',
        'threshold': threshold,
        'threshold_units': 'μmol/kg'
    }
    
    omz_mask.attrs = {
        'long_name': 'Oxygen minimum zone mask',
        'description': f'Boolean mask indicating where oxygen < {threshold} μmol/kg',
        'threshold': threshold,
        'threshold_units': 'μmol/kg'
    }
    
    return omz_thickness, omz_mask


@register_tools(
    name="compute_potential_density_anomaly",
    description="Compute potential density anomaly from in-situ temperature and practical salinity. "
                "Uses TEOS-10 standard via GSW library. Potential density is referenced to surface pressure (0 dbar).",
    args={
        "temp_data": "xr.DataArray, required. In-situ temperature in degrees Celsius.",
        "sal_data": "xr.DataArray, required. Practical salinity (PSS-78).",
        "ref_pressure": "float, optional. Reference pressure in dbar for potential density calculation. Default is 0 dbar (sigma0)."
    },
    returns="xr.DataArray. Potential density anomaly (σ_θ = ρ - 1000) in kg/m³, with same dimensions as input."
)
def compute_potential_density_anomaly(temp_data: xr.DataArray, sal_data: xr.DataArray, ref_pressure: float = 0.0) -> xr.DataArray:
    if temp_data.dims != sal_data.dims:
        raise ValueError(f"Temperature and salinity must have same dimensions. "
                        f"Got temp: {temp_data.dims}, sal: {sal_data.dims}")

    if 'depth' not in temp_data.coords:
        raise ValueError("Temperature and salinity data must have 'depth' coordinate")

    depth = temp_data['depth'].values
    if 'lat' in temp_data.coords:
        lat = temp_data['lat'].values
        if lat.ndim == 0:  # scalar
            lat = float(lat)
    else:
        lat = 0.0  # default latitude
        
    if 'lon' in temp_data.coords:
        lon = temp_data['lon'].values
        if lon.ndim == 0:  # scalar
            lon = float(lon)
    else:
        lon = 0.0
    
    pressure = gsw.p_from_z(-depth, lat)
    if 'depth' in temp_data.dims:
        depth_axis = temp_data.dims.index('depth')
        # Create shape for broadcasting
        shape = [1] * temp_data.ndim
        shape[depth_axis] = len(depth)
        pressure = pressure.reshape(shape)
    
    absolute_salinity = gsw.SA_from_SP(sal_data.values, pressure, lon, lat)
    conservative_temperature = gsw.CT_from_t(absolute_salinity, temp_data.values, pressure)
    sigma_theta = gsw.sigma0(absolute_salinity, conservative_temperature) if ref_pressure == 0.0 else gsw.pot_rho_t_exact(absolute_salinity, temp_data.values, pressure, ref_pressure) - 1000.0

    rho_data = xr.DataArray(
        sigma_theta,
        coords=temp_data.coords,
        dims=temp_data.dims,
        attrs={
            'long_name': f'Potential density anomaly (sigma_{int(ref_pressure)})',
            'units': 'kg/m³',
            'reference_pressure': f'{ref_pressure} dbar',
            'description': f'Potential density anomaly (rho - 1000 kg/m³) referenced to {ref_pressure} dbar',
            'method': 'TEOS-10 (GSW library)',
            'source_temperature': 'in-situ temperature',
            'source_salinity': 'practical salinity (PSS-78)'
        }
    )
    
    return rho_data


def _find_threshold_depth(data: xr.DataArray, depth: np.ndarray, threshold: float, ref_depth: float) -> xr.DataArray:
    """
    Helper function to find depth where data exceeds threshold.

    Args:
        - data (xr.DataArray): Data array with depth dimension (e.g., density difference)
        - depth (np.ndarray): Depth values, m
        - threshold (float): Threshold value
        - ref_depth (float): Reference depth, m

    Returns:
        - xr.DataArray: Depth (m, from existing discrete depth levels) where threshold is first exceeded
    """
    exceeds_threshold = data >= threshold

    depth_dim = data.dims.index('depth')

    # xarray -> numpy
    data_values = data.values
    threshold_mask = exceeds_threshold.values

    output_shape = list(data_values.shape)
    output_shape.pop(depth_dim)
    mld_values = np.full(output_shape, np.nan)

    # Get indices for iteration (all dimensions except depth)
    other_dims = list(range(data_values.ndim))
    other_dims.pop(depth_dim)
    for idx in np.ndindex(*[data_values.shape[i] for i in other_dims]):
        full_idx = list(idx)
        full_idx.insert(depth_dim, slice(None))
        full_idx = tuple(full_idx)

        mask_profile = threshold_mask[full_idx]

        exceeded_indices = np.where(mask_profile)[0]
        if len(exceeded_indices) > 0:
            first_idx = exceeded_indices[0]
            d_mld = depth[first_idx]

            if d_mld >= ref_depth:
                output_idx = list(idx)
                output_idx = tuple(output_idx)
                mld_values[output_idx] = d_mld

    mld_array = xr.DataArray(
        mld_values,
        coords={k: v for k, v in data.coords.items() if k != 'depth'},
        dims=[d for d in data.dims if d != 'depth']
    )

    return mld_array

@register_tools(
    name="compute_mixed_layer_depth",
    description="Calculate mixed layer depth (MLD) using density threshold criterion. "
                "MLD is defined as the depth where potential density exceeds the near-surface "
                "reference density by a specified threshold (default: 0.03 kg/m³).",
    args={
        "temp_data": "xr.DataArray, required. In-situ temperature in degrees Celsius.",
        "sal_data": "xr.DataArray, required. Practical salinity (PSS-78, unitless) with same dimensions as temp_data.",
        "ref_depth": "float, optional. Reference depth for surface density in meters. Default is 10 m.",
        "method": "str, optional. Method for MLD calculation: 'density' (default), 'temperature', or 'both'.",
        "density_threshold": "float, optional. Density difference threshold in kg/m³, used for method 'density'. Default is 0.03 kg/m³.",
        "temp_threshold": "float, optional. Temperature difference threshold in °C, used for method 'temperature'. Default is 0.2 °C."
    },
    returns="xr.DataArray. Mixed layer depth in meters with dimensions excluding depth."
)
def compute_mixed_layer_depth(temp_data: xr.DataArray, sal_data: xr.DataArray, 
                              ref_depth: float = 10.0, method: str = 'density', density_threshold: float = 0.03, temp_threshold: float = 0.2) -> xr.DataArray:
    if temp_data.dims != sal_data.dims:
        raise ValueError(f"Temperature and salinity must have same dimensions. "
                        f"Got temp: {temp_data.dims}, sal: {sal_data.dims}")
    
    if 'depth' not in temp_data.coords:
        raise ValueError("Data must have 'depth' coordinate")
    
    if method not in ['density', 'temperature', 'both']:
        raise ValueError(f"Method must be 'density', 'temperature', or 'both', got '{method}'")
    
    depth = temp_data['depth'].values
    if ref_depth < depth.min() or ref_depth > depth.max():
        raise ValueError(f"Reference depth {ref_depth}m is outside data depth range "
                        f"[{depth.min()}, {depth.max()}]m")
    
    sigma_theta = compute_potential_density_anomaly(temp_data, sal_data, ref_pressure=0.0)

    ref_idx = np.argmin(np.abs(depth - ref_depth))
    actual_ref_depth = depth[ref_idx]

    sigma_ref = sigma_theta.isel(depth=ref_idx)
    temp_ref = temp_data.isel(depth=ref_idx)

    if method in ['density', 'both']:
        delta_sigma = sigma_theta - sigma_ref
        mld_density = _find_threshold_depth(
            delta_sigma, 
            depth, 
            threshold=density_threshold,
            ref_depth=actual_ref_depth
        )
    if method in ['temperature', 'both']:
        delta_temp = np.abs(temp_data - temp_ref)
        mld_temperature = _find_threshold_depth(
            delta_temp,
            depth,
            threshold=temp_threshold,
            ref_depth=actual_ref_depth
        )
    
    if method == 'density':
        mld = mld_density
        method_description = f'Density threshold: Δσ_θ = {density_threshold} kg/m³'
    elif method == 'temperature':
        mld = mld_temperature
        method_description = f'Temperature threshold: ΔT = {temp_threshold} °C'
    else:  # both
        # Take minimum (shallowest) MLD from both criteria
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
    
    return mld_data


@register_tools(
    name="compute_nutricline_depth",
    description="Compute nutricline depth as the depth where vertical gradient (dN/dz) is maximum. ",
    args={
        "data": "xr.DataArray, required. Input nutrient concentration data with 'depth' dimension. "
                "Should have dimensions like (time, depth, lat, lon) or subset thereof.",
        "min_depth": "float, optional. Minimum depth (m) to search for nutricline. Default is 0.",
        "max_depth": "float, optional. Maximum depth (m) to search for nutricline. Default is None (search all depths)."
    },
    returns="xr.DataArray. Nutricline depth with same dimensions as input except 'depth' is removed. "
            "Units are meters. NaN where gradient is not positive or nutricline cannot be determined."
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
        # Create indexing for all dimensions
        mask_full = xr.DataArray(
            depth_mask,
            coords={'depth': gradient['depth']},
            dims=['depth']
        )
        gradient_masked = gradient_masked.where(mask_full, np.nan)

    max_idx = gradient_masked.argmax(dim='depth', skipna=True)
    nutricline_depth = gradient['depth'].isel(depth=max_idx)
    
    # Set NaN where gradient is not positive or all NaN
    max_gradient = gradient_masked.max(dim='depth', skipna=True)
    all_nan = gradient_masked.isnull().all(dim='depth')
    nutricline_depth = nutricline_depth.where(
        (max_gradient > 0) & (~all_nan),
        np.nan
    )
    
    nutricline_depth.attrs = {
        'long_name': f'Nutricline depth ({data.name or "nutrient"})',
        'standard_name': 'nutricline_depth',
        'units': 'meters',
        'description': 'Depth of maximum vertical gradient in nutrient concentration',
        'method': 'maximum_gradient',
        'positive': 'down',
        'source_variable': data.name or 'unknown'
    }
    nutricline_depth.name = f'{data.name}_nutricline_depth' if data.name else 'nutricline_depth'
    
    return nutricline_depth


@register_tools(
    name="compute_enso_nino_index",
    description="Compute ENSO Niño index for a specified region from monthly temperature data. "
                "Calculates area-weighted mean temperature anomaly for the selected Niño region.",
    args={
        "data": "xr.DataArray, required. Input temperature data with dimensions (time, depth, lat, lon).",
        "nino_region": "str, optional. Niño region to compute. Choose from 'nino1+2', 'nino3', 'nino3.4', 'nino4'. Default is 'nino3.4'.",
        "depth_level": "float, optional. Average across depth levels (in meters) that are shallower than depth_level in the data. Default is 10 (surface). ",
        "anomaly_freq": "str, optional. Frequency for climatology when computing anomaly. Choose from 'month', 'season', 'year'. Default is 'month'."
    },
    returns="xr.DataArray. Time series of Niño index (temperature anomaly in °C) for the specified region."
)
def compute_enso_nino_index(data: xr.DataArray, nino_region: str = "nino3.4", depth_level: float = 10.0, anomaly_freq: str = "month") -> xr.DataArray:
    # (lat_min, lat_max, lon_min, lon_max)
    nino_regions = {
        "nino1+2": (-10, 0, -90, -80),
        "nino3": (-5, 5, -150, -90),
        "nino3.4": (-5, 5, -170, -120),
        "nino4": (-5, 5, 160, -150)  # crosses dateline
    }

    if nino_region not in nino_regions:
        raise ValueError(f"nino_region must be one of {list(nino_regions.keys())}, got '{nino_region}'")

    lat_min, lat_max, lon_min, lon_max = nino_regions[nino_region]
    if lon_max < lon_min:  # Crosses dateline (e.g., Niño 4: 160 to -150)
        data_region1 = data.sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, 180))
        data_region2 = data.sel(lat=slice(lat_min, lat_max), lon=slice(-180, lon_max))
        data_region = xr.concat([data_region1, data_region2], dim='lon')
    else:
        data_region = data.sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))
    
    anomaly = compute_anomaly(data=data_region, var="temp", freq=anomaly_freq, method="subtract")

    if 'depth' in anomaly.dims:
        anomaly_surface = anomaly.sel(depth=slice(0, depth_level))
        # Average over depth dimension
        anomaly_surface = anomaly_surface.mean(dim='depth', skipna=True)
    else:
        anomaly_surface = anomaly
    
    nino_index = compute_basic_stat(data=anomaly_surface, op="area_weighted_mean", dim=["lat", "lon"])

    nino_index.attrs = {
        'long_name': f'{nino_region.upper()} Index',
        'units': '°C',
        'description': f'Area-weighted sea surface temperature anomaly for {nino_region.upper()} region',
        'region_bounds': f'lat: [{lat_min}, {lat_max}], lon: [{lon_min}, {lon_max}]',
        'depth_range': f'0-{depth_level}m averaged',
        'anomaly_frequency': anomaly_freq,
        'anomaly_method': 'subtract',
        'weighting': 'area weighting',
        'reference': 'WOA23 climatology'
    }
    nino_index.name = f"{nino_region}_index"
    
    return nino_index


@register_tools(
    name="compute_pdo_index",
    description="Compute Pacific Decadal Oscillation (PDO) index from sea surface temperature data. "
                "PDO is defined as the leading principal component of North Pacific SST anomalies (poleward of 20°N). "
                "The SST anomalies are computed by removing the global mean SST to isolate regional patterns and temperal trends.",
    args={
        "data": "xr.DataArray, required. Temperature data with dimensions (time, lat, lon, depth). "
                    "Must cover at least the North Pacific region (20°N-70°N, 110°E-100°W). ",
        "depth_level": "float, optional. Average across depth levels (in meters) that are shallower than depth_level in the data. Default is 10 (surface). ",
        "standardize": "bool, optional. Whether to standardize the PDO index to zero mean and unit variance. "
                       "Default is True."
    },
    returns="tuple[xr.DataArray, xr.DataArray, float]. Returns (pdo_index, eof_pattern, explained_variance) where: "
            "pdo_index is 1D time series of PDO values, eof_pattern is 2D spatial pattern (lat, lon) of the first EOF, "
            "explained_variance is the percentage of variance explained by the first EOF."
)
def compute_pdo_index(data: xr.DataArray, depth_level: float = 10.0, standardize: bool = True) -> tuple[xr.DataArray, xr.DataArray, float]:
    # Extract surface temperature
    if 'depth' in data.dims:
        sst_data = data.sel(depth=slice(0, depth_level))
        sst_data = sst_data.mean(dim='depth', skipna=True)
    
    # Subset to North Pacific region (20°N-70°N, 110°E-100°W)
    pacific_east = sst_data.sel(lat=slice(20, 70), lon=slice(110, 180))
    pacific_west = sst_data.sel(lat=slice(20, 70), lon=slice(-180, -100))
    north_pacific = xr.concat([pacific_east, pacific_west], dim='lon')
    
    if len(north_pacific.time) == 0 or len(north_pacific.lat) == 0 or len(north_pacific.lon) == 0:
        raise ValueError("No data in North Pacific region (20°N-70°N, 110°E-260°E)")
    
    # Remove global mean SST
    global_mean = compute_basic_stat(sst_data, op='area_weighted_mean', dim=['lat', 'lon'])
    north_pacific = north_pacific - global_mean
    
    # Get grid cell areas for North Pacific region
    lat_range = [float(north_pacific.lat.min()), float(north_pacific.lat.max())]
    lon_range = [float(north_pacific.lon.min()), float(north_pacific.lon.max())]
    areas = load_areas(lat_range, lon_range)
    areas = areas.sel(lat=north_pacific.lat, lon=north_pacific.lon, method="nearest")
    
    # Apply area weighting (square root for proper EOF analysis)
    weights = np.sqrt(areas.values)
    weighted_data = north_pacific * weights
    
    # Remove time mean to get anomalies
    weighted_anomalies = weighted_data - weighted_data.mean(dim='time')
    
    # Reshape data to (time, space) for EOF analysis
    n_times = len(weighted_anomalies.time)
    n_lat = len(weighted_anomalies.lat)
    n_lon = len(weighted_anomalies.lon)
    data_matrix = weighted_anomalies.values.reshape(n_times, n_lat * n_lon)
    
    # Remove any spatial points that are all NaN
    valid_points = ~np.isnan(data_matrix).all(axis=0)
    data_matrix_clean = data_matrix[:, valid_points]
    
    # Fill any remaining NaNs with 0 for PCA
    data_matrix_clean = np.nan_to_num(data_matrix_clean, nan=0.0)
    
    # Perform EOF analysis using PCA
    pca = PCA(n_components=1)
    pdo_index_values = pca.fit_transform(data_matrix_clean).flatten()
    
    # Get EOF pattern
    eof_pattern_1d = np.full(n_lat * n_lon, np.nan)
    eof_pattern_1d[valid_points] = pca.components_[0, :]
    eof_pattern_2d = eof_pattern_1d.reshape(n_lat, n_lon)
    
    # Unweight the EOF pattern for interpretation
    valid_mask = ~np.isnan(weights.values)
    eof_pattern_2d[valid_mask] = eof_pattern_2d[valid_mask] / weights.values[valid_mask]
    
    # Explained variance
    explained_var = float(pca.explained_variance_ratio_[0] * 100)
    
    if standardize:
        pdo_index_values = (pdo_index_values - np.nanmean(pdo_index_values)) / np.nanstd(pdo_index_values)
    
    pdo_index = xr.DataArray(
        pdo_index_values,
        coords={'time': north_pacific.time},
        dims=['time'],
        attrs={
            'long_name': 'Pacific Decadal Oscillation Index',
            'description': 'Leading principal component of North Pacific SST anomalies (20°N-70°N)',
            'units': 'standardized' if standardize else 'dimensionless',
            'method': 'EOF analysis of area-weighted SST anomalies',
            'region': 'North Pacific (20°N-70°N, 110°E-100°W)',
            'standardized': standardize
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
            'explained_variance_percent': explained_var
        }
    )
    
    return pdo_index, eof_pattern, explained_var