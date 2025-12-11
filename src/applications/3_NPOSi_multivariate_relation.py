from src.tools.data_loader import load_data_with_geolocator
from src.tools.basic_calculator import compute_basic_stat, fit_multivariate_relation

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

no3_path = "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_no3_sjtu_jingwei_v1_1.zarr"
po4_path = "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_po4_sjtu_jingwei_v1_1.zarr"
sio2_path = "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_si_sjtu_jingwei_v1_1.zarr"
o2_path = "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_o2_sjtu_jingwei_v1_1.zarr"

time_range = ["1960-01-01", "1965-12-31"]
geofeature_name = "water_body"
loc_name = "North Pacific Ocean"

no3_data = load_data_with_geolocator(
    zarr_path=no3_path,
    time_range=time_range,
    geofeature=geofeature_name,
    loc_name=loc_name
)
no3_data.name = "NO3"

po4_data = load_data_with_geolocator(
    zarr_path=po4_path,
    time_range=time_range,
    geofeature=geofeature_name,
    loc_name=loc_name
)
po4_data.name = "PO4"

sio2_data = load_data_with_geolocator(
    zarr_path=sio2_path,
    time_range=time_range,
    geofeature=geofeature_name,
    loc_name=loc_name
)
sio2_data.name = "SiO2"

o2_data = load_data_with_geolocator(
    zarr_path=o2_path,
    time_range=time_range,
    geofeature=geofeature_name,
    loc_name=loc_name
)
o2_data.name = "O2"

print(f"Data dimensions: {o2_data.dims}")
print(f"Data shape: {o2_data.shape}")
print(f"Number of depth layers: {len(o2_data.depth)}")

no3_mean = compute_basic_stat(no3_data, op="mean", dim=['lat', 'lon'])
po4_mean = compute_basic_stat(po4_data, op="mean", dim=['lat', 'lon'])
sio2_mean = compute_basic_stat(sio2_data, op="mean", dim=['lat', 'lon'])
o2_mean = compute_basic_stat(o2_data, op="mean", dim=['lat', 'lon'])

# Preserve variable names
no3_mean.name = "NO3"
po4_mean.name = "PO4"
sio2_mean.name = "SiO2"
o2_mean.name = "O2"

depths = o2_mean.depth.values
n_depths = len(depths)

# Store results
results_by_depth = {
    'depth': [],
    'no3_coef': [],
    'po4_coef': [],
    'sio2_coef': [],
    'intercept': [],
    'r2': [],
    'test_r2': [],
    'rmse': [],
    'n_samples': []
}

for i, depth_val in enumerate(depths):
    print(f"  Processing depth {depth_val:.1f} m ({i+1}/{n_depths})...")
    
    # Extract current depth layer data
    o2_layer = o2_mean.sel(depth=depth_val)
    no3_layer = no3_mean.sel(depth=depth_val)
    po4_layer = po4_mean.sel(depth=depth_val)
    sio2_layer = sio2_mean.sel(depth=depth_val)
    
    # Preserve variable names
    no3_layer.name = "NO3"
    po4_layer.name = "PO4"
    sio2_layer.name = "SiO2"
    
    try:
        # Perform multivariate linear regression
        reg_result = fit_multivariate_relation(
            target=o2_layer,
            predictors=[no3_layer, po4_layer, sio2_layer],
            method="linear"
        )
        
        # Extract results
        results_by_depth['depth'].append(depth_val)
        results_by_depth['no3_coef'].append(reg_result['coefficients'].get('NO3', np.nan))
        results_by_depth['po4_coef'].append(reg_result['coefficients'].get('PO4', np.nan))
        results_by_depth['sio2_coef'].append(reg_result['coefficients'].get('SiO2', np.nan))
        results_by_depth['intercept'].append(reg_result['coefficients'].get('intercept', np.nan))
        results_by_depth['r2'].append(reg_result['metrics']['r2'])
        results_by_depth['test_r2'].append(reg_result['metrics']['test_r2'])
        results_by_depth['rmse'].append(reg_result['metrics']['rmse'])
        results_by_depth['n_samples'].append(reg_result['metrics']['n_samples'])
        
    except Exception as e:
        print(f"    Regression failed at depth {depth_val:.1f} m: {e}")
        results_by_depth['depth'].append(depth_val)
        results_by_depth['no3_coef'].append(np.nan)
        results_by_depth['po4_coef'].append(np.nan)
        results_by_depth['sio2_coef'].append(np.nan)
        results_by_depth['intercept'].append(np.nan)
        results_by_depth['r2'].append(np.nan)
        results_by_depth['test_r2'].append(np.nan)
        results_by_depth['rmse'].append(np.nan)
        results_by_depth['n_samples'].append(0)

for key in results_by_depth:
    results_by_depth[key] = np.array(results_by_depth[key])

print("\nRegression analysis completed!")


fig, ax = plt.subplots(figsize=(10, 8))

depths_plot = results_by_depth['depth']

colors = {
    'NO3': '#e41a1c',
    'PO4': '#377eb8',
    'SiO2': '#4daf4a',
    'Intercept': '#984ea3'
}

ax.plot(results_by_depth['no3_coef'], depths_plot, 'o-', 
        color=colors['NO3'], label='NO₃⁻', linewidth=2, markersize=5)
ax.plot(results_by_depth['po4_coef'], depths_plot, 's-', 
        color=colors['PO4'], label='PO₄³⁻', linewidth=2, markersize=5)
ax.plot(results_by_depth['sio2_coef'], depths_plot, '^-', 
        color=colors['SiO2'], label='SiO₂', linewidth=2, markersize=5)
ax.plot(results_by_depth['intercept'], depths_plot, 'd-', 
        color=colors['Intercept'], label='Intercept (β₀)', linewidth=2, markersize=5)

ax.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.7)

ax.set_xlabel('Regression Coefficient', fontsize=12)
ax.set_ylabel('Depth (m)', fontsize=12)
ax.invert_yaxis()

ax.set_title('Regression Coefficients of Nutrients on Dissolved Oxygen\n'
             'O₂ = β₀ + β₁·NO₃⁻ + β₂·PO₄³⁻ + β₃·SiO₂',
             fontsize=13, fontweight='bold')
ax.legend(loc='lower right', fontsize=10, framealpha=0.9)
ax.grid(True, alpha=0.3)
ax.tick_params(axis='both', labelsize=10)

plt.tight_layout()

# Save figure
output_path = "/home/lorraine/jingwei_QA/src/applications/res/3_nutrient_o2_regression_coefficients.png"
plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')