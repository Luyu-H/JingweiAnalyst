from src.tools.data_loader import load_data_with_geolocator
from src.tools.basic_calculator import temporal_aggregate
from src.tools.scientific_calculator import compute_ratio

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd

ocean_basins = {
    "Southern Ocean": "SOUTHERN OCEAN",
    "Arctic Ocean": "Arctic Ocean", 
    "Pacific Ocean": ["North Pacific Ocean", "South Pacific Ocean"],
    "Atlantic Ocean": ["North Atlantic Ocean", "South Atlantic Ocean"],
    "Indian Ocean": "INDIAN OCEAN"
}

no3_path = "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_no3_sjtu_jingwei_v1_1.zarr"
po4_path = "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_po4_sjtu_jingwei_v1_1.zarr"

time_range = ["1975-01-01", "1985-12-31"]
#time_range = ["1975-01-01", "1975-12-31"]
depth_range = [0, 1000]

basin_np_ratios = {}

for basin, value in ocean_basins.items():
    print(f"Processing {basin}...")
    
    no3_data = load_data_with_geolocator(
        zarr_path=no3_path,
        time_range=time_range,
        depth_range=depth_range,
        geofeature="water_body",
        loc_name=value
    )
    
    po4_data = load_data_with_geolocator(
        zarr_path=po4_path,
        time_range=time_range,
        depth_range=depth_range,
        geofeature="water_body",
        loc_name=value
    )
    
    no3_clim = temporal_aggregate(no3_data, freq='climatology', op='mean')
    po4_clim = temporal_aggregate(po4_data, freq='climatology', op='mean')
    no3_jan_mean = no3_clim.sel(month=1)
    po4_jan_mean = po4_clim.sel(month=1)
    
    np_ratio = compute_ratio(no3_jan_mean, po4_jan_mean)
    
    # Flatten and remove NaN values for boxplot
    ratio_values = np_ratio.values.flatten()
    ratio_values = ratio_values[~np.isnan(ratio_values)]
    
    # Remove extreme outliers (keep values within reasonable range)
    ratio_values = ratio_values[(ratio_values > 0) & (ratio_values < 100)]
    
    basin_np_ratios[basin] = ratio_values
    
    print(f"  {basin}: n={len(ratio_values)}, "
          f"median={np.median(ratio_values):.2f}, "
          f"mean={np.mean(ratio_values):.2f}")


width_cm = 9
height_cm = 5
width_inch = width_cm / 2.54
height_inch = height_cm / 2.54
plt.rcParams['font.family'] = 'DejaVu Serif'
fig, ax = plt.subplots(figsize=(width_inch*2.2, height_inch*2.2))

# Prepare data for boxplot
basin_names = list(basin_np_ratios.keys())
basin_data = [basin_np_ratios[basin] for basin in basin_names]

# Create boxplot with custom styling
bp = ax.boxplot(
    basin_data,
    labels=basin_names,
    patch_artist=True,
    showfliers=True,
    flierprops=dict(marker='o', markersize=3, alpha=0.3),
    medianprops=dict(color='black', linewidth=2),
    whiskerprops=dict(linewidth=1.5),
    capprops=dict(linewidth=1.5)
)

# Define colors for each ocean basin
colors = ['#8dd3c7', '#80b1d3', '#fb8072', '#bebada', '#fdb462']

for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

ax.axhline(y=16, color='red', linestyle='--', linewidth=2, label='Redfield Ratio (16:1)')

ax.set_ylim(top=30)

ax.tick_params(axis='x', labelsize=13)
ax.tick_params(axis='y', labelsize=13)
ax.set_xlabel('Ocean Basin', fontsize=14)
ax.set_ylabel('N:P Ratio', fontsize=14)
ax.set_title('N:P Ratio Comparison Across Ocean Basins', fontsize=16)

ax.yaxis.grid(True, linestyle='--', alpha=0.3)
ax.set_axisbelow(True)

ax.legend(loc='upper right', fontsize=14)

# Rotate x-axis labels for better readability
plt.xticks(rotation=15, ha='right')

# for i, basin in enumerate(basin_names):
#     n = len(basin_np_ratios[basin])
#     median_val = np.median(basin_np_ratios[basin])
#     ax.text(i + 1, ax.get_ylim()[1] * 0.95, f'n={n}', 
#             ha='center', va='top', fontsize=9, color='gray')

plt.tight_layout()

output_path = "/home/lorraine/jingwei_QA/src/applications/res/4_np_ratio_ocean_basins_boxplot.png"
plt.savefig(output_path, dpi=450, facecolor='white')