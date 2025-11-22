from src.tools.data_loader import load_data_with_geolocator
from src.tools.basic_calculator import temporal_aggregate, compute_basic_stat

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

zarr_path = "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_so_sjtu_jingwei_v1_1.zarr"
time_range = ["1990-01-01", "2020-12-31"]
# Load data using geolocator to select Indian Ocean region
sal_data = load_data_with_geolocator(
    zarr_path=zarr_path,
    time_range=time_range,
    depth_range=(0, 50),  # Surface layer
    geofeature="water_body",
    loc_name="INDIAN OCEAN"
)

print(f"Data dimensions: {sal_data.dims}")
print(f"Data shape: {sal_data.shape}")
print(f"Time range: {sal_data.time.values[0]} to {sal_data.time.values[-1]}")
print(f"Depth range: {float(sal_data.depth.min())} to {float(sal_data.depth.max())} m")

surface_sal = compute_basic_stat(data=sal_data, op="mean", dim="depth")

print(f"Surface salinity dimensions: {surface_sal.dims}")
print(f"Surface salinity shape: {surface_sal.shape}")

indian_ocean_mean = compute_basic_stat(
    data=surface_sal,
    op="area_weighted_mean",
    dim=["lat", "lon"]
)

print(f"Indian Ocean mean salinity dimensions: {indian_ocean_mean.dims}")
print(f"Indian Ocean mean salinity shape: {indian_ocean_mean.shape}")

seasonal_cycle = temporal_aggregate(
    data=indian_ocean_mean,
    freq="climatology",
    op="mean"
)

print(f"Seasonal cycle dimensions: {seasonal_cycle.dims}")
print(f"Seasonal cycle shape: {seasonal_cycle.shape}")

months = seasonal_cycle.month.values
month_names = seasonal_cycle.month_name.values
salinity_values = seasonal_cycle.values

print("\n=== Monthly Climatological Salinity ===")
for i, (month, name, sal) in enumerate(zip(months, month_names, salinity_values)):
    print(f"{name}: {sal:.3f} PSU")


width_cm = 9
height_cm = 5
width_inch = width_cm / 2.54
height_inch = height_cm / 2.54
plt.rcParams['font.family'] = 'DejaVu Serif'
fig, ax = plt.subplots(figsize=(width_inch*2.2, height_inch*2.2))

# Color scheme: ocean blue gradient
line_color = '#1f77b4'
marker_color = '#fb8072'
fill_color = '#a6cee3'

# Plot the seasonal cycle line
ax.plot(
    months, salinity_values,
    color=line_color,
    linewidth=1.5,
    marker='o',
    markersize=6,
    markerfacecolor=marker_color,
    markeredgecolor='white',
    markeredgewidth=1,
    label='Monthly Mean Salinity',
    zorder=3,
    alpha=0.7
)

# Add subtle fill below the line
ax.fill_between(
    months, 
    np.min(salinity_values) - 0.02,
    salinity_values,
    color=fill_color,
    alpha=0.3,
    zorder=1
)

# Add data labels on each point
for i, (month, sal) in enumerate(zip(months, salinity_values)):
    ax.annotate(
        f'{sal:.3f}',
        (month, sal),
        textcoords="offset points",
        xytext=(0, 12),
        ha='center',
        fontsize=9,
        color='#333333'
    )

# Customize x-axis with month names
ax.set_xticks(months)
ax.set_xticklabels(month_names, fontsize=13)

ax.tick_params(axis='x', labelsize=13)
ax.tick_params(axis='y', labelsize=13)

# Labels and title
ax.set_xlabel('Month', fontsize=14)
ax.set_ylabel('Salinity (PSU)', fontsize=14)
ax.set_title(
    'Indian Ocean Surface Salinity Seasonal Cycle',
    fontsize=16
)

# Set y-axis limits with padding
y_min = np.min(salinity_values)
y_max = np.max(salinity_values)
y_range = y_max - y_min
ax.set_ylim(y_min - 0.1 * y_range, y_max + 0.15 * y_range)

# Add grid
ax.grid(True, which='major', alpha=0.3, zorder=0)
ax.set_axisbelow(True)

# ax.legend(loc='upper right', fontsize=14)

# Add statistics annotation
sal_mean = np.mean(salinity_values)
sal_range = y_max - y_min
stats_text = f'Annual Mean: {sal_mean:.2f} PSU\nSeasonal Range: {sal_range:.2f} PSU'
ax.text(
    0.05, 0.92,
    stats_text,
    transform=ax.transAxes,
    fontsize=13,
    verticalalignment='top',
    bbox=dict(
        boxstyle='round',
        facecolor='white',
        edgecolor='lightgray',
        linewidth=0.8,
        alpha=0.9
    )
)

# Spine customization
for spine in ['top', 'right']:
    ax.spines[spine].set_visible(False)
ax.spines['left'].set_linewidth(1.2)
ax.spines['bottom'].set_linewidth(1.2)

plt.tight_layout()

output_path = "/home/lorraine/jingwei_QA/src/applications/res/7_indian_ocean_salinity_seasonal_cycle.png"
plt.savefig(output_path, dpi=450, facecolor='white')