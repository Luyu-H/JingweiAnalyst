from src.tools.data_loader import load_regular_data
from src.tools.scientific_calculator import compute_enso_nino_index
from src.tools.basic_calculator import compute_basic_stat

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.dates as mdates
from matplotlib.colors import LinearSegmentedColormap

zarr_path = "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_thetao_sjtu_jingwei_v1_1.zarr"
time_range = ["1960-01-01", "2020-12-31"]
temp_data = load_regular_data(
    zarr_path=zarr_path,
    time_range=time_range,
    depth_range=None,  # Will be handled by compute_enso_nino_index
    lat_range=None,
    lon_range=None
)

print(f"Data dimensions: {temp_data.dims}")
print(f"Data shape: {temp_data.shape}")
print(f"Time range: {temp_data.time.values[0]} to {temp_data.time.values[-1]}")

nino34_index = compute_enso_nino_index(
    data=temp_data,
    nino_region="nino3.4",
    depth_level=10.0,  # Average 0-10m for SST
    anomaly_freq="month"
)
nino34_index = nino34_index - 273.15

print(f"\nNiño 3.4 Index calculated")
print(f"Data dimensions: {nino34_index.dims}")
print(f"Data shape: {nino34_index.shape}")
print(f"Index range: {float(nino34_index.min()):.2f} to {float(nino34_index.max()):.2f} °C")


# Standard thresholds:
# El Niño: Niño 3.4 index >= 0.5°C for 5 consecutive months
# La Niña: Niño 3.4 index <= -0.5°C for 5 consecutive months
el_nino_threshold = 0.5
la_nina_threshold = -0.5

# Get time coordinates and index values
times = nino34_index.time.values
values = nino34_index.values

# Create masks for potential events (simplified: based on threshold crossing)
el_nino_mask = values >= el_nino_threshold
la_nina_mask = values <= la_nina_threshold


width_cm = 9
height_cm = 5
width_inch = width_cm / 2.54
height_inch = height_cm / 2.54
plt.rcParams['font.family'] = 'DejaVu Serif'
fig, ax = plt.subplots(figsize=(width_inch*2.2, height_inch*2.2))

time_numeric = np.arange(len(times))

# Create filled areas for El Niño (red) and La Niña (blue) events
# Fill between for positive anomalies (El Niño)
ax.fill_between(
    times, 
    values, 
    el_nino_threshold,
    where=(values >= el_nino_threshold),
    color='#d6604d',
    alpha=0.6,
    label='El Niño (≥0.5°C)',
    interpolate=True
)

# Fill between for negative anomalies (La Niña)
ax.fill_between(
    times, 
    values, 
    la_nina_threshold,
    where=(values <= la_nina_threshold),
    color='#4393c3',
    alpha=0.6,
    label='La Niña (≤-0.5°C)',
    interpolate=True
)

# Plot the main line
ax.plot(times, values, color='#2C3E50', linewidth=1.2, zorder=5)

# Add threshold lines
ax.axhline(y=el_nino_threshold, color='#d6604d', linestyle='--', linewidth=1, alpha=0.7)
ax.axhline(y=la_nina_threshold, color='#4393c3', linestyle='--', linewidth=1, alpha=0.7)
ax.axhline(y=0, color='black', linestyle='-', linewidth=0.8, alpha=0.5)

# Add threshold labels
ax.text(times[-1], el_nino_threshold + 0.5, '+0.5°C', fontsize=9, color='#d6604d', 
        ha='right', va='bottom')
ax.text(times[-1], la_nina_threshold - 0.5, '-0.5°C', fontsize=9, color='#4393c3', 
        ha='right', va='top')

ax.set_xlabel('Year', fontsize=14)
ax.set_ylabel('Niño 3.4 Index (°C)', fontsize=14)
ax.set_title('Niño 3.4 Index (1960-2020)', fontsize=16, pad=15)

y_min = min(values.min(), -4)
y_max = max(values.max(), 4)
ax.set_ylim(y_min - 0.1, y_max + 0.1)

ax.xaxis.set_major_locator(mdates.YearLocator(10))
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax.xaxis.set_minor_locator(mdates.YearLocator(1))

ax.grid(True, which='major', alpha=0.3)

legend = ax.legend(loc='upper left', frameon=True, fontsize=13)
legend.get_frame().set_alpha(0.9)

ax.tick_params(axis='x', labelsize=13)
ax.tick_params(axis='y', labelsize=13)

# Add annotation for major events
# Find significant El Niño events (index > 1.5)
major_el_nino_idx = np.where(values > 1.5)[0]
if len(major_el_nino_idx) > 0:
    # Group consecutive indices to find event peaks
    peaks = []
    start = major_el_nino_idx[0]
    for i in range(1, len(major_el_nino_idx)):
        if major_el_nino_idx[i] - major_el_nino_idx[i-1] > 6:  # New event
            # Find peak in previous event
            event_idx = major_el_nino_idx[np.where(major_el_nino_idx < major_el_nino_idx[i])[0]]
            event_idx = event_idx[event_idx >= start]
            if len(event_idx) > 0:
                peak_idx = event_idx[np.argmax(values[event_idx])]
                peaks.append(peak_idx)
            start = major_el_nino_idx[i]
    # Don't forget the last event
    event_idx = major_el_nino_idx[major_el_nino_idx >= start]
    if len(event_idx) > 0:
        peak_idx = event_idx[np.argmax(values[event_idx])]
        peaks.append(peak_idx)
    
    # Annotate top 3 El Niño events
    peak_values = [(p, values[p]) for p in peaks]
    peak_values.sort(key=lambda x: x[1], reverse=True)
    
    for i, (peak_idx, peak_val) in enumerate(peak_values[:3]):
        year = str(times[peak_idx])[:4]
        ax.annotate(
            f'{year}',
            xy=(times[peak_idx], peak_val),
            xytext=(0, 15),
            textcoords='offset points',
            fontsize=10,
            color='#b2182b',
            ha='center',
            fontweight='bold'
        )

# Find significant La Niña events (index < -1.5)
major_la_nina_idx = np.where(values < -1.5)[0]
if len(major_la_nina_idx) > 0:
    peaks = []
    start = major_la_nina_idx[0]
    for i in range(1, len(major_la_nina_idx)):
        if major_la_nina_idx[i] - major_la_nina_idx[i-1] > 6:
            event_idx = major_la_nina_idx[np.where(major_la_nina_idx < major_la_nina_idx[i])[0]]
            event_idx = event_idx[event_idx >= start]
            if len(event_idx) > 0:
                peak_idx = event_idx[np.argmin(values[event_idx])]
                peaks.append(peak_idx)
            start = major_la_nina_idx[i]
    event_idx = major_la_nina_idx[major_la_nina_idx >= start]
    if len(event_idx) > 0:
        peak_idx = event_idx[np.argmin(values[event_idx])]
        peaks.append(peak_idx)
    
    # Annotate top 3 La Niña events
    peak_values = [(p, values[p]) for p in peaks]
    peak_values.sort(key=lambda x: x[1])
    
    for i, (peak_idx, peak_val) in enumerate(peak_values[:3]):
        year = str(times[peak_idx])[:4]
        ax.annotate(
            f'{year}',
            xy=(times[peak_idx], peak_val),
            xytext=(0, -15),
            textcoords='offset points',
            fontsize=10,
            color='#2166ac',
            ha='center',
            fontweight='bold'
        )

plt.tight_layout()

output_path = "/home/lorraine/jingwei_QA/src/applications/res/6_nino34_index_timeseries.png"
plt.savefig(output_path, dpi=450, facecolor='white')