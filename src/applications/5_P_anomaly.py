from src.tools.data_loader import load_regular_data
from src.tools.basic_calculator import temporal_aggregate, compute_basic_stat
from src.tools.scientific_calculator import compute_anomaly

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import matplotlib.dates as mdates

zarr_path = "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_po4_sjtu_jingwei_v1_1.zarr"

time_range = ["1995-01-01", "2000-12-31"]
po4_data = load_regular_data(
    zarr_path=zarr_path,
    time_range=time_range,
    depth_range=[0, 100],
    lat_range=[-5, 5],
    lon_range=[120, 290]
)

print(f"Original dimensions: {po4_data.dims}")
print(f"Original data shape: {po4_data.shape}")
print(f"Time range: {po4_data.time.values[0]} - {po4_data.time.values[-1]}")

po4_anomaly = compute_anomaly(
    data=po4_data,
    var="po4",
    freq="season",
    method="subtract"
)

print(f"Dimensions of the anomaly data: {po4_anomaly.dims}")
print(f"Shape of the anomaly data: {po4_anomaly.shape}")

po4_anomaly = compute_basic_stat(
    data=po4_anomaly,
    op="mean",
    dim="depth"
)

po4_anomaly_seasonal = temporal_aggregate(
    data=po4_anomaly,
    freq="season",
    op="mean"
)

print(f"Dimensions after seasonal aggregation: {po4_anomaly_seasonal.dims}")
print(f"Data shape after seasonal aggregation: {po4_anomaly_seasonal.shape}")

po4_hovmoller = compute_basic_stat(
    data=po4_anomaly_seasonal,
    op="mean",
    dim="lat"
)

print(f"Dimensions of Hovmöller data: {po4_hovmoller.dims}")
print(f"Shape of Hovmöller data: {po4_hovmoller.shape}")


width_cm = 9
height_cm = 5
width_inch = width_cm / 2.54
height_inch = height_cm / 2.54
plt.rcParams['font.family'] = 'DejaVu Serif'
fig, ax = plt.subplots(figsize=(width_inch*2.2, height_inch*2.2))

lon = po4_hovmoller.lon.values
time = po4_hovmoller.time.values
data_plot = po4_hovmoller.values.T

# 将 -180~180 格式的经度转换为连续的绘图坐标
# 负值经度加 360 转换为 0~360 格式用于绘图
lon_plot = np.where(lon < 0, lon + 360, lon)

LON, TIME = np.meshgrid(lon_plot, range(len(time)))

vmax = np.nanpercentile(np.abs(data_plot), 95)
vmin = -vmax

cmap = plt.cm.RdBu_r

norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

mesh = ax.pcolormesh(
    LON, TIME, data_plot,
    cmap=cmap,
    norm=norm,
    shading='auto'
)

cbar = plt.colorbar(mesh, ax=ax, orientation='vertical', pad=0.02, shrink=1)
cbar.set_label('Phosphate Anomaly (μmol/kg)', fontsize=13)
cbar.ax.tick_params(labelsize=13)

ax.set_xlabel('Longitude', fontsize=14)
ax.set_xlim([lon_plot.min(), 280])

# 设置经度刻度和标签，调整范围到280
lon_ticks = np.arange(120, 281, 20)  # 包含280
lon_labels = []
for l in lon_ticks:
    if l <= 180:
        lon_labels.append(f'{int(l)}°E')
    else:
        lon_labels.append(f'{int(360-l)}°W')
ax.set_xticks(lon_ticks)
ax.set_xticklabels(lon_labels)

ax.set_ylabel('Time', fontsize=14)

n_times = len(time)

if hasattr(po4_hovmoller, 'year') and hasattr(po4_hovmoller, 'season'):
    years = po4_hovmoller.year.values
    seasons = po4_hovmoller.season.values
    
    tick_positions = []
    tick_labels = []
    for i in range(0, n_times, 4):  # 每4个季节（一年）显示一个标签
        if i < n_times:
            tick_positions.append(i)
            tick_labels.append(str(int(years[i])))
    
    ax.set_yticks(tick_positions)
    ax.set_yticklabels(tick_labels)
else:
    n_ticks = min(15, n_times)
    tick_step = max(1, n_times // n_ticks)
    tick_positions = range(0, n_times, tick_step)
    ax.set_yticks(tick_positions)

ax.set_title(
    'Equatorial Pacific PO₄ Anomaly',
    fontsize=16
)

ax.tick_params(axis='x', labelsize=13)
ax.tick_params(axis='y', labelsize=13)

cs = ax.contour(
    LON, TIME, data_plot,
    levels=[0],
    colors='black',
    linewidths=1,
    linestyles='--'
)

ax.grid(True, which='major', alpha=0.3)

plt.tight_layout()

output_path = "/home/lorraine/jingwei_QA/src/applications/res/5_po4_anomaly_hovmoller.png"
plt.savefig(output_path, dpi=450, facecolor='white')