from src.tools.data_loader import load_regular_data
from src.tools.basic_calculator import compute_basic_stat
from src.tools.scientific_calculator import identify_oxygen_minimum_zone

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

zarr_path = "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_o2_sjtu_jingwei_v1_1.zarr"
time_range = ["2010-01-01", "2020-12-31"]
# time_range = ["2010-01-01", "2010-03-31"]
o2_data = load_regular_data(
    zarr_path=zarr_path,
    time_range=time_range,
    depth_range=None,
    lat_range=None,
    lon_range=None
)
print(f"Dimensions: {o2_data.dims}")
print(f"Data shape: {o2_data.shape}")

omz_thickness, omz_mask = identify_oxygen_minimum_zone(
    o2_data=o2_data,
    threshold=20.0
)
print(f"Dimensions of the OMZ data: {omz_thickness.dims}")
print(f"Shape of the OMZ thickness data: {omz_thickness.shape}")

omz_thickness_mean = compute_basic_stat(
    data=omz_thickness,
    op="mean",
    dim="time"
)
print(f"Range of the OMZ thickness: {float(omz_thickness_mean.min()):.1f} - {float(omz_thickness_mean.max()):.1f} m")


width_cm = 9
height_cm = 5
width_inch = width_cm / 2.54
height_inch = height_cm / 2.54
plt.rcParams['font.family'] = 'DejaVu Serif'
fig = plt.figure(figsize=(width_inch*2.2, height_inch*2.2))

# 使用GridSpec控制布局，不使用tight_layout
gs = gridspec.GridSpec(
    1, 1,
    left=0.05, right=0.85,  # 留出右边空间放colorbar inset
    top=0.95, bottom=0.05
)

ax = fig.add_subplot(gs[0], projection=ccrs.Robinson(central_longitude=180))
cax = inset_axes(
    ax,
    width="3%",       # 颜色条宽度（相对 ax）
    height="100%",    # 必须写100%，才能和地图等高
    loc="lower left",
    bbox_to_anchor=(1.05, 0.0, 1, 1),  # x偏移1.02，使其放在右侧
    bbox_transform=ax.transAxes,
    borderpad=0
)

ax.set_global()
ax.set_aspect(1.15)
ax.add_feature(cfeature.LAND, facecolor='lightgray', edgecolor='black', linewidth=0.5)
ax.add_feature(cfeature.COASTLINE, linewidth=0.5)

lon = omz_thickness_mean.lon.values
lat = omz_thickness_mean.lat.values
data_plot = omz_thickness_mean.values

colors = [
    '#f7fbff',
    '#deebf7',
    '#c6dbef',
    '#9ecae1',
    '#6baed6',
    '#4292c6',
    '#2171b5',
    '#08519c',
    '#08306b',
    '#4a1486',
    '#6a0d83',
]
cmap = LinearSegmentedColormap.from_list('omz_cmap', colors, N=256)

vmin = 0
vmax = np.nanpercentile(data_plot, 99)

mesh = ax.pcolormesh(
    lon, lat, data_plot,
    transform=ccrs.PlateCarree(),
    cmap=cmap,
    vmin=vmin,
    vmax=vmax,
    shading='auto'
)

# Colorbar使用独立的axes
cbar = plt.colorbar(mesh, cax=cax, orientation='vertical')
cbar.set_label('OMZ Thickness (m)', fontsize=13)
cbar.ax.tick_params(labelsize=13)

# 标题不加粗
ax.set_title(
    'Global OMZ Thickness Distribution',
    fontsize=16,
    pad=10
)

# 添加经纬度网格线和标注
gl = ax.gridlines(
    crs=ccrs.PlateCarree(),
    draw_labels=True,
    linewidth=0.5,
    color='gray',
    alpha=0.5,
    linestyle='--',
    x_inline=False,
    y_inline=False
)

# 设置经纬度标注位置
gl.top_labels = False
gl.right_labels = False
gl.left_labels = True
gl.bottom_labels = True

# 设置经纬度间隔（经度更稀疏）
gl.xlocator = plt.FixedLocator([-180, -90, 0, 90, 180])
gl.ylocator = plt.FixedLocator([-60, -30, 0, 30, 60])

# 设置标注字体大小
gl.xlabel_style = {'size': 13}
gl.ylabel_style = {'size': 13}

output_path = "/home/lorraine/jingwei_QA/src/applications/res/2_omz_thickness_distribution.png"
plt.savefig(output_path, dpi=450, facecolor='white')
plt.show()
