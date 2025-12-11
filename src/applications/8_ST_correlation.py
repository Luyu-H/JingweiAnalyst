from src.tools.data_loader import load_regular_data
from src.tools.basic_calculator import compute_basic_stat, analyze_correlation

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import TwoSlopeNorm

temp_zarr_path = "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_temp_sjtu_jingwei_v1_1.zarr"
sal_zarr_path = "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_so_sjtu_jingwei_v1_1.zarr"

# time_range = ["1980-01-01", "2000-12-31"]
time_range = ["1980-01-01", "1980-01-31"]

temp_data = load_regular_data(
    zarr_path=temp_zarr_path,
    time_range=time_range,
    depth_range=[0, 50],
    lat_range=None,
    lon_range=None
)

sal_data = load_regular_data(
    zarr_path=sal_zarr_path,
    time_range=time_range,
    depth_range=[0, 50],
    lat_range=None,
    lon_range=None
)

print(f"Temperature data dimensions: {temp_data.dims}")
print(f"Temperature data shape: {temp_data.shape}")
print(f"Salinity data dimensions: {sal_data.dims}")
print(f"Salinity data shape: {sal_data.shape}")

sst = compute_basic_stat(data=temp_data, op="mean", dim="depth")
sss = compute_basic_stat(data=sal_data, op="mean", dim="depth")

print(f"\nSST data dimensions: {sst.dims}")
print(f"SST data shape: {sst.shape}")
print(f"SSS data dimensions: {sss.dims}")
print(f"SSS data shape: {sss.shape}")
print(f"Time range: {sst.time.values[0]} to {sst.time.values[-1]}")

correlation_result = analyze_correlation(
    var1=sst,
    var2=sss,
    dim="time",
    method="pearson",
    lag=0,
    min_periods=12  # Require at least 12 months of data
)

corr = correlation_result["correlation"]
pvalue = correlation_result["pvalue"]
n_samples = correlation_result["n_samples"]

print(f"\nCorrelation data dimensions: {corr.dims}")
print(f"Correlation data shape: {corr.shape}")
print(f"Correlation range: {float(corr.min()):.3f} to {float(corr.max()):.3f}")


fig = plt.figure(figsize=(14, 8))
ax = plt.axes(projection=ccrs.Robinson(central_longitude=180))

ax.set_global()
ax.add_feature(cfeature.LAND, facecolor='lightgray', edgecolor='black', linewidth=0.5)
ax.add_feature(cfeature.COASTLINE, linewidth=0.5)

lon = corr.lon.values
lat = corr.lat.values
data_plot = corr.values

# Use diverging colormap (blue-white-red)
# Blue for negative correlation, red for positive correlation
cmap = plt.cm.RdBu_r  # Reversed RdBu colormap (red for positive)

# Use TwoSlopeNorm to ensure zero value corresponds to white
vmin = -1
vmax = 1
norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

# Plot filled contours
mesh = ax.pcolormesh(
    lon, lat, data_plot,
    transform=ccrs.PlateCarree(),
    cmap=cmap,
    norm=norm,
    shading='auto'
)

# Add colorbar
cbar = plt.colorbar(
    mesh, 
    ax=ax, 
    orientation='horizontal',
    pad=0.05,
    shrink=0.8,
    aspect=40,
    extend='neither'
)
cbar.set_label('Correlation Coefficient', fontsize=12)
cbar.ax.tick_params(labelsize=10)

# Set colorbar ticks
cbar.set_ticks([-1, -0.5, 0, 0.5, 1])

# Add title
ax.set_title(
    'SST-SSS Correlation',
    fontsize=14,
    fontweight='bold',
    pad=10
)

# Add gridlines
gl = ax.gridlines(
    draw_labels=False,
    linewidth=0.5,
    color='gray',
    alpha=0.5,
    linestyle='--'
)

plt.tight_layout()

output_path = "/home/lorraine/jingwei_QA/src/applications/res/8_sst_sss_correlation_map.png"
plt.savefig(output_path, dpi=450, facecolor='white')