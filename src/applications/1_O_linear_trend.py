from src.tools.data_loader import load_regular_data
from src.tools.basic_calculator import compute_basic_stat, fit_1D_trend, temporal_aggregate
from src.tools.scientific_calculator import convert_concentration_to_content

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

zarr_path = "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_o2_sjtu_jingwei_v1_1.zarr"
time_range = ["2000-01-01", "2020-12-31"]
o2_concentration = load_regular_data(
    zarr_path=zarr_path,
    time_range=time_range,
    depth_range=None,
    lat_range=None,
    lon_range=None
)
print(f"Dimensions: {o2_concentration.dims}")
print(f"Shapes: {o2_concentration.shape}")

o2_content = convert_concentration_to_content(
    concentration=o2_concentration,
    density=1025.0
)
print(f"Shape of the content data: {o2_content.shape}")

global_o2_timeseries = compute_basic_stat(
    data=o2_content,
    op="sum",
    dim=["depth", "lat", "lon"]
)
print(f"Length of the time series: {len(global_o2_timeseries)}")
print(f"Time range: {global_o2_timeseries.time.values[0]} - {global_o2_timeseries.time.values[-1]}")

y = global_o2_timeseries.values.astype(float)
time_values = pd.to_datetime(global_o2_timeseries.time.values)

# 使用temporal_aggregate计算气候态月平均
monthly_climatology = temporal_aggregate(
    data=global_o2_timeseries,
    freq='climatology',
    op='mean'
)
print(f"Monthly climatology shape: {monthly_climatology.shape}")

# 构建季节性分量：将气候态映射回原始时间序列
months = time_values.month
seasonal_component = np.array([monthly_climatology.sel(month=m).values for m in months])

# 去除季节性（anomaly = 原始值 - 该月气候态）
y_deseasonalized = y - seasonal_component

# 对去季节性后的数据进行标准化
y_deseas_mean = y_deseasonalized.mean()
y_deseas_std = y_deseasonalized.std()

# 如果去季节性后标准差很小，使用原始数据的标准差来避免数值问题
if y_deseas_std < 1e-10:
    y_deseas_std = y.std()

y_norm = (y_deseasonalized - y_deseas_mean) / y_deseas_std

global_o2_timeseries_norm = xr.DataArray(
    y_norm,
    coords=global_o2_timeseries.coords,
    dims=global_o2_timeseries.dims,
    attrs=global_o2_timeseries.attrs
)

trend_result = fit_1D_trend(
    data=global_o2_timeseries_norm,
    method="linear, a * x + b"
)

fitted_anomaly = trend_result["fitted_data"].values * y_deseas_std + y_deseas_mean
# 将趋势直线平移到原始数据的均值位置（不加回季节性，保持直线）
fitted_line = fitted_anomaly + y.mean() - y_deseas_mean
fitted_data = xr.DataArray(
    fitted_line,
    coords=global_o2_timeseries.coords,
    dims=global_o2_timeseries.dims,
    attrs={**global_o2_timeseries.attrs, "fitted_by": "linear_trend_line"}
)
equation = trend_result["equation"]
metrics = trend_result["metrics"]

y_anomaly_centered = y_deseasonalized - y_deseas_mean + y.mean()
ss_res = np.sum((y_anomaly_centered - fitted_line) ** 2)
ss_tot = np.sum((y_anomaly_centered - y_anomaly_centered.mean()) ** 2)
r2_trend = 1 - (ss_res / ss_tot)

print(f"\nFitting Result:")
print(f"Equation: {equation}")
print(f"R² (deseasonalized): {metrics['r2']:.4f}")
print(f"R² (original scale): {r2_trend:.4f}")
print(f"RMSE: {metrics['rmse']:.4e}")
print(f"MAE: {metrics['mae']:.4e}")


# Visualization
time_values = pd.to_datetime(global_o2_timeseries.time.values)

width_cm = 9
height_cm = 5
width_inch = width_cm / 2.54
height_inch = height_cm / 2.54
plt.rcParams['font.family'] = 'DejaVu Serif'
fig, ax = plt.subplots(figsize=(width_inch*2.2, height_inch*2.2))

n_months = len(time_values)

ax.plot(time_values, global_o2_timeseries.values, 
        '#2171b5', linewidth=1.5, alpha=0.7, label='Global O₂ Content')

ax.plot(time_values, fitted_data.values, 
        'r-', linewidth=2, label=f'Linear Trend (R²={metrics["r2"]:.3f})')

ax.tick_params(axis='x', labelsize=13)
ax.tick_params(axis='y', labelsize=13)
ax.set_xlabel('Year', fontsize=14)
ax.set_ylabel('Global O₂ Content (mol)', fontsize=14)
ax.set_title('Global Ocean DO Content Trend (2000-2020)', fontsize=16)

ax.xaxis.set_major_locator(mdates.YearLocator(5))  
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax.xaxis.set_minor_locator(mdates.YearLocator(1))

ax.grid(True, which='major', alpha=0.3)
ax.grid(True, which='minor', alpha=0.1)
ax.legend(loc='upper right', fontsize=14)

plt.tight_layout()

output_path = "/home/lorraine/jingwei_QA/src/applications/res/1_global_o2_trend_2000_2020.png"
plt.savefig(output_path, dpi=450)