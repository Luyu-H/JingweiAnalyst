**Evaluation.** The three CMIP6 historical simulations (CESM2, GFDL-ESM4, MPI-ESM1-2-LR) were evaluated against the WOA23 1971–2000 annual 0 m temperature climatology. Monthly `tos` data were loaded from `cmip6_historical` for 1971-01 through 2000-12 on each model’s native grid, and the 1971–2000 climatology was formed with `temporal_aggregate(freq="mean", op="mean")`. Native-grid global means were computed with `compute_basic_stat(..., op="area_weighted_mean", dim=["lat","lon"])` on each model grid, without applying the common WOA23 mask. The WOA23 global mean was computed with the same area-weighted statistic over all valid WOA23 0 m cells; the `landsea_01` ocean mask was not applied for that global value. For the model–WOA comparison, each model climatology was conservatively regridded to the WOA23 1° grid using `regrid_horizontal(target="woa23", method="conservative")`. The common comparison mask required valid values in all three regridded models, WOA23, and `landsea_01 > 1`; this yielded 40,830 common 1° cells. Biases, RMSE and correlations were computed on those common cells with area weighting. Metrics therefore characterize the regridded fields, not the native-grid climatologies.

**Key numerical results.**

| Quantity | CESM2 | GFDL-ESM4 | MPI-ESM1-2-LR |
|---|---:|---:|---:|
| Native-grid global mean SST (°C) | 18.52 | 18.13 | 17.58 |
| Common-cell mean SST (°C) | 18.48 | 18.12 | 17.57 |
| Mean bias vs WOA23 (°C) | +0.28 | −0.09 | −0.63 |
| RMSE vs WOA23 (°C) | 1.062 | 0.832 | 1.339 |
| Spatial Pearson correlation | 0.995 | 0.997 | 0.993 |

The WOA23 global mean over all valid 1° cells is **18.19 °C**; on the common comparison cells its area-weighted mean is **18.20 °C**.

**Ranking by absolute mean bias (smallest to largest):**  
**GFDL-ESM4 (|0.09| °C) < CESM2 (|0.28| °C) < MPI-ESM1-2-LR (|0.63| °C).**

**Spatial structure and interpretation.** All three models capture the first-order global SST pattern, as shown by correlations of 0.993–0.997. The remaining bias is strongly regional rather than uniform. Large positive and negative departures occur at high latitudes, especially near the Southern Ocean sea-ice margin, and in eastern equatorial/coastal upwelling regions and western boundary current zones. The bias fields reach local extrema of roughly −4.8 to +8.8 °C for CESM2, −4.5 to +7.4 °C for GFDL-ESM4, and −7.3 to +7.6 °C for MPI-ESM1-2-LR. CESM2 has a net warm global-area bias with localized strong warm anomalies. GFDL-ESM4 has the smallest global-area bias and the lowest RMSE, indicating substantial cancellation between regional warm and cold errors. MPI-ESM1-2-LR has a widespread cold bias, the largest RMSE, and the lowest correlation, so it is the least consistent with WOA23 on the common comparison grid.

**Visualizations.** Figure 1 is a 2×2 global map panel. The first three panels show the model-minus-WOA23 SST bias fields on the common comparison cells, using a shared diverging colour scale centred on zero. The fourth panel shows the WOA23 0 m reference SST on the same cells with a sequential colour scale. All panels are labelled, and the colour bars are labelled in degrees Celsius. Figure 2 is a bar chart of the four global mean SST values: CESM2, GFDL-ESM4, MPI-ESM1-2-LR on their native ocean grids and WOA23 over all valid 1° cells. Each bar is annotated with its value, and the y-axis is labelled in degrees Celsius.
