**Methodology used.** Monthly CMIP6 historical `tos` was loaded for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR over **1971-01 to 2000-12**, converted to **°C**, and averaged in time with `temporal_aggregate(freq='mean')` to form each model’s 1971–2000 SST climatology. Native-grid global means were computed with `compute_basic_stat(op='area_weighted_mean', dim=['lat','lon'])`; each model’s `areacello` was loaded and summed over valid climatology cells to record the total area. WOA23 `temperature` was loaded as the **annual 1971–2000 climate normal at 0 m depth** on the 1° grid, converted to °C, and its whole-grid global mean was computed analogously. WOA23 and each model climatology were compared after conservative horizontal regridding of the models to the WOA23 1° target (`regrid_horizontal(method='conservative')`). A single shared comparison mask was built with `build_comparison_mask(..., ocean_only=True, extra_mask=landsea_mask>1)` over all three regridded model fields plus WOA23 SST and applied to every field. This gave **40,830 comparison cells**. Metrics used `compute_difference_metrics(dim=None, weighting='auto')` and area-weighted Pearson correlation with `analyze_correlation(method='pearson', weighting='area', dim=['lat','lon'])`. All values below are °C unless noted.

**Native-grid global mean SST (Figure 2 quantities).**

| Field | Global mean SST (°C) | Valid area used (m²) |
|---|---:|---:|
| CESM2 | 18.52 | \(3.6051\times10^{14}\) |
| GFDL-ESM4 | 18.13 | \(3.6111\times10^{14}\) |
| MPI-ESM1-2-LR | 17.58 | \(3.6598\times10^{14}\) |
| WOA23, whole 1° grid | 18.19 | \(3.4933\times10^{14}\) |

These are native-grid/whole-WOA23 means and are **not** computed over the shared comparison cells.

**Shared-sample comparison with WOA23 on the 1° grid.** The WOA23 mean returned by `compute_difference_metrics` on the shared cells is **18.20°C**; this differs slightly from the whole-grid WOA23 mean of 18.19°C because the comparison is restricted to cells common to all fields.

| Model | Shared sample mean SST (°C) | WOA23 shared mean (°C) | Bias, model − WOA23 (°C) | RMSE (°C) | Area-weighted Pearson r | Largest single departure |
|---|---:|---:|---:|---:|---:|---|
| CESM2 | 18.48 | 18.20 | **+0.28** | **1.062** | **0.995** | 8.84°C at 43.5°N, 64.5°W |
| GFDL-ESM4 | 18.12 | 18.20 | **−0.09** | **0.832** | **0.997** | 7.38°C at 16.5°S, 11.5°E |
| MPI-ESM1-2-LR | 17.57 | 18.20 | **−0.63** | **1.339** | **0.993** | 7.61°C at 16.5°S, 11.5°E |

**Ranking by absolute mean bias, smallest to largest:**  
**GFDL-ESM4 (0.09°C) < CESM2 (0.28°C) < MPI-ESM1-2-LR (0.63°C)**. The ranking was determined from unrounded absolute biases.

**Spatial structure of the bias.** All three models capture the large-scale SST pattern very well, with area-weighted spatial correlations of 0.993–0.997. The biases are regional rather than uniform: each model has both warm and cold local departures. CESM2 has a positive shared-domain mean bias and a large positive local maximum near 43.5°N, 64.5°W; its difference field ranges from about −4.80 to +8.84°C. GFDL-ESM4 is closest to the WOA23 shared-cell mean, slightly cold overall, with local differences from about −4.49 to +7.38°C and its largest absolute departure near 16.5°S, 11.5°E. MPI-ESM1-2-LR is the coldest of the three on the comparison cells, with a stronger negative mean bias and local differences from about −7.33 to +7.61°C; its largest absolute departure is also near 16.5°S, 11.5°E. Thus, although the large-scale patterns are highly correlated with WOA23, the residual bias is spatially structured, with the largest local errors concentrated in specific regions rather than being a simple global offset.

**Visualizations produced.**  
- **Figure 1** is a 2×2 global map panel: the first three panels show CESM2, GFDL-ESM4 and MPI-ESM1-2-LR SST biases (`model − WOA23`) on the shared comparison cells, all using one shared diverging colour scale centred on zero and clipped at the 98th percentile of absolute bias across the three models; the fourth panel shows the WOA23 0 m reference SST on the same shared cells with a sequential colour scale (2nd–98th percentile limits). Every panel is labelled, and the colour bars are labelled in °C.  
- **Figure 2** is a bar chart of the four native-grid/whole-grid global mean SST values (CESM2, GFDL-ESM4, MPI-ESM1-2-LR and WOA23), each bar annotated with its value, with the y-axis labelled in °C. Its title notes that each bar is over its own native ocean grid and not over the shared comparison cells.

**Key insight.** GFDL-ESM4 performs best against WOA23 by mean bias and RMSE in this comparison, CESM2 is intermediate with a warm shared-sample bias, and MPI-ESM1-2-LR has the largest cold shared-sample bias and highest RMSE, while all models retain very high spatial pattern correlation with WOA23.
