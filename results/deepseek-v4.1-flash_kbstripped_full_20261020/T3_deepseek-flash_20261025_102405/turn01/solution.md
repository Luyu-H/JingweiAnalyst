The 1971–2000 CMIP6 historical SST climatologies reproduce the broad spatial pattern of the WOA23 0 m SST climatology very well, but their global-mean biases differ. On the common WOA23 1° comparison mask, **GFDL-ESM4 is closest to WOA23**, followed by **CESM2**, then **MPI-ESM1-2-LR**. CESM2 has a warm global-mean bias, while GFDL-ESM4 and MPI-ESM1-2-LR are cold-biased relative to WOA23.

**Methodology actually used**
- CMIP6 `tos` from `cmip6_historical` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, monthly from 1971-01 to 2000-12, loaded with `load_regular_data`.
- Values were kept only where finite and `|tos| < 1e10`; converted from the loaded unit to `degC` with `transform_unit`.
- A 1971–2000 climatology was formed with `temporal_aggregate(freq="mean", op="mean")`, followed by finite-value and physical-range filtering `(-5, 45) °C`.
- Native model `areacello` fields were loaded with `load_areas`; the native-grid mean SST was computed with `compute_basic_stat(op="area_weighted_mean", dim=["lat", "lon"])`. The reported native ocean area is the sum of `areacello` over valid SST cells.
- WOA23 temperature was loaded with `load_climatology_data(dataset="woa23", variable="temperature", time_slice="annual", period="1971-2000", depth_range=[0,0])`. Ocean cells were defined from `landsea_01 > 1`. WOA23 areas were loaded with `load_areas(dataset="woa23", variable="temperature")`.
- Each model climatology was regridded to the WOA23 1° grid with `regrid_horizontal(..., target="woa23", method="conservative")`.
- A single common comparison mask was built over WOA23 and all three regridded model fields using `build_comparison_mask(ocean_only=True)`, then applied with `.where(mask)` to every field before all shared-cell metrics.
- On the shared cells, means were computed with `compute_basic_stat(op="area_weighted_mean")`; bias and RMSE with `compute_difference_metrics(..., weighting="auto")`; and spatial Pearson correlation with `analyze_correlation(..., weighting="area")`.
- No lag analysis was performed; all comparisons are contemporaneous 1971–2000 climatology comparisons.

**Native-grid 1971–2000 global mean SST**
| Product | Native-grid mean SST (°C) | Approx. area covered (m²) |
|---|---:|---:|
| CESM2 | **18.52** | 3.605×10¹⁴ |
| GFDL-ESM4 | **18.13** | 3.611×10¹⁴ |
| MPI-ESM1-2-LR | **17.58** | 3.660×10¹⁴ |
| WOA23 | **18.19** | 3.493×10¹⁴ |

**Shared WOA23 1° comparison cells**
- Number of common comparison cells: **40,830**.
- Shared-cell mean SST: CESM2 **18.48 °C**, GFDL-ESM4 **18.12 °C**, MPI-ESM1-2-LR **17.57 °C**, WOA23 **18.20 °C**.

| Model | Bias, model − WOA23 (°C) | RMSE (°C) | Spatial Pearson r |
|---|---:|---:|---:|
| CESM2 | **+0.28** | **1.062** | **0.995** |
| GFDL-ESM4 | **−0.09** | **0.832** | **0.997** |
| MPI-ESM1-2-LR | **−0.63** | **1.339** | **0.993** |

**Ranking by smallest to largest absolute shared-cell bias**
1. **GFDL-ESM4** (|−0.086| ≈ 0.09 °C)
2. **CESM2** (0.280 °C)
3. **MPI-ESM1-2-LR** (0.631 °C)

**Spatial structure of the bias**
All three models capture the large-scale SST pattern strongly, with spatial correlations above 0.993. The remaining errors are regional rather than a uniform global offset: the bias arrays span roughly −4.8 to +8.8 °C for CESM2, −4.5 to +7.4 °C for GFDL-ESM4, and −7.3 to +7.6 °C for MPI-ESM1-2-LR on the shared cells. CESM2 has a net warm shared-cell bias and larger positive local excursions; GFDL-ESM4 has the smallest absolute mean bias and lowest RMSE; MPI-ESM1-2-LR has the strongest cold mean bias and largest RMSE and localized negative excursions. Thus, the global ranking is controlled by small mean offsets, while the spatial maps show substantial but highly correlated regional structure.

**Visualizations**
- **Figure 1** is a 2×2 global map: the first three panels show the shared-cell SST bias (`model − WOA23`) for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR using one common diverging `RdBu_r` colour scale centred on zero, with limits set from the 98th percentile of absolute bias. The fourth panel shows WOA23 0 m SST on the same shared cells with a sequential `viridis` colour scale. All panels are labelled and the colour bars are labelled in degrees Celsius.
- **Figure 2** is a bar chart of the four native-grid global mean SST values: CESM2, GFDL-ESM4, MPI-ESM1-2-LR and WOA23. Each bar is annotated with its value (and the code also annotates the corresponding area), and the y-axis is labelled in degrees Celsius.

**Caveats**
These are single-member coupled historical simulations, so the 1971–2000 climatological biases and ranking mix forced model error with internal variability, and event timing is not comparable between models and observations. Also, CMIP6 `tos` is a top-cell bulk temperature including under-ice values rather than a skin temperature; differing top-cell thicknesses and grid resolutions can contribute to inter-model differences, especially for narrow or polar features.
