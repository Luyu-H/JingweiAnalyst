**Methodology (as actually executed).** Monthly CMIP6 historical `tos` for CESM2, GFDL-ESM4, and MPI-ESM1-2-LR was loaded for 1971-01 to 2000-12. Where the loaded units were kelvin, the fields were converted to degrees Celsius with `transform_unit`; the 1971–2000 climatology was then formed with `temporal_aggregate(data=tos_c, freq="mean", op="mean")`. Each model’s climatology was kept on its original grid for the native-grid mean and area accounting. Native-grid means were computed with `compute_basic_stat(..., op="area_weighted_mean", dim=None)`; `areacello` was loaded and used to report the total area covered by valid climatology cells. WOA23 `temperature` was loaded as the 1971–2000 annual climatology at 0 m on the 1° grid, and WOA23 `landsea_01 > 1` was used as the ocean mask for its native-grid mean. The WOA23 mean is therefore an area-weighted mean over valid ocean cells, not an unweighted average over all 1° grid cells.

For the model–WOA comparison, each model climatology was conservatively regridded to the WOA23 1° grid using `regrid_horizontal(..., target="woa23", method="conservative")`. A common comparison mask was built over the three regridded model fields plus WOA23 using `build_comparison_mask(..., ocean_only=True)`; it contains **40,830 cells**. On those common cells, `compute_difference_metrics(..., weighting="auto")` was used for shared-cell means, bias (`delta_SST`), and RMSE, and `analyze_correlation(..., method="pearson", weighting="area")` was used for spatial Pearson correlation.

**Native-grid global mean SST, 1971–2000**

| Field | Area-weighted mean SST (°C) | Total ocean area covered (m²) |
|---|---:|---:|
| CESM2 | 18.52 | 3.6051 × 10¹⁴ |
| GFDL-ESM4 | 18.13 | 3.6111 × 10¹⁴ |
| MPI-ESM1-2-LR | 17.58 | 3.6598 × 10¹⁴ |
| WOA23 | 18.19 | 3.4933 × 10¹⁴ |

The WOA23 native-grid value, **18.19 °C**, is the area-weighted mean over valid WOA23 ocean cells with `landsea_01 > 1`. GFDL-ESM4 was processed on its 180×360 lat–lon archive grid, so its “native-grid” result refers to that archive grid rather than to a different-resolution source grid.

**Common-cell comparison on the WOA23 1° grid**

| Model | Model mean SST (°C) | WOA23 mean SST (°C) | `delta_SST` model − WOA23 (°C) | RMSE (°C) | Spatial Pearson correlation |
|---|---:|---:|---:|---:|---:|
| CESM2 | 18.48 | 18.20 | +0.28 | 1.062 | 0.995 |
| GFDL-ESM4 | 18.12 | 18.20 | −0.09 | 0.832 | 0.997 |
| MPI-ESM1-2-LR | 17.57 | 18.20 | −0.63 | 1.339 | 0.993 |

Ranking by smallest to largest absolute mean bias:

**GFDL-ESM4 (0.09 °C) < CESM2 (0.28 °C) < MPI-ESM1-2-LR (0.63 °C).**

**Interpretation.** GFDL-ESM4 has the smallest absolute mean bias and the lowest RMSE on the common comparison cells. CESM2 is slightly warm on average (+0.28 °C), while MPI-ESM1-2-LR is cold on average (−0.63 °C) and has the largest RMSE. All three models have very high spatial correlations with WOA23 (0.993–0.997), indicating that the large-scale SST pattern is well captured. The bias is not simply a uniform offset: the bias-field ranges and standard deviations are larger than the mean biases, with CESM2 ranging from about −4.80 to +8.84 °C, GFDL-ESM4 from about −4.49 to +7.38 °C, and MPI-ESM1-2-LR from about −7.33 to +7.61 °C. Thus the spatial structure consists of a correct broad-scale SST gradient with regionally varying warm/cold model biases; GFDL-ESM4 has the smallest bias amplitude, while MPI-ESM1-2-LR has the largest regional bias spread.

**Figures.** Figure 1 is a 2×2 map panel showing the CESM2, GFDL-ESM4, and MPI-ESM1-2-LR SST bias fields (model minus WOA23) on the common comparison cells, using one shared diverging colour scale centred on zero and clipped at the pooled 98th percentile of absolute bias. The fourth panel shows the WOA23 0 m SST reference on the same cells with a sequential colour scale. All panels are labelled and the colour bars carry °C units. Figure 2 is a bar chart of the four native-grid global mean SST values (CESM2, GFDL-ESM4, MPI-ESM1-2-LR, WOA23), with each bar annotated by its value and the y-axis labelled in °C. The native-grid means and the common-cell means are distinct quantities and are labelled separately.
