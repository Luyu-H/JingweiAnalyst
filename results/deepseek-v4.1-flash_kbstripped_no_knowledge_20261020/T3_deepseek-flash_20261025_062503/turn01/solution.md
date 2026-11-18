**Approach used by the executed code.** The CMIP6 historical `tos` fields for CESM2, GFDL-ESM4, and MPI-ESM1-2-LR were loaded with `load_regular_data(dataset="cmip6_historical", variable="tos", time_range=["1971-01","2000-12"])`. The monthly records were reduced with `temporal_aggregate(freq="mean", op="mean")` to a 1971–2000 long-term climatology. Native-grid global means were computed with `compute_basic_stat(op="area_weighted_mean", dim=["lat","lon"])`. The WOA23 field was loaded with `load_climatology_data(dataset="woa23", variable="temperature", time_slice="annual", period="1971-2000", depth_range=[0.0,0.0])`; its full-grid mean was computed by the same area-weighted operation over valid WOA23 cells. Each model climatology was regridded conservatively to the WOA23 1° grid with `regrid_horizontal(target="woa23", method="conservative")`. The comparison mask was created by `build_comparison_mask(fields=[...], ocean_only=True)` on the WOA23 grid, giving **40,830 common comparison cells**. On those cells, model and WOA23 means, model-minus-WOA23 biases, RMSE, and Pearson correlations were computed with the corresponding tools. Note: the executed code loaded a WOA23 land-sea mask but did not pass it into the comparison-mask construction; the actual mask came from `build_comparison_mask(..., ocean_only=True)`. The code also requested area-weighted statistics through `compute_basic_stat`, rather than explicitly calling `load_areas` in the shown execution.

**Key numerical results.** Values are in °C; `δSST = model − WOA23`.

| Product | Native-grid global mean | Comparison-cell mean | δSST | RMSE | Spatial Pearson r |
|---|---:|---:|---:|---:|---:|
| CESM2 | **18.52** | **18.48** | **+0.28** | **1.062** | **0.996** |
| GFDL-ESM4 | **18.13** | **18.12** | **−0.09** | **0.832** | **0.997** |
| MPI-ESM1-2-LR | **17.58** | **17.57** | **−0.63** | **1.339** | **0.995** |
| WOA23 | **18.19** | **18.20** | — | — | — |

WOA23’s full 1° valid-grid global mean is **18.19 °C**; on the common comparison cells its mean is **18.20 °C**. The ranking from smallest to largest absolute mean bias is:

**GFDL-ESM4 (|−0.09|) < CESM2 (|+0.28|) < MPI-ESM1-2-LR (|−0.63|).**

**Interpretation.** All three models reproduce the large-scale SST pattern very well, with spatial correlations of 0.995–0.997. GFDL-ESM4 is clearly the closest to WOA23 overall: smallest absolute mean bias, lowest RMSE, and highest correlation, with a very slight cold bias. CESM2 has a modest warm global-mean bias (+0.28 °C) and intermediate RMSE. MPI-ESM1-2-LR has the strongest cold mean bias (−0.63 °C), the largest RMSE, and the lowest correlation.

The bias is not a uniform offset. The mapped bias fields have local extrema of several °C: CESM2 ranges from −4.80 to +8.84 °C with bias standard deviation 1.03 °C; GFDL-ESM4 from −4.49 to +7.38 °C with std 0.82 °C; MPI-ESM1-2-LR from −7.33 to +7.61 °C with std 1.20 °C. Thus the high correlations indicate a correct broad zonal SST pattern, while the residual biases are regionally varying. CESM2 is warm on average, GFDL-ESM4 is nearly unbiased/slightly cold, and MPI-ESM1-2-LR is cold on average; however, all models show both warm and cold local biases. The bias spread is smallest for GFDL-ESM4 and largest for MPI-ESM1-2-LR.

**Figures.** Figure 1 is a 2×2 panel on the common WOA23 comparison cells: the three model-minus-WOA23 bias maps use a shared zero-centred diverging colour scale and a common °C label, while the WOA23 0 m SST reference panel uses a sequential colour scale in °C. This figure is designed to show the sign, spatial coherence, and regional structure of each model’s bias. Figure 2 is a bar chart of the four global-mean SST values—CESM2, GFDL-ESM4, MPI-ESM1-2-LR, and WOA23—each computed over its own native grid as requested, annotated with the numerical value and with the y-axis labelled in °C.

**Caveat.** Native-grid global means are not directly cell-comparable across models, and the comparison statistics apply only to the 40,830 common cells after conservative regridding; the WOA23 comparison-cell mean therefore differs slightly from its full-grid mean.
