The three CMIP6 historical simulations reproduce the WOA23 1971–2000 SST climatology well in spatial pattern, with small global-mean offsets but non-uniform regional biases. Using the executed workflow, `tos` was loaded for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR over 1971-01 to 2000-12, converted from K to °C, and temporally aggregated with `temporal_aggregate(freq="mean")` to form the 1971–2000 climatology. WOA23 temperature was loaded as the 1971–2000 annual climatology, depth range `[0,10]`, then the first depth level was selected; it was treated as °C. The WOA23 land–sea mask was loaded with `load_landsea_mask` and ocean cells were defined by `landsea_01 > 1`. Model climatologies were conservatively regridded to the WOA23 1° grid with `regrid_horizontal(target="woa23", method="conservative")`. A common comparison mask was built with `build_comparison_mask(..., ocean_only=True, extra_mask=ocean)`, restricted to cells valid for all three regridded models and WOA23; 40,830 common cells were used. Means, biases and RMSE were computed with `compute_difference_metrics`, and spatial Pearson correlations with `analyze_correlation(..., method="pearson", weighting="area")`. Global/native means were computed with `compute_basic_stat(..., op="area_weighted_mean", dim=["lat","lon"])`. The executed code does not show an explicit `load_areas`/`areacello` argument; the weighted values reported below are those returned by those `compute_basic_stat` calls.

Key results:

| Field | Native/global mean SST (°C) | Common-cell mean SST (°C) | ΔSST model − WOA23 (°C) | RMSE (°C) | Spatial Pearson r |
|---|---:|---:|---:|---:|---:|
| CESM2 | 18.52 | 18.48 | +0.28 | 1.062 | 0.995 |
| GFDL-ESM4 | 18.13 | 18.12 | −0.09 | 0.832 | 0.997 |
| MPI-ESM1-2-LR | 17.58 | 17.57 | −0.63 | 1.339 | 0.993 |
| WOA23 | 18.19 | 18.20 | — | — | — |

The WOA23 global mean over its valid 1° SST field, without applying the multi-model comparison mask, is 18.19°C. On the common comparison cells the WOA23 mean is 18.20°C.

Ranking from smallest to largest absolute mean bias:
1. GFDL-ESM4: |ΔSST| = 0.09°C  
2. CESM2: |ΔSST| = 0.28°C  
3. MPI-ESM1-2-LR: |ΔSST| = 0.63°C  

GFDL-ESM4 has the smallest global-mean bias and the lowest RMSE. CESM2 has a moderate warm mean bias, while MPI-ESM1-2-LR has the largest cold mean bias and largest RMSE. All models nevertheless capture the large-scale SST pattern very strongly, with correlations of 0.993–0.997. The bias fields are not uniform: the unweighted spread of the model-minus-WOA23 fields is about 1.03°C standard deviation for CESM2 (range roughly −4.80 to +8.84°C), 0.82°C for GFDL-ESM4 (−4.49 to +7.38°C), and 1.20°C for MPI-ESM1-2-LR (−7.33 to +7.61°C). Thus regional warm and cold biases are larger than the global-mean offsets and partly compensate in the global mean; the spatial maps are important for interpreting model performance beyond the global scalar metrics.

Figure 1 is a 2×2 panel of global maps on the common WOA23 comparison cells: the first three panels show the CESM2, GFDL-ESM4 and MPI-ESM1-2-LR bias fields (model minus WOA23) with a shared diverging color scale centered at zero, clipped at the 98th percentile of pooled absolute bias, enabling direct comparison of sign and magnitude. The fourth panel shows the WOA23 0 m reference SST on the same cells with a separate sequential color scale, clipped at the 1st–99th percentiles. All panels and colorbars are labelled, with SST/bias units in °C. The figure is designed to reveal the regional spatial structure and sign of the biases against the observed SST pattern.

Figure 2 is a bar chart of the four global-mean SST values over their own native/valid grids: CESM2 18.52°C, GFDL-ESM4 18.13°C, MPI-ESM1-2-LR 17.58°C and WOA23 18.19°C. Each bar is annotated with its value, and the y axis is labelled in °C. The chart summarizes the overall global-mean SST offsets before focusing on the common-cell comparison statistics and maps.
