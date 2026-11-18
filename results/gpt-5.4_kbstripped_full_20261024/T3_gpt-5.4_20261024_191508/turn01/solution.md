I evaluated 1971–2000 SST climatologies from the CMIP6 historical `tos` fields for CESM2, GFDL-ESM4, and MPI-ESM1-2-LR against the WOA23 annual 1971–2000 temperature climatology at 0 m depth.

Method used, exactly as run:
- For each CMIP6 model, monthly `tos` was loaded for 1971-01 to 2000-12 from `cmip6_historical`.
- Each model climatology was formed with `temporal_aggregate(data=..., freq="mean", op="mean")`.
- Native-grid global means were computed with `compute_basic_stat(op="area_weighted_mean", dim=["lat","lon"])`.
- WOA23 `temperature` was loaded for `period="1971-2000"`, `time_slice="annual"`, `depth_range=[0.0,0.0]`.
- The WOA23 land-sea mask was loaded and ocean was defined as `landsea_01 > 1`; that mask was applied before computing the WOA23 global mean.
- Each model climatology was then regridded to the WOA23 1° grid using `regrid_horizontal(..., target="woa23", method="conservative", skipna=True, na_thres=0.5)`.
- A single shared comparison mask was built across all three regridded models and WOA23 with `build_comparison_mask(..., ocean_only=True, reduce_other_dims="all")`.
- On that shared mask, area-weighted means, bias, and RMSE were computed using WOA23 grid-cell areas; spatial Pearson correlation was computed with area weighting.
- The common comparison sample contains 40,830 grid cells.

Direct answer
The three models all reproduce the large-scale global SST pattern of WOA23 very well, with spatial correlations above 0.99 on the shared WOA23 1° comparison grid. Among them, GFDL-ESM4 is the closest overall to WOA23 by mean bias, RMSE, and correlation. CESM2 is slightly too warm globally, while MPI-ESM1-2-LR is notably too cool.

Global mean SST on each source’s own grid
- CESM2: 18.52 °C
- GFDL-ESM4: 18.13 °C
- MPI-ESM1-2-LR: 17.58 °C
- WOA23 (0 m annual field, 1° ocean grid): 18.19 °C

Shared-grid comparison results
These are all computed on the same 40,830 WOA23-grid cells after conservative regridding and common masking.

- WOA23 mean on shared cells: 18.20 °C for all three comparisons

CESM2
- Model mean SST: 18.48 °C
- WOA23 mean SST: 18.20 °C
- Mean bias, model minus WOA23: +0.28 °C
- RMSE: 1.062 °C
- Spatial Pearson correlation: 0.995

GFDL-ESM4
- Model mean SST: 18.12 °C
- WOA23 mean SST: 18.20 °C
- Mean bias, model minus WOA23: -0.09 °C
- RMSE: 0.832 °C
- Spatial Pearson correlation: 0.997

MPI-ESM1-2-LR
- Model mean SST: 17.57 °C
- WOA23 mean SST: 18.20 °C
- Mean bias, model minus WOA23: -0.63 °C
- RMSE: 1.339 °C
- Spatial Pearson correlation: 0.993

Ranking by smallest to largest absolute mean bias
1. GFDL-ESM4: 0.09 °C
2. CESM2: 0.28 °C
3. MPI-ESM1-2-LR: 0.63 °C

Interpretation
- GFDL-ESM4 has the best overall agreement with WOA23: near-zero mean bias, lowest RMSE, and highest spatial correlation.
- CESM2 is also very close structurally, but with a modest warm bias.
- MPI-ESM1-2-LR captures the pattern well but is cooler overall and has the largest RMS error.

Spatial structure of the bias
The bias fields show that differences are not spatially uniform. All three models retain the observed broad SST pattern very closely, but regional departures remain:
- CESM2 shows a generally warm-leaning bias pattern overall, consistent with its positive global mean bias.
- GFDL-ESM4 shows the smallest-amplitude bias field overall, consistent with its lower RMSE and near-zero mean bias.
- MPI-ESM1-2-LR shows broader cool-biased structure, consistent with its negative mean bias and larger RMSE.
- The largest local errors are regional rather than basin-wide mean offsets, as shown by maximum absolute errors of about 7–9 °C despite small global mean biases.
- Because correlations are extremely high for all models, the main differences are amplitude and regional mismatch rather than failure to reproduce the global SST pattern.

Visualization description
Figure 1
- A 2×2 panel was generated.
- Three panels map `model - WOA23` SST bias for CESM2, GFDL-ESM4, and MPI-ESM1-2-LR on the shared masked WOA23 1° grid.
- Those three bias maps use one shared diverging color scale centered on zero, with color-bar label in °C.
- The fourth panel shows the masked WOA23 0 m annual SST field on the same grid, with a sequential color scale and color-bar label in °C.
- The plotting used percentile-based color scaling: the bias range was based on the 98th percentile of absolute bias values, and the WOA23 SST panel used the 2nd to 98th percentiles.

Figure 2
- A bar chart was generated showing the four own-grid global mean SST values: CESM2, GFDL-ESM4, MPI-ESM1-2-LR, and WOA23.
- Each bar is annotated with its value, and the y-axis is labeled in °C.

Main takeaways
- All three CMIP6 historical simulations reproduce observed large-scale SST geography very well for 1971–2000.
- GFDL-ESM4 is the best of the three against WOA23 by all reported shared-grid skill metrics.
- CESM2 is slightly too warm on average.
- MPI-ESM1-2-LR is systematically cooler and has the largest spatial error magnitude.
- The dominant model-observation differences are regional bias structures superimposed on otherwise very realistic global SST climatologies.
