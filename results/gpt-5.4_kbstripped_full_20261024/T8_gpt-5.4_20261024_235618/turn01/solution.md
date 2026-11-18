I evaluated the 2002 annual-mean ocean carbon system in three CMIP6 historical models—CESM2, GFDL-ESM4, and MPI-ESM1-2-LR—against the GLODAP mapped climatology, using:
- CMIP6 monthly `dissic`, `talk`, `ph`, `thetao`, and `so` for January–December 2002,
- CMIP6 monthly `dissic`, `thetao`, and `so` for January 1971–December 2000,
- GLODAP climatology `TCO2`, `TAlk`, and `pHtsinsitutp`.

Method actually used in the executed code
- The code formed simple time means from the monthly CMIP6 fields using `temporal_aggregate(..., freq="mean", op="mean")`.
- For DIC and alkalinity, model values were converted from per-volume to per-mass units with `transform_unit(..., to_unit="umol kg-1", density=...)`, where density came from `compute_potential_density_anomaly(..., ref_pressure=None)` applied to the model annual-mean or 1971–2000 mean `thetao` and `so`.
- GLODAP 0–50 m scalar means were computed directly from the loaded climatology with `compute_basic_stat(..., op="volume_weighted_mean")`.
- Model 0–50 m scalar means were also computed with `compute_basic_stat(..., op="volume_weighted_mean")` after selecting `depth=slice(0,50)`.
- For the 0–50 m DIC maps, each source was first collapsed vertically with `compute_basic_stat(..., op="mean", dim="depth")`.
- Model DIC maps were then conservatively regridded to the GLODAP grid using `regrid_horizontal(..., method="conservative", skipna=True, na_thres=0.5)`.
- A single shared comparison mask was built across GLODAP and all three regridded model maps with `build_comparison_mask(..., ocean_only=False, reduce_other_dims="all")`.
- Map bias, RMSE, and cell count were computed with `compute_difference_metrics(..., weighting="auto")`; spatial correlation was computed with `analyze_correlation(..., method="pearson", weighting="area")`.
- For the depth comparison, model and GLODAP DIC were vertically interpolated to exact target depths 0, 100, 500, 1000, and 2000 m with `regrid_vertical(..., method="linear", fill_surface=True)`, then model slices were conservatively regridded horizontally to the GLODAP grid before global mean bias was computed.

Main scalar 0–50 m results
GLODAP reference:
- DIC: 2036 µmol kg⁻¹
- Total alkalinity: 2304 µmol kg⁻¹
- pH (total scale): 8.080

Models, with difference defined as model minus GLODAP:

CESM2
- DIC: 1997 µmol kg⁻¹, bias = -39 µmol kg⁻¹
- Total alkalinity: 2264 µmol kg⁻¹, bias = -40 µmol kg⁻¹
- pH: 8.080, bias = -0.001

GFDL-ESM4
- DIC: 2059 µmol kg⁻¹, bias = +23 µmol kg⁻¹
- Total alkalinity: 2336 µmol kg⁻¹, bias = +31 µmol kg⁻¹
- pH: 8.088, bias = +0.008

MPI-ESM1-2-LR
- DIC: 1974 µmol kg⁻¹, bias = -61 µmol kg⁻¹
- Total alkalinity: 2218 µmol kg⁻¹, bias = -87 µmol kg⁻¹
- pH: 8.059, bias = -0.022

Interpretation of the scalar comparison
- CESM2 is low in both DIC and alkalinity relative to GLODAP, but its surface-layer mean pH is nearly unbiased.
- GFDL-ESM4 is high in DIC, alkalinity, and pH.
- MPI-ESM1-2-LR is low in all three quantities, with the largest negative alkalinity and pH errors of the three models.

Thus, the sign of the errors is internally consistent within every model:
- CESM2: DIC low, alkalinity low, pH low
- GFDL-ESM4: DIC high, alkalinity high, pH high
- MPI-ESM1-2-LR: DIC low, alkalinity low, pH low

0–50 m DIC map skill on the GLODAP 1° grid
All map metrics were computed on the common shared mask, with 41,227 valid comparison cells.

CESM2
- Mean bias: -32 µmol kg⁻¹
- RMSE: 57 µmol kg⁻¹
- Spatial correlation: 0.830

GFDL-ESM4
- Mean bias: +28 µmol kg⁻¹
- RMSE: 66 µmol kg⁻¹
- Spatial correlation: 0.714

MPI-ESM1-2-LR
- Mean bias: -55 µmol kg⁻¹
- RMSE: 80 µmol kg⁻¹
- Spatial correlation: 0.740

Interpretation of the map metrics
- CESM2 gives the best spatial reproduction of the 0–50 m DIC field: it has the lowest RMSE and the highest spatial correlation.
- GFDL-ESM4 has the smallest absolute global 0–50 m DIC scalar bias (+23 µmol kg⁻¹), but weaker spatial skill than CESM2.
- MPI-ESM1-2-LR performs worst on the 0–50 m DIC maps, with the largest negative mean bias and largest RMSE.

Depthwise global mean DIC bias relative to GLODAP
Biases are model minus GLODAP after vertical interpolation to exact depths and horizontal regridding to the GLODAP grid.

CESM2
- 0 m: -29 µmol kg⁻¹
- 100 m: -46 µmol kg⁻¹
- 500 m: -12 µmol kg⁻¹
- 1000 m: +2 µmol kg⁻¹
- 2000 m: -2 µmol kg⁻¹

GFDL-ESM4
- 0 m: +29 µmol kg⁻¹
- 100 m: +19 µmol kg⁻¹
- 500 m: -11 µmol kg⁻¹
- 1000 m: -1 µmol kg⁻¹
- 2000 m: -11 µmol kg⁻¹

MPI-ESM1-2-LR
- 0 m: -58 µmol kg⁻¹
- 100 m: -59 µmol kg⁻¹
- 500 m: -78 µmol kg⁻¹
- 1000 m: -67 µmol kg⁻¹
- 2000 m: -61 µmol kg⁻¹

Interpretation of the vertical structure
- CESM2 is too low near the surface and at 100 m, but its bias weakens strongly with depth and is near zero by 1000–2000 m.
- GFDL-ESM4 is too high near the surface, then becomes slightly low below 500 m; its deep biases remain modest.
- MPI-ESM1-2-LR shows a large negative DIC bias at all sampled depths, indicating a pervasive low-carbon mean state rather than only a surface-layer problem.

2002 minus 1971–2000 change in the global 0–50 m mean DIC
CESM2
- Change: +9 µmol kg⁻¹
- Ratio to 2002 DIC bias vs GLODAP: -0.24

GFDL-ESM4
- Change: +12 µmol kg⁻¹
- Ratio to 2002 DIC bias vs GLODAP: 0.53

MPI-ESM1-2-LR
- Change: +11 µmol kg⁻¹
- Ratio to 2002 DIC bias vs GLODAP: -0.18

Interpretation of the period-change ratios
- CESM2 and MPI-ESM1-2-LR: the 2002 minus 1971–2000 DIC increase is small compared with their 2002 mismatch to GLODAP. That implies most of the apparent 2002 bias should be read as mean-state model error rather than simply a consequence of temporal evolution between the historical mean and 2002.
- GFDL-ESM4: the period change is about half the size of the 2002 GLODAP mismatch (ratio 0.53), so part of the apparent 2002 difference could reflect temporal evolution, but a substantial mean-state bias still remains.

Which model reproduces the observed carbon system best?
Overall, GFDL-ESM4 reproduces the observed carbon system best in the combined sense used by the executed code. It had:
- the smallest absolute scalar DIC bias,
- relatively small alkalinity and pH biases,
- modest depthwise DIC biases compared with MPI-ESM1-2-LR,
- and the lowest overall composite score assembled from scalar biases, map bias/RMSE/correlation penalty, and mean absolute depth bias.

However, there is an important nuance:
- If the main priority is the spatial pattern of upper-ocean DIC, CESM2 performs best, with the highest map correlation (0.830) and lowest RMSE (57 µmol kg⁻¹).
- If the priority is overall agreement across DIC, alkalinity, pH, and vertical-structure metrics taken together, the executed code ranks GFDL-ESM4 as best overall.

What the generated figures show
Figure 1
- Three global maps in one row, one for each model.
- Each map shows the 0–50 m DIC bias relative to GLODAP after the model map was regridded to the GLODAP grid and masked to the shared valid comparison cells.
- All panels use the same diverging color scale centered on zero, so positive and negative biases are directly comparable across models.

Figure 2
- A single-panel vertical profile plot of global mean DIC bias versus depth.
- One line is drawn per model.
- Depth increases downward, bias is on the horizontal axis, a vertical zero-bias line is included, and markers are placed at 0, 100, 500, 1000, and 2000 m.
- Its purpose is to compare how each model’s DIC bias changes with depth after interpolation to common depths and regridding to the GLODAP grid.

Main scientific takeaways
- The three models separate cleanly by sign:
  - CESM2 and MPI-ESM1-2-LR are generally too low in DIC, alkalinity, and pH.
  - GFDL-ESM4 is generally too high in all three.
- MPI-ESM1-2-LR has the most systematic low-DIC bias, extending from the surface to 2000 m.
- CESM2 has strong spatial skill in upper-ocean DIC despite being too low in the global mean.
- GFDL-ESM4 is the best all-round match to GLODAP in the combined metrics used here, though it tends to overestimate the upper-ocean carbon system and has weaker spatial DIC pattern skill than CESM2.

Caveat from the executed workflow
- The DIC and alkalinity comparison was based on converting CMIP6 `dissic` and `talk` to µmol kg⁻¹ using density derived from model `thetao` and `so` with `ref_pressure=None`, exactly as done in the executed code.
