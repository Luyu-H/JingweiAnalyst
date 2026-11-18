I assessed sea surface salinity (SSS) in CESM2, GFDL-ESM4, and MPI-ESM1-2-LR using the executed workflow exactly as run.

Method used
- WOA23 “salinity” annual climatology for period 1971–2000 was loaded and converted to psu.
- WOA23 SSS was formed from the selected WOA23 depths returned by the depth-level loader for depth_range [0, 10], which were 0, 5, and 10 m, using thickness weights from the returned depth bounds: 2.5, 5.0, and 5.0 m. The code therefore computed a thickness-weighted mean across those three levels, then masked to ocean where WOA23 landsea > 1.
- For each CMIP6 model, monthly historical so for 1950-01 to 2014-12 and monthly ssp245 so for 2071-01 to 2100-12 were loaded and converted to psu.
- Model SSS was formed from the model depth levels returned for depth_range [0, 10], using the loaded CMIP6 volume field over those selected depths. The code computed the sum of salinity × volume over depth divided by the sum of volume over depth.
- Global means on native grids were computed with area_weighted_mean.
- For mean-state comparison, each model’s 1971–2000 climatology was regridded to the WOA23 grid with conservative regridding, a shared comparison mask was built from model and WOA23 fields, and RMSE and area-weighted spatial Pearson correlation were computed on that mask.
- For time series and trends, monthly global means were computed first, then aggregated to annual means. Linear trends were fit separately for historical 1950–2014 and ssp245 2071–2100, with confidence intervals and significance coming from the trend tool’s serial-correlation-adjusted effective sample size.
- End-of-century change was computed as the 2071–2100 mean minus the 1995–2014 mean. Global changes were computed on native grids. Basin changes were computed after conservative regridding to WOA23 and then selecting Atlantic, Indian, and Pacific basins with marginal seas included.

Direct answer

1) Mean-state bias against WOA23
WOA23 global mean 1971–2000 SSS: 34.61 psu.

Model global means and biases relative to WOA23:
- CESM2: 34.32 psu, bias -0.29 psu
- GFDL-ESM4: 34.56 psu, bias -0.05 psu
- MPI-ESM1-2-LR: 34.36 psu, bias -0.25 psu

Map-based skill on the WOA23 1° grid:
- CESM2: RMSE 1.07 psu; spatial r = 0.815; comparison cells = 40,947
- GFDL-ESM4: RMSE 0.96 psu; spatial r = 0.839; comparison cells = 41,023
- MPI-ESM1-2-LR: RMSE 1.06 psu; spatial r = 0.807; comparison cells = 40,899

Interpretation:
- All three models are slightly fresher than WOA23 in their global 1971–2000 mean.
- GFDL-ESM4 has the smallest global-mean bias and the best spatial skill of the three by both RMSE and correlation.
- CESM2 and MPI-ESM1-2-LR are similar in overall mean-state error and both are notably fresher than WOA23 globally.

2) Historical and end-of-century trends in global-mean SSS

Historical 1950–2014 global-mean SSS trend:
- CESM2: -0.00268 psu/decade
  - 95% CI: [-0.00693, 0.00157]
  - R² = 0.245
  - CI excludes zero? No
- GFDL-ESM4: -0.00110 psu/decade
  - 95% CI: [-0.00361, 0.00141]
  - R² = 0.064
  - CI excludes zero? No
- MPI-ESM1-2-LR: -0.000551 psu/decade
  - 95% CI: [-0.00504, 0.00393]
  - R² = 0.007
  - CI excludes zero? No

ssp245 2071–2100 global-mean SSS trend:
- CESM2: -0.0209 psu/decade
  - 95% CI: [-0.0255, -0.0164]
  - R² = 0.910
  - CI excludes zero? Yes
- GFDL-ESM4: -0.00534 psu/decade
  - 95% CI: [-0.0106, -0.000039]
  - R² = 0.324
  - CI excludes zero? Yes
- MPI-ESM1-2-LR: -0.0231 psu/decade
  - 95% CI: [-0.0513, 0.00504]
  - R² = 0.701
  - CI excludes zero? No

ssp245-to-historical slope ratio:
- CESM2: 7.8
  - Not meaningful, because the executed analysis flagged CESM2 historical and ssp245 as different local members.
- GFDL-ESM4: 4.9
  - Meaningful under the member-continuity note used in the code.
- MPI-ESM1-2-LR: 42.0
  - Meaningful under the member-continuity note used in the code, but the historical slope is extremely small, so the ratio is numerically unstable and should be interpreted cautiously.

Interpretation:
- None of the three models shows a statistically significant historical global-mean trend over 1950–2014 once serial correlation is accounted for.
- Under late-century ssp245, CESM2 and GFDL-ESM4 show significant freshening of global-mean SSS; MPI-ESM1-2-LR also freshens, but its adjusted 95% interval still overlaps zero.

3) End-of-century mean change under ssp245
Change defined by the code as 2071–2100 mean minus 1995–2014 mean, in psu.

Global ocean:
- CESM2: -0.12
- GFDL-ESM4: -0.01
- MPI-ESM1-2-LR: -0.09

Atlantic:
- CESM2: +0.10
- GFDL-ESM4: +0.19
- MPI-ESM1-2-LR: -0.04

Indian:
- CESM2: +0.06
- GFDL-ESM4: -0.00
- MPI-ESM1-2-LR: -0.05

Pacific:
- CESM2: -0.24
- GFDL-ESM4: -0.04
- MPI-ESM1-2-LR: -0.11

Experiment-member caveat from the executed workflow:
- CESM2: interpret cross-experiment changes with caution, because historical and ssp245 use different local members, so the change mixes forced response with internal-variability differences between members.
- GFDL-ESM4 and MPI-ESM1-2-LR: the code treated historical r1 and ssp245 r1 as continuous for interpretation.

Interpretation:
- All three models show a negative global mean delta SSS by late century, strongest in CESM2 and MPI-ESM1-2-LR.
- CESM2 and GFDL-ESM4 show Atlantic salinification despite global freshening, while the Pacific freshens in all three models.
- MPI-ESM1-2-LR differs from the other two by freshening in all three requested basins.

Key scientific patterns
- Mean state: GFDL-ESM4 agrees best with WOA23 overall, both in global mean and spatial pattern.
- Historical evolution: global-mean SSS changes are weak and not statistically distinct from zero for 1950–2014 in all three models after accounting for autocorrelation.
- Future evolution: late-century ssp245 tends toward fresher global-mean surface salinity, with strongest and clearest declines in CESM2, weaker but significant decline in GFDL-ESM4, and a large but not adjusted-significant decline in MPI-ESM1-2-LR.
- Basin contrast: the Atlantic tends to become saltier in CESM2 and GFDL-ESM4 while the Pacific becomes fresher, suggesting an intensified inter-basin salinity contrast in those models.

What the figures show
- Figure 1: three side-by-side global maps of 1971–2000 SSS bias, defined as model climatology minus WOA23 after conservative regridding to the WOA23 grid and masking to cells common to each model-observation pair. All panels share one diverging color scale centered on zero, with color bar in psu.
- Figure 2: three side-by-side global maps of cellwise historical SSS trend over 1950–2014, calculated from annual means on each model’s own grid, then regridded to a regular WOA23 longitude-latitude grid for display. A common diverging color scale is used, and grid cells significant at p < 0.05 from the corrected trend test are stippled.
- Figure 3: one panel with annual global-mean SSS for all three models, showing historical 1950–2014 and ssp245 2071–2100 as separate segments, with the 2015–2070 gap left blank. Each segment has its own fitted line, and model colors are consistent across figures.

Important caveat on methodology
- Although the original request specified full-cell-volume weighting for WOA23 levels within 0–10 m, the executed code used WOA23 depth-thickness weights derived from the returned depth bounds, not cell volume weights. For the CMIP6 models, the executed code did use volume weighting over the selected top levels. The summary above reports exactly what the code did.
