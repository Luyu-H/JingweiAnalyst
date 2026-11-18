# Plan

1. Open the raw files with xarray/numpy/pandas: `ersst` monthly `sst` for 1950-01 to 2014-12 on its native 2-degree grid; and `cmip6_historical` monthly `tos` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR over the same window, each retained on its own native/archive grid. Also load each CMIP6 model’s published `areacello` fixed field for weighting. Do not substitute `thetao` for `tos`.
- Note [cmip6-tos-bulk-top-cell]: use `tos` directly with the default `depth_level=None`; it is the bulk temperature of the top model cell, not a skin value.
- Note [cmip6-static-grid-fields]: use each model’s published `areacello` as the area weight, with zero/missing values over land, and report the origin recorded in `attrs['area_source']` when it feeds the regional average.
- Note [cmip6-grid-labels]: GFDL-ESM4 is on the archive `gr` 1x1 degree regridded grid with a regridded land-sea mask; CESM2 and MPI-ESM1-2-LR native horizontal resolutions differ by about 2.5x, so narrow features are not comparable in amplitude across grids.

2. For each of the four sources, compute the monthly Nino 3.4 regional mean SST over 5S-5N, 170W-120W on that source’s own grid. For CMIP6 use the model’s `areacello`; for ERSSTv5 on its regular 2-degree grid use cosine-latitude cell-area weights. Restrict weights to the same valid cells as the data and keep the operation NaN-aware. Do not regrid any source.
- Note [omip-weighted-means-and-integrals]: the Nino 3.4 average must be area-weighted with weights restricted to the same valid cells as the data; an unweighted mean over latitude/longitude is wrong on a curvilinear grid.
- Note [omip-native-grid-vs-regridded]: compute every index and statistic on each source’s own grid with its own cell areas; no regridding is needed, so no regrid method is quoted.

3. Build each source’s monthly Nino 3.4 anomaly series over 1950-01 to 2014-12 by subtracting that source’s own 1950-2014 monthly climatology, matching each month to the climatology of its own calendar month. Record and report the climatology and baseline used for each source.
- Note [anomaly-baselines-and-climate-indices]: do not use `climatology='woa23'`; each anomaly must come from the field’s own 1950-2014 monthly climatology, with calendar-month matching, and the exceedance metric counts months with |anomaly| > 1.0 degree Celsius, so persistence (`min_duration`) plays no role.

4. For each source, compute the population standard deviation (`ddof=0`) of the Nino 3.4 anomaly over all 780 months. For each model, compute amplitude bias as model standard deviation minus the ERSSTv5 standard deviation, amplitude ratio as model standard deviation divided by the ERSSTv5 standard deviation, and rank the three models from smallest to largest absolute amplitude bias. Also compute, for each source, the population standard deviation of the anomaly for each of the 12 calendar months, identify the calendar month of largest and smallest standard deviation, and compute the largest-to-smallest ratio. Compare only these statistics between models and observations.
- Note [cmip6-experiments-historical-ssp245]: coupled historical runs have their own internal variability, so compare only statistics and never individual years or event timing; with one member per model, differences in amplitude cannot be split into forced versus internal-variability parts.

5. For each source, count the months over 1950-01 to 2014-12 whose Nino 3.4 anomaly strictly exceeds 1.0 degree Celsius in absolute value, and compute the percentage of the 780 months to 1 decimal place, using the strict |anomaly| > 1.0 criterion from [anomaly-baselines-and-climate-indices].

6. Make Figure 1 using the 12 monthly standard deviations for ERSSTv5, CESM2, GFDL-ESM4 and MPI-ESM1-2-LR: calendar month January-December on the x axis, standard deviation in degrees Celsius on the y axis, four curves on one axis, with observations visually distinct and each curve’s peak month marked. Make Figure 2 as four stacked panels in the order ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, sharing one time axis over 1950-2014 and an identical y-axis range containing every anomaly series, with +1.0 and -1.0 degree Celsius levels marked and each panel annotated with its full-period standard deviation and exceedance percentage.

Other notes:
- Each source’s anomaly is relative to that source’s own 1950-2014 monthly climatology; no external climatology such as `woa23` is used.
- The model-minus-ERSSTv5 amplitude difference compares each model’s bulk top-cell `tos` with an observational SST product; `tos` is not a skin temperature.
- With a single historical member per model, model-versus-ERSSTv5 differences cannot be separated into forced response and internal variability.
- GFDL-ESM4 numbers are on the archive 1x1 degree `gr` grid; CESM2 and MPI-ESM1-2-LR are on grids whose horizontal resolutions differ by about 2.5x, so features narrower than the coarse grid are not comparable in amplitude.
- [model-gfdl-esm4-reported-biases] GFDL-ESM4’s documented slightly stronger ENSO with a longer period makes a positive amplitude bias for GFDL-ESM4 consistent with its documentation; this is context, not a target value.
