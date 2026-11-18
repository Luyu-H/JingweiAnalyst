# Plan

1. Open the raw files for ERSSTv5 monthly `sst` and CMIP6 historical monthly `tos` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, plus each CMIP6 model’s static `areacello`. Restrict every source to 1950-01 through 2014-12, giving 780 months. Keep each source on its native grid; do not interpolate or remap.
- Note [cmip6-grid-labels]: treat GFDL-ESM4’s `grid_label='gr'` record as the 180×360 1° archive grid and use its published `areacello` on that grid; state this grid difference.

2. For each source, form the Niño 3.4 monthly regional mean over 5S-5N, 170W-120W by selecting native-grid cell centres inside the box, handling longitude conventions correctly. Compute the regional mean as the area-weighted sum divided by the sum of weights, using only valid ocean cells in both numerator and denominator. For ERSSTv5 use cosine-latitude weights on its native 2° grid; for each CMIP6 model use its published `areacello`.
- Note [omip-weighted-means-and-integrals]: the regional mean must be area-weighted over valid cells only; do not use an unweighted latitude/longitude mean.
- Note [cmip6-static-grid-fields]: use each model’s published `areacello` as the horizontal weight, with zero/missing land cells excluded; state the area field used rather than rescaling.

3. For each source, compute its own 1950-2014 calendar-month climatology of the Niño 3.4 regional mean: for each of the 12 calendar months, take the plain multi-year mean across 1950-2014. Subtract the matching calendar-month climatology from each monthly regional mean to obtain that source’s Niño 3.4 anomaly time series. State the 1950-2014 monthly climatology baseline for every source.
- Note [anomaly-baselines-and-climate-indices]: each source uses its own 1950-2014 calendar-month climatology, not an observed or WOA climatology; a different baseline would add a monthly offset. Count months with absolute anomaly greater than 1.0 directly, without event persistence rules, and compute seasonal amplitudes per calendar month rather than DJF/MAM/JJA/SON.
- Note [monthly-means-are-time-averages]: the monthly `tos`/`sst` values are monthly time averages, so the standard deviations and threshold counts describe monthly variability only; the calendar-month climatology is a plain mean across years with no day weighting.

4. For each of the four sources, compute the population standard deviation, ddof=0, of its Niño 3.4 anomaly over all 780 months. Compute each model’s amplitude bias as model standard deviation minus ERSSTv5 standard deviation, and its amplitude ratio as model standard deviation divided by the ERSSTv5 standard deviation. Rank the three models from smallest to largest absolute amplitude bias, using unrounded values internally.
- Note [cmip6-experiments-historical]: the CMIP6 historical runs are coupled and contain internal variability, so a given model year is not expected to match the observed year; compare amplitude statistics over the full 1950-2014 window rather than year-by-year event timing.

5. For each source, compute the population standard deviation, ddof=0, of the anomaly separately for each of the 12 calendar months using the 65 anomaly values from 1950-2014 for that calendar month. Report all 48 values. For each source, identify the calendar month of largest standard deviation and its value, the calendar month of smallest standard deviation and its value, and the largest-to-smallest ratio. Compare each model’s peak month and peak amplitude/largest-to-smallest ratio with ERSSTv5 to judge whether it reproduces the observed peak month and the strength of seasonal phase locking.

6. Using the [anomaly-baselines-and-climate-indices] counting rule, count for each source the number of months whose Niño 3.4 anomaly has absolute value strictly greater than 1.0 degrees Celsius. Compute the percentage as that count divided by the 780 total months, reported to 1 decimal place.

7. Produce Figure 1 as a line plot of the 12 calendar-month standard deviations for ERSSTv5, CESM2, GFDL-ESM4 and MPI-ESM1-2-LR on one axis, with calendar month January to December on the x axis, standard deviation in degrees Celsius on the y axis, the observations visually distinct from the models, and each curve’s peak month marked. Produce Figure 2 as four stacked panels sharing one 1950-2014 time axis, ordered ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, each showing that source’s Niño 3.4 anomaly series with an identical y-axis range containing all series, horizontal plus and minus 1.0 degrees Celsius levels marked, and each panel annotated with its standard deviation and exceedance percentage.

Other notes:
- Compute all intermediate values unrounded; apply the requested rounding only when reporting: standard deviations, biases and ratios to 2 decimal places, percentages to 1 decimal place.
- All standard deviations requested here are population standard deviations with ddof=0, and the exceedance percentage denominator is the full 780-month window.
