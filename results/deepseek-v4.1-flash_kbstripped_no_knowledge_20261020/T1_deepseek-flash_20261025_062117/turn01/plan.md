# Plan

1. Use load_regular_data to load ERSSTv5 monthly `sst` for 1950-01 to 2014-12 on its native 2-degree grid, and separately load CMIP6 historical monthly `tos` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR over the same window on each model's native grid. Preserve each source's native cell-area information for weighting.
2. For each of the four sources, use compute_enso_nino_index to form the monthly Nino 3.4 SST anomaly series for 5S-5N, 170W-120W, area-weighted with that source's native cell areas. The anomaly is relative to the 1950-2014 monthly climatology, which is exactly the loaded 780-month window.
3. For each source, use compute_basic_stat to compute the population standard deviation (`ddof=0`) over all 780 monthly Nino 3.4 anomalies. Keep the unrounded values, then compute model-minus-ERSSTv5 amplitude bias and model/ERSSTv5 amplitude ratio, rank the three models by smallest to largest absolute amplitude bias, and round only for reporting.
4. For each source, use compute_basic_stat on each calendar-month subset (65 January values, 65 February values, etc.) to compute the 12 population standard deviations (`ddof=0`). Use locate_extreme on each 12-value seasonal-amplitude series to identify the calendar month and value of the maximum and minimum; compute the maximum/minimum ratio. Compare each model's peak month, peak value and ratio with ERSSTv5 to judge whether the observed peak month and seasonal phase-locking strength are reproduced.
5. For each source, use compute_exceedance_fraction by count twice on the Nino 3.4 anomaly series: once for anomalies > +1.0 degrees Celsius and once for anomalies < -1.0 degrees Celsius. Sum the two counts to get the number of months with |anomaly| > 1.0 degrees Celsius; divide by 780 and multiply by 100 for the percentage, rounded to 1 decimal place at reporting.
6. Produce Figure 1 as a line plot of calendar month (January to December) versus the 12 seasonal standard deviations for ERSSTv5, CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, with the observation visually distinct and each curve's peak month marked. Produce Figure 2 as four stacked time-series panels in the order ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, sharing a 1950-2014 x axis and an identical y-axis range covering all four anomaly series, with +/-1.0 degrees Celsius reference lines and each panel annotated with its all-month standard deviation and exceedance percentage.

Other notes:
- All Nino 3.4 anomalies use a 1950-2014 calendar-month climatology and exactly 780 months.
- Population standard deviations use `ddof=0`; each calendar-month standard deviation is over 65 years.
- Exceedance is strict |anomaly| > 1.0 degrees Celsius, not >=, and is counted out of 780 months.
- Nino 3.4 is 5S-5N, 170W-120W; each source is area-weighted on its native grid without regridding.
- ERSSTv5 is the 2-degree in-situ observational product; CMIP6 historical `tos` are the three named model runs.
- Report bias, ratio, ranking and all standard deviations after rounding to the requested precision, using unrounded intermediates.
