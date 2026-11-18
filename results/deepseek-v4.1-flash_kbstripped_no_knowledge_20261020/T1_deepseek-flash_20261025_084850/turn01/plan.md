# Plan

1. Load data with load_regular_data: from cmip6_historical, variable `tos`, models CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, monthly 1950-01 to 2014-12, each on its native grid; from ersst, variable `sst`, monthly over the same window on its native 2-degree grid.
2. For each of the four sources separately, use compute_enso_nino_index to form the monthly Nino 3.4 index: area-weighted SST anomaly over 5S-5N, 170W-120W, using each source's native cell areas and anomaly relative to that source's 1950-2014 monthly climatology. This produces 780 monthly anomalies per source.
3. Use compute_basic_stat with population standard deviation (ddof = 0) over the full 780-month time dimension for each Nino 3.4 anomaly series. Report the four full-period standard deviations in degrees Celsius to 2 decimal places.
4. From the full-period standard deviations, compute each model's amplitude bias as model std minus ERSSTv5 std and amplitude ratio as model std divided by ERSSTv5 std; report bias and ratio to 2 decimal places, and rank models from smallest to largest absolute amplitude bias.
5. For each source, group the 780 monthly anomalies by calendar month (January-December; 65 values per calendar month) and use compute_basic_stat with population standard deviation (ddof = 0) on each calendar-month group. Report all 48 calendar-month standard deviations in degrees Celsius to 2 decimal places.
6. For each source, identify the calendar month and value of maximum standard deviation, the calendar month and value of minimum standard deviation, and compute the largest-to-smallest ratio to 2 decimal places. Compare each model's peak month and peak amplitude/phase-locking strength with ERSSTv5's observed peak month and strength.
7. For each source, use compute_exceedance_fraction with the threshold absolute Nino 3.4 anomaly strictly greater than 1.0 degrees Celsius, counted over months, to obtain the number and percentage of the 780 months; report the count and percentage to 1 decimal place.
8. Use custom plotting code to make Figure 1: calendar month January-December on the x axis, population standard deviation in degrees Celsius on the y axis, four curves on one axis, ERSSTv5 visually distinct from the three models, and each curve's peak month marked.
9. Use custom plotting code to make Figure 2: four stacked panels sharing one time axis over 1950-2014, in the order ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, each showing that source's Nino 3.4 anomaly series on an identical y-axis range that contains every series, with +1.0 and -1.0 degrees Celsius levels marked and each panel annotated with its full-period standard deviation and exceedance percentage.
Other notes:
- Nino 3.4 is 5S-5N, 170W-120W, and anomalies are relative to each source's 1950-2014 monthly climatology.
- Use population standard deviation with ddof = 0 throughout.
- No horizontal regridding is performed; each index uses the source's native grid and native cell areas.
- The exceedance criterion is strict absolute anomaly greater than 1.0 degrees Celsius.
- 1950-2014 contains 780 months, with 65 values for each calendar month.
- Keep intermediate values unrounded; round only reported values to the requested precision.
- Figure 2 panels must share one y-axis range and use the stated panel order.
