# Plan

1. Load the monthly SST records over 1950-01 to 2014-12 with `load_regular_data`: ERSSTv5 variable `sst` on its native 2-degree grid, and CMIP6 historical variable `tos` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, each on its native model grid. Load each source's native grid-cell areas with `load_areas` for area weighting.
2. For each of the four sources (ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR), compute the monthly Nino 3.4 index with `compute_enso_nino_index`: area-weighted SST anomaly over 5S-5N, 170W-120W, using that source's native cell areas, with anomalies defined relative to the 1950-2014 monthly climatology of the same source. This yields four 780-month anomaly series.
3. For each source, compute the population standard deviation (`ddof=0`) of its Nino 3.4 anomaly over all 780 months with `compute_basic_stat` reducing over time. Report each value in degrees Celsius to 2 decimal places.
4. For each model, compute amplitude bias as model standard deviation minus ERSSTv5 standard deviation and amplitude ratio as model standard deviation divided by ERSSTv5 standard deviation, using unrounded intermediate values; report bias and ratio to 2 decimal places. Rank the three models from smallest to largest absolute amplitude bias.
5. For each source, split the monthly Nino 3.4 anomaly series by calendar month (January-December) and compute the population standard deviation (`ddof=0`) for each calendar month with `compute_basic_stat` reducing over the years in that month. Report all 48 seasonal-amplitude values to 2 decimal places.
6. For each source, identify the calendar month of largest seasonal standard deviation and its value, the calendar month of smallest seasonal standard deviation, and the largest-to-smallest ratio to 2 decimal places. Compare each model's peak month and largest-to-smallest ratio with ERSSTv5 to state whether it reproduces the observed peak month and strength of seasonal phase locking.
7. For each source, count the months whose Nino 3.4 anomaly strictly exceeds 1.0 degrees Celsius in absolute value and compute the percentage of the 780 months to 1 decimal place; this is an own-code tabulation from the monthly index series.
8. Produce Figure 1 as an own-code line plot of the seasonal amplitude cycle: calendar month January-December on the x axis, standard deviation in degrees Celsius on the y axis, four curves on one axis, ERSSTv5 visually distinct from the three models, and each curve's peak month marked.
9. Produce Figure 2 as an own-code four-panel stacked plot sharing one 1950-2014 time axis, ordered ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, each showing that source's Nino 3.4 anomaly series on an identical y-axis range containing every series, with +1.0 and -1.0 degrees Celsius levels marked and each panel annotated with its standard deviation and exceedance percentage.
Other notes:
- All anomalies are relative to each source's own 1950-2014 monthly climatology, not a common observational climatology.
- All reported standard deviations are population standard deviations with `ddof=0`.
- The Nino 3.4 region is exactly 5S-5N, 170W-120W, and area weighting uses each source's native grid-cell areas without horizontal regridding.
- ERSSTv5 is the in situ-only 2-degree observational product; CMIP6 results use the provided historical run for each model on its native grid.
- Exceedance uses strict absolute value greater than 1.0 degrees Celsius, with denominator 780 months.
- Seasonal phase-locking strength is summarized by the largest-to-smallest ratio of the 12 calendar-month standard deviations.
- Values are rounded only when reported, to the precision requested.
