# Plan

1. Open `cmip6_historical` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR: variable `tos`, monthly, 1950-01 to 2014-12, top model level, on each model's native grid; open each model's fixed `areacello`. Open `ersst`: variable `sst`, monthly, 1950-01 to 2014-12, on its native 2-degree grid. Identify the single historical member available for each model and the 780 matching months.
- Note [cmip6-tos-bulk-top-cell]: use the top `tos` level (`depth_level=None`) for the Nino 3.4 index; do not apply an observed (WOA23) climatology to model `tos`.
- Note [cmip6-grid-labels]: GFDL-ESM4 `tos` is `gr` (regridded to 1x1 degree), CESM2 and MPI-ESM1-2-LR are `gn`; resolution differs by about 2.5x, so features narrower than the coarse grid are not comparable in amplitude and amplitude differences partly reflect resolution.

2. For each source, select the native grid cells in the Nino 3.4 box 5S-5N, 170W-120W, handling longitude conventions and wrap, and for each month compute the area-weighted mean SST over those valid cells. Use each model's own `areacello`; for ERSSTv5 compute regular 2-degree cell areas from longitude spacing and cos(latitude), restricted to valid `sst` cells. Mask weights to the valid SST cells in each month and record the area covered and the source of the area field.
- Note [omip-weighted-means-and-integrals]: the 5S-5N, 170W-120W average must be area-weighted with each source's own cell areas; weights must be NaN-restricted to the same valid cells; the reported standard deviation is the population standard deviation (ddof=0) of the index series over all 780 months, not an unweighted map standard deviation.
- Note [cmip6-static-grid-fields]: on curvilinear native grids the model's own `areacello` is the only correct area weight (zero/missing over land); record and report the area-field origin, and state the area each average was taken over rather than rescaling one model to match another; total ocean area differs between models by about a per cent.

3. For each source, form the monthly Nino 3.4 regional-mean SST series, compute its own 1950-2014 calendar-month climatology (12 values, the plain mean of that calendar month over the 65 years), and subtract the matching calendar-month climatology to produce the 780-month Nino 3.4 anomaly series.
- Note [monthly-means-are-time-averages]: no day weights enter the monthly climatology; each calendar month is the plain mean of that month's values over the 1950-2014 baseline; only interannual/inter-month variance is recoverable from monthly means, and the 12 per-month standard deviations are variability across years, not sub-monthly.
- Note [anomaly-baselines-and-climate-indices]: each model and ERSSTv5 index must be its anomaly from its own 1950-2014 monthly climatology, with each calendar month matched to that month's climatology to remove the seasonal cycle; state that baseline; exceedance counts are months, not events, so count months with |anomaly| > 1.0 °C as specified and do not apply event detection/minimum-duration rules.

4. From the unrounded anomaly series, compute for each source: population standard deviation (ddof=0) over all 780 months; for each model, amplitude bias = model standard deviation minus ERSSTv5 standard deviation and amplitude ratio = model standard deviation divided by ERSSTv5 standard deviation, then rank models from smallest to largest absolute bias; the 12 calendar-month standard deviations (ddof=0 across the 65 years for each month); the maximum and minimum calendar-month standard deviations, their months, and the largest-to-smallest ratio; and the count and percentage (to 1 decimal place) of months with |anomaly| > 1.0 °C. Round only reported values: standard deviations, biases and ratios to 2 decimal places, percentages to 1 decimal place. Compare each model's peak calendar month, peak standard deviation and peak/min standard-deviation ratio with ERSSTv5 to judge phase locking.
- Later step uses [anomaly-baselines-and-climate-indices] inline for the strict exceedance count.

5. Produce Figure 1: line plot of the 12 calendar-month standard deviations with January-December on the x axis and degrees Celsius on the y axis, four curves on one axis, ERSSTv5 visually distinct from the three models, and each curve's peak calendar month marked.

6. Produce Figure 2: four stacked panels sharing one 1950-2014 time axis, in the order ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, each showing that source's monthly Nino 3.4 anomaly on an identical y-axis range containing every series; mark +1.0 and -1.0 °C and annotate each panel with its standard deviation (2 decimal places) and exceedance percentage (1 decimal place).

7. Assemble the results tables and figure captions: 4 overall standard deviations, 3 biases/ratios/rankings, 48 seasonal standard deviations, 4 peak/min/ratio/phase-locking summaries, and 4 exceedance counts/percentages; state that the comparison is statistical rather than event-by-event.
- Note [cmip6-experiments-historical-ssp245]: the historical runs are concentration-driven coupled simulations, so a given model year is not expected to match the observed year; compare statistics, not event timing, and caveat that one member per model cannot separate forced change from internal variability.

Other notes:
- The anomaly baseline for every source is its own 1950-2014 monthly climatology; no observed climatology is applied to model output.
- Strict exceedance means |anomaly| > 1.0 °C, not >=; counts are monthly counts out of 780, not events.
- Use unrounded intermediate values throughout and round only the reported precision.
- The exact Nino 3.4 cell set and area differ among sources because of native grids and coastlines; report the area each average was taken over.
- ERSSTv5 is a reconstructed observational product while CMIP6 historical output is coupled-model simulation, so agreement is assessed through distributions and statistics.
