# Plan

1. Load the monthly SST data for 1950-01 to 2014-12 from `cmip6_historical` for models CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, variable `tos`, and from `ersst`, variable `sst`, over the Niño 3.4 box (5S–5N, 170W–120W) using `load_regular_data`.
- Note [cmip6-tos-bulk-top-cell]: use `tos` as the model surface temperature field directly, with its native valid-data mask and no added sea-ice mask.
- Note [cmip6-grid-labels]: keep each model on its own archived native grid when loading and analysing the regional index.
- Note [missing-values-and-land-cells]: ensure the field mask comes from the loaded SST variable so missing/fill values are excluded before any weighted reduction.

2. Load the published horizontal cell areas for the three CMIP6 model grids with `load_areas`, matched to CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, and prepare ERSST area weights from its native latitude–longitude grid for the same Niño 3.4 subset.

3. Compute the monthly Niño 3.4 regional mean SST series for each source on its own grid using `compute_basic_stat` with area weighting over latitude and longitude.
- Note [omip-weighted-means-and-integrals]: form each monthly regional mean as an area-weighted mean over valid ocean cells only, with weights masked wherever the SST field is invalid.
- Note [omip-native-grid-vs-regridded]: do not regrid any source; compute the index on each source’s own grid with its own cell areas.

4. Convert each monthly regional-mean SST series to a Niño 3.4 anomaly series by removing that source’s own 1950–2014 monthly climatology: use `temporal_aggregate` to form the 12-month climatology and `compute_anomaly_vs_own_climatology` to subtract it month-by-month.
- Note [anomaly-baselines-and-climate-indices]: define anomalies relative to each source’s own 1950–2014 monthly climatology and compare month with the same calendar month only.
- Note [monthly-means-are-time-averages]: use the plain mean of each calendar month across years for the climatology; no extra day weighting is applied to these monthly mean inputs.

5. From each 780-month Niño 3.4 anomaly series, compute the population standard deviation (ddof = 0) over time, then compute each model’s amplitude bias (model minus ERSSTv5) and amplitude ratio (model divided by ERSSTv5), and rank the three models by absolute amplitude bias. Use `compute_basic_stat` for the standard deviations; do the bias, ratio and ranking from those scalar results.

6. Compute the population standard deviation (ddof = 0) separately for each calendar month across years for each source’s Niño 3.4 anomaly series, using `temporal_aggregate` to group by climatological month and `compute_basic_stat` for the within-month standard deviations. From the 12 monthly values per source, identify the month of maximum and minimum standard deviation, compute the max/min ratio, and assess whether each model reproduces the observed peak month and the observed strength of seasonal phase locking.

7. Count, for each source, the number of months in 1950–2014 whose Niño 3.4 anomaly has absolute value strictly greater than 1.0 °C, and convert to a percentage of 780 months. Use `compute_exceedance_fraction` on the 1D anomaly series for the fraction and derive counts from the same boolean exceedance test.

8. Prepare Figure 1 from the 12 calendar-month standard deviations for ERSSTv5, CESM2, GFDL-ESM4 and MPI-ESM1-2-LR: one line plot with observations visually distinct, calendar month on the x-axis, standard deviation (°C) on the y-axis, and each source’s peak month marked.

9. Prepare Figure 2 from the four Niño 3.4 anomaly time series: four stacked panels in the order ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, sharing the 1950–2014 time axis and a common y-axis range spanning all series, with ±1.0 °C reference lines and each panel annotated with that source’s full-period standard deviation and exceedance percentage.

Other notes:
- [monthly-means-are-time-averages] The reported variability is the variability of monthly mean Niño 3.4 anomalies: over all 780 months for the full-period metric, and within each calendar month across years for the seasonal phase-locking metric.
