# Plan

1. Load the ERSSTv5 monthly SST field from dataset `ersst`, variable `sst`, over 1950-01 to 2023-12 and the Niño 3.4 region box 5S–5N, 170W–120W, on the native grid, using `load_regular_data`.
- Note [anomaly-baselines-and-climate-indices]: Reconstruct Niño 3.4 as monthly SST anomalies over 5S–5N, 170W–120W relative to the stated monthly climatology for each baseline exactly, and for ONI-style event detection use the 3-month running mean with at least 5 consecutive overlapping seasons at or beyond ±0.5 °C on anomalies.

2. Compute the area-weighted monthly regional mean SST series for Niño 3.4 from the loaded ERSST field using `compute_basic_stat` with a spatial weighted mean over latitude and longitude, leaving time intact.

3. Build the two reconstructions from the same ERSST regional mean series:
   - `RECON_9120`: remove the series’ monthly climatology over 1991-01 to 2020-12 using `compute_anomaly_vs_own_climatology`.
   - `RECON_8110`: remove the series’ monthly climatology over 1981-01 to 2010-12 using `compute_anomaly_vs_own_climatology`.

4. Load the published CPC monthly Niño indices from dataset `cpc_indices` using `load_regular_data`:
   - product `ersst5_91-20`, Niño 3.4 anomaly column, 1950-01 to 2023-12 as the main reference;
   - product `ersst5_81-10`, Niño 3.4 anomaly column, 1950-01 to 2023-12 as the control;
   - product `oni`, Niño 3.4 ONI series, 1950-01 to 2023-12.

5. Compare `RECON_9120` against the CPC `ersst5_91-20` Niño 3.4 monthly anomaly series over 1950-01 to 2023-12:
   - use `analyze_correlation` to report the corrected Pearson correlation on the monthly series without detrending;
   - use `compute_difference_metrics` to report RMSE and the largest absolute deviation with its year and month.
- Note [serial-correlation-and-significance]: Treat the reconstruction and CPC series as two estimates of the same quantity, correlate the monthly series directly without detrending, and report the corrected Pearson correlation; use RMSE and largest absolute deviation as direct agreement metrics.

6. Compute the 3-month centred running mean of `RECON_9120` using `temporal_aggregate` with a 3-month running window, then apply `detect_threshold_events` separately for warm events (threshold +0.5 °C) and cool events (threshold -0.5 °C), each with the operational minimum duration of 5 consecutive overlapping 3-month seasons. From the warm-event output, report:
   - the number of El Niño events over 1950-2023;
   - the three events with the highest peak values, each with peak year-month and peak value.
   Use the cool-event output for the Figure 1 La Niña shading.

7. Apply the identical ONI-style procedure to the published CPC `oni` monthly series:
   - use `detect_threshold_events` for El Niño events at or above +0.5 °C with minimum duration 5 seasons;
   - report the count and the three highest peaks with peak year-month and value;
   - compare the event list from Step 6 with CPC ONI and explain any count difference in terms of differences between a running mean applied to reconstructed monthly Niño 3.4 anomalies and the published official ONI series definition/production.

8. Quantify the baseline-period offset in the reconstruction by computing `RECON_8110 / RECON_9120` difference with code outside the listed tools, then reduce it over 1950-01 to 2020-12 with `compute_basic_stat` using the unweighted mean over time to report the mean reconstructed baseline offset.

9. Quantify the published baseline-period offset by differencing the CPC `ersst5_81-10` and `ersst5_91-20` Niño 3.4 anomaly series, then use `compute_basic_stat` with the unweighted mean over time for 1950-01 to 2020-12 to report the mean published baseline offset and compare it directly with the reconstructed offset.

10. Prepare the requested visualizations from the computed series:
   - Figure 1: two-row shared-x time series over 1950-2023, upper row overlaying `RECON_9120`, CPC `ersst5_91-20`, the 3-month running mean, ±0.5 °C thresholds, and warm/cool event shading from Steps 6; lower row plotting `RECON_9120 - CPC`.
   - Figure 2: square scatter plot of `RECON_9120` (y) against CPC `ersst5_91-20` (x), one point per month, identical axis ranges/ticks, 1:1 line, and annotations for correlation and RMSE.

Other notes:
- State the climatology baseline explicitly for each reconstructed or published anomaly series, because the 1981–2010 and 1991–2020 versions differ by baseline offset rather than by a different variability definition.
- Report all requested metrics on the exact windows named in the question: 1950-01 to 2023-12 for agreement and event analysis, and 1950-01 to 2020-12 for mean baseline offsets.
