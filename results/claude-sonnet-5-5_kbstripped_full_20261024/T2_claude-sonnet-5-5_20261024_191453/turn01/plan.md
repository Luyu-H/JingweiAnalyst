# Plan

1. Load the ERSSTv5 monthly SST (`ersst`, variable `sst`) with `load_regular_data`. Use the native 2-degree grid, 1950-01 to 2023-12, over the Nino 3.4 box (5S-5N, 170W-120W).

2. Rebuild the index with `compute_enso_nino_index` on ERSSTv5's own climatology.
   - Build RECON_9120 with the Nino 3.4 region and `base_period` 1991-2020.
   - Build RECON_8110 with the same code and `base_period` 1981-2010, over 1950-01 to 2023-12.
   - Note [anomaly-baselines-and-climate-indices]: match each calendar month to its own climatology. State the climatology and baseline used for each series. Compare RECON_9120 with `ersst5_91-20` and RECON_8110 with `ersst5_81-10`, since the CPC products are not interchangeable.

3. Load the published CPC series from `cpc_indices` with `load_regular_data`.
   - Use `product` = `ersst5_91-20` (main reference) and `ersst5_81-10` (control), with `variable` set to the Nino 3.4 anomaly column.
   - Load `oni` as well, which is the official 3-month running index.
   - Cover 1950-01 to 2023-12.

4. Compare RECON_9120 with CPC `ersst5_91-20`.
   - Use `analyze_correlation` with `detrend=False` and no lead-lag. Report the Pearson r to 6 decimals.
   - Use `compute_difference_metrics` for the RMSE (4 decimals) and the largest absolute deviation (4 decimals), including the year and month where it occurs.
   - Note [common-comparison-mask]: use the same months for all three statistics, taken where both series are defined.

5. Detect El Nino events and build the 3-month running mean.
   - Apply `temporal_aggregate` with a 3-month centred running window to RECON_9120. Use the `oni` series as published, since it is already a 3-month running index.
   - Apply `detect_threshold_events` with a threshold of +0.5 °C and `min_duration` of 5 consecutive seasons to RECON_9120's running mean, and separately to the CPC ONI series with the identical criterion.
   - Report the El Nino count for each, and the three events with the highest peaks (peak year-month and value to 2 decimals).
   - Run `detect_threshold_events` with -0.5 °C and `min_duration` 5 to find the La Nina episodes. These are needed only for the shading in Figure 1.

6. Compute the baseline offsets over 1950-01 to 2020-12.
   - Reconstructed offset: compute RECON_8110 minus RECON_9120 and take the unweighted mean over the months using `compute_basic_stat` (3 decimals).
   - Published offset: compute CPC `ersst5_81-10` minus CPC `ersst5_91-20` and take the unweighted mean over the same months (3 decimals).
   - Compare the two with `compute_difference_metrics` and state how closely they agree.

7. Figures.
   - Figure 1: two rows sharing a time axis over 1950-2023. The upper row overlays the RECON_9120 and CPC 1991-2020 monthly anomalies and the 3-month running mean, with ±0.5 °C threshold lines. It shades between the threshold and the running-mean curve during the detected El Nino episodes (warm) and La Nina episodes (cool). The lower row shows RECON_9120 minus CPC on its own vertical scale. Label both vertical axes in °C.
   - Figure 2: a square scatter of RECON_9120 (y) against CPC 1991-2020 (x), one point per month. Use identical axis ranges and ticks, a 1:1 line, and annotate the correlation and RMSE.

Other notes:
- [anomaly-baselines-and-climate-indices] The baseline offset is one constant per calendar month. It is not a reconstruction error, and the reconstructed and published offsets should agree.
- [anomaly-baselines-and-climate-indices] The ONI event count can differ from RECON_9120's. The CPC ONI uses centred 30-year base periods applied to successive 5-year spans, while RECON_9120 uses one fixed 1991-2020 baseline. Explain this in the answer.
- [serial-correlation-and-significance] The correlation is computed without detrending, since the shared trend is part of the agreement being checked. Do not quote a p-value from the naive monthly sample size, because the effective sample size is much smaller.
- Offsets cover 1950-01 to 2020-12 only. The correlation, RMSE and maximum deviation cover 1950-01 to 2023-12.
