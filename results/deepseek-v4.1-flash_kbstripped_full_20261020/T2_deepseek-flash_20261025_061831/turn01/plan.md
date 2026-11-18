# Plan

1. Load the ERSSTv5 monthly SST field (`ersst`, variable `sst`) for 1950-01 to 2023-12 restricted to 5°S–5°N, 170°W–120°W on its native 2° grid with `load_regular_data`. In the same step load the three published CPC series with `load_regular_data`, using `product` for the series and `variable` for the region (Nino 3.4): `ersst5_91-20` (main reference), `ersst5_81-10` (control) and `oni` (official 3-month running index).

2. Turn the loaded ERSSTv5 regional SST field into an anomaly field with `compute_anomaly_vs_own_climatology`, base period 1991-2020 (the RECON_9120 anomaly field over 1950-2023).
   - Note [anomaly-baselines-and-climate-indices]: the anomaly must come from the field's own 1991-2020 monthly climatology, each month differenced against its own calendar month; record the climatology and the `base_period` the tool writes into the result attrs.

3. Reduce the anomaly field over latitude and longitude with `compute_basic_stat` (area-weighted spatial mean) to obtain RECON_9120, the monthly Nino 3.4 index series, 1950-01 to 2023-12.

4. Repeat steps 2–3 on the same ERSSTv5 regional SST field but with `compute_anomaly_vs_own_climatology` using `base_period` 1981-2010, then `compute_basic_stat` area-weighted spatial mean, to obtain RECON_8110, the same index on the 1981-2010 baseline.

5. Score RECON_9120 against the CPC 1991-2020 Nino 3.4 anomaly series over 1950-01 to 2023-12: `analyze_correlation` for the Pearson correlation to 6 decimal places, and `compute_difference_metrics` for the RMSE to 4 decimal places and the largest absolute deviation to 4 decimal places together with the year and month where it occurs.
   - Note [serial-correlation-and-significance]: correlate the reconstruction with the CPC 1991-2020 series as they are with `detrend=False` — two estimates of the same quantity whose shared trend is part of the agreement being checked — and quote any significance from the effective sample size the tools apply by default, reporting the unadjusted values alongside rather than the naive monthly step count.

6. Compute the 3-month centred running mean of RECON_9120 with `temporal_aggregate` (running window of 3), then apply `detect_threshold_events` at +0.5 °C (events that cross the threshold and stay across it) and keep only episodes lasting at least 5 consecutive overlapping 3-month seasons → the El Nino count over 1950-2023 and each event's peak date/value; take the three highest peaks and report each peak value to 2 decimal places. Apply the same `detect_threshold_events` run on the 3-month running mean at −0.5 °C with the same ≥5-season rule to obtain the La Nina events needed for the Figure 1 shading.

7. Apply the identical criterion directly to the CPC ONI series (`oni` is already the 3-month running index) with `detect_threshold_events` at +0.5 °C, keeping episodes of at least 5 consecutive seasons → the count of El Nino events and the three highest peaks with their year/month and peak values to 2 decimal places, for comparison with step 6.

8. Form the reconstruction's baseline offset as RECON_8110 minus RECON_9120 and take its unweighted mean over the months 1950-01 to 2020-12 with `compute_basic_stat` (mean in time, unweighted); form the published offset as CPC 1981-2010 minus CPC 1991-2020 Nino 3.4 anomaly and take the same unweighted monthly mean over the same months; report both to 3 decimal places and state how closely they agree. Per [anomaly-baselines-and-climate-indices], the reconstructed difference is a baseline effect (one constant offset per calendar month, variability identical), not a reconstruction error, and likewise for the published offset — both are reported as unweighted means over the months.

9. Figure 1 with own plotting code (no tool covers figures): a two-row figure sharing a time axis over 1950-2023. Upper row: RECON_9120 and the CPC 1991-2020 monthly anomalies overlaid, the 3-month running mean, the +0.5 and −0.5 °C thresholds, and shading between the threshold and the running-mean curve — warm during the step 6 El Nino events, cool during the step 6 La Nina events. Lower row: RECON_9120 minus CPC on its own vertical scale. Both vertical axes labelled in degrees Celsius.

10. Figure 2 with own plotting code: square scatter of RECON_9120 (y) against the CPC 1991-2020 series (x), one point per month, identical axis ranges and ticks, a 1:1 reference line, and the step 5 correlation and RMSE annotated.

Other notes:
- The two El Nino counts need not agree: the detection rule is identical, but RECON_9120 uses one fixed 1991-2020 monthly climatology while the CPC index carries its own base-period convention, so events whose 3-month running mean hovers at the +0.5 °C threshold can qualify in one series and not the other; the count difference therefore comes from such borderline events, not from the criterion.
- The reconstruction is built on the ERSSTv5 native 2° grid as loaded; no regridding is applied, because the comparison partners are index series rather than gridded fields.
- Report r to 6 decimal places, RMSE and largest absolute deviation to 4 decimal places (with year and month), event peaks to 2 decimal places and the two baseline offsets to 3 decimal places; intermediate values are carried unrounded.
