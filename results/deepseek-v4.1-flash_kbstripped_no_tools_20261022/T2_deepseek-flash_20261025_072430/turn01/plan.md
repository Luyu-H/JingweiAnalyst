# Plan

1. Open `ersst` monthly `sst` files for 1950-01 to 2023-12 on the native 2° grid; open `cpc_indices` products `ersst5_91-20`, `ersst5_81-10`, and `oni`, selecting the Niño 3.4 anomaly column from each CPC Niño product. Normalize time to monthly timestamps and handle the ERSST longitude convention so 170W-120W is included. No CMIP models are used.

2. Build the reconstructed Niño 3.4 indices with one function applied twice. For each month, select the Niño 3.4 box 5S-5N, 170W-120W; for each calendar month compute the grid-cell SST climatology over 1991-2020 for RECON_9120 and over 1981-2010 for RECON_8110; subtract the corresponding calendar-month climatology to form anomalies; then area-weight the anomaly field over the box to a monthly regional series for 1950-01 to 2023-12. Store baseline and climatology details in the result attrs.
- Note [anomaly-baselines-and-climate-indices]: Build RECON_9120 as SST minus the climatology of its own calendar month over 1991-2020, and build RECON_8110 the same way with the 1981-2010 baseline; record climatology and baseline in the result attrs; later count events from the running mean by persistence, not by eye.
- Note [omip-weighted-means-and-integrals]: The Niño 3.4 regional mean must be area-weighted using cosine-latitude cell areas on the regular 2° grid, with weights restricted to the same valid non-missing cells; report the area mean itself, not an unweighted cell count as an area share.

3. Align RECON_9120 with the CPC `ersst5_91-20` Niño 3.4 anomaly over the common months 1950-01 to 2023-12. Compute Pearson correlation and report it to 6 decimals, RMSE = sqrt(mean((RECON minus CPC)^2)) and report it to 4 decimals, and the maximum absolute deviation with its year and month and report it to 4 decimals. Keep full precision until reporting.

4. Compute RECON_8110 minus RECON_9120 for every month in 1950-01 to 2020-12 and take the unweighted mean over those months; report to 3 decimals. This is a time mean of monthly regional series, not an additional area weighting.

5. Detect El Niño and La Niña events in RECON_9120 using the persistence rule from [anomaly-baselines-and-climate-indices]: form the 3-month centred running mean of the monthly anomaly, classify a warm event as a contiguous run of at least 5 monthly running-mean values at or above +0.5°C and a cold event as a contiguous run of at least 5 values at or below -0.5°C, and use a `detect_threshold_events`-style contiguous-run check rather than visual inspection. Count El Niño events over 1950-2023; for each El Niño event record its peak running-mean value and the month/year of that peak; sort by peak value and report the three largest peaks to 2 decimals. Keep the La Niña intervals for Figure 1 shading.

6. Apply the identical event criterion to the official CPC `oni` series, which is already a 3-month running index: do not smooth it again. Count El Niño events over 1950-2023 and report the three largest peak values with their peak month/year to 2 decimals. Compare the RECON_9120 and ONI event sets and threshold-crossing months to explain the count difference, noting that the official ONI uses its operational base-period convention rather than the fixed 1991-2020 climatology used for RECON_9120, with resulting marginal threshold changes.

7. Compute the published baseline offset: align the CPC `ersst5_81-10` and `ersst5_91-20` Niño 3.4 anomaly columns over 1950-01 to 2020-12, take the unweighted mean of CPC 1981-2010 minus CPC 1991-2020, and report it to 3 decimals. Compare it with the RECON_8110 minus RECON_9120 offset from step 4; report both means and their difference as the measure of agreement.

8. Figure 1: two rows sharing a 1950-2023 time axis. Upper row: plot RECON_9120 and CPC 1991-2020 monthly anomalies, the 3-month centred running mean of RECON_9120, and +0.5 and -0.5°C threshold lines; shade between the +0.5°C threshold and the running-mean curve for RECON_9120 El Niño event intervals, and between the -0.5°C threshold and the curve for detected La Niña intervals. Lower row: RECON_9120 minus CPC 1991-2020 on its own vertical scale. Label both vertical axes in °C.

9. Figure 2: square scatter of RECON_9120 on the y-axis against CPC 1991-2020 on the x-axis, one point per common month; use identical x and y ranges and ticks, add the 1:1 reference line, and annotate the Pearson correlation and RMSE from step 3. Label both axes in °C.

Other notes:
- The reconstructed index uses a fixed 1991-2020 monthly climatology; the official ONI uses its operational base-period convention, so marginal events near the ±0.5°C threshold can be counted differently.
- CPC `oni` is already a 3-month running series; applying a second 3-month running mean would be wrong.
- Baseline offsets are unweighted means over months, while the Niño 3.4 regional means themselves are area-weighted over valid grid cells.
- Report only final metrics at the requested precision; use unrounded intermediate values in all calculations.
- The comparison must use common monthly timestamps between the ERSSTv5 reconstruction and each CPC product.
