**Methodology (as executed).** ERSSTv5 monthly `sst` was loaded for 1950-01–2023-12 in the latitude range [-6, 6] and longitude range [-172, -118]. CPC `cpc_indices` NINO3.4 series were loaded for products `ersst5_91-20`, `ersst5_81-10`, and `oni` over the same period. RECON_9120 and RECON_8110 were built with `compute_enso_nino_index` using `nino_region="nino3.4"`, `anomaly_freq="month"`, `climatology="self"`, and base periods 1991-01–2020-12 and 1981-01–2010-12 respectively. Monthly series were aligned by inner join and `dropna`; comparison used 1950-01–2023-12, and baseline offsets used 1950-01–2020-12. Pearson correlation used `analyze_correlation(method="pearson")`; RMSE and maximum absolute deviation used `compute_difference_metrics`; the 3-month running mean used `temporal_aggregate(freq="rolling", window=3, center=True, op="mean")`. Events were detected with `detect_threshold_events` at +0.5 °C, direction above, minimum duration 5 overlapping seasons; La Niña shading used -0.5 °C, direction below, same duration. CPC ONI events were detected directly on the loaded ONI series with the same +0.5 °C / 5-season rule. Baseline offsets were computed from monthly differences and averaged with `compute_basic_stat(op="mean")`.

**RECON_9120 vs CPC 1991–2020 Niño 3.4, 1950-01–2023-12 (N=888):**
- Pearson correlation: **0.999972** (6 dp)
- RMSE: **0.0075 °C** (4 dp)
- Largest absolute deviation: **0.0226 °C** at **1992-01** (4 dp)

The reconstruction is essentially identical to CPC’s published 1991–2020 Niño 3.4 anomaly. The residual RMSE is below 0.01 °C, and the largest monthly departure is only about 0.023 °C, consistent with official processing/rounding differences.

**El Niño events in RECON_9120** using the 3-month centered running mean and at least 5 consecutive seasons ≥ +0.5 °C:
- Count: **16**
- Three highest peaks:
  - **2015-12: 2.64 °C**
  - **1997-12: 2.35 °C**
  - **1983-01: 2.15 °C**

For context, the same detection with threshold -0.5 °C found **25** La Niña events for shading.

**El Niño events in the CPC ONI** using the identical +0.5 °C / ≥5-season rule:
- Count: **22**
- Three highest peaks:
  - **2015-12: 2.59 °C**
  - **1997-12: 2.37 °C**
  - **1983-01: 2.14 °C**

The event counts differ (16 for RECON_9120 vs 22 for ONI) because the operational ONI is already a 3-month running mean and uses NOAA’s operational evolving 30-year base periods updated every five years, together with official processing/rounding. RECON_9120 instead uses one fixed 1991–2020 monthly climatology. Events near the ±0.5 °C threshold can therefore be classified differently.

**Baseline offsets, mean over 1950-01–2020-12 (N=852):**
- Reconstructed offset, RECON_8110 − RECON_9120: **0.060 °C** (3 dp)
- Published CPC offset, `ersst5_81-10` − `ersst5_91-20`: **0.058 °C** (3 dp)
- CPC minus reconstructed offset: **−0.002 °C** (3 dp)

Thus the reconstructed and published baseline offsets agree to within about **0.002 °C**.

**Figures.** Figure 1 is a two-panel time-series plot sharing the 1950–2023 axis. The upper panel overlays monthly RECON_9120 and CPC 1991–2020 Niño 3.4 anomalies, the RECON 3-month centered running mean, ±0.5 °C threshold lines, and warm/cool shading between the threshold and the running-mean curve for detected El Niño/La Niña episodes. The lower panel shows RECON_9120 minus CPC on its own vertical scale; both y-axes are labeled in °C. Figure 2 is a square monthly scatter of RECON_9120 (y) against CPC 1991–2020 (x), with identical axis ranges and ticks, a 1:1 reference line, and annotated Pearson correlation and RMSE.
