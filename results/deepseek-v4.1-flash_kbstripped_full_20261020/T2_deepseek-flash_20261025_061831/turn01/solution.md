**Outcome.** The Nino 3.4 index was successfully rebuilt from the ERSSTv5 monthly SST field over 5°S–5°N, 170°W–120°W for 1950-01 to 2023-12 on the native 2° grid, without regridding. The reconstruction was compared with the CPC Nino 3.4 anomaly series on the 1991–2020 baseline, event statistics were computed, and the 1981–2010 versus 1991–2020 baseline offset was reproduced.

**Methodology as executed.**
- ERSSTv5 monthly SST (`dataset="ersst"`, `variable="sst"`) was loaded for 1950-01 to 2023-12, restricted to lat −5 to 5 and lon −170 to −120.
- CPC indices were loaded with `load_regular_data` using variable `NINO3.4`, `field="anomaly"`, and products `ersst5_91-20` (main reference), `ersst5_81-10` (control), and `oni` (already a 3-month running index).
- Reconstructed anomalies were produced with `compute_anomaly_vs_own_climatology` using `freq="month"`, `method="subtract"`:
  - RECON_9120: `base_period=["1991-01", "2020-12"]`
  - RECON_8110: `base_period=["1981-01", "2010-12"]`
- The anomaly fields were reduced to monthly index series with `compute_basic_stat`, `op="area_weighted_mean"`, `dim=["lat", "lon"]`.
- RECON_9120 and CPC 1991–2020 were aligned with `xr.align(..., join="inner")`.
- Correlation used `analyze_correlation` with `method="pearson"` and `detrend=False`; difference metrics used `compute_difference_metrics` with RECON_9120 as data and CPC 1991–2020 as reference.
- The 3-month centred running mean of RECON_9120 used `temporal_aggregate` with `freq="rolling"`, `window=3`, `center=True`, `op="mean"`.
- El Niño/La Niña event detection used `detect_threshold_events` with `min_duration=5`: +0.5 °C and `direction="above"` for El Niño, −0.5 °C and `direction="below"` for La Niña. The same +0.5 °C rule was applied directly to the CPC ONI.
- Baseline offsets were formed as RECON_8110 − RECON_9120 and CPC 1981–2010 − CPC 1991–2020 after common alignment, then averaged over 1950-01 to 2020-12 with `compute_basic_stat`, `op="mean"`, `dim="time"` (unweighted in time).

**Key numerical results.**
- Pearson correlation between RECON_9120 and CPC 1991–2020 over 1950-01 to 2023-12: **r = 0.999972**.
- Significance: p-value reported as **0.0** using the effective sample size, and **0.0** unadjusted. Effective sample size **n_eff = 47.41**; actual monthly sample count **n = 888**.
- RMSE: **0.0075 °C**.
- Largest absolute deviation: **0.0226 °C**, occurring in **1992-01**.
- Mean difference/reconstruction bias: **−0.0030 °C**.

**El Niño event counts and top three peaks.**
- RECON_9120, using its 3-month centred running mean and the ≥5-season rule: **16 El Niño events**.
  - **2015-12: 2.64 °C**
  - **1997-12: 2.35 °C**
  - **1983-01: 2.15 °C**
- CPC ONI, using the identical +0.5 °C and ≥5-season rule: **22 El Niño events**.
  - **2015-12: 2.59 °C**
  - **1997-12: 2.37 °C**
  - **1983-01: 2.14 °C**

The counts differ even though the detection criterion is identical. RECON_9120 uses one fixed 1991–2020 monthly climatology, whereas the CPC index carries its own base-period convention. Events whose 3-month running mean is close to the +0.5 °C threshold can therefore pass in one series and fail in the other. The difference is due to borderline events and baseline convention, not to a different threshold rule. For the Figure 1 shading, La Niña detection at −0.5 °C with the same ≥5-season rule found **25 La Niña events**.

**Baseline offsets over 1950-01 to 2020-12.**
- Reconstructed offset, mean(RECON_8110 − RECON_9120): **0.060 °C**.
- Published offset, mean(CPC 1981–2010 − CPC 1991–2020): **0.058 °C**.
- Difference: **0.002 °C** (unrounded 0.0018 °C). The reconstructed and published baseline offsets therefore agree very closely.

**Visualizations.**
- **Figure 1** is a two-row figure sharing a 1950–2023 time axis. The upper panel overlays RECON_9120 and the CPC 1991–2020 monthly anomalies, the RECON 3-month running mean, the +0.5 °C and −0.5 °C thresholds, and shading between the threshold and the running-mean curve for detected El Niño events (warm shading, ≥5 seasons at/above +0.5 °C) and La Niña events (cool shading, ≥5 seasons at/below −0.5 °C). The lower panel shows RECON_9120 minus CPC on its own vertical scale. Both vertical axes are labelled in °C.
- **Figure 2** is a square monthly scatter plot of RECON_9120 (y-axis) against CPC 1991–2020 (x-axis), one point per month. It uses identical axis ranges and tick intervals, a 1:1 reference line, and annotates the correlation and RMSE.

**Interpretation.** The reconstruction is essentially indistinguishable from the CPC 1991–2020 Nino 3.4 anomaly at monthly resolution: r exceeds 0.99997, RMSE is below 0.008 °C, and the largest single-month deviation is only about 0.023 °C. The event counts differ because of threshold-borderline events under different baseline conventions, while the baseline-offset calculation reproduces the published 1981–2010 versus 1991–2020 offset to within about 0.002 °C.
