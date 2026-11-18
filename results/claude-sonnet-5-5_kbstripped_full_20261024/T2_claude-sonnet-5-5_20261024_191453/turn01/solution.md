## Direct answer

The Niño 3.4 index rebuilt from ERSSTv5 reproduces the CPC 1991–2020 series almost exactly (r = 0.999972, RMSE 0.0075 °C). The same code also reproduces the CPC offset between the 1981–2010 and 1991–2020 baselines to within about 0.002 °C. The El Niño count differs between the reconstruction (16) and the CPC ONI (22), while the three largest peaks agree closely.

## Methodology (as executed)

- **Data.** `load_regular_data` loaded ERSSTv5 `sst` (`ersst`) for 1950-01 to 2023-12 on the native 2° grid. The load box was 10S–10N, 175W–115W, giving 888 months × 11 lat × 30 lon points, all valid.
- **Index.** `compute_enso_nino_index` was called with `nino_region="nino3.4"`, `anomaly_freq="month"` and `climatology="self"`. Anomalies are therefore relative to ERSSTv5's own per-calendar-month climatology over the stated base period.
  - RECON_9120 used `base_period` 1991-01 to 2020-12.
  - RECON_8110 used `base_period` 1981-01 to 2010-12.
  - Both cover 1950-01 to 2023-12. I did not inspect the function's internals, so I make no claim about how it weights grid cells when averaging over the region.
- **Reference series.** CPC series were loaded from `cpc_indices` with variable `NINO3.4` and `field="anomaly"`, for products `ersst5_91-20`, `ersst5_81-10` and `oni`.
  - Time stamps were normalised to the first of the month.
  - The 1981–2010 product has only 852 months (through 2020-12).
- **Comparison.** RECON_9120 and CPC 1991–2020 were aligned with an inner join and restricted to months where both are valid. All 888 months were retained, and the same months were used for the correlation, RMSE and maximum deviation.
  - `analyze_correlation` ran with `detrend=False` and no lead-lag.
  - `compute_difference_metrics` is unweighted and its difference is RECON minus CPC.
- **Events.** `temporal_aggregate` with `freq="rolling"`, `window=3`, `center=True` gave the RECON_9120 3-month running mean, which has NaN at the first and last month (886 valid values).
  - The CPC `oni` series was used as published.
  - `detect_threshold_events` ran with +0.5 °C, `direction="above"` and `min_duration=5` on both series.
  - La Niña episodes were detected on the RECON running mean at −0.5 °C, `direction="below"`, `min_duration=5`, and used only for shading.
- **Offsets.** RECON_8110 − RECON_9120 and CPC 81–10 − CPC 91–20 were computed over 1950-01 to 2020-12 (852 months). Their unweighted means came from `compute_basic_stat`, and the two series were compared with `compute_difference_metrics`.
- **Not done.** RECON_8110 was not compared directly with CPC `ersst5_81-10` in the executed code. That CPC series was used only to form the published offset.

## Key results

**RECON_9120 vs CPC 1991–2020 (1950-01 to 2023-12, n = 888)**

| Statistic | Value |
|---|---|
| Pearson r | 0.999972 |
| RMSE | 0.0075 °C |
| Mean bias (RECON − CPC) | −0.0030 °C |
| Largest absolute deviation | 0.0226 °C, at 1992-01 (RECON lower than CPC) |

- No p-value is quoted. The effective sample size from the autocorrelation adjustment is about 47 against 888 nominal months, so a naive monthly p-value would be meaningless.
- The residuals are small but somewhat larger than CPC's two-decimal rounding alone would produce (maximum 0.0226 °C against a rounding limit of about 0.005 °C). The code did not investigate the source of the remaining differences.

**El Niño events (≥ +0.5 °C for at least 5 consecutive seasons, 1950–2023)**

| Series | Count | Three highest peaks |
|---|---|---|
| RECON_9120 (3-month centred running mean) | 16 | 2015-12 (2.64 °C), 1997-12 (2.35 °C), 1983-01 (2.15 °C) |
| CPC ONI (as published) | 22 | 2015-12 (2.59 °C), 1997-12 (2.37 °C), 1983-01 (2.14 °C) |

- The same three events top both lists, with peak values within about 0.05 °C.
- The printed event tables are truncated in the output. In the visible portion, the ONI list contains events that do not appear in the RECON list, such as 1951–52, 1953, 1968–69 and 1976–77. It also splits the late-1970s warm period into two events (1976–77 and 1977–78).
- **Why the counts differ.** The CPC ONI uses centred 30-year base periods applied to successive 5-year spans, so its baseline changes through the record. RECON_9120 uses one fixed 1991–2020 baseline.
  - A moving baseline shifts values across the ±0.5 °C threshold and can lengthen, shorten, merge or split episodes. This matters most in the early decades, where the two baselines differ most.
  - The mean of the published ONI over the window (+0.013 °C) is much higher than the mean of RECON_9120 (−0.185 °C), which is consistent with this.
  - This is the stated explanation. The code did not run a controlled experiment to isolate the baseline effect from other differences, such as CPC's rounding to two decimals.

**Baseline offset (1981–2010 minus 1991–2020), 1950-01 to 2020-12, unweighted mean over 852 months**

| Quantity | Value |
|---|---|
| Reconstructed offset, RECON_8110 − RECON_9120 | 0.060 °C |
| Published offset, CPC 81–10 − CPC 91–20 | 0.058 °C |
| Difference of the two means | 0.002 °C |
| RMSE between the two monthly offset series | 0.0046 °C |
| Maximum absolute difference | 0.0283 °C, at 2020-12 |

- The reconstructed offset ranges from about −0.007 to +0.137 °C across months. It is a calendar-month-dependent constant, not a reconstruction error, and the reconstruction reproduces the published value closely.
- The published values are rounded to 0.01 °C, which limits agreement. The 2020-12 maximum difference exceeds that rounding limit, and the code did not investigate its cause.

## Figures

- **Figure 1** has two rows sharing a time axis from 1950 to 2023.
  - The upper row overlays the CPC 1991–2020 and RECON_9120 monthly anomalies. It also shows the RECON_9120 3-month centred running mean, dashed lines at ±0.5 °C, and a zero line. Between the threshold and the running-mean curve it shades red for the detected El Niño episodes and blue for the detected La Niña episodes, both from the RECON running mean. The CPC ONI is not drawn.
  - The lower row plots RECON_9120 − CPC on its own vertical scale.
  - Both y-axes are labelled in °C.
- **Figure 2** is a square scatter plot with one point per month. CPC 1991–2020 is on the x-axis and RECON_9120 on the y-axis. Both axes share the same range (rounded outward to 0.5 °C) and the same ticks, with equal aspect and a dashed 1:1 line. The correlation, RMSE and n are annotated.

## Insights

- Over the full record the reconstruction tracks the published CPC 1991–2020 series to within about 0.02 °C, so the choice of reconstruction method is not a meaningful source of difference.
- The difference between the 1981–2010 and 1991–2020 products is a small, nearly constant baseline shift of about 0.06 °C. This is much smaller than ENSO variability, with a standard deviation of about 0.87 °C.
- The El Niño event count is sensitive to the baseline convention, even though the strongest events are insensitive to it.
