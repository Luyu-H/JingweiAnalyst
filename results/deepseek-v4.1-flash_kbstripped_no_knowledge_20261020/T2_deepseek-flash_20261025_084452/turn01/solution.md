## Summary

The Nino 3.4 index was rebuilt from the raw ERSSTv5 SST field and compared with the CPC published Nino 3.4 indices. The reconstruction is extremely close to CPC on the 1991–2020 baseline.

### Methodology used

- **Raw SST:** `load_regular_data` was called with `dataset="ersst"`, `variable="sst"`, `field="total"`, `time_range=["1950-01","2023-12"]`, `lat_range=[-5.0,5.0]`, and `lon_range=[-170.0,-120.0]`.
- **Regional mean:** `compute_basic_stat` with `op="area_weighted_mean"` and `dim=["lat","lon"]` produced the monthly area-weighted Nino 3.4 SST series.
- **Climatologies:** `temporal_aggregate` with `freq="climatology"` and `op="mean"` was applied to the regional-mean SST over **1991-01 to 2020-12** and **1981-01 to 2010-12**, giving 12 monthly climatological values for each baseline.
- **Reconstructed anomalies:** `RECON_9120` and `RECON_8110` were formed by subtracting the matching calendar-month climatology from the full 1950–2023 regional-mean SST series using `groupby("time.month")`. No intermediate rounding was applied.
- **CPC indices:** `load_regular_data` loaded `dataset="cpc_indices"`, `variable="NINO3.4"` for products `"ersst5_91-20"`, `"ersst5_81-10"`, and `"oni"`, all over 1950–2023.
- **Comparison:** `RECON_9120` and CPC 1991–2020 were inner-aligned over their common months. Pearson correlation used `analyze_correlation(..., dim="time", method="pearson")`. Difference statistics used `compute_difference_metrics(..., dim="time", weighting="none")`. The largest absolute deviation was taken from the absolute difference field, not the signed maximum. No lag analysis was performed; all comparisons are contemporaneous monthly values.
- **Event detection:** The 3-month centred running mean of `RECON_9120` was computed with `temporal_aggregate(freq="rolling", window=3, center=True, op="mean")`. Events were detected with `detect_threshold_events` at `+0.5 °C` for El Niño and `-0.5 °C` for La Niña, with `min_duration=5`. The same El Niño detector was applied directly to CPC ONI.
- **Baseline offsets:** The reconstructed offset was the unweighted time mean of `RECON_8110 - RECON_9120` over 1950-01 to 2020-12 using `compute_basic_stat(op="mean", dim="time")`. The published offset used CPC 1981–2010 minus CPC 1991–2020 over the same period after inner alignment.

### Key numerical results

- **Pearson correlation** between `RECON_9120` and CPC 1991–2020 Nino 3.4 anomaly, 1950-01 to 2023-12: **0.999972**
- **RMSE:** **0.0075 °C**
- **Largest absolute deviation:** **0.0226 °C**, occurring in **January 1992**
- The mean difference (bias, not requested but computed) was **-0.0030 °C**, confirming a very small overall offset.

### El Niño counts and largest peaks

Using the operational ONI criterion — at least 5 consecutive overlapping 3-month seasons with the centred running mean at or above +0.5 °C:

- **RECON_9120 El Niño count:** **16**
- Three largest RECON_9120 peaks:
  1. **December 2015**, **2.64 °C**
  2. **December 1997**, **2.35 °C**
  3. **January 1983**, **2.15 °C**

- **CPC ONI El Niño count:** **22**
- Three largest CPC ONI peaks:
  1. **December 2015**, **2.59 °C**
  2. **December 1997**, **2.37 °C**
  3. **January 1983**, **2.14 °C**

The same three strongest events appear in both series, but the counts differ. CPC ONI is the official 3-month running mean of CPC’s own ERSSTv5 Nino 3.4 anomaly product, computed with CPC’s base-period and processing choices. `RECON_9120` is rebuilt here from raw ERSSTv5 SST using an area-weighted average over 5°S–5°N, 170°W–120°W and a fixed 1991–2020 monthly climatology. Even under the same +0.5 °C and 5-season rule, small differences in base period, spatial averaging, and monthly anomaly construction can move a marginal event across the threshold or change its duration.

### Baseline offsets

- **Reconstructed offset** (`RECON_8110 - RECON_9120`), unweighted mean over 1950-01 to 2020-12: **0.060 °C**
- **Published CPC offset** (CPC 1981–2010 minus CPC 1991–2020), same months: **0.058 °C**
- **Difference:** **0.002 °C**

The reconstructed and published baseline offsets agree to within **0.002 °C**. The positive sign means the older 1981–2010 baseline is cooler than the newer 1991–2020 baseline in this comparison.

### Visualizations

- **Figure 1:** A two-row figure sharing a 1950–2023 time axis. The upper row overlays `RECON_9120`, CPC 1991–2020 monthly anomalies, the `RECON_9120` 3-month centred running mean, and the +0.5 °C and −0.5 °C thresholds. Warm shading fills between the +0.5 °C threshold and the running mean during detected El Niño events; cool shading fills between the −0.5 °C threshold and the running mean during detected La Niña events. The lower row shows `RECON_9120 - CPC_9120` on its own vertical scale. Both vertical axes are labelled in °C. The La Niña detector found 25 events for shading.
- **Figure 2:** A square scatter plot of monthly `RECON_9120` (y) against CPC 1991–2020 Nino 3.4 anomaly (x), with identical axis ranges and ticks, a 1:1 reference line, and annotation of the Pearson correlation and RMSE.

### Interpretation

The reconstruction reproduces the CPC 1991–2020 Nino 3.4 anomaly almost exactly at monthly resolution, with correlation **0.999972**, RMSE **0.0075 °C**, and maximum absolute deviation only **0.0226 °C**. The event analysis shows the same three strongest El Niño peaks in both `RECON_9120` and CPC ONI, but the total count differs (**16 vs 22**) because of CPC’s own ONI base-period and processing choices versus the fixed-baseline, area-weighted reconstruction used here. The mean baseline offset is also reproduced closely: **0.060 °C** reconstructed versus **0.058 °C** published, a difference of only **0.002 °C**.
