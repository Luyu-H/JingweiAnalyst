**Method summary**
- Loaded ERSSTv5 monthly SST for 1950-01–2023-12, restricted to latitude −5° to 5° and longitude −170° to −120°, in °C.
- Built the reconstructed index by subtracting the field’s own monthly climatology over 1991-01–2020-12 (`method="subtract"`), then reducing over latitude/longitude with an area-weighted mean. This gives RECON_9120. The identical calculation with the 1981-01–2010-12 monthly climatology gives RECON_8110.
- Loaded CPC NINO3.4 anomaly products `ersst5_91-20`, `ersst5_81-10`, and `oni` for 1950-01–2023-12. The `ersst5_81-10` product returned 852 valid months (1950-01–2020-12); the baseline-offset calculation was restricted to 1950-01–2020-12 as required.
- Compared series using Pearson correlation over time and `compute_difference_metrics` with `weighting="none"` (unweighted over months). Event detection used the 3-month centered running mean for RECON_9120 (`window=3`, `center=True`) and used ONI directly, with threshold +0.5 °C or −0.5 °C and `min_duration=5` overlapping seasons. Baseline offsets are unweighted monthly means of differences.

**Reconstruction vs CPC 1991-2020**
- Pearson correlation: **0.999972** (6 d.p.; 888 valid monthly pairs).
- RMSE: **0.0075 °C** (4 d.p.).
- Largest absolute deviation: **0.0226 °C**, occurring at **1992-01**.
- Mean bias (RECON_9120 minus CPC 1991-2020): **−0.0030 °C**.

The reconstruction is extremely close to the CPC 1991-2020 series, with sub-hundredth-degree RMSE and maximum absolute deviation of about 0.023 °C.

**El Niño events and top peaks**
- RECON_9120, using its 3-month centered mean with at least five consecutive seasons ≥ +0.5 °C: **16 El Niño events**.
  - Top three peaks:
    1. **2015-12: 2.64 °C**
    2. **1997-12: 2.35 °C**
    3. **1983-01: 2.15 °C**
- CPC ONI, using the identical criterion directly on the ONI series: **22 El Niño events**.
  - Top three peaks:
    1. **2015-12: 2.59 °C**
    2. **1997-12: 2.37 °C**
    3. **1983-01: 2.14 °C**

The same three peak dates appear at the top of both lists, but the counts differ (16 vs 22) and the peak values differ slightly. The code’s explanation is that CPC’s published index uses its own product version, quality control, baseline, and spatial averaging; small differences near the +0.5 °C threshold or near the five-season duration can change whether a run is classified as an event.

For the Figure 1 shading, detection at −0.5 °C with the same five-season criterion found **25 La Niña events** in RECON_9120.

**Baseline offset check, 1950-01–2020-12**
- Reconstructed offset, RECON_8110 minus RECON_9120, unweighted monthly mean: **0.060 °C** (3 d.p.).
- Published CPC offset, CPC 1981-2010 minus CPC 1991-2020, same window: **0.058 °C**.
- Difference, reconstructed minus published: **+0.002 °C** (absolute difference **0.002 °C**).

The reconstructed and published baseline offsets agree to about **0.002 °C**; the reconstruction slightly overestimates the published baseline shift.

**Figures**
- **Figure 1** is a two-row figure sharing the 1950–2023 time axis. The upper panel overlays RECON_9120 monthly anomalies, CPC 1991-2020 monthly anomalies, and the RECON_9120 3-month centered running mean. It also shows +0.5 °C and −0.5 °C thresholds and shades between the threshold and the running-mean curve for detected El Niño events (warm shading) and La Niña events (cool shading). The lower panel plots RECON_9120 minus CPC on its own vertical scale. Both vertical axes are labeled in °C, and the lower axis is labeled by year.
- **Figure 2** is a square scatter plot of RECON_9120 (y-axis) against CPC 1991-2020 (x-axis), one point per month. It uses identical x- and y-axis ranges and ticks determined from the pooled data, includes a 1:1 reference line, and annotates the Pearson correlation, RMSE, and maximum absolute deviation.
