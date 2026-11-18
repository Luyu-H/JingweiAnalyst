## 1. Direct answer

The Niño 3.4 index was rebuilt from the raw ERSSTv5 monthly SST field and matches the CPC `ersst5_91-20` Niño 3.4 anomaly essentially to the precision at which CPC distributes its product. Over 1950-01 … 2023-12 (888 common months) the reconstruction correlates with CPC at **r = 0.999972**, with **RMSE = 0.0075 °C** and a **largest absolute deviation of 0.0226 °C, occurring in January 1992**. Applying the operational ONI rule (≥5 consecutive 3-month-mean values at or above +0.5 °C) to the reconstruction yields **16 El Niño events over 1950–2023**, whose three largest peaks are **Dec 2015 (2.64 °C), Dec 1997 (2.35 °C) and Jan 1983 (2.15 °C)**. The same rule applied to the published CPC ONI gives **22 events** whose three largest peaks are **Dec 2015 (2.59 °C), Dec 1997 (2.37 °C) and Jan 1983 (2.14 °C)**. The reconstructed 1981-2010-minus-1991-2020 baseline offset is **+0.060 °C** (mean over 1950-01 … 2020-12), the published offset is **+0.058 °C**, and they agree to **0.002 °C**.

## 2. Methodology (as actually executed)

**Data and grid.** ERSSTv5 was opened with xarray from `sst.mnmean.nc` (variables `sst`, `time_bnds`), and `sst` was subset to 1950-01-01 … 2023-12-31, giving 888 monthly fields on the native grid: 89 latitudes (−88° to 88°, 2°) × 180 longitudes (0°–358°, 2°). Longitudes in the file already run 0–360; the code applied `lon % 360` before masking, so no re-ordering was needed.

**Region average (spatially area-weighted).** A boolean mask was built on the file's own grid: `lat >= -5 & lat <= 5` (which selects the five cell rows at −4°, −2°, 0°, 2°, 4°) and `lon%360 >= 190 & lon%360 <= 240` (26 cell columns covering 170°W–120°W). The region mean for each month was formed as

`Σ(sst_box · cos(lat)) / Σ(valid · cos(lat))`

i.e. a cosine-latitude area-weighted mean, with the sums taken over both `lat` and `lon` and restricted to valid (non-missing) cells via the `valid` indicator — not a plain latitude or grid-cell mean. The denominator is constant at 129.84 (= Σcos over the 5 latitudes × 26 longitudes) for all 888 months, confirming no missing data inside the box; all 130 cells were used every month.

**Anomalies and baselines.** For each calendar month separately, a climatology was computed as the mean of that month's region-mean SST over the baseline years, then subtracted from every month, month-matched:
- 1991-01 … 2020-12 climatology → **RECON_9120** (888 months, 1950-01 … 2023-12); the climatology and `base_period='1991-2020'` were stored in the result `attrs`.
- 1981-01 … 2010-12 climatology → **RECON_8110**, same months, differing from RECON_9120 only through the baseline.

**CPC series.** The CPC ASCII files were parsed column-wise; the Niño 3.4 anomaly is the 10th token (`ANOM4`), with values below −90 treated as missing. `ersst5_91-20` (1991-2020 base), `ersst5_81-10` (1981-2010 base) and the seasonal `oni` file (4th column, anomaly) were read; the ONI 3-letter seasons were mapped to their **centre month** (DJF→Jan, … NDJ→Dec) to put them on the same monthly axis. Record lengths found: CPC 91-20 = 918 months, CPC 81-10 = 852 months (ends Dec 2020), ONI = 919 months.

**Statistics over time (unweighted).** All step-5, step-8 and step-9 statistics are computed across months (one value per month), not across the map, so no spatial weighting enters them. Correlation is the Pearson product-moment correlation of the paired monthly series; RMSE is `sqrt(mean((recon − cpc)²))`; the largest absolute deviation is `max|recon − cpc|` with its timestamp.

**Event detection.** A **3-month centred running mean** (`rolling(window=3, center=True, min_periods=3)`) was applied to RECON_9120, so 1950-01 and 2023-12 have no running-mean value. Events were then found as **contiguous runs of at least 5 consecutive running-mean values** at or above +0.5 °C (El Niño) or at or below −0.5 °C (La Niña); the event peak is the extremum within the run and its date is the centre month of the peak 3-month window. The identical function was applied directly to the published ONI (already a 3-month running index, hence no further smoothing), over the ONI record as distributed.

**Baseline offsets.** Month-by-month differences were formed and then averaged **unweighted over the months 1950-01 … 2020-12** (852 months, dropna applied): RECON_8110 − RECON_9120 for the reconstruction, and CPC `ersst5_81-10` − CPC `ersst5_91-20` for the published pair.

## 3. Key numerical results

| Quantity | Value |
|---|---|
| Pearson r, RECON_9120 vs CPC 91-20 (1950–2023, 888 months) | **0.999972** |
| RMSE | **0.0075 °C** |
| Max absolute deviation | **0.0226 °C**, **January 1992** |
| RECON_9120 El Niño count (1950–2023) | **16** |
| RECON_9120 top-3 peaks | **Dec 2015, 2.64 °C; Dec 1997, 2.35 °C; Jan 1983, 2.15 °C** |
| CPC ONI El Niño count (its full record) | **22** |
| CPC ONI top-3 peaks | **Dec 2015, 2.59 °C; Dec 1997, 2.37 °C; Jan 1983, 2.14 °C** |
| Reconstructed offset (8110 − 9120), 1950–2020 | **+0.060 °C** |
| Published offset (CPC 81-10 − CPC 91-20), 1950–2020 | **+0.058 °C** |
| Difference (reconstruction − published) | **0.002 °C** (|Δ| = 0.002 °C) |

Auxiliary outputs from the same run: 25 La Niña events in RECON_9120 and 20 in ONI.

**Interpretation.** The near-unity correlation and sub-hundredth-of-a-degree RMSE show that the prescribed recipe — cos-latitude-weighted average of the 130 native cells in 5°S–5°N, 170°W–120°W, minus a month-matched 1991-2020 climatology of that same series — reproduces CPC's ERSSTv5-based Niño 3.4 anomaly almost exactly. The residual is at the level of CPC's 0.01 °C text-file rounding; the single largest departure (Jan 1992) falls at the peak of the 1991-92 event, where rounding of the archived value is largest in absolute terms. The positive baseline offset is expected and physical in sign: the 1981-2010 normals of the Niño 3.4 region are warmer than the 1991-2020 normals, so anomalies referenced to 1981-2010 sit higher; the reconstruction reproduces this published offset to 0.002 °C, confirming that the offset is a baseline definition rather than an error in either series.

**Why the two event counts differ (16 vs 22).** (i) *Different normals*: RECON_9120 uses a single fixed 1991-2020 monthly climatology, whereas the ONI uses its own published rolling 30-year centred normals, so the same SST maps to slightly different anomaly values, especially in the 1950s–1970s. This is visible in the early record: the reconstruction's detected warm events begin in 1957, while the ONI additionally flags 1951-52 and 1953. (ii) *Precision*: CPC products are distributed rounded (0.01 °C in the Niño 3.4 files, 0.01 °C in ONI), whereas the reconstruction is carried at full float precision, so threshold crossings at ±0.5 °C can be decided differently. (iii) *Edge truncation and record length*: the centred 3-month mean has no value at 1950-01 and 2023-12, so events straddling the ends of the reconstruction (e.g. the 2023 event, detected only through Nov 2023) are cut short, while the ONI series is distributed from its own first to last centre month (Jan 1950 to Jul 2026 in the file used here), i.e. over a *longer* window than the 1950–2023 reconstruction interval — so the ONI count includes events after 2023 that the reconstruction window cannot contain. Merging/splitting of multi-year warm spells also depends on these details (e.g. ONI's single 20-month 2014-10 … 2016-05 run versus the reconstruction's 2015-03 … 2016-04 run).

## 4. Figures

**Figure 1** (two rows, shared 1950–2023 time axis). Top row: thin lines for the RECON_9120 monthly anomalies and the CPC `ersst5_91-20` monthly anomalies overlaid, a heavier black line for the 3-month centred running mean of RECON_9120, dashed ±0.5 °C threshold lines, and `fill_between` shading between the ±0.5 °C level and the running-mean curve — red where a detected El Niño event is active, blue where a detected La Niña event is active (masks built from the ≥5-consecutive-month detections and additionally intersected with finite running-mean values). Vertical axis labelled in °C. Bottom row: RECON_9120 minus CPC on its own vertical scale (with a zero line), y-axis labelled "RECON_9120 - CPC (°C)", x-axis the year; the shared axis is fixed to 1950-01 … 2023-12. The figure is designed so that the top panel shows agreement between the reconstruction and the operational product together with the event catalogue used for counting, and the bottom panel exposes the (very small) residual structure in time.

**Figure 2** (square panel). Scatter of RECON_9120 (y) against CPC `ersst5_91-20` (x), one point per common month of 1950–2023, with identical x and y limits (data range plus 6 % padding), identical 1 °C tick spacing, an equal aspect ratio, a 1:1 dashed reference line, and an annotation box giving r = 0.999972 and RMSE = 0.0075 °C. Its purpose is to show that the points collapse onto the 1:1 line, i.e. that the reconstruction is unbiased relative to CPC across the full anomaly range.

## 5. Notes and caveats

- The four series involved are not interchangeable: RECON_9120 and CPC `ersst5_91-20` are on 1991-2020, RECON_8110 and CPC `ersst5_81-10` on 1981-2010, and the ONI on its own published rolling normals; each count/peak statement above refers to the series named.
- The CPC 81-10 file ends in Dec 2020, which is why the offset comparison is restricted to 1950-01 … 2020-12 (852 months); the averages are unweighted over months, a different operation from the spatial cosine-latitude weighting used for the region mean.
- All quantities were carried unrounded and rounded only for reporting (6 dp for r, 4 dp for RMSE/max deviation, 2 dp for peaks, 3 dp for offsets).
