**Direct answer.** Rebuilding Niño 3.4 from the raw ERSSTv5 field reproduces the CPC 1991–2020 Niño 3.4 anomaly series almost exactly: Pearson \(r = 0.999972\), RMSE \(= 0.0075^\circ\)C, and maximum absolute deviation \(= 0.0226^\circ\)C at **1992-01** (signed RECON minus CPC \(= -0.0226^\circ\)C). The reconstructed 1981–2010 minus 1991–2020 baseline offset is **0.060°C**, versus the published CPC offset **0.058°C**, a difference of **0.002°C**.

**Methodology from the executed code.** The code opened ERSSTv5 from `/data1/lorraine-20419/ocean_data/ERSSTv5/sst.mnmean.nc`, extracted `sst` for 1950-01-01 through 2023-12-31, and selected the Niño 3.4 region with latitude slice `5` to `-5` and longitude slice `190` to `240` on the 0–360° grid. Monthly climatologies were formed with `groupby("time.month").mean(dim="time")` over **1991–01-01 to 2020-12-31** for RECON_9120 and **1981–01-01 to 2010-12-31** for RECON_8110. Anomalies were formed by subtracting the calendar-month climatology from each month. The regional average used \(\cos(\mathrm{latitude})\) area weights and a finite-value mask, with the denominator restricted to valid cells. No lag analysis was applied; all comparisons are same-month pairings. No detrending was applied.

The CPC monthly files `ersst5.nino.mth.91-20.ascii` and `ersst5.nino.mth.81-10.ascii` were read from the `NINO3.4_ANOM` column. The CPC ONI file `oni.ascii.txt` was read from the `ANOM` column, mapped from seasonal labels to centre months (DJF→January, etc.), and then used directly. The ONI was **not smoothed a second time**: the +0.5°C threshold was applied directly to the parsed ONI series.

For event detection on RECON_9120, the monthly series was smoothed with a centred 3-month running mean using `rolling(window=3, center=True, min_periods=3)`, leaving the first and last months undefined. El Niño events were contiguous runs of at least 5 months with running-mean values \(\ge +0.5^\circ\)C; La Niña events used \(\le -0.5^\circ\)C. Event peaks are the maximum/minimum of the centred 3-month running mean inside each run, reported at the centre month. For ONI, the identical persistence rule was applied directly to the ONI values.

**Validation against CPC 1991–2020.** Over 1950-01 to 2023-12 with \(n=888\) monthly pairs:
- Pearson correlation: **0.999972**
- RMSE: **0.0075°C**
- Largest absolute deviation: **0.0226°C**
- Month of largest absolute deviation: **1992-01**
- Signed RECON minus CPC at that month: **\(-0.0226^\circ\)C** (reconstruction slightly below CPC)

**El Niño events in RECON_9120, 1950–2023.** Count: **16** events under the 5-season persistence rule. Three highest peaks:
1. **2015-12: 2.64°C**
2. **1997-12: 2.35°C**
3. **1983-01: 2.15°C**

The code also detected **25** La Niña events for the shading in Figure 1.

**El Niño events in CPC ONI, 1950–2023.** Count: **22** events, applying the threshold directly to the ONI series. Three highest peaks:
1. **2015-12: 2.59°C**
2. **1997-12: 2.37°C**
3. **1983-01: 2.14°C**

The counts differ because RECON_9120 uses a fixed 1991–2020 monthly climatology and its monthly anomalies were further smoothed with a centred 3-month running mean before applying the persistence criterion. CPC ONI is already the official 3-month running mean and uses CPC’s rolling 30-year base periods, updated every five years, rather than a fixed 1991–2020 baseline. These baseline and smoothing-definition differences change how many 3-month seasons cross \(\pm 0.5^\circ\)C, so the event counts and peak values need not match exactly.

**Baseline offsets over 1950-01 to 2020-12.** The reconstructed offset, RECON_8110 − RECON_9120, has an unweighted monthly mean of **0.060°C**. The published CPC offset, CPC 1981–2010 − CPC 1991–2020, has an unweighted monthly mean of **0.058°C**. Their difference is **0.002°C**. The offsets are not exactly constant month to month; the month-by-month offset difference has mean **0.0018°C**, standard deviation **0.0042°C**, maximum absolute difference **0.0283°C**, and correlation **0.996445**. This month-dependent offset reflects the difference between the 1981–2010 and 1991–2020 climatological baselines, not an error in either series.

**Figures.**  
- **Figure 1** is a two-row figure sharing the 1950–2024 time axis. The upper row overlays the monthly RECON_9120 and CPC 1991–2020 anomalies, their centred 3-month running means, and the +0.5/−0.5°C thresholds, with red shading between the +0.5°C threshold and the RECON running mean during detected RECON El Niño events and blue shading between the −0.5°C threshold and the RECON running mean during detected La Niña events. The lower row plots RECON_9120 minus CPC on its own vertical scale, with a zero line. Both vertical axes are labelled in °C.  
- **Figure 2** is a square scatter plot of monthly RECON_9120 (y) against CPC 1991–2020 (x), with identical axis limits and ticks, a 1:1 reference line, and annotation of the step-5 correlation, RMSE, and number of paired months.

**Interpretation.** The reconstruction is numerically almost indistinguishable from the published CPC 1991–2020 monthly Niño 3.4 anomaly at the monthly-pair level. The baseline-offset calculation also closely reproduces CPC’s published 1981–2010 versus 1991–2020 offset. The event-count difference is therefore attributable to the baseline/smoothing conventions used for event definition, not to a failure of the SST reconstruction.
