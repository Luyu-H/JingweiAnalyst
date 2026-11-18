The ERSSTv5 Niño 3.4 reconstruction closely reproduces the CPC 1991–2020 monthly Niño 3.4 anomaly series.

**Methodology actually used**
- ERSSTv5 monthly SST was opened for 1950-01 to 2023-12 on the native 2° grid.
- Latitude was sorted ascending and longitude converted to 0–360°. The Niño 3.4 box was selected as 5°S–5°N, 190°E–240°E, equivalent to 170°W–120°W.
- For each calendar month, a grid-cell climatology was computed over the chosen baseline: 1991–2020 for RECON_9120 and 1981–2010 for RECON_8110. The anomaly was formed by subtracting that calendar-month climatology.
- The regional monthly anomaly was computed with cosine-latitude area weights, restricted to valid non-missing cells.
- CPC Niño 3.4 anomaly columns were read from the `ersst5_91-20` and `ersst5_81-10` files. The ONI file was read by mapping each 3-month season to its centre month; the ONI series was not smoothed again.
- Comparisons used inner-joined common monthly timestamps. No lag was applied.
- Event detection used contiguous runs of at least 5 values of the centred 3-month running mean at or above +0.5°C for El Niño, and at or below −0.5°C for La Niña.

**Key numerical results**
- RECON_9120 versus CPC 1991–2020 Niño 3.4 anomaly, 1950-01 to 2023-12:
  - Pearson correlation = **0.999972**
  - RMSE = **0.0075°C**
  - Maximum absolute deviation = **0.0226°C**, occurring in **1992-01**

- El Niño events under the operational ONI criterion:
  - RECON_9120: **16 events**
  - Three largest RECON_9120 peaks:
    - **2015-12: 2.64°C**
    - **1997-12: 2.35°C**
    - **1983-01: 2.15°C**
  - CPC ONI: **22 events**
  - Three largest ONI peaks:
    - **2015-12: 2.59°C**
    - **1997-12: 2.37°C**
    - **1983-01: 2.14°C**

The same top three peak events appear in both series, but the event counts differ: RECON_9120 identifies 16 El Niño events, whereas ONI identifies 22. The executed comparison also found 172 threshold-crossing months at or above +0.5°C for RECON_9120 versus 226 for ONI. The reason is that RECON_9120 uses a fixed 1991–2020 calendar-month climatology, while the official ONI uses its operational base-period convention and is distributed rounded to 0.1°C. These choices alter marginal crossings near ±0.5°C and therefore change how many 5-season contiguous runs are counted. The ONI is already a 3-month running index, so no second smoothing was applied.

- Baseline offsets over 1950-01 to 2020-12:
  - Reconstructed offset, mean(RECON_8110 − RECON_9120) = **0.060°C**
  - Published CPC offset, mean(CPC 1981–2010 − CPC 1991–2020) = **0.058°C**
  - Published minus reconstructed difference = **−0.002°C**

Thus the reconstructed baseline offset agrees with the published CPC baseline offset to within **0.002°C**.

**Visualizations**
- Figure 1 uses two rows sharing a 1950–2023 time axis. The upper row is designed to overlay the RECON_9120 monthly anomaly, the CPC 1991–2020 monthly anomaly, the centred 3-month running mean of RECON_9120, and ±0.5°C threshold lines. It shades between the +0.5°C threshold and the running-mean curve during detected El Niño intervals, and between the −0.5°C threshold and the curve during detected La Niña intervals. The lower row shows RECON_9120 minus CPC 1991–2020 on its own vertical scale. Both vertical axes are labelled in °C.
- Figure 2 is a square scatter plot of RECON_9120 on the y-axis against CPC 1991–2020 on the x-axis, one point per common month, with identical axis ranges and ticks, a 1:1 reference line, and annotation of the Pearson correlation and RMSE. Both axes are labelled in °C.

Overall, the reconstruction is essentially the same as CPC’s 1991–2020 Niño 3.4 anomaly series, with very small numerical differences. The event-count difference relative to ONI is consistent with the different baseline and rounding conventions used for the official 3-month ONI series.
