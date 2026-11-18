The Nino 3.4 index was reconstructed from ERSSTv5 by averaging SST over 5°S–5°N, 170°W–120°W and subtracting the 1991–2020 monthly climatology.  The reconstruction (RECON_9120) was then compared with the CPC Nino 3.4 anomaly on the same 1991–2020 baseline over 1950-01 to 2023-12.

**Monthly validation against CPC 1991–2020 Nino 3.4**
- Pearson correlation: **0.999971** (to 6 decimal places)
- RMSE: **0.0075 °C** (to 4 decimal places)
- Largest absolute deviation: **0.0229 °C**, occurring in **1992-01**

This is an essentially exact reproduction of the CPC monthly Nino 3.4 anomaly, with residuals at the level of a few hundredths of a degree Celsius.

**El Nino events under the ONI criterion**

Using at least five consecutive overlapping 3-month seasons with the 3-month running mean ≥ +0.5 °C:

- RECON_9120: **16 El Nino events** over 1950–2023.
  - Three highest peaks:
    1. **2015-12, 2.64 °C**
    2. **1997-12, 2.35 °C**
    3. **1983-01, 2.15 °C**

- CPC ONI: **22 El Nino events** over the same period.
  - Three highest peaks:
    1. **2015-12, 2.59 °C**
    2. **1997-12, 2.37 °C**
    3. **1983-01, 2.14 °C**

The two series identify the same three strongest events, with small peak-amplitude differences.  The event-count difference arises because RECON_9120 uses a **fixed 1991–2020 monthly climatology** and unrounded ERSSTv5 area means, whereas the official CPC ONI uses its operational processing: a **rolling 30-year base period updated every five years**, CPC quality-control/missing-data treatment, and values rounded to 0.1 °C.  Near the ±0.5 °C threshold, those small baseline and rounding differences can lengthen, shorten, split, or merge marginal 3-month runs, producing 16 detected events in RECON_9120 versus 22 in the official ONI series.

**Baseline offset between 1981–2010 and 1991–2020**

Mean over 1950-01 to 2020-12 of (1981–2010 anomaly − 1991–2020 anomaly):

- Reconstructed offset, RECON_8110 − RECON_9120: **0.060 °C**
- Published CPC offset, CPC 1981–2010 − CPC 1991–2020: **0.058 °C**

The two offsets agree to **0.0018 °C** in absolute value; i.e. the reconstructed baseline offset is within about **0.002 °C** of the published CPC offset.  The positive sign means that anomalies referenced to the cooler 1981–2010 climatology are slightly higher, on average, than those referenced to 1991–2020.

**Figures generated**
- **Figure 1**: two-panel time-series figure for 1950–2023.  The upper panel overlays RECON_9120 and CPC 1991–2020 monthly anomalies, their 3-month running means, the ±0.5 °C thresholds, and warm/cool shading for detected El Nino/La Nina events.  The lower panel shows RECON_9120 − CPC on its own vertical scale, labelled in °C.
- **Figure 2**: square scatter plot of RECON_9120 (y) versus CPC 1991–2020 (x), with one point per month, identical axis ranges and ticks, a 1:1 reference line, and the correlation and RMSE annotated.
