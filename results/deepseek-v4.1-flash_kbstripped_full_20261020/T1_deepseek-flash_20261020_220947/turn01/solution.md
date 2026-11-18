**Direct answer.** Over 1950-01 to 2014-12, the three CMIP6 historical runs bracket the observed ERSSTv5 Niño 3.4 amplitude. MPI-ESM1-2-LR is closest to observations, GFDL-ESM4 is modestly stronger, and CESM2 is substantially stronger. Ranking of absolute amplitude bias from smallest to largest is **MPI-ESM1-2-LR < GFDL-ESM4 < CESM2**. Seasonal phase locking is only partly reproduced: GFDL-ESM4 and MPI-ESM1-2-LR peak in December like ERSSTv5, while CESM2 peaks in January; all three have weaker maximum-to-minimum seasonal amplitude ratios than ERSSTv5.

**Methodology actually used in the executed code.**
- The code loaded ERSSTv5 `sst` and CMIP6 historical `tos` for CESM2, GFDL-ESM4, and MPI-ESM1-2-LR using `load_regular_data`, over **1950-01 to 2014-12**, latitude **-5 to 5**, longitude **-170 to -120**.
- For each source, it called `compute_enso_nino_index` with `nino_region="nino3.4"`, `anomaly_freq="month"`, `climatology="self"`, and `base_period=time_range`. Thus each 780-month Niño 3.4 anomaly series uses that source’s own 1950–2014 calendar-month climatology.
- Overall and calendar-month standard deviations were computed with `compute_basic_stat(..., op="std", dim="time")`, using the 780 monthly anomaly values and the requested population-standard-deviation convention. Calendar-month statistics were obtained by manually selecting months 1–12 and applying the same standard-deviation operation.
- Exceedance counts use individual months with absolute Niño 3.4 anomaly **strictly greater than 1.0 °C**; no persistence-defined event detection or lag analysis was performed.
- Bias and ratio were computed from the unrounded overall standard deviations. The displayed code does not explicitly call `load_areas` and does not show a separate regridding step; the regional indices are those returned by `compute_enso_nino_index` from the loaded fields.

**Overall Niño 3.4 amplitude, exceedance, and seasonal-locking summary.**

| Source | Overall std (°C) | Bias vs ERSSTv5 (°C) | Ratio vs ERSSTv5 | Months with \|anomaly\|>1.0 °C | Percent | Peak month | Peak std (°C) | Min month | Min std (°C) | Max/min ratio |
|---|---:|---:|---:|---:|---:|---|---:|---|---:|---:|
| ERSSTv5 | 0.85 | — | — | 170 | 21.8% | Dec | 1.12 | Apr | 0.60 | 1.88 |
| CESM2 | 1.29 | +0.44 | 1.52 | 360 | 46.2% | Jan | 1.62 | Jun | 0.94 | 1.72 |
| GFDL-ESM4 | 1.00 | +0.15 | 1.18 | 237 | 30.4% | Dec | 1.25 | Jul | 0.75 | 1.67 |
| MPI-ESM1-2-LR | 0.85 | +0.01 | 1.01 | 199 | 25.5% | Dec | 1.03 | Jul | 0.71 | 1.45 |

The four overall standard deviations to two decimals are **ERSSTv5 0.85 °C, CESM2 1.29 °C, GFDL-ESM4 1.00 °C, MPI-ESM1-2-LR 0.85 °C**. The amplitude ranking by absolute bias is **MPI-ESM1-2-LR (smallest), then GFDL-ESM4, then CESM2 (largest)**.

**Calendar-month standard deviations (°C), January–December.**

| Month | ERSSTv5 | CESM2 | GFDL-ESM4 | MPI-ESM1-2-LR |
|---|---:|---:|---:|---:|
| Jan | 1.08 | 1.62 | 1.13 | 0.96 |
| Feb | 0.88 | 1.43 | 1.06 | 0.91 |
| Mar | 0.67 | 1.20 | 0.96 | 0.78 |
| Apr | 0.60 | 1.07 | 0.91 | 0.82 |
| May | 0.61 | 1.01 | 0.84 | 0.86 |
| Jun | 0.61 | 0.94 | 0.78 | 0.78 |
| Jul | 0.68 | 0.98 | 0.75 | 0.71 |
| Aug | 0.76 | 1.09 | 0.85 | 0.74 |
| Sep | 0.85 | 1.24 | 0.98 | 0.77 |
| Oct | 0.99 | 1.40 | 1.11 | 0.86 |
| Nov | 1.08 | 1.56 | 1.21 | 0.97 |
| Dec | 1.12 | 1.62 | 1.25 | 1.03 |

**Seasonal phase-locking interpretation.**
- ERSSTv5 peaks in **December (1.12 °C)** and has its minimum in **April (0.60 °C)**, with a max/min ratio of **1.88**.
- CESM2 peaks in **January (1.62 °C)** with a ratio of **1.72**; it does **not** reproduce the observed peak month and has weaker seasonal phase locking than ERSSTv5.
- GFDL-ESM4 peaks in **December (1.25 °C)** with a ratio of **1.67**; it reproduces the observed peak month but underestimates the observed seasonal amplitude ratio.
- MPI-ESM1-2-LR peaks in **December (1.03 °C)** with a ratio of **1.45**; it also reproduces the peak month but has the weakest seasonal phase locking of the four sources.

**Exceedance and amplitude interpretation.** All three models exceed the observed frequency of months with |Niño 3.4 anomaly| > 1.0 °C. CESM2 stands out with 360 months (46.2%), followed by GFDL-ESM4 with 237 months (30.4%), and MPI-ESM1-2-LR with 199 months (25.5%), compared with ERSSTv5’s 170 months (21.8%). This is consistent with the overall amplitude ranking: CESM2 has the largest positive amplitude bias, GFDL-ESM4 is intermediate, and MPI-ESM1-2-LR is nearly unbiased in overall standard deviation.

**Visualizations produced.**
- **Figure 1** is a line plot of the seasonal cycle of Niño 3.4 amplitude: calendar month January–December on the x-axis and standard deviation (°C) on the y-axis. Four curves are shown on one axis, with ERSSTv5 visually distinct in black and the three models in different colors/line styles. Each curve’s peak month is marked and identified in the legend.
- **Figure 2** consists of four vertically stacked panels in the order ERSSTv5, CESM2, GFDL-ESM4, and MPI-ESM1-2-LR, sharing one 1950–2014 time axis. Each panel shows that source’s Niño 3.4 anomaly series on an identical y-axis range containing all series, with 0 °C and ±1.0 °C levels marked. Each panel is annotated with the overall standard deviation and the percentage of months exceeding |1.0| °C.

**Important caveats.** The executed comparison uses one historical member per model, so model-to-model amplitude differences can include internal variability and cannot be attributed solely to forced differences. Also, the displayed code does not explicitly load published cell areas or regrid the fields; the reported indices are those produced directly by the `load_regular_data` and `compute_enso_nino_index` calls shown. Any additional native-grid or top-cell-definition effects are therefore not separately quantified in the executed code.
