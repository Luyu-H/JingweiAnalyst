**Direct answer.** The three CMIP6 historical runs differ substantially in simulated ENSO amplitude. Using the monthly Niño 3.4 anomaly over 1950-01 to 2014-12, ERSSTv5 has a standard deviation of **0.85 °C**. MPI-ESM1-2-LR is closest (**0.85 °C**, bias **+0.01 °C**, ratio **1.01**), GFDL-ESM4 is somewhat high (**1.00 °C**, bias **+0.15 °C**, ratio **1.18**), and CESM2 is strongly high (**1.29 °C**, bias **+0.44 °C**, ratio **1.52**). Ranking by smallest to largest absolute amplitude bias: **MPI-ESM1-2-LR, GFDL-ESM4, CESM2**.

**Methodology actually used.** The analysis loaded ERSSTv5 monthly `sst` via `dataset='ersst'` and CMIP6 historical monthly `tos` via `dataset='cmip6_historical'` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR. All sources were loaded for **1950-01 to 2014-12** on their native grids, first over latitude **-10 to 10°** and longitude **-180 to -100°**. The Niño 3.4 index was then formed with `compute_enso_nino_index` using `nino_region='nino3.4'`, `anomaly_freq='month'`, `climatology='self'`, and `base_period=['1950-01','2014-12']`, yielding exactly **780** monthly values per source. Anomalies are in **°C** relative to the 1950–2014 calendar-month climatology. No regridding was done by the executed code, and no lag analysis was performed. Standard deviations were computed with `compute_basic_stat(op='std')` over all 780 months and over each calendar-month subset of 65 values; the reported values are the requested population standard deviations. Strict exceedance was counted with `compute_exceedance_fraction`: anomalies **> +1.0 °C**, anomalies **< -1.0 °C**, and checked against an absolute-threshold count.

**All-month Niño 3.4 anomaly standard deviation (1950–2014):**

| Source | Std. dev. (°C) |
|---|---:|
| ERSSTv5 | 0.85 |
| CESM2 | 1.29 |
| GFDL-ESM4 | 1.00 |
| MPI-ESM1-2-LR | 0.85 |

**Model amplitude comparison against ERSSTv5:**

| Model | Std. dev. (°C) | Bias, model − ERSSTv5 (°C) | Ratio, model / ERSSTv5 | Rank by absolute bias |
|---|---:|---:|---:|---:|
| MPI-ESM1-2-LR | 0.85 | +0.01 | 1.01 | 1 |
| GFDL-ESM4 | 1.00 | +0.15 | 1.18 | 2 |
| CESM2 | 1.29 | +0.44 | 1.52 | 3 |

**Calendar-month standard deviations of the Niño 3.4 anomaly (°C):**

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

**Seasonal amplitude peaks and phase locking:**

| Source | Peak month | Peak std. dev. (°C) | Minimum month | Minimum std. dev. (°C) | Max/min ratio | Reproduces observed peak month? | Ratio vs ERSSTv5 max/min |
|---|---:|---:|---:|---:|---:|---|---:|
| ERSSTv5 | Dec | 1.12 | Apr | 0.60 | 1.88 | — | 1.00 |
| CESM2 | Jan | 1.62 | Jun | 0.94 | 1.72 | No | 0.92 |
| GFDL-ESM4 | Dec | 1.25 | Jul | 0.75 | 1.67 | Yes | 0.89 |
| MPI-ESM1-2-LR | Dec | 1.03 | Jul | 0.71 | 1.45 | Yes | 0.77 |

The observed seasonal phase locking peaks in **December** with a max/min ratio of **1.88**. GFDL-ESM4 and MPI-ESM1-2-LR both reproduce the December peak month but with weaker phase-locking strength, ratios **1.67** and **1.45**, respectively. CESM2 peaks in **January** rather than December and has a ratio of **1.72**, so it does not reproduce the observed peak month, though its seasonal contrast is closer to observed than the other two models.

**Strict exceedance of |Niño 3.4 anomaly| > 1.0 °C:**

| Source | Months > +1.0 °C | Months < -1.0 °C | Total | Percentage of 780 months |
|---|---:|---:|---:|---:|
| ERSSTv5 | 85 | 85 | 170 | 21.8% |
| CESM2 | 176 | 184 | 360 | 46.2% |
| GFDL-ESM4 | 124 | 113 | 237 | 30.4% |
| MPI-ESM1-2-LR | 107 | 92 | 199 | 25.5% |

All three models exceed the observed frequency of strong Niño 3.4 anomalies. CESM2 has the largest overestimate: **46.2%** of months versus **21.8%** in ERSSTv5. GFDL-ESM4 is intermediate at **30.4%**, and MPI-ESM1-2-LR is closest at **25.5%**, despite having an all-month standard deviation nearly identical to the observations. Thus amplitude in the total standard deviation and the frequency of large excursions are related but not identical diagnostics.

**Visualizations.** Figure 1 is the seasonal-cycle amplitude plot: calendar month January–December on the x-axis and the 12 calendar-month Niño 3.4 standard deviations in °C on the y-axis, with four curves on one axis. ERSSTv5 is visually distinct using a thick black line with markers, while the three models use separate line styles/colors; each curve’s peak month is marked with an open circle. Figure 2 consists of four stacked time-series panels in the order ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, sharing the 1950–2014 time axis. Each panel plots that source’s Niño 3.4 anomaly series on the same y-axis range, includes ±1.0 °C reference lines, and annotates the panel with its all-month standard deviation and strict exceedance percentage/count.
