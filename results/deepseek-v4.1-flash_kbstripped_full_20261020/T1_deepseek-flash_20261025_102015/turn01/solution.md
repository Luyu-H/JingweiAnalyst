The computation succeeded. ENSO amplitude was measured as the population standard deviation (ddof = 0) of the monthly Niño 3.4 SST anomaly over 1950-01 to 2014-12 (780 months). The Niño 3.4 region used in the data loading was 5°S–5°N, 170°W–120°W. The monthly index was produced with `compute_enso_nino_index` using `nino_region="nino3.4"`, `anomaly_freq="month"`, `climatology="self"`, and `base_period=["1950-01", "2014-12"]`, so each source used its own 1950–2014 monthly climatology. The code filtered invalid/fill SST values with `abs(SST) < 1e10`. Overall and calendar-month standard deviations were computed with `compute_basic_stat(op="std", dim="time")`. Exceedances were counted directly from the monthly index as months with `abs(Niño 3.4 anomaly) > 1.0 °C`, then expressed as a percentage of 780 months. No lagged analysis was performed.

Direct answer: CESM2 strongly overestimates ENSO amplitude, GFDL-ESM4 moderately overestimates it, and MPI-ESM1-2-LR is closest to ERSSTv5. Ranked from smallest to largest absolute amplitude bias: MPI-ESM1-2-LR, GFDL-ESM4, CESM2.

Overall Niño 3.4 anomaly standard deviation:

| Source | Std. dev. (°C) |
|---|---:|
| ERSSTv5 | 0.85 |
| CESM2 | 1.29 |
| GFDL-ESM4 | 1.00 |
| MPI-ESM1-2-LR | 0.85 |

Amplitude bias and ratio relative to ERSSTv5:

| Model | Bias (°C) | Ratio | Absolute-bias rank |
|---|---:|---:|---:|
| CESM2 | +0.44 | 1.52 | 3 |
| GFDL-ESM4 | +0.15 | 1.18 | 2 |
| MPI-ESM1-2-LR | +0.01 | 1.01 | 1 |

Calendar-month standard deviations of the Niño 3.4 anomaly (°C):

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

Seasonal amplitude peak/minimum and ratio:

| Source | Peak month | Peak std. dev. (°C) | Minimum month | Minimum std. dev. (°C) | Max/min ratio |
|---|---|---:|---|---:|---:|
| ERSSTv5 | Dec | 1.12 | Apr | 0.60 | 1.88 |
| CESM2 | Jan | 1.62 | Jun | 0.94 | 1.72 |
| GFDL-ESM4 | Dec | 1.25 | Jul | 0.75 | 1.67 |
| MPI-ESM1-2-LR | Dec | 1.03 | Jul | 0.71 | 1.45 |

Phase-locking comparison with ERSSTv5: ERSSTv5 peaks in December. GFDL-ESM4 and MPI-ESM1-2-LR also peak in December, so they reproduce the observed peak month. CESM2 peaks in January and does not reproduce the observed peak month. Peak strength relative to ERSSTv5 is 1.45 for CESM2, 1.11 for GFDL-ESM4, and 0.92 for MPI-ESM1-2-LR. The max/min ratio is smaller than observed in all three models, indicating weaker seasonal amplitude contrast than in ERSSTv5.

Strict exceedances of |Niño 3.4 anomaly| > 1.0 °C:

| Source | Exceedance months | Percentage of 780 months |
|---|---:|---:|
| ERSSTv5 | 170 | 21.8% |
| CESM2 | 360 | 46.2% |
| GFDL-ESM4 | 237 | 30.4% |
| MPI-ESM1-2-LR | 199 | 25.5% |

Figure 1 shows the seasonal cycle of ENSO amplitude: calendar months January–December on the x-axis, standard deviation of the Niño 3.4 anomaly in °C on the y-axis, and four curves for ERSSTv5, CESM2, GFDL-ESM4, and MPI-ESM1-2-LR on one axis. ERSSTv5 is visually distinguished, and each curve’s peak month is marked.

Figure 2 shows four vertically stacked panels sharing one 1950–2014 time axis, ordered ERSSTv5, CESM2, GFDL-ESM4, and MPI-ESM1-2-LR. Each panel plots that source’s monthly Niño 3.4 anomaly on an identical y-axis range spanning all four series. The +1.0 °C and −1.0 °C levels are marked with dashed lines, and each panel is annotated with its standard deviation and exceedance percentage.

Key patterns: CESM2 has the largest amplitude bias, the largest exceedance frequency, and a January peak rather than the observed December peak. GFDL-ESM4 is closer in overall amplitude and reproduces the December peak, but has more exceedances than observed. MPI-ESM1-2-LR has almost the same overall amplitude as ERSSTv5 and reproduces the December peak, but still shows a slightly higher exceedance percentage and weaker seasonal amplitude contrast. The comparison is statistical; because the historical runs are not year-synchronized with the observed record, no event-by-event matching was performed.
