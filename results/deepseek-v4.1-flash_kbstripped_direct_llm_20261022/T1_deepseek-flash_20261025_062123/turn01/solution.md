I assessed ENSO amplitude using the monthly Niño 3.4 index over 5°S–5°N, 170°W–120°W, for 1950-01 to 2014-12. Each source was first averaged over its native grid with appropriate area weights, then anomalies were computed relative to each source’s own 1950–2014 monthly climatology. Standard deviations are population standard deviations (ddof=0).

### Overall Niño 3.4 amplitude

| Source | Niño 3.4 anomaly SD (°C) | Amplitude bias vs ERSSTv5 (°C) | Amplitude ratio vs ERSSTv5 |
|---|---:|---:|---:|
| ERSSTv5 | 0.85 | — | 1.00 |
| CESM2 | 1.29 | +0.44 | 1.52 |
| GFDL-ESM4 | 1.00 | +0.15 | 1.18 |
| MPI-ESM1-2-LR | 0.85 | +0.01 | 1.01 |

MPI-ESM1-2-LR is essentially unbiased in overall amplitude: its exact SD is 0.8541°C, bias is +0.0051°C, and ratio is 1.006, which round to 0.85°C, +0.01°C, and 1.01. CESM2 strongly overestimates ENSO amplitude; GFDL-ESM4 also overestimates, but less so.

**Ranking by smallest to largest absolute amplitude bias:**
1. MPI-ESM1-2-LR (+0.01°C)
2. GFDL-ESM4 (+0.15°C)
3. CESM2 (+0.44°C)

### Seasonal cycle of Niño 3.4 amplitude
Standard deviation of Niño 3.4 anomalies for each calendar month, in °C:

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

### Peak and minimum seasonal amplitude

| Source | Peak month | Peak SD (°C) | Minimum month | Minimum SD (°C) | Peak/min ratio |
|---|---|---:|---|---:|---:|
| ERSSTv5 | Dec | 1.12 | Apr | 0.60 | 1.88 |
| CESM2 | Jan | 1.62 | Jun | 0.94 | 1.72 |
| GFDL-ESM4 | Dec | 1.25 | Jul | 0.75 | 1.67 |
| MPI-ESM1-2-LR | Dec | 1.03 | Jul | 0.71 | 1.45 |

Observed ENSO phase locking peaks in December with a peak/minimum ratio of 1.88. CESM2 does **not** reproduce the observed peak month: it peaks in January. Its peak amplitude is about 1.45× the observed peak, while its peak/minimum ratio is about 0.92× the observed ratio, so its seasonal contrast is slightly weaker than observed despite overestimating overall amplitude.

GFDL-ESM4 reproduces the observed December peak. Its peak amplitude is about 1.11× observed, and its peak/minimum ratio is about 0.89× observed, indicating a close but slightly weaker seasonal phase locking.

MPI-ESM1-2-LR also reproduces the observed December peak, but its peak amplitude is about 0.92× observed and its peak/minimum ratio is only about 0.77× observed, so it clearly underestimates the strength of seasonal phase locking despite having nearly unbiased overall amplitude.

### Frequency of strong Niño 3.4 anomalies
Number and percentage of months with |Niño 3.4 anomaly| strictly greater than 1.0°C:

| Source | Months exceeding ±1.0°C | Percentage |
|---|---:|---:|
| ERSSTv5 | 170 / 780 | 21.8% |
| CESM2 | 360 / 780 | 46.2% |
| GFDL-ESM4 | 237 / 780 | 30.4% |
| MPI-ESM1-2-LR | 199 / 780 | 25.5% |

All three models produce more months exceeding ±1.0°C than observed. CESM2 has by far the largest overestimate; MPI-ESM1-2-LR is closest to the observed exceedance frequency, but still slightly high.

### Visualization summary
Figure 1 shows the seasonal cycle of Niño 3.4 amplitude by calendar month for all four sources, with the ERSSTv5 curve visually distinct and peak months marked. Figure 2 shows four stacked time-series panels in the order ERSSTv5, CESM2, GFDL-ESM4, and MPI-ESM1-2-LR on a common y-axis range, with ±1.0°C reference lines and each panel annotated with its SD and exceedance percentage.

Overall, MPI-ESM1-2-LR best matches the observed total Niño 3.4 amplitude, followed by GFDL-ESM4, while CESM2 substantially overestimates ENSO amplitude. However, MPI-ESM1-2-LR underestimates the seasonal phase-locking strength. GFDL-ESM4 captures the observed December peak and has a relatively realistic peak amplitude. CESM2 has the strongest ENSO amplitude and the largest excess of extreme months, and it peaks in January rather than the observed December.
