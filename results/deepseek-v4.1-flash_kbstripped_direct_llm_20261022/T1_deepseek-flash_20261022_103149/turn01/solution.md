Using the monthly Niño 3.4 SST anomaly over 5°S–5°N, 170°W–120°W, computed relative to the 1950–2014 monthly climatology with population standard deviations (ddof=0):

### 1. Overall Niño 3.4 amplitude and model bias

| Source | Std. dev. (°C, 1950–2014) | Bias vs ERSSTv5 (°C) | Ratio vs ERSSTv5 |
|---|---:|---:|---:|
| ERSSTv5 | 0.85 | — | — |
| CESM2 | 1.29 | +0.44 | 1.52 |
| GFDL-ESM4 | 1.00 | +0.15 | 1.18 |
| MPI-ESM1-2-LR | 0.85 | +0.01 | 1.01 |

Ranking from smallest to largest absolute amplitude bias:

1. **MPI-ESM1-2-LR** (+0.01°C)
2. **GFDL-ESM4** (+0.15°C)
3. **CESM2** (+0.44°C)

MPI-ESM1-2-LR is essentially unbiased in total Niño 3.4 amplitude, GFDL-ESM4 is moderately too strong, and CESM2 strongly overestimates Niño 3.4 amplitude.

### 2. Seasonal cycle of Niño 3.4 amplitude by calendar month

Standard deviations in °C, by month:

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

### 3. Seasonal phase locking

| Source | Peak month | Peak std. dev. (°C) | Minimum month | Minimum std. dev. (°C) | Peak/min ratio |
|---|---:|---:|---:|---:|---:|
| ERSSTv5 | Dec | 1.12 | Apr | 0.60 | 1.88 |
| CESM2 | Jan | 1.62 | Jun | 0.94 | 1.72 |
| GFDL-ESM4 | Dec | 1.25 | Jul | 0.75 | 1.67 |
| MPI-ESM1-2-LR | Dec | 1.03 | Jul | 0.71 | 1.45 |

Phase-locking assessment relative to ERSSTv5:

- **CESM2** does **not** reproduce the observed peak month: it peaks in **January**, not December. Its peak-to-minimum ratio is slightly weaker than observed, 1.72 vs 1.88.
- **GFDL-ESM4** reproduces the observed **December** peak month, but its seasonal amplitude ratio is weaker than observed, 1.67 vs 1.88.
- **MPI-ESM1-2-LR** also reproduces the observed **December** peak month, but its seasonal phase locking is substantially weaker, with a peak/min ratio of 1.45 vs 1.88.

Thus, GFDL-ESM4 and MPI-ESM1-2-LR capture the December peak month but underestimate the observed strength of seasonal phase locking. CESM2 misses the peak month and has too large a background Niño 3.4 variance.

### 4. Frequency of strong Niño 3.4 anomalies

Number and percentage of months with |Niño 3.4 anomaly| strictly > 1.0°C:

| Source | Count | Percentage |
|---|---:|---:|
| ERSSTv5 | 170 | 21.8% |
| CESM2 | 360 | 46.2% |
| GFDL-ESM4 | 237 | 30.4% |
| MPI-ESM1-2-LR | 199 | 25.5% |

All three models produce too many months exceeding the ±1.0°C threshold relative to observations. This is most severe in CESM2, where the exceedance frequency is more than double the observed value, consistent with its large positive amplitude bias. MPI-ESM1-2-LR is closest to the observed exceedance frequency.

### Overall assessment

- **MPI-ESM1-2-LR** has the best total Niño 3.4 amplitude, with an almost perfect standard deviation ratio (1.01), but it underestimates seasonal phase locking.
- **GFDL-ESM4** has a modest positive amplitude bias and reproduces the December peak, but its seasonal modulation is weaker than observed.
- **CESM2** strongly overestimates Niño 3.4 amplitude and produces far too many large anomalies. Its seasonal peak occurs in January rather than December.

Figure 1 shows the four seasonal-amplitude curves with each source’s peak month marked. Figure 2 shows the stacked 1950–2014 Niño 3.4 anomaly time series on a common y-axis with ±1.0°C reference lines and each panel annotated with its standard deviation and exceedance percentage.

Caveat: these standard deviations are computed from monthly-climatology anomalies without detrending. Therefore, long-term model drift can inflate the apparent amplitude and exceedance statistics, especially for CESM2.
