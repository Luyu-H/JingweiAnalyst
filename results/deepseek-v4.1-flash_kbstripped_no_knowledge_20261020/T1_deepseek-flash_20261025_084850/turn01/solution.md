## Methodology (as implemented)

The code loaded monthly SST for 1950-01 through 2014-12:

- CMIP6 historical `tos` for CESM2, GFDL-ESM4, and MPI-ESM1-2-LR from `cmip6_historical`;
- ERSSTv5 `sst` from `ersst`.

Data were loaded over a wider box of 10°S–10°N, 180°W–110°W, but the ENSO index itself was computed with `compute_enso_nino_index(..., nino_region="nino3.4", anomaly_freq="month", climatology="self", base_period=time_range)`. Thus the Nino 3.4 region was 5°S–5°N, 170°W–120°W, anomalies were relative to each source’s own 1950–2014 monthly climatology, and no horizontal regridding was applied: each source was processed on its native grid through the index function. This produced 780 monthly Nino 3.4 anomalies per source.

Full-period and calendar-month standard deviations were computed with `compute_basic_stat(..., op="std", dim="time")`. The analysis convention is population standard deviation, `ddof=0`. Calendar-month statistics were formed by grouping on `time.dt.month`, giving 65 values per calendar month. Exceedances were computed with `compute_exceedance_fraction(data=nino[s], threshold=1.0, comparison="abs_greater", dim=None, weighting="none")`, so the criterion was strictly `|Nino3.4 anomaly| > 1.0 °C` over the 780 months.

## Full-period ENSO amplitude and model biases

Full-period Nino 3.4 anomaly standard deviation:

| Source | Std. dev. (°C) |
|---|---:|
| ERSSTv5 | 0.85 |
| CESM2 | 1.29 |
| GFDL-ESM4 | 1.00 |
| MPI-ESM1-2-LR | 0.85 |

Amplitude bias and ratio relative to ERSSTv5:

| Model | Bias, model − ERSSTv5 (°C) | Ratio, model / ERSSTv5 |
|---|---:|---:|
| CESM2 | +0.44 | 1.52 |
| GFDL-ESM4 | +0.15 | 1.18 |
| MPI-ESM1-2-LR | +0.01 | 1.01 |

Ranking from smallest to largest absolute amplitude bias:

**MPI-ESM1-2-LR < GFDL-ESM4 < CESM2**

Thus CESM2 strongly overestimates Nino 3.4 amplitude, GFDL-ESM4 overestimates moderately, and MPI-ESM1-2-LR is closest to the observed amplitude.

## Calendar-month standard deviations

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

Seasonal amplitude extremes:

| Source | Peak month | Peak std. dev. (°C) | Trough month | Trough std. dev. (°C) | Peak/trough ratio |
|---|---|---:|---|---:|---:|
| ERSSTv5 | Dec | 1.12 | Apr | 0.60 | 1.88 |
| CESM2 | Jan | 1.62 | Jun | 0.94 | 1.72 |
| GFDL-ESM4 | Dec | 1.25 | Jul | 0.75 | 1.67 |
| MPI-ESM1-2-LR | Dec | 1.03 | Jul | 0.71 | 1.45 |

Observed ERSSTv5 peaks in December with a peak/trough ratio of 1.88. GFDL-ESM4 and MPI-ESM1-2-LR also peak in December, while CESM2 peaks in January. All three models have lower peak/trough ratios than the observations, indicating weaker simulated seasonal phase locking/seasonal contrast than in ERSSTv5.

## Exceedance of |Nino 3.4 anomaly| > 1.0 °C

| Source | Months exceeding | Percentage of 780 months |
|---|---:|---:|
| ERSSTv5 | 170 | 21.8% |
| CESM2 | 360 | 46.2% |
| GFDL-ESM4 | 237 | 30.4% |
| MPI-ESM1-2-LR | 199 | 25.5% |

CESM2 has by far the most frequent strong anomalies, consistent with its high amplitude. MPI-ESM1-2-LR has an exceedance rate closest to the observed value.

## Figures

**Figure 1** is a line plot of the seasonal cycle of ENSO amplitude: calendar month January–December on the x-axis and population standard deviation of the Nino 3.4 anomaly in °C on the y-axis. Four curves are drawn on one axis, with ERSSTv5 plotted distinctly from the three models. Each curve’s peak calendar month is marked and annotated. The figure is designed to compare the annual cycle of Nino 3.4 amplitude and phase locking among the observations and the three CMIP6 runs.

**Figure 2** contains four vertically stacked panels sharing one time axis from 1950 to 2014, ordered ERSSTv5, CESM2, GFDL-ESM4, and MPI-ESM1-2-LR. Each panel shows that source’s monthly Nino 3.4 anomaly series. All panels use an identical y-axis range containing every series. The +1.0 °C and −1.0 °C levels are marked, and |anomaly| > 1.0 °C periods are shaded. Each panel is annotated with the full-period standard deviation and the exceedance count/percentage, allowing direct comparison of amplitude and event frequency across sources.

## Key interpretation

The ENSO amplitude ranking relative to ERSSTv5 is clear: MPI-ESM1-2-LR is nearly unbiased in full-period standard deviation (ratio 1.01), GFDL-ESM4 is moderately too strong (1.18), and CESM2 is substantially too strong (1.52). The exceedance statistics follow the same order, with CESM2 producing strong |anomaly| > 1.0 °C events far more often than observed. Although GFDL-ESM4 and MPI-ESM1-2-LR reproduce the observed December peak month, all models show weaker seasonal phase locking than ERSSTv5, as measured by their smaller peak-to-trough ratios. CESM2 additionally shifts the peak month to January.
