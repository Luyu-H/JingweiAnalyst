Direct answer: ERSSTv5 has a Nino 3.4 anomaly standard deviation of 0.85 °C over 1950–2014. The three CMIP6 historical runs give CESM2 = 1.29 °C, GFDL-ESM4 = 1.00 °C, and MPI-ESM1-2-LR = 0.85 °C. Thus MPI-ESM1-2-LR has the smallest total amplitude bias, GFDL-ESM4 is intermediate, and CESM2 strongly overamplifies ENSO. Seasonal phase locking is only partly reproduced: GFDL-ESM4 and MPI-ESM1-2-LR peak in December like ERSSTv5, while CESM2 peaks in January; all three models have weaker largest-to-smallest seasonal amplitude ratios than the observations.

Methodology from the executed code: ERSSTv5 `sst` and CMIP6 `tos` for CESM2, GFDL-ESM4, and MPI-ESM1-2-LR were loaded with `load_regular_data` for 1950-01 through 2014-12, with `field='total'`, latitude range [-10, 10], and longitude range [-180, -110]. Each source remained on its native grid; no horizontal regridding is present in the executed code. The Nino 3.4 index was computed with `compute_enso_nino_index` using `nino_region='nino3.4'`, `depth_level=None`, `anomaly_freq='month'`, `climatology='self'`, and `base_period=time_range`, so each source’s anomalies are relative to its own 1950–2014 monthly climatology. The code does not show a separate `load_areas` call; regional averaging is delegated to `compute_enso_nino_index`. Standard deviations were computed with `compute_basic_stat(..., op='std', dim='time')` over all 780 months and separately for each calendar month. Exceedance counts are own-code counts of months with `abs(anomaly) > 1.0 °C`, divided by 780.

Overall amplitude and exceedance:

| Source | N months | Nino 3.4 std (°C) | |anomaly| > 1.0 count | Percentage |
|---|---:|---:|---:|---:|
| ERSSTv5 | 780 | 0.85 | 170 | 21.8% |
| CESM2 | 780 | 1.29 | 360 | 46.2% |
| GFDL-ESM4 | 780 | 1.00 | 237 | 30.4% |
| MPI-ESM1-2-LR | 780 | 0.85 | 199 | 25.5% |

Model amplitude bias and ratio relative to ERSSTv5:

| Model | Bias (°C) | Ratio | Absolute-bias rank |
|---|---:|---:|---:|
| CESM2 | +0.44 | 1.52 | 3 |
| GFDL-ESM4 | +0.15 | 1.18 | 2 |
| MPI-ESM1-2-LR | +0.01 | 1.01 | 1 |

Ranking from smallest to largest absolute amplitude bias: MPI-ESM1-2-LR, GFDL-ESM4, CESM2.

Seasonal (calendar-month) Nino 3.4 anomaly standard deviations, °C:

| Source | Jan | Feb | Mar | Apr | May | Jun | Jul | Aug | Sep | Oct | Nov | Dec |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ERSSTv5 | 1.08 | 0.88 | 0.67 | 0.60 | 0.61 | 0.61 | 0.68 | 0.76 | 0.85 | 0.99 | 1.08 | 1.12 |
| CESM2 | 1.62 | 1.43 | 1.20 | 1.07 | 1.01 | 0.94 | 0.98 | 1.09 | 1.24 | 1.40 | 1.56 | 1.62 |
| GFDL-ESM4 | 1.13 | 1.06 | 0.96 | 0.91 | 0.84 | 0.78 | 0.75 | 0.85 | 0.98 | 1.11 | 1.21 | 1.25 |
| MPI-ESM1-2-LR | 0.96 | 0.91 | 0.78 | 0.82 | 0.86 | 0.78 | 0.71 | 0.74 | 0.77 | 0.86 | 0.97 | 1.03 |

Seasonal phase-locking summary and comparison with ERSSTv5:

| Source | Peak month | Peak std (°C) | Minimum month | Minimum std (°C) | Max/min ratio | Peak month matches ERSSTv5? | Strength relative to ERSSTv5 |
|---|---:|---:|---|---:|---:|---|---|
| ERSSTv5 | Dec | 1.12 | Apr | 0.60 | 1.88 | — | — |
| CESM2 | Jan | 1.62 | Jun | 0.94 | 1.72 | No | Weaker ratio (1.72 vs 1.88) |
| GFDL-ESM4 | Dec | 1.25 | Jul | 0.75 | 1.67 | Yes | Weaker ratio (1.67 vs 1.88) |
| MPI-ESM1-2-LR | Dec | 1.03 | Jul | 0.71 | 1.45 | Yes | Weaker ratio (1.45 vs 1.88) |

Thus, GFDL-ESM4 and MPI-ESM1-2-LR reproduce the observed December peak month, while CESM2 peaks one month later in January. None reproduces the observed phase-locking strength in a strict sense: all model max/min ratios are below the ERSSTv5 ratio of 1.88. GFDL-ESM4 is the closest of the two December-peaking models on ratio, whereas MPI-ESM1-2-LR is closest in total amplitude but has the weakest seasonal contrast. CESM2 has the largest overall amplitude and the highest exceedance rate, and its ratio is closer to observed than GFDL/MPI, but its peak occurs in January rather than December.

Visualizations: Figure 1 is an own-code line plot of the seasonal amplitude cycle, with January–December on the x-axis, Nino 3.4 anomaly standard deviation (°C) on the y-axis, and four curves on one axis. ERSSTv5 is plotted as a thick solid black line and the three models use different colors and line styles; each curve’s maximum seasonal standard deviation is marked. This figure is intended to compare the annual cycle of ENSO amplitude and seasonal phase locking across sources.

Figure 2 is an own-code four-panel stacked time-series plot with a shared 1950–2014 time axis. The panels are ordered ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR. Each panel shows that source’s monthly Nino 3.4 anomaly series on an identical y-axis range spanning all four series, with horizontal lines at +1.0 °C, 0 °C, and -1.0 °C. Each panel is annotated with its overall standard deviation and its exceedance count/percentage, so the figure compares event amplitude, threshold exceedance, and temporal evolution across the four sources.

Key insight: total ENSO amplitude is best captured by MPI-ESM1-2-LR, but its seasonal phase locking is too weak; GFDL-ESM4 has a moderate positive amplitude bias and correct peak month but also underestimates the seasonal amplitude ratio; CESM2 has a large positive amplitude bias, many more months exceeding ±1 °C, and a one-month shift in its peak seasonal amplitude. Because each source uses its own 1950–2014 monthly climatology, these results describe simulated variability amplitude, exceedance behavior, and seasonal phase locking rather than mean-state SST biases.
