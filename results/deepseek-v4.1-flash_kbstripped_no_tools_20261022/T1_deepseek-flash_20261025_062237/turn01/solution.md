**Methodology used.** The analysis used monthly sea surface temperature over **1950-01 through 2014-12** (780 months). For each source, the Niño 3.4 regional mean was computed over **5°S–5°N, 170°W–120°W** (190°E–240°E) from native-grid cell centres inside the box, with no interpolation or remapping. The regional mean was an area-weighted sum divided by the sum of weights over valid ocean cells only. ERSSTv5 used cosine-latitude weights on its native 2° grid. CESM2 and MPI-ESM1-2-LR used their published `areacello` on their `gn` grids; CESM2 was further masked with `sftof > 0`. GFDL-ESM4 was processed on its `grid_label='gr'` archive grid using its published `areacello` (a grid difference from the `gn` models), without remapping. All three CMIP6 runs used variant `r1i1p1f1`.

Each source’s Niño 3.4 anomaly was formed relative to **its own 1950–2014 calendar-month climatology**, i.e. the plain 65-year mean for each calendar month; no common observed climatology was applied to models. Monthly `sst`/`tos` values are monthly means, so the climatology has no day weighting. All standard deviations are population standard deviations (`ddof=0`). Exceedance counts use the strict rule `|anomaly| > 1.0°C`, with denominator 780 months and no event-persistence rule.

**Overall Niño 3.4 amplitude and model comparison.**

| Source | Overall SD (°C) | Bias vs ERSSTv5 (°C) | Ratio vs ERSSTv5 | Months with |anomaly| > 1.0°C | Percentage |
|---|---:|---:|---:|---:|---:|
| ERSSTv5 | **0.85** | — | — | **170** | **21.8%** |
| CESM2 | **1.29** | **+0.44** | **1.52** | **360** | **46.2%** |
| GFDL-ESM4 | **1.00** | **+0.15** | **1.18** | **237** | **30.4%** |
| MPI-ESM1-2-LR | **0.85** | **+0.01** | **1.01** | **199** | **25.5%** |

Ranking of models by smallest to largest absolute amplitude bias: **MPI-ESM1-2-LR (0.01°C), GFDL-ESM4 (0.15°C), CESM2 (0.44°C)** using unrounded values internally.

**Monthly Niño 3.4 anomaly standard deviations by calendar month (°C, ddof=0).**

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

**Seasonal peak, minimum, and phase-locking ratio.**

| Source | Peak month | Peak SD (°C) | Minimum month | Minimum SD (°C) | Peak/min ratio |
|---|---|---:|---|---:|---:|
| ERSSTv5 | **December** | **1.12** | **April** | **0.60** | **1.88** |
| CESM2 | **January** | **1.62** | **June** | **0.94** | **1.72** |
| GFDL-ESM4 | **December** | **1.25** | **July** | **0.75** | **1.67** |
| MPI-ESM1-2-LR | **December** | **1.03** | **July** | **0.71** | **1.45** |

ERSSTv5 peaks in **December** with a **1.88** peak-to-minimum ratio. CESM2’s nominal peak is **January** (with December nearly tied at 1.62°C), so it does not reproduce the observed December peak month, and its phase-locking ratio (1.72) is somewhat weaker despite much larger overall amplitude. GFDL-ESM4 and MPI-ESM1-2-LR both reproduce the **December** peak month, but their phase-locking ratios (1.67 and 1.45, respectively) are weaker than observed. MPI-ESM1-2-LR has the closest overall amplitude to ERSSTv5, while CESM2 is strongly overamplified.

**Visualizations.** Figure 1 shows the seasonal cycle of amplitude as a line plot with calendar month January–December on the x-axis and Niño 3.4 anomaly SD on the y-axis. ERSSTv5 is plotted as a visually distinct black solid curve; the three models are dashed colored curves. Each curve’s peak month is marked with a star and labelled. Figure 2 shows four stacked panels sharing the 1950–2014 time axis, ordered ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR. Each panel displays that source’s Niño 3.4 anomaly series on a common y-axis range containing all series, with horizontal ±1.0°C reference lines and annotated overall SD and exceedance percentage.

**Key insights.** ERSSTv5 has a Niño 3.4 monthly anomaly SD of **0.85°C** and a strong December peak with April minimum. MPI-ESM1-2-LR reproduces the overall amplitude almost exactly (**0.85°C**, +0.01°C bias), and GFDL-ESM4 is only moderately high (**1.00°C**, +0.15°C). CESM2 substantially overestimates monthly Niño 3.4 variability (**1.29°C**, +52%) and has far more months exceeding ±1.0°C (**46.2%** vs **21.8%** observed). All models show weaker seasonal phase locking than ERSSTv5, judged by the peak-to-minimum ratio, even when the December peak month is reproduced.
