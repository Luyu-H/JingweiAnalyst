## Direct answer

The three CMIP6 historical runs bracket the observed ENSO amplitude. Over 1950-01 to 2014-12, the ERSSTv5 Nino 3.4 anomaly has a population standard deviation of **0.85 degC**. **MPI-ESM1-2-LR** is closest to the observed amplitude (**0.85 degC**, bias **+0.01 degC**, ratio **1.01**), followed by **GFDL-ESM4** (**1.00 degC**, bias **+0.15 degC**, ratio **1.18**), and **CESM2** has the largest positive amplitude bias (**1.29 degC**, bias **+0.44 degC**, ratio **1.52**). The ranking from smallest to largest absolute amplitude bias is therefore:

**MPI-ESM1-2-LR < GFDL-ESM4 < CESM2**.

All four sources show December or January seasonal amplitude maxima, but no model fully matches both the observed peak month and the observed strength of seasonal phase locking as measured by the maximum-to-minimum monthly standard-deviation ratio.

## Methodology used by the executed code

The code computed the monthly Nino 3.4 index for ERSSTv5 and for CMIP6 historical `tos` from CESM2, GFDL-ESM4, and MPI-ESM1-2-LR. The period was **1950-01 to 2014-12**, giving **780 months**. The region was **5S-5N, 170W-120W**, implemented as latitude between **-5 and +5 degrees** and longitude modulo 360 between **190 and 240 degrees east**. No regridding was performed. No lagged analysis or lag range was computed; all statistics are contemporaneous over the full period.

For each CMIP6 model, the regional average used the model’s published `areacello` cell areas. For CESM2, the code additionally masked the `areacello` weights using `sftof > 0` and non-null `sftof`. For ERSSTv5, the weights were cosine-latitude weights on the native 2-degree grid, multiplied by ones over longitude. For every source, the weights were restricted to cells where the SST data were valid; the numerator was the sum of SST times weight and the denominator was the sum of valid weights. The CMIP6 output used `tos` directly, not `thetao`.

Each source’s Nino 3.4 anomaly was formed by subtracting that source’s own **1950-2014 monthly climatology**, matched by calendar month. No external climatology such as `woa23` was used. Anomalies are in degrees Celsius. Full-period and calendar-month standard deviations are population standard deviations with `ddof=0`. Exceedance counts use the strict criterion **|anomaly| > 1.0 degC** over all 780 months.

Caveats from the code and setup: each model has one historical member (`r1i1p1f1`), so model-observational differences cannot be split into forced response and internal variability. GFDL-ESM4 is on the archive `gr` 1x1-degree grid, while CESM2 and MPI-ESM1-2-LR are on native grids whose horizontal resolutions differ; no regridding was applied. `tos` is the bulk temperature of the top model cell, not a skin temperature.

## Key numerical results

### Full-period Nino 3.4 amplitude and exceedance

| Source | Full-period std (degC) | Peak month | Peak monthly std (degC) | Min month | Min monthly std (degC) | Max/min ratio | Months with \|anomaly\| > 1.0 degC | Percentage |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| ERSSTv5 | 0.85 | Dec | 1.12 | Apr | 0.60 | 1.88 | 170 | 21.8% |
| CESM2 | 1.29 | Jan | 1.62 | Jun | 0.94 | 1.72 | 360 | 46.2% |
| GFDL-ESM4 | 1.00 | Dec | 1.25 | Jul | 0.75 | 1.67 | 237 | 30.4% |
| MPI-ESM1-2-LR | 0.85 | Dec | 1.03 | Jul | 0.71 | 1.45 | 199 | 25.5% |

### Amplitude bias and ratio against ERSSTv5

Using the unrounded standard deviations from the computation:

| Model | Model std (degC) | Bias vs ERSSTv5 (degC) | Ratio vs ERSSTv5 | Absolute bias (degC) | Rank by absolute bias |
|---|---:|---:|---:|---:|---:|
| MPI-ESM1-2-LR | 0.85 | +0.01 | 1.01 | 0.005 | 1 |
| GFDL-ESM4 | 1.00 | +0.15 | 1.18 | 0.150 | 2 |
| CESM2 | 1.29 | +0.44 | 1.52 | 0.438 | 3 |

The ranking from smallest to largest absolute amplitude bias is **MPI-ESM1-2-LR < GFDL-ESM4 < CESM2**.

### Monthly standard deviations of the Nino 3.4 anomaly

Values are degC, rounded to two decimals.

| Source | Jan | Feb | Mar | Apr | May | Jun | Jul | Aug | Sep | Oct | Nov | Dec |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ERSSTv5 | 1.08 | 0.88 | 0.67 | 0.60 | 0.61 | 0.61 | 0.68 | 0.76 | 0.85 | 0.99 | 1.08 | 1.12 |
| CESM2 | 1.62 | 1.43 | 1.20 | 1.07 | 1.01 | 0.94 | 0.98 | 1.09 | 1.24 | 1.40 | 1.56 | 1.62 |
| GFDL-ESM4 | 1.13 | 1.06 | 0.96 | 0.91 | 0.84 | 0.78 | 0.75 | 0.85 | 0.98 | 1.11 | 1.21 | 1.25 |
| MPI-ESM1-2-LR | 0.96 | 0.91 | 0.78 | 0.82 | 0.86 | 0.78 | 0.71 | 0.74 | 0.77 | 0.86 | 0.97 | 1.03 |

### Seasonal phase-locking comparison

ERSSTv5 peaks in **December**, with monthly standard deviation **1.12 degC** and a maximum-to-minimum ratio of **1.88** (minimum in April, **0.60 degC**).

- **CESM2**: peak month is **January**, not December. Its peak monthly standard deviation is **1.62 degC** and its ratio is **1.72**. It does not reproduce the observed peak month; its peak amplitude and full-period amplitude are larger than observed, but its phase-locking ratio is smaller.
- **GFDL-ESM4**: peak month is **December**, matching the observed peak month. Its peak monthly standard deviation is **1.25 degC** and its ratio is **1.67**. It matches the peak month but has a smaller maximum-to-minimum ratio than observed.
- **MPI-ESM1-2-LR**: peak month is **December**, matching the observed peak month. Its peak monthly standard deviation is **1.03 degC** and its ratio is **1.45**. It matches the peak month but shows weaker seasonal amplitude contrast than observed.

Thus, GFDL-ESM4 and MPI-ESM1-2-LR reproduce the observed peak calendar month, while CESM2 peaks one month later. All three models have a lower largest-to-smallest monthly standard-deviation ratio than ERSSTv5, indicating weaker seasonal phase locking by this metric.

### Baseline monthly climatologies

The baseline is each source’s own 1950-2014 monthly climatology, matched by calendar month. Values are degC, rounded to two decimals.

| Source | Jan | Feb | Mar | Apr | May | Jun | Jul | Aug | Sep | Oct | Nov | Dec |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ERSSTv5 | 26.38 | 26.59 | 27.09 | 27.54 | 27.65 | 27.46 | 27.03 | 26.64 | 26.53 | 26.49 | 26.48 | 26.43 |
| CESM2 | 27.25 | 27.47 | 27.76 | 28.17 | 28.45 | 28.27 | 27.79 | 27.07 | 26.56 | 26.55 | 26.83 | 27.09 |
| GFDL-ESM4 | 26.09 | 26.12 | 26.23 | 26.81 | 27.44 | 27.31 | 26.84 | 26.00 | 25.45 | 25.45 | 25.61 | 25.79 |
| MPI-ESM1-2-LR | 25.38 | 25.65 | 26.17 | 26.80 | 27.07 | 26.75 | 25.96 | 25.20 | 25.13 | 25.27 | 25.38 | 25.38 |

## Visualization descriptions

**Figure 1** is a line plot of the seasonal cycle of Amplitude: calendar month January-December on the x-axis and the monthly Nino 3.4 anomaly standard deviation in degC on the y-axis. Four curves are shown on one axis, with ERSSTv5 visually distinct from the three models. Each curve’s peak month is marked. The intended logic is to compare the monthly amplitude cycle and the timing of maximum amplitude across observations and models.

**Figure 2** consists of four stacked panels in the order ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR. The panels share one time axis over 1950-2014 and an identical y-axis range that contains all four anomaly series. Each panel shows that source’s monthly Nino 3.4 anomaly series, with +1.0 and -1.0 degC reference levels marked and the panel annotated with the full-period standard deviation and exceedance count/percentage. The intended logic is to display the anomaly time series on a common scale while separating the sources into panels, allowing amplitude and exceedance behavior to be compared without implying event-by-event correspondence.

## Interpretation and insights

- **MPI-ESM1-2-LR** has near-observed full-period ENSO amplitude in Nino 3.4 and the smallest absolute amplitude bias. Its exceedance percentage, 25.5%, is also closest to the observed 21.8%, though it is still somewhat higher.
- **GFDL-ESM4** has a modest positive amplitude bias, with a standard deviation of 1.00 degC versus 0.85 degC observed. Its exceedance percentage is higher than observed, and it matches the observed December peak month but has weaker maximum-to-minimum seasonal contrast.
- **CESM2** clearly overestimates Nino 3.4 amplitude in this single historical member: its full-period standard deviation is 1.29 degC, its amplitude ratio is 1.52, and 46.2% of months exceed |1.0| degC. Its seasonal peak occurs in January rather than December.
- All models have positive full-period amplitude bias in this comparison. However, because each model has only one member, these differences cannot be cleanly attributed to model formulation versus internal variability.
- The seasonal phase-locking ratio is lower in all models than in ERSSTv5 by this metric, suggesting weaker seasonal amplitude contrast even where the peak month is correctly reproduced.
- Comparisons are statistical only and should not be interpreted as matching individual ENSO events or years. Different native/archive grids and the use of bulk top-cell `tos` versus an observational SST product also affect detailed amplitude comparisons.
