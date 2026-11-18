## Direct answer

For the Nino 3.4 region (5°S–5°N, 170°W–120°W) over **1950-01 to 2014-12 (780 months)**, using each source’s **own 1950–2014 monthly climatology** and area weighting on its native grid, the observed ERSSTv5 Nino 3.4 anomaly standard deviation is **0.85 °C**. The three CMIP6 historical runs give:

- **MPI-ESM1-2-LR: 0.85 °C**, bias **+0.01 °C**, ratio **1.01** — best overall amplitude match, smallest absolute bias.
- **GFDL-ESM4: 1.00 °C**, bias **+0.15 °C**, ratio **1.18** — slightly high amplitude.
- **CESM2: 1.29 °C**, bias **+0.44 °C**, ratio **1.52** — clearly high amplitude.

Ranking by smallest absolute amplitude bias: **MPI-ESM1-2-LR, GFDL-ESM4, CESM2**.

Seasonal phase locking is strongest in the observations: ERSSTv5 peaks in **December at 1.12 °C**, with a minimum in **April at 0.60 °C** and peak/min ratio **1.88**. GFDL-ESM4 and MPI-ESM1-2-LR both peak in December, but with weaker peak/min ratios (**1.67** and **1.45**). CESM2 peaks in **January**, not December, with the largest absolute seasonal peak (**1.62 °C**) but a weaker peak/min ratio (**1.72**) than ERSSTv5.

## Methodology used by the executed code

- Sources analysed: **ERSSTv5**, **CESM2**, **GFDL-ESM4**, and **MPI-ESM1-2-LR**.
- Time window: **1950-01 through 2014-12**, exactly **780 monthly values**.
- CMIP6 input: historical `tos` files were opened with `xr.open_mfdataset`, and the single historical variant present in the paths was used. The code found one variant for each model: **r1i1p1f1** for CESM2, GFDL-ESM4, and MPI-ESM1-2-LR. If a vertical dimension was present, the first level was selected.
- Area weighting: each source was averaged with its own cell areas. For CESM2 and MPI-ESM1-2-LR the `areacello` field was additionally masked by `sftof > 0`; for GFDL-ESM4 the area source was `areacello` directly. For ERSSTv5, regular 2-degree spherical cell areas were computed from latitude spacing and longitude spacing with **R = 6371 km**.
- Region: latitude from **−5° to +5°**, longitude **190°E to 240°E** on a 0–360 longitude axis, equivalent to **170°W–120°W**. Weighted means were computed as \(\sum(\mathrm{SST}\times\mathrm{area}) / \sum(\mathrm{area})\), with weights NaN-restricted to valid SST cells in each month.
- Baseline: each source’s anomaly was defined from **its own 1950–2014 calendar-month climatology**. The climatology for each calendar month is the plain mean of that month over the 65 years. No observed climatology or WOA climatology was applied to the models.
- Standard deviations: **population standard deviations**, `ddof=0`. Seasonal standard deviations were computed across the 65 years for each calendar month.
- Exceedance: strict threshold **|Nino 3.4 anomaly| > 1.0 °C**, counted as months out of 780, with no event-detection or minimum-duration rule.
- No lagged correlation or event-by-event comparison was performed; all metrics are contemporaneous monthly statistics.

The average Nino 3.4 box area over which each index was computed differs by source:
- ERSSTv5: **6.42 × 10⁶ km²**
- CESM2: **6.35 × 10⁶ km²**
- GFDL-ESM4: **6.19 × 10⁶ km²**
- MPI-ESM1-2-LR: **6.21 × 10⁶ km²**

The area series were constant through time after validity masking in these runs.

## Overall amplitude and exceedance statistics

| Source | Variant | Nino 3.4 anomaly std (°C) | Months with \|anomaly\| > 1.0 °C | Percent |
|---|---:|---:|---:|---:|
| ERSSTv5 | — | **0.85** | **170 / 780** | **21.8%** |
| CESM2 | r1i1p1f1 | **1.29** | **360 / 780** | **46.2%** |
| GFDL-ESM4 | r1i1p1f1 | **1.00** | **237 / 780** | **30.4%** |
| MPI-ESM1-2-LR | r1i1p1f1 | **0.85** | **199 / 780** | **25.5%** |

Amplitude bias and ratio relative to ERSSTv5:

| Model | Bias, model − ERSSTv5 (°C) | Ratio, model / ERSSTv5 | Rank by smallest absolute bias |
|---|---:|---:|---:|
| MPI-ESM1-2-LR | **+0.01** | **1.01** | **1** |
| GFDL-ESM4 | **+0.15** | **1.18** | **2** |
| CESM2 | **+0.44** | **1.52** | **3** |

Thus MPI-ESM1-2-LR essentially matches the observed overall ENSO/Nino 3.4 anomaly amplitude, GFDL-ESM4 is modestly stronger than observed, and CESM2 is substantially stronger. CESM2 also has far more months exceeding the ±1.0 °C threshold. The exceedance percentages should be read as monthly frequency statistics, not as event counts.

## Calendar-month standard deviations (°C)

| Source | Jan | Feb | Mar | Apr | May | Jun | Jul | Aug | Sep | Oct | Nov | Dec |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ERSSTv5 | 1.08 | 0.88 | 0.67 | 0.60 | 0.61 | 0.61 | 0.68 | 0.76 | 0.85 | 0.99 | 1.08 | **1.12** |
| CESM2 | **1.62** | 1.43 | 1.20 | 1.07 | 1.01 | 0.94 | 0.98 | 1.09 | 1.24 | 1.40 | 1.56 | **1.62** |
| GFDL-ESM4 | 1.13 | 1.06 | 0.96 | 0.91 | 0.84 | 0.78 | **0.75** | 0.85 | 0.98 | 1.11 | 1.21 | **1.25** |
| MPI-ESM1-2-LR | 0.96 | 0.91 | 0.78 | 0.82 | 0.86 | 0.78 | **0.71** | 0.74 | 0.77 | 0.86 | 0.97 | **1.03** |

For CESM2 the maximum seasonal standard deviation is shared/effectively at January and December when rounded to 2 decimal places, but the unrounded maximum used by the code was **January**. For the other sources, the peak month is **December**.

## Peak, minimum, and phase-locking assessment

| Source | Peak month | Peak std (°C) | Min month | Min std (°C) | Peak/min ratio | Reproduces observed peak month? | Peak-strength ratio vs ERSSTv5 | Peak/min ratio relative to ERSSTv5 |
|---|---|---:|---|---:|---:|---|---:|---:|
| ERSSTv5 | **December** | **1.12** | **April** | **0.60** | **1.88** | — | — | — |
| CESM2 | **January** | **1.62** | **June** | **0.94** | **1.72** | **No** | **1.45** | **0.92** |
| GFDL-ESM4 | **December** | **1.25** | **July** | **0.75** | **1.67** | **Yes** | **1.11** | **0.89** |
| MPI-ESM1-2-LR | **December** | **1.03** | **July** | **0.71** | **1.45** | **Yes** | **0.92** | **0.77** |

Interpretation:
- **CESM2** does not reproduce the observed December peak; its strongest seasonal amplitude occurs in January. It has a much stronger absolute December/January seasonal amplitude than observed, but its peak-to-minimum seasonal contrast is slightly weaker than observed.
- **GFDL-ESM4** reproduces the December peak month. Its December amplitude is about 11% stronger than observed, while its peak/min seasonal contrast is weaker than observed.
- **MPI-ESM1-2-LR** reproduces the December peak month, but its December amplitude is slightly weaker than observed and its seasonal peak/min contrast is substantially weaker than observed.

Overall, MPI-ESM1-2-LR matches the overall ENSO amplitude most closely, GFDL-ESM4 is intermediate, and CESM2 simulates an ENSO that is too energetic. All three models show weaker seasonal phase locking, measured as the seasonal peak/min ratio, than ERSSTv5.

## Visualization description

- **Figure 1** shows the seasonal cycle of Nino 3.4 amplitude. The x-axis is calendar month from January to December; the y-axis is the standard deviation of the monthly Nino 3.4 anomaly in °C. Four curves are plotted on the same axes: ERSSTv5 in black with a thicker line, and CESM2, GFDL-ESM4, and MPI-ESM1-2-LR in distinct colors. Each curve’s peak calendar month is marked with a star and a month label. The figure is intended to diagnose the annual cycle of ENSO amplitude and the seasonal phase locking.
- **Figure 2** shows four vertically stacked panels sharing the same 1950–2014 time axis, ordered **ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR**. Each panel plots that source’s monthly Nino 3.4 anomaly with an identical y-axis range containing all series. Dashed horizontal lines mark **+1.0 °C** and **−1.0 °C**, and a thin zero line is included. Each panel is annotated with that source’s overall standard deviation and its strict exceedance percentage. The figure is intended to compare time-domain amplitude, exceedance frequency, and visual clustering of warm/cold anomaly months under each source’s own climatology.

## Important insights and caveats

The comparison is statistical rather than event-by-event. The CMIP6 historical runs are coupled simulations and were not expected to reproduce the observed year-to-year timing of individual ENSO events. Only one historical member per model was available, so internal variability cannot be separated from forced behaviour with these data. In addition, each source is on its own native grid with its own area field, and the Nino 3.4 box area differs by source; this affects small amplitude differences and means that fine-scale features are not strictly comparable across resolutions. The reported standard deviations, biases, ratios, and exceedance percentages are based on unrounded anomaly series, with rounding applied only for reporting.
