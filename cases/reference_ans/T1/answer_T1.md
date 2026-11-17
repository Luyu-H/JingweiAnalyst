# Reference answer - T1: multi-model assessment of simulated ENSO amplitude

Window 1950-01 to 2014-12, 780 months. Nino 3.4 index area-weighted over 5S-5N / 170W-120W
on each source's own native grid, anomalies against that source's own 1950-2014 monthly
climatology. Standard deviations use ddof = 0.

## Step 2. Overall amplitude

| Source | std (degC) | bias vs obs (degC) | ratio to obs |
|---|---|---|---|
| ERSSTv5 (observed) | 0.85 | - | - |
| CESM2 | 1.29 | +0.44 | 1.52 |
| GFDL-ESM4 | 1.00 | +0.15 | 1.18 |
| MPI-ESM1-2-LR | 0.85 | +0.01 | 1.01 |

Ranking by |amplitude bias|, closest first: **MPI-ESM1-2-LR < GFDL-ESM4 < CESM2**.

## Step 3. Seasonal phase locking

Standard deviation within each calendar month (degC):

| Source | Jan | Feb | Mar | Apr | May | Jun | Jul | Aug | Sep | Oct | Nov | Dec |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ERSSTv5 | 1.08 | 0.88 | 0.67 | 0.60 | 0.61 | 0.61 | 0.68 | 0.76 | 0.85 | 0.99 | 1.08 | 1.12 |
| CESM2 | 1.62 | 1.43 | 1.20 | 1.07 | 1.01 | 0.94 | 0.98 | 1.09 | 1.24 | 1.40 | 1.56 | 1.62 |
| GFDL-ESM4 | 1.13 | 1.06 | 0.96 | 0.91 | 0.84 | 0.78 | 0.75 | 0.85 | 0.98 | 1.11 | 1.21 | 1.25 |
| MPI-ESM1-2-LR | 0.96 | 0.91 | 0.78 | 0.82 | 0.86 | 0.78 | 0.71 | 0.74 | 0.77 | 0.86 | 0.97 | 1.03 |

| Source | peak month | peak (degC) | trough month | trough (degC) | peak/trough |
|---|---|---|---|---|---|
| ERSSTv5 | Dec | 1.12 | Apr | 0.60 | 1.88 |
| CESM2 | Jan | 1.62 | Jun | 0.94 | 1.72 |
| GFDL-ESM4 | Dec | 1.25 | Jul | 0.75 | 1.67 |
| MPI-ESM1-2-LR | Dec | 1.03 | Jul | 0.71 | 1.45 |

## Step 4. Months with |anomaly| > 1.0 degC

| Source | months | percentage |
|---|---|---|
| ERSSTv5 | 170 / 780 | 21.8% |
| CESM2 | 360 / 780 | 46.2% |
| GFDL-ESM4 | 237 / 780 | 30.4% |
| MPI-ESM1-2-LR | 199 / 780 | 25.5% |

## Figures

- `fig1_seasonal_amplitude.png` - the 12-month amplitude cycle, four curves on one axis,
  each source's peak month ringed.
- `fig2_nino34_panels.png` - the four series stacked on one identical vertical scale (+/-4.00 degC),
  with the +/-1.0 degC levels marked and each panel annotated with its std and exceedance.
