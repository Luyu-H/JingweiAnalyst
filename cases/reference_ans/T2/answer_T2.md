# Reference answer - T2: Nino 3.4 index reconstruction validation

Window 1950-01 to 2023-12, 888 months. Index rebuilt from ERSSTv5 monthly SST, area-weighted
over 5S-5N / 170W-120W on the native 2-degree grid, anomalies against a 1991-2020 monthly
climatology from the same field.

## Step 2. Agreement with the CPC 1991-2020 series

| Quantity | Value |
|---|---|
| Pearson r | 0.999972 |
| RMSE | 0.0075 degC |
| Largest absolute deviation | 0.0226 degC, at 1992-01 |
| Mean bias (reconstruction - CPC) | -0.0030 degC |

CPC publishes the series to 2 decimal places, and the reconstruction sits well inside that
last digit, so the correlation is quoted to 6 places and the errors to 4; at the 4 and 3
places of the original draft these would read 1.0000 and 0.007 and say nothing.

## Step 3. El Nino events in the reconstruction, ONI criterion

**16 El Nino events** (and 25 La Nina events) over 1950-2023.

| Peak | Peak value (degC) | Event span | Seasons |
|---|---|---|---|
| 2015-12 | 2.64 | 2015-03 – 2016-04 | 14 |
| 1997-12 | 2.35 | 1997-05 – 1998-04 | 12 |
| 1983-01 | 2.15 | 1982-07 – 1983-05 | 11 |

## Step 4. The same criterion on the CPC ONI

**22 El Nino events** (and 20 La Nina events) over the same window.

| Peak | Peak value (degC) | Event span | Seasons |
|---|---|---|---|
| 2015-12 | 2.59 | 2014-10 – 2016-05 | 20 |
| 1997-12 | 2.37 | 1997-05 – 1998-04 | 12 |
| 1983-01 | 2.14 | 1982-05 – 1983-06 | 14 |

The counts differ (16 against 22) because the ONI re-centres its
baseline every five years and therefore carries no long-term warming trend, whereas the
reconstruction is referenced to a single fixed 1991-2020 climatology. Against that fixed
baseline the 1950s-1980s read cold, and the early El Nino events fall short of +0.5 degC.
The three strongest events are the same three in both, at the same peak months.

## Steps 5-6. Baseline offset, 1981-2010 relative to 1991-2020

Both means taken over 1950-01 to 2020-12 (852 months), where the CPC 1981-2010 series ends.

| Source | Mean offset (degC) |
|---|---|
| Reconstructed, RECON_8110 - RECON_9120 | +0.060 |
| Published, CPC 81-10 - CPC 91-20 | +0.058 |

The two agree to 0.0018 degC. The reconstruction
therefore reproduces CPC's baseline arithmetic, not merely the variability of the index.

## Figures

- `fig1_series_overlay.png` - the two monthly series overlaid with the 3-month running mean,
  warm/cool shading over each detected event, and the reconstruction-minus-CPC residual below.
- `fig2_scatter.png` - reconstruction against CPC on identical square axes with a 1:1 line.
