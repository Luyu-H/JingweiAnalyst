# Reference answer - T5: sea surface salinity bias, trend and scenario change

SSS is the 0-10 m volume-weighted mean salinity throughout, on each source's own grid.
Models are CESM2, GFDL-ESM4 and MPI-ESM1-2-LR; observations are WOA23 1971-2000.

## Step 2. Mean state, 1971-2000

| Source | global mean SSS (psu) | bias vs WOA23 (psu) |
|---|---|---|
| WOA23 (observed) | 34.63 | - |
| CESM2 | 34.32 | -0.32 |
| GFDL-ESM4 | 34.56 | -0.07 |
| MPI-ESM1-2-LR | 34.36 | -0.27 |

## Step 3. Pattern agreement on the shared 1-degree mask (40829 cells)

| Model | RMSE (psu) | spatial correlation |
|---|---|---|
| CESM2 | 0.77 | 0.903 |
| GFDL-ESM4 | 0.60 | 0.926 |
| MPI-ESM1-2-LR | 0.76 | 0.893 |

## Step 4. Global mean SSS trend

| Model | Period | trend (psu/decade) | 95% CI | R2 | CI excludes zero |
|---|---|---|---|---|---|
| CESM2 | historical 1950-2014 | -0.00268 | [-0.00693, 0.00157] | 0.245 | NO |
| CESM2 | ssp245 2071-2100 | -0.0209 | [-0.0255, -0.0164] | 0.910 | yes |
| GFDL-ESM4 | historical 1950-2014 | -0.0011 | [-0.00361, 0.00141] | 0.064 | NO |
| GFDL-ESM4 | ssp245 2071-2100 | -0.00534 | [-0.0106, -3.86e-05] | 0.324 | yes |
| MPI-ESM1-2-LR | historical 1950-2014 | -0.000551 | [-0.00504, 0.00393] | 0.007 | NO |
| MPI-ESM1-2-LR | ssp245 2071-2100 | -0.0231 | [-0.0513, 0.00504] | 0.701 | NO |

| Model | ssp245 / historical slope ratio | interpretable? |
|---|---|---|
| CESM2 | 7.8 | NO - denominator not distinguishable from zero |
| GFDL-ESM4 | 4.9 | NO - denominator not distinguishable from zero |
| MPI-ESM1-2-LR | 42.0 | NO - denominator not distinguishable from zero |

**Every historical trend has a confidence interval that includes zero**, so none of the three
is distinguishable from no trend at all. The ssp245 trends are an order of magnitude larger
and two of the three exclude zero. The ratio is therefore reported but is NOT a meaningful
multiple for any model: its denominator is a number consistent with zero, and the ratio's
magnitude is set by how close to zero that denominator happens to fall.

## Step 5. Scenario change, ssp245 2071-2100 minus historical 1995-2014

Marginal seas are folded into their parent ocean. The Southern and Arctic Oceans are separate
basins in the WOA23 mask and are not part of the Atlantic, Indian or Pacific numbers.

Ensemble members: CESM2 historical r1i1p1f1, ssp245 r10i1p1f1; GFDL-ESM4 historical r1i1p1f1, ssp245 r1i1p1f1; MPI-ESM1-2-LR historical r1i1p1f1, ssp245 r1i1p1f1. CESM2 uses different members in the two experiments (the ssp245 run continues a historical member that is not in the archive), so its delta SSS contains the difference in internal variability between two members as well as the forced change.

| Model | Global | Atlantic | Indian | Pacific |
|---|---|---|---|---|
| CESM2 | -0.12 | +0.11 | +0.06 | -0.24 |
| GFDL-ESM4 | -0.01 | +0.19 | -0.00 | -0.04 |
| MPI-ESM1-2-LR | -0.09 | -0.04 | -0.05 | -0.11 |

## Figures

- `fig1_sss_bias_maps.png` - model minus WOA23 SSS bias, three panels on one shared scale.
- `fig2_sss_trend_maps.png` - historical 1950-2014 SSS trend in psu/decade, p < 0.05 stippled.
- `fig3_global_mean_timeseries.png` - global mean SSS, the two experiments drawn as separate
  segments with the 2015-2070 gap left empty, each with its own fitted line.
