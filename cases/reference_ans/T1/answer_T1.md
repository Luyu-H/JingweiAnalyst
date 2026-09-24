# Reference answer - T1: Multi-model bias of the global SST climatology

Period 1971-01 to 2000-12. Reference: WOA23 temperature, 1971-2000 climate normal, annual, objectively analysed, 0 m.
Common comparison mask: 40830 cells on the WOA23 1-degree grid.

## Step 2-3. Global mean SST, each product on its own native ocean grid

| Product | global mean SST (degC) |
|---|---|
| CESM2 | 18.52 |
| GFDL-ESM4 | 18.13 |
| MPI-ESM1-2-LR | 17.58 |
| WOA23 | 18.19 |

These four numbers each cover that product's own ocean mask, so their differences are NOT
delta_SST. delta_SST is defined on the common mask below.

## Step 5-6. On the common comparison mask

The mean SST column is given to 3 decimals so that delta_SST can be checked against it. The
graded value of delta_SST is the 2-decimal one, and of RMSE and the correlation the 3-decimal
one; rounding the means to 2 decimals first would not reproduce them.

| Model | mean SST (degC) | delta_SST (degC) | RMSE (degC) | pattern correlation |
|---|---|---|---|---|
| CESM2 | 18.485 | +0.28 | 1.062 | 0.995 |
| GFDL-ESM4 | 18.118 | -0.09 | 0.832 | 0.997 |
| MPI-ESM1-2-LR | 17.573 | -0.63 | 1.339 | 0.993 |
| WOA23 (reference) | 18.204 | - | - | - |

Ranking by |delta_SST|, smallest first: **GFDL-ESM4 < CESM2 < MPI-ESM1-2-LR**.

delta_SST at full precision: CESM2 +0.2802, GFDL-ESM4 -0.0863, MPI-ESM1-2-LR -0.6308. GFDL-ESM4's value sits close to the -0.085 rounding boundary, so grade it against the
full-precision number in `answer_T1.json` rather than the rounded string.

## Figures

- `fig1_bias_maps.png` - 2x2 global maps: three bias fields on one shared diverging scale (+/-3.5 degC) plus the WOA23 reference field.
- `fig2_global_mean_bars.png` - bar chart of the four native-grid global means, each annotated with its value.
