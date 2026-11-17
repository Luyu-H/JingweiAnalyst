# Reference answer - T6: nitracline depth, multi-model assessment

Nitracline = depth of the maximum vertical nitrate gradient, searched within the top 1000 m; profiles are loaded to 1500 m so the deepest layer searched
has the level below it and its gradient is not one-sided.
The climatological profile is formed first and the nitracline taken from it, because WOA23
nitrate exists only as one all-period annual climatology at these depths.
Models cover 1971-2000; WOA23 nitrate is a 1965-2022 composite - the periods genuinely differ.

## Step 2. Global mean nitracline depth

| Source | global mean (m) | bias vs WOA23 (m) |
|---|---|---|
| WOA23 (observed) | 166 | - |
| CESM2 | 147 | -18 |
| GFDL-ESM4 | 134 | -32 |
| MPI-ESM1-2-LR | 132 | -34 |

## Step 3. Basin means

| Source | Pacific (m) | Atlantic (m) | Indian (m) |
|---|---|---|---|
| WOA23 | 201 | 165 | 180 |
| CESM2 | 163 | 110 | 166 |
| GFDL-ESM4 | 160 | 91 | 170 |
| MPI-ESM1-2-LR | 184 | 112 | 121 |

| Model | Pacific bias | Atlantic bias | Indian bias | ordering matches WOA23 |
|---|---|---|---|---|
| CESM2 | -38 | -56 | -14 | NO |
| GFDL-ESM4 | -41 | -74 | -10 | NO |
| MPI-ESM1-2-LR | -17 | -54 | -58 | yes |

Observed ordering, shallowest first: **Atlantic < Indian < Pacific**.
Marginal seas are folded into their parent ocean. The Southern and Arctic Oceans are separate basins in the WOA23 mask and are not part of these numbers.

## Step 4. Pattern agreement on the shared 1-degree mask (40623 cells)

| Model | spatial correlation | RMSE (m) |
|---|---|---|
| CESM2 | 0.404 | 162.9 |
| GFDL-ESM4 | 0.499 | 161.3 |
| MPI-ESM1-2-LR | 0.598 | 151.7 |

Closest to observations by RMSE: **MPI-ESM1-2-LR**.

## Step 5. What the comparison shows

- All three models place the nitracline too shallow globally, by 18 to 34 m. CESM2 and GFDL-ESM4 have their largest basin shortfalls in the Atlantic (56 and 74 m); MPI-ESM1-2-LR's largest shortfall is in the Indian Ocean (58 m, versus 54 m in the Atlantic).
- Only MPI-ESM1-2-LR reproduces the observed ordering of the three basins (Atlantic < Indian < Pacific); CESM2 and GFDL-ESM4 both put the Indian Ocean deepest.
- Spatial correlations are low (0.40 to 0.60) and the RMSE is 152 to 163 m, far worse than a temperature-based diagnostic would give. Part of that is the diagnostic itself: where the profile is smooth the depth of the maximum gradient jumps between levels, and WOA23 has 47 levels in the top 1000 m against the models' 19 to 40, so the observed map is the patchiest of the four.
- The equatorial Pacific section exposes a large CESM2 bias that the depth diagnostic alone does not show: below about 300 m in the eastern basin CESM2's nitrate falls to roughly a fifth of the observed value (8.5 against 40.2 micromoles per kilogram at 787 m, 110W), the signature of excessive water-column denitrification in an over-extended eastern tropical Pacific oxygen minimum zone. GFDL-ESM4 and MPI-ESM1-2-LR do not show it.

## Figures

- `fig1_nitracline_maps.png` - 2x2 maps of nitracline depth, one shared scale, models
  regridded to 1 degree for display only.
- `fig2_equatorial_pacific_sections.png` - 2S-2N nitrate sections from 120E to 80W in the
  top 1000 m, one shared scale in micromoles per kilogram, each source's nitracline overlaid.
