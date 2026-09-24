# Reference answer - T7: intermediate-water N:P against Redfield

Depth range 200-1000 m. Models 1971-2000; WOA23 nitrate and phosphate are the
1965-2022 all-period composite, the only climatology published for the nutrients.
All numbers on each source's native grid, weighted by its own cell volumes.

## Step 2. Mean N:P, two different averages

| Source | mean of cell ratios | ratio of means |
|---|---|---|
| WOA23 | 14.15 | 14.00 |
| CESM2 | 11.83 | 10.93 |
| GFDL-ESM4 | 14.13 | 13.65 |
| MPI-ESM1-2-LR | 12.28 | 11.78 |

The first gives every parcel of water the same say and is the headline value; the second is
pulled towards the richest water. They are different quantities, not two estimates of one.

## Step 3. NO3 = a x PO4 + b, unweighted OLS over all valid cells

| Source | a | b (umol/kg) | R2 |
|---|---|---|---|
| WOA23 | 12.90 | +1.88 | 0.908 |
| CESM2 | 4.99 | +10.52 | 0.238 |
| GFDL-ESM4 | 11.30 | +4.05 | 0.817 |
| MPI-ESM1-2-LR | 8.81 | +5.34 | 0.680 |

## Step 4. Volume inside the Redfield band (14 < N:P < 18)

| Source | volume within band |
|---|---|
| WOA23 | 55.0% |
| CESM2 | 17.2% |
| GFDL-ESM4 | 74.6% |
| MPI-ESM1-2-LR | 0.5% |

## Step 5. Basin structure and the denitrification signal

| Source | Pacific | Atlantic | Indian | N:P<12 global | N:P<12 Pacific | N:P<12 Indian |
|---|---|---|---|---|---|---|
| WOA23 | 13.57 | 15.55 | 13.73 | 6.2% | 8.1% | 10.3% |
| CESM2 | 10.58 | 13.55 | 11.18 | 36.7% | 59.4% | 35.8% |
| GFDL-ESM4 | 13.39 | 15.07 | 14.72 | 7.3% | 14.9% | 0.5% |
| MPI-ESM1-2-LR | 11.31 | 13.15 | 12.60 | 22.5% | 39.6% | 15.0% |

Marginal seas are folded into their parent ocean. The Southern and Arctic Oceans are separate basins in the WOA23 mask and are not part of these numbers.

## What the comparison shows

- Observed intermediate-water N:P is 14.15, clearly BELOW the canonical Redfield 16, and the observed regression slope is 12.90 with a positive intercept of +1.88 micromoles per kilogram. The 16:1 ratio is an uptake ratio for organic matter, not the ratio of the dissolved pool in aged intermediate water, so this is the expected result rather than an error.
- The two averages of step 2 differ for every source, which is why both are reported: the volume-weighted mean of the cell ratios gives each parcel the same say, while the ratio of the volume-weighted means is pulled towards the nutrient-rich water.
- The regression of step 3 weights cells equally and the mean of step 2 weights them by volume, so the slope and the mean ratio answer different questions and do not coincide.
- All three models place the low-N:P water in the right places - the eastern tropical Pacific, the Arabian Sea and the Bay of Bengal - so the denitrification signal is reproduced in pattern. All three overdo it in amount: 6.2 per cent of the observed 200-1000 m volume sits below N:P 12, against CESM2 36.7, GFDL-ESM4 7.3, MPI-ESM1-2-LR 22.5 per cent in the models.
- CESM2 is the outlier. Its nitrate-phosphate relationship has largely broken down: the slope is 4.99 against an observed 12.90, and R-squared falls to 0.238 against 0.908. The scatter shows why: a large branch of its water runs to high phosphate with nitrate collapsing towards zero, and 59.4 per cent of its Pacific intermediate water falls below N:P 12. This is the same excessive eastern-tropical-Pacific denitrification that case T6 finds as a deep nitrate deficit in the same model.
- GFDL-ESM4 reproduces the observed mean N:P almost exactly (14.13 against 14.15) and keeps the tightest nitrate-phosphate relationship of the three (R-squared 0.817). It holds MORE of its volume inside the Redfield band than the observations do (74.6 against 55.0 per cent), so it is not merely closer to observations but closer to Redfield than the real ocean is.
- MPI-ESM1-2-LR has almost no water in the Redfield band at all (0.5 per cent), despite a mean N:P of 12.28: its distribution is tightly clustered below 14 rather than spread across the band. A mean close to the observed one can hide a distribution that is not.
- Every source puts the Atlantic highest and the Pacific lowest, so the inter-basin ordering of N:P is reproduced by all three models even where the absolute values are not.
- Part of the spread is model structure rather than circulation. GFDL-ESM4's COBALTv2 uses fixed N:P ratios for each plankton type and MPI-ESM1-2-LR's HAMOCC6 a prescribed Redfield composition of organic matter, so in those two models dissolved N:P can depart from the uptake ratio only through the nitrogen cycle - N2 fixation, denitrification, deposition - and transport; CESM2's MARBL has variable C:P stoichiometry (Dunne et al. 2020; Mauritsen et al. 2019; Danabasoglu et al. 2020).

## Figures

- `fig1_np_scatter.png` - nitrate against phosphate for the four sources on identical axes,
  each with its own fit and the Redfield 16:1 line.
- `fig2_np_departure_maps.png` - 200-1000 m volume-weighted N:P minus 16, one shared
  diverging scale, models regridded for display only.
