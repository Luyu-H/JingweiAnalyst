The computation succeeded for the global mean state, global trends, and global change, but the WOA23 basin-mask reader did not work in this run: `basin_woa` had `n_valid = 0`, so all Atlantic/Indian/Pacific basin changes came out as `NaN`. Those basin values therefore cannot be reported from this execution and should be recomputed after fixing the `basinmask_01.msk` parsing. The global and map results below are unaffected.

### 1. Mean-state SSS and bias against WOA23, 1971–2000

WOA23 global mean SSS over the ocean is **34.63 psu**.

| Model | Global mean SSS (psu) | Bias vs WOA23 (psu) | RMSE vs WOA23 (psu) | Spatial correlation |
|---|---:|---:|---:|---:|
| CESM2 | 34.32 | −0.32 | 1.07 | 0.874 |
| GFDL-ESM4 | 34.56 | −0.07 | 1.02 | 0.882 |
| MPI-ESM1-2-LR | 34.36 | −0.28 | 1.21 | 0.841 |

All three models are slightly **fresher than WOA23** in the global mean. GFDL-ESM4 has the smallest global-mean bias and the best RMSE/correlation. CESM2 and MPI-ESM1-2-LR are fresher by about −0.32 and −0.28 psu, respectively. The spatial correlations are all high (0.84–0.88), so the large-scale SSS pattern is reasonably represented, but the RMSE of about 1 psu indicates substantial local biases. The bias maps contain some very large local values, especially in coastal/regridded cells, even though the global-mean biases are small.

### 2. Linear trends of annual global-mean SSS

Historical period: 1950–2014. SSP245 period: 2071–2100. Slopes are in psu per decade; confidence intervals are 95%.

| Model | Historical slope | 95% CI | R² | Excludes zero? | SSP245 slope | 95% CI | R² | Excludes zero? | SSP/hist ratio | Ratio meaningful? |
|---|---:|---:|---:|---|---:|---:|---:|---|---:|---|
| CESM2 | −0.00268 | [−0.00386, −0.00149] | 0.245 | Yes | −0.0209 | [−0.0235, −0.0184] | 0.910 | Yes | 7.8 | Yes* |
| GFDL-ESM4 | −0.00110 | [−0.00215, −0.0000376] | 0.0636 | Yes | −0.00534 | [−0.00832, −0.00236] | 0.325 | Yes | 4.9 | Yes |
| MPI-ESM1-2-LR | −0.000687 | [−0.00235, 0.000977] | 0.0107 | No | −0.0229 | [−0.0289, −0.0168] | 0.684 | Yes | 33.3 | No |

\* CESM2’s ratio is statistically nonzero, but it is affected by the ensemble-member caveat described below.

Historically, all three models show weak global freshening. The trend is statistically significant for CESM2 and GFDL-ESM4, but only marginally for GFDL-ESM4 given its very low R²; for MPI-ESM1-2-LR the historical confidence interval includes zero, so its historical trend is not robust. Over 2071–2100 under SSP245, all three models show significant freshening trends, with CESM2 and MPI-ESM1-2-LR freshening much faster than GFDL-ESM4. The ratio is not meaningful for MPI-ESM1-2-LR because the historical slope is statistically indistinguishable from zero. For CESM2 the ratio should also be interpreted cautiously because the historical and SSP245 experiments use different ensemble members.

### 3. End-of-century SSS change: SSP245 2071–2100 minus historical 1995–2014

Values are area-weighted mean ΔSSS in psu, rounded to two decimals. The first global value uses the WOA comparison grid; the native-grid value is given in parentheses.

| Model | Global ΔSSS | Atlantic ΔSSS | Indian ΔSSS | Pacific ΔSSS | Caveat |
|---|---:|---:|---:|---:|---|
| CESM2 | −0.11 (−0.12 native) | NaN | NaN | NaN | historical uses r1i1p1f1; SSP245 uses r10i1p1f1 |
| GFDL-ESM4 | −0.00 (−0.01 native; unrounded −0.0049) | NaN | NaN | NaN | same member r1i1p1f1 |
| MPI-ESM1-2-LR | −0.08 (−0.09 native) | NaN | NaN | NaN | same member r1i1p1f1 |

The global ocean freshens by about **−0.11 psu for CESM2** and **−0.08 psu for MPI-ESM1-2-LR**, while **GFDL-ESM4 is essentially unchanged** (slight freshening of about −0.005 psu on the WOA grid, −0.008 psu on its native grid). The basin breakdown is missing because the basin mask was not loaded correctly. In a corrected run, the WOA23 basin mask would need to be parsed with its surface level and basin codes so that marginal seas are folded into their parent ocean and the Southern and Arctic basins are excluded.

### 4. Ensemble-member caveat

The main caveat concerns CESM2: `so`, `volcello`, and `areacello` are from **r1i1p1f1** in the historical experiment but from **r10i1p1f1** in SSP245. Therefore the CESM2 ΔSSS and the SSP/historical slope ratio conflate the scenario change with an ensemble-member difference. GFDL-ESM4 and MPI-ESM1-2-LR use r1i1p1f1 consistently in both experiments, so their changes are not affected by this particular issue.

### 5. Overall assessment

- **Mean state:** all models are slightly too fresh globally; GFDL-ESM4 is closest to WOA23, while CESM2 and MPI-ESM1-2-LR are fresher by about 0.3 psu. Spatial correlations are high but RMSE is about 1 psu, indicating regionally important SSS biases.
- **Historical trend:** global freshening is weak and not significant in MPI-ESM1-2-LR; CESM2 and GFDL-ESM4 show statistically significant but small freshening.
- **SSP245 end-century:** all models freshen over 2071–2100, with much stronger trends in CESM2 and MPI-ESM1-2-LR than in GFDL-ESM4. The global mean change by the end of the century is modest in these experiments: about −0.1 psu or less.
- **Basin means:** could not be computed in this run because the WOA23 basin mask returned no valid basin values. The global conclusions are robust, but the Atlantic/Indian/Pacific breakdown still needs to be recomputed after fixing the mask parser.
- **Figures:** Figure 1 shows the shared-scale model-minus-WOA SSS bias maps; Figure 2 shows the 1950–2014 cell-wise annual-mean SSS trend with p < 0.05 stippling; Figure 3 shows the global-mean SSS time series with the historical and SSP245 segments separated by the 2015–2070 gap and with fitted lines.
