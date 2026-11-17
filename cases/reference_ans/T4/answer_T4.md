# Reference answer - T4: MLD operator against each model's own mlotst

1971-2000 March and September climatologies. Potential density threshold 0.03 kg/m3
referenced to each model's shallowest level. MLD computed month by month and then averaged.
Every number on the model's native grid, area-weighted by its own areacello.

## Steps 3-5. Bias, error, correlation and area of disagreement

| Model | Month | bias (m) | RMSD (m) | correlation | area with abs diff > 20 m | shared area |
|---|---|---|---|---|---|---|
| CESM2 | March | -4.5 | 7.5 | 0.998 | 2.7% | 96.4% of ocean |
| CESM2 | September | -2.6 | 5.1 | 0.998 | 0.2% | 97.9% of ocean |
| GFDL-ESM4 | March | -9.6 | 21.4 | 0.945 | 6.3% | 97.2% of ocean |
| GFDL-ESM4 | September | -8.0 | 14.3 | 0.987 | 7.5% | 98.5% of ocean |
| MPI-ESM1-2-LR | March | -2.5 | 6.8 | 0.995 | 1.5% | 96.8% of ocean |
| MPI-ESM1-2-LR | September | -2.4 | 14.3 | 0.983 | 3.5% | 97.6% of ocean |

The bias is the area-weighted mean of (operator MLD - mlotst); a negative value means the
operator puts the mixed layer shallower than the model's own diagnostic. The shared-area
column is the ocean area where both fields are defined in all 30 years, as a percentage of
the model's total ocean area.

## Step 6. Consistency ranking

| Model | mean RMSD over both months (m) |
|---|---|
| CESM2 | 6.3 |
| MPI-ESM1-2-LR | 10.6 |
| GFDL-ESM4 | 17.9 |

Most internally consistent: **CESM2**. Least: **GFDL-ESM4**.

## Why the two are not expected to agree exactly

- Both fields come from the same simulation, so the residual is a property of the two diagnostics. The protocol computes mlotst online at every time step, from a surface parcel moved to each level's own pressure, and archives the monthly mean of those instantaneous depths (Griffies et al. 2016, H24). The operator takes a sigma0 difference on the monthly-mean profile. A monthly mean is itself an averaged profile, so the smoothing that makes a 30-year-mean profile too shallow also makes the monthly one shallower than the mean of the instantaneous depths: the bias is expected to be negative, and it is negative in 6 of the 6 model-and-month combinations.
- GFDL-ESM4 computes mlotst on its native 75 hybrid layers with a 2 m top layer, but its archived thetao and so were remapped to 35 z levels and regridded to 1 degree, so the operator sees a coarser profile than the model's own diagnostic did. Part of its larger RMSD is therefore the archive's resolution, and the ranking measures consistency between the archived fields, not the quality of the model's mixed layer.
- Threshold check, 1971-03: RMSD against mlotst at 0.03 and at 0.125 kg/m3 is CESM2 17.7 and 86.3 m; GFDL-ESM4 38.6 and 90.1 m; MPI-ESM1-2-LR 22.7 and 60.3 m. The protocol's 0.03 kg/m3 is the one mlotst is consistent with.

## Figures

- `fig1_difference_maps.png` - 3x2 maps of the difference field, one shared diverging scale
  (+/-30 m), regridded to 1 degree for display only.
- `fig2_scatter_panels.png` - 3x2 scatter of operator MLD against mlotst on the native grid,
  logarithmic axes with identical ranges, dashed y = x.
