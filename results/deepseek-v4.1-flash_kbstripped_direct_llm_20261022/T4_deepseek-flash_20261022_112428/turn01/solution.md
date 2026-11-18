The offline operator MLD reproduces the published `mlotst` well in the gross spatial pattern—spatial correlations are \(\ge 0.947\) in all six cases—but it does not reproduce it exactly. The comparison is an internal consistency check of the diagnostic; it is not a test of the model physics.

Area-weighted statistics for the 1971–2000 March/September climatologies, using only cells where both fields are defined in all 30 instances of that month:

| Model | Month | Mean diff (operator − mlotst, m) | RMSE (m) | Correlation | Area with \|diff\| > 20 m (%) | Common area (% of total ocean) |
|---|---:|---:|---:|---:|---:|---:|
| CESM2 | March | -2.1 | 5.6 | 0.998 | 1.0 | 96.4 |
| CESM2 | September | -0.3 | 4.0 | 0.998 | 0.0 | 97.8 |
| GFDL-ESM4 | March | -6.9 | 20.2 | 0.947 | 5.8 | 97.2 |
| GFDL-ESM4 | September | -5.7 | 13.0 | 0.987 | 6.9 | 98.5 |
| MPI-ESM1-2-LR | March | -0.7 | 6.1 | 0.996 | 1.3 | 96.7 |
| MPI-ESM1-2-LR | September | -0.3 | 12.1 | 0.988 | 3.4 | 97.5 |

The sign convention is operator MLD minus published `mlotst`, so negative values mean the offline operator is, on average, shallower than the published mixed-layer depth. The common-area coverage is high, about 96–98% of each model’s total ocean area, so the statistics are representative of nearly the full ocean.

**Ranking from most to least internally consistent**, using the lowest mean area-weighted RMSE over March and September as the criterion:

1. **CESM2** — mean RMSE \(4.8\) m
2. **MPI-ESM1-2-LR** — mean RMSE \(9.1\) m
3. **GFDL-ESM4** — mean RMSE \(16.6\) m

CESM2 is clearly the most self-consistent: its correlations are essentially 0.998, its September bias is negligible, and only about 1% or less of the common area has \(|diff|>20\) m. MPI-ESM1-2-LR is intermediate: small mean biases and high correlations, but larger RMSE in September and locally larger residuals. GFDL-ESM4 is the least internally consistent: it has the largest negative bias in both months (about -6 to -7 m), the lowest March correlation, the largest RMSE, and the largest area fraction exceeding 20 m (about 6–7%).

The disagreement is not globally uniform. It concentrates in regions where the mixed layer is deep or the vertical grid is comparatively coarse: high-latitude deep-convection regions such as the North Atlantic subpolar/Labrador/Nordic seas, the Southern Ocean and Antarctic marginal seas, and some marginal or coastal seas. In those regions, small differences in the criterion or in vertical interpolation can produce large depth differences. In the open subtropical and tropical ocean, the two fields agree much more closely, which explains the very high correlations despite local large residuals.

The operator and `mlotst` are not expected to agree exactly even when both are based on the same nominal density criterion. Possible reasons include differences in the reference depth and density threshold implementation, use of potential vs practical temperature/salinity and different density formulae, online versus offline computation, vertical interpolation and remapping, treatment of the surface layer and 10 m reference, masking of shallow or marginal cells, and the loss of sub-monthly nonlinear information when computing the operator from monthly mean `thetao` and `so`. Thus the residual disagreement measures the sensitivity of the mixed-layer diagnostic, not an inconsistency in the model physics.
