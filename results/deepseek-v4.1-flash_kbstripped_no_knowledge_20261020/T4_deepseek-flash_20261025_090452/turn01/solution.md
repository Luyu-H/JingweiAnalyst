## Direct answer

The offline operator MLD does **largely reproduce** each model’s published `mlotst` in spatial pattern, with very high spatial correlations (0.937–0.998), but it does **not** agree exactly. The operator MLD is systematically **shallower** than the published `mlotst` in every model/month combination, with area-weighted biases from about **−2.4 to −7.2 m** and RMSDs from **5.1 to 21.2 m**. The largest local disagreements occur in isolated high-latitude/deep-convection and marginal-sea cells.

Ranking of internal consistency, using the mean of the March and September area-weighted RMSD over the common mask as the stated criterion, from most to least consistent:

1. **CESM2** — mean RMSD ≈ **6.3 m**
2. **MPI-ESM1-2-LR** — mean RMSD ≈ **10.6 m**
3. **GFDL-ESM4** — mean RMSD ≈ **17.2 m**

## Methodology actually used

- Loaded CMIP6 historical monthly `thetao`, `so`, and `mlotst` on Omon, plus fixed `areacello`, for CESM2, GFDL-ESM4, and MPI-ESM1-2-LR over **1971-01 to 2000-12** using `load_regular_data`.
- Selected **March** and **September** separately, giving 30 monthly instances of each month.
- For each model and month, computed the offline operator MLD by calling `compute_mixed_layer_depth` directly on the model’s `thetao` and `so`, with:
  - `ref_depth = 10.0 m`
  - `method = "density"`
  - `density_threshold = 0.03 kg m⁻³`
  - crossing interpolated between bracketing levels by the tool.
- A common mask was built directly on each model’s native grid as the conjunction of non-null operator MLD and non-null published `mlotst` in **all 30 instances** of that month.
- Applied the mask to the monthly series, then used `temporal_aggregate(... freq="mean")` to obtain the 1971–2000 March and September climatologies.
- Computed the difference field as **operator MLD minus published mlotst** from the two climatologies.
- Metrics:
  - bias, RMSD, MAE, centered RMSE, and maximum absolute error via `compute_difference_metrics(... weighting="area")`;
  - spatial Pearson correlation via `analyze_correlation(... method="pearson")`;
  - percentage of common-mask area with `|difference| > 20 m` via `compute_exceedance_fraction(... threshold=20.0, comparison="abs_greater", weighting="area")`;
  - locations of largest absolute departures via `locate_extreme(... op="abs_max")`;
  - coverage via `areacello`, as common-mask area divided by the area where `mlotst` is ever defined, times 100.
- No lag analysis was performed.
- For display only, difference fields were regridded with `regrid_horizontal(... target="woa23", method="bilinear")`; all metrics were computed on native cells.

## Key numerical results

| Model | Month | Bias, operator − mlotst (m) | RMSD (m) | Spatial r | Area with \|diff\| > 20 m (%) | Common-mask coverage (% of model ocean area) |
|---|---:|---:|---:|---:|---:|---:|
| CESM2 | March | −4.5 | 7.5 | 0.998 | 2.7 | 96.4 |
| CESM2 | September | −2.6 | 5.1 | 0.998 | 0.2 | 97.9 |
| GFDL-ESM4 | March | −7.2 | 21.2 | 0.937 | 5.9 | 97.1 |
| GFDL-ESM4 | September | −5.8 | 13.2 | 0.984 | 6.9 | 98.4 |
| MPI-ESM1-2-LR | March | −2.5 | 6.8 | 0.994 | 1.5 | 96.8 |
| MPI-ESM1-2-LR | September | −2.4 | 14.3 | 0.969 | 3.5 | 97.6 |

All six cases have negative bias, meaning the offline operator MLD is on average shallower than the published `mlotst`. The common mask covers more than 96% of the model ocean area in all cases, so the comparison is nearly global.

Largest absolute departures:

- **CESM2 March:** −372.9 m near 59.71°N, 58.59°W.
- **CESM2 September:** −51.2 m near 75.48°S, 56.31°W.
- **GFDL-ESM4 March:** +854.9 m near 42.5°N, 130.5°E.
- **GFDL-ESM4 September:** +429.7 m near 37.5°S, 139.5°E.
- **MPI-ESM1-2-LR March:** −433.7 m near 53.96°N, 51.79°W.
- **MPI-ESM1-2-LR September:** −697.8 m near 72.18°S, 139.08°W.

These extremes indicate that disagreement is concentrated in localized high-latitude/deep-convection and marginal-sea regions rather than being uniform over the ocean. GFDL-ESM4 has the broadest area exceeding the 20 m threshold and the largest RMSD, while CESM2 has the smallest overall RMSD and the highest correlations.

## Interpretation

The high spatial correlations show that the offline diagnostic captures the broad structure of the published mixed layer depth well. The systematic negative bias and localized large differences show that it is not a bitwise reproduction of the online `mlotst` diagnostic. This is expected: `mlotst` is computed online within each model on its native vertical grid and algorithm, whereas the operator MLD is computed offline from archived `thetao` and `so`, with interpolation between levels, a 10 m reference depth, a 0.03 kg m⁻³ density threshold, and possible differences in surface/reference-level treatment, vertical discretisation, density-criterion implementation, and post-processing/archiving.

## Figures

- **Figure 1** shows a 3-by-2 panel of global maps of operator MLD minus published `mlotst`, with rows CESM2, GFDL-ESM4, and MPI-ESM1-2-LR and columns March and September. The fields were regridded to a regular longitude–latitude display grid only. All panels use one shared diverging `RdBu_r` colour scale centred on zero, clipped for display at the pooled 98th percentile of absolute difference. Each panel is annotated with its bias, RMSD, and correlation.
- **Figure 2** shows a 3-by-2 panel of scatter plots in the same layout, operator MLD on the x-axis and published `mlotst` on the y-axis, one point per common-mask ocean cell on the native grid. The axes are logarithmic with identical ranges and ticks in every panel, and each panel includes a dashed y = x reference line and annotations of bias, RMSD, and correlation. Only finite positive values were used for the logarithmic axes.

Both figures were generated successfully.
