## 1. Answer to the question

Yes — an offline mixed-layer-depth operator driven only by each model’s own monthly `thetao` and `so` reproduces the same model’s published `mlotst` to within a few metres in most of the ocean, but the level of internal consistency is model-dependent. Over the 1971–2000 March and September climatologies, all six model-month pairs show a **negative bias** (operator MLD shallower than `mlotst`) and **very high spatial correlation (r = 0.941–0.998)**, so the diagnostic is broadly consistent with the published field. However, RMSE ranges from 5.1 m (CESM2, September) to 21.2 m (GFDL-ESM4, March), and the ranking is:

**CESM2 (6.3 m mean RMSE) > MPI-ESM1-2-LR (10.6 m) > GFDL-ESM4 (17.2 m)**, from most to least internally consistent.

## 2. Methodology (exactly as executed)

- **Data and period.** CMIP6 `historical` monthly `thetao`, `so` and `mlotst` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR were read with `load_regular_data` over 1971-01 to 2000-12. Each model’s published cell areas were read with `load_areas(dataset="cmip6_historical", model=m, variable="thetao")`; the model’s total ocean area is the sum of `areacello`.
- **Operator MLD.** For each model, and separately for March (month 3) and September (month 9), the code selected the 30 corresponding months and called `compute_mixed_layer_depth` on potential temperature and salinity with `method="density"`, `density_threshold=0.03` kg m⁻³, `ref_depth=10.0` m and `temperature_type="potential"` (TEOS-10 potential density anomaly), interpolating the crossing depth between bracketing levels.
- **Common valid mask.** The operator MLD and `mlotst` were aligned (`join="inner"`), and a cell was retained only if **both** fields were non-missing in **all 30** instances of that calendar month — `(mld_op.notnull() & ml_sel.notnull()).all("time")`. The mask therefore differs by model and by month.
- **Climatology.** `temporal_aggregate(..., freq="mean")` was applied to both masked monthly fields, giving one March and one September 1971–2000 climatology per model. Provenance attributes from the native fields were carried onto the climatologies so that the weighting tools could use the native grid.
- **Statistics (native grid, area-weighted).** `compute_difference_metrics(..., weighting="area")` gave bias, RMSE, and the largest absolute error with its location (obtained from `loc_of_max` on the native grid); `analyze_correlation(..., method="pearson", weighting="area")` gave the spatial Pearson correlation; `compute_exceedance_fraction(..., threshold=20.0, comparison="abs_greater", weighting="area")` gave the area-weighted percentage of |operator − `mlotst`| > 20 m. All metrics were computed on the common-valid cells only.
- **Valid-area percentage.** The binary mask was multiplied by `areacello`, summed, and divided by the model’s total ocean area; reported as a percentage.
- **Ranking.** The mean of the March and September RMSE was used as the primary criterion (lower = more consistent).
- **Regridding** was used only for the map display (`regrid_horizontal(diff, target="woa23", method="bilinear")`); no statistic used the regridded field.

## 3. Key numerical results

| Model | Month | Bias (m) | RMSE (m) | r | \|Δ\|>20 m (% area) | Valid area (% ocean) | Max \|Δ\| (m) at (lat, lon) |
|---|---|---|---|---|---|---|---|
| CESM2 | March | −4.5 | 7.5 | 0.998 | 2.7 | 96.4 | 372.9 at (59.7, −58.6) |
| CESM2 | September | −2.6 | 5.1 | 0.998 | 0.2 | 97.9 | 51.2 at (−75.5, −56.3) |
| GFDL-ESM4 | March | −7.2 | 21.2 | 0.941 | 5.9 | 97.1 | 854.9 at (42.5, 130.5) |
| GFDL-ESM4 | September | −5.8 | 13.2 | 0.986 | 6.9 | 98.4 | 429.7 at (−37.5, 139.5) |
| MPI-ESM1-2-LR | March | −2.5 | 6.8 | 0.995 | 1.5 | 96.8 | 433.7 at (54.0, −51.8) |
| MPI-ESM1-2-LR | September | −2.4 | 14.3 | 0.983 | 3.5 | 97.6 | 697.8 at (−72.2, −139.1) |

**Interpretation.**
- The sign of every bias is negative: the offline operator systematically yields a **shallower** mixed layer than the model’s published `mlotst`, by 2.4–7.2 m on area average.
- Correlations are extremely high everywhere (≥ 0.94), so the spatial pattern of the operator MLD tracks `mlotst` well; the disagreement is a small-scale, amplitude-level issue rather than a structural one.
- CESM2 is clearly the most internally consistent: RMSE 5.1–7.5 m, and in September essentially no area exceeds the 20 m threshold (0.2 %). GFDL-ESM4 is the least consistent: RMSE 13.2–21.2 m, correlation as low as 0.941, and 5.9–6.9 % of the valid area exceeding 20 m.
- MPI-ESM1-2-LR sits in between: it has the smallest biases (−2.4 and −2.5 m) but degrades strongly in September (RMSE 14.3 m, 697.8 m maximum departure near 72°S, 139°W).
- The valid area is high for every model-month (96.4–98.4 % of the total ocean area), so the comparison is not driven by a small or unusual subset of cells.

**Where the disagreement concentrates.**
- CESM2: the largest departures are in the subpolar North Atlantic / Labrador Sea region (≈ 60° N, 59° W in March; ≈ 75° S, 56° W in September), i.e. deep-convection and high-latitude mixed-layer regimes. Exceedance is 2.7 % of area in March but only 0.2 % in September.
- GFDL-ESM4: the largest March departure is in the North Pacific near 42.5° N, 130.5° E (854.9 m), and the September maximum near 37.5° S, 139.5° E; exceedance is the largest of the three models and is slightly higher in September (6.9 %) than March (5.9 %).
- MPI-ESM1-2-LR: largest departures near 54° N, 52° W in March and near 72° S, 139° W in September, again high-latitude/subpolar mixed layers, with exceedance rising from 1.5 % (March) to 3.5 % (September).

**Why exact agreement is not expected.** The operator is an offline diagnostic from monthly `thetao`/`so` with a 0.03 kg m⁻³ potential-density increase relative to a 10 m reference, interpolated between the model’s vertical levels. The published `mlotst` is the model’s native mixed-layer diagnostic, computed online with its own reference level, its own threshold implementation, its own vertical grid and interpolation, and its own sea-ice and missing-data conventions. Monthly averaging also aliases sub-monthly mixed-layer variability, and the common-valid mask removes cells that are not defined in all 30 instances. All of these produce a small, systematically negative offset in the operator relative to `mlotst`, largest where vertical gradients and mixed-layer variability are strongest.

## 4. Visualization description

**Figure 1 (3 × 2 difference maps).** Rows are CESM2, GFDL-ESM4 and MPI-ESM1-2-LR; columns are March and September. Each panel shows operator MLD minus `mlotst` for the 1971–2000 climatology, regridded bilinearly to a common 1° longitude–latitude frame (`target="woa23"`, 180 × 360) purely for display. All six panels use a single diverging `RdBu_r` colour scale centred on zero, with the limit set to the 98th percentile of |difference| pooled across all six panels and a colourbar extended at both ends. Longitude limits are chosen according to the convention of the regridded field (0–360 or −180–180), latitude limits −80° to 90°, with one shared colourbar labelled "Operator MLD − mlotst (m)". The figure is designed to reveal where the offline operator and the native diagnostic diverge, with the shared zero-centred scale making the sign of the disagreement comparable across models and seasons.

**Figure 2 (3 × 2 scatter plots).** Same layout and ordering. Each panel scatters operator MLD (x) against published `mlotst` (y), one point per common-valid ocean cell, on identical logarithmic axes and identical major tick locations across all panels (ticks derived from the 0.5th and 99.5th percentiles of the pooled positive values, rounded to decades). A dashed black y = x reference line marks perfect agreement, and each panel is annotated with its area-weighted bias, RMSE and Pearson correlation. The figure is designed to show whether departures from the 1:1 line are systematic (offset/bias) or scatter, and whether they grow at large mixed-layer depths.

## 5. Insights

- The diagnostic is self-consistent at the pattern level for all three models (r ≥ 0.94), so an offline density-threshold operator is a reasonable proxy for `mlotst`, but not a substitute at the metre level.
- Bias is uniformly negative and small (≤ 7.2 m), consistent with the operator using a fixed 10 m reference and linear interpolation on a coarse monthly grid, which tends to place the threshold crossing slightly shallower than the native online diagnostic.
- Consistency is not simply a function of horizontal resolution or cell count: CESM2 (82 360/84 146 valid cells) is most consistent, while GFDL-ESM4 (~42 000 cells) is least consistent, and MPI-ESM1-2-LR has the fewest valid cells (~33 000–35 000) yet ranks second.
- The seasonal contrast is model-specific: CESM2 and MPI-ESM1-2-LR degrade in September (Southern Hemisphere winter/spring deep mixed layers), whereas GFDL-ESM4 is worst in March. Disagreement is consistently concentrated in deep-convection and high-latitude/subpolar regions, exactly where mixed-layer depth is large and vertical resolution and threshold placement matter most.
- Valid-area fractions of 96–98 % indicate the comparison is robust and not an artefact of aggressive masking.
