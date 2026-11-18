## 1. Direct answer

For all three models, the MLD computed offline from each model's own monthly-mean `thetao` and `so` reproduces the published `mlotst` closely in pattern. Spatial correlations are 0.945–0.998. The operator MLD is systematically shallower than `mlotst`, with a negative mean difference in all six model-month cases. Agreement is best for CESM2 and worst for GFDL-ESM4. The remaining differences concentrate in deep-convection and sea-ice-edge regions.

## 2. What was actually done

- **Data:** `load_regular_data` from `cmip6_historical` for `thetao`, `so` and `mlotst`, monthly, 1971-01 to 2000-12, global. Each model was kept on its own grid, and metrics were never computed on a regridded field. If `mlotst` had a depth dimension, its first level was taken. If `so` was on a different vertical grid from `thetao`, it was interpolated to the `thetao` depths.
- **Reference depth:** each model's shallowest `thetao` level from `load_depth_levels`.
  - CESM2: 5.0 m
  - GFDL-ESM4: 2.5 m
  - MPI-ESM1-2-LR: 6.0 m
- **Areas:** `load_areas`, which returns each model's own published `areacello`. Total ocean areas are 3.605×10¹⁴ m² (CESM2), 3.611×10¹⁴ m² (GFDL-ESM4) and 3.613×10¹⁴ m² (MPI-ESM1-2-LR).
- **Operator MLD:** for each model, the 30 Marches and 30 Septembers were selected. `compute_mixed_layer_depth` was called with `method="density"`, `density_threshold=0.03` (kg m⁻³), `ref_depth` set to the model's shallowest level, and `temperature_type="auto"`. The density and salinity handling is done inside that tool; the script does not convert salinity itself. This was applied to each monthly-mean profile, so **monthly-mean T and S were used and the MLD was computed first**. The 30 MLDs were then averaged with `temporal_aggregate(freq="mean")`, a plain mean over years. `mlotst`, the time mean of instantaneous diagnostics, was averaged over the same 30 months in the same way. Values ≥1e10 (unmasked fill) were set to NaN in both fields.
- **No-crossing columns:** I take them to be NaN in the tool's output, so they drop out of the mask. The code defines the "area removed" as cells that are ocean at the surface in all 30 months but have a NaN operator MLD in at least one month. The script has no per-model handling of no-crossing encodings (such as 0 m in CESM2). The scatter step found no non-positive values among shared-mask cells (0 dropped).
- **Shared mask:** a cell is kept only if both the operator MLD and `mlotst` are non-NaN in all 30 instances of that month. This was built directly from the two fields and applied with `.where(mask)`; `build_comparison_mask` was not called. Because both fields are on the same grid within a model, no regridding was needed.
- **Metrics:**
  - `compute_difference_metrics(operator, mlotst)` gives bias and RMSE.
  - `analyze_correlation(..., weighting="area")` gives the spatial correlation.
  - `compute_exceedance_fraction(diff, threshold=20, comparison="abs_greater")` gives the >20 m share.
  - Area weights were not passed explicitly in these calls, so the weighting is handled inside the tools.
  - `locate_extreme(op="abs_max")` finds the largest difference.
  - `select_basin` with Atlantic, Pacific, Indian, Southern and Arctic gives a basin breakdown.

## 3. Key numerical results (operator MLD − `mlotst`, 1971–2000 climatology, shared mask)

| Model | Month | Mean diff (m) | RMSE (m) | Correlation | Area with \|diff\| > 20 m (% of shared area) | Shared area (% of model ocean area) |
|---|---|---|---|---|---|---|
| CESM2 | March | −4.5 | 7.5 | 0.998 | 2.7 | 96.4 |
| CESM2 | September | −2.6 | 5.1 | 0.998 | 0.2 | 97.9 |
| GFDL-ESM4 | March | −9.6 | 21.4 | 0.945 | 6.3 | 97.2 |
| GFDL-ESM4 | September | −8.0 | 14.3 | 0.987 | 7.5 | 98.5 |
| MPI-ESM1-2-LR | March | −2.5 | 6.8 | 0.995 | 1.5 | 96.8 |
| MPI-ESM1-2-LR | September | −2.4 | 14.3 | 0.983 | 3.5 | 97.6 |

**Mask sizes** (March / September, in native-grid cells):

| Model | March cells | September cells |
|---|---|---|
| CESM2 | 82,360 | 84,146 |
| GFDL-ESM4 | 42,122 | 43,259 |
| MPI-ESM1-2-LR | 32,978 | 34,894 |

The shared mask covers 96–98.5% of each model's ocean area, so the two products agree closely on where data are defined.

**Area removed by excluding no-crossing columns** (% of model ocean area, March / September):
- CESM2: 3.6 / 2.1
- GFDL-ESM4: 2.8 / 1.5
- MPI-ESM1-2-LR: 3.2 / 2.4

For example, CESM2 March removes about 1.29×10¹³ m², and the shared area plus the removed area account for essentially all of the ocean area.

**Mean fields (March / September):**

| Model | Operator MLD (m) | `mlotst` (m) |
|---|---|---|
| CESM2 | 56.1 / 71.6 | 60.6 / 74.2 |
| GFDL-ESM4 | 32.0 / 46.6 | 41.5 / 54.6 |
| MPI-ESM1-2-LR | 44.3 / 69.6 | 46.9 / 72.0 |

## 4. Ranking and where the disagreement sits

**Criterion:** the lowest mean of the March and September area-weighted RMSE, with correlation and the >20 m area share as supporting criteria.

| Rank | Model | Mean RMSE (m) | Mean correlation | Mean >20 m share (%) |
|---|---|---|---|---|
| 1 | CESM2 | 6.3 | 0.998 | 1.4 |
| 2 | MPI-ESM1-2-LR | 10.6 | 0.989 | 2.5 |
| 3 | GFDL-ESM4 | 17.9 | 0.966 | 6.9 |

All three criteria give the same order. In March alone, MPI-ESM1-2-LR (6.8 m) is marginally lower than CESM2 (7.5 m). MPI's September RMSE (14.3 m) is much larger, which is why it ranks second overall.

**Where the differences concentrate** (from `locate_extreme` and the basin tables):
- **North Atlantic and Labrador Sea in March** (deep convection):
  - CESM2 has its largest difference (−373 m) at about 59.7°N, 58.6°W. The Atlantic basin RMSE is 12.3 m, with 11% of the area above 20 m.
  - MPI's March maximum (−434 m) is at about 54°N, 51.8°W.
  - GFDL-ESM4's March maximum (+814 m) is at 54.5°N, 10.5°W. Its Atlantic RMSE is 34.3 m, with 15.6% of the area above 20 m.
- **Arctic in March:** GFDL-ESM4 has an RMSE of 25.9 m with 10.9% of the area above 20 m. MPI has an RMSE of 13.9 m with 5.4% above 20 m.
- **Southern Ocean in September:**
  - MPI has an RMSE of 37.4 m and a bias of −14.0 m, with 13.7% of the area above 20 m. Its largest difference (−698 m) is at about 72°S, 139°W.
  - GFDL-ESM4 has an RMSE of 20.0 m, with 12.1% of the area above 20 m.
  - GFDL-ESM4's largest September difference (+432 m) is at about 37.5°S, 139.5°E.
  - CESM2's Southern Ocean September RMSE is only 1.9 m.
- **GFDL-ESM4 Indian Ocean in September:** RMSE of 18.1 m, with 12.5% of the area above 20 m.

**Why exact agreement is not expected:**
- The operator MLD was computed from monthly-mean `thetao` and `so`, whereas `mlotst` is the time mean of instantaneous diagnostics. The threshold criterion is nonlinear, so the MLD of the mean profile is not the mean of the MLDs.
- Vertical interpolation and the reference depth differ.
- The equation of state may differ.
- No-crossing columns are encoded differently across models.

The consistently negative bias is plausibly linked to the monthly-mean averaging, but that was not tested separately.

## 5. Visualizations

- **Figure 1** is a 3×2 grid of difference maps, with rows CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, and columns March and September.
  - Each native-grid difference field (shared mask applied) was conservatively regridded to the `woa23` 1° longitude-latitude grid with `regrid_horizontal`, for display only.
  - The maps are drawn with `pcolormesh` on a common frame (−180 to 180°E, −90 to 90°N, equal aspect, grey background for missing cells). They are plain longitude-latitude axes, not a map projection.
  - All panels share one `RdBu_r` colour scale centred on zero, with limits set to ± the 98th percentile of the pooled |difference|. The colourbar extends at both ends.
  - Panel titles give each case's bias and RMSE, computed on the native grid.
- **Figure 2** is a 3×2 grid of scatter plots in the same layout, with operator MLD on x and `mlotst` on y.
  - Each point is one shared-mask cell on the native grid.
  - All axes are logarithmic with range 1–10⁴ m, ticks at 1, 10, 100, 1000 and 10000, and equal aspect.
  - Each panel has a dashed y = x line and an annotation of bias, RMSE and correlation.

## 6. Caveats

- The metrics are area-weighted means on the native grid of each model. Cell counts are not comparable across models because the grids differ.
- Figure 1 uses the regridded field for display only. GFDL-ESM4 is already on a regular 1° grid, so its displayed field is effectively the native one.
- No-crossing columns are excluded only to the extent that the tool returns NaN for them. The script has no model-specific treatment of 0 m or full-column-depth encodings, so any such values that fall inside the shared mask would still be counted.
