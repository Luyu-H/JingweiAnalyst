## Direct answer

The offline operator MLD does **not** reproduce the published `mlotst` exactly in any of the six model/month cases, but the level of internal consistency is high for CESM2, intermediate for MPI-ESM1-2-LR, and poorest for GFDL-ESM4. All six cases have a negative area-weighted mean difference, meaning the operator MLD is generally **shallower** than the published `mlotst`. Spatial correlations are high (0.983–0.998), so the main spatial pattern is captured; the disagreements are mainly in bias, root-mean-square difference, and localized high-latitude/midlatitude regions.

## Methodology actually used

- CMIP6 historical monthly `thetao`, `so`, `mlotst`, and fixed `areacello` were opened for CESM2, GFDL-ESM4, and MPI-ESM1-2-LR, all variant `r1i1p1f1`, for **1971-01 through 2000-12**.
- CESM2 area was additionally masked with `sftof > 0`.
- For each model and each of the 30 Marches and 30 Septembers, the operator MLD was computed from that month’s `thetao`/`so` profile using:
  - pressure from model depth and latitude;
  - practical salinity converted to Absolute Salinity;
  - potential temperature converted to Conservative Temperature;
  - potential density anomaly `sigma0`;
  - shallowest model level as reference;
  - threshold `sigma0 - sigma0_ref = 0.03 kg m-3`;
  - linear interpolation of the crossing depth between bracketing levels;
  - missing operator MLD when the threshold was never reached.
- The operator MLD was computed on each monthly profile first, then averaged over the 30 instances. It was not computed from the 30-year mean temperature and salinity.
- For each model/month, the all-time-valid mask required both operator MLD and published `mlotst` to be finite in **all 30** instances of that calendar month. The climatologies and comparisons were made on that common mask.
- Metrics were computed on each model’s native grid with `areacello` area weights. The comparison mask was built with `build_comparison_mask(..., ocean_only=False)` and applied with `.where(mask)`.
- Area-weighted Pearson correlation used `areacello` weights. The exceedance percentage is the area with `|D| > 20 m` divided by common-mask area. Coverage is common-mask area divided by total model ocean area.
- No required model variable was missing; the computation succeeded.

## Numerical results

Here `D = operator MLD - published mlotst`; negative values mean the operator MLD is shallower.

| Model | Month | Mean D (m) | RMSE (m) | Correlation | Area with \|D\|>20 m (%) | Common-mask coverage (%) | Common-mask cells |
|---|---:|---:|---:|---:|---:|---:|---:|
| CESM2 | March | -4.5 | 7.5 | 0.998 | 2.7 | 96.4 | 82,360 |
| CESM2 | September | -2.6 | 5.1 | 0.998 | 0.2 | 97.9 | 84,146 |
| GFDL-ESM4 | March | -9.6 | 21.4 | 0.945 | 6.3 | 97.2 | 42,122 |
| GFDL-ESM4 | September | -8.0 | 14.3 | 0.987 | 7.5 | 98.5 | 43,259 |
| MPI-ESM1-2-LR | March | -2.5 | 6.8 | 0.995 | 1.6 | 96.8 | 33,260 |
| MPI-ESM1-2-LR | September | -2.4 | 14.3 | 0.983 | 3.5 | 97.6 | 35,179 |

### Ranking of internal consistency

Primary criterion: average of the March and September area-weighted RMSE; lower RMSE means more internally consistent.

| Rank | Model | Mean March/September RMSE (m) | Mean absolute bias (m) | Mean \|D\|>20 m area (%) |
|---:|---|---:|---:|---:|
| 1 | CESM2 | 6.3 | 3.5 | 1.4 |
| 2 | MPI-ESM1-2-LR | 10.5 | 2.5 | 2.5 |
| 3 | GFDL-ESM4 | 17.9 | 8.8 | 6.9 |

CESM2 is the most internally consistent by RMSE and by the exceedance area, although MPI-ESM1-2-LR has a slightly lower mean absolute bias. GFDL-ESM4 is clearly the least consistent, with the largest RMSE, largest absolute bias, and largest area exceeding 20 m.

## Regional concentration of disagreement

The code summarized disagreement by latitude bands using area-weighted means and RMSEs.

- **CESM2 March:** largest disagreement is in the Northern Hemisphere subtropics/midlatitudes, especially 20–40°N (mean -11.9 m, RMSE 14.2 m) and 40–60°N (RMSE 12.7 m). September disagreement is weaker and more Southern Hemisphere/tropical, with the largest band RMSE at 40–20°S (mean -6.5 m, RMSE 8.6 m).
- **GFDL-ESM4 March:** strong disagreement in the Northern Hemisphere, especially 20–40°N (mean -16.9 m, RMSE 24.1 m), 40–60°N (RMSE 60.7 m), and 60–90°N (RMSE 34.6 m). September disagreement is larger in the high southern latitudes and Southern Hemisphere subtropics, e.g. 90–60°S (mean -8.7 m, RMSE 21.1 m), 60–40°S (RMSE 22.2 m), and 40–20°S (mean -13.5 m, RMSE 18.7 m).
- **MPI-ESM1-2-LR March:** disagreement concentrates in the Northern Hemisphere mid-to-high latitudes, with RMSE 11.3 m at 20–40°N, 14.0 m at 40–60°N, and 12.9 m at 60–90°N. September is dominated by the far Southern Ocean: 90–60°S has mean -22.5 m and RMSE 52.1 m, with additional Southern Hemisphere midlatitude disagreement.

## Non-crossing excluded-area report

The code reports the area where at least one of the 30 monthly profiles had a valid surface but never reached the 0.03 kg m-3 threshold, and therefore produced missing operator MLD.

| Model | Month | Non-crossing area (% of total ocean) | Native cells |
|---|---:|---:|---:|
| CESM2 | March | 3.6 | 3,736 |
| CESM2 | September | 2.1 | 1,950 |
| GFDL-ESM4 | March | 2.8 | 2,442 |
| GFDL-ESM4 | September | 1.5 | 1,305 |
| MPI-ESM1-2-LR | March | 3.2 | 3,531 |
| MPI-ESM1-2-LR | September | 2.4 | 1,612 |

These excluded non-crossing areas differ by model and month, and the common-mask coverage also differs. The ranking should therefore not be read as a comparison over identical ocean areas.

## Why the operator and `mlotst` are not expected to agree exactly

Even though both quantities come from the same simulation, the offline operator differs from the online published `mlotst` diagnostic because:
- the reference level and detailed density criterion may differ;
- the operator uses a `sigma0` difference from the shallowest model level, whereas the model diagnostic may use a different density formulation or a displaced-parcel criterion at local pressure;
- the operator is computed on individual monthly profiles and then averaged, so it is not the same as the model’s online time-mean MLD;
- this time-averaging order tends to make the operator MLD shallower, consistent with the negative mean differences found here;
- non-crossing columns are set to missing, and their area differs by model/month;
- `thetao` masks, `so` masks, and `mlotst` masks can differ, so a common all-time-valid mask is required.

These factors are diagnostic properties, not model physics errors.

## Visualization description

Two figures were generated.

**Figure 1** is a 3-by-2 panel of maps of `D = operator MLD - mlotst`, with rows CESM2, GFDL-ESM4, and MPI-ESM1-2-LR and columns March and September. Each native curvilinear model field was conservatively remapped with area weights to a regular 2° longitude-latitude grid for display only. All six panels use a common geographic frame, a single shared diverging colour scale centred on zero with symmetric limits clipped at the pooled 98th percentile of `|D|`, coastlines/ocean-mask contours, and one shared colorbar. All numerical metrics were computed on the native grids, not on the display grid.

**Figure 2** is a 3-by-2 scatter panel in the same model/month layout, with operator MLD on the x-axis and published `mlotst` on the y-axis, one point per common-mask ocean cell. Both axes are logarithmic with identical limits and tick locations in every panel, and a dashed `y = x` reference line is drawn. Each panel is annotated with its area-weighted bias, RMSE, and correlation. Cells with zero or negative values cannot appear on logarithmic axes; those cells were omitted from the scatter display only, with the omitted count annotated, while all metrics used the full common mask.
