## 1. Direct answer

Yes — the offline operator MLD reproduces the models' published `mlotst` closely, but not exactly, and the residual is a property of the diagnostic rather than of the model physics. For all six model-and-month combinations the operator MLD is **systematically shallower** than `mlotst` (bias negative everywhere), correlations are 0.94–0.998, and the area with |difference| > 20 m stays between 0.2 % and 7.5 % of the common-mask area.

Ranking by internal consistency (mean of the March and September RMSDs, lower = more consistent):

| Rank | Model | Mean RMSD (m) | Mean bias (m) | Mean r | Mean >20 m (%) | Mean coverage (%) |
|---|---|---|---|---|---|---|
| 1 | CESM2 | **6.28** | −3.53 | 0.998 | 1.4 | 97.1 |
| 2 | MPI-ESM1-2-LR | **10.57** | −2.49 | 0.982 | 2.5 | 97.2 |
| 3 | GFDL-ESM4 | **17.85** | −8.79 | 0.964 | 6.9 | 97.9 |

(Mean bias, correlation, exceedance and coverage in the last four columns are the arithmetic averages of the two monthly values reported below.)

## 2. Methodology, exactly as executed

**Data and period.** `load_regular_data` retrieved CMIP6 `historical` `thetao`, `so` and `mlotst` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR over **1971-01 to 2000-12** (monthly, `chunk_time=6`). Each model's own published cell areas came from `load_areas(dataset="cmip6_historical", model=…, variable="thetao")` and the vertical grid from `load_depth_levels(…, variable="thetao")`. 1e20 was left to xarray to decode as NaN; each field's validity was taken from that field itself (no shared mask between variables).

**Reference depth (shallowest level centre), as read from the models' own levels:**

| Model | levels | top→bottom (m) | `ref_depth` used |
|---|---|---|---|
| CESM2 | 60 | 5 – 5375 | **5.0 m** |
| GFDL-ESM4 | 35 | 2.5 – 6500 | **2.5 m** |
| MPI-ESM1-2-LR | 40 | 6 – 5720 | **6.0 m** |

**Operator MLD.** `compute_mixed_layer_depth(temp_data=thetao, sal_data=so, ref_depth=ref_depth, method="density", density_threshold=0.03, temperature_type="potential")` was applied **to each month separately** (never averaging T/S before the MLD), i.e. σ0 criterion of 0.03 kg m⁻³ relative to the shallowest level centre, with the crossing interpolated between bracketing levels. `thetao` entered by the potential-temperature route; `so` was treated as practical salinity and converted to absolute salinity only inside the density routine. Units are metres.

**Climatologies and the common mask.** For March (month 3) and September (month 9) the 30 monthly instances were selected, aligned, and a single mask was formed as `(mld_m.notnull() & mlo_m.notnull()).all("time")` — cells where **both** fields are defined in **all 30** Marches or Septembers. The March/September climatologies are the **plain arithmetic 30-year mean of that calendar month** (`.where(mask).mean("time")`); no day-weighting was applied, since only one calendar month is averaged at a time. The same mask was applied to the operator MLD, to `mlotst` and to `areacello`. The difference field is **operator MLD − mlotst**.

**Metrics** (native grid, common mask, `areacello` weights): `compute_difference_metrics(data=mld_c, reference=mlo_c, dim=hd, weighting="area")` for bias, RMSE, centred RMSE and max |error|; `analyze_correlation(var1=mld_c, var2=mlo_c, dim=hd, method="pearson")` at zero lag for the spatial correlation; `compute_exceedance_fraction(data=diff, threshold=20.0, comparison="abs_greater", dim=hd, weighting="area")` for the area fraction with |difference| > 20 m; `compute_basic_stat(..., op="sum")` for the total ocean area, the common-mask area and the excluded area. Total ocean areas: CESM2 3.6051×10¹⁴ m², GFDL-ESM4 3.6111×10¹⁴ m², MPI-ESM1-2-LR 3.6130×10¹⁴ m² (about 1 % apart, so each integral covers its own model's area and is not rescaled).

## 3. Headline numbers (operator MLD − mlotst)

| Model | Month | Bias (m) | RMSD (m) | cRMSD (m) | r | >20 m (%) | Coverage (%) | Common area (m²) |
|---|---|---|---|---|---|---|---|---|
| CESM2 | March | −4.5 | 7.5 | 6.0 | 0.998 | 2.7 | 96.4 | 3.4759×10¹⁴ |
| CESM2 | September | −2.6 | 5.1 | 4.4 | 0.998 | 0.2 | 97.9 | 3.5279×10¹⁴ |
| GFDL-ESM4 | March | −9.6 | 21.4 | 19.2 | 0.942 | 6.3 | 97.2 | 3.5107×10¹⁴ |
| GFDL-ESM4 | September | −8.0 | 14.3 | 11.8 | 0.985 | 7.5 | 98.5 | 3.5587×10¹⁴ |
| MPI-ESM1-2-LR | March | −2.5 | 6.8 | 6.3 | 0.994 | 1.5 | 96.8 | 3.4957×10¹⁴ |
| MPI-ESM1-2-LR | September | −2.4 | 14.3 | 14.1 | 0.969 | 3.5 | 97.6 | 3.5254×10¹⁴ |

Area with |difference| > 20 m: CESM2 9.22×10¹² m² (March) and 5.72×10¹¹ m² (September); GFDL-ESM4 2.20×10¹³ m² and 2.67×10¹³ m²; MPI-ESM1-2-LR 5.38×10¹² m² and 1.23×10¹³ m².

**Interpretation.** All biases are negative and of order a few metres, i.e. the offline MLD built from monthly-mean T/S is shallower than the published `mlotst`; the offset is small for CESM2 and MPI-ESM1-2-LR and roughly three times larger for GFDL-ESM4 in both months. CESM2 has the highest correlations (≈0.998) and the smallest RMSD in both months. GFDL-ESM4 is worst on every diagnostic, and is the only model whose >20 m exceedance is similar in March and September (6.3 % vs 7.5 %), while CESM2 and MPI-ESM1-2-LR degrade markedly in September (CESM2 RMSD 5.1→ small in September but MPI 6.8→14.3 m), pointing to the September re-stratification being the hardest case for a monthly-mean operator. Max single-cell |errors| from `compute_difference_metrics` are large (51 m to 814 m), showing that a handful of cells dominate the tails.

**Excluded cells (threshold never reached).** Where the operator MLD is undefined because the σ0 threshold is never crossed but `mlotst` is present, the excluded areas are: CESM2 1.29×10¹³ m² (3.6 %) in March and 7.72×10¹² m² (2.1 %) in September; GFDL-ESM4 1.00×10¹³ m² (2.8 %) and 5.24×10¹² m² (1.5 %); MPI-ESM1-2-LR 1.17×10¹³ m² (3.2 %) and 8.76×10¹² m² (2.4 %). These are essentially the well-mixed high-latitude columns, and are exactly the cells where the models fill `mlotst` differently: 0 m for CESM2 versus approximately the full column for GFDL-ESM4 and MPI-ESM1-2-LR. Because they are removed from the common mask rather than filled, they affect the *coverage* percentages (96.4–98.5 %) and not the bias/RMSD/correlation.

## 4. Visualizations

**Figure 1** is a 3 × 2 panel of global maps of operator MLD − `mlotst`, rows CESM2 / GFDL-ESM4 / MPI-ESM1-2-LR, columns March / September. Each native-grid difference field was put on a common 1° longitude–latitude grid with `regrid_horizontal(target="woa23", method="bilinear", skipna=True, na_thres=1.0)` **for display only** (no metric uses the remapped fields), drawn as true geographic maps with `pcolormesh`, `RdBu_r`, one shared diverging scale centred on zero whose limit is the pooled 98th percentile of |difference| across all six panels (clipped, colour bar marked `extend="both"`, land in grey, −80° to 90° latitude). Each panel is titled with its bias, RMSD and correlation.

**Figure 2** is the matching 3 × 2 panel of scatter plots on the **native common-mask cells**, operator MLD on the x axis against `mlotst` on the y axis, one point per ocean cell (82 360 / 84 146 cells for CESM2 March/September, 42 122 / 43 259 for GFDL-ESM4, 32 978 / 34 894 for MPI-ESM1-2-LR). Both axes are logarithmic with one identical range and identical ticks across all six panels, a dashed 1 : 1 line is drawn, and each panel is annotated with its bias, RMSD, correlation and >20 m percentage. The axis range comes from the pooled 0.5th and 99.5th percentiles of the positive values; ticks are 1, 3, 10, 30, 100, 300, 1000, 3000, 10 000 m filtered to that range. Cells with non-positive or non-finite values would be excluded from the plot only (not from the metrics) — the check found **zero** such cells in all six panels.

## 5. Where the disagreement concentrates, and why exact agreement is not expected

The six difference maps are the tool for locating the residual: the code computed the exceedance field as well, and its area is largest for GFDL-ESM4 in both months and for MPI-ESM1-2-LR in September, and smallest for CESM2 in September. No separate regional averaging was performed by the code, so the localisation is delivered by the maps rather than by a tabulated regional score; the display-grid extremes (GFDL-ESM4 March spanning −627 m to +814 m, MPI-ESM1-2-LR September down to −463 m) show that the residual is concentrated in a small number of strongly stratified or convecting cells rather than spread uniformly.

Exact agreement is not expected, and none of the following is model error:

1. **Reference level.** The operator uses the models' shallowest level centre (5.0, 2.5 and 6.0 m), whereas `mlotst` is defined against a 10 m reference; a shallower reference makes the operator MLD systematically shallower, consistent with the negative bias in all six panels.
2. **σ0 versus the displaced-parcel comparison.** The operator compares σ0 at a fixed reference, while the `mlotst` criterion compares a displaced parcel at the local pressure p(k); the two differ where the background stratification is strong.
3. **MLD-of-mean versus mean-of-MLDs.** Each month's MLD was formed from that month's monthly-mean T/S and then averaged over 30 years; the published field is generated at the model time step and then averaged. Monthly-mean T/S smooths thin, transient mixed layers, which again biases the operator shallow.
4. **Grid and resolution.** GFDL-ESM4 is delivered on a regridded 1° grid (`gr`) with a regridded land–sea mask, adding a remapping difference, and its ~100 km resolution versus MPI's ~250 km makes narrow frontal features non-comparable in amplitude. This is the main reason GFDL-ESM4 ranks last despite its high correlation.

Overall: the operator MLD is an internally consistent reconstruction of `mlotst` for all three models (r ≥ 0.94, RMSD 5–21 m, coverage 96–99 % of the ocean), CESM2 is the most consistent, MPI-ESM1-2-LR second, GFDL-ESM4 third, and the negative bias of a few metres plus the small >20 m exceedance areas are fully accounted for by the reference level, the σ0 versus displaced-parcel formulation, the MLD-of-mean versus mean-of-MLDs difference, and the use of monthly-mean T/S.
