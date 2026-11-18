## 1. Direct answer

Yes — the offline operator reproduces each model's published `mlotst` closely in the spatial sense (area-weighted Pearson r = 0.945–0.998 for all six model-month combinations), but it is **systematically shallower**: every one of the six bias values (operator MLD − `mlotst`) is negative, between −2.4 m and −9.6 m. The disagreement is therefore a property of the diagnostic construction, not of model physics, and it behaves as the plan anticipated: the operator (shallowest-level reference, displaced-parcel Δρ = 0.03 kg m⁻³ criterion at local pressure, mean of 30 per-profile MLDs) sits above the published field in every case.

## 2. Methodology actually executed

- **Data**: CMIP6 `historical` monthly `thetao`, `so`, `mlotst` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, **1971-01 to 2000-12** (360 monthly records each). Grids as published: CESM2 `gn` (nlat 384, nlon 320), GFDL-ESM4 `gr` (regridded 1°, lat 180 × lon 360, 35 z levels), MPI-ESM1-2-LR `gn` (j 220, i 254). All metrics kept on the native grid with each model's own `areacello`.
- **Geometry**: `load_depth_levels` gave each model's levels/bounds; the reference depth used for the surface parcel is the **shallowest level**: CESM2 5.0 m, GFDL-ESM4 2.5 m, MPI-ESM1-2-LR 6.0 m.
- **Operator MLD**: `compute_mixed_layer_depth` with `method='density'`, `density_threshold=0.03` (δB_crit = 0.0003 m s⁻²), `temperature_type='auto'` (so `thetao` was used as potential temperature), computed **per monthly profile before any averaging**; the first threshold crossing was interpolated between bracketing levels. Non-crossing columns were excluded (counted and reported separately).
- **Averaging**: `temporal_aggregate(..., freq='climatology', op='mean')` was applied to the full 12-month field, and the March / September member was then selected. Averaging was therefore over the **30 monthly means of that calendar month** (no day weighting, no averaging of profiles), i.e. exactly the "mean of MLDs," not the "MLD of a mean."
- **Comparison mask**: one shared mask per model and month, `(op.notnull() & mlotst.notnull()).all(dim='time')`, retaining cells valid in **all 30** Marches or Septembers; applied with `.where(mask)` before any metric. Coverage is reported as retained `areacello` area divided by the model's total ocean area.
- **Metrics**: `compute_difference_metrics(dim=None)` (area-weighted by default) for bias and RMSD; `analyze_correlation(dim=['lat','lon'], method='pearson', weighting='area')`; `compute_exceedance_fraction(threshold=20, comparison='abs_greater', weighting='area')` for the >20 m share (area-based, missing cells out of the denominator).
- **Display**: Figure 1 fields were regridded with `regrid_horizontal(target='1deg', method='conservative', skipna=True)` — display only.

## 3. Key numerical results

| Model (grid) | Month | Bias (m) | RMSD (m) | r | \|Δ\|>20 m (% area) | Covered (% ocean) | Retained cells | Non-crossing area (any / all, 10⁶ km²) |
|---|---|---|---|---|---|---|---|---|
| CESM2 (gn) | Mar | **−4.5** | **7.5** | **0.998** | 2.7 (9.22×10⁶ km²) | 96.4 | 82 360 | 12.92 / 3.50 |
| CESM2 (gn) | Sep | −2.6 | 5.1 | 0.998 | 0.2 (0.57×10⁶ km²) | 97.9 | 84 146 | 7.72 / 1.88 |
| GFDL-ESM4 (gr) | Mar | −9.6 | 21.4 | 0.945 | 6.3 (22.01×10⁶ km²) | 97.2 | 42 122 | 10.04 / 0.99 |
| GFDL-ESM4 (gr) | Sep | −8.0 | 14.3 | 0.987 | 7.5 (26.71×10⁶ km²) | 98.5 | 43 259 | 5.24 / 0.84 |
| MPI-ESM1-2-LR (gn) | Mar | −2.5 | 6.8 | 0.995 | 1.5 (5.38×10⁶ km²) | 96.8 | 32 978 | 11.73 / 0.98 |
| MPI-ESM1-2-LR (gn) | Sep | −2.4 | 14.3 | 0.983 | 3.5 (12.30×10⁶ km²) | 97.6 | 34 894 | 8.76 / 1.10 |

Total ocean areas: CESM2 360.51×10⁶ km², GFDL-ESM4 361.11×10⁶ km², MPI-ESM1-2-LR 361.30×10⁶ km². Retained areas are 347.6–355.9×10⁶ km² in every case, so ~96–99 % of each model's ocean is compared, and the excluded portion is dominated by columns where the operator finds no crossing (`op_missing_any`), roughly 5–13×10⁶ km².

**Ranking (most → least internally consistent)**, criterion = area-weighted RMSD averaged over the March and September panels, with |bias| and r as secondary descriptors:
1. **CESM2** — mean RMSD 6.3 m, |bias| 2.5/4.5 m, r 0.998.
2. **MPI-ESM1-2-LR** — mean RMSD 10.6 m, |bias| ~2.5 m, r 0.983–0.995.
3. **GFDL-ESM4** — mean RMSD 17.9 m, |bias| 8.0–9.6 m, r 0.945–0.987, and the largest >20 m area in both months.

Note that MPI ranks second by RMSD despite the smallest bias: its September panel carries a heavy RMSD penalty from a few large outliers (minimum difference −698 m) while remaining almost unbiased. GFDL-ESM4 is worst on every measure — largest bias, largest RMSD, lowest correlation, largest exceedance area — and is also the only model on `gr`, i.e. the 1°/35-level regridded product used as its own published grid.

## 4. Where the disagreement concentrates and why exact agreement is not expected

The difference fields are strongly non-Gaussian: CESM2 March spans −373 to +166 m with an area-weighted mean of −4.1 m and a cell-mean of similar sign, GFDL-ESM4 March spans −627 to +814 m, MPI-ESM1-2-LR September reaches −698 m. The area-weighted RMSD (7.5, 21.4, 6.8 m) is far larger than |bias|, so the error is concentrated in a small area of very large local mismatch rather than being a uniform offset — the signature of high-latitude deep-convection sites, sea-ice margins and the Southern Ocean, plus shallow shelves where the vertical resolution is coarse and the criterion can be crossed or missed between adjacent levels. Figure 1 is the diagnostic for localising these.

Three differences from the published product explain the systematic shallow bias, and none of them is model error:
1. **Reference level** — the operator displaces a parcel from the model's shallowest level (2.5 m GFDL, 5 m CESM2, 6 m MPI), whereas `mlotst` is computed from its own reference (nominally 10 m). Starting shallower makes the criterion easier to satisfy, giving a shallower depth.
2. **Density route** — the operator uses a displaced-parcel comparison at local pressure with potential temperature and Practical-Salinity-converted-to-Absolute-Salinity input, while `mlotst` in these models is a σ_t (potential-density) diagnostic.
3. **Time averaging** — the operator averages 30 per-profile MLDs, which is systematically deeper than the MLD of the 30-year mean profile would be, but still cannot offset 1 and 2; the published field's own averaging order differs.

In addition, the non-crossing columns (CESM2 often returns a 0 m / no-crossing column, GFDL and MPI fill nearly the whole column) are excluded from the comparison and are reported separately above rather than folded into the metrics.

## 5. Visualizations

- **Figure 1** — a 3×2 panel of true geographic maps of operator MLD minus `mlotst` for the 1971–2000 March and September climatologies, rows CESM2 / GFDL-ESM4 / MPI-ESM1-2-LR and columns March / September. Each native curvilinear or regular field was conservatively regridded to a common 1° longitude–latitude frame for display only, and all six panels share one diverging `RdBu_r` colour scale centred on zero, fixed at ±29.95 m (the 98th percentile of the pooled absolute differences, so that the bulk of the distribution is resolved and a handful of ±600–800 m outliers saturate). Panel titles carry the model, grid label, month and its native-grid bias, RMSD and r.
- **Figure 2** — a 3×2 panel of cell-wise scatter plots, operator MLD on x and `mlotst` on y, one point per retained ocean cell, on identical logarithmic axes over a common range derived from the pooled positive values (ticks at 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000 as permitted by the range), with a dashed y = x line and each panel annotated with its native-grid bias, RMSD and r. This directly visualises the sign of the bias (points displaced above or below y = x) and whether disagreement is a uniform offset or confined to a tail of deep values.

## 6. Important insights

- Internal consistency of the diagnostic is high in the spatial sense for all three models: r ≥ 0.945 and typically ≈ 0.99, so the operator and `mlotst` describe the same structures; the disagreement is in amplitude.
- The bias is always negative and small at low and mid latitudes but locally reaches hundreds of metres, so RMSD is dominated by tail events while the area-weighted bias is dominated by the well-behaved majority — the two must be read together.
- Seasonal contrast matters: CESM2 and MPI agree best in September in bias but CESM2 also has the smallest >20 m area in September (0.2 %); MPI's September RMSD (14.3 m) is twice its March value (6.8 m) despite an unchanged bias, i.e. a few very deep September columns.
- GFDL-ESM4's metrics should be interpreted with its `gr` grid label in mind; part of its larger RMSD may reflect the 1°/35-level product it publishes rather than the native model grid.
- The residual is a diagnostic property: reference level, displaced-parcel vs σ_t, and averaging order of profiles vs profiles of the mean. None of these should be read as a model physics deficiency.
