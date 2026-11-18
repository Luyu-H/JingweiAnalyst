## 1. Direct answer

Three CMIP6 historical simulations (CESM2, GFDL-ESM4, MPI-ESM1-2-LR) were compared with the WOA23 1971–2000 annual 0 m temperature climatology. All three reproduce the broad-scale SST pattern very well (area-weighted spatial correlation with WOA23 of 0.993–0.996 over 40,455 shared ocean cells), but they differ in their globally averaged thermal offset:

- **GFDL-ESM4 is the least biased** (mean bias −0.08 °C, RMSE 0.837 °C, r = 0.996).
- **CESM2 is warm** (+0.29 °C, RMSE 1.063 °C, r = 0.995).
- **MPI-ESM1-2-LR is cold** (−0.63 °C, RMSE 1.342 °C, r = 0.993) and is the most biased of the three.

Ranking from smallest to largest absolute bias: **GFDL-ESM4 (0.08 °C) < CESM2 (0.29 °C) < MPI-ESM1-2-LR (0.63 °C)**.

## 2. Methodology actually executed

**Reference (WOA23).** The file `woa23_decav71A0_t00_01.nc` was opened without time decoding; the objectively analysed field `t_an` (the only `*_an` variable) was selected, the leading `time` and `depth` indices dropped, and the field transposed to (`lat`, `lon`) on the native 1° grid (180 × 360, latitudes −89.5…89.5, longitudes −179.5…179.5). The land–sea mask was read from `landsea_01.msk` as comma-separated text, using the last column (bottom standard level); ocean cells were defined as mask value > 1 and reindexed to the WOA lat/lon coordinates. Units are degrees Celsius as published.

**WOA23 global mean.** Computed over the full 1° WOA ocean domain (mask > 1 AND finite SST), weighted by cell area = (R⊕π/180)² · cos(lat) with R⊕ = 6,371,000 m, i.e. cos(latitude) weights. Unrounded value kept.

**CMIP6 input.** Monthly `tos` from `tos_Omon_{model}_historical_*.nc` was opened with `open_mfdataset(combine="by_coords", use_cftime=True)`, sorted by time, and sliced from the first 1971-01 to the last 2000-12 (360 months). Missing values (1e20) were replaced with NaN via `where(tos < 1e19)` **before** any arithmetic. The climatology is the plain time mean of the monthly `tos` over those 360 months. `areacello` was taken from the single `areacello_Ofx_{model}_historical_*.nc` file, squeezed, and transposed to match the climatology's dimension order.

**Native-grid global means.** Σ(tos·areacello)/Σ(areacello) using each model's own areas, restricted to cells where `tos` and `areacello` are both finite and area > 0. The total ocean area used for each mean is reported below; means were not rescaled to a common area.

**Remapping to the WOA grid.** CESM2 (gx1v7 ≈ 100 km) and MPI-ESM1-2-LR (≈ 250 km) native curvilinear fields were remapped with a conservative, area-weighted scheme implemented by supersampling each source cell with 8 × 8 interior points, assigning each point's share of the source cell area to the target 1° cell that contains it, and forming area-weighted target means. GFDL-ESM4 `tos` is already on its regular 1° `gr` grid and was instead longitude-shifted to [−180, 180), sorted, and bilinearly interpolated to the WOA grid (no conservative remap for that model). Remapped fields were used only for maps, RMSE, bias and correlation — never for the model global means.

**Shared comparison mask.** WOA ocean cells (mask > 1) that are finite in WOA23 SST **and** in all three remapped model fields: **40,455 cells**. This mask was applied to WOA23 and to all three models before bias, RMSE and correlation.

**Metrics on the shared mask.** With cos(latitude) weights on the WOA grid: weighted means, bias = model − WOA23, area-weighted RMSE = √(Σw·Δ²/Σw), and area-weighted Pearson correlation using weighted covariance and variances (cov/√(vx·vy)). Unrounded biases were used for ranking; rounding only at reporting.

## 3. Key numerical results

**Native-grid global means (each over its own ocean grid):**

| Field | Mean SST (°C) | Ocean area (m²) |
|---|---|---|
| CESM2 | **18.52** | 3.605 × 10¹⁴ |
| GFDL-ESM4 | **18.13** | 3.611 × 10¹⁴ |
| MPI-ESM1-2-LR | **17.58** | 3.660 × 10¹⁴ |
| WOA23 (full 1° ocean domain) | **18.23** | 3.495 × 10¹⁴ |

The model domains differ from the WOA domain by a few per cent (≈ 3.5–4.7 % larger), reflecting differing coastline representation and the WOA mask being restricted to cells with valid analysed values.

**Shared-cell statistics (40,455 cells; WOA23 mean on those cells = 18.27 °C):**

| Model | Mean SST (°C) | Bias ΔSST (°C) | RMSE (°C) | Correlation |
|---|---|---|---|---|
| CESM2 | 18.56 | **+0.29** | **1.063** | **0.995** |
| GFDL-ESM4 | 18.19 | **−0.08** | **0.837** | **0.996** |
| MPI-ESM1-2-LR | 17.64 | **−0.63** | **1.342** | **0.993** |

Ranking by |bias|: **GFDL-ESM4 (0.08 °C) → CESM2 (0.29 °C) → MPI-ESM1-2-LR (0.63 °C)**.

## 4. Figures

**Figure 1 (2 × 2 global maps, shared comparison cells).** Three bias panels — CESM2 − WOA23, GFDL-ESM4 − WOA23, MPI-ESM1-2-LR − WOA23 — drawn with `pcolormesh` on the WOA 1° lat/lon grid, all using the same diverging colormap (`RdBu_r`) with a symmetric limit ±vmax, where vmax is the 98th percentile of the pooled absolute bias values (with `extend="both"`), so the three panels are directly comparable and centred on zero. The fourth panel shows WOA23 SST on the same cells with a sequential colormap (`viridis`). Every panel is labelled, axes are labelled Longitude/Latitude, each colour bar carries a °C unit label (SST bias °C; SST °C), and the figure has a common title.

**Figure 2 (bar chart).** The four native-domain global means (CESM2, GFDL-ESM4, MPI-ESM1-2-LR, WOA23) as bars, each annotated with its value to two decimals, y-axis labelled SST (°C), with the axis starting at 0 to make the ≈ 1 °C spread visible.

Both figures were generated successfully.

## 5. Interpretation and insights

- **Large-scale pattern captured; residual biases are spatially structured.** Spatial correlations of ≈ 0.99 indicate that all three models place the main SST gradients correctly. However, the area-weighted RMSE (0.84–1.34 °C) is several times larger than the absolute global-mean bias (0.08–0.63 °C), showing that regional biases of opposite sign largely cancel in the global average. The bias maps in Figure 1 are therefore the informative diagnostic, and should be read as a spatial field rather than a single offset; the shared ±vmax scale puts the three models on equal footing.
- **The models bracket the reference.** CESM2 is warm (+0.29 °C), MPI-ESM1-2-LR cold (−0.63 °C), and GFDL-ESM4 nearly unbiased (−0.08 °C); the spread across models (~0.9 °C between warmest and coldest bias) is an inter-model uncertainty rather than a data-processing artefact.
- **Model ranking is consistent between metrics:** GFDL-ESM4 has both the smallest bias and smallest RMSE and highest correlation, MPI-ESM1-2-LR the largest of all three — so the ranking is not sensitive to which of the three diagnostics is used.
- **Two different "global means" must not be conflated.** The native-grid means (18.52 / 18.13 / 17.58 / 18.23 °C) come from each model's own ocean mask and areas; the shared-cell means (18.56 / 18.19 / 17.64 / 18.27 °C) come from the common 40,455-cell mask. Their differences are not the same quantity as the shared-cell bias, and biases were computed only on the shared mask.
- **Caveats from the set-up.** Each simulation is a single-member coupled historical run, so internal variability is present and no year-to-year or event-level agreement with WOA23 should be expected — only climatological statistics. `tos` is the bulk temperature of the top model cell including under sea ice, and the top layers differ (CESM2 0–10 m, GFDL-ESM4 0–2 m, MPI-ESM1-2-LR 0–12 m), so "SST" is not exactly the same layer across models. Also, coarse grids (especially MPI-ESM1-2-LR, ≈ 250 km) cannot represent narrow features at the amplitude seen in a 1° analysis, and remapping was performed with a conservative area-weighted supersampling scheme for CESM2/MPI-ESM1-2-LR but with bilinear interpolation for the already-1° GFDL-ESM4 field.
