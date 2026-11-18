## 1. Direct answer

On the shared 1° ocean sample, GFDL-ESM4 has the smallest absolute mean SST bias against WOA23 (−0.09 °C). CESM2 is next (+0.28 °C) and MPI-ESM1-2-LR has the largest (−0.63 °C). The ranking from smallest to largest |ΔSST| is **GFDL-ESM4 < CESM2 < MPI-ESM1-2-LR**.

GFDL-ESM4 also has the lowest RMSE and the highest spatial correlation. All three models reproduce the large-scale SST pattern very well (r ≥ 0.993). Most of the error sits in regional features. The largest departures are in the North Atlantic (CESM2, MPI-ESM1-2-LR), the south-east Atlantic off Angola/Namibia (GFDL-ESM4, MPI-ESM1-2-LR), and the Persian Gulf (GFDL-ESM4).

## 2. What the code did

- **Data.**
  - `tos` was loaded for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR for 1971-01 to 2000-12, each on its own grid.
  - `temporal_aggregate(freq="mean", op="mean")` gave one 2-D long-term mean per model. The plan calls this the day-weighted mean. I did not check the weighting separately.
  - WOA23 `temperature` was loaded as the 1971-2000 annual field at 0 m (`depth_range=[0,1]`, first level selected) on the 1° grid.
- **Fill values.** The maximum absolute value of each climatology was checked and asserted to be below 1e10. It was about 31 °C for every model, so no 1e20 fill values entered any statistic. The same check was applied to WOA23. Sea-ice cells were not removed.
- **Native-grid means.**
  - `compute_basic_stat(op="area_weighted_mean", dim=["lat","lon"])` was applied to each model climatology and to WOA23.
  - The code did not call `load_areas` itself. The area weights were resolved inside the tool.
  - The code recorded "areacello" for every model. That string is also the fallback in the code if the tool returned no weight attribute, so it does not independently confirm the weight source.
  - The ~1 % difference in total ocean area between models was not computed in this run, and the areas were not rescaled.
- **Regridding.** `regrid_horizontal(..., target="woa23", method="conservative")` was used for all three models. GFDL-ESM4 is already on a 1° `gr` grid, and its regridded statistics are identical to its native ones (same 44,564 valid cells).
- **Common mask.**
  - `build_comparison_mask` was called once on the three regridded models plus WOA23, with `ocean_only=True`.
  - The WOA23 `landsea_01` test (value > 1) removed 0 cells.
  - The mask was applied with `.where(mask)` to all fields before any metric.
- **Scores.**
  - `compute_difference_metrics` (area-weighted by default) gave mean model, mean WOA23, bias, RMSE and the location of the largest departure.
  - `analyze_correlation(weighting="area")` gave the spatial Pearson correlation.
  - `locate_extreme` gave the largest positive and negative differences.
  - Ranking used the unrounded absolute bias.

## 3. Key numerical results

**Native-grid global means (area-weighted, each on its own ocean grid; Figure 2 values)**

| Product | Global mean SST (°C) |
|---|---|
| CESM2 (`gn`) | 18.52 |
| GFDL-ESM4 (`gr`, 1° archive grid) | 18.13 |
| MPI-ESM1-2-LR (`gn`) | 17.58 |
| WOA23 (whole 1° grid, valid cells) | 18.19 |

GFDL-ESM4's "native" grid is a 1° regridded grid, not a tripolar grid.

The WOA23 field's plain unweighted mean over valid cells is about 13.85 °C. That figure is not the area-weighted global mean and is not used in the results.

**Shared 1° comparison sample**

- The mask has **40,830 cells**, out of 41,088 valid WOA23 cells (the smallest input).
- Valid cells per field before masking were 42,576 (CESM2), 44,564 (GFDL-ESM4), 42,720 (MPI-ESM1-2-LR) and 41,088 (WOA23).
- The land-sea test removed none of the 258 dropped cells. They were dropped because at least one field lacked data there.
- These means come from `compute_difference_metrics` on the shared sample. They differ from the native-grid means above and are not interchangeable with them.

| Model | Mean model SST (°C) | Mean WOA23 SST (°C) | Bias (°C) | RMSE (°C) | Spatial r |
|---|---|---|---|---|---|
| CESM2 | 18.48 | 18.20 | +0.28 | 1.062 | 0.995 |
| GFDL-ESM4 | 18.12 | 18.20 | −0.09 | 0.832 | 0.997 |
| MPI-ESM1-2-LR | 17.57 | 18.20 | −0.63 | 1.339 | 0.993 |

**Ranking by |bias|, smallest to largest:** GFDL-ESM4 (0.086), CESM2 (0.280), MPI-ESM1-2-LR (0.631).

## 4. Spatial structure of the bias (extremes from `locate_extreme`)

Coordinates are cell centres, with longitude in −180 to 180.

| Model | Largest positive bias | Largest negative bias |
|---|---|---|
| CESM2 | +8.84 °C at 43.5°N, 64.5°W (north-west Atlantic, Nova Scotia/Gulf of Maine shelf region) | −4.80 °C at 47.5°N, 39.5°W (North Atlantic) |
| GFDL-ESM4 | +7.38 °C at 16.5°S, 11.5°E (south-east Atlantic, Angola/Namibia coast) | −4.49 °C at 24.5°N, 52.5°E (Persian Gulf) |
| MPI-ESM1-2-LR | +7.61 °C at 16.5°S, 11.5°E (south-east Atlantic) | −7.33 °C at 44.5°N, 41.5°W (North Atlantic) |

- **CESM2:** a warm bias along the north-west Atlantic shelf paired with a cold bias east of it in the North Atlantic. This pattern is consistent with a displaced Gulf Stream/North Atlantic Current. The mean bias is positive.
- **MPI-ESM1-2-LR:** a deep cold bias in the central North Atlantic, the largest in magnitude of any model. It also has the same south-east Atlantic warm bias as GFDL-ESM4. Its negative mean bias (−0.63 °C) is the largest of the three.
- **GFDL-ESM4:** the near-zero mean bias comes from compensation between a warm south-east Atlantic bias and a cold bias in the Persian Gulf. Its RMSE is the lowest because its errors are the smallest and least widespread.
- **South-east Atlantic:** the warm bias at 16.5°S, 11.5°E is the largest positive error in two of the three models. This is the Angola/Benguela margin, where coarse models tend to be too warm. That reading is an interpretation, not a diagnostic from the code.

The code reported only the extremes, not region-by-region statistics. Broader patterns would have to be read from the maps.

## 5. Figures

- **Figure 1** is a 2×2 panel of global maps on the common mask.
  - Three panels show the bias (model minus WOA23) for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR. They use the `RdBu_r` colormap with one shared colour scale, symmetric about zero. The limit is ±(98th percentile of the pooled |bias|), with out-of-range values clipped (`extend="both"`). The titles include each model's mean bias.
  - The fourth panel shows the WOA23 0 m annual SST on the same masked cells with the sequential `viridis` colormap.
  - The axes are labelled, and the colourbars are labelled "SST bias (°C)" and "SST (°C)".
- **Figure 2** is a bar chart of the four native-grid global means (CESM2, GFDL-ESM4, MPI-ESM1-2-LR, WOA23). Each bar is annotated with its value to 2 decimals, and the y axis is labelled "Global mean SST (°C)".

## 6. Caveats

- **Top-cell thickness.** `tos` is a top-cell bulk temperature, with the top cell about 10 m thick in CESM2, 2 m in GFDL-ESM4 and 12 m in MPI-ESM1-2-LR. It is compared like-for-like with WOA23 at 0 m. This depth mismatch is not corrected and adds some uncertainty to the small mean biases.
- **Sea ice.** Under-ice SST of about −1.8 °C is retained in all means.
- **Regridding.** RMSE and correlation depend on the conservative regridding to 1°.
- **Native versus shared means.** Native-grid means and shared-mask means differ and must not be mixed. For example, CESM2 minus WOA23 is +0.33 °C in the native means (Figure 2) but +0.28 °C on the shared sample.
- **Rounding.** Values were rounded only for reporting.
