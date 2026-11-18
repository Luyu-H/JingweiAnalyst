# Plan

1. Open the raw files for WOA23 and CMIP6 historical. From WOA23 take `temperature` for the 1971–2000 climate normal, annual time slice, 0 m depth, on the 1-degree grid, and take `landsea_01`; define WOA ocean cells as `landsea_01 > 1`. From CMIP6 historical take monthly `tos` for 1971-01 to 2000-12 and `areacello` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR. Decode CF metadata and explicitly replace the missing value 1e20 with NaN before any arithmetic.
- Note [missing-values-and-land-cells]: CMIP6 land and below-sea-floor cells use 1e20; mask these out before sums/means so they cannot dominate the statistics.

2. For each CMIP6 model, compute the 1971–2000 climatological SST map on its native grid by averaging the monthly `tos` values over the 360 months from 1971-01 to 2000-12. GFDL-ESM4 `tos` is already on its `gr` 180×360 1-degree grid, so its native map/mask are the regridded ones; CESM2 and MPI-ESM1-2-LR are on their native `gn` grids.
- Note [cmip6-grid-labels]: GFDL-ESM4 `tos` is `gr` already remapped to 180×360 1-degree; CESM2 is gx1v7 (~100 km) and MPI-ESM1-2-LR is ~250 km, so the coarse grid cannot represent narrow features at the same amplitude.

3. Compute each model’s native-grid global mean SST using its own `areacello`, restricted to cells where `tos` is valid: Σ(`tos`×`areacello`)/Σ`areacello`, with weights masked to the same valid cells. Also compute the WOA23 global mean over the full 1-degree WOA ocean domain (`landsea_01 > 1`), not restricted to the later comparison cells, using cos(latitude) cell-area weights. Report each mean in °C to 2 decimal places and state the ocean area each mean was taken over. Keep unrounded values for later ranking.
- Note [omip-weighted-means-and-integrals]: global mean SST is area-weighted over the surface field, Σ`tos`·`areacello`/Σ`areacello`, with weights masked to valid cells; a plain mean over latitude/longitude is wrong.
- Note [cmip6-static-grid-fields]: use each model’s own `areacello`; total ocean area differs between models by about a per cent because each resolves the coastline differently, so state the area each mean was taken over rather than rescaling one model to another.
- Note [cmip6-tos-bulk-top-cell]: `tos` is the bulk temperature of the top model cell, including under sea ice; the top cells differ (CESM2 0–10 m, GFDL-ESM4 0–2 m, MPI-ESM1-2-LR 0–12 m), so “SST” is not exactly the same layer across models.

4. Remap each model’s 1971–2000 SST climatology to the WOA23 1-degree grid using conservative area-weighted remapping. Use the native fields for the model global means in step 3 and the remapped fields only for the comparison metrics and maps. State that conservative remapping was used.
- Note [omip-native-grid-vs-regridded]: compute each global mean on the model’s own grid with its own `areacello`; regrid only for the maps, RMSE and pattern correlation on the common mask, using conservative remapping.

5. On the WOA23 1-degree grid, build one shared comparison mask: WOA ocean cells (`landsea_01 > 1`) that are also finite in WOA23 SST and in all three remapped model SST fields. Apply `.where(mask)` to WOA23 and to all three model fields before computing bias, RMSE and correlation. Report the number of cells in this mask.
- Note [common-comparison-mask]: put all three models and WOA23 on the WOA23 1-degree grid, use one mask with `ocean_only=True`, and apply it to every field; report the shared-sample `mean_data`/`mean_reference`, and state that each model’s own-grid global mean is a different quantity whose difference from another global mean is not the shared-cell bias.

6. On the shared mask, using WOA-grid cos(latitude) cell-area weights, compute for each model: the weighted mean model SST, the weighted mean WOA23 SST over the same cells, the mean bias delta_SST = model minus WOA23, the area-weighted RMSE, and the area-weighted spatial Pearson correlation with WOA23. Use the unrounded biases for ranking and round only when reporting: means and biases to 2 decimal places, RMSE and correlation to 3 decimal places. Rank the models from smallest to largest absolute delta_SST.

7. Produce Figure 1 as a 2-by-2 panel of global maps on the shared comparison cells: CESM2 minus WOA23, GFDL-ESM4 minus WOA23, MPI-ESM1-2-LR minus WOA23, and WOA23 SST. Use one shared diverging colour scale centred on zero for the three bias panels, with a symmetric limit set by the largest absolute bias among them; use a sequential colour scale for the WOA23 reference panel. Label every panel and every colour bar in °C.

8. Produce Figure 2 as a bar chart of the four native-domain global mean SST values from steps 3: CESM2, GFDL-ESM4, MPI-ESM1-2-LR and WOA23, each over its own native ocean grid. Annotate each bar with its value and label the y axis in °C.

Other notes:
- [cmip6-experiments-historical-ssp245] These are coupled single-member historical runs with internal variability; a 1971–2000 climatology is a legitimate statistic, but no match of individual years or event timing with WOA23 should be expected.
