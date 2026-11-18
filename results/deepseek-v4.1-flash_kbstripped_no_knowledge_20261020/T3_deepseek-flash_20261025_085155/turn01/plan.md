# Plan

1. Load CMIP6 historical `tos` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, monthly 1971-01 to 2000-12, on each model’s native grid using `load_regular_data`; convert from K to degrees Celsius with `transform_unit`. Load each model’s native `areacello` cell areas with `load_areas`. For each model, form the 1971-2000 SST climatology with `temporal_aggregate` using the day-weighted long-term mean.
2. For each model, compute the native-ocean global mean of its 1971-2000 SST climatology with `compute_basic_stat`, using `areacello` as spatial weights and restricting to its native valid ocean cells.
3. Load WOA23 `temperature`, 1971-2000 climate normal, annual time slice, 0 m depth, with `load_climatology_data`; load the WOA23 1-degree cell areas with `load_areas`. Compute the WOA23 global mean over the whole 1-degree grid with `compute_basic_stat`, area-weighted by the WOA23 1-degree cell areas.
4. Load the WOA23 1-degree `landsea_01` mask with `load_landsea_mask`; treat cells with value greater than 1 as ocean, as specified.
5. Regrid each model’s 1971-2000 SST climatology horizontally to the WOA23 1-degree grid with `regrid_horizontal`. Use `build_comparison_mask` on the regridded CESM2, GFDL-ESM4 and MPI-ESM1-2-LR SST fields together with WOA23 SST, restricted to WOA ocean cells (`landsea_01 > 1`) and cells valid for all fields.
6. On the common comparison mask, compute the area-weighted mean SST of each regridded model and of WOA23 with `compute_basic_stat`, using WOA23 1-degree cell areas as weights. For each model, use `compute_difference_metrics` with model as the first field and WOA23 as the reference to obtain mean bias `delta_SST` (model minus WOA23) and RMSE. Use `analyze_correlation` to obtain the spatial Pearson correlation between each model SST and WOA23 over the same cells.
7. Form each model-minus-WOA23 bias field on the common comparison cells for mapping. Rank CESM2, GFDL-ESM4 and MPI-ESM1-2-LR from smallest to largest absolute mean `delta_SST`. Use the bias maps to describe broad spatial structure, such as coherent regional warm or cold bias patterns.
8. Produce Figure 1 as a 2-by-2 panel: the CESM2, GFDL-ESM4 and MPI-ESM1-2-LR bias fields on the comparison cells, all using one shared diverging colour scale centred on zero, plus WOA23 0 m SST on the same cells using a sequential colour scale; label every panel and colourbar in degrees Celsius.
9. Produce Figure 2 as a bar chart of the four global mean SST values from steps 2 and 3, each over its own native grid, annotated with its value and with the y axis labelled in degrees Celsius.

Other notes:
- CMIP6 `tos` is converted from K to degrees Celsius before all climatology, bias and mean calculations; WOA23 temperature is already in degrees Celsius.
- Native-grid model global means use each model’s own `areacello` and native ocean cells; the WOA23 global mean uses WOA23 1-degree areas over the whole grid, not the land-sea comparison mask.
- The model-WOA23 comparison is performed on the WOA23 1-degree grid after horizontal regridding, using only WOA ocean cells (`landsea_01 > 1`) that are valid for all three models and WOA23.
- Bias is defined as model minus WOA23; the model ranking uses absolute mean bias.
- RMSE and spatial Pearson correlation are computed over the same common comparison cells as the comparison means and biases.
- Figure 1 uses a shared diverging scale centred at zero for the three model biases and a separate sequential scale for the WOA23 reference field; all colour scales are labelled in degrees Celsius.
