# Plan

1. Load CMIP6 historical monthly `thetao`, `so` and `mlotst` for models CESM2, GFDL-ESM4 and MPI-ESM1-2-LR over 1971-01 to 2000-12, and load each model’s `areacello` and vertical level metadata. Tools: `load_regular_data`, `load_volumes_for_cmip6`, `load_depth_levels`.
- Note [missing-values-and-land-cells]: treat exposed fill values as missing immediately and build validity from each loaded field itself, not from an external mask.

2. From each model’s monthly `thetao` and `so`, compute monthly offline operator MLD on the native model grid using the same density-threshold definition as `mlotst`, with threshold-crossing depth linearly interpolated between bracketing levels. Use the model’s shallowest level as `ref_depth`, and exclude cells/months where the threshold is never reached from later operator-versus-`mlotst` comparison. Tools: `compute_potential_density_anomaly`, `compute_mixed_layer_depth`.
- Note [cmip6-mlotst-definition]: compute the offline operator monthly first, then compare only its March/September climatologies against published `mlotst` climatologies, interpreting residuals as diagnostic differences rather than model error.
- Note [cmip6-thetao-potential-temperature]: use `thetao` as potential temperature in the density/MLD calculation.

3. For each model, isolate March and September monthly fields for both operator MLD and published `mlotst`, then enforce the strict shared comparison mask separately for March and September: keep only cells where both fields are defined in all 30 March instances, or all 30 September instances, before any climatology is formed. Then compute the 1971-2000 March and September climatologies as plain means across the 30 monthly values. Tools: `temporal_aggregate`.
- Note [monthly-means-are-time-averages]: compute MLD from each monthly `thetao`/`so` field first and average the 30 same-calendar-month MLDs; do not compute MLD from a month-mean profile.
- Note [common-comparison-mask]: build and apply the all-30-month shared valid-data mask before climatologies and all subsequent metrics.

4. On each model’s native archived grid, compute operator-minus-`mlotst` difference maps for March and September climatologies, then calculate area-weighted mean difference and area-weighted RMSD over the shared mask using each model’s `areacello`. Also compute total valid comparison area and express it as a percentage of total ocean area from `areacello`. Tools: `compute_difference_metrics`, `compute_basic_stat`.
- Note [omip-native-grid-vs-regridded]: keep all summary metrics on the native model grid; any horizontal regridding is only for the display maps.
- Note [omip-weighted-means-and-integrals]: use `areacello` weights on the shared valid cells for map means, RMSD, exceedance area, and valid-area percentage.

5. Compute the area-weighted percentage of the shared comparison area where `|operator MLD - mlotst| > 20 m` for each model and month, using the same native-grid shared mask and `areacello` weights. Tool: `compute_exceedance_fraction`.

6. Compute the spatial Pearson correlation between the March climatology pair and the September climatology pair for each model on the shared native-grid mask. Tool: `analyze_correlation`.

7. Identify where disagreement concentrates by locating the largest positive/negative departures and summarizing the geographical patterns in the native-grid difference fields; use these diagnostics together with the bias, RMSD, correlation, and exceedance-area results to rank the three models from most to least internally consistent, stating the ranking criterion explicitly. Tools: `locate_extreme`, `compute_difference_metrics`.

8. Regrid each March and September climatological difference map from each model’s native grid to a regular longitude-latitude grid for Figure 1 display only, then plot the 3-by-2 panel with one shared diverging colour scale centred on zero. Tool: `regrid_horizontal`.

9. Prepare Figure 2 scatter-plot inputs from the native-grid March and September climatology pairs on the shared mask, one point per ocean cell, and annotate each panel with the already computed bias, RMSD and correlation; use identical logarithmic axes, identical limits/ticks, and a dashed `y=x` line across all six panels. Tools: `analyze_correlation`, `compute_difference_metrics`.

Other notes:
- [cmip6-mlotst-definition] Published `mlotst` and the offline operator are not expected to match exactly because `mlotst` is based on the model top level, a displaced-parcel density test at local pressure, and a monthly mean of instantaneous MLDs, whereas the offline operator uses monthly `thetao`/`so`; the offline result is therefore expected to be systematically shallower in many places.
