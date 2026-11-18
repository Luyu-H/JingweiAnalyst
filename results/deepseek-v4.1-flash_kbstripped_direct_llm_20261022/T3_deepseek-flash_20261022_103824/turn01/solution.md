The comparison as executed is not usable: the WOA23 land–sea mask reader failed. In the output, `mask_arr` contains zero valid values, so `woa_ocean` is all `False`; consequently `woa_global_mean`, all `comp_stats`/`comp_stats_weighted`, and the common‑cell comparison are `NaN`. Figure 1 is therefore masked out everywhere, and Figure 2 has a `NaN` WOA23 bar. The requested WOA23 global mean, model–WOA comparison means, RMSEs, correlations, ranking, and spatial bias maps cannot be obtained from this run as it stands.

The only reliable numbers printed by the run are the native-grid, area-weighted global means of the three CMIP6 `tos` climatologies over each model’s own ocean grid using its published `areacello` (and `sftof` for CESM2):

- CESM2: **18.52 °C**
- GFDL-ESM4: **18.13 °C**
- MPI-ESM1-2-LR: **17.58 °C**

For diagnosis, if the WOA mask failure is bypassed and the finite WOA `t_an` values are used directly as the ocean mask (41,088 cells), the recovered unweighted statistics on the WOA 1° grid are:

- WOA23 mean on those cells: **13.85 °C**
- CESM2 mean: **14.00 °C**, bias **+0.15 °C**, RMSE **1.054 °C**
- GFDL-ESM4 mean: **13.77 °C**, bias **−0.07 °C**, RMSE **0.828 °C**
- MPI-ESM1-2-LR mean: **13.34 °C**, bias **−0.50 °C**, RMSE **1.325 °C**
- Ranking by absolute mean bias, smallest to largest: **GFDL-ESM4 < CESM2 < MPI-ESM1-2-LR**

The spatial correlation with WOA23 and the area-weighted WOA23 global mean cannot be recovered from the failed output, and the generated bias maps do not show a valid spatial pattern because the common mask is empty. To complete the evaluation, the WOA `landsea_01.msk` file must be parsed correctly (or, more simply, the ocean mask can be taken as `np.isfinite(woa_sst)` since WOA23 `t_an` is already `NaN` over land). Then set `common_mask = woa_ocean & np.isfinite(model_on_woa[model])`, recompute the area-weighted WOA global mean, the common-cell statistics, and regenerate Figure 1 with a non-empty mask.
