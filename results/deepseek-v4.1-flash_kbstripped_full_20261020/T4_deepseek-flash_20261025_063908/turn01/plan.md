# Plan

1. Load the model data with `load_regular_data` from `cmip6_historical`: variables `thetao`, `so`, `mlotst` and the static `areacello`, for models CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, monthly, 1971-01 to 2000-12 (all three models' historical runs cover this window; a variable a model does not publish comes back as missing and must be reported as such).

2. Load each model's vertical coordinate with `load_depth_levels` (level centres and layer boundaries) for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR.
   - Note [cmip6-mlotst-definition]: compute the MLD as a displaced surface parcel compared with the local parcel at each level's in-situ pressure (δB = −g(ρ_displaced − ρ_local)/ρ_local), taking the first level where δB > 0.0003 m s⁻² (Δρ_crit = 0.03 kg m⁻³), linearly interpolated between the bracketing levels, with `ref_depth` set to the model's shallowest level rather than 10 m; compute the MLD of each monthly profile and average the MLDs afterwards (never average `thetao`/`so` first); exclude columns that never reach the threshold, report the area excluded, and state the three definitional differences (reference level, σ0-difference versus displaced-parcel comparison, averaging order) so the residual against `mlotst` is not read as model error; expect the offline MLD to be shallower than `mlotst`.

3. Compute density from temperature and salinity with `compute_potential_density_anomaly` (TEOS-10), per model, per month, per level.
   - Note [cmip6-thetao-potential-temperature]: `thetao` is potential temperature referenced to the surface, so the calculation must follow the potential-temperature route (`temperature_kind()` reads this from the metadata); the in-situ route biases density increasingly with depth, and no archive carries `bigthetao`.
   - Note [cmip6-so-practical-salinity]: `so` stores practical-salinity numbers (about 35) despite the "0.001" units — do not multiply by 0.001 or 1000; convert S_P to S_A with the location-dependent function only at the TEOS-10 call that expects S_A.

4. Compute the operator MLD with `compute_mixed_layer_depth` from the density field, using a 0.03 kg m⁻³ threshold, `ref_depth` = each model's shallowest level, and the interpolated crossing described in step 2, for every monthly profile of every model.

5. Build the March and September climatologies of both fields with `temporal_aggregate` (climatological-monthly), then select March and September, separately for the operator MLD and for the published `mlotst`.
   - Note [monthly-means-are-time-averages]: the diagnostic is applied first and averaged second — f(monthly mean) ≠ mean of f, so the MLD of an averaged profile is never used; the March (September) climatology is the plain mean of that calendar month's 30 values, with no day weighting for a single-calendar-month climatology. State this order in the answer.

6. Build, with `build_comparison_mask`, one mask per model and calendar month over the 30 monthly instances of that month, keeping only cells where both the operator MLD and `mlotst` are defined in all of the 30 Marches (or all 30 Septembers), and apply that single mask to both climatologies before any statistic.
   - Note [common-comparison-mask]: a bias, RMSE or pattern correlation exists only where both fields have data, so put both fields on the model's own grid, build one mask over both and apply it to both; WOA23's land-sea mask stays out of this comparison because WOA23 is not involved.

7. Score each of the six model–month combinations with `compute_difference_metrics` on the masked fields: mean difference (operator MLD minus `mlotst`), root-mean-square difference, and the single largest departure with its location, weighted in space by the model's own cell areas.
   - Note [cmip6-static-grid-fields]: use each model's published `areacello` as the area weight and state the total ocean area each percentage is taken over (total ocean area differs between models by about a per cent) rather than rescaling one model to another.

8. Compute the spatial Pearson correlation between the two masked climatologies for each of the six combinations with `analyze_correlation`.

9. Compute, with `compute_exceedance_fraction` applied to |operator MLD − `mlotst`| with a 20 m threshold and area weighting by `areacello`, the percentage of the common area where the absolute difference exceeds 20 m; and compute with `compute_basic_stat` (area-weighted, using `areacello`) the area of the common mask as a percentage of each model's total ocean area (the coverage percentage).

10. Rank the three models from most to least internally consistent for March and for September, using the area-weighted RMSE against `mlotst` as the primary criterion with the spatial correlation and the mean difference reported alongside; report both calendar months and state the criterion explicitly in the answer.

11. Figure 1: with `regrid_horizontal` put each of the six difference fields on one common regular longitude–latitude frame for display, drawn as six true geographic maps (rows CESM2, GFDL-ESM4, MPI-ESM1-2-LR; columns March, September) with one shared diverging colour scale centred on zero.
   - Note [omip-native-grid-vs-regridded]: regrid only for the maps and the cell-by-cell comparison, using 'conservative' when going to a coarser grid and stating the method; take the means, the area percentages and the correlations on each model's own grid with its own cell areas, and do not let the display grid's areas enter the area-weighted percentages.

12. Figure 2: scatter plots for the six combinations on the native masked grids, operator MLD (x) against `mlotst` (y), one point per ocean cell, logarithmic axes with identical range and identical ticks in every panel, a dashed y = x reference line, and each panel annotated with its bias, RMSE and correlation from steps 7 and 8.
   - Note [cmip6-grid-labels]: GFDL-ESM4 is on the `gr` grid, already remapped from 75 hybrid layers to 35 z levels (conserve_order1), an extra resolution difference; CESM2 (~100 km) and MPI-ESM1-2-LR (~250 km) differ by about 2.5×, so features narrower than the coarse grid are not comparable in amplitude.

Other notes:
- Both fields come from the same simulation, so every disagreement is a property of the diagnostic (definition and discretisation) and of the vertical grid, not of the model's physics; the answer must say so explicitly.
- The disagreement is expected to concentrate where the mixed layer is deep and the vertical spacing near its base is coarse — the Southern Ocean, the North Atlantic and Nordic Seas deep-convection regions, and the Arctic — which the Figure 1 maps must be read for.
- The ranking is a ranking of the internal consistency of the diagnostic, not of model quality.
- Report values rounded only at the point of reporting: mean and RMS differences and the two percentages to 1 decimal place, correlations to 3 decimal places; carry unrounded intermediate fields through all steps.
