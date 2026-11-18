# Plan

1. Load the CMIP6 historical monthly fields with `load_regular_data`: `thetao`, `so` and `mlotst` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, 1971-01 to 2000-12 (360 monthly records per model). Keep each model on its published grid and record the grid label.
   - Note [cmip6-grid-labels]: GFDL-ESM4 is on `gr` (regridded to 1°, 35 z levels); CESM2 and MPI-ESM1-2-LR are on `gn`; state each model's grid wherever its numbers are quoted.

2. Load each model's geometry with `load_depth_levels` (level values and layer boundaries, needed for the interpolated crossing and for `ref_depth`) and `load_areas` (the model's own published `areacello`, in m²) for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR.

3. Compute the operator MLD per model and per monthly instance with `compute_mixed_layer_depth`, before any time averaging, from that model's own `thetao` and `so`: use the model's shallowest level as `ref_depth`, the Δρ_crit = 0.03 kg m⁻³ criterion (δB_crit = 0.0003 m s⁻²) applied to a displaced surface parcel compared at local pressure, and interpolate the depth between the bracketing levels where the criterion is first crossed. Do not use the KPP boundary-layer depth. Exclude, and count the area of, columns where the threshold is never reached (CESM2 often returns a 0 m column; GFDL/MPI fill roughly the full column), and report that excluded area.
   - Note [cmip6-mlotst-definition]: the operator is the interpolated first-crossing depth for Δρ_crit = 0.03 kg m⁻³ with `ref_depth` = the model's shallowest level; exclude non-crossing columns and report their area.
   - Note [cmip6-thetao-potential-temperature]: `thetao` is potential temperature referenced to the surface; use the potential-temperature route read by `temperature_kind()` in the TEOS-10 density calculation, never the in-situ route.
   - Note [cmip6-so-practical-salinity]: `so` is already Practical Salinity of order 35 with units "0.001" — apply no ×0.001 or ×1000 — and convert S_P to S_A before calling the gsw density routines.

4. For each model, build the March and the September comparison sample: take the 30 March instances of the operator-MLD field and of `mlotst`, and the 30 September instances of each, and call `build_comparison_mask` for each calendar month across its 30 instances so that a cell is retained only where both fields are defined in all 30 instances of that month; apply the mask with `.where(mask)` to both fields, and report the retained cell count and area. Then use `temporal_aggregate` (climatological-monthly, i.e. a plain mean over the 30 years) to form the March and September climatologies of both fields.
   - Note [common-comparison-mask]: one shared mask per model and calendar month, built with `ocean_only=False` (WOA23 is not part of this comparison) and applied with `.where(mask)` before any metric; report retained cells/area beside the metrics.
   - Note [monthly-means-are-time-averages]: compute the MLD of each monthly profile first and average the 30 March (or 30 September) MLDs afterwards — the MLD of an averaged profile is systematically shallower; averaging the 30 monthly means of one calendar month is a plain mean over years with no day weights.

5. For each of the six model-and-month combinations, compute the bias (operator MLD − `mlotst`) and the root-mean-square difference of the two climatologies with `compute_difference_metrics`, area-weighted by the model's own `areacello`.
   - Note [omip-weighted-means-and-integrals]: every map statistic must be area-weighted — `compute_difference_metrics` is area-weighted by default, and `analyze_correlation` needs `weighting='area'`; the >20 m share and the covered-area share are shares of area, not of cells, so use `compute_exceedance_fraction` for the former.

6. Compute the spatial Pearson correlation between the operator-MLD and `mlotst` climatologies for each of the six combinations with `analyze_correlation` and `weighting='area'`.

7. Compute, for each of the six combinations, the area fraction with |operator MLD − `mlotst`| > 20 m using `compute_exceedance_fraction` (area-based, missing cells excluded from the denominator), and the covered-area fraction using `compute_basic_stat` with the model's `areacello`: sum of the retained-mask area divided by the model's total ocean area. Report the percentage and its area.

8. Rank the three models from most to least internally consistent, stating the criterion used (e.g. the area-weighted RMSD averaged over the March and September panels, with |bias| and correlation as secondary descriptors), and report where the disagreement concentrates (for example high-latitude deep-convection and sea-ice-margin regions, the Southern Ocean, and shallow shelves where the vertical resolution is coarse or the threshold is never crossed).

9. Produce Figure 1 by regridding the six difference maps to a common regular longitude–latitude display grid with `regrid_horizontal` (display only), drawing them as true geographic maps in a 3×2 panel (rows CESM2, GFDL-ESM4, MPI-ESM1-2-LR; columns March, September) with one shared diverging colour scale centred on zero. Produce Figure 2 as a 3×2 panel of cell-wise scatter plots (operator MLD on x, `mlotst` on y) on identical logarithmic axes with ticks in every panel, a dashed y = x line, and each panel annotated with its bias, RMSD and correlation from steps 5–6.
   - Note [omip-native-grid-vs-regridded]: compute the bias, RMSD and correlation on each model's own grid with its own `areacello`, and regrid to a regular longitude–latitude grid for the Figure 1 maps only, stating the method ('conservative' when going coarser, 'bilinear'/'nearest_s2d' fallbacks).

Other notes:
- The residual operator-minus-`mlotst` difference is a property of the diagnostic, not of model physics: `mlotst` is a single published field computed with its own reference level (nominally 10 m), its own density route and its own time averaging, whereas the operator uses the model's shallowest level, a displaced-parcel comparison at local pressure and a mean of per-profile MLDs; the operator is therefore expected to be the shallower of the two.
- Report the three differences explicitly (reference level, displaced-parcel vs σ0, time averaging of profiles vs averaging of MLD) instead of interpreting the residual as model error.
- Report the non-crossing/excluded columns and their area for every model and month, separately from the metrics.
- Report each model's grid label beside its numbers, and keep the metrics on the native grid.
