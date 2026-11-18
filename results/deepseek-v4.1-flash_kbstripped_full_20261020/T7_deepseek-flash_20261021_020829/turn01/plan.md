# Plan

1. **Observations (WOA23).** Load `woa23` with `load_climatology_data`: `nitrate` and `phosphate`, the annual objectively analysed climatology (all-period composite), all standard levels covering 200–1000 m, global. Also load the WOA23 depth levels and layer boundaries (`load_depth_levels`) and the grid-cell areas (`load_areas`) — the WOA grid has no partial bottom cells, so its own level boundaries give the layer thicknesses.
   - Note [common-comparison-mask]: build one shared comparison mask with `build_comparison_mask` (WOA23 as reference, `ocean_only=True`, the default) and apply it with `.where(mask)` before any cell-by-cell regression, ratio or metric, reporting the mask's cell counts; a mask does not carry between levels, so within 200–1000 m take the valid cells at each depth from the fields themselves, not from a single surface mask.

2. **Models (CMIP6 historical).** Load `cmip6_historical` with `load_regular_data`: `no3` and `po4` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, monthly, 1971-01 to 2000-12. Load each model's `volcello` with `load_volumes_for_cmip6` and its depth levels with `load_depth_levels`. Variables a model does not publish are reported as missing.

3. **Model climatologies.** Reduce each model's monthly 1971–2000 `no3`/`po4` to a 1971–2000 climatology with `temporal_aggregate` (day-weighted long-term mean), giving one annual-mean field per model and variable on that model's own grid.

4. **Layer definition (200–1000 m).** On each of the four sources' own vertical levels, keep the levels whose depths fall in 200–1000 m. Do not vertically interpolate: layer quantities use each product's own levels and thicknesses ([omip-native-grid-vs-regridded] applies).
   - Note [cmip6-static-grid-fields]: the 200–1000 m layer thicknesses and volume weights come from each model's `volcello` / `areacello` (cell-thickness form, which respects partial bottom cells), not from `lev_bnds` differences; use each model's own `areacello`/`volcello` and report `attrs['area_source']` and `attrs['volume_source']`, never rescaling one model's area to another's.

5. **Cell-by-cell N:P.** Compute N:P at every cell of the 200–1000 m layer for all four sources with `compute_ratio` (nitrate / phosphate), then form the depth-mean over the layer at each horizontal cell using the layer thicknesses/`volcello` as vertical weights (step 4), and build the depth-integrated volume weights for the layer.
   - Note [omip-weighted-means-and-integrals]: all layer means are volume-weighted, Σ x·volcello / Σ volcello over valid cells (`op='volume_weighted_mean'`); the mean of the cell-by-cell N:P and the ratio of the two volume-weighted means are different quantities and both are reported, clearly labelled; volume percentages come from `compute_exceedance_fraction`, not from a share of cells, with missing cells left out of the denominator.

6. **Headline stoichiometry, global.** With `compute_basic_stat` (`volume_weighted_mean` over the whole 200–1000 m layer, global domain) report for each of the four sources: (a) the mean of the cell-by-cell N:P ratio (headline), and (b) the ratio of the volume-weighted mean nitrate to the volume-weighted mean phosphate. Both to 2 decimal places.

7. **Concentrations for the regression.** Express the model nitrate and phosphate of the 200–1000 m layer in micromoles per kilogram at the constant density of 1025 kg/m³ prescribed by the question (µmol kg⁻¹ = mol m⁻³ × 10⁶ / 1025); WOA23 nutrient fields are micromoles per litre, taken as micromoles per kilogram.
   - Note [cmip6-bgc-units-per-volume]: `no3`/`po4` are archived per volume (mol m⁻³), so convert with the ×10⁶/ρ factor and never relabel a ×1000 prefix change as "µmol/kg"; N:P is invariant under a common density factor, so the mean-of-ratio and ratio-of-means values need no density — only the regression concentrations do.

8. **Regression.** For each of the four sources, stack all valid cells of the 200–1000 m layer (level × latitude × longitude, cells weighted equally) and fit NO3 = a·PO4 + b with `fit_multivariate_relation`; report slope a (2 dp), intercept b in µmol kg⁻¹ (2 dp) and R² (3 dp), with the diagnostics the tool returns.

9. **Volume beyond thresholds.** With `compute_exceedance_fraction` (volume-based, NaN-aware, returning the m³ behind each %) compute for each source: the percentage of water volume with N:P strictly between 14 and 18 (global), and the percentage below 12 (global, and within the Pacific and the Indian Ocean). 1 decimal place.

10. **Basin means.** Load the WOA23 `basinmask_01` with `load_basin_mask` and, with `select_basin`, keep the Pacific, Atlantic and Indian Oceans in turn, blanking the rest, for the 200–1000 m N:P fields of all four sources; then repeat the volume-weighted headline mean of step 6 within each basin (2 dp).
    - Note [basin-and-region-definitions]: use the WOA23 `basinmask_01` with `include_marginal_seas=True` so marginal seas fold into their parent ocean, and `depth_dependent=False` (surface labels, what a named ocean usually means); state that the OMIP 10-region Ofx `basin` field is absent locally so a WOA labelling was used, and that the Southern and Arctic Oceans are excluded here as the question defines.

11. **Model-versus-observation diagnosis.** Compare each model's headline N:P, basin means and low-N:P volume fractions with WOA23, and state for each source whether its stoichiometry sits at, above or below Redfield (16:1), and whether each model reproduces the observed low-N:P volume and its basin distribution. Report the shared-mask cell counts from step 1 here for the like-for-like assessment [common-comparison-mask].

12. **Figure 1.** 4-panel scatter of nitrate (y) against phosphate (x) in µmol kg⁻¹ over 200–1000 m, panels ordered WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, identical axis ranges and ticks, each panel with its own fitted line from step 8 and the Redfield 16:1 line through the origin, visually distinguishable, annotated with a, b and R².
   - Note [omip-native-grid-vs-regridded]: every number in steps 5–11 is computed on each source's own grid with its own areas and volumes; regrid only for the Figure 2 maps (display only), stating the method — 'conservative' onto the coarser regular WOA grid for the fields, 'nearest_s2d' for masks — and do not vertically interpolate for the 200–1000 m layer.

13. **Figure 2.** 2×2 panel of global maps of the 200–1000 m mean N:P minus 16 for the four sources, each on a common regular longitude–latitude frame (curvilinear model fields put on that grid with `regrid_horizontal`, conservative, display only), on one shared diverging colour scale centred on zero.

Other notes:
- [cmip6-bgc-units-per-volume] The question prescribes a constant ρ = 1025 kg/m³, whereas the protocol asks for in-situ density (`compute_potential_density_anomaly`, `ref_pressure=None`); the constant is followed and the ~1% error it introduces in the µmol kg⁻¹ concentrations, and hence in the regression slope/intercept, is stated.
- WOA23 nutrients are micromoles per litre and are treated as micromoles per kilogram (difference <2% for ρ ≈ 1024–1028 kg/m³); only the model fields need the density factor.
- [cmip6-grid-labels] GFDL-ESM4 is `gr` (already remapped to 180×360, conserve_order1), so its "native" numbers sit on the 1° archive grid with a regridded land–sea mask and only approximately conservative volumes; CESM2 (~100 km) and MPI-ESM1-2-LR (~250 km) differ in resolution, so narrow features are not comparable in amplitude.
- [omip-native-grid-vs-regridded] is followed for the numbers: the headline per-source statistics, regressions and exceedance fractions use each source's own grid, own cells and own volumes, and the shared comparison mask of step 1 is reported for cell counts and for the model-vs-WOA diagnosis rather than replacing the per-source numbers.
- The Southern and Arctic Oceans are excluded from the basin results; the three reported basins are Pacific, Atlantic and Indian.
