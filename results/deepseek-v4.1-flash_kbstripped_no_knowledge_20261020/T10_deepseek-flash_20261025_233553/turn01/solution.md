## 1. Direct answer

The CMIP6 models reproduce the observed oxygen-minimum zone only poorly and disagree about its future. Using the WOA23 all-period oxygen climatology, water with O₂ < 20 µmol kg⁻¹ occupies 14.2 × 10¹⁵ m³ and water with O₂ < 60 µmol kg⁻¹ occupies 80.3 × 10¹⁵ m³. In the 1971–2000 historical means both models overestimate the OMZ: GFDL‑ESM4 gives 43.1 and 105 × 10¹⁵ m³, MPI‑ESM1‑2‑LR gives 92.0 and 212 × 10¹⁵ m³. GFDL‑ESM4 is clearly the better of the two at both thresholds (about 1.3× the observed <60 volume versus about 2.6× for MPI). Under ssp245 (2071–2100) the models do not agree: GFDL‑ESM4 expands at 60 µmol kg⁻¹ (+4.41 %) while contracting at 20 µmol kg⁻¹ (−4.06 %); MPI‑ESM1‑2‑LR contracts at both (−3.03 % at 20, −0.63 % at 60). CESM2, which has no historical O₂, has an ssp245 <20 volume of 60.9 and <60 volume of 165 × 10¹⁵ m³. The WOA23 OMZ at 60 µmol kg⁻¹ is strongly respiration‑dominated: volume‑weighted AOU = 271.3 µmol kg⁻¹ and saturation = 11.3 %.

## 2. Methodology as actually executed

**Observations (WOA23).** `load_climatology_data(dataset="woa23", variable=...)` was called for `oxygen`, `percent_oxygen_saturation` and `apparent_oxygen_utilization` with no time-selection argument, so the all-period composite on the native 102-level, 1° grid (180 lat × 360 lon) was used; no separate 1971–2000 oxygen normal was requested. Depth levels and bounds came from `load_depth_levels(dataset="woa23", variable="oxygen")`, and layer thickness `dz` was formed as the difference between the lower and upper bound of each layer. Areas came from `load_areas(dataset="woa23")` and WOA23 cell volumes were built as `area × dz`. `load_landsea_mask()` was called and the mask object was retained, but the volume product itself was formed from the area returned by `load_areas` alone; no separate additional multiplication by the land–sea mask was performed in the executed volume calculation.

**CMIP6.** `load_regular_data` was used with `dataset="cmip6_historical"` (GFDL‑ESM4, MPI‑ESM1‑2‑LR) over `time_range=["1971-01", "2000-12"]` and `dataset="cmip6_ssp245"` (CESM2, GFDL‑ESM4, MPI‑ESM1‑2‑LR) over `["2071-01", "2100-12"]`, for `o2`, `thetao` and `so`. Each monthly series was reduced with `temporal_aggregate(freq="mean", op="mean")`. Potential density anomaly was computed with `compute_potential_density_anomaly(temp_data=thm, sal_data=sam, ref_pressure=None)`, and `transform_unit` converted the mean `o2` from mol m⁻³ to µmol kg⁻¹ using that density anomaly as the `density` argument. The unconverted mol m⁻³ mean was kept for the inventory. Model volumes were taken from `load_volumes_for_cmip6`, areas from `load_areas`, and level bounds from `load_depth_levels(..., variable="o2")`.

**Volume diagnostics.** OMZ volumes were obtained from `compute_exceedance_fraction(field, threshold=T, comparison="less", dim=None, weighting="volume")` for T = 20 and 60 µmol kg⁻¹, using each product’s own cell volumes (WOA23: area × layer thickness; CMIP6: `volcello`), reading `exceeding_amount` in m³ and dividing by 10¹⁵. Thresholds are strict.

**Thickness.** `identify_oxygen_minimum_zone(field, threshold=60.0)[0]` produced the column-integrated thickness field at 60 µmol kg⁻¹, and `locate_extreme(..., op="max")` returned the maximum and its coordinates.

**Change.** Percent change = 100 × (ssp245 volume − historical volume)/historical volume, computed on unrounded m³.

**Inventory.** `convert_concentration_to_content(o2_mol_mean, density=1025.0)` was applied to the original mol m⁻³ mean and summed globally with `compute_basic_stat(op="sum", dim=None)`.

**Basins.** `select_basin(field, basin=B, include_marginal_seas=True, depth_dependent=False)` was applied directly to the WOA23, GFDL‑ESM4 and MPI‑ESM1‑2‑LR historical fields (the WOA23 basin mask was not regridded in the executed code), followed by the same volume-weighted `compute_exceedance_fraction` at 60 µmol kg⁻¹. Southern and Arctic were not requested.

**AOU/saturation.** `compute_basic_stat(..., op="volume_weighted_mean", dim=None)` on `aou_woa` and `po2_woa` where `o2_woa < 60`.

## 3. Key numbers

OMZ volume (10¹⁵ m³, strict thresholds):

| Source | < 20 | < 60 | max thickness at <60 |
|---|---|---|---|
| WOA23 | 14.2 | 80.3 | 3125.0 m at 7.50°N, 121.50°E |
| GFDL‑ESM4 historical | 43.1 | 105 | 4343.2 m at 13.50°N, −93.50° |
| MPI‑ESM1‑2‑LR historical | 92.0 | 212 | 4960.0 m at −9.39°, −81.22° |
| CESM2 ssp245 | 60.9 | 165 | 3631.9 m at 7.45°N, 121.44°E |
| GFDL‑ESM4 ssp245 | 41.3 | 110 | 4318.2 m at 13.50°N, −93.50° |
| MPI‑ESM1‑2‑LR ssp245 | 89.2 | 211 | 4950.0 m at −9.39°, −81.22° |

Percentage change, ssp245 − historical (positive = expansion): GFDL‑ESM4 −4.06 % (<20) and +4.41 % (<60); MPI‑ESM1‑2‑LR −3.03 % (<20) and −0.63 % (<60).

Total dissolved oxygen inventory (mol): GFDL‑ESM4 historical 2.20 × 10¹⁷, MPI‑ESM1‑2‑LR historical 2.23 × 10¹⁷, CESM2 ssp245 2.02 × 10¹⁷, GFDL‑ESM4 ssp245 2.15 × 10¹⁷, MPI‑ESM1‑2‑LR ssp245 2.17 × 10¹⁷.

Basin OMZ volume at <60 µmol kg⁻¹ (10¹⁵ m³), Pacific / Indian / Atlantic: WOA23 67.4 / 11.6 / 1.29 (basin totals 642 / 235 / 300); GFDL‑ESM4 historical 92.8 / 5.43 / 6.39 (627 / 227 / 289); MPI‑ESM1‑2‑LR historical 180 / 19.2 / 11.3 (627 / 226 / 289).

WOA23 inside the <60 OMZ: volume-weighted AOU 271.3 µmol kg⁻¹, saturation 11.3 %.

## 4. Interpretation

- **Model skill.** At both thresholds the observed volume lies far below both models, but GFDL‑ESM4 is much closer (105 vs 80.3 × 10¹⁵ m³ at 60 µmol kg⁻¹, ratio ≈ 1.31; 43.1 vs 14.2 at 20, ratio ≈ 3.0) than MPI‑ESM1‑2‑LR (212 vs 80.3, ratio ≈ 2.6; 92.0 vs 14.2, ratio ≈ 6.5). Both models overestimate the volume of suboxic water, so the near-20 µmol kg⁻¹ core is overestimated more severely than the broader 60 µmol kg⁻¹ envelope in relative terms.
- **Disagreement on the scenario change.** The two models agree only in sign at the 20 µmol kg⁻¹ threshold (both contract, −3.0 to −4.1 %), and disagree at 60 µmol kg⁻¹: GFDL‑ESM4 expands (+4.41 %) while MPI‑ESM1‑2‑LR contracts slightly (−0.63 %). The disagreement is not merely one of magnitude: the direction of the change in the broader OMZ is opposite. Even the largest change (+4.41 %) is only about 4.6 × 10¹⁵ m³, far smaller than the inter-model bias in the historical volume (about 107 × 10¹⁵ m³ at 60 µmol kg⁻¹), so the scenario signal is small compared with model structural error.
- **Basin structure.** In WOA23, and in both models, the Pacific dominates the <60 volume; MPI‑ESM1‑2‑LR overestimates all three basins, while GFDL‑ESM4 is close to observed in the Pacific but has too much Atlantic OMZ (6.39 vs 1.29) and too little Indian OMZ (5.43 vs 11.6).
- **Mechanism in the observations.** The volume-weighted AOU of 271.3 µmol kg⁻¹ inside the <60 zone, combined with a saturation of only 11.3 %, means the water has lost nearly 90 % of its saturation capacity and accumulated a large respiration signal. This points to oxygen consumption by remineralisation, i.e. a respiration-driven (and advectively sustained) regime rather than a purely poorly-ventilated, recently isolated water mass, though the two cannot be fully separated with these tracers alone.
- **Inventory.** Global oxygen inventories differ between models by only a few percent (2.20 vs 2.23 × 10¹⁷ mol in the historical runs) and both decrease under ssp245 (to 2.15 and 2.17 × 10¹⁷ mol, i.e. −2.4 % and −2.3 %), a small change relative to the OMZ volume disagreements.

## 5. Figures

**Figure 1** shows three global maps in one row of the historical OMZ thickness at 60 µmol kg⁻¹, for GFDL‑ESM4, MPI‑ESM1‑2‑LR and WOA23. The model thickness fields were first mapped onto the WOA23 1° grid with `regrid_horizontal(target="woa23", method="bilinear")`, then drawn as pcolormesh with a common viridis scale running from 0 to the 98th percentile of the pooled values across the three panels, so the panels are directly comparable and the highest 2 % of thicknesses saturate the colour bar.

**Figure 2** shows three panels in one row: the ssp245-minus-historical difference of the regridded thickness field for GFDL‑ESM4 and MPI‑ESM1‑2‑LR on a shared RdBu_r diverging scale symmetric about zero (clipped at the 98th percentile of |change|, blue = thinning, red = thickening), and, as the third panel, the CESM2 ssp245 thickness alone on its own viridis sequential scale. The CESM2 panel carries an explicit label stating that no change can be formed because the model publishes no historical O₂.

## 6. Caveats

- The thickness DataArrays returned by `identify_oxygen_minimum_zone` carry the parent field’s `units` attribute ("micromoles_per_kilogram" / "umol kg-1") even though the values are column-integrated metres; the maxima quoted above are the metre values (Fig. 1/2 colour bars are labelled in metres).
- WOA23 oxygen, saturation and AOU are the all-period composite; the 1971–2000 oxygen normal was not separately available through the call made.
- GFDL‑ESM4’s O₂ fields contain slightly negative concentrations (minimum −27.2 µmol kg⁻¹ historical, −21.1 µmol kg⁻¹ ssp245), which are counted as below both thresholds.
- Basins exclude Southern and Arctic, so basin volumes sum to less than the global totals (e.g. WOA23 80.3 vs 67.4 + 11.6 + 1.29 = 80.3, the remaining volume falling in the excluded basins).
- Inventories used a constant reference density of 1025 kg m⁻³ and the time-mean mol m⁻³ fields, so they are the model-native concentrations combined with `volcello`, not the µmol kg⁻¹ fields used for the thresholds.
