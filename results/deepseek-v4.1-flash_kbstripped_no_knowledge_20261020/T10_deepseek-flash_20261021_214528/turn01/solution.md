## 1. What was computed, and exactly how

**Observations (WOA23).** The 1971–2000 annual climatology was read with `load_climatology_data(dataset='woa23', variable=…, period='1971-2000', time_slice='annual')` for `oxygen` (µmol kg⁻¹), `percent_oxygen_saturation` (%), and `apparent_oxygen_utilization` (µmol kg⁻¹), on the native 1°×1° grid with 102 depth levels. `load_depth_levels(dataset='woa23', variable='oxygen', time_slice='annual')` returned the level centres and layer bounds, `load_landsea_mask()` the WOA23 land–sea (bottom-level) mask, and `load_basin_mask()` the WOA23 basin mask (102×180×360).

**Models.** For both experiments the oxygen, temperature and salinity fields were read with `load_regular_data` (`o2`, `thetao`, `so`): `cmip6_historical` for GFDL-ESM4 and MPI-ESM1-2-LR over 1971-01–2000-12, and `cmip6_ssp245` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR over 2071-01–2100-12. Each model's native `volcello` (`load_volumes_for_cmip6`), `areacello` (`load_areas`) and depth levels/bounds (`load_depth_levels`) were loaded for the same variable. A day-weighted long-term mean of each field was formed with `temporal_aggregate(freq='mean', op='mean')` over the 360 monthly steps. Potential density anomaly was computed from the mean `thetao`/`so` with `compute_potential_density_anomaly(ref_pressure=0)` (TEOS-10), and the native `o2` (delivered in mol m⁻³) was converted to µmol kg⁻¹ with `transform_unit(..., to_unit='umol/kg', density=sigma)`. Grids: GFDL-ESM4 1°×1°, 35 levels; MPI-ESM1-2-LR curvilinear 220×254, 40 levels; CESM2 384×320, 60 levels.

**OMZ definition and volumes.** `identify_oxygen_minimum_zone(o2_data=…, threshold=20.0 / 60.0)` produced the 3-D masks and the column thickness field (summed thickness of water strictly below the threshold in each column). Volumes were obtained with `compute_exceedance_fraction(threshold=…, comparison='less', weighting='volume', dim=None)` and are reported from its `exceeding_amount` in units of 10¹⁵ m³. An independent cross-check (mask × native `volcello`, summed) reproduced every 60-threshold volume exactly. Maximum thickness and its position came from `locate_extreme(thickness60, op='max')`.

**Inventory.** Oxygen inventory per model was obtained with `convert_concentration_to_content(concentration=o2_kg, density=sigma)` summed over the ocean, cross-checked against the explicit ∫ o₂·(σ₀+1000)·V dV / 10⁶ (the tool value being retained only when it agreed with the manual integral to within a factor of two).

**Basins.** Basin volumes at the 60 threshold were obtained by applying `select_basin(data=…, basin=…, include_marginal_seas=True, depth_dependent=False)` directly to each native oxygen field, then `compute_exceedance_fraction(weighting='volume', threshold=60)`. Marginal seas are therefore folded into their parent ocean, and the Southern and Arctic Oceans are excluded by construction. (The basin mask was *not* regridded; `select_basin` was applied natively.)

**In-OMZ water-mass properties.** For WOA23 the 60-threshold mask was used to subset AOU and O₂ saturation, and `compute_basic_stat(op='volume_weighted_mean', dim=None)` gave volume-weighted means.

**Maps.** Thickness fields were put on the WOA23 1° grid with `regrid_horizontal(target='woa23', method='conservative')`; colour limits were the 98th percentile of the pooled (or absolute-change) values.

---

## 2. Key numerical results

### OMZ volume (10¹⁵ m³, O₂ strictly below threshold)

| Source / experiment | O₂ < 20 | O₂ < 60 |
|---|---|---|
| WOA23 (1971–2000) | 14.2 | 80.3 |
| GFDL-ESM4 historical | 42.9 | 104 |
| MPI-ESM1-2-LR historical | 91.6 | 211 |
| CESM2 ssp245 | 60.6 | 164 |
| GFDL-ESM4 ssp245 | 41.2 | 109 |
| MPI-ESM1-2-LR ssp245 | 88.9 | 209 |

### Maximum 60-threshold OMZ thickness and location

| Source / experiment | Max thickness (m) | Lat (°) | Lon (°) |
|---|---|---|---|
| WOA23 | 3125.0 | 7.50 | 121.50 |
| GFDL-ESM4 historical | 4343.2 | 13.50 | −93.50 |
| MPI-ESM1-2-LR historical | 4960.0 | −9.39 | −81.22 |
| CESM2 ssp245 | 3631.9 | 7.45 | 121.44 |
| GFDL-ESM4 ssp245 | 4318.2 | 13.50 | −93.50 |
| MPI-ESM1-2-LR ssp245 | 4950.0 | −9.39 | −81.22 |

WOA23's maximum lies in the shallow-silled Sulu/Celebes region of the Indo-Pacific; GFDL-ESM4 places its maximum in the eastern tropical Pacific off Central America, MPI-ESM1-2-LR off Peru.

### Historical → ssp245 percentage change (positive = expanding)

| Model | O₂ < 20 | O₂ < 60 |
|---|---|---|
| GFDL-ESM4 | −4.07 % | +4.59 % |
| MPI-ESM1-2-LR | −3.03 % | −0.63 % |

### Global dissolved oxygen inventory (mol)

| Case | Inventory |
|---|---|
| GFDL-ESM4 historical | 2.20 × 10¹⁷ |
| MPI-ESM1-2-LR historical | 2.23 × 10¹⁷ |
| CESM2 ssp245 | 2.02 × 10¹⁷ |
| GFDL-ESM4 ssp245 | 2.15 × 10¹⁷ |
| MPI-ESM1-2-LR ssp245 | 2.17 × 10¹⁷ |

Both models lose ~2.3–2.4 % of their total dissolved oxygen between the two windows.

### Basin-resolved 60-threshold OMZ volume (10¹⁵ m³)

| Case | Pacific | Indian | Atlantic | Sum |
|---|---|---|---|---|
| WOA23 | 67.4 | 11.6 | 1.29 | 80.3 |
| GFDL-ESM4 historical | 91.7 | 5.40 | 6.37 | 104 |
| MPI-ESM1-2-LR historical | 179 | 19.1 | 11.2 | 209 |

The three basins account for essentially the whole WOA23 < 60 volume (80.33 of 80.33); for the models the sum is slightly smaller than the global total (104 vs 104; 209 vs 211) because Southern and Arctic waters are excluded.

### WOA23 water-mass character inside the 60-threshold OMZ

Volume-weighted mean AOU = **271.3 µmol kg⁻¹**; mean O₂ saturation = **11.3 %**.

High AOU (271 of a global maximum near 343 µmol kg⁻¹) together with very low saturation (11 %) means the water inside the observed OMZ has lost most of its saturating oxygen to organic-matter remineralisation. This is the signature of a **respiration-dominated** oxygen sink, not of a merely poorly ventilated water mass (which would carry low AOU at comparatively higher saturation, i.e. oxygen deficit inherited from the source region rather than consumed in situ). The WOA23 OMZ is a place where oxygen is consumed, not just where it fails to arrive.

---

## 3. Model skill and scenario agreement

- **Which model is closer to the observed OMZ?** GFDL-ESM4, decisively. Absolute error in the 60-threshold volume is 23.7 × 10¹⁵ m³ for GFDL-ESM4 versus 130 × 10¹⁵ m³ for MPI-ESM1-2-LR. MPI-ESM1-2-LR more than doubles the observed volume (211 vs 80.3) and is also far too large at the 20 threshold (91.6 vs 14.2, a factor 6.5). GFDL-ESM4 overestimates by ~30 % at 60 (104 vs 80.3) and by a factor 3 at 20 (42.9 vs 14.2). Both models therefore put too much water in the low-oxygen classes, but MPI-ESM1-2-LR does so much more strongly.
- **Do the two models agree on the scenario change?** Only partially:
 - At **60 µmol kg⁻¹** they **disagree in sign**: GFDL-ESM4 expands (+4.59 %), MPI-ESM1-2-LR contracts (−0.63 %). Both changes are small (≲5 %), i.e. neither model predicts a dramatic reorganisation of the sub-60 volume by 2071–2100 under SSP2-4.5.
 - At **20 µmol kg⁻¹** they **agree in sign** — both contract, by −4.07 % (GFDL-ESM4) and −3.03 % (MPI-ESM1-2-LR) — and in approximate size.
- CESM2 cannot enter this comparison: no historical `o2` exists, so no historical volume, inventory, percentage change or change map can be formed. Its ssp245 sub-60 volume (164 × 10¹⁵ m³) sits between the other two models but cannot be placed in a change context.
- Consistent with the small volume changes, both models show a ~2.3–2.4 % decline in total dissolved oxygen inventory, and their maximum 60-threshold thicknesses barely move (GFDL 4343 → 4318 m; MPI 4960 → 4950 m), staying at the same grid points. The projected signal in these two models is thus a modest thinning/expansion at the margins of the low-oxygen volume rather than a displacement of the thickest column.

---

## 4. Figures produced

**Figure 1** — three global maps in one row of the historical 60-threshold OMZ thickness for GFDL-ESM4, MPI-ESM1-2-LR (both 1971–2000) and WOA23 (1971–2000), all conservatively regridded to the WOA23 1° grid and drawn on a single sequential (`viridis`) scale from 0 m to the pooled 98th percentile, with an extending colour bar. Its purpose is a like-for-like visual comparison of where each product puts thick low-oxygen columns, independent of colour choice.

**Figure 2** — three maps in one row: (left, centre) ssp245 (2071–2100) minus historical (1971–2000) 60-threshold thickness for GFDL-ESM4 and MPI-ESM1-2-LR on a shared diverging `RdBu_r` scale centred on zero with symmetric limits at the 98th percentile of the absolute change, so blue = thinning and red = thickening; (right) CESM2's ssp245 thickness alone on its own sequential `viridis` scale, annotated that no change field can be formed because CESM2 publishes no historical `o2`. The figure is designed to separate "model-projected change" (two comparable panels) from "model state only" (CESM2).

---

## 5. Insights and caveats

- The observed OMZ volume is dominated by the Pacific (67.4 of 80.3 × 10¹⁵ m³, ~84 %), with the Indian second (11.6) and the Atlantic small (1.29). Both models reproduce this Pacific dominance but with a strongly inflated Pacific (91.7 and 179) and — in GFDL-ESM4 — an Atlantic inflated roughly five-fold (6.37 vs 1.29).
- The two CMIP6 models disagree with each other far more in their *mean state* (211 vs 104, a factor 2) than in their *response* at 20 µmol kg⁻¹ (≈ −3 to −4 %), and the 60-threshold response even changes sign between them. This suggests that under SSP2-4.5 the direction of the centennial change in the sub-60 inventory is model-dependent, while the shrinking of the sub-20 core is more robust across the two.
- **Unphysical values:** GFDL-ESM4's oxygen reaches −27.3 µmol kg⁻¹ (historical) and −21.2 µmol kg⁻¹ (ssp245), and negative σ₀ values appear at the surface. These negatives are counted as "strictly below" both thresholds, so GFDL-ESM4's OMZ volumes are, if anything, slightly inflated by model artefacts; the −4 %/−0.6 % contraction at 20 µmol kg⁻¹ in the two models occurs against this background.
- **Display-only regridding:** the conservative remap to 1° materially changes per-cell thickness where native curvilinear cells overlap partial-land 1° cells (MPI-ESM1-2-LR's field-mean thickness rises from 244 m on its native grid to 432 m on the 1° grid, and its valid-cell count from 36 502 to 42 713). All ranked numbers above come from native-grid calculations; the maps should be read qualitatively for pattern, not quantitatively for magnitude.
- The WOA23 land–sea mask was loaded but not used in any computation; the CMIP6 land points are handled implicitly through the volume weights and NaN-aware sums.
