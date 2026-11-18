## 1. What was done, and what it shows

**Direct answer.** Both CMIP6 models that have a historical `o2` simulation over-predict the observed volume of water below both OMZ thresholds, and they disagree on the sign of the 21st-century change at the 60 µmol kg⁻¹ threshold: GFDL-ESM4 expands (+4.41 %) while MPI-ESM1-2-LR contracts slightly (−0.63 %). GFDL-ESM4 is the better of the two in volumetric terms, mainly because MPI-ESM1-2-LR is roughly 2.6× too voluminous at 60 µmol kg⁻¹ and 6.5× too voluminous at 20 µmol kg⁻¹. The observed WOA23 sub-60 µmol kg⁻¹ zone has a volume-weighted AOU of 271.3 µmol kg⁻¹ and a saturation of 11.3 %, i.e. it is strongly respiration-loaded, not merely under-ventilated.

## 2. Methodology (exactly as executed)

**Data and periods.**
- `load_regular_data` for `cmip6_historical` (`o2`, `thetao`, `so`) for GFDL-ESM4 and MPI-ESM1-2-LR, monthly, 1971-01 to 2000-12.
- `load_regular_data` for `cmip6_ssp245` (`o2`, `thetao`, `so`) for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, monthly, 2071-01 to 2100-12.
- WOA23 `oxygen`, `apparent_oxygen_utilization` and `percent_oxygen_saturation` via `load_climatology_data(dataset="woa23", …, time_slice="annual")`. The WOA23 oxygen reference therefore is the available **all-period composite**; no separate 1971-2000 oxygen normal was assumed (its metadata describe oxygen/nutrients as an all-period composite).

**Unit conversion and time averaging (order of operations).**
- In-situ density from `compute_potential_density_anomaly(temp_data=thetao, sal_data=so, ref_pressure=None, temperature_type="potential")`; the returned field is a density anomaly read from its long_name as sigma + 1000 (cell-by-cell ranges, e.g. 995.0–1057.3 kg m⁻³ for GFDL-ESM4 historical), and `transform_unit(from_unit="mol m-3", to_unit="umol kg-1", density=rho)` converts at each cell's own density. The tool's `density_assumption` attribute was carried.
- Monthly fields were then combined with `temporal_aggregate(freq="mean", op="mean")` over each 30-year period (calendar-month-length weighted), i.e. thresholds were applied to **period-mean converted oxygen**, not to monthly thresholded fields. Period means of the raw `o2` (mol m⁻³) were formed the same way for the inventory.

**OMZ diagnostics.**
- `identify_oxygen_minimum_zone(o2_mean, threshold=20.0/60.0)` on the period-mean µmol kg⁻¹ field gives column thickness (sum of volcello/areacello over sub-threshold cells) and a 3-D mask.
- `compute_exceedance_fraction(…, comparison="less", dim=None, weighting="volume")` gives the volume (m³) below each threshold; divided by 1e15 for the reported units.
- `locate_extreme(thick60, op="max")` gives the maximum thickness and its latitude/longitude.
- `convert_concentration_to_content(concentration=o2_raw_mean)` with each model's `volcello`, summed with `skipna=True`, gives the dissolved-oxygen inventory in moles (volume basis, no density factor).
- Basins: `select_basin(…, basin=b, include_marginal_seas=True, depth_dependent=False)` for Pacific, Indian and Atlantic only, applied **before** `compute_exceedance_fraction`; Southern and Arctic are separate basins and excluded. (The Ofx `basin` field is absent locally, so OMIP's 10-region labels could not be used.)
- WOA23 AOU/saturation: `woa_aou.where(woa_mask60)` and `woa_sat.where(woa_mask60)`, then `compute_basic_stat(op="volume_weighted_mean", dim=["depth","lat","lon"])`.
- Model–observation comparison: `regrid_horizontal(…, target="woa23", method="conservative", skipna=True)` on the historical 60-threshold thickness fields, one `build_comparison_mask(fields=[...], ocean_only=True)`, `.where(mask)` on every field, then `compute_difference_metrics(data=…, reference=woa_thick60_m, dim=["lat","lon"], weighting="area")`.
- Change: 100 × (ssp245 − historical)/historical on unrounded volumes, for GFDL-ESM4 and MPI-ESM1-2-LR only. No CESM2 change was formed: CESM2 historical `o2` is not published, and its ssp245 is r10 whose parent historical r10 is unavailable locally (historical archive is r1).

## 3. Key numerical results

**Volume below threshold (1e15 m³) and maximum 60-threshold thickness**

| Source | vol <20 | vol <60 | max thickness (m) | lat, lon |
|---|---|---|---|---|
| WOA23 | 14.2 | 80.3 | 3125.0 | 7.50, 121.50 |
| GFDL-ESM4 historical | 43.1 | 105 | 4343.2 | 13.50, −93.50 |
| MPI-ESM1-2-LR historical | 92.0 | 212 | 4960.0 | −9.39, −81.22 |
| CESM2 ssp245 | 60.9 | 165 | 3631.9 | 7.45, 121.44 |
| GFDL-ESM4 ssp245 | 41.3 | 110 | 4318.2 | 13.50, −93.50 |
| MPI-ESM1-2-LR ssp245 | 89.2 | 211 | 4950.0 | −9.39, −81.22 |

Observed fractions of ocean volume are 1.03 % (<20) and 5.85 % (<60). WOA23's thickest sub-60 column is in the western Pacific warm-pool/Indonesian region; both models' maxima sit in the eastern tropical Pacific, where their sub-60 layers fill nearly the whole column.

**Percentage change historical → ssp245 (positive = expanding)**

| Model | <20 µmol kg⁻¹ | <60 µmol kg⁻¹ |
|---|---|---|
| GFDL-ESM4 | −4.06 % | +4.41 % |
| MPI-ESM1-2-LR | −3.03 % | −0.63 % |

**Dissolved oxygen inventory (mol)**

| Source | Inventory |
|---|---|
| GFDL-ESM4 historical | 2.20e17 |
| MPI-ESM1-2-LR historical | 2.23e17 |
| CESM2 ssp245 | 2.02e17 |
| GFDL-ESM4 ssp245 | 2.15e17 |
| MPI-ESM1-2-LR ssp245 | 2.17e17 |

Both models lose oxygen between the two periods (GFDL −2.4 %, MPI −2.3 % of inventory, from the tabulated moles), even where the sub-60 volume change is of opposite sign.

**Basin OMZ volume at 60 µmol kg⁻¹ (1e15 m³)**

| Source | Pacific | Indian | Atlantic |
|---|---|---|---|
| WOA23 | 67.4 | 11.6 | 1.29 |
| GFDL-ESM4 historical | 92.8 | 5.43 | 6.39 |
| MPI-ESM1-2-LR historical | 180 | 19.2 | 11.3 |

Basin totals are smaller than the global totals because the Southern and Arctic basins are excluded by construction. The models place too much sub-60 water in the Pacific (GFDL +38 %, MPI +167 % versus WOA23) and too much in the Atlantic (GFDL 5.0×, MPI 8.8× observed), while GFDL-ESM4 under-fills the Indian Ocean (5.43 vs 11.6).

**WOA23 water quality inside the sub-60 zone:** AOU = **271.3 µmol kg⁻¹**, percent saturation = **11.3 %**. AOU is high and saturation is low, so the zone is consistent with respiration-driven oxygen consumption on top of restricted ventilation, not with poor ventilation alone. Caveat: the model OMIP `o2` values may use a different solubility/pressure convention than the WOA23 AOU and saturation products, so observed-vs-model AOU/saturation differences include a definitional component; AOU is reported in µmol kg⁻¹.

**Model vs WOA23 thickness metrics over the shared sample (n = 40 893 cells)**

| Model | bias (m) | RMSE (m) | MAE (m) | max abs error (m) | at lat, lon | mean_data (m) | mean_reference (m) |
|---|---|---|---|---|---|---|---|
| GFDL-ESM4 | 68.61 | 493.15 | 208.75 | 2875.66 | 13.50, −93.50 | 299.61 | 231.00 |
| MPI-ESM1-2-LR | 373.07 | 853.44 | 394.85 | 3946.61 | −20.50, −84.50 | 604.07 | 231.00 |

Both models are biased thick relative to the observed 231.00 m shared-sample mean; MPI-ESM1-2-LR is about 2.6× too thick, GFDL-ESM4 about 1.3×.

**Which model is better, and do they agree?** GFDL-ESM4 reproduces the observed OMZ volume better on both thresholds (mean absolute volumetric error 117.5 % vs 356.6 % for MPI-ESM1-2-LR; at 60 µmol kg⁻¹ the errors are 30.9 % and 164.2 %). The models **do not agree** on the scenario change: at 60 µmol kg⁻¹ GFDL-ESM4 expands by +4.41 % while MPI-ESM1-2-LR contracts by −0.63 %, and both shrink at 20 µmol kg⁻¹ (−4.06 % and −3.03 %). The disagreement at 60 µmol kg⁻¹ is a sign disagreement, and in both cases the change is small compared with the models' bias against the observed volume.

## 4. Figures

**Figure 1** — three global maps in one row of the historical (1971-2000) 60 µmol kg⁻¹ OMZ thickness for GFDL-ESM4, MPI-ESM1-2-LR and WOA23, all three on one shared sequential (viridis) colour scale in metres, from 0 to the pooled 98th percentile so the rare very thick columns do not compress the rest of the range. The intent is a like-for-like visual test of where each model puts thick sub-60 water compared with the observed climatology; both model fields and the WOA23 field were regridded to the WOA23 grid with conservative regridding so the panels are spatially comparable.

**Figure 2** — three global maps in one row. Panels 1 and 2 show the ssp245-minus-historical change in 60 µmol kg⁻¹ OMZ thickness for GFDL-ESM4 and MPI-ESM1-2-LR on one shared diverging (RdBu_r) scale centred on zero, with the symmetric limit set from the pooled 98th percentile of |change|; the shared zero-centred scale is what makes the expansion/contraction contrast between the two models readable. Panel 3 shows CESM2's ssp245 60 µmol kg⁻¹ thickness on its own sequential scale and is explicitly labelled that no change can be formed, because CESM2 has no usable historical counterpart (no published historical `o2`; its ssp245 r10 parent historical r10 is not in the local archive).

## 5. Insights and caveats

- Thresholding was applied to period-mean converted oxygen, not to monthly exceedance, which matters because sub-threshold occurrence is nonlinear; the comparison with the observational climatology is therefore on the same averaging footing.
- CMIP6 `o2` was converted monthly from mol m⁻³ to µmol kg⁻¹ with cell-by-cell in-situ TEOS-10 density (only ±2 % about 1035 kg m⁻³), and the tool's `density_assumption` string states this; the inventory instead uses the raw volume-basis field with `volcello`, so no density assumption enters it.
- GFDL-ESM4's converted field contains small negative oxygen values (minimum −27.2 µmol kg⁻¹ historical, −21.1 ssp245). These are retained and counted as sub-threshold, which inflates its <20 µmol kg⁻¹ volume somewhat; MPI-ESM1-2-LR stays positive (minimum 0.045 µmol kg⁻¹).
- Both models over-state OMZ volume, especially MPI-ESM1-2-LR and especially at the stricter 20 µmol kg⁻¹ threshold, and both put the thickest sub-60 columns in the eastern tropical Pacific rather than the observed western Pacific/Indonesian maximum.
- The two models bracketing the sign of the 60 µmol kg⁻¹ change, combined with a −2.3 to −2.4 % inventory decline in both, supports the view that total ocean oxygen loss is robust while the volume of severely depleted water is model-dependent.
- WOA23 reference volumes rest on the all-period oxygen climatology; the AOU/saturation comparison carries an additional definitional (solubility/pressure convention) uncertainty.
