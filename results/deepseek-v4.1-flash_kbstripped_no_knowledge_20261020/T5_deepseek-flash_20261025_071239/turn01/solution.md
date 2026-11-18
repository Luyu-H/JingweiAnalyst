## 1. Bottom line

The workflow ran to completion and produced all three figures, but **the absolute salinity numbers are not physical**. The CMIP6 `so` fields were already stored as grams per kilogram (numerically ≈ psu: global mean 34.1, range 6.9–48.2), yet carried the attribute `units = "0.001"`. The helper `to_psu()` matched that attribute and multiplied the fields by 1000, so every model SSS field is **1000× too large**, while the WOA23 field (units `"1"`, sample check > 1) was left unscaled. Consequently:

- Model global means and their bias against WOA23, the bias maps (Fig. 1) and the RMSE are internally inconsistent (model ~3.4×10⁴ "psu" vs. WOA23 34.63 psu) and the RMSE is **not recoverable** from the printed value.
- Quantities that are invariant to a constant positive scaling remain valid: the spatial Pearson correlations, the ssp245/historical slope ratio, and the *sign/significance* of trends. Quantities that are linear in the model field alone (global means, global and basin ΔSSS, trend slopes) can be brought back to psu by dividing by 1000; I give those rescaled values below, clearly flagged.
- The figures are plotted on the erroneous scale (Fig. 1 colour bar in "thousands of psu", Fig. 3 y-axis ≈34 000 instead of ≈34.3).

## 2. Methodology exactly as executed

**Data and periods.** `load_regular_data` for `cmip6_historical` `so`, monthly 1950‑01 to 2014‑12 (780 months), and `cmip6_ssp245` `so`, monthly 2071‑01 to 2100‑12 (360 months), for CESM2, GFDL‑ESM4 and MPI‑ESM1‑2‑LR, both loaded with `depth_range=[0, 10]`. Per model: `load_volumes_for_cmip6(variable='so')` for volcello, `load_areas(variable='so')` for areacello, `load_depth_levels(variable='so')`. Recorded variants: historical r1i1p1f1 for all three; ssp245 r1i1p1f1 for GFDL‑ESM4 and MPI‑ESM1‑2‑LR but **r10i1p1f1 for CESM2**.

**Levels retained by the 0–10 m filter.** CESM2: 1 level (5 m); GFDL‑ESM4: 2 levels (2.5 and 10 m); MPI‑ESM1‑2‑LR: 1 level (6 m); WOA23: 3 levels (0, 5, 10 m) with bounds 0–2.5, 2.5–7.5, 7.5–12.5 m, i.e. full-cell thicknesses 2.5, 5, 5 m (the 10 m cell extends to 12.5 m).

**SSS definition.** A depth-weighted mean over the retained levels, with *full cell volume* weights, no overlap weighting with the 0–10 m slab: models use Σ(so·volcello)/Σ(volcello) over valid levels (GFDL‑ESM4 volcello is time-varying, 1980/1032 time steps); WOA23 uses Σ(S·A·Δz)/Σ(A·Δz) with `woa_area` and the level thickness from `load_depth_levels`. Same construction for both, but as noted the model fields were scaled ×1000 by `to_psu()`.

**WOA23.** `load_climatology_data(dataset='woa23', variable='salinity', period='1971-2000', time_slice='annual')`; `load_landsea_mask()` (ocean where value > 1, used only via the mask/comparison cells); `load_basin_mask(depth_levels='0m')`.

**1971–2000 mean state.** `temporal_aggregate(sss_hist, freq='mean', op='mean')` over 1971‑01 to 2000‑12, then `compute_basic_stat(op='area_weighted_mean', dim=['lat','lon'])`; WOA23 global mean by the same call on the SSS field. (areacello/1° areas were loaded but **not** passed explicitly to `compute_basic_stat`; for CESM2 the dims are `nlat/nlon` while `dim=['lat','lon']` was requested, a further silent-mismatch risk.)

**Bias metrics.** Model climatology → `regrid_horizontal(target='woa23', method='conservative')` → `build_comparison_mask([model, woa], ocean_only=True)` → `compute_difference_metrics(..., weighting='area')` for bias/RMSE and `analyze_correlation(..., method='pearson', weighting='area')` for r.

**Trends.** Monthly global mean via `compute_basic_stat(op='area_weighted_mean', dim=['lat','lon'])`, `temporal_aggregate(freq='year')`, then `fit_1D_trend(method='linear')`; the reported 95 % CIs use an autocorrelation-adjusted effective sample size (n_eff ≈ 5–15 of 30–65), slopes converted to psu/decade by ×10.

**Trend maps.** `temporal_aggregate(freq='year')` then `fit_trend_map(dim=year, per_decade=True, alpha=0.05)`; slope regridded conservatively to the WOA23 1° grid for display, significance regridded with `nearest_s2d` and thresholded at > 0.5.

**ΔSSS.** ssp245 2071–2100 mean minus historical 1995–2014 mean. Global: area-weighted mean of the difference on the native model grid. Basins: each period-mean regridded conservatively to WOA23 1°, differenced, then `select_basin(basin=..., include_marginal_seas=True, depth_dependent=False)` for Atlantic, Indian, Pacific and area-weighted means. Southern and Arctic are excluded.

## 3. Key numerical results

**1971–2000 global-mean SSS (psu; model values divided by 1000 to undo the scaling error)**

| | global mean | minus WOA23 |
|---|---|---|
| WOA23 | 34.63 | — |
| CESM2 | 34.32 | −0.32 |
| GFDL‑ESM4 | 34.56 | −0.07 |
| MPI‑ESM1‑2‑LR | 34.36 | −0.27 |

Raw (scaled) output was 34315.64, 34560.23, 34359.29 and 34.63 for WOA23, giving the meaningless "biases" 34281.01, 34525.60, 34324.66. All three models are slightly fresher than the WOA23 normal.

**Spatial metrics on the common WOA23 1° ocean/valid cells** (correlations are unaffected by the scaling and are usable):
CESM2 r = 0.901; GFDL‑ESM4 r = 0.924; MPI‑ESM1‑2‑LR r = 0.891. RMSE as printed (34364.26, 34632.87, 34405.81) is **invalid** — the reference was unscaled — and cannot be rescued algebraically; it must be recomputed after fixing units.

**Historical 1950–2014 global-mean trend (psu/decade, rescaled ÷1000; the printed values were −2.68, −1.10, −0.551 "psu/decade")**

| model | slope (psu/dec) | 95 % CI | R² | CI excludes 0 |
|---|---|---|---|---|
| CESM2 | −2.68×10⁻³ | (−6.93×10⁻³, +1.57×10⁻³) | 0.245 | no (p = 0.18) |
| GFDL‑ESM4 | −1.10×10⁻³ | (−3.61×10⁻³, +1.41×10⁻³) | 0.064 | no (p = 0.36) |
| MPI‑ESM1‑2‑LR | −5.51×10⁻⁴ | (−5.04×10⁻³, +3.93×10⁻³) | 0.007 | no (p = 0.79) |

**ssp245 2071–2100 trend (rescaled from −20.9, −5.34, −23.1)**

| model | slope (psu/dec) | 95 % CI | R² | CI excludes 0 |
|---|---|---|---|---|
| CESM2 | −2.09×10⁻² | (−2.55×10⁻², −1.64×10⁻²) | 0.910 | yes (p ≈ 1×10⁻⁶) |
| GFDL‑ESM4 | −5.34×10⁻³ | (−1.06×10⁻², −3.9×10⁻⁵) | 0.324 | yes, marginal (p = 0.049) |
| MPI‑ESM1‑2‑LR | −2.31×10⁻² | (−5.13×10⁻², +5.04×10⁻³) | 0.701 | no (p = 0.079) |

**Ratio ssp245/historical slope** (scale-invariant, so these are valid as printed): CESM2 7.8, GFDL‑ESM4 4.9, MPI‑ESM1‑2‑LR 42.0. **The ratio is not meaningful**: the two slopes cover different periods and mean states, are fitted to different forcing, carry independent uncertainties, and here the historical slopes are individually indistinguishable from zero — MPI's 42.0 is an artefact of dividing by a near-zero denominator.

**ΔSSS, ssp245 2071–2100 minus historical 1995–2014 (psu, rescaled ÷1000)**

| model | global | Atlantic | Indian | Pacific |
|---|---|---|---|---|
| CESM2 | −0.12 | +0.10 | +0.06 | −0.24 |
| GFDL‑ESM4 | −0.01 | +0.19 | −0.00 | −0.04 |
| MPI‑ESM1‑2‑LR | −0.09 | −0.04 | −0.05 | −0.11 |

Global values are not the area-weighted sum of the three basins: the global mean includes the Southern and Arctic Oceans (excluded from the basin analysis) and is taken on the native grid rather than after regridding.

**Ensemble-member caveat (as recorded by the code).** GFDL‑ESM4 and MPI‑ESM1‑2‑LR use r1i1p1f1 for both experiments; **CESM2 uses r1i1p1f1 historically but r10i1p1f1 for ssp245**. CESM2's ΔSSS and slope ratio therefore mix member-to-member variability with the forced response, and its basin pattern (Atlantic/Indian salinification with Pacific freshening) should not be read as a forced signal.

## 4. Visualizations

- **Figure 1** — three global panels (one per model) of the 1971–2000 SSS bias (model climatology minus WOA23) on the common comparison cells after conservative regridding to the WOA23 1° grid. A single diverging `RdBu_r` scale centred on zero is shared by all panels, clipped at the 98th percentile of |bias|, with a colour bar labelled in psu. Logic: whether each model is too salty or too fresh relative to the 1971–2000 WOA23 normal, and where the errors are structured. Because of the scaling error the fields are ~10³ psu in magnitude and the panel content is dominated by the offset rather than the physical pattern.
- **Figure 2** — three panels of the cell-by-cell linear trend of annual-mean SSS over 1950–2014 from `fit_trend_map(per_decade=True)`, conservatively regridded to a common regular longitude–latitude frame for display only, one shared diverging scale centred on zero in psu/decade, clipped at the 98th percentile, with p < 0.05 cells stippled (contourf hatching of the regridded significance field). Logic: spatial structure of historical surface-salinity change and where it is statistically distinguishable from zero.
- **Figure 3** — one panel: annual global-mean SSS series for all three models, historical 1950–2014 and ssp245 2071–2100 as separate segments with 2015–2070 left empty (shaded grey band marks the gap), each segment overlaid with its own `fit_1D_trend` fitted line (dashed for historical, dotted for ssp245), model colours kept consistent across all three figures, y-axis intended in psu (actually in the scaled, erroneous units).

## 5. Patterns, and what must be fixed

- All three models are slightly **fresher than the WOA23 1971–2000 normal** in the global mean (−0.32, −0.07, −0.27 psu after rescaling), and all reproduce the large-scale salinity pattern well (r = 0.89–0.92). CESM2 and MPI have larger global offsets than GFDL‑ESM4.
- **Every fitted global-mean trend is negative**, in both experiments: the historical global surface ocean is weakly freshening in these single members (none significant once autocorrelation is accounted for), and the ssp245 segment freshens faster in CESM2 and MPI (2×10⁻² psu/decade, i.e. ~0.06 psu over 30 years). GFDL‑ESM4's ssp245 trend is significant only marginally.
- **Ranking is consistent**: GFDL‑ESM4 is the least biased globally, has the highest pattern correlation, the smallest ΔSSS, and the weakest trends; MPI‑ESM1‑2‑LR and CESM2 are noisier (MPI historical R² = 0.007).
- **Required fixes before these numbers can be quoted as science:** (i) remove the ×1000 conversion for `so` (values are already g kg⁻¹ ≈ psu) while keeping the WOA23 field as is — this is what corrupts the global means, the bias maps, the RMSE, the ΔSSS magnitudes and the Fig. 1/Fig. 3 axes; (ii) recompute RMSE, which is not recoverable from the printed value; (iii) pass areacello / 1° areas explicitly as weights and reconcile CESM2's `nlat/nlon` dims with the `dim=['lat','lon']` argument; (iv) note that CESM2's ΔSSS and slope ratio mix experiment members (r1 vs r10). Correlations and slope ratios already reported stand as computed.
