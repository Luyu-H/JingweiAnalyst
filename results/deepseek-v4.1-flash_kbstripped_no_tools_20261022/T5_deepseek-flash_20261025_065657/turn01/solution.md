## 1. What was done and what it shows

**Data and SSS definition.** Monthly CMIP6 `so` was read for CESM2, GFDL‑ESM4 and MPI‑ESM1‑2‑LR over 1950‑01–2014‑12 (historical) and 2071‑01–2100‑12 (ssp245), together with each experiment's `areacello`, `volcello` and `sftof`. For every model month and for the WOA23 1971–2000 annual climate normal, SSS was defined as the **full-cell-volume-weighted mean of the salinity levels whose depths lie within 0–10 m inclusive** (WOA23: the 0, 5 and 10 m levels, each with its full 5 m layer thickness times the 1° cell area computed from the latitude/longitude cell edges; models: their own `volcello`), with no partial-overlap weighting. WOA23 was restricted to `landsea_01 > 1`. Both datasets are practical salinity, so model minus WOA differences are taken directly in psu with no rescaling.

**Comparison.** The WOA23 1971–2000 SSS was built on its 180×360 1° grid; each model's 1971–2000 mean SSS (time mean of monthly fields) was regridded to that grid with a linear `griddata` interpolation (nearest fill where linear failed), as was its `areacello`. A single common mask (`woa_ocean & isfinite(woa_sss)`) gave **41 088 comparison cells**, applied to all fields before any metric. Metrics are area-weighted on that shared sample, with the regridded model area as weights; the WOA value quoted alongside is therefore the WOA mean *over the model-weighted shared sample*, which is why it is not identical to the WOA global mean. Bias, RMSE and the (area-weighted) spatial Pearson correlation come from that same sample.

**Trends.** Monthly SSS was area-weighted to a global mean with each experiment's `areacello`, then reduced to annual means (65 years for 1950–2014, 30 years for 2071–2100). A linear fit was made separately to each segment with a lag‑1 effective-sample-size correction; the corrected slope, its 95 % CI, R² and the CI-excludes-zero test are reported. The two segments were never fitted as one series and the 2015–2070 gap was left empty.

**ΔSSS.** ssp245 2071–2100 mean minus historical 1995–2014 mean, regridded to 1°, using area-weighted means over valid cells; basins are the WOA23 `basinmask_01` with marginal seas folded in (Atlantic = 1, 6, 7, 8, 11; Pacific = 2; Indian = 3, 9, 10; Southern = 5 and Arctic = 4 excluded from basin reporting but included in the global).

## 2. Key numerical results

**Mean state, 1971–2000 (psu)**

| | model (native grid, area-weighted) | diff. vs WOA | model on shared sample | WOA on shared sample | bias | RMSE | corr |
|---|---|---|---|---|---|---|---|
| WOA23 | **34.63** | — | — | — | — | — | — |
| CESM2 | 34.32 | −0.31 | 34.20 | 34.43 | **−0.23** | **1.34** | **0.767** |
| GFDL‑ESM4 | 34.56 | −0.07 | 34.64 | 34.66 | **−0.02** | **0.77** | **0.889** |
| MPI‑ESM1‑2‑LR | 34.36 | −0.27 | 34.36 | 34.59 | **−0.23** | **1.06** | **0.805** |

All three models are slightly fresh relative to WOA23. GFDL‑ESM4 is by far the closest in the mean, RMSE and pattern (RMSE 0.77 psu, r = 0.889); CESM2 is the freshest and least correlated (RMSE 1.34 psu, r = 0.767). A caveat: the native-grid global means and the shared-sample means differ by up to ≈0.11 psu (largest for CESM2), because the model climatology and its area were both interpolated to 1° and evaluated on the WOA mask; the bias/RMSE/correlation quoted are the shared-sample values. GFDL‑ESM4 is a `gr` model, so its native numbers are already archive-grid numbers, whereas CESM2 (~100 km) and MPI (~250 km) cannot reproduce features narrower than their grids.

**Linear trends of the annual global mean SSS (psu per decade, ESS-corrected, 95 % CI)**

| Model | Historical 1950–2014 | R² | Excl. 0 | ssp245 2071–2100 | R² | Excl. 0 | Ratio |
|---|---|---|---|---|---|---|---|
| CESM2 | −2.68×10⁻³ [−7.47×10⁻³, 2.12×10⁻³] | 0.245 | no (p = 0.22) | **−2.09×10⁻²** [−2.55×10⁻², −1.64×10⁻²] | 0.910 | **yes** (p ≈ 1×10⁻⁶) | 7.8 |
| GFDL‑ESM4 | −1.10×10⁻³ [−3.61×10⁻³, 1.41×10⁻³] | 0.064 | no (p = 0.36) | −5.34×10⁻³ [−1.07×10⁻², 3.82×10⁻⁵] | 0.325 | no (p = 0.051) | 4.9 |
| MPI‑ESM1‑2‑LR | −6.87×10⁻⁴ [−5.48×10⁻³, 4.10×10⁻³] | 0.011 | no (p = 0.76) | −2.29×10⁻² [−5.08×10⁻², 5.07×10⁻³] | 0.684 | no (p = 0.081) | 33.3 |

Effective sample sizes were small (n_eff = 7.8, 15.1, 11.6 historical; 12.3, 12.1, 5.1 for ssp245), reflecting strong serial correlation, so the CIs are wide. The **historical trends are not significant for any model** — all three intervals straddle zero, and the negative point estimates are a fraction of the ssp245 slopes. Only **CESM2's ssp245 trend is significant** (−0.021 psu/decade, R² = 0.91). The **ssp245/historical slope ratios (7.8, 4.9, 33.3) are not meaningful**: the historical denominator is statistically indistinguishable from zero, and both numerator and denominator are short, internally variable segments; they should not be interpreted as an acceleration factor. For CESM2 an additional contamination applies (see caveats).

**ΔSSS, ssp245 2071–2100 minus historical 1995–2014 (psu)**

| Model | Global | Atlantic | Indian | Pacific |
|---|---|---|---|---|
| CESM2 | −0.11 | −0.10 | −0.01 | −0.21 |
| GFDL‑ESM4 | −0.01 | +0.07 | −0.04 | −0.04 |
| MPI‑ESM1‑2‑LR | −0.10 | −0.06 | −0.06 | −0.12 |

End-of-century changes are small and mostly negative (freshening), largest in the Pacific for CESM2 and MPI. GFDL‑ESM4 is nearly neutral globally, and is the only model with a (weakly) **saline** Atlantic. These are basin contrasts, not a strong global signal.

**Member-continuity caveat.** CESM2's historical run is **r1i1p1f1** but its ssp245 run is **r10i1p1f1**, so its ΔSSS and slope ratio mix two different members and therefore include internal-variability differences as well as forcing. GFDL‑ESM4 and MPI‑ESM1‑2‑LR pair **r1i1p1f1** with **r1i1p1f1**, i.e. continuous members.

## 3. Visualizations

- **Figure 1** (three panels, one row, psu): maps of the 1971–2000 SSS bias (model climatology minus WOA23) on the 41 088 common-mask cells, with one shared diverging scale centred on zero (RdBu_r, clipped symmetrically at the pooled 98th percentile of |bias|) and a horizontal colour bar labelled in psu; each panel title carries that model's colour. Intended to show, on a like-for-like mask, where each model is too salty or too fresh relative to the WOA23 climate normal.
- **Figure 2** (three panels, one row): the cell-by-cell linear trend of annual-mean SSS over 1950–2014 in psu per decade, after the fitted-field and p-value fields were put on a common regular 1° longitude–latitude frame for display only (linear interpolation for the trend, nearest for the p-value field), masked to the common comparison mask. One shared diverging scale centred on zero, clipped at the pooled 98th percentile of |trend|; cells with p < 0.05 (ESS-corrected) are stippled over the trend field. Intended to separate the magnitude and sign of local trends from where they are statistically robust.
- **Figure 3** (single panel, psu): annual global-mean SSS for all three models, historical 1950–2014 and ssp245 2071–2100 drawn as two separate segments with the 2015–2070 gap empty, each segment overlaid with its own linear fit (dashed, same model colour), y axis labelled in psu, models coloured consistently with Figures 1–2.

## 4. Insights and caveats

1. **Mean-state skill is model dependent:** GFDL‑ESM4 is close to WOA23 in the mean (bias −0.02 psu) and in pattern (r = 0.89); CESM2 and MPI‑ESM1‑2‑LR are both about 0.23 psu fresh on the shared sample with RMSE ≈ 1.1–1.3 psu and r ≈ 0.77–0.80. Part of the CESM2/MPI RMSE is genuinely resolution-limited structure (≈100 km and ≈250 km grids) rather than a large-scale bias.
2. **Trends are hard to detect historically:** after ESS correction the 1950–2014 interval always includes zero, while the ssp245 segment has a steeper negative slope but is significant only for CESM2. The slope ratio is therefore uninformative.
3. **ΔSSS is small and basin-dependent:** the global freshening is 0.01–0.11 psu, with the Pacific contributing most; the Atlantic is near zero or, for GFDL, weakly saline. Because these are 20–30-year means and CESM2 mixes members (r1 historical, r10 ssp245), basin values should be read as indicative rather than precise.
4. **Consistency of convention:** all model–WOA and model–model comparisons use the same 0–10 m full-volume weighting, the same 41 088-cell mask and the same area weighting, so differences between models are not artefacts of differing masks or weighting; however the regridding to 1° (linear with nearest fill) is itself a source of about 0.1 psu in the global means, as seen in the CESM2 native (34.32) versus shared-sample (34.20) values.
