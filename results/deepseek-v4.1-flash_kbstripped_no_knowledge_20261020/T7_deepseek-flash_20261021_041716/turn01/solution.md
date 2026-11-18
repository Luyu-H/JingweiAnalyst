## 1. What the executed code actually did

**Sources and period.** Three CMIP6 historical models (CESM2, GFDL‑ESM4, MPI‑ESM1‑2‑LR) were loaded with `load_regular_data` for `no3` and `po4` over 1971‑01 to 2000‑12 and a 200–1000 m `depth_range`, then reduced to a long-term mean with `temporal_aggregate(freq="mean", op="mean")` (day-weighted). WOA23 `nitrate` and `phosphate` were taken from the annual objectively-analysed climatology with `load_climatology_data(time_slice="annual", depth_range=[200,1000])` and used directly; no averaging over time was applied. The WOA23 product is an all-period climatology, so the model (1971–2000) and observed fields are **not synchronous** — a caveat that qualifies every comparison below.

**Units.** Model concentrations were converted with `transform_unit(..., to_unit="umol/kg", density=1.025)` i.e. the prescribed constant density of 1025 kg m⁻³. WOA23 values are already µmol kg⁻¹.

**Weighting.** All means and exceedance fractions used `compute_basic_stat(op="volume_weighted_mean")` / `compute_exceedance_fraction(weighting="volume")`, i.e. volume weighting on each source's own grid. No explicit `volcello`/`cell-area` array was passed in the executed code; the weights came from the statistic's internal volume treatment. The regression was the only unweighted step (each valid 3‑D cell counted once; because model grids differ in resolution, cell counts — 1.56 M valid CESM2 cells vs 0.33 M MPI — weight models differently in that fit).

**Metrics.** Cell-by-cell N:P via `compute_ratio` on cells where both tracers were valid and PO₄ > 0; the 14–18 band via `compute_exceedance_fraction(threshold=2.0, comparison="abs_less", reference=16.0)` (i.e. |N:P−16| < 2, the strict 14–18 interval); the low-N:P fraction via `threshold=12.0, comparison="less"`. Basins via `select_basin(basin=..., include_marginal_seas=True)` for Pacific, Atlantic, Indian only.

## 2. Headline results (200–1000 m)

| Source | Mean cell-wise N:P | NO₃/PO₄ ratio-of-means | mean NO₃ | mean PO₄ | slope a | intercept b | R² | % vol 14–18 | % vol <12 |
|---|---|---|---|---|---|---|---|---|---|
| WOA23 | **14.15** | 14.00 | 27.49 | 1.96 | 13.04 | 1.69 | 0.908 | 55.0 | 6.2 |
| CESM2 | **11.83** | 10.93 | 21.31 | 1.95 | 4.89 | 10.84 | 0.238 | 17.2 | 36.7 |
| GFDL‑ESM4 | **14.13** | 13.65 | 24.47 | 1.79 | 11.37 | 4.02 | 0.817 | 74.6 | 7.3 |
| MPI‑ESM1‑2‑LR | **12.28** | 11.77 | 20.38 | 1.73 | 8.59 | 5.64 | 0.678 | 0.5 | 22.5 |

(Concentrations in µmol kg⁻¹; b in µmol kg⁻¹.)

## 3. Basin means (headline cell-wise N:P) and low-N:P volumes

| Source | Pacific | Atlantic | Indian | % vol <12 Pacific | % vol <12 Indian |
|---|---|---|---|---|---|
| WOA23 | 13.57 | 15.55 | 13.73 | 8.1 | 10.3 |
| CESM2 | 10.58 | 13.55 | 11.18 | 59.4 | 35.8 |
| GFDL‑ESM4 | 13.39 | 15.07 | 14.72 | 14.9 | 0.5 |
| MPI‑ESM1‑2‑LR | 11.29 | 13.15 | 12.60 | 40.0 | 15.0 |

## 4. Direct answer: where each source sits relative to Redfield, and who reproduces the denitrification signal

- **WOA23 — below Redfield (14.15 vs 16).** The observed intermediate water is sub-Redfield, with a genuine low-N:P tail: 6.2 % of the global 200–1000 m volume has N:P < 12, concentrated in the Pacific (8.1 %) and Indian (10.3 %), and only 55 % of the volume lies in the 14–18 band. Basin ranking (Atlantic 15.55 > Indian 13.73 ≈ Pacific 13.57) is the expected imprint of fixed-nitrogen loss in the Pacific/Indian versus the better-ventilated Atlantic.
- **GFDL‑ESM4 — at Redfield (14.13).** It is the only model matching the observed global mean and the observed global low-N:P volume (7.3 % vs 6.2 %), and its Atlantic (15.07) and Pacific (13.39) means are close to WOA23. However it **does not reproduce the basin distribution**: Pacific low-N:P volume is nearly double the observed (14.9 % vs 8.1 %) while the Indian is essentially absent (0.5 % vs 10.3 %). Its internal spread is also too narrow — 74.6 % of volume in the 14–18 band versus 55 % observed — i.e. the model reproduces the mean but not the observed variance/denitrification geography.
- **CESM2 — well below Redfield (11.83).** The low-N:P signal is grossly over-expressed: 36.7 % of global volume below 12 (six times observed), 59.4 % in the Pacific and 35.8 % in the Indian. The weak, shallow slope (a = 4.89) with a large offset (b = 10.84) and R² = 0.238 shows nitrate is almost decoupled from phosphate — a nitrate pool drawn down to a near-constant low level rather than a spatially structured denitrification pattern. Basin values are all ~2–3 units too low.
- **MPI‑ESM1‑2‑LR — below Redfield (12.28).** Low-N:P volume is ~3.5× observed globally (22.5 %), with the excess in both the Pacific (40.0 %) and Indian (15.0 %); the basin ordering (Pacific 11.29 < Indian 12.60 < Atlantic 13.15) is qualitatively right but every basin is displaced downward. Almost no water falls in the 14–18 band (0.5 %), so the model's N:P distribution is centred on ~12 rather than straddling Redfield with a low-N:P tail.

**Verdict.** No model reproduces both the magnitude and the spatial structure of the observed low-N:P denitrification signal. GFDL‑ESM4 gets the global mean and global low-N:P fraction right but puts the denitrification in the wrong basins (Pacific instead of Indian) and under-represents the spread. CESM2 and MPI‑ESM1‑2‑LR both sit below Redfield and reproduce the *sign* of the signal (Pacific low-N:P > global > Indian), but with strongly exaggerated amplitude and, in CESM2's case, a nitrate–phosphate relationship that is almost decoupled (slope far below Redfield, large positive intercept).

A consistent mechanistic reading: modelled mean phosphate is close to observations (1.73–1.96 µmol kg⁻¹), whereas modelled mean nitrate is 20.4–24.5 versus 27.5 µmol kg⁻¹ observed. The low modelled N:P is therefore driven by **nitrate deficiency**, not phosphate excess, and the depressed regression slopes (4.9–11.4 versus 13.0) indicate that this deficiency is not organised like the observed denitrification pattern.

## 5. Figures

**Figure 1 (4-panel NO₃ vs PO₄ scatter, order WOA23, CESM2, GFDL‑ESM4, MPI‑ESM1‑2‑LR).** Each panel is a density-preserving random subsample (up to 20 000 cells, seed 0) of the valid 200–1000 m PO₄/NO₃ pairs, plotted with identical axis limits and tick sets derived from the pooled 99.5th percentile of all four sources, so panels are directly comparable. Each panel overlays the source's own equal-weight OLS fit (solid) and the Redfield 16:1 line through the origin (dashed), plus a box with a, b, R² and the headline mean N:P; a legend distinguishes cells, fit and Redfield line. The figure's logic is to expose how tightly nitrate tracks phosphate along a Redfield-like ray, and where the cloud departs from that ray.

**Figure 2 (2×2 global maps of 200–1000 m mean N:P − 16).** For each source, N:P is averaged over the depth dimension with volume weighting to give one value per water column, 16 is subtracted, and the field is drawn as a true geographic map on a common longitude–latitude frame. Model fields are bilinearly regridded to a 1° grid (`regrid_horizontal`, display only); WOA23 is left on its native 1° grid. All four panels share one diverging colour scale centred on zero, clipped at the 98th percentile of the pooled |anomaly| (with a floor of 0.5) so the Redfield-relative sign and spatial pattern are comparable across sources.

## 6. Important insights and caveats

- The observed N:P is not Redfield in intermediate water: at 14.15 with 6.2 % of volume below 12 and more than half the volume in 14–18, the observation is a sub-Redfield distribution with a distinct low tail, not a narrow 16:1 mode.
- The gap between the cell-wise mean and the ratio-of-means (14.15 vs 14.00 for WOA23; 11.83 vs 10.93 for CESM2) is small but systematic and larger in the models, reflecting skewed, low-ratio tails that the ratio-of-means damps.
- R² and slope together are a useful diagnostic: WOA23 (a = 13.04, R² = 0.908) is a tightly coupled near-Redfield system slightly offset by denitrification; GFDL‑ESM4 (11.37, 0.817) is similar but too homogeneous; MPI (8.59, 0.678) and especially CESM2 (4.89, 0.238) have progressively more nitrate variability unexplained by phosphate.
- Non-physical values exist in the raw model fields (CESM2 nitrate minimum −1.24 µmol kg⁻¹); ratios were only formed where both tracers were valid and PO₄ > 0, and none of the reported statistics were masked or repaired beyond that.
- WOA23's per-column anomaly field shows extreme positive values (up to ~119) at a few cells where near-zero phosphate inflates the ratio; the pooled 98th-percentile clipping in Figure 2 keeps those from swamping the shared colour scale.
- Because the regression is unweighted per cell, models with finer grids contribute more cells to their own fit, so cross-model intercept/slope comparisons mix resolution effects with genuine stoichiometric differences; the volume-weighted means and volume fractions do not have this issue.
- Finally, the model fields are 1971–2000 climatologies while WOA23 is an all-period analysis, so part of any mean offset may reflect the different sampling eras rather than model error alone.
