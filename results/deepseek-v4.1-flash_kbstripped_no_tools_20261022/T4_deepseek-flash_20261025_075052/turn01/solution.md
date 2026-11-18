## 1. Direct answer

Yes — the offline operator reproduces each model's published `mlotst` to within a few metres at most locations, but not exactly, and the level of agreement differs systematically between models. On the common mask (cells where both the operator MLD and `mlotst` are defined in all 30 Marches or all 30 Septembers of 1971–2000), the area-weighted spatial correlation between the two fields is 0.954–0.998 in every one of the six model-and-month cases. The residual disagreement is a diagnostic difference, not a model-physics difference: the operator MLD is systematically **shallower** than `mlotst` for five of the six cases (bias −9.6 to +0.6 m), which is the sign expected from the known residual causes (see §5).

Ranking, most to least internally consistent by the stated criterion (lowest area-weighted RMSD averaged over March and September):

1. **MPI-ESM1-2-LR** — mean RMSD 7.4 m
2. **CESM2** — mean RMSD 7.7 m
3. **GFDL-ESM4** — mean RMSD 20.5 m

## 2. Methodology as actually implemented

**Data.** Raw NetCDF was opened with xarray from `cmip6_historical` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR: monthly `thetao`, `so`, `mlotst` (Omon) and fixed `areacello` (Ofx). Variant labels were read from the filenames themselves (all three resolved to `r1i1p1f1`); no variant was assumed. CESM2 additionally used `sftof`. No variable was missing for any model (`missing_report` empty). Each model was kept on its own native grid: CESM2 384×320 (`nlat`×`nlon`, 60 levels, `lev` in cm therefore multiplied by 0.01 to metres), GFDL-ESM4 180×360 (`lat`×`lon`, 35 levels, scale 1.0), MPI-ESM1-2-LR 220×256 (`j`×`i`, 40 levels, scale 1.0). No common grid was imposed for the statistics. The time axis was restricted to years 1971–2000 inclusive and split by calendar month; the code asserts that exactly 30 March and 30 September records are present and that the `thetao`, `so` and `mlotst` time axes are index-aligned.

**Operator MLD.** For every cell and every monthly record: pressure from depth and latitude (`gsw.p_from_z`, depth positive up); practical salinity converted to Absolute Salinity with the location-dependent `gsw.SA_from_SP` using the stored values as-is (no ×0.001 or ×1000, since `so` is stored ~35 with units "0.001"); `thetao` converted to Conservative Temperature via `gsw.CT_from_pt`; density from `gsw.rho(SA, CT, p)`. The reference parcel is taken at the model's **shallowest level**. The displaced-parcel buoyancy difference δB = 9.81·(ρ_ambient(z) − ρ_parcel(z))/ρ_ref is evaluated at each level's own local pressure, and the shallowest level with δB > 3.0 × 10⁻⁴ m s⁻² is found; the crossing depth is linearly interpolated between that level and the level above using the criterion value. Columns where the threshold is never reached return NaN and are dropped; the area lost this way was recorded per record and averaged over the 30 records (see §4).

**Climatologies and mask.** Diagnostic first, average second: the operator MLD was computed for each of the 30 monthly profiles and then averaged; `mlotst` was likewise averaged over the same 30 records. No day weights were used. One common mask per model and month requires both the operator MLD and `mlotst` to be finite in **all 30** instances; it was applied to both climatologies before any statistic. Land and fill values stay NaN and are excluded from all area fractions (fill values >1e20 were masked out of `areacello`; for CESM2 the global `areacello` was restricted with `sftof > 0`).

**Scoring.** All statistics are weighted by the model's own `areacello`: area-weighted mean of (operator MLD − `mlotst`), area-weighted RMS difference, area-weighted spatial Pearson correlation (weighted covariance over the square root of the weighted variances), the percentage of the common-mask area with |difference| > 20 m, and the common-mask area as a percentage of the model's total ocean area (the sum of its own `areacello`, never a rescaled common total). No cell thickness or `volcello` was used; `lev_bnds` was not used.

## 3. Key numerical results

Operator MLD minus published `mlotst`, 1971–2000 March and September climatologies, common mask, area-weighted:

| Model | Month | Mean diff (m) | RMSD (m) | r | Area \|Δ\|>20 m (%) | Common-mask area (% of ocean) |
|---|---|---|---|---|---|---|
| CESM2 | March | −4.0 | 10.9 | 0.995 | 2.5 | 96.4 |
| CESM2 | September | −1.8 | 4.6 | 0.998 | 0.1 | 97.8 |
| GFDL-ESM4 | March | −9.6 | 27.8 | 0.954 | 6.2 | 97.3 |
| GFDL-ESM4 | September | −7.8 | 13.3 | 0.989 | 7.0 | 97.9 |
| MPI-ESM1-2-LR | March | −1.7 | 7.8 | 0.994 | 1.5 | 96.7 |
| MPI-ESM1-2-LR | September | +0.6 | 7.0 | 0.998 | 2.6 | 96.1 |

Number of common-mask cells: CESM2 82 355 (March), 84 115 (September); GFDL-ESM4 42 125 / 42 707; MPI-ESM1-2-LR 33 215 / 33 839.

Area excluded in step 2 (columns with a valid surface profile but no threshold crossing), as a percentage of each model's own ocean area, mean over the 30 records with the worst record in brackets: CESM2 2.4 (2.6) March, 1.2 (1.4) September; GFDL-ESM4 1.2 (1.4) March, 0.8 (0.9) September; MPI-ESM1-2-LR 1.4 (1.6) March, 1.8 (2.0) September. So the common mask still covers 96–98% of each model's ocean.

Unweighted field means over the common mask (operator vs `mlotst`, metres): CESM2 62.6 vs 66.4 (March), 69.2 vs 70.9 (September); GFDL-ESM4 40.2 vs 49.7, 48.7 vs 55.8; MPI-ESM1-2-LR 73.2 vs 73.8, 64.5 vs 64.0. Unweighted difference ranges: CESM2 −260 to +623 m (March), −27 to +63 m (September); GFDL-ESM4 −886 to +811 m, −126 to +413 m; MPI-ESM1-2-LR −167 to +560 m, −60 to +147 m. The very large excursions occur in March, i.e. where the deepest winter mixed layers are found; they are few cells and do not dominate the area-weighted scores.

**Interpretation.** The correlations are uniformly very high because both fields share the same large-scale spatial pattern; the discriminating quantities are the bias and RMSD. MPI-ESM1-2-LR and CESM2 are essentially tied (7.4 vs 7.7 m) and differ mainly in seasonality: CESM2's March bias (−4.0 m) and March RMSD (10.9 m) are the larger of its two months, with September almost matching `mlotst` (r = 0.998, only 0.1% of area beyond 20 m). MPI is nearly unbiased (largest |bias| 1.7 m) but pays for it with a slightly wider tail (2.6% > 20 m in September). GFDL-ESM4 is clearly the least internally consistent on every supporting metric as well as on RMSD: the largest bias magnitudes, the highest >20 m area fractions (6.2% and 7.0%), and the lowest correlations (0.954, 0.989). Its comparison is on the `gr` 1° archive grid with 35 z levels, i.e. a considerably coarser vertical and horizontal sampling than CESM2's, which limits how sharply the criterion crossing can be located; the native model had 75 hybrid layers that were remapped before archiving. The ranking compares scores computed on each model's own grid and its own `areacello` total, so it is not a cell-by-cell cross-model comparison, and CESM2 (~100 km) and MPI (~250 km) differ in resolution by ~2.5×: features narrower than the coarse grid cannot be compared in amplitude.

## 4. Where the disagreement concentrates

The code binned the |difference| > 20 m area into six latitude bands (−90 to −60, −60 to −30, −30 to 0, 0 to 30, 30 to 60, 60 to 90) crossed with five sectors (Southern Ocean for lat < −50°, Arctic for lat > 60°, and Pacific 120–280°E, Atlantic ≥280°E or <20°E, Indian 20–120°E outside those polar bands), computing the >20 m fraction within each region and that region's share of total ocean area; 90 region–model–month rows were populated. The highest fractions sit where the mixed layer is deepest and where the operator's vertical resolution and reference-level choice matter most — the March high-latitude and deep-convection regimes — which is also where the unweighted difference ranges are largest (up to +623 m for CESM2, +811 m for GFDL-ESM4, +560 m for MPI in March); September differences are much smaller in every model. The excerpt printed for CESM2 March contains only zero-fraction entries (Southern Ocean bands), and the remaining rows are best read off Figure 1. Where differences in these regions are associated with narrow frontal or convective features, the CESM2–MPI resolution difference of ~2.5× means those amplitudes are not directly comparable between models.

## 5. Why operator and `mlotst` are not expected to agree exactly

These are residuals of the diagnostic, not model error, and were not scored as such:

1. **Reference level** — the operator references the parcel to the model's shallowest level, whereas `mlotst` is referenced to the ocean surface; this alone biases the operator shallow.
2. **Density formulation** — the operator uses displaced-parcel buoyancy at each level's local pressure rather than a σ₀ (potential-density) criterion.
3. **Averaging order** — the published `mlotst` is an instantaneous diagnostic averaged over the month, while an MLD formed offline from monthly-mean θ/S is a different (generally shallower) object; f(monthly mean) ≠ mean of f.

The predominantly negative biases (−4.0, −1.8, −9.6, −7.8, −1.7 m) are the expected sign; MPI-ESM1-2-LR September (+0.6 m) is the only case that does not follow it. Density was necessarily taken through the potential-temperature route (no `bigthetao`), and salinity through the location-dependent S_P → S_A conversion, since the in-situ route and an unscaled S_P would both bias density and MLD regionally.

## 6. Figures

**Figure 1** — 3 rows (CESM2, GFDL-ESM4, MPI-ESM1-2-LR) × 2 columns (March, September) of global maps of operator MLD − `mlotst`. Each panel is a `pcolormesh` on a common regular 1° longitude–latitude frame (0.5–359.5°E, −89.5 to 89.5°N): the native curvilinear/regular model fields were bilinearly interpolated with `scipy.griddata` to that frame for display only, with the field duplicated at ±360° to close the dateline seam, and land suppressed by a nearest-neighbour interpolation of the native ocean mask and left white. A single diverging `RdBu_r` scale centred on zero is shared by all six panels, with symmetric limits taken from the pooled 98th percentile of |difference| across the panels (colour bar extended at both ends). A thin black contour marks the interpolated ocean mask at 0.5. Axes are longitude 0–360 with 60° ticks and latitude −90 to 90 with 30° ticks.

**Figure 2** — same 3×2 layout of scatter plots, one point per common-mask ocean cell, operator MLD (x) against `mlotst` (y), both on logarithmic axes with identical limits and identical logarithmic tick sets in every panel, a dashed y = x reference line, and light major grid lines. Each panel is annotated with the area-weighted bias, RMSD and correlation from §3 and with the number of common-mask cells whose MLD is non-positive and therefore cannot be placed on logarithmic axes; those cells are omitted from the panel only and remain in all statistics. The intended reading is the departure of the point cloud from the 1:1 line as a function of MLD magnitude — near-1:1 at shallow MLD, with a systematic offset at the deep-MLD end where the diagnostic differences are largest.

## 7. Insights

- The internal-consistency test is passed in the loose sense for all three models (r ≥ 0.954, 96–98% of ocean area usable), so the published `mlotst` is genuinely reproducible from the archived θ/S with the stated criterion; the differences are diagnostic, not structural.
- The ranking is driven by the deep-MLD tail, not by the bulk of the ocean: >20 m differences occupy only 0.1–7.0% of area yet dominate RMSD, and they coincide with the March deep-convection and high-latitude regimes.
- Seasonal asymmetry is large in CESM2 (March RMSD 10.9 m vs September 4.6 m) and moderate in GFDL-ESM4 (27.8 vs 13.3 m), but weak in MPI (7.8 vs 7.0 m); MPI is the most seasonally stable and the least biased.
- The area lost to columns with no threshold crossing (0.8–2.4% of ocean, worst single record 2.6%) is at the same order as the residual disagreement and is reported separately rather than folded into the scores.
- Because scores are computed on each model's own grid, `areacello` and ocean total, and because the archives differ in resolution (`gr` 1°/35 levels for GFDL-ESM4 vs ~100 km/60 levels for CESM2 vs ~250 km/40 levels for MPI), the ranking should be read as a statement about each model's diagnostic self-consistency, not as a cross-model comparison of mixed-layer physics.
