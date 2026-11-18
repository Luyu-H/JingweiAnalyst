## 1. Direct answer

Yes — the test can be run, and the answer is that a density-criterion MLD computed offline from each model's own `thetao`/`so` largely reproduces the model's published `mlotst`, but **not exactly**: for all six model-and-month combinations the operator MLD is systematically *shallower* than `mlotst` (all biases negative), with weighted RMS differences ranging from 5.1 m to 31.7 m and spatial correlations of 0.938–0.998. Internal consistency, ranked by the mean of the March and September area-weighted RMS difference (smallest = most consistent), is:

1. **CESM2** — mean RMS 6.3 m (7.5 m March, 5.1 m September), r = 0.998/0.998, bias −4.5/−2.6 m
2. **MPI-ESM1-2-LR** — mean RMS 19.2 m (10.8 m March, 27.5 m September), r = 0.995/0.956, bias −0.9/−8.0 m
3. **GFDL-ESM4** — mean RMS 25.5 m (31.7 m March, 19.2 m September), r = 0.938/0.985, bias −11.4/−7.5 m

The residual disagreements are definitional (reference level, displaced-parcel density at local pressure vs. σ0, and MLD-of-monthly-means vs. instantaneous MLD), not model physics errors.

## 2. Methodology as actually executed

**Data and period.** Files were read from `/data1/lorraine-20419/benchmark_cmip6/historical/<model>/<var>/`. For each variable the variant present in the filenames was identified and the most common one retained (all variables resolved to `r1i1p1f1`). `thetao` and `so` (Omon) were opened with `xr.open_mfdataset(combine="by_coords", use_cftime=True, decode_times=True, chunks={"time":1})`; `mlotst` (Omon) and the Ofx fields `areacello` (and `sftof` for CESM2) likewise. The time axis was restricted to `dt.year` 1971–2000 inclusive and the code verified 360 months; `so` and `mlotst` were aligned to the `thetao` time axis with `.sel(time=...)`. March and September indices were selected by `dt.month == 3` and `== 9`, each checked to contain exactly 30 months.

**Grids and levels (no regridding at this stage).** CESM2: curvilinear `nlat`×`nlon` = 384×320, `lev` scaled by 0.01 (60 levels, 5–5375 m), masked by `sftof > 0`. GFDL-ESM4: regular `lat`×`lon` = 180×360, `lev` scale 1.0 (35 levels, 2.5–6500 m). MPI-ESM1-2-LR: curvilinear `j`×`i` = 220×256, `lev` scale 1.0 (40 levels, 6–5720 m). Level order was sorted ascending (`argsort`) and converted to metres below the surface. `areacello` was transposed to (y, x); values with magnitude > 1e10 were set to NaN; for CESM2 cells with `sftof == 0` were additionally set to NaN. Total ocean area = `nansum(areacello)` over the surviving wet cells (CESM2 3.605×10¹⁴ m², GFDL-ESM4 1.453×10¹⁴ m², MPI-ESM1-2-LR 7.846×10¹³ m²).

**Operator MLD.** For every monthly profile and each horizontal column: depth z = −depth_m, p = `gsw.p_from_z(z, lat)`, longitude wrapped to [−180, 180]; `so` treated as practical salinity (units "0.001") and converted with `gsw.SA_from_SP(so, p, lon, lat)`; `thetao` treated as potential temperature and converted along the potential-temperature route with `gsw.CT_from_pt(SA, thetao)`; potential density `gsw.rho(SA, CT, p_ref)` referenced to the pressure of the *shallowest level in each column* (`p_ref = p[0]`), i.e. effectively the surface, not 10 m. The first level where ρ exceeds ρ(surface) + 0.03 kg m⁻³ was located with `argmax` on the boolean condition, and the crossing depth linearly interpolated between that level and the one above (if the threshold is already exceeded at the first level, that level's depth is used). Columns never reaching the threshold, or with missing surface T/S, were left NaN. Fill values were additionally scrubbed with `|x| > 1e10 → NaN` before use.

**Climatology and shared mask.** For each model and month the 30 operator-MLD fields and the 30 `mlotst` fields were stacked; `build_comparison_mask` returns True only where **both** fields are finite in **all 30** instances. The mask was applied to both fields, then the 30 values were averaged with `nanmean` (each month's MLD computed first, then averaged — never averaging T/S first). No day weights were used.

**Metrics.** All statistics used `areacello` as weights on the native grid: weighted bias (operator − `mlotst`), weighted RMS, weighted spatial Pearson correlation (weighted covariance over weighted variances), and the area fraction with |difference| > 20 m expressed as a share of the shared-mask weighted area; the shared-mask area is also reported as a percentage of the model's total ocean area. Cells with non-finite `areacello` (e.g. decoded fill values) drop out of the weighted sums, so the number of cells entering the metrics can be smaller than the mask cell count (GFDL-ESM4: 23 649/24 777 of 42 122/43 259 mask cells; MPI: 18 905/20 783 of 33 260/35 177; CESM2: all mask cells used).

**Display regridding.** Only for Figure 1 (and the cell-by-cell map view) was a common 1° longitude–latitude grid used: source cell centres were binned into 1°×1° target cells and averaged with `areacello` as weight (labelled in the code as a conservative area-weighted cell-centre regrid). Area means, counts and all Step-4 metrics remain on each model's native grid.

## 3. Key numerical results

| Model | Month | Bias (m) | RMS (m) | r | \|Δ\|>20 m (%) | Shared area (% of ocean) | Mask cells |
|---|---|---|---|---|---|---|---|
| CESM2 | March | −4.5 | 7.5 | 0.998 | 2.7 | 96.4 | 82 360 |
| CESM2 | September | −2.6 | 5.1 | 0.998 | 0.2 | 97.9 | 84 146 |
| GFDL-ESM4 | March | −11.4 | 31.7 | 0.938 | 10.6 | 94.1 | 42 122 |
| GFDL-ESM4 | September | −7.5 | 19.2 | 0.985 | 15.0 | 97.3 | 43 259 |
| MPI-ESM1-2-LR | March | −0.9 | 10.8 | 0.995 | 3.7 | 91.4 | 33 260 |
| MPI-ESM1-2-LR | September | −8.0 | 27.5 | 0.956 | 7.5 | 94.1 | 35 177 |

Shared-mask areas: CESM2 3.476×10¹⁴ (March) and 3.528×10¹⁴ m² (September); GFDL-ESM4 1.367×10¹⁴/1.414×10¹⁴ m²; MPI-ESM1-2-LR 7.168×10¹³/7.385×10¹³ m². The 91–98 % coverage confirms that the comparison is essentially global, with the shortfall coming from columns where the 0.03 kg m⁻³ threshold is never reached or where a field is missing in at least one of the 30 months.

Local extremes are large in the deep-convection cells (GFDL-ESM4 March difference spans −626 to +814 m; MPI September reaches −685 m; CESM2 March −372 to +162 m), i.e. a small density error in a weakly stratified, deeply mixed column translates into a very large depth error. Away from those cells the fields track each other almost one-to-one (r ≥ 0.938 everywhere).

**Where the disagreement concentrates** (area-weighted exceedance fractions from the latitude-band and basin diagnostics):
- **CESM2** — March: North Atlantic 29.5 %, band 30–40°N 26.6 %, 40–50°N 8.1 %, 50–60°N 6.6 %; essentially nothing in the Southern Ocean or tropics. September: near-zero everywhere (largest band value 0.7 % at 10–20°N).
- **GFDL-ESM4** — March: North Atlantic 57.0 %, bands 30–40°N 60.0 %, 60–70°N 42.1 %, 50–60°N 33.2 %, 40–50°N 29.1 %, North Pacific 19.8 %, and 9.3 % in 80–70°S. September: Southern Ocean 20.1 %, bands −40 to −30°S 41.5 %, −50 to −40°S 31.8 %, −80 to −70°S 23.4 %, −60 to −50°S 13.8 %.
- **MPI-ESM1-2-LR** — March: North Atlantic 14.5 %, bands 60–70°N 27.2 %, 50–60°N 21.3 %, 30–40°N 8.4 %. September: Southern Ocean 13.4 %, bands −80 to −70°S 50.5 %, −70 to −60°S 24.8 %.

So the exceedance is a winter deep-mixing signal: the North Atlantic subpolar/Nordic seas and the 30–40°N band in March, and the Southern Ocean/Antarctic bands in September (GFDL-ESM4 and MPI).

## 4. What the figures show

**Figure 1** (`visualization...turn1_1.png`) is the requested 3×2 panel of global maps of operator MLD minus `mlotst`: rows CESM2, GFDL-ESM4, MPI-ESM1-2-LR; columns March, September. Each model's curvilinear fields were placed on the common 1° longitude–latitude frame (0.5–359.5°E, −89.5–89.5°N) with the conservative area-weighted cell-centre regrid described above, purely for display. All six panels share one diverging `RdBu_r` colour scale centred on zero and extending over ±36.6 m (the 98th percentile of |difference| pooled over the six regridded maps), with a single vertical colourbar labelled "MLD difference (m)", flat shading, gridlines, 60° longitude ticks and 30° latitude ticks, and equal aspect.

**Figure 2** (`visualization...turn1_2.png`) is the matching 3×2 panel of scatter plots, operator MLD (x) against published `mlotst` (y), one point per native shared-mask ocean cell (82 360/84 146 for CESM2; 23 649/24 777 for GFDL-ESM4; 18 905/20 783 for MPI-ESM1-2-LR), using only cells where both climatologies are positive. Axes are logarithmic with identical limits (0.1–10 000 m) and identical decade tick labels in every panel, a dashed y = x reference line, dotted gridlines at all minor and major ticks, and each panel's title annotated with the area-weighted bias, RMS difference and correlation from the native-grid metrics table. Both figures carry a suptitle stating the 1971–2000 baseline, the display regridding method, and the common scatter/colour limits.

## 5. Insights and patterns

- **Sign is consistent and expected.** Every bias is negative: the operator MLD is shallower than `mlotst` in all six combinations. This is the signature of the definitional differences flagged in the analysis plan — the reference is the shallowest model level (effectively the surface) rather than the 10 m often quoted, the criterion is applied to a displaced-parcel density at local pressure rather than a σ0 difference, and the operator is a mean of 30 monthly MLDs rather than a climatological MLD (or a monthly-mean MLD from instantaneous fields). The operator is never expected to reproduce `mlotst` exactly.
- **CESM2 is by far the most internally consistent** (r ≈ 0.998, RMS 5–7 m, exceedance ≤ 2.7 %), and its disagreement is confined almost entirely to the March North Atlantic deep-convection region, where the mixed layer is deep and a 0.03 kg m⁻³ threshold crossing is sensitive to the density formulation.
- **GFDL-ESM4 is the least consistent**, with the largest bias (−11.4 m in March) and RMS (31.7 m in March), and exceedance concentrated in the March North Atlantic (57 % of that basin's area) and in the September Southern Ocean (20 %). The GFDL-ESM4 profile has only 35 levels with coarse near-surface spacing (2.5, 10, 20, 32.5 m …), so the interpolated crossing depth is more strongly quantised than in the finer CESM2 (60 levels) and MPI (40 levels) columns; this is a plausible contributor to the larger RMS in addition to the definitional terms.
- **MPI-ESM1-2-LR shows a strong seasonal asymmetry**: nearly unbiased in March (−0.9 m, RMS 10.8 m, r 0.995) but a −8.0 m bias and 27.5 m RMS in September, driven by the Antarctic/Southern Ocean bands (−80 to −70°S: 50.5 % of area exceeding 20 m).
- **Correlations are high everywhere (≥ 0.938) and RMS is concentrated in a small area fraction** (0.2–15 %), i.e. the two fields agree in pattern almost everywhere and disagree in magnitude only where the mixed layer is deep or the vertical grid is coarse. Any use of `mlotst` as a validation target should therefore expect a shallow bias of a few metres in the mean, with local errors of tens to hundreds of metres in winter deep-mixing regions.
- **Caveat on the area percentages:** the shared-area percentage is relative to `nansum(areacello)` after the code's own masking (`|area| > 1e10 → NaN`, plus the `sftof` wet-cell mask for CESM2). The resulting total ocean areas differ substantially between models (3.605×10¹⁴, 1.453×10¹⁴, 7.846×10¹³ m²), so the "% of ocean" column should be compared across models with that denominator definition in mind; the bias, RMS, correlation and exceedance fractions are unaffected because they are normalised by the same weighted shared area.
