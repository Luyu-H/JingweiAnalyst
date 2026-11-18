## 1. Answer to the question

All four sources have a 200–1000 m nitrate-to-phosphate stoichiometry **below Redfield (16:1)**; none sits at or above it. The observed WOA23 headline mean N:P is **14.15** (mean of the cell-by-cell ratio) / **14.00** (ratio of the volume-weighted means), and the OLS relation NO3 = 13.04·PO4 + 1.69 µmol kg⁻¹ (R² = 0.908) — i.e. a slope well below 16 plus a positive nitrate intercept, which is the signature of nitrate excess at low phosphate (the denitrification/low-N:P signal). Only **6.2 %** of the observed 200–1000 m volume is below N:P = 12, but that low-N:P volume is strongly basin-specific (Pacific 8.1 %, Indian 10.3 %, with the Atlantic mean highest at 15.55 and the Pacific lowest at 13.57).

The three CMIP6 models split into two behaviours:

* **GFDL-ESM4 is the only model that reproduces the observed global stoichiometry** (headline 14.13, ratio-of-means 13.65, slope 11.37, b = 4.02, R² = 0.817) and it gets the Pacific essentially right (13.39 vs 13.57) — but it **over-fills the 14–18 band** (74.6 % vs 55.0 % of volume) and **fails on the low-N:P volume**: Pacific below-12 is about twice the observation (14.9 % vs 8.1 %), the Indian is almost empty (0.5 % vs 10.3 %), and the Atlantic is too low (15.07 vs 15.55).
* **CESM2 and MPI-ESM1-2-LR are far too low.** CESM2: 11.83 / 10.93, with a pathological linear relation (slope 4.89, intercept 10.84, R² = 0.238) and **36.7 %** of volume below 12 (Pacific 59.4 %, Indian 35.8 %); MPI: 12.28 / 11.77, slope 8.59, R² = 0.678, **22.5 %** below 12 (Pacific 40.0 %, Indian 15.0 %) and essentially no volume in 14–18 (0.5 %).

So **yes in sign, no in amplitude and no in basin fidelity**: all models produce a sub-Redfield, Pacific-depleted signal, but CESM2 and MPI exaggerate the low-N:P volume by a factor of ~4–6 and shift the whole distribution downward, while GFDL matches the mean yet concentrates the water in the Redfield band and misplaces the low-N:P volume (too much Pacific, almost none Indian). None reproduces the observed combination of a 14–18-dominated distribution *plus* a small (≈6 %) low-N:P tail confined to the Pacific and Indian.

---

## 2. Methodology actually executed

**Observations.** WOA23 `nitrate` and `phosphate` loaded with `load_climatology_data(dataset="woa23", variable=..., time_slice="annual", depth_range=[200, 1000])` — the annual objectively analysed all-period composite, on the 180 × 360 WOA grid, 23 standard levels from 200 m to 1000 m (187.5–1025 m bounds). Depth levels/bounds from `load_depth_levels(dataset="woa23", variable="nitrate", ...)`; cell areas from `load_areas(dataset="woa23")`. WOA has no partial bottom cells, so layer thicknesses were taken as the differences of its own level bounds (`woa_dz` = 25 m down to 500 m, then 37.5–50 m), giving a 200–1000 m layer volume of **2.7217 × 10¹⁷ m³** (Σ area·dz over cells where nitrate is valid). WOA nutrients are µmol kg⁻¹ as loaded and were used unchanged.

**Models.** `load_regular_data(dataset="cmip6_historical", variable in {no3, po4}, model, time_range=["1971-01","2000-12"], depth_range=[200, 1000])` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR. Monthly fields were reduced to a single 1971–2000 climatological mean per model/variable with `temporal_aggregate(freq="mean", op="mean")`. Volumes came from `load_volumes_for_cmip6(dataset="cmip6_historical", model, depth_range=[200,1000], variable="no3")` (first time slice taken when present), carrying `volume_source = "volcello published by cmip6_historical"` and `area_source = None` for all three models — i.e. the cell-thickness form of `volcello`, no `lev_bnds` differencing and no rescaling between models. Layer volumes: CESM2 2.7667 × 10¹⁷ m³, GFDL-ESM4 2.8655 × 10¹⁷ m³, MPI-ESM1-2-LR 2.7235 × 10¹⁷ m³. Levels retained in 200–1000 m on each model's own grid: **23 (WOA23), 20 (CESM2), 10 (GFDL-ESM4), 11 (MPI-ESM1-2-LR)**; no vertical interpolation was performed anywhere.

**N:P and layer means.** Cell-by-cell N:P from `compute_ratio(no3, po4)`; the depth mean at each horizontal cell and the global statistics both via `compute_basic_stat(op="volume_weighted_mean")` over depth + horizontal dimensions. Two different global quantities were reported: (a) the volume-weighted mean of the cell-by-cell ratio (headline) and (b) the ratio of the volume-weighted mean nitrate to the volume-weighted mean phosphate.

**Regression.** All valid cells of the 200–1000 m layer were stacked (level × latitude × longitude) and weighted equally; NO3 = a·PO4 + b fitted with `fit_multivariate_relation(target=NO3, predictors=[PO4], method="linear", validation="random")`, with the OLS coefficients and R² also recomputed directly by `np.polyfit` and reported (identical values). Model concentrations were converted from mol m⁻³ to µmol kg⁻¹ with `transform_unit(from_unit="mol m-3", to_unit="umol/kg", density=1.025)`, i.e. the prescribed constant ρ = 1025 kg m⁻³ (≈ ×10⁶/1025); WOA23 was taken as already µmol kg⁻¹.

**Volume fractions.** `compute_exceedance_fraction(..., dim=None, weighting="volume")`: for 14 < N:P < 18 the call used `threshold=2.0, comparison="abs_less", reference=16.0` (|N:P − 16| < 2); for below 12, `threshold=12.0, comparison="less", reference=0.0`. Missing cells were excluded from the denominator.

**Basins.** `load_basin_mask(depth_levels="0m")` returned the WOA23 basin mask on the 180 × 360 grid (values 1–56); `select_basin(..., include_marginal_seas=True, depth_dependent=False)` isolated Pacific, Atlantic and Indian in turn (marginal seas folded into their parent ocean, surface labels); the same volume-weighted headline mean was then recomputed within each basin, and the below-12 fraction recomputed within the Pacific and Indian. The Southern and Arctic Oceans were not reported.

**Comparison mask.** `build_comparison_mask(..., ocean_only=True)` was applied to the four *200–1000 m layer-mean maps regridded conservatively to the WOA23 1° grid*, and only its cell count was used for the like-for-like diagnosis: **38 167 common ocean cells**. Important caveat: the mask was **not** applied with `.where(mask)` to the native-grid fields before the regression, ratio or exceedance metrics — all headline numbers above are on each source's own grid, own valid cells and own volumes, and the mask enters only in this diagnostic role.

**Caveats stated.** (i) Constant ρ = 1025 kg m⁻³ instead of in-situ density introduces ~1 % error in the µmol kg⁻¹ concentrations and hence in the slope/intercept (N:P ratios are unaffected). (ii) WOA µmol L⁻¹ treated as µmol kg⁻¹ (<2 % effect). (iii) GFDL-ESM4 is the `gr` 1° archive grid with an already-remapped land–sea mask and only approximately conservative volumes; CESM2 (~100 km) and MPI (~250 km) differ in resolution, so narrow features are not comparable in amplitude. (iv) The three sources sample the 200–1000 m band with 10–23 levels, so the vertical weighting differs between them.

---

## 3. Key numerical results

| | WOA23 | CESM2 | GFDL-ESM4 | MPI-ESM1-2-LR |
|---|---|---|---|---|
| Levels in 200–1000 m | 23 | 20 | 10 | 11 |
| Layer volume (m³) | 2.722e17 | 2.767e17 | 2.866e17 | 2.724e17 |
| **Mean of cell-by-cell N:P** (headline) | **14.15** | **11.83** | **14.13** | **12.28** |
| Ratio of volume-weighted means | 14.00 | 10.93 | 13.65 | 11.77 |
| Slope a | 13.04 | 4.89 | 11.37 | 8.59 |
| Intercept b (µmol kg⁻¹) | 1.69 | 10.84 | 4.02 | 5.64 |
| R² | 0.908 | 0.238 | 0.817 | 0.678 |
| % volume 14 < N:P < 18 | 55.0 | 17.2 | 74.6 | 0.5 |
| % volume N:P < 12 (global) | 6.2 | 36.7 | 7.3 | 22.5 |
| Pacific mean | 13.57 | 10.58 | 13.39 | 11.29 |
| Atlantic mean | 15.55 | 13.55 | 15.07 | 13.15 |
| Indian mean | 13.73 | 11.18 | 14.72 | 12.60 |
| % < 12, Pacific | 8.1 | 59.4 | 14.9 | 40.0 |
| % < 12, Indian | 10.3 | 35.8 | 0.5 | 15.0 |
| Regression n (valid cells) | 859 059 | 1 558 515 | 394 259 | 330 410 |
| Position vs Redfield | below | below | below | below |

Global volume-weighted mean concentrations: WOA23 NO3 27.49, PO4 1.964 µmol kg⁻¹; CESM2 0.02184 / 0.001999 mol m⁻³; GFDL-ESM4 0.02508 / 0.001838 mol m⁻³; MPI 0.02089 / 0.001774 mol m⁻³.

---

## 4. Figures produced

**Figure 1** — 2 × 2 panel of NO3 (y) versus PO4 (x) in µmol kg⁻¹ for the 200–1000 m layer, panels ordered WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, with identical axis limits and ticks derived from the pooled 99.5th percentile across all four sources. Each panel shows 12 000 randomly sub-sampled cells (seed 0, alpha-blended, for display only), that source's own OLS fit line (orange) with a=, b=, R²= annotated in a text box, and the Redfield 16:1 line through the origin (dashed grey), so the two lines are visually distinguishable. The panel title carries the source name and R²; the figure footnote states that model concentrations were converted at ρ = 1025 kg m⁻³.

**Figure 2** — 2 × 2 panel of global maps of the 200–1000 m mean N:P minus 16 for the four sources, all on the common 1° WOA23 longitude–latitude frame (curvilinear model fields regridded conservatively with `regrid_horizontal(..., target="woa23", method="conservative")`, display only), longitude sorted to −180…180, latitude limited to −80…90, with an equator line. A single shared diverging `RdBu_r` colour scale is used, centred on zero via `TwoSlopeNorm` and clipped symmetrically at the 98th percentile of |anomaly| (minimum 1.0), masked land in grey; one common colour bar is drawn beside the panels and each title repeats the source's global headline mean.

---

## 5. Insights and patterns

1. **The observation is not Redfield at the cell level.** The mean of the ratio (14.15) is higher than the ratio of the means (14.00), and the regression slope is 13.04 with a +1.69 µmol kg⁻¹ intercept — the classic sub-Redfield/positive-intercept geometry of waters influenced by denitrification, which produces N:P below 16 where phosphate is low.
2. **Basin ordering in the observation is Atlantic > Indian > Pacific** (15.55, 13.73, 13.57); the low-N:P tail (<12) is 6.2 % globally but 8.1 % in the Pacific and 10.3 % in the Indian. All four sources and all three models preserve the Atlantic-highest ordering, so the basin *pattern* is robust even where amplitudes are not.
3. **GFDL-ESM4 matches the level but not the shape.** Its headline (14.13) is almost identical to WOA23 (14.15) and its Pacific mean is close, yet it is much more tightly clustered around Redfield (74.6 % in 14–18 vs 55.0 %) and its low-N:P water is concentrated in the Pacific with almost none in the Indian — the opposite of the observed Indian contribution.
4. **CESM2 and MPI-ESM1-2-LR are biased low by 2–3 N:P units and over-produce low-N:P water by roughly a factor of 4–6** (36.7 % and 22.5 % below 12 vs 6.2 %), concentrated in the Pacific (59.4 % and 40.0 %). CESM2's very low slope (4.89) and large intercept (10.84) with R² = 0.238 show that its NO3–PO4 relationship is not a simple linear Redfield-like stoichiometry at all, so the "low N:P" there is a mean-state nutrient bias rather than a denitrification-like offset.
5. **The 38 167-cell shared comparison mask** confirms that the four layer fields can be compared on the same 1° ocean domain; the differences above are therefore not caused by differing valid-area coverage, though the amplitudes of narrow features remain affected by the differing native resolutions and by GFDL-ESM4's archive-grid (`gr`) volumes.
