## Direct answer

The 200–1000 m nitrate-to-phosphate stoichiometry was computed for WOA23 and the three CMIP6 historical models using their native grids and volume weights for the mean/percentage diagnostics, and an equal-weight-per-cell OLS regression for the NO3 = a·PO4 + b fit. All four sources sit **below Redfield (16:1)** in their headline mean N:P. **GFDL-ESM4 is closest to WOA23** in both mean N:P (14.13 versus 14.15) and global low-N:P volume fraction (7.3% versus 6.2% below N:P = 12). **CESM2 and MPI-ESM1-2-LR strongly overestimate the global low-N:P volume**, and CESM2 in particular has a much weaker NO3–PO4 relationship. The basin-resolved diagnostics could not be produced because the WOA basin-mask reader returned zero valid mask cells; therefore the Pacific/Atlantic/Indian means and Pacific/Indian below-12 percentages are unavailable, not model results.

## Methodology as executed

- **Models:** CMIP6 historical monthly `no3` and `po4` from CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, restricted to 1971–2000.
- **Model climatology:** calendar-month-length weighted monthly climatology using `days_in_month`; the cell-by-cell N:P ratio was formed from the monthly climatological NO3 and PO4 fields and then averaged over the 12 months. The reported concentration means and regression use the days-in-month-weighted 1971–2000 annual mean concentrations.
- **WOA23:** annual objectively analysed `nitrate` and `phosphate` from the single all-time-period composite, used as distributed.
- **Depth layer:** model and WOA levels whose nominal depth lies within 200–1000 m, stacked into one set of cells.
- **Volume weights:** each model’s own `volcello`; GFDL-ESM4 `volcello` was taken from the 1971–2000 Omon time mean, while CESM2 and MPI-ESM1-2-LR used the Ofx fixed field. WOA23 volume was computed from WOA depth-bound thickness times 1° spherical cell area.
- **Units:** model `no3`/`po4` were converted from mol m⁻³ to µmol kg⁻¹ with a constant density of 1025 kg m⁻³: multiplication by 1.0×10⁶/1025 = 975.609756. WOA23 is already in µmol kg⁻¹. N:P and volume percentages are scale-invariant; the density assumption affects the µmol kg⁻¹ regression concentrations.
- **Weighting distinction:** headline mean N:P and exceedance percentages are volume-weighted; the OLS regression is deliberately unweighted per cell.
- **Regression:** ordinary least squares NO3 = a·PO4 + b over all valid 200–1000 m cells, each cell weighted equally.
- **Time mismatch:** WOA23 is an all-period objectively analysed composite, not a 1971–2000 normal; the models are 1971–2000 climatologies.

## Key numerical results

Volume-weighted diagnostics for the 200–1000 m layer. “Mean cell-by-cell N:P” is the headline value; “ratio of means” is the ratio of the volume-weighted mean NO3 to the volume-weighted mean PO4. OLS is unweighted per cell.

| Source | Mean cell-by-cell N:P | Ratio of means | OLS slope a | OLS intercept b (µmol kg⁻¹) | R² | Volume with 14 < N:P < 18 | Global volume with N:P < 12 | Total valid volume (m³) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| WOA23 | 14.15 | 14.00 | 13.04 | 1.69 | 0.908 | 54.9% (1.495×10¹⁷ m³) | 6.2% (1.689×10¹⁶ m³) | 2.723×10¹⁷ |
| CESM2 | 11.83 | 10.93 | 4.89 | 10.84 | 0.238 | 17.2% (4.748×10¹⁶ m³) | 36.7% (1.016×10¹⁷ m³) | 2.767×10¹⁷ |
| GFDL-ESM4 | 14.13 | 13.65 | 11.37 | 4.02 | 0.817 | 74.6% (2.137×10¹⁷ m³) | 7.3% (2.083×10¹⁶ m³) | 2.865×10¹⁷ |
| MPI-ESM1-2-LR | 12.28 | 11.78 | 8.62 | 5.60 | 0.680 | 0.5% (1.485×10¹⁵ m³) | 22.5% (6.212×10¹⁶ m³) | 2.761×10¹⁷ |

Volume sources: WOA23 from WOA depth bounds × 1° cell area; CESM2 and MPI-ESM1-2-LR from `volcello` (Ofx); GFDL-ESM4 from `volcello` (Omon 1971–2000 time mean).

## Basin results

The basin-resolved calculations were attempted using `basinmask_01.msk` and the WOA surface labels, with marginal seas folded into their parent ocean and Southern/Arctic Oceans excluded. However, the executed mask reader produced **zero valid basin-labelled cells** (`basin_raw` and `basin_parent_2d` had `n_valid = 0`). Consequently, for every source the Pacific, Atlantic and Indian mean N:P values and the Pacific/Indian below-12 percentages are **NaN / unavailable**. The global results above are unaffected because they do not use the basin mask. The basin distribution of the low-N:P signal therefore cannot be assessed from this run.

## Interpretation

- **All sources are below Redfield in the headline mean N:P.** WOA23 is 14.15, GFDL-ESM4 is 14.13, MPI-ESM1-2-LR is 12.28 and CESM2 is 11.83. The fitted OLS slopes are also all below 16, though the intercepts are positive.
- **GFDL-ESM4 best reproduces the observed global mean and low-N:P volume.** Its global N:P < 12 fraction is 7.3%, close to WOA23’s 6.2%, and its mean N:P is nearly identical to WOA23. It does, however, place a larger fraction of volume in the 14–18 band (74.6% versus 54.9%).
- **CESM2 and MPI-ESM1-2-LR overestimate the global low-N:P volume.** CESM2 has 36.7% of volume below N:P = 12 and MPI-ESM1-2-LR has 22.5%, compared with 6.2% in WOA23. Both also have very small fractions in the 14–18 band (CESM2 17.2%, MPI 0.5%), and their mean N:P is substantially below Redfield.
- **CESM2 has the weakest NO3–PO4 coupling in this layer**, with R² = 0.238 and a low slope of 4.89. GFDL-ESM4 and WOA23 have much stronger relationships (R² = 0.817 and 0.908 respectively).
- The cross-model N:P spread should be read partly as structural: fixed N:P per plankton type in GFDL-ESM4, prescribed Redfield organic matter in MPI-ESM1-2-LR, and variable C:P in CESM2. A feature not detailed in a model paper should not be read as absent from the model.
- Because the basin mask failed, the question of whether models reproduce the **basin distribution** of the observed low-N:P signal remains unanswered in this execution.

## Figures generated

**Figure 1** is a 4-panel scatter of NO3 versus PO4 in µmol kg⁻¹ over 200–1000 m, ordered WOA23, CESM2, GFDL-ESM4 and MPI-ESM1-2-LR. All panels share identical axis limits and ticks. Each panel shows its own OLS fit and the Redfield 16:1 line through the origin as a visually distinct line, and is annotated with the fitted slope `a`, intercept `b` and R². The scatter is subsampled to at most 20,000 points per panel for display, while the regression uses all valid cells.

**Figure 2** is a 2×2 panel of global maps of the 200–1000 m mean N:P minus 16 for the four sources. The native model fields were conservatively regridded to a common 1° longitude–latitude grid for display only, using area-weighted conservative binning. All four panels use one shared diverging colour scale centred on zero, with the scale clipped at the 98th percentile of the pooled absolute anomaly (approximately ±9.87). The map is diagnostic/display only; all numerical diagnostics remain on each source’s native grid with its own volumes.
