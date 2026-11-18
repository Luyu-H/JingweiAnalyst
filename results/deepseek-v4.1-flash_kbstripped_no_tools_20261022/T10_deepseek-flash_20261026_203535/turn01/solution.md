## 1. What the run actually did (methodology, as executed)

**Reference (WOA23).** The three 1971–2000 climate-normal objectively-analysed fields were opened from `woa23_*_decav71A0_**00_01.nc`: oxygen (`o_an`, µmol kg⁻¹), apparent oxygen utilization (`A_an`, µmol kg⁻¹) and percent oxygen saturation (`O_an`, %), then aligned onto one depth–lat–lon grid (102 levels, 1°, −179.5…179.5). Geometry was rebuilt by the code, not read from the files: cell areas from a spherical formula (R = 6 371 000 m, mean Δlat = Δlon = 1°), layer thickness from midpoints between the 102 standard levels with the top bound at 0 m and the bottom bound extrapolated to 5550 m, and an ocean-level mask derived from `landsea_01.msk` bottom levels (level kept where lev ≤ bottom level). Basin labels came from `basinmask_01.msk`, reduced to the surface level (`depth_dependent=False`) and **folded into parents**: codes 1/2/3 kept as Atlantic/Pacific/Indian, 6–8 → Atlantic, 9–10 → Indian, other shelves assigned by their centroid latitude/longitude (>60°N → Arctic, <−50° → Southern, 110–290°E → Pacific, otherwise Atlantic); Arctic (4) and Southern (5) were then excluded from all basin totals.

**CMIP6.** Monthly `o2`, `thetao`, `so` were opened per model/variable (one ensemble variant chosen = the variant with the most files) and averaged over calendar years **1971–2000** (historical: GFDL-ESM4, MPI-ESM1-2-LR) and **2071–2100** (ssp245: CESM2, GFDL-ESM4, MPI-ESM1-2-LR) by simple time mean. Density was computed with a UNESCO-1983 in-situ polynomial in T, S and hydrostatic pressure P = 1.013 + 1.0055·z (CESM2 levels scaled ×0.01). `volcello` (GFDL time-varying, averaged over the same window), `areacello` and `deptho` were taken on each model's own grid; column thickness = volcello/areacello. **Oxygen unit conversion used the executed expression `o2 × 1000 / ρ`, not `o2 × 1e6 / ρ`** — see the caveat below.

**Threshold operation.** Volume was accumulated over cells where the field is strictly below the threshold, with field and weights jointly masked (NaN/zero-volume cells excluded from both numerator and denominator). Inventories used native mol m⁻³ × volcello with the same masking. AOU and percent saturation were volume-weighted over WOA23 cells with o2 < 60 µmol kg⁻¹ only.

## 2. Key numbers

**Threshold volumes (10¹⁵ m³)**

| source | < 20 | < 60 |
|---|---|---|
| WOA23 1971–2000 | 14.2 | 80.5 |
| GFDL-ESM4 historical | 1 330 | 1 330 |
| MPI-ESM1-2-LR historical | 1 360 | 1 360 |
| CESM2 ssp245 | 1 330 | 1 330 |
| GFDL-ESM4 ssp245 | 1 330 | 1 330 |
| MPI-ESM1-2-LR ssp245 | 1 360 | 1 360 |

**Maximum <60 µmol kg⁻¹ column thickness.** WOA23 3125.0 m at 7.50°, 121.50°; GFDL-ESM4 historical 6486.5 m at 31.50°, 142.50° (ssp245 6486.6 m, same point); MPI-ESM1-2-LR historical 6020.0 m at 27.73°, 301.58° (ssp245 identical); CESM2 ssp245 5500.0 m at −41.82°, 212.56°.

**Dissolved-oxygen inventory (mol).** GFDL-ESM4 historical 2.20×10¹⁷; MPI-ESM1-2-LR historical 2.25×10¹⁷; CESM2 ssp245 2.02×10¹⁷; GFDL-ESM4 ssp245 2.15×10¹⁷; MPI-ESM1-2-LR ssp245 2.20×10¹⁷.

**Basin <60 µmol kg⁻¹ volumes (10¹⁵ m³).** WOA23: Pacific 67.5, Indian 11.6, Atlantic 1.33. GFDL-ESM4 historical: 630, 227, 290. MPI-ESM1-2-LR historical: 649, 227, 291.

**Scenario change (historical → ssp245, positive = expansion).** GFDL-ESM4 +0.00 % at both 20 and 60 (+7.7×10⁻⁴ %); MPI-ESM1-2-LR 0.00 % at both.

**Model–WOA differences (own-grid integrals, 10¹⁵ m³).** GFDL-ESM4 Δ₂₀ = 1321, Δ₆₀ = 1254; MPI-ESM1-2-LR Δ₂₀ = 1345, Δ₆₀ = 1279. Closer historical model (sum of absolute differences 2575 vs 2624): **GFDL-ESM4**.

**WOA23 OMZ water-mass character (<60 µmol kg⁻¹, volume-weighted).** AOU = 271.3 µmol kg⁻¹; saturation = 11.3 %. That is near-complete oxygen utilisation, i.e. the diagnosed zone is dominated by **respiratory consumption** (high AOU, saturation barely above one-tenth of equilibrium), with the persistence of such a large volume requiring sluggish ventilation to let the AOU accumulate.

## 3. Important caveat on the model numbers

The executed conversion was `o2 × 1000/ρ` where mol m⁻³ → µmol kg⁻¹ requires `×1e6/ρ`. As a consequence the "converted" model fields span −0.026 to 0.49, i.e. they remain on the ~mol m⁻³ scale and every wet cell falls below both 20 and 60 µmol kg⁻¹. This is why (i) the <20 and <60 volumes are identical for every model, (ii) each model value equals its own total ocean volume (≈1.33–1.36×10¹⁸ m³, and GFDL/MPI totals differ only by grid/volume differences), (iii) model maximum thicknesses equal the full water column (6020–6487 m), (iv) basin totals are whole-basin volumes, and (v) the scenario changes are numerical zeros. **Those five families of numbers are therefore degenerate and not physical OMZ estimates.** The WOA23 columns (threshold volumes, AOU/saturation, basin volumes, thickness maximum), the inventories, and the grid/volume bookkeeping are unaffected by the conversion and stand as computed. Because the model thresholds never bite, the "closer model" comparison is dominated by the difference in total ocean volume, not by OMZ skill, and must not be read as a model-quality ranking.

Also note that the code never applied a common `build_comparison_mask`; every total is an own-grid integral, so model-vs-WOA differences mix basin-grid, mask and volume differences.

## 4. Figures

- **Figure 1** (one row, three panels): column thickness of water below 60 µmol kg⁻¹ for GFDL-ESM4 historical, MPI-ESM1-2-LR historical and WOA23, on one shared sequential (`viridis`) scale from 0 to the pooled 98th percentile of the three fields, with a single colour bar in metres. Given the conversion caveat, the two model panels show essentially the wet-column thickness (bathymetry-shaped) rather than an oxygen minimum zone, whereas the WOA23 panel shows the genuine observed sub-60 µmol kg⁻¹ layer.
- **Figure 2** (one row, three panels): the two change panels (ssp245 2071–2100 minus historical 1971–2000 thickness at the 60 threshold) for GFDL-ESM4 and MPI-ESM1-2-LR on a shared zero-centred diverging (`RdBu_r`) scale, and CESM2's ssp245 thickness on its own sequential scale, labelled explicitly as ssp245-only because no CESM2 historical `o2`/member exists to form a difference.

## 5. Insights

1. Only the WOA23 branch is currently interpretable: a global sub-60 µmol kg⁻¹ volume of 80.5×10¹⁵ m³, of which 84 % is Pacific (67.5), 14 % Indian (11.6) and <2 % Atlantic (1.33), with a maximum column of 3125 m in the western equatorial Pacific (7.5°N, 121.5°E) and a volume-weighted AOU/saturation signature (271.3 µmol kg⁻¹, 11.3 %) that points to respiration-supported oxygen depletion.
2. The oxygen inventories (2.0–2.25×10¹⁷ mol) are the only model quantities in the run that are dimensionally consistent; they decrease from historical to ssp245 in both models (GFDL 2.20→2.15, MPI 2.25→2.20 ×10¹⁷ mol), a ~2 % loss, which is the direction expected under warming but is several times smaller than the degenerate 0.00 % OMZ-volume change would suggest.
3. Both models would agree in sign only trivially here; no statement about the size of the OMZ expansion can be carried out of this run. Re-running the model thresholds with `o2 × 1e6/ρ` (and optionally a common comparison mask for the model-vs-WOA step) is required before the "which model is closer" and "do the models agree" questions can be answered scientifically.
