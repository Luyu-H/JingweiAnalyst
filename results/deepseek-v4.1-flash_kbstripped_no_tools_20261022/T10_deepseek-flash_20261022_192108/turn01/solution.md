## 1. Direct answer

Six source/experiment combinations were evaluated on their own native grids with their own `volcello`, using period‑mean oxygen converted to µmol kg⁻¹ (WOA23 needs no conversion; CMIP6 `o2` in mol m⁻³ was divided by an EOS‑80 σ₀ computed from the period‑mean `thetao`/`so`). The WOA23 1971–2000 climatology gives an OMZ volume of **14.2 ×10¹⁵ m³ below 20 µmol kg⁻¹** and **80.5 ×10¹⁵ m³ below 60 µmol kg⁻¹**. Both historical CMIP6 models overestimate these volumes: GFDL‑ESM4 gives **42.9** and **104 ×10¹⁵ m³**, MPI‑ESM1‑2‑LR gives **91.6** and **212 ×10¹⁵ m³**; GFDL‑ESM4 is the better of the two at both thresholds. For the ssp245 2071–2100 window the models do **not** agree on the sign of the change at 60 µmol kg⁻¹: GFDL‑ESM4 expands by **+4.59 %** while MPI‑ESM1‑2‑LR contracts slightly by **−0.71 %**; at 20 µmol kg⁻¹ both contract (**−4.07 %** and **−3.03 %**). CESM2 contributes only an ssp245 snapshot (r10i1p1f1) and no change can be formed — and its volume diagnostic is additionally unusable (see §6).

## 2. Methodology as executed

**Observations.** WOA23 `oxygen`, `percent_oxygen_saturation` and `apparent_oxygen_utilization` were read from the `1971-2000_climate_normal` `decav71A0` files (variable names ending `_an`), on a 1° grid (lat −89.5…89.5, lon −179.5…179.5) with **102 depth levels, 0–5500 m, in metres**. The `oxygen` attribute units are `micromoles_per_kilogram`, so no density conversion was applied to WOA. Cell volumes were built from spherical geometry (R = 6 371 000 m, sin(lat)-weighted areas) with layer thicknesses from midpoints between standard levels (top bound 0, bottom 5550 m).

**WOA masks.** `landsea_01.msk` and `basinmask_01.msk` were parsed as plain tables, the 1° index rebuilt from the lat/lon columns, and the third column used as the surface value (`depth_dependent = False`, one label applied through the whole column).

**CMIP6.** GFDL‑ESM4 (180×360, 35 levels) and MPI‑ESM1‑2‑LR (220×256, 40 levels, `j`/`i`) historical 1971‑01–2000‑12; CESM2 (384×320, 60 levels, `nlat`/`nlon`, r10i1p1f1), GFDL‑ESM4 and MPI‑ESM1‑2‑LR ssp245 2071‑01–2100‑12 (both r1i1p1f1). No CESM2 historical `o2` exists and none was substituted.

**Weighting.** Monthly fields were combined into a single period mean by **month-length weighting**: each file was masked to the target years, multiplied by the number of days to the next time step (last month's length repeated), summed in chunks, and divided by the total weight — i.e. average-then-diagnose, not threshold-then-average, and not an equal-weight per calendar month mean.

**Units.** `o2` in mol m⁻³ → µmol kg⁻¹ via `o2 × 1e6/ρ`, with ρ the EOS‑80 potential density at p = 0 (`ref_pressure=None`) from the period-mean `thetao`/`so`. The recorded `density_assumption` for every CMIP6 combination is: *"monthly-mean thetao/so -> EOS-80 sigma0 (ref_pressure=None); umol/kg = mol/m3 * 1e6/rho"*. The inventory is the exception: it stays native mol m⁻³, multiplied by `volcello`.

**OMZ volume and thickness.** Masks `o2 < 20` and `o2 < 60` µmol kg⁻¹ were formed per level; volume is `volcello` summed over masked cells (cell counts reported). Thickness per column is the sum of `volcello/areacello` over levels below 60 µmol kg⁻¹ (partial bottom cells retained), and the maximum was located with `nanargmax` on that 2‑D field, longitude wrapped to (−180, 180].

**Basins.** Since the Ofx `basin` field is absent locally, basins come from `basinmask_01` only. Codes 10 and 11 were identified as Southern (mean lat < −55°) and Arctic (mean lat > 65°); the three largest remaining codes were named by their longitude weighting (150–285°E Pacific, 40–120°E Indian, 290–20°E Atlantic) — code 2 → Pacific, 1 → Atlantic, 3 → Indian — and all smaller codes were folded into the parent oceans by iterative nearest-neighbour dilation. Atlantic, Pacific and Indian only were summed; Southern and Arctic excluded.

**Figures.** Native-grid fields were interpolated to 1° by nearest-neighbour `griddata` for display only; maps span 0–360° and 80°S–90°N.

## 3. Key numerical results

**OMZ volume (10¹⁵ m³) and cell counts**

| Source | < 20 µmol kg⁻¹ | cells | < 60 µmol kg⁻¹ | cells |
|---|---|---|---|---|
| WOA23 1971–2000 | 14.2 | 34 250 | 80.5 | 180 370 |
| GFDL‑ESM4 historical | 42.9 | 1 086 136 | 104 | 1 131 345 |
| MPI‑ESM1‑2‑LR historical | 91.6 | 1 196 749 | 212 | 1 230 955 |
| GFDL‑ESM4 ssp245 | 41.2 | 1 084 781 | 109 | 1 133 527 |
| MPI‑ESM1‑2‑LR ssp245 | 88.9 | 1 194 689 | 211 | 1 230 906 |
| CESM2 ssp245 | 1510* | 3 357 732 | 1620* | 3 726 659 |

\*CESM2 values are contaminated — see §6.

**Maximum 60 µmol kg⁻¹ OMZ thickness**

| Source | Max thickness (m) | Lat | Lon |
|---|---|---|---|
| WOA23 | 3125.0 | 7.50 | 121.50 |
| GFDL‑ESM4 historical | 4343.2 | 13.50 | −93.50 |
| MPI‑ESM1‑2‑LR historical | 4960.0 | −9.39 | −81.22 |
| GFDL‑ESM4 ssp245 | 4318.2 | 13.50 | −93.50 |
| MPI‑ESM1‑2‑LR ssp245 | 4950.0 | −9.39 | −81.22 |
| CESM2 ssp245 | 5500.0* | 43.29 | 26.55 |

**Scenario change (ssp245 2071‑2100 minus historical 1971‑2000, % of historical volume; + = expansion)**

- GFDL‑ESM4: **−4.07 %** at 20 µmol kg⁻¹, **+4.59 %** at 60 µmol kg⁻¹.
- MPI‑ESM1‑2‑LR: **−3.03 %** at 20 µmol kg⁻¹, **−0.71 %** at 60 µmol kg⁻¹.

**Total dissolved O₂ inventory (mol)**

- GFDL‑ESM4 historical **2.20 ×10¹⁷**, ssp245 **2.15 ×10¹⁷**.
- MPI‑ESM1‑2‑LR historical **2.25 ×10¹⁷**, ssp245 **2.20 ×10¹⁷**.
- CESM2 ssp245 **2.02 ×10¹⁷**.

**Basin split at 60 µmol kg⁻¹ (10¹⁵ m³, marginal seas folded in, Southern/Arctic excluded)**

| Source | Pacific | Indian | Atlantic |
|---|---|---|---|
| WOA23 | 67.5 | 11.7 | 0.779 |
| GFDL‑ESM4 historical | 92.0 | 5.60 | 6.34 |
| MPI‑ESM1‑2‑LR historical | 181 | 19.7 | 11.2 |

**WOA23 within its own 60 µmol kg⁻¹ OMZ (volume-weighted):** AOU **271.3 µmol kg⁻¹**, oxygen saturation **11.3 %**.

**Which model is better / do they agree?** GFDL‑ESM4 reproduces the observed volume better (combined absolute misfit 52.3 ×10¹⁵ m³ versus 209 ×10¹⁵ m³ for MPI‑ESM1‑2‑LR); both overestimate, MPI by roughly a factor of 2.6 at 60 µmol kg⁻¹. The two models **do not agree on the sign at 60 µmol kg⁻¹** (one expands, one contracts) and the magnitudes differ by more than a factor of six; at 20 µmol kg⁻¹ they agree in sign (both shrink) with similar magnitude.

## 4. Figures produced

**Figure 1** — a single row of three global maps of historical OMZ thickness at 60 µmol kg⁻¹ (GFDL‑ESM4, MPI‑ESM1‑2‑LR, WOA23), one shared sequential (`viridis`) scale anchored at 0 and clipped at the pooled 98th percentile of the three fields (annotated in the title), each panel on its own native grid, 0–360° longitude and 80°S–90°N. The design intent is a directly comparable visual test of whether the models place thick suboxic columns where WOA23 observes them.

**Figure 2** — three panels: the change in 60 µmol kg⁻¹ OMZ thickness (ssp245 2071‑2100 minus historical 1971‑2000) for GFDL‑ESM4 and MPI‑ESM1‑2‑LR on a shared diverging `RdBu_r` scale centred on zero and clipped at the 98th percentile of |change| (both panels share one colour bar), plus CESM2's ssp245 OMZ thickness on its own sequential `viridis` scale (clipped at its own 98th percentile), annotated in-panel that no change can be formed because ssp245 is r10i1p1f1 and no r10i1p1f1 historical `o2` exists. The logic is a common-scale red/blue field for the two models whose historical and ssp245 members are a continuous trajectory, versus a standalone state map for the model that has no counterpart.

## 5. Interpretation

- The observed 60 µmol kg⁻¹ zone carries a very high AOU (271 µmol kg⁻¹) together with very low saturation (11 %): oxygen has been consumed in situ, so the zone is **respiration-driven**, not merely an old, poorly ventilated water mass. This is the correct reading of the pair, and it is the sense in which WOA's AOU/percent-saturation fields are informative (they rest on WOA's own solubility and pressure conventions, so they are not strictly comparable to the models' OMIP-derived oxygen).
- Both models place the thickest suboxic columns in the **eastern tropical Pacific** (13.5°N/93.5°W in GFDL‑ESM4; 9.4°S/81.2°W in MPI‑ESM1‑2‑LR), where WOA23's maximum sits in the **western Pacific/Indonesian seas** (7.5°N/121.5°E). WOA23's thinner maximum (3125 m vs 4343 and 4960 m) is consistent with its much smaller total volume.
- The Pacific dominates the OMZ in all three sources (67.5, 92.0, 181 ×10¹⁵ m³), with the Indian second and the Atlantic far smaller in WOA23 (0.78) than in either model (6.3, 11.2) — the models spread low oxygen into the Atlantic much more than observed.
- Inventories differ by only ~10 % between models (2.02–2.25 ×10¹⁷ mol) while OMZ volumes differ by factors of two to three, i.e. the models differ far more in how the oxygen is distributed than in how much there is.

## 6. Caveats that the numbers force

- **CESM2 ssp245 volumes are not physical.** Its `volcello` contains finite values in every one of the 7 372 800 cells (no land NaN mask), and its `o2` is 0 over land, so the `< 20` and `< 60` sums count the entire water column everywhere: 1510 and 1620 ×10¹⁵ m³, its maximum "thickness" is 5500 m (the full model depth) at 43.29°N/26.55°E, and its mean column thickness is 3352 m. These values exceed plausible global OMZ volumes and must be read as an artefact of the missing ocean mask in the volume sum, not as a modelled anoxic ocean. Its inventory (2.02 ×10¹⁷ mol) is unaffected, because land cells contribute zero oxygen to that integral.
- Volumes are integrals over **each model's own wet domain and its own `volcello`**, with no interpolation to a common grid; total ocean area and bathymetry differ by about a per cent between models, so part of every difference is a domain difference. Each source is evaluated on its own grid against the WOA23 reference.
- The historical (1971–2000) and ssp245 (2071–2100) windows are different 30‑year periods, so the reported changes are **period differences, not trends**. They are only meaningful for GFDL‑ESM4 and MPI‑ESM1‑2‑LR, whose r1i1p1f1 members form one continuous trajectory.
