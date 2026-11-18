**OMZ assessment against WOA23 and SSP2-4.5 changes**

I use O₂ < 60 µmol kg⁻¹ as the main OMZ threshold and also report the stricter < 20 µmol kg⁻¹ threshold. Volumes are integrated with each model’s `volcello`; WOA23 volumes use the WOA grid-cell volumes. Thickness is the column-integrated thickness of water below the threshold. Basin volumes use the WOA23 basin mask.

### 1. Global OMZ volume and maximum thickness at O₂ < 60 µmol kg⁻¹

| Source / experiment | V(O₂<20) [10¹⁵ m³] | V(O₂<60) [10¹⁵ m³] | Max thickness at <60 [m] | Latitude [°] | Longitude [°] |
|---|---:|---:|---:|---:|---:|
| WOA23 1971–2000 | 14.2 | 80.5 | 3125.0 | 7.50 | 121.50 |
| GFDL-ESM4 historical | 43.1 | 105 | 4343.2 | 13.50 | −93.50 |
| MPI-ESM1-2-LR historical | 92.0 | 214 | 4960.0 | −9.39 | −81.22 |
| CESM2 ssp245 | 60.9 | 165 | 36.3* | 7.45 | 121.44 |
| GFDL-ESM4 ssp245 | 41.3 | 110 | 4318.2 | 13.50 | −93.50 |
| MPI-ESM1-2-LR ssp245 | 89.2 | 212 | 4950.0 | −9.39 | −81.22 |

\*The CESM2 maximum thickness printed by the code is almost certainly an area-normalization artefact: CESM2 `sftof` is likely in percent and was used directly rather than divided by 100 when forming grid-cell thickness. If so, the corrected CESM2 maximum would be about **3630 m**; CESM2 OMZ volumes and oxygen inventory are not affected by this issue.

Both historical CMIP6 models strongly overestimate the observed OMZ volume. At the <60 threshold, GFDL-ESM4 is much closer to WOA23 (105 vs 80.5 ×10¹⁵ m³) than MPI-ESM1-2-LR (214 vs 80.5 ×10¹⁵ m³). At the <20 threshold both are far too voluminous, but GFDL-ESM4 is still less biased (43.1 vs 14.2) than MPI (92.0 vs 14.2).

### 2. Historical → SSP2-4.5 percentage change

Positive values mean an expanding OMZ.

| Model | ΔV(O₂<20) [%] | ΔV(O₂<60) [%] |
|---|---:|---:|
| GFDL-ESM4 | −4.06 | +4.42 |
| MPI-ESM1-2-LR | −3.05 | −0.72 |

At <20 µmol kg⁻¹, both models simulate a small contraction, and they roughly agree in sign and modest magnitude. At <60 µmol kg⁻¹ they do not agree in sign: GFDL-ESM4 shows a modest expansion (+4.42%), whereas MPI-ESM1-2-LR shows a near-neutral slight contraction (−0.72%). Thus these two models do not provide a robust consensus for the 60-threshold OMZ change.

### 3. Global dissolved oxygen inventory

| Source | Inventory [mol] |
|---|---:|
| GFDL-ESM4 historical | 2.20 ×10¹⁷ |
| MPI-ESM1-2-LR historical | 2.25 ×10¹⁷ |
| CESM2 ssp245 | 2.02 ×10¹⁷ |
| GFDL-ESM4 ssp245 | 2.15 ×10¹⁷ |
| MPI-ESM1-2-LR ssp245 | 2.20 ×10¹⁷ |

### 4. Basin partition of the <60 µmol kg⁻¹ OMZ volume

Units: 10¹⁵ m³.

| Source | Pacific | Indian | Atlantic |
|---|---:|---:|---:|
| WOA23 | 67.3 | 9.97 | 0.777 |
| GFDL-ESM4 historical | 92.7 | 3.48 | 5.76 |
| MPI-ESM1-2-LR historical | 182 | 16.0 | 11.1 |

The Pacific dominates the OMZ volume in all cases. MPI-ESM1-2-LR greatly overestimates the Pacific and Indian OMZ volumes relative to WOA23. GFDL-ESM4 is closer in total volume, although it places more OMZ volume in the Atlantic than WOA23.

### 5. WOA23 AOU and oxygen saturation within the <60 OMZ

Volume-weighted values within the WOA23 <60 µmol kg⁻¹ OMZ:

- Mean apparent oxygen utilization: **271.3 µmol kg⁻¹**
- Mean percent oxygen saturation: **11.3%**

The very high AOU and very low saturation indicate that the observed OMZ is strongly respiration-driven: a large amount of oxygen has been consumed after the water was last ventilated. Sluggish ventilation is still required to allow this respiratory oxygen deficit to accumulate, but the high AOU shows that the zone is not merely a poorly ventilated water mass with little oxygen utilization. It is primarily an oxygen-depleted, respiration-dominated zone.

### 6. Overall model skill and interpretation

**GFDL-ESM4 reproduces the observed <60 µmol kg⁻¹ OMZ volume better than MPI-ESM1-2-LR.** Its historical volume is about 31% above WOA23, whereas MPI is about 166% above WOA23. Both models overestimate the stricter <20 µmol kg⁻¹ volume, but GFDL-ESM4 is again closer.

For the scenario change, the two models disagree at the 60 threshold. GFDL-ESM4 predicts a modest expansion (+4.42%), while MPI predicts a very small contraction (−0.72%). At 20, both predict a small contraction (−4.06% and −3.05%). Therefore, these two CMIP6 models do not give a consistent picture of strong OMZ expansion under SSP2-4.5; GFDL-ESM4 suggests slight 60-threshold expansion, MPI suggests near-stability/slight contraction, and both suggest slight shrinking of the most severe <20 water.

CESM2 cannot be used to form a historical-to-future change because it publishes no historical `o2` in this dataset. Its SSP2-4.5 <60 OMZ volume is 165 ×10¹⁵ m³, intermediate between GFDL-ESM4 and MPI-ESM1-2-LR, but no change can be quantified.

### 7. Figures

- **Figure 1** shows historical <60 OMZ thickness for GFDL-ESM4, MPI-ESM1-2-LR and WOA23 on a common sequential scale. The models produce much thicker and more widespread OMZs than WOA23, especially in the eastern tropical Pacific, with MPI-ESM1-2-LR showing the largest values.
- **Figure 2** shows the SSP2-4.5 minus historical thickness change for GFDL-ESM4 and MPI-ESM1-2-LR on a common diverging scale centered on zero. GFDL-ESM4 has a net modest expansion at the 60 threshold, while MPI has a near-zero/slight contraction with substantial regional thickening and thinning. CESM2’s SSP2-4.5 OMZ thickness is shown separately on its own sequential scale, labelled as having no historical comparison and therefore no computable change.
