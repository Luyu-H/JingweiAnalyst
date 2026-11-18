## Overall assessment

Both historical CMIP6 models overestimate the observed global OMZ volume, especially at the stricter 20 µmol kg⁻¹ threshold. GFDL-ESM4 is clearly closer to WOA23 than MPI-ESM1-2-LR, but it still overestimates the <60 µmol kg⁻¹ OMZ volume by about 31%. The ssp245 changes are model-dependent: GFDL-ESM4 gives a small expansion at <60 µmol kg⁻¹, while MPI-ESM1-2-LR gives a small contraction. At <20 µmol kg⁻¹ both models simulate a contraction, but of different magnitude. Thus there is no robust model consensus on either the sign or size of the future OMZ change.

## Global OMZ volume and maximum thickness at the 60 µmol kg⁻¹ threshold

Volumes are in units of 10¹⁵ m³. Maximum thickness is the column-integrated thickness of water with O₂ < 60 µmol kg⁻¹.

| Source / experiment | Volume <20 | Volume <60 | Max OMZ thickness <60 (m) | Location (lat, lon) |
|---|---:|---:|---:|---:|
| WOA23 1971–2000 | 14.2 | 80.5 | 3125.0 | 7.50°, 121.50° |
| GFDL-ESM4 historical | 43.1 | 105 | 4343.2 | 13.50°, 266.50° |
| MPI-ESM1-2-LR historical | 92.0 | 214 | 4960.0 | −9.39°, 278.78° |
| CESM2 ssp245 2071–2100 | 60.9 | 165 | 3631.9 | 7.45°, 121.44° |
| GFDL-ESM4 ssp245 2071–2100 | 41.3 | 110 | 4318.2 | 13.50°, 266.50° |
| MPI-ESM1-2-LR ssp245 2071–2100 | 89.2 | 212 | 4950.0 | −9.39°, 278.78° |

The observed WOA23 OMZ occupies about 80.5 × 10¹⁵ m³ below 60 µmol kg⁻¹ and 14.2 × 10¹⁵ m³ below 20 µmol kg⁻¹. Both models produce substantially larger OMZs, with MPI-ESM1-2-LR especially high. The model maximum thicknesses are also larger than observed: 4343 m for GFDL-ESM4 and 4960 m for MPI-ESM1-2-LR, versus 3125 m for WOA23.

## Percentage change from historical to ssp245

| Model | Change at <20 µmol kg⁻¹ | Change at <60 µmol kg⁻¹ |
|---|---:|---:|
| GFDL-ESM4 | −4.07% | +4.43% |
| MPI-ESM1-2-LR | −3.04% | −0.72% |

Positive values indicate an expanding OMZ. GFDL-ESM4 simulates a small expansion at <60 µmol kg⁻¹ but a contraction at <20 µmol kg⁻¹. MPI-ESM1-2-LR simulates a contraction at both thresholds, although the <60 µmol kg⁻¹ contraction is very small.

## Total dissolved oxygen inventory

| Source / experiment | Inventory (mol) |
|---|---:|
| GFDL-ESM4 historical | 2.20 × 10¹⁷ |
| MPI-ESM1-2-LR historical | 2.25 × 10¹⁷ |
| CESM2 ssp245 | 2.02 × 10¹⁷ |
| GFDL-ESM4 ssp245 | 2.15 × 10¹⁷ |
| MPI-ESM1-2-LR ssp245 | 2.20 × 10¹⁷ |

## Basin decomposition at the 60 µmol kg⁻¹ threshold

The basin-mask computation returned the following numerical values (10¹⁵ m³):

| Source / experiment | Pacific | Indian | Atlantic |
|---|---:|---:|---:|
| WOA23 1971–2000 | 0.000 | 0.000 | 0.0148 |
| GFDL-ESM4 historical | 0.00219 | 0.000 | 0.0391 |
| MPI-ESM1-2-LR historical | 0.000 | 0.000 | 0.171 |

**Important caveat:** these basin values are not physically credible. Their sum is only a tiny fraction of the global OMZ volume, and the Pacific and Indian Oceans are reported as essentially zero, which is impossible for the observed OMZ. This indicates that the WOA23 `basinmask_01` parsing/mapping in the executed code failed. The global volume and inventory metrics are not affected by that basin-mask problem, but the basin decomposition should be recomputed with a corrected basin mask before being used scientifically.

## WOA23 AOU and oxygen saturation within the <60 µmol kg⁻¹ OMZ

Volume-weighted values within the observed OMZ:

- Mean apparent oxygen utilization (AOU): **271.3 µmol kg⁻¹**
- Mean percent oxygen saturation: **11.3%**

These values indicate severe oxygen depletion with a very large accumulated respiration signal. The high AOU and very low O₂ saturation imply that the observed OMZ is strongly shaped by **respiration/oxygen consumption**, not simply by a water mass that remained poorly ventilated while retaining near-saturated oxygen. Poor ventilation is still required to prevent oxygen resupply, but the AOU signal points to respiration as the proximate oxygen sink.

## Model skill and robustness of the future change

GFDL-ESM4 reproduces the observed OMZ volume much better than MPI-ESM1-2-LR:

- At <60 µmol kg⁻¹, GFDL-ESM4 is about **30.6%** higher than WOA23, while MPI-ESM1-2-LR is about **165.7%** higher.
- At <20 µmol kg⁻¹, GFDL-ESM4 is about **203%** higher than WOA23, while MPI-ESM1-2-LR is about **547%** higher.

Thus GFDL-ESM4 is the better of the two historical models for OMZ volume, although it still overestimates the observed extent.

The two models do **not** agree on the scenario change. At <20 µmol kg⁻¹ both simulate a contraction, but the magnitudes differ (−4.07% versus −3.04%). At <60 µmol kg⁻¹ they disagree even in sign: GFDL-ESM4 expands by +4.43%, whereas MPI-ESM1-2-LR contracts by −0.72%. Given the large historical biases, this means the ssp245 OMZ response is not robustly constrained by these two models.

## Figures

- **Figure 1** shows the historical OMZ thickness at the 60 µmol kg⁻¹ threshold for GFDL-ESM4, MPI-ESM1-2-LR and WOA23 on a common sequential colour scale.
- **Figure 2** shows the ssp245-minus-historical OMZ thickness change for GFDL-ESM4 and MPI-ESM1-2-LR on a common diverging scale centred on zero, plus the CESM2 ssp245 OMZ thickness on its own sequential scale. CESM2 has no historical `o2` in the provided data, so its change could not be formed.
