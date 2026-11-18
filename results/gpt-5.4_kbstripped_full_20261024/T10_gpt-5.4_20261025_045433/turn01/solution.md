Using the executed workflow, I assessed OMZs from WOA23 and CMIP6 model period means, with model oxygen first converted from native CMIP6 concentration units to \( \mu\mathrm{mol}\,\mathrm{kg}^{-1} \) using density computed from the period-mean temperature and salinity fields. Historical model means were formed from monthly data over 1971-01 to 2000-12 for GFDL-ESM4 and MPI-ESM1-2-LR; ssp245 means were formed over 2071-01 to 2100-12 for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR. WOA23 used the 1971-2000 annual climatology directly. OMZs were defined by oxygen strictly below 20 and strictly below 60 \( \mu\mathrm{mol}\,\mathrm{kg}^{-1} \). Global OMZ volumes were computed on each source’s own grid from cell volumes; OMZ thickness was the column-summed thickness of waters below threshold. For WOA23 basin volumes at the 60 threshold, the WOA basin mask was applied directly. For the two historical models, 60-threshold basin volumes were derived after horizontally regridding OMZ volume fields to the WOA23 grid and then applying the WOA basin mask. WOA23 mean apparent oxygen utilization and percent oxygen saturation were computed within the WOA23 OMZ at the 60 threshold using volume-weighted means.

## 1) Main scientific answer

WOA23 shows a substantial global OMZ, but both CMIP6 historical models overestimate its volume, especially MPI-ESM1-2-LR. Of the two historical models, **GFDL-ESM4 reproduces the observed global OMZ volume much better** than MPI-ESM1-2-LR. Under ssp245, the two models **agree that the most severe OMZ core** (\(<20\)) **shrinks slightly**, but they **do not agree at the broader \(<60\) threshold**: GFDL-ESM4 projects a modest expansion, while MPI-ESM1-2-LR projects a very slight contraction. CESM2, available only for ssp245, produces a large late-century OMZ comparable to the other models’ biased-high historical/ssp245 states.

## 2) Key numerical results

### Global OMZ volume and maximum 60-threshold thickness

All OMZ volumes below are in \(10^{15}\,\mathrm{m}^3\).

| Source | OMZ volume \(<20\) | OMZ volume \(<60\) | Max thickness at \(<60\) (m) | Location |
|---|---:|---:|---:|---|
| WOA23 | 14.2 | 80.5 | 3125.0 | 7.50°N, 121.50°E |
| GFDL-ESM4 historical | 43.1 | 105 | 4343.2 | 13.50°N, 93.50°W |
| MPI-ESM1-2-LR historical | 92.0 | 212 | 4960.0 | 9.39°S, 81.22°W |
| CESM2 ssp245 | 60.9 | 165 | 3631.9 | 7.45°N, 121.44°E |
| GFDL-ESM4 ssp245 | 41.3 | 110 | 4318.2 | 13.50°N, 93.50°W |
| MPI-ESM1-2-LR ssp245 | 89.2 | 211 | 4950.0 | 9.39°S, 81.22°W |

### Historical-to-ssp245 change in OMZ volume

Positive means expansion.

| Model | Change in \(<20\) OMZ | Change in \(<60\) OMZ |
|---|---:|---:|
| GFDL-ESM4 | -4.06% | +4.41% |
| MPI-ESM1-2-LR | -3.03% | -0.63% |

### Total dissolved oxygen inventory

Inventories are in moles.

| Source | Inventory |
|---|---:|
| GFDL-ESM4 historical | \(2.20 \times 10^{17}\) |
| MPI-ESM1-2-LR historical | \(2.23 \times 10^{17}\) |
| CESM2 ssp245 | \(2.02 \times 10^{17}\) |
| GFDL-ESM4 ssp245 | \(2.15 \times 10^{17}\) |
| MPI-ESM1-2-LR ssp245 | \(2.17 \times 10^{17}\) |

The model oxygen inventory declines from historical to ssp245 in both GFDL-ESM4 and MPI-ESM1-2-LR, consistent with a less oxygenated future ocean in bulk inventory terms, even though the \(<60\) OMZ response differs by model.

### Basin OMZ volume at the 60 threshold

Pacific, Indian and Atlantic only; Southern and Arctic excluded. Volumes are in \(10^{15}\,\mathrm{m}^3\).

| Source | Pacific | Indian | Atlantic |
|---|---:|---:|---:|
| WOA23 | 67.5 | 11.7 | 1.29 |
| GFDL-ESM4 historical | 92.8 | 5.43 | 6.39 |
| MPI-ESM1-2-LR historical | 390 | 40.7 | 16.4 |

Interpretation:
- Observationally, the OMZ is overwhelmingly **Pacific-dominated**, with a secondary Indian Ocean contribution and only a small Atlantic contribution.
- GFDL-ESM4 captures that broad ranking, but still overstates the Pacific and Atlantic OMZs and understates the Indian share.
- MPI-ESM1-2-LR is far too large in all three basins, especially the Pacific.

### WOA23 OMZ biogeochemical context at the 60 threshold

Within the WOA23 OMZ (\(<60\)):
- Mean apparent oxygen utilization: **271.3 \( \mu\mathrm{mol}\,\mathrm{kg}^{-1} \)**
- Mean percent oxygen saturation: **11.3%**

These values indicate **very high oxygen consumption relative to saturation and extremely low saturation itself**, implying an OMZ maintained primarily by **respiration acting on poorly ventilated waters**, rather than by weak solubility effects alone.

## 3) Which model performs better, and do the models agree on future change?

### Model fidelity against WOA23
Using absolute error in global OMZ volume at both thresholds:
- **GFDL-ESM4** combined absolute error: **53.6 \(\times 10^{15}\,\mathrm{m}^3\)**
- **MPI-ESM1-2-LR** combined absolute error: **210 \(\times 10^{15}\,\mathrm{m}^3\)**

So **GFDL-ESM4 is clearly closer to WOA23**.

### Agreement on scenario response
- At \(<20\), both models show a **small contraction** of the most oxygen-depleted waters.
- At \(<60\), they **disagree on sign**:
  - GFDL-ESM4: **expansion** (+4.41%)
  - MPI-ESM1-2-LR: **slight contraction** (-0.63%)

Thus, the two models **do not agree on the sign or size of the broader OMZ change** under ssp245.

## 4) Interpretation of the pattern

A clear pattern emerges:
1. **Observed OMZ structure** is large but not as extensive as in either historical model.
2. **GFDL-ESM4** is biased high, but still within a more plausible range relative to WOA23.
3. **MPI-ESM1-2-LR** strongly overestimates OMZ extent and thickness, especially in the Pacific.
4. In future projections, both models lose global oxygen inventory, but that does **not translate identically into threshold-based OMZ expansion**:
   - The severe OMZ core (\(<20\)) shrinks slightly in both.
   - The broader hypoxic envelope (\(<60\)) is model-dependent.

This implies that projected OMZ change is sensitive not just to overall deoxygenation, but to how each model redistributes oxygen around threshold values and across depth and region.

## 5) Visualization summary

Two figures were generated.

### Figure 1
A one-row, three-panel global map comparison of **historical OMZ thickness at the 60 \( \mu\mathrm{mol}\,\mathrm{kg}^{-1} \) threshold** for:
- GFDL-ESM4 historical,
- MPI-ESM1-2-LR historical,
- WOA23.

The model thickness fields were horizontally regridded to the WOA23 grid for display. All three panels use one shared sequential colour scale in metres, clipped at the pooled 98th percentile of the displayed thickness values.

### Figure 2
A one-row, three-panel global map figure showing:
- GFDL-ESM4 OMZ thickness change at the 60 threshold: ssp245 minus historical,
- MPI-ESM1-2-LR OMZ thickness change at the 60 threshold: ssp245 minus historical,
- CESM2 ssp245 OMZ thickness itself, since no historical oxygen field was available for CESM2.

The first two change panels were computed on native grids, then regridded to the WOA23 grid for display, and plotted on a shared diverging colour scale centred on zero. The CESM2 panel uses its own sequential scale and is explicitly labelled as a thickness map rather than a change map.

## 6) Important caveat from the executed workflow

For MPI-ESM1-2-LR, the global maxima of 60-threshold thickness were successfully located, but the coordinates came through as native-grid latitude/longitude metadata rather than standard `lat`/`lon` fields in the final reporting table. The actual reported maximum locations are:
- historical: **9.39°S, 81.22°W**
- ssp245: **9.39°S, 81.22°W**

## Bottom line

Relative to WOA23, both CMIP6 models overestimate OMZ volume, with **GFDL-ESM4 performing substantially better** than **MPI-ESM1-2-LR**. The observed OMZ is concentrated mainly in the Pacific and is characterized by **high AOU and very low oxygen saturation**, consistent with **respiration in poorly ventilated waters**. Under ssp245, both models project a **small reduction in the most extreme low-oxygen volume** (\(<20\)), but they **diverge for the broader \(<60\) OMZ**, so the sign of future large-scale OMZ change is not robust across these two models.
