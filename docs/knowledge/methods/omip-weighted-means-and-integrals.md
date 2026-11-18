---
id: omip-weighted-means-and-integrals
title: "How OMIP defines global and regional means and integrals of ocean fields"
type: method
summary: "OMIP global means are volume-weighted (Boussinesq) or mass-weighted (non-Boussinesq) over all ocean cells; surface means are area-weighted over the top level; integrals use every grid point in the region with no subsampling."
keywords: [global mean, 全球平均, volume weighted, 体积加权, area weighted, 面积加权, mass weighted, 质量加权, spatial integral, 空间积分, regional mean, 区域平均, thetaoga, tosga, soga, sosga, Boussinesq, weighted RMSE, pattern correlation, 空间相关, area fraction, 面积占比, mean of ratios, ratio of means, 比值平均]
datasets: [cmip6_historical, cmip6_ssp245]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [thetao, so, tos, areacello, volcello]
tools: [compute_basic_stat, load_areas, load_volumes_for_cmip6, select_basin, compute_difference_metrics, analyze_correlation, compute_exceedance_fraction, compute_ratio]
severity: should
related: [cmip6-static-grid-fields, omip-native-grid-vs-regridded, cf-cell-methods]
sources:
  - citation: "Griffies, S. M., et al. (2016). OMIP contribution to CMIP6: experimental and diagnostic protocol for the physical component of the Ocean Model Intercomparison Project. Geosci. Model Dev., 9, 3231–3296."
    locator: "Appendix A2 (p. 3239); Appendix C1 (p. 3244); Appendix G3 (p. 3251); Appendix F3, Eq. F13 (pp. 3250–3251); Appendix H12.1, Eqs. H28–H29 (p. 3259); H13, Eq. H34 (p. 3260); H17–H18, Eqs. H35–H37 (p. 3261)"
    url: https://doi.org/10.5194/gmd-9-3231-2016
  - citation: "Mauritsen, T., et al. (2019). Developments in the MPI-M Earth System Model version 1.2 (MPI-ESM1.2). J. Adv. Model. Earth Syst., 11, 998–1038."
    locator: "Section 4 (p. 1010): MPIOM hydrostatic and Boussinesq"
    url: https://doi.org/10.1029/2018MS001400
---

## Key points

- Spatial integrals "are computed using all grid points within the region or
  line. There is no sub-sampling." [OMIP-phys A2, C1]
- Global mean of a 3-D tracer T [OMIP-phys H12.1]:
  - Boussinesq model, **volume-weighted**: ⟨T⟩ = Σ T dA dz / Σ dA dz (Eq. H29)
  - non-Boussinesq model, **mass-weighted**: ⟨T⟩ = Σ ρ T dA dz / Σ ρ dA dz
    (Eq. H28)
  The same forms hold for salinity (Eqs. H35–H36).
- Global mean SST and SSS are **area-weighted over the top level**:
  tosga = Σ θ(k=1) dA / Σ dA (Eq. H34, likewise sosga, Eq. H37).
  [OMIP-phys H13, H18]
- Because the terms of these averages vary in time, OMIP computes the official
  global means online each time step (products of time-dependent terms).
  [OMIP-phys H12.1, H17]
- MPIOM (MPI-ESM1.2) is hydrostatic and Boussinesq. [MPI §4]

## Implications for this project

- Offline, with monthly fields: a global mean of `thetao` is
  Σ thetao · volcello / Σ volcello over valid ocean cells. A surface mean is
  Σ tos · areacello / Σ areacello. Unweighted `mean()` over (y, x) on a
  curvilinear grid is wrong. On a regular grid, unweighted means over latitude
  are also wrong.
- The weight must be restricted to the **same cells** as the valid data
  (NaN-aware weights). Otherwise the denominator counts land or missing cells.
- Volume weighting and mass weighting differ only through density variations.
  In-situ density varies by no more than about 2 % around 1035 kg m-3 over
  most of the ocean (OMIP-phys G3). For a Boussinesq model the volume-weighted
  form is the protocol form. The papers read here confirm the Boussinesq formulation only
  for MPIOM. For the other two models, state the assumption.
- A mean over a region follows the same rule, with the region mask applied on
  the model's own grid.
- Use a weighted mean for any physical average. A plain mean is right only
  when every element represents the same amount of ocean, which is rarely
  true: cells shrink towards the poles, and layer thicknesses vary by more than
  an order of magnitude with depth. With the tools: `compute_basic_stat` with
  `op='area_weighted_mean'` for a horizontal mean, and
  `op='volume_weighted_mean'` for any mean whose dimensions include depth (the
  mean temperature or concentration of a depth range). Its plain `op='mean'`
  weights every element equally.
- The same holds for every statistic taken over a map, not only for means: an
  unweighted RMSE, bias or pattern correlation counts a polar cell as heavily
  as a tropical one many times its size. `compute_difference_metrics` weights
  by area by default; `analyze_correlation` needs `weighting='area'` for a
  pattern correlation between two maps.
- A share of the ocean beyond a threshold is a share of its area (or volume),
  not of its cells, and missing cells must stay out of the denominator.
  `compute_exceedance_fraction` does both and returns the amount behind the
  percentage (m2 or m3), which belongs next to the percentage when the question
  is how much ocean.
- The mean of a ratio is not the ratio of the means. Averaging a cell-by-cell
  ratio such as N:P weights every cell's composition by its volume alone; the
  ratio of the two volume-weighted means weights each cell by its content and
  is pulled towards the richest water. Both are defensible, so state which one
  is reported. `compute_ratio` gives the cell-by-cell ratio; the ratio of means
  is the weighted mean of each field, divided.

## Common mistakes

- Averaging a depth profile level by level with equal weight. Level thickness
  varies from 10 m to 250 m or more.
- Dividing an area-weighted sum by the total area of the region including
  land.

## Source excerpts

> "In a Boussinesq model, the mean is computed as the volume weighted mean"
> (OMIP-phys, H12.1)
