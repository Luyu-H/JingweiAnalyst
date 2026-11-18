---
id: omip-native-grid-vs-regridded
title: "When to work on the native grid and when regridding is acceptable"
type: method
summary: "CMIP6 protocol puts diagnostics on native grids for conservation; regridding changes the land-sea mask and breaks conservation even with \"conservative\" methods, so compute integrals and budgets natively and regrid only for maps and pattern comparisons."
keywords: [native grid, 原生网格, regrid, 插值, 重网格, conservative regridding, 守恒插值, bilinear, land-sea mask, 海陆掩膜, pattern correlation, difference map, 差值图, WOA 1 degree grid, zig-zag, vertical interpolation, 垂直插值, fixed depths, 固定深度, nutricline, nitracline, 营养盐跃层]
datasets: [cmip6_historical, cmip6_ssp245, woa23, glodap]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [thetao, so, o2, no3, areacello]
tools: [regrid_horizontal, regrid_vertical, build_comparison_mask, compute_difference_metrics, compute_nutricline_depth]
severity: should
related: [cmip6-grid-labels, velocity-components-are-grid-relative, omip-weighted-means-and-integrals]
sources:
  - citation: "Griffies, S. M., et al. (2016). OMIP contribution to CMIP6: experimental and diagnostic protocol for the physical component of the Ocean Model Intercomparison Project. Geosci. Model Dev., 9, 3231–3296."
    locator: "Appendix A3–A3.4 (pp. 3239–3240); Appendix A4 (p. 3240); Appendix C2–C2.1 (pp. 3244–3245); Appendix C4 (pp. 3245–3246)"
    url: https://doi.org/10.5194/gmd-9-3231-2016
---

## Key points

- Two aims conflict: spherical-grid output is easy to map and difference,
  while native-grid output preserves conservation for budgets. "regridding …
  generally involves a land–sea mask distinct from the model native mask. When
  changing the land–sea mask it is difficult to ensure conservation even when
  using 'conservative' regridding methods." [OMIP-phys A3]
- The loss of conservation "may be of minor consequence for many metric-based
  analysis purposes. However, it sacrifices the needs of budget analyses."
  Regridding "merely for visualization purposes and/or to produce difference
  maps" can ignore non-conservation. [OMIP-phys A3, A3.2]
- Structured-grid models archive on their native grid. The standard spherical
  targets are the 1° or ¼° grids of the World Ocean Atlas. [OMIP-phys A3.3]
- In the vertical, models on z, z\*, p or p\* need no regridding. Coarsening to
  the 33 WOA standard levels is allowed, but "care should be given to models
  that do not share the same 5500 m depth as the World Ocean Atlas".
  [OMIP-phys A4, C2]
- For transports across sections (overturning, passages, heat transport) on
  non-spherical grids, compute along a native-grid "zig-zag" path that
  approximates the section rather than regridding fluxes. [OMIP-phys C4]

## Implications for this project

- A good pattern for model-versus-observation work: compute every **number**
  (global and basin means, inventories, integrals) on each source's native grid
  with its own weights. Regrid to the WOA 1° grid only for **maps, pattern
  correlation and RMSE on a common mask**. Say which numbers come from which
  grid. The exception is an analysis whose result depends on the grid and is
  compared between sources, such as an EOF mode with its explained variance
  and its regression pattern (see `eof-modes-and-pdo-index`): put every
  source, the observations included, on one common grid before the analysis,
  and take its numbers there.
- After regridding a field for a cell-by-cell comparison, build one shared
  valid-data mask for all sources before computing pattern statistics. Coastal
  cells otherwise differ between the sources (see `common-comparison-mask`,
  which also says when a pair's own shared cells are the sample instead, as
  for the patterns of a variability mode).
- Regrid only for the cell-by-cell comparison. An average, a trend or an
  integral belongs on the field's own grid with its own cell areas, because
  regridding replaces those with the target's.
- Method in `regrid_horizontal`: 'conservative' is the right choice going to a
  coarser grid and before any area-weighted step; 'bilinear' is the fallback
  when a model's corner mesh is unusable; 'nearest_s2d' is for masks and other
  categorical fields. A regridded RMSE or pattern correlation depends on the
  method, so state the method whenever one is quoted.
- Vertical interpolation onto another level set changes layer-integrated
  quantities. Integrate thickness-weighted quantities (inventories, OMZ
  thickness) on the original levels.
- Values at fixed depths (a bias "at 100, 500 and 1000 m") need each product
  interpolated to those depths (`regrid_vertical`), not the nearest level
  selected: vertical grids differ between products and rarely line up, so each
  product would answer at a different depth, and the products differ most
  where the gradients are steepest. This interpolation is for comparison, not
  for budgets.
- Diagnostics located at a vertical extremum, such as a nutricline (the depth of
  the maximum vertical gradient), are quantised by the vertical grid. Where the
  profile is smooth, in the subtropical gyres especially, the maximum can jump
  between levels, so read such a map together with the gradient itself.

## Common mistakes

- Computing a global inventory after regridding to 1° and comparing it with
  the native-grid value as if they should agree.
- Regridding velocity components as scalars (see
  `velocity-components-are-grid-relative`).

## Source excerpts

> "The absence of conservation may be of minor consequence for many
> metric-based analysis purposes. However, it sacrifices the needs of budget
> analyses." (OMIP-phys, A3)
