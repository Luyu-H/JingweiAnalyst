---
id: cmip6-grid-labels
title: "CMIP6 grid labels gn / gr and what they imply for the three local models"
type: dataset
summary: "gn means the model's native grid and gr regridded to the provider's preferred target grid; CESM2 and MPI-ESM1-2-LR are gn (curvilinear, native levels), GFDL-ESM4 is gr (conservatively regridded from 0.5° tripolar hybrid layers to 1°×1°, 35 z levels), which affects masks, conservation and which operations are exact."
keywords: [grid_label, 网格标签, gn, gr, gr1, native grid, 原生网格, regridded, 重网格, nominal_resolution, tripolar, displaced pole, bipolar, curvilinear, 曲线网格]
datasets: [cmip6_historical, cmip6_ssp245]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [areacello, volcello]
tools: [regrid_horizontal, load_areas, load_volumes_for_cmip6]
severity: should
related: [omip-native-grid-vs-regridded, model-gfdl-esm4-configuration, velocity-components-are-grid-relative]
sources:
  - citation: "WCRP-CMIP CMIP6_CVs, CMIP6_grid_label.json, CV collection version 6.2.60.0."
    locator: "grid_label gn, gr, gr1, gm, gnz"
    url: https://github.com/WCRP-CMIP/CMIP6_CVs
  - citation: "Griffies, S. M., et al. (2016). OMIP contribution to CMIP6: experimental and diagnostic protocol for the physical component of the Ocean Model Intercomparison Project. Geosci. Model Dev., 9, 3231–3296."
    locator: "Appendix A3.3–A3.4 (p. 3240); Appendix B1.2 (p. 3243)"
    url: https://doi.org/10.5194/gmd-9-3231-2016
---

## Key points

- Grid labels [CVs grid_label]:
  - `gn`: "data reported on a model's native grid"
  - `gr`: "regridded data reported on the data provider's preferred target grid"
  - `gr1`, `gr2`, …: other regridded grids. `gnz`: zonal mean on the native
    latitude grid. `gm`: global mean.
- OMIP: structured-grid models should archive on the native grid. The same
  variable on native and spherical grids keeps the same CMOR name, and the
  two "can never coexist in the same file". GFDL described its plan to coarsen
  ocean output to a 1° spherical grid to reduce archive size. [OMIP-phys A3.3–A3.4]
- Models with a moving vertical coordinate (z\*, isopycnal, hybrid) have
  **time-dependent cell thickness and volume**. [OMIP-phys B1.2]

## Checked in the local files

`grid_label` and `grid` global attributes of `thetao`, read on 2026-09-24:

| model | grid_label | grid attribute | nominal_resolution |
|---|---|---|---|
| CESM2 | gn | "native gx1v7 displaced pole grid (384x320 latxlon)" | 100 km |
| GFDL-ESM4 | **gr** | "ocean data regridded from tripolar - nominal 0.5 deg latitude/longitude to 180,360; interpolation method: conserve_order1" | 1x1 degree |
| MPI-ESM1-2-LR | gn | "gn" | 250 km |

## Implications for this project

- GFDL-ESM4 values are already remapped, so "native-grid" numbers for GFDL
  are numbers on its 1° archive grid. The land–sea mask is the regridded one,
  and budgets are only approximately conservative. The native hybrid layers
  (75) were remapped to 35 z levels, which is why `volcello` is monthly for
  GFDL-ESM4.
- Resolution differs by a factor of about 2.5 between CESM2 (100 km) and MPI
  (250 km). Features narrower than the coarse grid (boundary currents,
  equatorial undercurrent, coastal upwelling) are not comparable in amplitude.
- Within one model, variables can sit at different staggered points on the
  same gn grid (see `velocity-components-are-grid-relative`).

## Source excerpts

> "gr": "regridded data reported on the data provider's preferred target grid"
> (CMIP6_grid_label.json)
