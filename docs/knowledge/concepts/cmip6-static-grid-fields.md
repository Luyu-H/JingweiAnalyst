---
id: cmip6-static-grid-fields
title: "Meaning of the CMIP6 ocean grid fields areacello, volcello, deptho, sftof (and the absent thkcello/masscello)"
type: concept
summary: "areacello is tracer-cell horizontal area, volcello tracer-cell volume (extensive, \"area: sum\"), deptho the model sea-floor depth below the geoid for tracer cells, sftof the ocean percentage; layer thickness can be recovered as volcello/areacello."
keywords: [areacello, volcello, deptho, sftof, thkcello, masscello, rhozero, cell area, 网格面积, cell volume, 网格体积, bathymetry, 水深, sea area fraction, 海洋面积占比, layer thickness, 层厚, weights, 权重]
datasets: [cmip6_historical, cmip6_ssp245]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [areacello, volcello, deptho, sftof, thkcello, masscello]
tools: [load_areas, load_volumes_for_cmip6, load_depth_levels, identify_oxygen_minimum_zone]
severity: should
related: [omip-weighted-means-and-integrals, cf-cell-methods, cmip6-grid-labels]
sources:
  - citation: "PCMDI cmip6-cmor-tables, Tables/CMIP6_Ofx.json and CMIP6_Omon.json, data_specs_version 01.00.33 (18 Nov 2020)."
    locator: "variable_entry areacello, volcello, deptho, sftof, thkcello, masscello"
    url: https://github.com/PCMDI/cmip6-cmor-tables
  - citation: "Griffies, S. M., et al. (2016). OMIP contribution to CMIP6: experimental and diagnostic protocol for the physical component of the Ocean Model Intercomparison Project. Geosci. Model Dev., 9, 3231–3296."
    locator: "Appendix B (pp. 3242–3244); Appendix G3, G5, G7 (pp. 3251–3253); Table G1 (p. 3252)"
    url: https://doi.org/10.5194/gmd-9-3231-2016
---

## Key points

- **areacello**: "Horizontal area of ocean grid cells" (m2, cell_methods
  `area: sum`). This is the area of a **tracer** cell, and it "should be zero or
  missing over land". All CMIP6 ocean models have fixed horizontal grids, so
  the area is static. [CMOR Ofx; OMIP-phys B, G7]
- **volcello**: "Ocean Grid-Cell Volume" (m3, `area: sum`). This is an
  extensive field. It exists as a fixed field (Ofx) and as a time-varying one
  (Omon). [CMOR Ofx, Omon]
- **deptho**: "sea floor depth for present day relative to z=0 geoid.
  Reported as missing for land grid cells." It is the bottom "from the
  perspective of the model tracer fields". [CMOR Ofx; OMIP-phys G5]
- **sftof**: "Percentage of horizontal area occupied by ocean" (%). [CMOR Ofx]
- **thkcello** (cell thickness, m) and **masscello** (mass per area, kg m-2)
  are part of the request. For models with static cells,
  masscello = rhozero × thkcello, where rhozero is the model's Boussinesq
  reference density. That density is not the same in all models: OMIP
  suggests 1035 kg m-3, and early GFDL models used 1000 kg m-3. Where cell
  thickness varies in time, both are monthly fields. [CMOR Ofx/Omon; OMIP-phys
  B1–B2, G3]
- For full cells, thickness can be taken as the difference of the layer
  depth bounds. That fails at partial top and bottom cells. [OMIP-phys B1.1]

## Implications for this project

- None of the three models' local archive has `thkcello`, `masscello` or `rhozero`.
  The layer thickness that respects partial bottom cells is
  **volcello / areacello**. The `lev_bnds` difference is right only for full
  cells. Checked on the local historical files (January 2000, cells where
  `thetao` is valid):
  - **GFDL-ESM4**: area × layer-bound thickness gives 1.53e18 m3 against
    1.33e18 m3 from `volcello`, and 14 % of cells differ by more than 1 %.
    Its remapped z levels carry values in bottom cells that are only partly
    water.
  - **MPI-ESM1-2-LR**: partial bottom steps differ by up to ±50 % in 3.3 % of
    cells.
  - **CESM2**: the two agree within 0.6 %.
  `cell_thickness_like` in `src/tools/data_loader.py` returns volcello /
  areacello for CMIP6 fields. `identify_oxygen_minimum_zone` and the depth
  averages use it.
- `deptho` is model bathymetry, not true bathymetry. Deciding "is this level
  above the model sea floor" should use `deptho` or the fill pattern of the
  field itself, not an external topography.
- GFDL-ESM4's `volcello` varies in time because its native vertical coordinate
  moves (see `cmip6-grid-labels` and `configs/datasets.yaml`). The loaders
  already handle the time axis.
- On a curvilinear grid the model's own `areacello` is the only correct area
  weight, and the published `volcello` alone accounts for the partial cells at
  the sea floor. `load_areas` and `load_volumes_for_cmip6` use them where the
  model publishes them and record the origin in `attrs['area_source']` and
  `attrs['volume_source']`, which is worth reporting when they feed an
  integral.
- Total ocean area differs between models by a per cent or so, because a
  coarser grid resolves the coastline differently. That is a property of the
  model, not an error to correct: state the area each integral was taken over
  rather than rescaling one model to match another.

## Source excerpts

> "masscello = rhozero*thickcello, where thickcello is static cell thickness and
> rhozero is constant Boussinesq reference density. More generally, masscello
> is time dependent and reported as part of Omon." (CMOR Ofx, masscello)
