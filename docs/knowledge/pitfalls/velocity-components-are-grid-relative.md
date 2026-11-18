---
id: velocity-components-are-grid-relative
title: "CMIP6 uo/vo/umo/vmo are grid-relative (x-ward/y-ward), not eastward/northward, on curvilinear grids"
type: pitfall
summary: "uo/vo are \"x-ward/y-ward\" components along native grid lines (standard names sea_water_x/y_velocity) and sit at staggered C-grid positions; on the CESM2 displaced-pole and MPI bipolar grids they must be rotated before being treated as zonal/meridional, and OMIP advises against regridding them."
keywords: [uo, vo, wo, umo, vmo, wmo, velocity, 流速, current, 洋流, zonal, meridional, 纬向, 经向, x-ward, y-ward, rotation, 旋转, curvilinear, 曲线网格, C-grid, staggered, 交错网格, vector regridding, 矢量插值, upward]
datasets: [cmip6_historical, cmip6_ssp245]
models: [CESM2, MPI-ESM1-2-LR, GFDL-ESM4]
variables: [uo, vo, wo, umo, vmo, wmo]
tools: [regrid_horizontal, load_regular_data]
severity: must
related: [overturning-from-mass-transport, omip-native-grid-vs-regridded, model-cesm2-configuration, model-mpi-esm1-2-lr-configuration]
sources:
  - citation: "PCMDI cmip6-cmor-tables, Tables/CMIP6_Omon.json, data_specs_version 01.00.33 (18 Nov 2020)."
    locator: "variable_entry uo, vo, wo, umo, vmo, wmo"
    url: https://github.com/PCMDI/cmip6-cmor-tables
  - citation: "Griffies, S. M., et al. (2016). OMIP contribution to CMIP6: experimental and diagnostic protocol for the physical component of the Ocean Model Intercomparison Project. Geosci. Model Dev., 9, 3231–3296."
    locator: "Appendix C3 (p. 3245); Appendix I2–I5 and Table I1 (pp. 3264–3266)"
    url: https://doi.org/10.5194/gmd-9-3231-2016
  - citation: "Danabasoglu, G., et al. (2020). The Community Earth System Model Version 2 (CESM2). J. Adv. Model. Earth Syst., 12, e2019MS001916."
    locator: "Section 2.2 (p. 4): grid North Pole displaced into Greenland"
    url: https://doi.org/10.1029/2019MS001916
  - citation: "Mauritsen, T., et al. (2019). Developments in the MPI-M Earth System Model version 1.2 (MPI-ESM1.2). J. Adv. Model. Earth Syst., 11, 998–1038."
    locator: "Section 4 (pp. 1010–1011): Arakawa-C grid, bipolar GR1.5 with poles under Greenland and Antarctica"
    url: https://doi.org/10.1029/2018MS001400
  - citation: "WCRP-CMIP CMIP6_CVs, CMIP6_grid_label.json, CV collection version 6.2.60.0."
    locator: "grid_label gn, gr"
    url: https://github.com/WCRP-CMIP/CMIP6_CVs
---

## Key points

- CMOR standard names are `sea_water_x_velocity` and `sea_water_y_velocity`,
  with comments "Prognostic **x-ward** / **y-ward** velocity component resolved
  by the model". `umo`/`vmo` are "X-ward / Y-ward mass transport from resolved
  and parameterized advective transport" in kg s-1. [CMOR Omon]
- `wo` is `upward_sea_water_velocity`, **positive upward**. `wmo` is the upward
  mass transport (kg s-1, `area: sum`). [CMOR Omon; OMIP-phys I4–I5]
- OMIP recommends native-grid sampling of vector components. Transforming
  discrete vectors is ambiguous because "the commonly used C-grid has
  horizontal velocity components sitting at distinct spatial positions".
  OMIP therefore "generally recommend[s] against regridding of vector components
  onto a spherical grid". Regridded vector fields "can have spurious divergence
  and curl". [OMIP-phys C3, Table I1]
- `uo`/`vo` hold the resolved (Eulerian) velocity only. They exclude
  parameterised eddy-induced velocity, which `umo`/`vmo` include. [OMIP-phys I1–I3]
- Grids: CESM2 (POP2) has its grid North Pole displaced into Greenland. MPIOM
  GR1.5 is bipolar, with poles under Greenland and Antarctica, on an Arakawa-C
  grid. GFDL-ESM4's archived fields are regridded to a regular 1° grid (`gr`).
  [CESM2 §2.2; MPI §4; CVs grid_label]

## Implications for this project

- On CESM2 and MPI-ESM1-2-LR, the x/y directions depart from east/north
  (increasingly toward the displaced poles, i.e. the northern North Atlantic
  and Arctic). Before plotting currents, computing a zonal-mean zonal
  velocity or a northward transport, rotate (u, v) with the local grid angle
  derived from the 2-D lat/lon of adjacent cells. State that the rotation was
  applied.
- `uo`/`vo` sit at **velocity points, not tracer points**. Checked on the local
  historical files: the 2-D lat/lon stored with `uo`/`vo` differ from those of
  `thetao` on the same (j, i) index.
  - CESM2: `uo` and `vo` share one set of points, offset by a median of 0.22°
    in latitude and 0.56° in longitude, i.e. cell corners.
  - MPI-ESM1-2-LR: `uo` is offset mainly in longitude (median 0.62°) and `vo`
    mainly in latitude (median 0.32°), which is the C-grid face staggering.

  Masks, `areacello` and region selections built for tracers therefore do not
  apply directly to velocities. Use the velocity variable's own lat/lon, or
  average the velocities to tracer points first.
- For regridding or plotting, first rotate to east/north on the native grid,
  then interpolate each component (higher than linear order, per OMIP). Do not
  compute divergence or vorticity from the regridded field.
- Tool behaviour: the loader tags `uo`/`vo`/`umo`/`vmo` with
  `grid_point="velocity"` and `vector_component="x"/"y"`.
  - `regrid_horizontal` refuses these components on CESM2 and MPI-ESM1-2-LR.
  - `areas_like` and the area-weighted means still work, but record that
    tracer-cell areas were applied half a cell away. Those areas are adequate
    as averaging weights, not as face widths for transports.
  - No tool rotates vectors. Rotation must be done explicitly.

## Common mistakes

- Treating CESM2 `vo` as northward velocity in the North Atlantic.
- Regridding `uo` and `vo` independently with a conservative scalar method.

## Source excerpts

> "We therefore generally recommend against regridding of vector components
> onto a spherical grid." (OMIP-phys, C3)
