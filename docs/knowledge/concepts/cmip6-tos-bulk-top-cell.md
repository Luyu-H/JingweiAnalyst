---
id: cmip6-tos-bulk-top-cell
title: "CMIP6 tos is the bulk temperature of the top model cell, including under sea ice"
type: concept
summary: "tos is the upper-boundary (top-cell bulk) potential temperature of the liquid ocean, also under sea ice and ice shelves; it is not a skin temperature, and top-cell thickness differs by model."
keywords: [tos, SST, sea surface temperature, 海表温度, top level, 表层, skin temperature, sea ice, 海冰, bulk, thetao first level, 表层温度替代]
datasets: [cmip6_historical, cmip6_ssp245, ersst, woa23]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [tos, thetao]
tools: [load_regular_data, compute_enso_nino_index, compute_pdo_index, compute_anomaly_vs_woa23, compute_anomaly_vs_own_climatology, regrid_vertical]
severity: should
related: [cmip6-thetao-potential-temperature, derive-unavailable-diagnostics]
sources:
  - citation: "Griffies, S. M., et al. (2016). OMIP contribution to CMIP6: experimental and diagnostic protocol for the physical component of the Ocean Model Intercomparison Project. Geosci. Model Dev., 9, 3231–3296."
    locator: "Appendix H13 (p. 3260)"
    url: https://doi.org/10.5194/gmd-9-3231-2016
  - citation: "PCMDI cmip6-cmor-tables, Tables/CMIP6_Omon.json, data_specs_version 01.00.33 (18 Nov 2020)."
    locator: "variable_entry tos"
    url: https://github.com/PCMDI/cmip6-cmor-tables
  - citation: "WCRP-CMIP CMIP6_CVs, CMIP6_source_id.json, CV collection version 6.2.60.0."
    locator: "source_id CESM2 / GFDL-ESM4 / MPI-ESM1-2-LR → model_component.ocean.description"
    url: https://github.com/WCRP-CMIP/CMIP6_CVs
  - citation: "Danabasoglu, G., et al. (2020). The Community Earth System Model Version 2 (CESM2). J. Adv. Model. Earth Syst., 12, e2019MS001916."
    locator: "Section 7, p. 15 (SST defined as 0–10 m average)"
    url: https://doi.org/10.1029/2019MS001916
---

## Key points

- `tos` is the "temperature of upper boundary of the liquid ocean, including
  temperatures below sea-ice and floating ice shelves" (CMOR). Under sea ice it
  is the water temperature just below the ice, not the ice-surface temperature.
  [CMOR Omon; OMIP-phys H13]
- Model surface tracers are **bulk values averaged over the top grid cell**,
  which is generally at least a metre thick. They are not skin values. [OMIP-phys H13]
- "In nearly all ocean climate models, the surface temperature is given by the
  prognostic temperature in the upper-most ocean grid cell." [OMIP-phys H13]
- `tos` is surface potential temperature, which equals surface in-situ
  temperature and differs from surface Conservative Temperature. [OMIP-phys H13]
- Top-cell thickness of the three native ocean grids, according to the CMIP6 source
  registry [CVs source_id]:

  | model | ocean component | top grid cell |
  |---|---|---|
  | CESM2 | POP2, 60 levels | 0–10 m |
  | GFDL-ESM4 | GFDL-OM4p5 (MOM6), 75 levels | 0–2 m |
  | MPI-ESM1-2-LR | MPIOM1.63, 40 levels | 0–12 m |

  The CESM2 paper itself evaluates SST as "0‐ to 10‐m average in the ocean"
  against ERSSTv5. [CESM2 §7]

## Implications for this project

- Using the first level of `thetao` as a substitute for `tos` (needed for
  ssp245, where `tos` is missing; see `configs/datasets.yaml`) reproduces the
  model's own SST definition, because tos is the top-cell prognostic
  temperature. Minor differences are possible where the archived vertical grid
  differs from the native grid. GFDL-ESM4's archived `thetao` has been
  regridded to 35 z levels, the first centred at 2.5 m. Checked on GFDL-ESM4,
  1990–1999: a Niño3.4 index from the top level of `thetao` matches the one
  from `tos` to an RMS of 0.005 °C. A 0–10 m thickness-weighted average
  (two levels) departs by 0.031 °C. `compute_enso_nino_index` and
  `compute_pdo_index` therefore take the top level by default
  (`depth_level=None`), the like-for-like substitute where `tos` is missing.
  A `depth_level` of 10 m instead averages two levels in GFDL-ESM4 but one in
  CESM2 and MPI-ESM1-2-LR, so it is not the same quantity across models.
- A model's first level stands for its surface. A value asked for at 0 m,
  above every model's shallowest level centre, is that shallowest level's value
  (`regrid_vertical` with `fill_surface=True`).
- An ENSO index from a model is its anomaly from its **own** climatology.
  Subtracting an observed climatology puts the model's mean-state bias into
  the index, which is then not an ENSO index (`compute_enso_nino_index` with
  `climatology='woa23'` on model output does exactly that).
- Sea-ice-covered cells carry the under-ice water temperature (near freezing).
  They are not missing. A global-mean SST therefore includes polar values near
  −1.8 °C. Observational products that also report under-ice SST near freezing
  are comparable. If a product masks sea ice, the ocean masks must be made
  consistent before comparing.

## Source excerpts

> "Surface tracer fields produced from a climate model generally do not
> correspond to skin properties. Rather, they are bulk properties averaged over
> the top grid cell, which is generally no less than a metre thick."
> (OMIP-phys, H13)
