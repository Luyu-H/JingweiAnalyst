---
id: model-cesm2-configuration
title: "CESM2 ocean (POP2) and biogeochemistry (MARBL) configuration"
type: dataset
summary: "CESM2 uses POP2 on a nominal 1° displaced-pole grid (pole in Greenland, 0.27° meridional at the equator), 60 levels with 10 m spacing in the top 160 m, KPP mixing, and MARBL with N/P/Si/Fe co-limitation, three explicit phytoplankton groups and variable C:P."
keywords: [CESM2, POP2, MARBL, NCAR, displaced pole, 位移极点网格, gx1v7, KPP, ocean model, 海洋模式, biogeochemistry, 生物地球化学, phytoplankton, 浮游植物, variable stoichiometry, Caspian]
datasets: [cmip6_historical, cmip6_ssp245]
models: [CESM2]
variables: [thetao, so, no3, po4, si, dfe, chl, dissic, talk, ph, dissoc]
tools: []
severity: info
related: [model-cesm2-reported-biases, bgc-model-structure-comparison, velocity-components-are-grid-relative, historical-ssp245-member-continuity]
sources:
  - citation: "Danabasoglu, G., et al. (2020). The Community Earth System Model Version 2 (CESM2). J. Adv. Model. Earth Syst., 12, e2019MS001916."
    locator: "Section 2.2 (pp. 4–5); Section 4 (pp. 9–10)"
    url: https://doi.org/10.1029/2019MS001916
  - citation: "WCRP-CMIP CMIP6_CVs, CMIP6_source_id.json, CV collection version 6.2.60.0."
    locator: "source_id CESM2 → model_component.ocean / ocnBgchem"
    url: https://github.com/WCRP-CMIP/CMIP6_CVs
---

## Key points

**Ocean physics (POP2)** [CESM2 §2.2; CVs]
- A level-coordinate model (POP2), with KPP vertical mixing via CVMix, Langmuir
  mixing and a salinity-dependent freezing point. Shortwave absorption uses
  prognostic chlorophyll.
- Grid: spherical in the Southern Hemisphere. In the Northern Hemisphere the grid North Pole
  is **displaced into Greenland**. Nominal 1° with a uniform 1.125° zonal
  spacing. Meridional spacing is 0.27° at the equator, rising to 0.53° at 32°S
  and constant further south. In the Northern Hemisphere high latitudes it ranges from about 0.38°
  (north-western Atlantic) to 0.64° (north-western Pacific).
- 60 levels, maximum depth 5500 m. Spacing is a uniform **10 m in the upper
  160 m**, increasing to 250 m by 3500 m. The CMIP6 registry gives "320x384
  longitude/latitude; 60 levels; top grid cell 0-10 m".
- The **Caspian Sea is treated as a lake** in the land model and is not part of
  the ocean model.

**Biogeochemistry (MARBL, the updated BEC)** [CESM2 §2.2]
- Multiple nutrient co-limitation (**N, P, Si, Fe**). Three explicit
  phytoplankton functional groups (diatoms, diazotrophs, pico/nano) and one
  implicit group (calcifiers), with one zooplankton group.
- Fully prognostic carbonate chemistry. Sinking particles are implicit, with
  ballasting by dust, CaCO3 and Si.
- Updates relative to CESM1: sub-grid light variation, **variable C:P
  stoichiometry**, sea-floor burial balancing riverine input, semilabile **and
  refractory** DOM pools, prognostic NH3 emission and an explicit iron-binding
  ligand tracer.
- Riverine nutrients, carbon and alkalinity come from data sets. Inorganic N and P
  vary over 1900–2005 and are held constant outside that period.

**Initial state** [CESM2 §4]
- Ocean T/S come from PHC2, integrated about 1070 years before the PI control.
  BGC tracers were spun up offline with repeated forcing and a Newton–Krylov
  solver.
- The 11 historical members start from piControl years 501–871.

## Implications for this project

- Curvilinear grid: never select by lat/lon value (see `configs/datasets.yaml`).
  Velocities are grid-relative (see `velocity-components-are-grid-relative`).
- Fine equatorial meridional resolution (0.27°) means equatorial features are
  sharper than in the 1° observational grids. Regridding smooths them.
- The model domain has no Caspian Sea. Observational masks that include it
  must be reconciled.

## Source excerpts

> "The model uses spherical coordinates in the Southern Hemisphere. In the
> Northern Hemisphere, the grid North Pole is displaced into Greenland."
> (CESM2, §2.2)
