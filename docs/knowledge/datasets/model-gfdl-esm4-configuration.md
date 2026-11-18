---
id: model-gfdl-esm4-configuration
title: "GFDL-ESM4 ocean (MOM6) and biogeochemistry (COBALTv2) configuration"
type: dataset
summary: "GFDL-ESM4 (described in the paper as ESM4.1) runs MOM6 at nominal 0.5° on a tripolar C-grid with 75 hybrid z*-isopycnal layers and a 2 m top layer; its CMIP6 ocean output is conservatively regridded to a 1°×1° grid (gr); COBALTv2 has 33 tracers, fixed N:P per plankton type and mocsy carbonate chemistry."
keywords: [GFDL-ESM4, ESM4.1, MOM6, OM4, COBALTv2, tripolar, 三极网格, hybrid coordinate, 混合坐标, isopycnal, regridded, gr, conserve_order1, mocsy, fixed N:P, 固定化学计量, ocean model, 海洋模式]
datasets: [cmip6_historical, cmip6_ssp245]
models: [GFDL-ESM4]
variables: [thetao, so, o2, no3, po4, dissic, talk, ph, volcello]
tools: []
severity: info
related: [model-gfdl-esm4-reported-biases, cmip6-grid-labels, bgc-model-structure-comparison, cmip6-carbonate-system-definitions]
sources:
  - citation: "Dunne, J. P., et al. (2020). The GFDL Earth System Model Version 4.1 (GFDL-ESM 4.1): Overall coupled model description and simulation characteristics. J. Adv. Model. Earth Syst., 12, e2019MS002015."
    locator: "Sections 2.5–2.8 (pp. 7–8); Section 4.1 (p. 14); data citation https://doi.org/10.22033/ESGF/CMIP6.1407"
    url: https://doi.org/10.1029/2019MS002015
  - citation: "WCRP-CMIP CMIP6_CVs, CMIP6_source_id.json, CV collection version 6.2.60.0."
    locator: "source_id GFDL-ESM4 → model_component.ocean / ocnBgchem"
    url: https://github.com/WCRP-CMIP/CMIP6_CVs
---

## Key points

**Naming.** The paper describes "ESM4.1". Its CMIP6 output is published as
source_id `GFDL-ESM4` (data DOI 10.22033/ESGF/CMIP6.1407). [GFDL; CVs]

**Ocean physics (MOM6 / OM4p5)** [GFDL §2.5; CVs]
- Nominal 1/2° (720 × 576) on an Arakawa C grid. Meridional resolution is
  refined to 0.26° between 5°N and 5°S. The grid is tripolar north of 65°N
  (nominal 0.21°) and fixed at 0.19° south of 68°S.
- **75 hybrid depth–isopycnal layers**. The CMIP6 registry gives "top grid cell 0-2 m".
- A mesoscale eddy (MEKE) parameterisation and enhanced submesoscale
  restratification were needed at 0.5° to fix mid-latitude SST cold biases
  (list of development challenges, p. 10).
- Tracers are conserved to numerical precision. Salt drift over 900 years of
  piControl is about 1 part in 10^6. [GFDL §2.8]

**Archived grid (checked in the local files).** The global `grid` attribute
reads "ocean data regridded from tripolar - nominal 0.5 deg latitude/longitude
to 180,360; interpolation method: conserve_order1", with grid_label `gr`. The
vertical axis is 35 z levels (see `DATA_METADATA.md`).

**Biogeochemistry (COBALTv2)** [GFDL §2.6]
- 33 tracers covering C, alkalinity, O2, N, P, Fe, Si, calcite and lithogenic
  material. The plankton are small, large (diatom and non-diatom) and
  diazotrophic phytoplankton, three zooplankton groups, bacteria and three
  DOM labilities.
- "COBALTv2's reliance on **fixed nitrogen to phosphorus ratios** for plankton
  function types" is noted as the main simplification relative to TOPAZ.
- Carbonate chemistry uses **mocsy 2.0** (Lee et al. 2010 boron, Millero 2010
  K1/K2, Dickson & Riley 1979 Kf). The improved solver removes the unrealistically
  high pH that earlier GFDL models produced near river and ice-melt freshwater.
- O2 gas exchange uses a fixed atmospheric O2 concentration. [GFDL §2.7]

**Spin-up.** The ocean was initialised from observations, followed by a
400-year prototype spin-up and a 900-year final piControl. The historical run
began after about 500 years. [GFDL §4.1, §4.11]

## Implications for this project

- The archived fields are **already regridded**, horizontally (first-order
  conservative, to 1°) and vertically (from 75 hybrid layers to 35 z levels).
  Their land–sea mask is the regridded one. They can be selected by lat/lon
  directly, but budgets and integrals from `gr` fields are not exactly
  conservative (see `omip-native-grid-vs-regridded`).
- `volcello` varies in time because native layer thicknesses move with the
  hybrid coordinate. Loaders handle this (see `configs/datasets.yaml`).

## Source excerpts

> "ESM4.1 uses the same 75 vertical hybrid depth‐isopycnal layers as CM4.0."
> (GFDL, §2.5)
