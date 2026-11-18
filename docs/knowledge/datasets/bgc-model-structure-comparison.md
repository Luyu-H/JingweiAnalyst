---
id: bgc-model-structure-comparison
title: "Side-by-side structure of the three ocean biogeochemistry models (MARBL, COBALTv2, HAMOCC6)"
type: dataset
summary: "The three models differ in nutrients carried, plankton groups, stoichiometry (CESM2 variable C:P, GFDL fixed N:P per plankton type, MPI fixed Redfield), DOM pools and N2 fixation, which shapes cross-model differences in nutrients, O2, chlorophyll and DOC independently of physics."
keywords: [biogeochemistry, 生物地球化学, MARBL, COBALTv2, HAMOCC6, stoichiometry, 化学计量, Redfield, N:P, C:P, nitrogen fixation, 固氮, plankton functional types, 浮游植物功能群, iron, 铁, DOM, model comparison, 多模式比较]
datasets: [cmip6_historical, cmip6_ssp245]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [no3, po4, si, dfe, o2, chl, phycos, dissoc, dissic, talk, ph]
tools: [compute_ratio, fit_multivariate_relation]
severity: info
related: [model-cesm2-configuration, model-gfdl-esm4-configuration, model-mpi-esm1-2-lr-configuration, cmip6-bgc-variable-scope]
sources:
  - citation: "Danabasoglu, G., et al. (2020). The Community Earth System Model Version 2 (CESM2). J. Adv. Model. Earth Syst., 12, e2019MS001916."
    locator: "Section 2.2 (pp. 4–5)"
    url: https://doi.org/10.1029/2019MS001916
  - citation: "Dunne, J. P., et al. (2020). The GFDL Earth System Model Version 4.1 (GFDL-ESM 4.1). J. Adv. Model. Earth Syst., 12, e2019MS002015."
    locator: "Sections 2.6–2.7 (pp. 7–8)"
    url: https://doi.org/10.1029/2019MS002015
  - citation: "Mauritsen, T., et al. (2019). Developments in the MPI-M Earth System Model version 1.2 (MPI-ESM1.2). J. Adv. Model. Earth Syst., 11, 998–1038."
    locator: "Section 5 (pp. 1011–1015)"
    url: https://doi.org/10.1029/2018MS001400
---

## Key points

| aspect | CESM2 (MARBL) | GFDL-ESM4 (COBALTv2) | MPI-ESM1-2-LR (HAMOCC6) |
|---|---|---|---|
| nutrients / limitation | N, P, Si, Fe co-limitation | C, alkalinity, O2, N, P, Fe, Si cycles (33 tracers) | limitation scheme not detailed in paper; nitrate, phosphate, iron and silicate appear (Fe limitation of cyanobacteria, weathering inputs) |
| phytoplankton | diatoms, diazotrophs, pico/nano (explicit); calcifiers (implicit) | small, large (diatom, non-diatom), diazotrophs | bulk phytoplankton + prognostic cyanobacteria |
| zooplankton | 1 group | 3 groups (+ bacteria) | – (not detailed in paper) |
| stoichiometry | **variable C:P** | **fixed N:P per plankton type** | **prescribed Redfield** in organic matter |
| DOM | semilabile and refractory pools | 3 labilities (incl. semilabile) | one combined pool |
| N2 fixation | diazotroph group | diazotroph group | prognostic cyanobacteria, ~40°S–40°N, ~82 Tg N yr-1 |
| O2 and remineralisation | – (not detailed in paper) | – (not detailed in paper) | Monod O2 limitation (10 µmol L-1), anaerobic on NO3, N2O, SO4 |
| sinking particles | implicit, mineral ballasting | temperature-dependent remineralisation of sinking organic matter | settling 3.5 → 80 m d-1 with depth |
| carbonate chemistry | fully prognostic | mocsy 2.0 | total pH scale, Dickson constants, P+Si alkalinity |
| iron | dust + black carbon deposition, ligand tracer | interactive dust/iron, iceberg iron | dust climatology (Mahowald et al. 2005) |

[CESM2 §2.2; GFDL §2.6–2.7; MPI §5]

## Implications for this project

- Cross-model spread in nutrient ratios, chlorophyll or DOC reflects model
  **structure** as well as circulation. For example, a fixed-Redfield model
  can only produce N:P departures through N-cycle and transport processes.
  Mention this when interpreting spread.
- "Not detailed in paper" means the description paper does not state it. It
  does not mean the feature is absent. Do not infer absence.

## Source excerpts

> "The most noteworthy simplification relative to TOPAZ is COBALTv2's
> reliance on fixed nitrogen to phosphorus ratios for plankton function
> types." (GFDL, §2.6)
