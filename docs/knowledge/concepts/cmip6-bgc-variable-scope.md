---
id: cmip6-bgc-variable-scope
title: "What CMIP6 chl, phycos, dissoc and dfe actually include"
type: concept
summary: "chl sums all modelled phytoplankton groups, phycos is surface-only phytoplankton carbon, dissoc counts only explicitly modelled DOC pools (usually not the ~40 µM refractory pool), and dfe includes free and complexed dissolved iron."
keywords: [chl, chlorophyll, 叶绿素, phycos, phytoplankton carbon, 浮游植物碳, dissoc, DOC, dissolved organic carbon, 溶解有机碳, refractory DOC, 难降解, dfe, dissolved iron, 溶解铁, ligand]
datasets: [cmip6_historical, cmip6_ssp245]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [chl, chlos, phycos, dissoc, intdoc, dfe]
tools: []
severity: should
related: [cmip6-bgc-units-per-volume, cmip6-integrated-carbon-content, bgc-model-structure-comparison]
sources:
  - citation: "PCMDI cmip6-cmor-tables, Tables/CMIP6_Omon.json, data_specs_version 01.00.33 (18 Nov 2020)."
    locator: "variable_entry chl, chlos, phycos, dissoc, intdoc, dfe"
    url: https://github.com/PCMDI/cmip6-cmor-tables
  - citation: "Orr, J. C., et al. (2017). Biogeochemical protocols and diagnostics for the CMIP6 Ocean Model Intercomparison Project (OMIP). Geosci. Model Dev., 10, 2169–2199."
    locator: "Table 6 footnotes d–e (p. 2184); Table 11 footnote d (p. 2189)"
    url: https://doi.org/10.5194/gmd-10-2169-2017
  - citation: "Danabasoglu, G., et al. (2020). The Community Earth System Model Version 2 (CESM2). J. Adv. Model. Earth Syst., 12, e2019MS001916."
    locator: "Section 2.2 (p. 5), MARBL DOM pools"
    url: https://doi.org/10.1029/2019MS001916
  - citation: "Dunne, J. P., et al. (2020). The GFDL Earth System Model Version 4.1 (GFDL-ESM 4.1). J. Adv. Model. Earth Syst., 12, e2019MS002015."
    locator: "Section 2.6 (p. 7), COBALTv2 DOM labilities"
    url: https://doi.org/10.1029/2019MS002015
  - citation: "Mauritsen, T., et al. (2019). Developments in the MPI-M Earth System Model version 1.2 (MPI-ESM1.2). J. Adv. Model. Earth Syst., 11, 998–1038."
    locator: "Section 5.1 (p. 1013), DOM pools combined"
    url: https://doi.org/10.1029/2018MS001400
---

## Key points

- **chl**: "Sum of chlorophyll from all phytoplankton group concentrations. In
  most models this is equal to chldiat+chlmisc", in kg m-3 on the full
  3-D grid. [CMOR Omon; OMIP-BGC Table 6 fn e]
- **phycos**: "sum of phytoplankton organic carbon component concentrations
  at the sea surface", in mol m-3 on (lon, lat, time). It is surface only.
  [CMOR Omon]
- **dissoc**: "Sum of dissolved carbon component concentrations explicitly
  represented (i.e. not ~40 uM refractory unless explicit)". `intdoc` likewise
  covers "explicit pools only". [CMOR Omon; OMIP-BGC Table 11 fn d]
- **dfe**: "Dissolved iron in sea water, including both Fe2+ and Fe3+ ions (but
  not particulate detrital iron)". OMIP adds that modelled dissolved iron
  "includes all simulated dissolved species, both free and organically
  complexed". [CMOR Omon; OMIP-BGC Table 6 fn d]

## Implications for this project

- **Model DOC is not total DOC.** Observed bulk DOC includes a large refractory
  background. Model `dissoc` counts only the pools the model carries. The
  models differ here: the CESM2 paper states MARBL has "both semilabile and
  refractory dissolved organic material pools", while GFDL COBALTv2 has "three
  labilities of dissolved organic matter" and MPI HAMOCC6 merged its DOM into
  one pool (see the model configuration fragments). Cross-model DOC
  differences partly reflect which pools exist. Model-versus-observation DOC
  differences are dominated by this definitional gap.
- `chl` and `phycos` are both totals over each model's phytoplankton groups,
  and the number and type of groups differ between models. Differences therefore
  include model structure, not only simulated biomass.
- In CESM2, `chl` covers only the top 150 m (15 levels; see
  `configs/datasets.yaml`). Any depth integral of CESM2 chlorophyll is an
  integral over 0–150 m.
- In CESM2 historical, `phycos` comes from a different ensemble member (r6)
  than the other variables (see `configs/datasets.yaml`).

## Source excerpts

> "Sum of dissolved carbon component concentrations explicitly represented
> (i.e. not ~40 uM refractory unless explicit)" (CMOR Omon, dissoc)
