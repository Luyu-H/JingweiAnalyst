---
id: overturning-from-mass-transport
title: "Meridional overturning streamfunction (AMOC) from CMIP6 transports"
type: method
summary: "The overturning streamfunction is Ψ(y,z) = −∫∫ ρ v† dz dx in kg s-1, reported for Atlantic–Arctic, Indo-Pacific and global basins; msftmz is not in the local archive, so it must be built from vmo (MPI only) or vo, along native grid rows, with mass Sv = 1e9 kg s-1."
keywords: [AMOC, overturning, 经向翻转环流, streamfunction, 流函数, msftmz, msftyz, vmo, vo, Sverdrup, Sv, RAPID, 26.5N, NADW, AABW, sigma2000, density space]
datasets: [cmip6_historical, cmip6_ssp245]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [msftmz, msftyz, vmo, vo, umo, uo]
tools: [compute_basic_stat, select_basin, load_basin_mask]
severity: should
related: [velocity-components-are-grid-relative, derive-unavailable-diagnostics, basin-and-region-definitions, model-cesm2-reported-biases, model-gfdl-esm4-reported-biases]
sources:
  - citation: "Griffies, S. M., et al. (2016). OMIP contribution to CMIP6: experimental and diagnostic protocol for the physical component of the Ocean Model Intercomparison Project. Geosci. Model Dev., 9, 3231–3296."
    locator: "Appendix C4 (pp. 3245–3246); Appendix I1, Eqs. I1–I4 (p. 3264); I3 (p. 3266); I6, Eq. I6 (pp. 3266–3267)"
    url: https://doi.org/10.5194/gmd-9-3231-2016
  - citation: "PCMDI cmip6-cmor-tables, Tables/CMIP6_Omon.json and CMIP6_coordinate.json, data_specs_version 01.00.33 (18 Nov 2020)."
    locator: "variable_entry msftmz, msftyz, vmo, umo; axis_entry basin, rho"
    url: https://github.com/PCMDI/cmip6-cmor-tables
  - citation: "Danabasoglu, G., et al. (2020). The Community Earth System Model Version 2 (CESM2). J. Adv. Model. Earth Syst., 12, e2019MS001916."
    locator: "Section 9 (pp. 22–23)"
    url: https://doi.org/10.1029/2019MS001916
  - citation: "Dunne, J. P., et al. (2020). The GFDL Earth System Model Version 4.1 (GFDL-ESM 4.1). J. Adv. Model. Earth Syst., 12, e2019MS002015."
    locator: "Section 4.1 (p. 14); Section 4.10, Fig. 21 caption (pp. 30–31)"
    url: https://doi.org/10.1029/2019MS002015
---

## Key points

- Definition: Ψ(y, s, t) = − ∫_{xa}^{xb} dx ∫_{−H}^{z(s)} ρ v† dz, where v†
  is the **residual-mean** meridional velocity (resolved plus parameterised
  eddy-induced). s is depth or pressure, or potential density σ2000. [OMIP-phys I6, Eq. I6]
- Transports are archived as **mass transport, kg s-1**. 1 mass Sv = 10^9 kg s-1
  and 1 volume Sv = 10^6 m3 s-1. [OMIP-phys I1, Eqs. I3–I4]
- The three physically meaningful basins are **Atlantic–Arctic, Indo-Pacific
  and global**. "To separate the Indian and Pacific oceans is not sensible,
  since there is no meridional boundary separating these basins."
  The CMOR `basin` axis values are `atlantic_arctic_ocean`,
  `indian_pacific_ocean` and `global_ocean`. [OMIP-phys I6; CMOR coordinate]
- On non-spherical grids, integrate along native grid lines that approximate
  latitude circles (a "zig-zag" path) and report the result against latitude.
  The ŷ-ward native streamfunction `msftyz` "closely approximates the spherical
  meridional streamfunction" on tripolar and displaced-pole grids.
  [OMIP-phys C4, I6]
- `vmo`/`umo` = mass transport through cell faces from "resolved and
  parameterized advective transport". `vo`/`uo` are the resolved velocity only.
  [CMOR Omon; OMIP-phys I2–I3]
- Reported model values:
  - CESM2, 1995–2014, Eulerian-mean (resolved) AMOC: NADW maximum about 24 Sv.
    RAPID at 26.5°N gives 17.0 ± 3.3 Sv (April 2004 to February 2017). The CESM2
    maximum at 26.5°N is 2–3 Sv above RAPID. AABW cell about 0.5 Sv against an
    observed 1–5 Sv. [CESM2 §9]
  - GFDL-ESM4.1 piControl: maximum North Atlantic overturning 21.6 ± 1.4 Sv.
    In the historical run the AMOC at 24°N is shallower than observed and ends
    at 3000 m. [GFDL §4.1, §4.10]
  - GFDL computes its overturning "just from the Eulerian mean transport", and
    measures the strength as the maximum minus the surface value. [GFDL Fig. 21
    caption]

## Implications for this project

- `msftmz`/`msftyz` are **not** in the local archive. Of the transports, only
  MPI-ESM1-2-LR has `vmo`. CESM2 has `vo` (resolved only), and GFDL-ESM4 has
  neither (see `configs/datasets.yaml`). An AMOC comparison can therefore cover
  at most MPI (from `vmo`) and CESM2 (from `vo` × face area × ρ_o).
- `vo` on a curvilinear grid is **y-ward, not northward** (see
  `velocity-components-are-grid-relative`). Sum along a grid row (constant j),
  mask the Atlantic–Arctic basin, integrate cumulatively in the vertical and
  label rows with their mean latitude. Converting `vo` needs the face widths,
  which are not in the archive; approximating them from cell geometry adds
  error, which should be stated.
- A streamfunction from `vo` is Eulerian-mean. One from `vmo` includes
  parameterised eddies. Say which one was used.

## Source excerpts

> "Atlantic–Arctic, Indian–Pacific, and World Ocean are the only three
> physically relevant partitions available." (OMIP-phys, I6)
