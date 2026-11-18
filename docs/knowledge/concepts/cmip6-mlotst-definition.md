---
id: cmip6-mlotst-definition
title: "How CMIP6 mlotst is defined, and why an offline MLD will not match it exactly"
type: concept
summary: "mlotst uses a surface-referenced sigma-t / buoyancy criterion (ΔB = 0.0003 m s-2 ≈ Δρ 0.03 kg m-3) evaluated at local pressure, is meant to be computed online every time step, and is then time-averaged over the month."
keywords: [mlotst, mixed layer depth, 混合层深度, MLD, sigma-t, density threshold, 密度阈值, 0.03, buoyancy, reference depth, 参考深度, de Boyer Montégut, Levitus, mlotstmax, mlotstmin]
datasets: [cmip6_historical, seanoe_mld, woa23]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [mlotst, mlotstmax, mlotstmin, thetao, so]
tools: [compute_mixed_layer_depth, compute_potential_density_anomaly]
severity: must
related: [monthly-means-are-time-averages, cmip6-thetao-potential-temperature, cmip6-so-practical-salinity]
sources:
  - citation: "Griffies, S. M., et al. (2016). OMIP contribution to CMIP6: experimental and diagnostic protocol for the physical component of the Ocean Model Intercomparison Project. Geosci. Model Dev., 9, 3231–3296."
    locator: "Appendix H24 (pp. 3262–3263); Appendix A1 (p. 3239)"
    url: https://doi.org/10.5194/gmd-9-3231-2016
  - citation: "PCMDI cmip6-cmor-tables, Tables/CMIP6_Omon.json, data_specs_version 01.00.33 (18 Nov 2020)."
    locator: "variable_entry mlotst, mlotstmax"
    url: https://github.com/PCMDI/cmip6-cmor-tables
  - citation: "Dunne, J. P., et al. (2020). The GFDL Earth System Model Version 4.1 (GFDL-ESM 4.1). J. Adv. Model. Earth Syst., 12, e2019MS002015."
    locator: "Section 4.9, Figs. 16–17 captions (pp. 25–26)"
    url: https://doi.org/10.1029/2019MS002015
---

## Key points

- CMOR: `mlotst` = "Ocean Mixed Layer Thickness Defined by Sigma T", standard
  name `ocean_mixed_layer_thickness_defined_by_sigma_t`, comment "Sigma T is
  potential density referenced to ocean surface", cell_methods `time: mean`.
  `mlotstmax` and `mlotstmin` are the monthly maximum and minimum of
  instantaneous MLD. [CMOR Omon; OMIP-phys A1]
- The OMIP criterion follows Levitus (1982). A **surface parcel** (S and Θ of
  level k = 1) is displaced down to level k **at that level's in-situ pressure
  p(k)**. Its density is compared with the local density at p(k):
  - ρ_displaced = ρ[S(1), Θ(1), p(k)] and ρ_local = ρ[S(k), Θ(k), p(k)]
  - δB = −g (ρ_displaced − ρ_local) / ρ_local
  - The MLD is the first depth where δB > ΔB_crit = 0.0003 m s-2, linearly
    interpolated between that depth and the one above it. [OMIP-phys H24.2,
    eqs. H39–H42]
- With g = 9.8 m s-2 and ρ ≈ 1035 kg m-3, ΔB_crit corresponds to
  **Δρ_crit = 0.03 kg m-3**, the de Boyer Montégut et al. (2004) value. Some
  studies use 0.125 kg m-3, which gives deeper mixed layers. The choice "should be
  compatible across models and observations". [OMIP-phys H24.2]
- The MLD asked for is **not** the KPP boundary-layer depth, and not the bulk
  mixed layer of isopycnal models. [OMIP-phys H24.2]
- OMIP recommends computing `mlotst` **online, every model time step**, "to
  avoid aliasing". The monthly value is the time mean of those instantaneous
  depths. [OMIP-phys H24.2, A1]
- Observational MLD climatologies often use a 10 m reference depth. The GFDL paper,
  for example, evaluates its model against WOA13 with "the 0.03 density threshold
  relative to the upper 10 m". [GFDL §4.9]

## Implications for this project

- To compare an offline MLD with a model's own `mlotst`:
  - pass `ref_depth` equal to the model's shallowest level to
    `compute_mixed_layer_depth`, since `mlotst` is referenced to the ocean
    surface, not to 10 m. That removes the largest definitional difference,
    but not all of it: the result comes close to `mlotst`, never identical.
  - compute the MLD of each monthly profile and then average the monthly MLDs.
    Never average temperature and salinity first: the MLD of an averaged
    profile is not the average of the MLDs and is systematically shallower,
    because averaging smooths the pycnocline away.
  - leave out the columns that never reach the threshold (below) and report how
    much area that excludes.
- An MLD recomputed offline from monthly-mean `thetao`/`so` differs from the
  model's own `mlotst` for three reasons, even with the same 0.03 threshold:
  1. The reference level: OMIP uses the model's top level, while the default
     `ref_depth` of `compute_mixed_layer_depth` is 10 m.
  2. The density form: OMIP uses a displaced-parcel comparison at local pressure,
     the offline operator a σ0 difference.
  3. Time averaging: OMIP averages instantaneous MLDs, the offline operator uses
     the MLD of averaged profiles.
  State these differences when comparing the two. Do not treat the residual
  as model error.
- Because of point 3, an offline MLD from monthly means is expected to come
  out **shallower** than `mlotst`, even with the same threshold and reference
  level. A monthly mean is itself an averaged profile.
- **Columns where the threshold is never reached** are encoded differently
  in the three models' `mlotst`. This was checked on March 1990 against cells
  where the offline operator finds no crossing: GFDL-ESM4 and MPI-ESM1-2-LR
  report roughly the full column depth (median `mlotst/deptho` 0.96 and 0.85),
  while CESM2 reports 0 m in 1077 of 2486 such cells, mostly on Arctic
  shelves. Exclude these columns from any operator-versus-`mlotst` comparison.
  No single fill value treats all three models fairly.
- GFDL-ESM4 computes `mlotst` on its native 75 hybrid layers. The archived
  `thetao`/`so` it would be compared with were remapped to 35 z levels, which
  gives that model an extra resolution difference.
- The defaults of `compute_mixed_layer_depth` (0.03 kg m-3 from 10 m) are the
  de Boyer Montégut criterion that observational MLD climatologies use, so an
  MLD computed with them is directly comparable with those climatologies.
  Comparisons with SEANOE (0.03 kg m-3, 5 m or 10 m reference; see
  `configs/datasets.yaml`) mix criterion and reference-depth differences with
  real bias. Say which ones apply.

## Common mistakes

- Calling `mlotst` a "0.03 kg m-3 from 10 m" MLD.
- Comparing a 0.125 kg m-3 MLD with `mlotst`.

## Source excerpts

> "we recommend computing mlotst online in order to avoid aliasing."
> (OMIP-phys, H24.2)

> "Notably, the boundary layer depth is not what is asked for with the OMIP
> mixed layer depth diagnostic." (OMIP-phys, H24.2)
