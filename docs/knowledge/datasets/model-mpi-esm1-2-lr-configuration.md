---
id: model-mpi-esm1-2-lr-configuration
title: "MPI-ESM1-2-LR ocean (MPIOM1.6) and biogeochemistry (HAMOCC6) configuration"
type: dataset
summary: "MPI-ESM1-2-LR runs MPIOM1.6.3 on the bipolar GR1.5 grid (poles under Greenland and Antarctica) with 40 z-levels (20 in the top 700 m), Boussinesq and hydrostatic; HAMOCC6 uses a fixed Redfield composition, prognostic cyanobacteria N2 fixation, O2-dependent remineralisation, and OMIP carbonate chemistry."
keywords: [MPI-ESM1-2-LR, MPI-ESM1.2, MPIOM, HAMOCC6, GR1.5, bipolar grid, 双极网格, Boussinesq, Redfield, 雷德菲尔德比, nitrogen fixation, 固氮, cyanobacteria, denitrification, 反硝化, remineralization, 再矿化, total pH scale, tuning, 调参]
datasets: [cmip6_historical, cmip6_ssp245]
models: [MPI-ESM1-2-LR]
variables: [thetao, so, o2, no3, po4, dissic, talk, ph, uo, vo, umo, vmo, wmo]
tools: []
severity: info
related: [bgc-model-structure-comparison, velocity-components-are-grid-relative, cmip6-carbonate-system-definitions, model-drift-in-deep-trends]
sources:
  - citation: "Mauritsen, T., et al. (2019). Developments in the MPI-M Earth System Model version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 998–1038."
    locator: "Section 2 and Table 1 (pp. 1000–1001); Section 4 (pp. 1010–1011); Section 5 (pp. 1011–1015); Sections 8.1–8.2 (pp. 1021–1022)"
    url: https://doi.org/10.1029/2018MS001400
  - citation: "WCRP-CMIP CMIP6_CVs, CMIP6_source_id.json, CV collection version 6.2.60.0."
    locator: "source_id MPI-ESM1-2-LR → model_component.ocean / ocnBgchem"
    url: https://github.com/WCRP-CMIP/CMIP6_CVs
---

## Key points

**Ocean physics (MPIOM1.6.3)** [MPI §4, Table 1; CVs]
- Arakawa-C grid, **z levels, hydrostatic and Boussinesq**. Isopycnal (Redi)
  mixing and GM eddy transport. Vertical mixing combines the
  Pacanowski–Philander Richardson-number scheme with wind-driven mixed-layer
  mixing.
- LR uses the **bipolar GR1.5 grid (~1.5°, ~150 km)**, "one grid pole under
  Greenland and one under Antarctica". This gives higher resolution in the
  deep-water formation regions and the Greenland–Scotland overflows.
- **40 levels, the first 20 within the top 700 m**. Partial-step bottom
  topography. The registry gives "256 x 220 longitude/latitude; 40 levels; top grid cell 0-12 m".
- A vertical-viscosity coding bug was fixed **in LR (not HR)**. The fix
  considerably reduced SST biases in upwelling regions and the Southern Ocean,
  and cooled the model by 0.2–0.3 K, which required retuning.
- The atmosphere and ocean couple **daily** in LR (hourly in HR). [Table 1]

**Biogeochemistry (HAMOCC6)** [MPI §5]
- **Prognostic cyanobacteria** (Trichodesmium-like) provide N2 fixation,
  restricted to roughly 40°S–40°N by temperature. Global N2 fixation is about
  82 Tg N yr-1 in LR, and export about 6.7 Gt C yr-1.
- The remineralisation composition follows a **prescribed Redfield ratio** in
  organic material. Remineralisation proceeds on O2, nitrate, N2O and sulphate,
  with a Monod-type O2 limitation (half-saturation 10 µmol O2 L-1) that makes
  the aerobic–anaerobic transition gradual. H2S is a new tracer in LR only.
- Detritus sinks at 3.5 m d-1 in the upper 100 m, increasing linearly to
  80 m d-1 at 6000 m (Martin-curve fit). Plankton dynamics now act over the
  whole water column, not only the top 100 m.
- OMIP-compliant carbonate chemistry: **total pH scale**, Dickson constants,
  and alkalinity including phosphate and silicate. Gas exchange follows
  Wanninkhof (2014). N deposition is added to NO3.
- The two DOM pools were merged into one.

**Spin-up and tuning** [MPI §8.1–8.2]
- The cumulative ocean spin-up exceeds 20,000 years. The model was cooled to
  about 13.6 °C global mean for acceptable Arctic sea ice.
- The historical warming was an **explicit tuning target**. The paper warns that
  the close match to observed warming "should not be taken as a sign of
  quality". The control runs have no background volcanic aerosol, so historical
  runs cool slightly when volcanoes are switched on.

## Implications for this project

- Curvilinear grid: velocities are grid-relative, and `wo`/`wmo` sit on 41
  interface levels (see `configs/datasets.yaml`).
- MPI is the only local model with `umo`/`vmo`/`wmo`. Mass transports are
  available for overturning analysis.
- Because organic matter has a fixed Redfield composition, N:P departures in
  MPI dissolved nutrients come from N-cycle processes (N2 fixation,
  denitrification, deposition) and circulation, not from variable plankton
  stoichiometry. This is an inference from the model structure.

## Source excerpts

> "The bipolar grids (GR1.5 in MPI-ESM1.2-LR and GR3.0 in MPI-ESM1.2-CR)
> feature one grid pole under Greenland and one under Antarctica." (MPI, §4)
