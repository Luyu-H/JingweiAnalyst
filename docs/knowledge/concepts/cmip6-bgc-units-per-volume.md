---
id: cmip6-bgc-units-per-volume
title: "CMIP6 biogeochemical tracers are per unit volume (mol m-3, kg m-3)"
type: concept
summary: "All CMIP6/OMIP BGC tracer concentrations are stored per unit volume (mol m-3; chlorophyll kg m-3) because the models conserve volume; observational products are per unit mass, and the scale factors to the usual units differ by variable."
keywords: [units, 单位, mol m-3, umol/kg, μmol kg-1, per volume, per mass, 体积浓度, 质量浓度, unit conversion, 单位换算, chl, kg m-3, mg m-3, nitrate, phosphate, oxygen, iron, nM]
datasets: [cmip6_historical, cmip6_ssp245, woa23, glodap]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [o2, o2sat, no3, po4, si, dfe, dissic, dissoc, talk, chl, phycos]
tools: [transform_unit, convert_concentration_to_content, compute_potential_density_anomaly, identify_oxygen_minimum_zone, compute_ratio]
severity: must
related: [cmip6-carbonate-system-definitions, cmip6-bgc-variable-scope, cmip6-integrated-carbon-content]
sources:
  - citation: "Orr, J. C., et al. (2017). Biogeochemical protocols and diagnostics for the CMIP6 Ocean Model Intercomparison Project (OMIP). Geosci. Model Dev., 10, 2169–2199."
    locator: "Section 2.4, Eq. 1 (p. 2177); Tables 6, 10, 12 (pp. 2184, 2188, 2190)"
    url: https://doi.org/10.5194/gmd-10-2169-2017
  - citation: "PCMDI cmip6-cmor-tables, Tables/CMIP6_Omon.json, data_specs_version 01.00.33 (18 Nov 2020)."
    locator: "variable_entry o2, no3, po4, si, dfe, chl, phycos (units, comment)"
    url: https://github.com/PCMDI/cmip6-cmor-tables
---

## Key points

- In the OMIP tracer equation "conservation of volume is assumed … and
  standard units of mol m−3 are used for all tracers". [OMIP-BGC §2.4]
- The CMOR comments define "mole concentration" as "number of moles per unit
  volume, also called 'molarity'". [CMOR Omon, o2/no3/po4/si]
- Units as archived [CMOR Omon; OMIP-BGC Tables 6, 10, 12]:

  | variable | CMOR units | common reporting unit | factor (volume basis) |
  |---|---|---|---|
  | o2, o2sat, no3, po4, si, dissic, dissoc, talk, phycos | mol m-3 | mmol m-3 = µmol L-1 | × 1e3 |
  | dfe | mol m-3 | nmol L-1 (nM) | × 1e6 |
  | chl | kg m-3 | mg m-3 (≈ µg L-1) | × 1e6 |
  | ph | 1 | – | – |

  The factors in the last column are unit arithmetic. They change the prefix
  but keep the volume basis.

## Implications for this project

- Numbers in mmol m-3 are **still per volume**. WOA23 and GLODAP report
  µmol kg-1 (per mass). Changing the prefix does not make the two comparable.
  Only a density conversion does: µmol kg-1 = mmol m-3 × 1000 / ρ, where ρ is in
  kg m-3. `configs/datasets.yaml` (`cross_dataset_notes.units`) records this
  and notes that no loader converts automatically.
- The density error only cancels in quantities that are invariant under scaling:
  ratios of two tracers converted with the same ρ, or depths of extrema and
  gradients. Absolute concentrations, inventories against observations and
  thresholds on a mass basis all need the conversion. The scale-invariant
  quantities need no conversion at all and are best taken from the native
  field: a constant ρ cancels exactly, while the in-situ density, which grows
  with depth, rescales each level by a different factor and so itself shifts
  a depth located from the shape of a profile (a nutricline, the depth of a
  maximum gradient) by a little.
- Inventories (mol) are volume integrals of the per-volume field: Σ c · volcello.
  No density is needed for these on the model side.
- Two concentrations on different bases cannot be compared, differenced or
  evaluated against each other until one has been converted. The conversion
  needs a density, and the one that mass-per-volume arithmetic needs is the
  in-situ density, `compute_potential_density_anomaly` with `ref_pressure=None`
  (a sigma at a fixed reference pressure is not). `transform_unit` defaults to a
  constant density, the mean surface value, which errs by about a per cent in
  the deep ocean and in fresh polar water, enough to rival the differences a
  model-observation comparison sets out to measure. Pass the in-situ density
  field for anything below the surface that needs converting (not for the
  scale-invariant quantities above), and state which density was used (the
  result's `attrs['density_assumption']`).
- For inventories prefer the volume basis: a model's native mol m-3 field
  integrates exactly with `convert_concentration_to_content`, with no density
  in the arithmetic. Handing it a µmol kg-1 field without an explicit density
  applies the tool's constant default instead, which undoes any in-situ density
  used to make that field and biases the inventory by about one per cent.

## Common mistakes

- Multiplying mol m-3 by 1000 and labelling the result "µmol/kg".
- Applying a µmol kg-1 threshold (for example a hypoxia threshold defined in
  µmol kg-1) to a mol m-3 field without conversion: every cell then falls below
  it, and the whole ocean is marked as depleted.
- Reading chlorophyll in kg m-3 as mg m-3. The two differ by a factor of 1e6.

## Source excerpts

> "Conservation of volume is assumed in Eq. (1) and standard units of mol m−3
> are used for all tracers." (OMIP-BGC, §2.4)
