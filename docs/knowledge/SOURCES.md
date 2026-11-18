# Sources

Each fragment's `sources` block repeats the full citation. This file lists every
source once, records the version that was read, and identifies the parts used
by the current knowledge fragments. 

---

## OMIP-phys: physical ocean diagnostics protocol

Griffies, S. M., Danabasoglu, G., Durack, P. J., et al. (2016). OMIP contribution
to CMIP6: experimental and diagnostic protocol for the physical component of the
Ocean Model Intercomparison Project. *Geoscientific Model Development*, 9,
3231–3296. https://doi.org/10.5194/gmd-9-3231-2016 (CC BY 3.0)

Local copy: `docs/full/gmd-9-3231-2016.pdf`

- **Used:** Appendix A1–A4 and A7 (sampling, averaging, native and regridded output, missing
  fields), B (cell area, volume and mass), C1–C4 (integration, vertical coordinates,
  vector fields and zig-zag transports), D2–D5 and E (temperature, salinity, heat and
  salt content), F3 (time-dependent thickness weighting), G3–G7 and Table G1
  (reference density, heat capacity, bathymetry, basin masks and cell area), H10–H13
  and H16–H18 (temperature, global means, bottom properties and surface salinity),
  H21–H22 (buoyancy frequency and locally referenced potential density), H24–H25 and
  the Table H1 caption (mixed layer depth, squared diagnostics and variance), I1–I6
  and Table I1 (velocity, mass transport and overturning).
- **Knowledge fragments:**
  [basin-and-region-definitions](methods/basin-and-region-definitions.md),
  [cmip6-grid-labels](datasets/cmip6-grid-labels.md),
  [cmip6-mlotst-definition](concepts/cmip6-mlotst-definition.md),
  [cmip6-so-practical-salinity](concepts/cmip6-so-practical-salinity.md),
  [cmip6-static-grid-fields](concepts/cmip6-static-grid-fields.md),
  [cmip6-thetao-potential-temperature](concepts/cmip6-thetao-potential-temperature.md),
  [cmip6-tos-bulk-top-cell](concepts/cmip6-tos-bulk-top-cell.md),
  [common-comparison-mask](methods/common-comparison-mask.md),
  [derive-unavailable-diagnostics](methods/derive-unavailable-diagnostics.md),
  [missing-values-and-land-cells](pitfalls/missing-values-and-land-cells.md),
  [model-drift-in-deep-trends](pitfalls/model-drift-in-deep-trends.md),
  [monthly-means-are-time-averages](pitfalls/monthly-means-are-time-averages.md),
  [ocean-heat-and-salt-content](methods/ocean-heat-and-salt-content.md),
  [omip-native-grid-vs-regridded](methods/omip-native-grid-vs-regridded.md),
  [omip-weighted-means-and-integrals](methods/omip-weighted-means-and-integrals.md),
  [overturning-from-mass-transport](methods/overturning-from-mass-transport.md),
  [velocity-components-are-grid-relative](pitfalls/velocity-components-are-grid-relative.md)

## OMIP-BGC: biogeochemical protocol

Orr, J. C., Najjar, R. G., Aumont, O., et al. (2017). Biogeochemical protocols and
diagnostics for the CMIP6 Ocean Model Intercomparison Project (OMIP).
*Geoscientific Model Development*, 10, 2169–2199.
https://doi.org/10.5194/gmd-10-2169-2017 (CC BY 3.0)

Local copy: `docs/full/gmd-10-2169-2017.pdf`

- **Used:** Section 2 introduction (forced OMIP versus coupled historical runs and carbonate chemistry recommendations), 2.1.2 (oxygen saturation convention), 2.4 and Eq. (1)
  (volume conservation and per-volume tracer units), 2.5.3 and Eq. (16) (oxygen
  saturation and atmospheric pressure), 2.6 and Eqs. (32)–(38) (carbonate chemistry,
  alkalinity terms and total pH scale), Tables 6, 7, 10–12 and 15, including the cited
  footnotes (tracer scope, units, integrated carbon, oxygen saturation and oxygen
  minima).
- **Knowledge fragments:**
  [cmip6-bgc-units-per-volume](concepts/cmip6-bgc-units-per-volume.md),
  [cmip6-bgc-variable-scope](concepts/cmip6-bgc-variable-scope.md),
  [cmip6-carbonate-system-definitions](concepts/cmip6-carbonate-system-definitions.md),
  [cmip6-experiments-historical-ssp245](datasets/cmip6-experiments-historical-ssp245.md),
  [cmip6-integrated-carbon-content](concepts/cmip6-integrated-carbon-content.md),
  [cmip6-o2sat-and-derived-oxygen](concepts/cmip6-o2sat-and-derived-oxygen.md)

## CMOR tables: CMIP6 data request variable definitions

PCMDI. *cmip6-cmor-tables*, `Tables/CMIP6_Omon.json`, `Tables/CMIP6_Ofx.json`,
`Tables/CMIP6_coordinate.json`. data_specs_version 01.00.33, table_date
18 November 2020. https://github.com/PCMDI/cmip6-cmor-tables

Local copies: `docs/full/CMIP6_Omon.json`, `docs/full/CMIP6_Ofx.json`,
`docs/full/CMIP6_coordinate.json`

- **Used:** The Omon/Ofx headers (missing values); variable definitions, units, comments,
  dimensions, cell_methods and cell_measures for thetao, bigthetao, so, sos, sob, tos,
  tob, mlotst, mlotstmax, mlotstmin, uo, vo, wo, umo, vmo, wmo, dissic, dissoc, talk,
  ph, o2, o2sat, o2min, zo2min, no3, po4, si, dfe, chl, chlos, phycos, intdic, intdoc,
  intpoc, areacello, volcello, deptho, sftof, thkcello, masscello, basin, msftmz,
  msftyz, obvfsq, zos, spco2, fgco2, intpp, epc100 and hfds. Coordinate entries
  include basin, oline and rho, and the ocean vertical coordinates used by those
  definitions. These include unavailable variables whose definitions support the
  derivation and availability discussion in the fragments.
- **Knowledge fragments:**
  [basin-and-region-definitions](methods/basin-and-region-definitions.md),
  [cf-cell-methods](concepts/cf-cell-methods.md),
  [cmip6-bgc-units-per-volume](concepts/cmip6-bgc-units-per-volume.md),
  [cmip6-bgc-variable-scope](concepts/cmip6-bgc-variable-scope.md),
  [cmip6-carbonate-system-definitions](concepts/cmip6-carbonate-system-definitions.md),
  [cmip6-integrated-carbon-content](concepts/cmip6-integrated-carbon-content.md),
  [cmip6-mlotst-definition](concepts/cmip6-mlotst-definition.md),
  [cmip6-o2sat-and-derived-oxygen](concepts/cmip6-o2sat-and-derived-oxygen.md),
  [cmip6-so-practical-salinity](concepts/cmip6-so-practical-salinity.md),
  [cmip6-static-grid-fields](concepts/cmip6-static-grid-fields.md),
  [cmip6-thetao-potential-temperature](concepts/cmip6-thetao-potential-temperature.md),
  [cmip6-tos-bulk-top-cell](concepts/cmip6-tos-bulk-top-cell.md),
  [derive-unavailable-diagnostics](methods/derive-unavailable-diagnostics.md),
  [missing-values-and-land-cells](pitfalls/missing-values-and-land-cells.md),
  [overturning-from-mass-transport](methods/overturning-from-mass-transport.md),
  [velocity-components-are-grid-relative](pitfalls/velocity-components-are-grid-relative.md)

## CMIP6 CVs: controlled vocabularies

WCRP-CMIP. *CMIP6_CVs*, `CMIP6_grid_label.json`, `CMIP6_experiment_id.json`,
`CMIP6_source_id.json`. CV_collection_version 6.2.60.0 (modified 2025-09-04).
https://github.com/WCRP-CMIP/CMIP6_CVs

Local copies: `docs/full/CMIP6_grid_label.json`, `docs/full/CMIP6_experiment_id.json`,
`docs/full/CMIP6_source_id.json`

- **Used:** Grid labels gn, gr, gr1/gr2, gm and gnz; experiment entries historical, ssp245,
  piControl and esm-hist, including their parent-experiment metadata; source entries
  for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, especially ocean and ocnBgchem components,
  grid descriptions and top-cell thicknesses. Local variant labels and archived-grid
  attributes reported in the fragments come from local files, rather than from these
  controlled vocabularies.
- **Knowledge fragments:**
  [cmip6-experiments-historical-ssp245](datasets/cmip6-experiments-historical-ssp245.md),
  [cmip6-grid-labels](datasets/cmip6-grid-labels.md),
  [cmip6-tos-bulk-top-cell](concepts/cmip6-tos-bulk-top-cell.md),
  [historical-ssp245-member-continuity](pitfalls/historical-ssp245-member-continuity.md),
  [model-cesm2-configuration](datasets/model-cesm2-configuration.md),
  [model-gfdl-esm4-configuration](datasets/model-gfdl-esm4-configuration.md),
  [model-mpi-esm1-2-lr-configuration](datasets/model-mpi-esm1-2-lr-configuration.md),
  [velocity-components-are-grid-relative](pitfalls/velocity-components-are-grid-relative.md)

## CESM2 model description

Danabasoglu, G., Lamarque, J.-F., Bacmeister, J., et al. (2020). The Community
Earth System Model Version 2 (CESM2). *Journal of Advances in Modeling Earth
Systems*, 12, e2019MS001916. https://doi.org/10.1029/2019MS001916

Local copy: `docs/full/J Adv Model Earth Syst - 2020 - Danabasoglu - The Community Earth System Model Version 2  CESM2.pdf`

- **Used:** Section 2.2 (POP2 grid, levels, mixing, Caspian Sea treatment and MARBL structure,
  including DOM pools, stoichiometry and carbonate chemistry), Section 4
  (initialisation, spin-up and historical member start years), Section 6 (control-run
  energy imbalance and ocean drift), Section 7 (SST definition and biases), Section
  8.4 and Figure 14 (ENSO composites), and Section 9 (Eulerian-mean AMOC, AABW and
  deep ventilation).
- **Knowledge fragments:**
  [bgc-model-structure-comparison](datasets/bgc-model-structure-comparison.md),
  [cmip6-bgc-variable-scope](concepts/cmip6-bgc-variable-scope.md),
  [cmip6-carbonate-system-definitions](concepts/cmip6-carbonate-system-definitions.md),
  [cmip6-tos-bulk-top-cell](concepts/cmip6-tos-bulk-top-cell.md),
  [historical-ssp245-member-continuity](pitfalls/historical-ssp245-member-continuity.md),
  [model-cesm2-configuration](datasets/model-cesm2-configuration.md),
  [model-cesm2-reported-biases](datasets/model-cesm2-reported-biases.md),
  [model-drift-in-deep-trends](pitfalls/model-drift-in-deep-trends.md),
  [overturning-from-mass-transport](methods/overturning-from-mass-transport.md),
  [velocity-components-are-grid-relative](pitfalls/velocity-components-are-grid-relative.md)

## GFDL-ESM4 model description

Dunne, J. P., Horowitz, L. W., Adcroft, A. J., et al. (2020). The GFDL Earth
System Model Version 4.1 (GFDL-ESM 4.1): Overall coupled model description and
simulation characteristics. *Journal of Advances in Modeling Earth Systems*, 12,
e2019MS002015. https://doi.org/10.1029/2019MS002015

Local copy: `docs/full/J Adv Model Earth Syst - 2020 - Dunne - The GFDL Earth System Model Version 4 1  GFDL‐ESM 4 1   Overall Coupled Model.pdf`

The paper calls the model ESM4.1. Its CMIP6 output is published under
source_id GFDL-ESM4 (the paper's data citation is
https://doi.org/10.22033/ESGF/CMIP6.1407).

- **Used:** Sections 2.5–2.8 (MOM6 ocean grid and levels, COBALTv2 structure and DOM pools, mocsy
  constants, oxygen exchange and conservation), the development challenges on p. 10,
  Section 4.1 (spin-up, control drift and AMOC), Section 4.9 and Figures 16–17
  (surface biases and the MLD evaluation criterion), Section 4.10 and Figure 21
  (interior biases and Eulerian overturning), Section 4.11 (biogeochemistry, spin-up
  and remineralisation biases), Section 4.12 (ENSO), and the CMIP6 data citation
  identifying source_id GFDL-ESM4.
- **Knowledge fragments:**
  [bgc-model-structure-comparison](datasets/bgc-model-structure-comparison.md),
  [cmip6-bgc-variable-scope](concepts/cmip6-bgc-variable-scope.md),
  [cmip6-carbonate-system-definitions](concepts/cmip6-carbonate-system-definitions.md),
  [cmip6-mlotst-definition](concepts/cmip6-mlotst-definition.md),
  [model-drift-in-deep-trends](pitfalls/model-drift-in-deep-trends.md),
  [model-gfdl-esm4-configuration](datasets/model-gfdl-esm4-configuration.md),
  [model-gfdl-esm4-reported-biases](datasets/model-gfdl-esm4-reported-biases.md),
  [overturning-from-mass-transport](methods/overturning-from-mass-transport.md)

## MPI-ESM1.2 model description

Mauritsen, T., Bader, J., Becker, T., et al. (2019). Developments in the MPI-M
Earth System Model version 1.2 (MPI-ESM1.2) and its response to increasing CO2.
*Journal of Advances in Modeling Earth Systems*, 11, 998–1038.
https://doi.org/10.1029/2018MS001400

Local copy: `docs/full/J Adv Model Earth Syst - 2019 - Mauritsen - Developments in the MPI‐M Earth System Model version 1 2  MPI‐ESM1 2  and Its.pdf`

- **Used:** Section 2 and Table 1 (model configurations and coupling frequency), Section 4
  (MPIOM1.6, hydrostatic/Boussinesq formulation, bipolar C grid, vertical levels,
  mixing and the LR viscosity fix), Section 5, including 5.1–5.2 (HAMOCC6, DOM,
  nitrogen cycling, sinking particles and carbonate chemistry), and Sections 8.1–8.2
  (spin-up, tuning and historical experiments).
- **Knowledge fragments:**
  [bgc-model-structure-comparison](datasets/bgc-model-structure-comparison.md),
  [cmip6-bgc-variable-scope](concepts/cmip6-bgc-variable-scope.md),
  [cmip6-carbonate-system-definitions](concepts/cmip6-carbonate-system-definitions.md),
  [model-drift-in-deep-trends](pitfalls/model-drift-in-deep-trends.md),
  [model-mpi-esm1-2-lr-configuration](datasets/model-mpi-esm1-2-lr-configuration.md),
  [omip-weighted-means-and-integrals](methods/omip-weighted-means-and-integrals.md),
  [velocity-components-are-grid-relative](pitfalls/velocity-components-are-grid-relative.md)

## CF-1.7: NetCDF Climate and Forecast Metadata Conventions

Eaton, B., Gregory, J., Drach, B., et al. (2017). *NetCDF Climate and Forecast
(CF) Metadata Conventions*, Version 1.7.
https://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/cf-conventions.html

Local copy: `docs/full/cf-conventions-1.7.html` (downloaded 2026-09-24)

- **Used:** Section 2.5.1 (missing data, valid ranges and masking before scaling), Section 6.1.1
  (geographic region names), Sections 7.2–7.4 (cell measures, cell-method order,
  intensive/extensive quantities, sea-area restrictions and climatological
  statistics), and Appendix E (cell methods and units).
- **Knowledge fragments:**
  [basin-and-region-definitions](methods/basin-and-region-definitions.md),
  [cf-cell-methods](concepts/cf-cell-methods.md),
  [missing-values-and-land-cells](pitfalls/missing-values-and-land-cells.md)

## Bretherton et al. 1999: effective sample size

Bretherton, C. S., Widmann, M., Dymnikov, V. P., Wallace, J. M., and Bladé, I.
(1999). The effective number of spatial degrees of freedom of a time-varying
field. *Journal of Climate*, 12, 1990–2009.
https://doi.org/10.1175/1520-0442(1999)012<1990:TENOSD>2.0.CO;2

Local copy:

- `docs/full/bretherton-et-al-1999-effective-spatial-degrees-of-freedom.pdf`
  Downloaded from: https://atmos.washington.edu/wallace/PDFs/1999%20-%20The%20effective%20number%20of%20spatial%20degrees%20of%20freedom%20of%20a%20time-varying%20field.pdf

- **Used:** Section 5, pp. 2003–2004 (effective sample size for serial correlation), including
  Eqs. (30)–(31), the red-noise order-1 sample-mean formula, the two-series
  correlation formula, and uncertainty when estimating the covariance structure from a
  finite record.
- **Knowledge fragments:**
  [serial-correlation-and-significance](methods/serial-correlation-and-significance.md)

## Lauvset et al. 2016: GLODAPv2 mapped climatology

Lauvset, S. K., Key, R. M., Olsen, A., et al. (2016). A new global interior
ocean mapped climatology: the 1° × 1° GLODAP version 2. *Earth System Science
Data*, 8, 325–340. https://doi.org/10.5194/essd-8-325-2016 (CC BY 3.0)

Local copy:

- `docs/full/essd-8-325-2016.pdf`
  Downloaded from: https://essd.copernicus.org/articles/8/325/2016/essd-8-325-2016.pdf

- **Used:** The abstract (GLODAPv2.2016b global 1° × 1° mapped climatologies on 33 standard depth
  surfaces). The local WOA-mask overlap counts in the fragment are local checks,
  rather than results reported in this paper.
- **Knowledge fragments:**
  [common-comparison-mask](methods/common-comparison-mask.md)

## CPC: changes to the Oceanic Niño Index

NOAA Climate Prediction Center. *Changes to the Oceanic Niño Index (ONI)*.
https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_change.shtml

Local copy:

- `docs/full/cpc-ONI_change.html`
  Downloaded from: https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_change.shtml

- **Used:** The opening explanation of the ONI base-period strategy: Niño-3.4 warming since 1950,
  the reason for updating climatologies, and centred 30-year base periods applied to
  successive five-year spans. The fragment's local-product comparisons and tool
  instructions are project-specific interpretations.
- **Knowledge fragments:**
  [anomaly-baselines-and-climate-indices](methods/anomaly-baselines-and-climate-indices.md)

## WOA23 product documentation

Reagan, J. R., Garcia, H. E., Boyer, T. P., et al. (2024). *World Ocean Atlas
2023: Product Documentation*. NOAA National Centers for Environmental
Information, Silver Spring, MD.
https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DOCUMENTATION/WOA23_Product_Documentation.pdf

Local copy:

- `docs/full/WOA23_Product_Documentation.pdf`
  Downloaded from: https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DOCUMENTATION/WOA23_Product_Documentation.pdf

- **Used:** The file-naming convention and averaging-period codes on p. 13: 13–16 denote Northern
  Hemisphere winter January–March, spring April–June, summer July–September and autumn
  October–December, respectively. These support the distinction between WOA calendar
  quarters and meteorological DJF/MAM/JJA/SON seasons.
- **Knowledge fragments:**
  [anomaly-baselines-and-climate-indices](methods/anomaly-baselines-and-climate-indices.md)

## NCAR Climate Data Guide: PDO and EOF analysis

National Center for Atmospheric Research Staff (Eds.). *The Climate Data Guide:
Pacific Decadal Oscillation (PDO): Definition and Indices*.
https://climatedataguide.ucar.edu/climate-data/pacific-decadal-oscillation-pdo-definition-and-indices

National Center for Atmospheric Research Staff (Eds.). *The Climate Data Guide:
Empirical Orthogonal Function (EOF) Analysis and Rotated EOF Analysis*.
https://climatedataguide.ucar.edu/climate-data-tools-and-analysis/empirical-orthogonal-function-eof-analysis-and-rotated-eof-analysis


Local copies:

- **Used:** The PDO guide's definition and expert guidance (leading North Pacific SST-anomaly EOF,
  removal of the seasonal cycle and global-mean SST anomaly, positive-phase sign
  convention, and regression onto a normalised PC); the EOF guide's spatial weighting,
  including sqrt(cos(latitude)), and domain dependence of EOF patterns. The fragment's
  implementation-specific instructions refer to the project's own tools.
- **Knowledge fragments:**
  [eof-modes-and-pdo-index](methods/eof-modes-and-pdo-index.md)
