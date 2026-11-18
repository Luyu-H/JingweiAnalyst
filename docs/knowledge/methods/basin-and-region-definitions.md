---
id: basin-and-region-definitions
title: "Ocean basin partitions in CMIP6/OMIP and CF standard region names"
type: method
summary: "OMIP's recommended 10-region tracer mask separates the Southern Ocean, marginal seas and the Arctic from the Atlantic/Pacific/Indian basins, while overturning uses Atlantic–Arctic / Indo-Pacific / global; basin statistics depend on which partition is used and must say so."
keywords: [basin mask, 海盆掩膜, region, 区域, southern ocean, 南大洋, Atlantic, Pacific, Indian, Arctic, marginal seas, 边缘海, Mediterranean, flag_values, basin partition, 海盆划分]
datasets: [cmip6_historical, cmip6_ssp245, woa23]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [basin]
tools: [select_basin, load_basin_mask, build_comparison_mask]
severity: should
related: [overturning-from-mass-transport, omip-native-grid-vs-regridded]
sources:
  - citation: "Griffies, S. M., et al. (2016). OMIP contribution to CMIP6: experimental and diagnostic protocol for the physical component of the Ocean Model Intercomparison Project. Geosci. Model Dev., 9, 3231–3296."
    locator: "Appendix G6 (p. 3252); Table G1 (p. 3252); Appendix I6 (p. 3266)"
    url: https://doi.org/10.5194/gmd-9-3231-2016
  - citation: "PCMDI cmip6-cmor-tables, Tables/CMIP6_Ofx.json and CMIP6_coordinate.json, data_specs_version 01.00.33 (18 Nov 2020)."
    locator: "Ofx variable_entry basin; coordinate axis_entry basin, oline"
    url: https://github.com/PCMDI/cmip6-cmor-tables
  - citation: "Eaton, B., et al. (2017). NetCDF Climate and Forecast (CF) Metadata Conventions, Version 1.7."
    locator: "Section 6.1.1 Geographic Regions"
    url: https://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/cf-conventions.html
---

## Key points

- OMIP recommends a tracer region mask with CF standard region names and
  `flag_values = 1..10` [OMIP-phys G6]:

  | flag | region | flag | region |
  |---|---|---|---|
  | 1 | southern_ocean | 6 | mediterranean_sea |
  | 2 | atlantic_ocean | 7 | black_sea |
  | 3 | pacific_ocean | 8 | hudson_bay |
  | 4 | arctic_ocean | 9 | baltic_sea |
  | 5 | indian_ocean | 10 | red_sea |

  Some enclosed seas "may be missing from the model". [OMIP-phys G6]
- The CMOR Ofx `basin` variable (standard_name `region`) carries this mask on
  the model grid. It is requested on both the native and the spherical grid.
  [CMOR Ofx; OMIP-phys Table G1]
- For **overturning and meridional transports** the partition is different:
  `atlantic_arctic_ocean`, `indian_pacific_ocean`, `global_ocean`.
  [CMOR coordinate `basin`; OMIP-phys I6]
- The CMOR `oline` axis lists the standard passages (e.g. drake_passage,
  indonesian_throughflow, bering_strait, fram_strait, denmark_strait,
  florida_bahamas_strait, pacific_equatorial_undercurrent). [CMOR coordinate]
- CF: regions with complex boundaries are identified by a label variable with
  standard_name `region`, using the standardised region names. [CF 6.1.1]

## Implications for this project

- The Ofx `basin` field is not in the local archive. Basin statistics here use
  the WOA23 basin mask carried onto each grid (`select_basin`,
  `load_basin_mask`). That mask follows a WOA labelling, not the OMIP 10-region
  one. Report which mask was used, whether the Southern Ocean is a separate
  basin or split among the three oceans, and whether marginal seas are
  folded into their parent ocean. Each of these choices changes basin means.
- In the WOA23 mask each marginal sea carries its own code, so an ocean taken
  as its main code alone silently excludes a large part of what the name
  usually means. `select_basin` folds the marginal seas into the ocean they
  open onto by default (`include_marginal_seas=True`); either choice is
  defensible as long as it is stated.
- A marginal sea's label stops at its sill, and below that depth the water
  belongs to the open ocean, so the mask labels each standard level
  separately. The surface labels applied through the whole water column
  (`depth_dependent=False`) are what a named ocean usually means; each level's
  own labels (`depth_dependent=True`) are what a water-mass analysis wants.
- Model domains differ: the CESM2 paper notes that the Caspian Sea is treated
  as a lake in the land model and "is no longer included in the ocean model"
  (see `model-cesm2-configuration`). Marginal seas present in an observational
  mask may therefore have no model counterpart.

## Common mistakes

- Reporting "Indian Ocean" overturning. OMIP says the Indian and Pacific
  basins cannot be separated for overturning.
- Comparing basin means computed with different basin masks for model and
  observations.

## Source excerpts

> "This analysis generally involves the use of a mask that partitions the
> model grid into ocean basins (some enclosed seas may be missing from the
> model)." (OMIP-phys, G6)
