The requested sea-surface-salinity assessment could not be completed because the execution failed before any numerical results or figures were produced.

The failure occurred while applying the WOA23 land–sea mask:

`landsea_da.sel(lat=woa_lat, lon=woa_lon, method="nearest")`

with the error:

`"not all values found in index 'lat'"`

This is a coordinate-alignment problem, not a scientific result. It typically occurs when the target WOA latitude/longitude coordinates are outside the mask coordinate range or use a different longitude convention, so the nearest-neighbour `.sel()` lookup cannot find an index for every requested point. Consequently, no values were generated for:

- global mean 1971–2000 SSS or model biases versus WOA23;
- RMSE or spatial correlation against WOA23;
- historical or ssp245 global-mean SSS trends, confidence intervals, R², or slope ratios;
- basin-scale ΔSSS for the Atlantic, Indian, and Pacific;
- Figures 1–3.

To obtain the requested assessment, the WOA mask must first be robustly mapped onto the WOA23 salinity grid, e.g. by explicitly aligning longitude conventions and using a nearest-neighbour lookup such as `cKDTree`, or by clipping/aligning coordinates before using `reindex(..., method="nearest")` or `interp`. After fixing this mask-alignment step, the same workflow can be re-run to produce the summary table and the three requested figures.
