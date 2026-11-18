## What was actually produced — and the one thing that broke

The run "succeeded" in the sense that every script statement executed and both figures were written, but the **WOA23 land–sea mask was never read**. The parser opened `/data1/.../masks/landsea_01.msk` and returned an **empty table (0 rows)** — i.e. every line failed the `len(parts) >= 3` or the float-conversion test. Consequently:

- `mask_arr` has **0 valid cells** → `woa_ocean = (mask_arr > 1) & finite = all False`
- `woa_global_mean = NaN`, `n_common = 0`
- All comparison statistics (`model_mean`, `woa_mean`, `bias`, `RMSE`, `corr`) are **NaN**, and `ranking` fell back to the input order rather than being ordered by |bias|.

So of the requested numbers, **only the three native‑grid model means are trustworthy**. Everything that depends on the WOA mask (WOA global mean, comparison‑cell means, bias, RMSE, correlation, ranking, and the bias maps in Figure 1) is undefined in this run.

## Numbers that are valid: native‑grid global mean SST, 1971–2000 (°C)

Area‑weighted with each model's own published `areacello` over its own native ocean grid:

| Model | Global mean SST (°C) |
|---|---|
| CESM2 | **18.52** |
| GFDL‑ESM4 | **18.13** |
| MPI‑ESM1‑2‑LR | **17.58** |
| WOA23 | **undefined (NaN)** |

The three CMIP6 values are physically sensible (17.6–18.5 °C) and their spread of ~0.9 °C is typical of CMIP6 historical SST climatologies. The WOA23 bar is missing entirely from Figure 2 because its value is NaN (the bar has zero height and no annotation).

*Diagnostic only, not the requested quantity:* the unweighted mean of the finite values of the WOA23 0 m annual field is **13.85 °C** (41 088 finite cells of 64 800). This is **not** a substitute for the requested mask‑based global mean — it is unweighted and not restricted to the `landsea_01` ocean cells — so it should not be reported as the WOA23 reference.

## Numbers that could not be produced

On the WOA23 1‑degree comparison cells, none of the following exist in this run:

- mean SST of each model and of WOA23 on the common cells → NaN
- mean bias ΔSST (model − WOA23) → NaN
- RMSE and spatial Pearson correlation → NaN
- ranking from smallest to largest |ΔSST| → not meaningful (returned as CESM2 < GFDL‑ESM4 < MPI‑ESM1‑2‑LR only because the sort key was NaN for all three)

Note also that the regridding step itself mostly worked: CESM2 produced 45 710, GFDL‑ESM4 47 623 and MPI‑ESM1‑2‑LR 45 262 finite cells on the 1‑degree grid. The intersection with WOA23 was empty **only** because `woa_ocean` was empty — not because of the regridding or the distance filter.

## Figures

**Figure 1** (2×2 global maps): the three bias panels (CESM2, GFDL‑ESM4, MPI‑ESM1‑2‑LR) and the WOA23 reference panel are **blank**, because `bias_fields` and `woa_field` were masked with the all‑False `common` array. Colour bars and axis labels render, but there is no field drawn and therefore no bias structure to interpret. The shared diverging scale is meaningless here (its 99th percentile was taken over an empty array).

**Figure 2** (bar chart): shows the three model bars annotated 18.52, 18.13 and 17.58; the WOA23 bar is absent (NaN).

## Why the bias structure cannot be described

The requested spatial description — e.g. the characteristic warm cold‑tongue / equatorial Pacific bias, the warm high‑latitude Southern Ocean and North Atlantic biases typical of CMIP6 `tos` — **cannot be read off this run**, since no bias field was computed. Any such description would have to come from a corrected run, not from these outputs.

## Fix required before the comparison can be reported

The `landsea_01.msk` file is evidently **not** in the `lat lon value` whitespace format the parser assumed (WOA `.msk` files are fixed‑format/row‑oriented, and the naive `split()` + `float()` route yields zero records). Options:

1. Parse the mask with its actual fixed‑width/row format (read the header record, then slice fixed columns per row), or
2. Use the WOA23 1‑degree land–sea mask distributed as a netCDF/`sftof`‑style field instead of the `.msk` file, or
3. As a stopgap, take the ocean domain as `np.isfinite(woa_arr)` from the `t_an` field itself (41 088 cells) — this reproduces a plausible ocean domain but is **not** the `landsea_01 > 1` criterion the task specifies, so label it clearly if used.

Once `mask_arr` is populated, the same script will produce the WOA23 global mean, the common‑cell means, ΔSST, RMSE, correlation, the ranking, and a populated Figure 1 and Figure 2 without further changes.
