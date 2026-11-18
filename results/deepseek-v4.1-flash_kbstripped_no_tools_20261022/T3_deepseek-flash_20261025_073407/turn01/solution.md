## 1. Direct answer

All three CMIP6 historical 1971–2000 SST climatologies reproduce the large-scale WOA23 pattern almost perfectly (spatial r ≥ 0.99), but they differ in their globally averaged offset. On the shared 1° comparison mask, CESM2 is warm (+0.30 °C), GFDL‑ESM4 is nearly unbiased (−0.04 °C), and MPI‑ESM1‑2‑LR is cold (−0.55 °C). Ranked from smallest to largest absolute bias: **GFDL‑ESM4 (0.04 °C) < CESM2 (0.30 °C) < MPI‑ESM1‑2‑LR (0.55 °C)**. RMSE follows the same ordering (0.894, 1.156, 1.415 °C), so GFDL‑ESM4 is the best of the three both in the mean and in the pointwise sense, and MPI‑ESM1‑2‑LR the worst.

## 2. Methodology exactly as executed

**WOA23 reference.** Opened `woa23_decav71A0_t00_01.nc` with `decode_times=False`, took the objectively analysed field `t_an` (the only variable ending in `_an`), the single annual time slice (`isel(time=0)`) and the 0 m level (`isel(depth=0)`), verified `depth = 0.0 m`. That field is the 1971–2000 annual climate normal on the 1° grid (180 × 360, lon −179.5…179.5, lat −89.5…89.5). A `|t| < 1e10` guard removed undecoded fills. The land–sea mask was read from the text file `landsea_01.msk` (64800 rows, columns lat, lon, value), binned to the 1° grid by nearest-index rounding, and ocean defined as `landsea > 1`. The requested **WOA23 global mean is the plain unweighted arithmetic mean of all non‑missing 0 m cells** (41088 valid cells; the `landsea > 1` count also equals 41088).

**CMIP6 model data.** `tos` and the Ofx `areacello` were opened for CESM2, GFDL‑ESM1‑2‑LR — i.e. CESM2, GFDL‑ESM4 and MPI‑ESM1‑2‑LR — all `r1i1p1f1` (variant parsed from the filenames, not assumed). Monthly data were restricted to years 1971–2000 (360 months for every model, confirmed), masked with `|tos| < 1e19`, and reduced to a cell-wise time mean (`mean(dim='time', skipna=True)`). Each model's global mean is the area-weighted sum Σ(tos·areacello)/Σ(areacello) over cells where both `tos` and `areacello` are finite and positive — each model on its own native grid with its own published cell areas, no rescaling between models. Units °C throughout.

**Regridding.** Conservative first-order donor-cell remapping onto the WOA23 1° grid: every native cell centre is binned into the 1° target cell containing it and donates its full `value × areacello`; the target value is the area-weighted mean of its donors. Nothing else was regridded — the native-grid means were computed before, and independently of, this step.

**Comparison mask.** On the 1° grid: `landsea > 1` AND valid `woa_sst` AND valid regridded fields of all three models → **21 961 cells, 186.35 × 10⁶ km²**. All fields were `.where(mask)`-ed before any metric.

**Metrics.** Area weighting proportional to cos(latitude) for the shared-mask statistics: mean model SST, mean WOA23 over the *same* cells, ΔSST = mean(model) − mean(WOA23), RMSE = sqrt(Σw·d²/Σw), and weighted Pearson r about the weighted means.

## 3. Key numerical results

| Quantity | CESM2 | GFDL‑ESM4 | MPI‑ESM1‑2‑LR | WOA23 |
|---|---|---|---|---|
| Native-grid area-weighted global mean (°C) | **18.52** | **18.13** | **17.58** | — |
| Native ocean area Σ`areacello` (10¹² m²) | 3.6051 | 3.6111 | 3.6598 | — |
| Whole-1°-grid unweighted mean (°C) | — | — | — | **13.85** |
| Shared-mask mean model SST (°C) | 17.03 | 16.70 | 16.19 | — |
| Shared-mask mean WOA23 SST (°C) | 16.73 | 16.73 | 16.73 | 16.73 |
| ΔSST model − WOA23 (°C) | **+0.30** | **−0.04** | **−0.55** | — |
| Area-weighted RMSE (°C) | **1.156** | **0.894** | **1.415** | — |
| Area-weighted spatial r | **0.994** | **0.996** | **0.992** | — |
| Bias-field min / max (°C) | −4.88 / +8.61 | −4.41 / +7.39 | −7.45 / +7.68 |
| Bias-field std over mask (°C) | 1.14 | 0.90 | 1.34 |

Ranking by |ΔSST| (unrounded): GFDL‑ESM4 (0.036) < CESM2 (0.303) < MPI‑ESM1‑2‑LR (0.545).

Two numbers deserve emphasis. First, the WOA23 bar (13.85 °C) is a *different quantity* from the shared-mask WOA23 mean (16.73 °C): the former is the unweighted mean over the full 1° grid including all valid ocean cells, the latter is a cos(lat)-weighted mean over the 21 961 cells where all four products are simultaneously valid; the four bars in Figure 2 are therefore not strictly comparable. Second, the native-grid model means (18.52, 18.13, 17.58 °C) are much warmer than the shared-mask means (17.03, 16.70, 16.19 °C), because the native means cover each model's own full ocean (including high-latitude ocean weighted by genuine cell area) whereas the shared mask is restricted to cells where even the sparsest product (MPI, 26 313 valid regridded cells) is defined. This is exactly the native-grid-vs-regridded discrepancy the plan warns about, and the two must not be expected to agree.

## 4. Interpretation

The models capture the spatial pattern of SST essentially perfectly (r = 0.992–0.996), so the ranking is governed by the mean offset and by regionally concentrated errors. GFDL‑ESM4 has both the smallest global offset (−0.04 °C) and the smallest dispersion of bias (std 0.90 °C, RMSE 0.894 °C), i.e. a nearly unbiased field with modest regional structure. CESM2 is warm on average by 0.30 °C with larger excursions (std 1.14, RMSE 1.156). MPI‑ESM1‑2‑LR is cold by 0.55 °C and has by far the widest bias range (−7.45 to +7.68 °C, std 1.34, RMSE 1.415); its poor RMSE despite a modest mean bias is a signature of compensating regional errors, not of a uniform cold offset. The three-model spread of ~0.85 °C in the shared-mask bias is comparable to the inter-model spread one expects from a single historical member, and the common-mask restriction (21 961 of ~41 000 WOA ocean cells, 186 × 10⁶ km²) reflects disagreement about where the ocean is at the sea-ice margins and around coastlines as much as it reflects model skill.

## 5. Figures produced

**Figure 1** is a 2 × 2 pcolormesh panel on the WOA23 1° grid, all four panels using the same lat/lon axes and the same shared comparison mask. The three panels CESM2, GFDL‑ESM4 and MPI‑ESM1‑2‑LR show (model − WOA23) computed from the conservatively regridded climatologies, drawn with the diverging `RdBu_r` colormap on a symmetric scale ±v_max where v_max is the pooled 98th percentile of |bias| over the three models (annotated in the suptitle), so that colours are directly comparable between models and the sign is read off the zero-centred scale; each panel title also carries that model's ΔSST, RMSE and r. The fourth panel shows the WOA23 0 m 1971–2000 SST on the same cells with the sequential `viridis` colormap, limited to its 2nd–98th percentiles, as the reference field. Every panel is labelled with longitude/latitude and every colour bar carries a °C unit label ("SST bias (degC)" and "SST (degC)") and triangular `extend` ends. The suptitle records the v_max clipping and the mask cell count.

**Figure 2** is a bar chart of four numbers, one per product, each bar annotated with its value to 2 decimal places and the y axis labelled in °C: the three native-grid, own-`areacello`-weighted model means and the WOA23 whole-1°-grid unweighted mean. A caption note states explicitly that the model bars are area-weighted over each native ocean grid while the WOA23 bar is the unweighted whole-grid mean, so the four bars represent slightly different quantities.

Describing the spatial structure seen in Figure 1: the bias is not a uniform offset but has a clear latitudinal/regional organisation. The largest amplitudes are positive and occur in the high-latitude Southern Ocean and Arctic/sea-ice-margin cells (up to about +7.7 to +8.6 °C in CESM2 and MPI), where the models' bulk top-cell `tos` reports near-freezing under-ice water while WOA23's 0 m analysis differs in its sea-ice treatment, and in eastern-boundary upwelling regions. Strong negative biases (down to −4.4 to −7.4 °C) appear in the same high-latitude and western-boundary/current-extension zones. The tropical and subtropical open ocean, which dominates the area and hence the mean, is reproduced within a few tenths of a degree by all three models — that is why r is so high and the global ΔSST so small compared with the local extremes. MPI's bias field is the noisiest and most sign-alternating at high latitudes, consistent with its larger std and RMSE.

## 6. Caveats from the executed analysis

- `tos` is bulk top-cell temperature (CESM2 0–10 m, GFDL‑ESM4 0–2 m, MPI‑ESM1‑2‑LR 0–12 m), not skin SST, so the comparison with a 0 m WOA23 analysis is not strictly like-for-like; sea-ice-covered cells are not missing in `tos` and carry near-freezing values.
- GFDL‑ESM4's "native" `tos` is on the `gr` 1° archive grid with the regridded land–sea mask, so its native-grid number is an archive-grid number.
- A single `r1i1p1f1` historical member cannot separate forced response from internal variability, and a model year is not expected to match the observed year; the 1971–2000 climatology bias carries that caveat.
- The native-grid means and the shared-mask/regridded means are different quantities by construction; only the latter enter the bias, RMSE and correlation statistics and Figure 1.
