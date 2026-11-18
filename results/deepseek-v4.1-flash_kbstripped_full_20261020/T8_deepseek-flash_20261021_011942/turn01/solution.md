## 1. What was done (as executed)

**Data.** Monthly CMIP6 `historical` output from CESM2, GFDL-ESM4 and MPI-ESM1-2-LR: `dissic`, `talk`, `ph`, `thetao`, `so` for calendar year 2002 (averaged over the 12 monthly values into one annual-mean state), plus `dissic`, `thetao`, `so` for 1971-01 to 2000-12 (loaded with `depth_range=[0,50]` and averaged over time into one 1971–2000 mean state). Grid weights were taken from each model's own published `areacello` and `volcello` (`area_source`/`volume_source` = "… published by cmip6_historical" for all three models; the code never used `thkcello`/`masscello`/`rhozero`).

**Reference.** `load_climatology_data(dataset="glodap")` for `TCO2`, `TAlk` and `pHtsinsitutp` — GLODAPv2.2016b mapped climatology, 33 standard levels, no time dimension. Model `ph` (total scale, in situ) was compared only with `pHtsinsitutp`; pH was not rescaled.

**Units.** Model `dissic`/`talk` are mol m⁻³. They were converted with `transform_unit` to mmol m⁻³ and then divided by a cell-by-cell absolute in-situ density, ρ = `compute_potential_density_anomaly(thetao, so, ref_pressure=None)` + 1000, to give µmol kg⁻¹. Density assumption recorded for all three models: *"cell-by-cell in-situ density (sigma computed with TEOS-10, ref_pressure=None, +1000); umol kg⁻¹ = mmol m⁻³ × 1000 / ρ"*. `so` was used as loaded, without rescaling (its "0.001" units string was not applied as a multiplier).

**Weighting.** Depth-spanning means (all 0–50 m means) used `compute_basic_stat(op="volume_weighted_mean")` restricted to the levels each column has in 0–50 m; horizontal metrics used area weighting on the GLODAP 1° grid. No unweighted `mean()` over (y, x) was taken. Native-grid and regridded numbers are reported separately below.

## 2. Direct answer

**No single model wins everywhere, but CESM2 is the best overall reproduction of the observed carbon system.** It has by far the smallest pH error (−0.001), the smallest alkalinity error of the two low-bias models in relative terms and the best spatial pattern of surface DIC (r = 0.828, RMSE 60 µmol kg⁻¹). GFDL-ESM4 is closest on the surface DIC and alkalinity *inventory* (+23 and +31 µmol kg⁻¹ vs GLODAP) but has the worst spatial pattern (r = 0.709) and a larger pH error (+0.008). MPI-ESM1-2-LR is worst on every carbon metric (DIC −61, TAlk −87 µmol kg⁻¹, pH −0.022, RMSE 83 µmol kg⁻¹, r = 0.736).

**Sign consistency: yes.** Each model's three errors have the same sign — CESM2 DIC/TAlk/pH all negative, GFDL-ESM4 all positive, MPI-ESM1-2-LR all negative — the physically expected behaviour for a carbon system coupled through carbonate chemistry (low DIC + low TAlk ⇒ low pH).

## 3. Key numbers

### 3.1 Global 0–50 m means (each source on its own native grid, volume-weighted; µmol kg⁻¹ / pH)

| Source | DIC | Δ vs GLODAP | TAlk | Δ vs GLODAP | pH | Δ vs GLODAP |
|---|---|---|---|---|---|---|
| CESM2 | 1997 | −39 | 2264 | −40 | 8.080 | −0.001 |
| GFDL-ESM4 | 2059 | +23 | 2336 | +31 | 8.088 | +0.008 |
| MPI-ESM1-2-LR | 1974 | −61 | 2218 | −87 | 8.059 | −0.022 |
| GLODAP | 2036 | — | 2304 | — | 8.080 | — |

These averages each cover that source's *own* valid 0–50 m cells (a different set per source), so they are not directly the same as the shared-mask biases below.

### 3.2 0–50 m DIC map, conservative regrid to the GLODAP 1° grid, single shared mask (41 227 cells, `ocean_only=False`)

| Model | Mean bias (µmol kg⁻¹) | RMSE (µmol kg⁻¹) | Spatial correlation |
|---|---|---|---|
| CESM2 | −37 | 60 | 0.828 |
| GFDL-ESM4 | +25 | 65 | 0.709 |
| MPI-ESM1-2-LR | −59 | 83 | 0.736 |

On this mask the model means are 1997.5, 2059.3 and 1974.7 vs a GLODAP mean of 2034.0 µmol kg⁻¹; the biases (−36.5, +25.3, −59.4) are close to, but not identical with, the native-grid differences of §3.1.

### 3.3 Vertical profile of DIC bias (model − GLODAP, area-weighted on the 1° grid, both defined)

| Depth | CESM2 | GFDL-ESM4 | MPI-ESM1-2-LR | cells (CESM2/GFDL/MPI) |
|---|---|---|---|---|
| 0 m | −29 | +29 | −58 | 40 921 / 41 301 / 40 828 |
| 100 m | −46 | +19 | −59 | 38 318 / 38 534 / 38 355 |
| 500 m | −12 | −11 | −78 | 36 370 / 36 561 / 36 294 |
| 1000 m | +2 | −1 | −67 | 35 727 / 35 975 / 35 688 |
| 2000 m | −2 | −11 | −61 | 33 041 / 33 795 / 33 109 |

CESM2's and GFDL-ESM4's biases converge to near zero below ~500 m (a near-cancellation of the surface offset by the deeper ocean), whereas MPI-ESM1-2-LR carries a persistent ≈ −60 to −80 µmol kg⁻¹ low bias over the whole water column — a whole-ocean carbon inventory deficit, not a surface-only problem.

### 3.4 Anthropogenic change and the bias/change ratio

| Model | 2002 minus 1971–2000 (0–50 m DIC, native grid) | Ratio to the step-5 DIC bias |
|---|---|---|
| CESM2 | +9 µmol kg⁻¹ | −0.24 |
| GFDL-ESM4 | +12 µmol kg⁻¹ | +0.53 |
| MPI-ESM1-2-LR | +11 µmol kg⁻¹ | −0.18 |

**Implication for reading the bias.** In every model the 1971–2000 → 2002 change (+9 to +12 µmol kg⁻¹) is only 18–53 % of the magnitude of that model's DIC bias against GLODAP. The bias is therefore *larger* than the entire accumulated anthropogenic DIC signal of the last three decades of the 20th century: it cannot be interpreted as a timing/lag error (i.e. a model simply "running a few years ahead or behind" the observations), because that would require a temporal offset of roughly two to five decades. The DIC offsets are genuine systematic model field errors (inventory, saturation state and/or biological/ventilation parameterisations), not drift artefacts. The sign of the ratio also differs between models, so no common phase offset can reconcile them with GLODAP.

## 4. Figures

**Figure 1** — three panels in one row (one per model) of the 0–50 m DIC bias (model − GLODAP) on the 41 227-cell shared GLODAP-grid mask, drawn with a single diverging `RdBu_r` colour scale centred on zero and shared across all three panels, symmetric and clipped at the pooled 98th percentile of |bias| so the large high-latitude outliers do not flatten the open-ocean pattern; one common colourbar labelled in µmol kg⁻¹. Each title carries that model's map bias, RMSE and correlation. The intended reading is the *spatial structure* of the surface carbon error under a common scale and a common cell set.

**Figure 2** — a single panel showing the vertical profile of the global mean DIC bias: depth on the vertical axis increasing downward, bias in µmol kg⁻¹ on the horizontal axis, one line per model with markers at 0, 100, 500, 1000 and 2000 m and a vertical zero reference line. It is designed to show how each model's bias is distributed through the water column rather than at the surface only.

## 5. Insights and caveats

- The three models bracket the observation on the surface carbon inventory: GFDL-ESM4 is high (+23 DIC, +31 TAlk), CESM2 mildly low (−39, −40), MPI-ESM1-2-LR strongly low (−61, −87). Because DIC and TAlk errors travel together, the DIC/TAlk *difference* (the carbonate-system buffer) is better preserved than the absolute inventory.
- pH is a much more forgiving diagnostic: despite ±40–90 µmol kg⁻¹ DIC/TAlk errors, CESM2's surface pH is within 0.001 of GLODAP and even MPI-ESM1-2-LR's is only 0.022 low — consistent with pH being set by the DIC/TAlk ratio rather than by either absolute quantity. Ranking models on pH alone would therefore be misleading.
- Interpretation of the pH result depends on the model's carbonate constants. The code compared the models' own reported `ph` on the total scale in situ, but the CESM2 metadata as loaded does not document which carbonate dissociation constants its `ph` diagnostic uses — this is an **unverified assumption** and should be flagged when attributing CESM2's near-zero pH difference to model skill; the same caveat applies to GFDL-ESM4 and MPI-ESM1-2-LR.
- Provenance caveat: the loaded attributes show MPI-ESM1-2-LR `ph` comes from variant `r10i1p1f1`, whereas its `dissic`/`talk`/`thetao`/`so` come from `r1i1p1f1`; the MPI pH comparison is therefore not a strictly single-member, single-trajectory comparison of the carbon system, and its pH error should be read with that in mind. CESM2 and GFDL-ESM4 are internally `r1i1p1f1` throughout.
- The 0–50 m global means and the map/profile metrics rest on different cell sets (native valid cells vs the shared GLODAP-grid mask, `ocean_only=False` because the reference is GLODAP, not WOA23); the model-minus-GLODAP difference of §3.1 and the shared-mask bias of §3.2 differ by 1–2 µmol kg⁻¹ and must not be quoted interchangeably.
- GLODAPv2 has no time dimension, so its single TCO2/TAlk/pH climatology is the fixed reference against which each model's 2002 annual mean is scored; a 2002 model state is being compared with a multi-decadal observational composite, which is part of why the near-zero deep biases of CESM2 and GFDL-ESM4 should not be over-interpreted as exact skill.
