The nitracline was defined as the depth of the maximum vertical nitrate gradient within the top 1000 m, using the level just below 1000 m for the deepest gradient. CMIP6 values are 1971–2000 climatologies from each model’s native grid; global model means use each model’s own `areacello`, and WOA23 global means use cosine-latitude weighting. Basin means use the WOA23 basin mask with marginal seas folded into their parent oceans; Southern and Arctic Oceans are excluded.

Negative differences below mean the simulated nitracline is shallower than WOA23; positive means deeper.

### Global mean nitracline depth

| Source | Global mean depth (m) | Difference from WOA23 (m) |
|---|---:|---:|
| WOA23 | 156 | — |
| CESM2 | 140 | -16 |
| GFDL-ESM4 | 117 | -40 |
| MPI-ESM1-2-LR | 124 | -32 |

All three models place the global-mean nitracline too shallow relative to WOA23. CESM2 is closest in the global mean, although MPI-ESM1-2-LR has the best spatial statistics.

### Basin-mean nitracline depth

| Source | Pacific (m) | Atlantic (m) | Indian (m) |
|---|---:|---:|---:|
| WOA23 | 196 | 141 | 112 |
| CESM2 | 153 | 96 | 161 |
| GFDL-ESM4 | 138 | 72 | 127 |
| MPI-ESM1-2-LR | 173 | 90 | 83 |

Differences from WOA23:

| Model | Pacific difference (m) | Atlantic difference (m) | Indian difference (m) |
|---|---:|---:|---:|
| CESM2 | -43 | -46 | +49 |
| GFDL-ESM4 | -58 | -69 | +15 |
| MPI-ESM1-2-LR | -24 | -51 | -29 |

The models are generally too shallow in the Pacific and Atlantic. The Indian Ocean bias differs strongly among models: CESM2 and GFDL-ESM4 place the Indian nitracline too deep, whereas MPI-ESM1-2-LR places it too shallow.

### Spatial comparison with WOA23 on the 1° grid

| Model | Pearson correlation | RMSE (m) |
|---|---:|---:|
| CESM2 | 0.356 | 150.9 |
| GFDL-ESM4 | 0.457 | 148.4 |
| MPI-ESM1-2-LR | 0.581 | 139.2 |

### Best model and basin ordering

MPI-ESM1-2-LR reproduces the observed nitracline best overall: it has the highest spatial correlation with WOA23 (0.581) and the lowest RMSE (139.2 m). However, its absolute nitracline depths are still too shallow in all three basins.

The observed WOA23 basin ordering is:

**Pacific > Atlantic > Indian**

Only **MPI-ESM1-2-LR** reproduces this ordering exactly.  
CESM2 gives **Indian > Pacific > Atlantic**, so it does not match.  
GFDL-ESM4 gives **Pacific > Indian > Atlantic**, so it also does not match.

The generated figures are consistent with these diagnostics: Figure 1 shows the global nitracline-depth maps on a common 0–1000 m scale with deeper values visually heavier, and Figure 2 shows equatorial Pacific nitrate sections averaged over 2°S–2°N from 120°E to 80°W, with the nitracline depth overlaid in each panel.
