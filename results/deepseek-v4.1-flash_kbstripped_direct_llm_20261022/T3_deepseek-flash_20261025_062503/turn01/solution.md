The 1971–2000 SST climatologies were computed from monthly CMIP6 historical `tos` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, then area-weighted and binned to the WOA23 1° grid for comparison with WOA23 at 0 m. The comparison uses the 21,961 common ocean cells where all three models and WOA23 are valid.

**Numerical results**

| Product | Native-grid global mean SST (°C) | Mean on common WOA grid (°C) | Bias vs WOA23 (°C) | RMSE (°C) | Pearson correlation |
|---|---:|---:|---:|---:|---:|
| CESM2 | 18.52 | 13.18 | +0.19 | 1.157 | 0.995 |
| GFDL-ESM4 | 18.13 | 12.98 | −0.02 | 0.905 | 0.997 |
| MPI-ESM1-2-LR | 17.58 | 12.61 | −0.39 | 1.399 | 0.993 |
| WOA23 | 13.85 | 13.00 | — | — | — |

The WOA23 value of 13.85 °C is the simple mean over valid ocean cells on the 1° grid. For context, the cosine-latitude area-weighted WOA23 mean is 18.23 °C, so the large WOA–model offset in Fig. 2 is mainly a native-grid weighting/definition artifact. The like-for-like comparison is on the common WOA grid, where the models range from 12.61 to 13.18 °C against WOA23 at 13.00 °C.

**Ranking by absolute mean bias, smallest to largest**
1. GFDL-ESM4: |ΔSST| = 0.02 °C
2. CESM2: |ΔSST| = 0.19 °C
3. MPI-ESM1-2-LR: |ΔSST| = 0.39 °C

**Spatial structure of the bias**
All three models reproduce the broad observed SST pattern very well, with spatial correlations of 0.993–0.997. GFDL-ESM4 is clearly the best overall: it has the smallest mean bias, the lowest RMSE (0.905 °C) and the highest correlation. CESM2 has a weak warm mean bias (+0.19 °C) and slightly larger RMSE (1.157 °C). MPI-ESM1-2-LR has the strongest cool mean bias (−0.39 °C), the largest bias variability and the largest RMSE (1.399 °C).

The bias fields are spatially coherent rather than random. The largest local errors occur in dynamically active and coastal/high-latitude regions—especially high-latitude/sea-ice margins and eastern-boundary/equatorial upwelling and current-extension regions—while much of the subtropical open ocean is close to WOA23. The models show alternating positive and negative local bias patches, so their global mean biases can be small even when local biases are large: the model bias fields reach approximately +8.6/−4.9 °C for CESM2, +7.4/−4.4 °C for GFDL-ESM4, and +7.7/−7.4 °C for MPI-ESM1-2-LR. MPI-ESM1-2-LR has the strongest spatial heterogeneity, consistent with its largest RMSE.

Figure 1 displays the three model-minus-WOA23 bias maps on the common comparison cells with a shared diverging colour scale, plus the WOA23 0 m SST reference with a sequential scale. Figure 2 shows the four global mean SST values with each bar annotated in °C.
