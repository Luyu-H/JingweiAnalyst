**Direct answer.** On the common WOA23 1° comparison mask, all three CMIP6 historical simulations reproduce the large-scale SST pattern very well (area-weighted spatial correlations 0.992–0.996). GFDL-ESM4 has the smallest absolute mean SST bias (−0.04 °C) and the lowest RMSE (0.894 °C). CESM2 has a warm mean bias of +0.30 °C, and MPI-ESM1-2-LR has the largest cold mean bias (−0.55 °C), largest RMSE (1.415 °C), and lowest pattern correlation (0.992). The ranking by smallest to largest absolute mean bias is GFDL-ESM4, CESM2, MPI-ESM1-2-LR.

**Method as executed.**
- CMIP6 historical monthly `tos` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR was read from the raw netCDF files, concatenated along time when needed, and CF-masked so `_FillValue`/`missing_value` and values ≥1e19 become NaN. Data were converted to °C if their units indicated Kelvin.
- The 1971–2000 climatology was formed by selecting all monthly time steps whose calendar year lies in 1971–2000 and taking the time mean over those 360 months.
- Each model’s native-grid global mean was computed as Σ(`tos_clim`·`areacello`)/Σ`areacello` with `areacello` restricted to cells where `tos_clim` is valid and area > 0.
- WOA23 was read from `woa23_decav71A0_t00_01.nc` using `t_an`, the first annual time slice, and the nearest 0 m depth level. The WOA23 land–sea mask was read from `landsea_01.msk` and placed on the 180×360 1° grid; ocean was defined as `landsea > 1`.
- The WOA23 whole-grid SST mean was computed as a cos(latitude)-weighted mean over valid WOA23 temperature cells that satisfy the ocean mask.
- For the model–WOA comparison, each model’s 1971–2000 climatology was conservatively binned to the WOA23 1° grid: each native cell centre was assigned to its containing 1° cell, and the target value was the `areacello`-weighted mean of contributing native cells.
- The common comparison mask required `landsea > 1`, valid WOA23 temperature, and valid regridded SST for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR simultaneously. It contains **21,961** 1° cells. All shared-mask means, biases, RMSEs and correlations were computed on this mask using cos(latitude) area weights.

**Global means on each product’s own grid.**
- CESM2: **18.52 °C**
- GFDL-ESM4: **18.13 °C**
- MPI-ESM1-2-LR: **17.58 °C**
- WOA23 0 m annual field, cos-latitude-weighted over the valid 1° ocean grid: **18.23 °C**

These native-grid means are not expected to equal the shared-mask means because the model means use each model’s own grid, own `areacello`, and own land–sea mask, while the shared-mask statistics use a common 1° mask and conservatively binned fields.

**Common-mask comparison on WOA23 1° grid (21,961 cells).**
The WOA23 shared-mask mean is **16.73 °C**.

| Model | Shared-mask mean SST (°C) | Mean bias model−WOA23 (°C) | Area-weighted RMSE (°C) | Area-weighted spatial correlation |
|---|---:|---:|---:|---:|
| CESM2 | 17.03 | +0.30 | 1.156 | 0.994 |
| GFDL-ESM4 | 16.70 | −0.04 | 0.894 | 0.996 |
| MPI-ESM1-2-LR | 16.19 | −0.55 | 1.415 | 0.992 |

**Ranking by absolute mean bias, smallest to largest.**
1. **GFDL-ESM4** (|−0.036| °C)
2. **CESM2** (|+0.303| °C)
3. **MPI-ESM1-2-LR** (|−0.545| °C)

**Spatial structure of the bias.** The bias fields are spatially heterogeneous rather than a uniform offset. The unweighted map ranges are approximately −4.88 to +8.61 °C for CESM2, −4.41 to +7.39 °C for GFDL-ESM4, and −7.45 to +7.68 °C for MPI-ESM1-2-LR, with map standard deviations of about 1.14, 0.90 and 1.34 °C, respectively. Thus the small area-weighted global-mean biases arise from cancellation between regional warm and cold biases, with local SST biases of several °C. All models still have very high spatial correlations with WOA23 because the dominant large-scale SST pattern is captured; the RMSE differences show that the local amplitude and placement of the biases differ substantially between models. MPI-ESM1-2-LR has the largest mapped bias amplitude and the largest RMSE. The comparison is on conservatively binned 1° fields, so the model-native resolution differences can affect the apparent amplitude of fine-scale features.

**Figure 1.** A 2×2 map panel on the shared WOA23 comparison cells. The first three panels show CESM2, GFDL-ESM4 and MPI-ESM1-2-LR SST biases (model minus WOA23) using one shared diverging `RdBu_r` colour scale centred on zero, clipped at the 98th percentile of the pooled absolute bias. The fourth panel shows the WOA23 1971–2000 0 m SST reference on the same cells using a sequential `viridis` colour scale, clipped to its 2nd–98th percentile range. Each panel is labelled by model/reference and axis, and each colour bar is labelled in °C.

**Figure 2.** A bar chart of the four global mean SST values: CESM2 (18.52 °C), GFDL-ESM4 (18.13 °C), MPI-ESM1-2-LR (17.58 °C) and WOA23 (18.23 °C). Each bar is annotated with its two-decimal value and the y axis is labelled “Mean SST (°C)”. The three model bars are on their own native ocean grids with their own `areacello`; the WOA23 bar is over the valid 1° grid.

**Interpretation and caveats.** GFDL-ESM4 performs best on the common mask for mean bias, RMSE and pattern correlation. CESM2 has a modest warm global-mean bias but higher local RMSE than GFDL-ESM4. MPI-ESM1-2-LR has the strongest cold global-mean bias and the largest RMSE. The native-grid global means and the shared-mask means differ by construction: WOA23’s whole-grid mean (18.23 °C) is higher than its shared-mask mean (16.73 °C), and the model native-grid means are likewise not directly comparable to the shared-mask means. Finally, `tos` is a top-cell bulk temperature and can include under-ice values near −1.8 °C rather than missing, whereas WOA23 sea-ice-covered cells are ice-free-water values; this can contribute to high-latitude differences in the maps and RMSE.
