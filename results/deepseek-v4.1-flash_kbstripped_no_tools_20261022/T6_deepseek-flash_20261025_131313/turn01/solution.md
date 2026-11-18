**Bottom line.** On their native grids, all three CMIP6 models place the global-mean nitracline shallower than WOA23. WOA23 gives **156 m**; CESM2 **140 m**, GFDL-ESM4 **117 m**, and MPI-ESM1-2-LR **124 m**. On the common WOA 1° comparison mask, the model-minus-WOA differences are **CESM2 −14 m**, **GFDL-ESM4 −38 m**, and **MPI-ESM1-2-LR −32 m**. Using the stated criterion of highest area-weighted pattern correlation and lowest RMSE, **MPI-ESM1-2-LR reproduces the observed nitracline best** (r = 0.566, RMSE = 151.7 m), although CESM2 has the smallest global mean bias. CESM2 and MPI-ESM1-2-LR reproduce WOA23’s basin ordering; GFDL-ESM4 does not.

**Methodology as executed.**  
- **WOA23 nitrate:** used the objectively analysed all-time-period annual climatology `n_an` from `woa23_all_n00_01.nc`, on the 1° grid (180×360, 102 depth levels), in µmol kg⁻¹. The WOA `landsea_01.msk` mask was used with ocean defined as value > 1; the WOA `basinmask_01.msk` surface value column was used for basin labels. Southern and Arctic were excluded. No Ofx basin field was used.  
- **CMIP6 data:** historical monthly `no3`, `thetao`, `so`, and `areacello` for CESM2, GFDL-ESM4, and MPI-ESM1-2-LR over **1971-01-01 to 2000-12-31**, each on its native horizontal grid and native vertical levels. The time mean was a custom **month-length-weighted** mean using `days_in_month`. Native grids used were CESM2 384×320×60, GFDL-ESM4 180×360×35, and MPI-ESM1-2-LR 220×256×40.  
- **Nitracline definition:** adjacent-level finite differences on each native vertical grid; only pairs whose upper level is at or above 1000 m were searched, with the next deeper level allowed to close the deepest pair. The nitracline depth is the depth of the upper level of the maximum vertical gradient. No vertical interpolation was used to locate the nitracline. Models were located on their native archived **mol m⁻³** `no3`; WOA23 was located on its native **µmol kg⁻¹** field.  
- **Units/density:** model `no3` was converted to µmol kg⁻¹ for sections and shared-unit plotting using UNESCO 1983 in-situ density from time-mean `thetao` and `so`, with pressure approximated as **1.0055 dbar per metre depth**. The recorded density assumption is: “UNESCO 1983 in-situ density from time-mean thetao and so; pressure approximated as 1.0055 dbar per metre depth.”  
- **Means:** native global and basin means are area-weighted, using WOA cos(latitude) cell areas for WOA23 and each model’s own `areacello` for models, with weights masked to the same valid cells. Basin labels were mapped to each model grid by nearest neighbour in 3-D longitude/latitude space.  
- **Comparison:** model nitracline fields were conservatively mapped to the WOA 1° grid by areacello-weighted binning with nearest filling. CESM2 went from 384×320 to 180×360 (coarser in latitude, finer in longitude); GFDL-ESM4 was already 180×360; MPI-ESM1-2-LR went from 220×256 to 180×360 (finer target). A single common mask was applied to all four fields, giving **41,071 comparison cells**. Correlations, RMSE, and differences were area-weighted on that shared sample; reported model-minus-WOA differences use the shared-sample means, not differences of native-grid means.  
- **Basin caveat:** the executed classifier retained only one final WOA label per basin: Atlantic = 5, Pacific = 12, Indian = 56, with Southern = 10 and Arctic = 11 excluded. It did **not** fully fold every marginal-sea label into its parent ocean. The basin means/differences/order below therefore refer to those selected WOA labels, not a complete marginal-sea-folded basin definition.

**Global nitracline depth**
| Source | Native-grid global mean (m) | Shared-sample mean (m) | Shared-sample difference vs WOA23 (m) |
|---|---:|---:|---:|
| WOA23 | 156 | 156.3 | — |
| CESM2 | 140 | 142.7 | −13.6 (−14) |
| GFDL-ESM4 | 117 | 118.5 | −37.8 (−38) |
| MPI-ESM1-2-LR | 124 | 124.4 | −31.9 (−32) |

**Native-grid basin mean nitracline depth (m)**
| Source | Pacific | Atlantic | Indian |
|---|---:|---:|---:|
| WOA23 | 78 | 51 | 75 |
| CESM2 | 82 | 31 | 71 |
| GFDL-ESM4 | 42 | 30 | 58 |
| MPI-ESM1-2-LR | 98 | 31 | 59 |

**Shared-sample basin differences, model minus WOA23 (m)**
| Model | Pacific | Atlantic | Indian |
|---|---:|---:|---:|
| CESM2 | +6 | −9 | −2 |
| GFDL-ESM4 | −38 | −13 | −15 |
| MPI-ESM1-2-LR | +21 | −15 | −15 |

**Pattern metrics on the common WOA 1° mask (41,071 cells)**
| Model | Area-weighted Pearson r | RMSE (m) |
|---|---:|---:|
| CESM2 | 0.371 | 159.7 |
| GFDL-ESM4 | 0.454 | 161.4 |
| MPI-ESM1-2-LR | 0.566 | 151.7 |

**Basin ordering, shallow to deep**
- WOA23: Atlantic (51 m) < Indian (75 m) < Pacific (78 m).
- CESM2: Atlantic (31 m) < Indian (71 m) < Pacific (82 m) — matches WOA23.
- GFDL-ESM4: Atlantic (30 m) < Pacific (42 m) < Indian (58 m) — does not match WOA23.
- MPI-ESM1-2-LR: Atlantic (31 m) < Indian (59 m) < Pacific (98 m) — matches WOA23.

**Visualizations.** Figure 1 is a 2×2 global map panel of nitracline depth for WOA23, CESM2, GFDL-ESM4, and MPI-ESM1-2-LR on the common-masked WOA 1° grid, using one shared sequential colour scale (0 to the pooled 98th percentile, with deeper values visually heavier) and land masked. Figure 2 is a four-panel longitude–depth section of nitrate through the equatorial Pacific (2°S–2°N mean), covering 120°E–80°W and the top 1000 m with depth increasing downward, in the order WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR. The model sections are in µmol kg⁻¹ after density conversion and are display-interpolated to a common longitude grid and WOA depth levels; all panels share one concentration colour scale. Each panel overlays its own nitracline depth line obtained with the same gradient rule from the displayed 2°S–2°N-mean section. These display interpolations do not affect the metrics.

**Interpretation.** The models capture only part of the observed nitracline spatial structure: correlations are modest (0.371–0.566) and RMSEs are large (151.7–161.4 m). MPI-ESM1-2-LR is best by the specified correlation/RMSE criterion, but it has a substantial shallow global bias (−32 m), a deep Pacific bias (+21 m), and shallow Atlantic/Indian biases (−15 m each). CESM2 has the smallest global bias (−14 m) and the best basin-order match with relatively small Pacific/Indian biases, but its spatial correlation is lower and RMSE higher. GFDL-ESM4 has the largest shallow global bias (−38 m), the largest RMSE, and reverses the Pacific–Indian ordering. Native global means are deeper than the three-basin means because the global calculation includes the Southern and Arctic Oceans, which are excluded from the Pacific/Atlantic/Indian basin means. Nitracline depths are quantized by each product’s native vertical grid, so the maps should be read together with the underlying vertical gradient. The WOA23 field is an all-time-period composite, not a 1971–2000 normal, so a period mismatch is inherent in the comparison. The basin results should be interpreted with the above label-selection caveat.
