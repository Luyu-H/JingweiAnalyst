import numpy as np
import matplotlib.pyplot as plt

COL = {"CESM2": "#1f77b4", "GFDL-ESM4": "#d62728", "MPI-ESM1-2-LR": "#2ca02c"}

# ------------------------------------------------------------------ Figure 1
fig1, axes1 = plt.subplots(1, 3, figsize=(16, 4.0), constrained_layout=True)
pooled = np.concatenate([np.asarray(diff_maps[m].values).ravel() for m in MODELS])
vmax1 = float(np.nanpercentile(np.abs(pooled), 98))
for ax, m in zip(axes1, MODELS):
    f = diff_maps[m]
    pm = ax.pcolormesh(f["lon"].values, f["lat"].values, f.values,
                       cmap="RdBu_r", vmin=-vmax1, vmax=vmax1, shading="auto")
    ax.set_title("%s: bias vs WOA23\nr = %.3f, RMSE = %.2f psu"
                 % (m, corr_woa[m], rmse_woa[m]), fontsize=10)
    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
fig1.colorbar(pm, ax=axes1, label="SSS bias (psu)", extend="both")
fig1.suptitle("Figure 1 - 1971-2000 SSS bias against WOA23 on the common comparison cells "
              "(colours clipped at the 98th percentile of |bias|)", fontsize=11)

# ------------------------------------------------------------------ Figure 2
lon2d, lat2d = np.meshgrid(trend_map_slope[MODELS[0]]["lon"].values,
                           trend_map_slope[MODELS[0]]["lat"].values)
pooled_t = np.concatenate([np.asarray(trend_map_slope[m].values).ravel() for m in MODELS])
vmax2 = float(np.nanpercentile(np.abs(pooled_t), 98))

fig2, axes2 = plt.subplots(1, 3, figsize=(16, 4.0), constrained_layout=True)
for ax, m in zip(axes2, MODELS):
    f = trend_map_slope[m]
    pm2 = ax.pcolormesh(f["lon"].values, f["lat"].values, f.values,
                        cmap="RdBu_r", vmin=-vmax2, vmax=vmax2, shading="auto")
    sig = np.asarray(trend_map_sig[m].values) > 0.5
    ax.plot(lon2d[sig], lat2d[sig], ".", color="k", markersize=0.6, alpha=0.6)
    ax.set_title("%s: 1950-2014 SSS trend\n(stipple: p < 0.05)" % m, fontsize=10)
    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
fig2.colorbar(pm2, ax=axes2, label="SSS trend (psu/decade)", extend="both")
fig2.suptitle("Figure 2 - cell-by-cell linear trend of annual-mean SSS, 1950-2014 "
              "(colours clipped at the 98th percentile of |trend|)", fontsize=11)

# ------------------------------------------------------------------ Figure 3
fig3, ax3 = plt.subplots(figsize=(11, 5), constrained_layout=True)
for m in MODELS:
    h = gm_hist_year[m]
    s = gm_ssp_year[m]
    ax3.plot(h["year"].values, np.asarray(h.values), color=COL[m], lw=1.4,
             label=m)
    ax3.plot(s["year"].values, np.asarray(s.values), color=COL[m], lw=1.4)
    fh = trend_hist[m]["fitted_data"]
    fs = trend_ssp[m]["fitted_data"]
    ax3.plot(fh["year"].values, np.asarray(fh.values), color=COL[m], ls="--", lw=1.6)
    ax3.plot(fs["year"].values, np.asarray(fs.values), color=COL[m], ls="--", lw=1.6)

ax3.axvspan(2015, 2070, color="0.9", zorder=0)
ax3.text(2042, ax3.get_ylim()[1], "no data\n2015-2070", ha="center", va="top", fontsize=8)
ax3.set_xlabel("year")
ax3.set_ylabel("global mean SSS (psu)")
ax3.set_title("Figure 3 - annual global-mean SSS, historical 1950-2014 and ssp245 2071-2100\n"
              "dashed lines: fitted linear trends (per decade slopes in the printed table)",
              fontsize=10)
ax3.legend(loc="best", fontsize=9, frameon=False)
ax3.grid(alpha=0.3)

plt.show()
