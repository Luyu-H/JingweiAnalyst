# =============================================================================
# Figure 1: 0-50 m DIC bias maps (model - GLODAP) on the shared comparison mask
# =============================================================================
lon1d = bias_maps[models[0]].lon.values
lat1d = bias_maps[models[0]].lat.values

pooled = np.concatenate([np.abs(bias_maps[m].values).ravel() for m in models])
vmax = float(np.nanpercentile(pooled, 98))
if not np.isfinite(vmax) or vmax == 0:
    vmax = float(np.nanmax(pooled))

fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.6), constrained_layout=True)
for ax, m in zip(axes, models):
    d = bias_maps[m]
    pm = ax.pcolormesh(lon1d, lat1d, d.values, cmap="RdBu_r",
                       vmin=-vmax, vmax=vmax, shading="auto")
    ax.set_facecolor("0.85")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-80, 90)
    ax.set_xticks([-180, -120, -60, 0, 60, 120, 180])
    ax.set_yticks([-60, -30, 0, 30, 60])
    ax.set_title("%s\nbias %.0f, RMSE %.0f, r = %.3f" % (
        m, map_metrics[m]["bias"], map_metrics[m]["rmse"], map_metrics[m]["correlation"]),
        fontsize=10)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

fig.colorbar(pm, ax=axes, label="DIC bias (µmol kg$^{-1}$)", extend="both")
fig.suptitle("0–50 m DIC bias vs GLODAPv2.2016b (model − GLODAP), shared diverging scale "
             "clipped at the 98th percentile of |bias|; common mask: %d cells" % mask_n,
             fontsize=11)

# =============================================================================
# Figure 2: vertical profile of the global mean DIC bias (model - GLODAP)
# =============================================================================
fig2, ax2 = plt.subplots(figsize=(6.4, 5.2), constrained_layout=True)
markers = {"CESM2": "o", "GFDL-ESM4": "s", "MPI-ESM1-2-LR": "^"}
for m in models:
    x = [vert_bias[m][d] for d in profile_depths]
    ax2.plot(x, profile_depths, marker=markers.get(m, "o"), ms=7, lw=1.8, label=m)

ax2.axvline(0.0, color="k", lw=0.9, ls="--")
ax2.set_yticks(profile_depths)
ax2.set_yticklabels([str(d) for d in profile_depths])
ax2.invert_yaxis()
ax2.set_xlabel("DIC bias (µmol kg$^{-1}$)")
ax2.set_ylabel("Depth (m)")
ax2.set_title("Global mean DIC bias vs GLODAP (model − GLODAP)\non the GLODAP 1° grid, "
              "pairwise valid cells at each depth", fontsize=10)
ax2.grid(alpha=0.3)
ax2.legend(loc="best")
