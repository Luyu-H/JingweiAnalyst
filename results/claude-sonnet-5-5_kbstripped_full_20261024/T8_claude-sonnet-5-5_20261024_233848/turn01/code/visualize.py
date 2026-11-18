fig, axes = plt.subplots(1, 3, figsize=(18, 4.8), constrained_layout=True, sharey=True)
pooled = np.concatenate([np.abs(bias_maps[m].values[np.isfinite(bias_maps[m].values)]).ravel() for m in models])
vmax = float(np.nanpercentile(pooled, 98))
im = None
for ax, m in zip(axes, models):
    b = bias_maps[m]
    im = ax.pcolormesh(b["lon"].values, b["lat"].values, b.values, cmap="RdBu_r",
                       vmin=-vmax, vmax=vmax, shading="auto")
    ax.set_title(f"{m}\nmean bias {map_metrics[m]['bias']:.0f}, RMSE {map_metrics[m]['rmse']:.0f}")
    ax.set_xlabel("Longitude (°E)")
axes[0].set_ylabel("Latitude (°N)")
cb = fig.colorbar(im, ax=axes, orientation="vertical", extend="both", shrink=0.9, pad=0.02)
cb.set_label("DIC bias (µmol kg-1)")
fig.suptitle("0-50 m DIC bias vs GLODAPv2 (model 2002 minus GLODAP), common cells; colour scale clipped at 98th percentile of |bias|")

fig2, ax2 = plt.subplots(figsize=(6.5, 7), constrained_layout=True)
colors = {"CESM2": "tab:blue", "GFDL-ESM4": "tab:orange", "MPI-ESM1-2-LR": "tab:green"}
for m in models:
    ax2.plot(prof_bias[m], PROF, "-o", color=colors[m], label=m, markersize=6)
ax2.axvline(0, color="k", lw=1)
for d in PROF:
    ax2.axhline(d, color="grey", lw=0.5, ls=":")
ax2.set_yticks(PROF)
ax2.invert_yaxis()
ax2.set_xlabel("Global mean DIC bias (µmol kg-1)")
ax2.set_ylabel("Depth (m)")
ax2.set_title("Global mean DIC bias profile (model 2002 minus GLODAP)")
ax2.legend(loc="best")
