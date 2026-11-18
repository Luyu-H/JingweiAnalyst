# ---------------- Figure 1: bias maps
pooled = np.concatenate([np.abs(bias_maps[m].values).ravel() for m in models])
vmax1 = float(np.nanpercentile(pooled, 98))
fig1, axes = plt.subplots(1, 3, figsize=(18, 4.5), constrained_layout=True, sharey=True)
for ax, m in zip(axes, models):
    b = bias_maps[m]
    pm = ax.pcolormesh(b["lon"], b["lat"], b.values, cmap="RdBu_r", vmin=-vmax1, vmax=vmax1, shading="auto")
    ax.set_title(f"{m}\nbias {bias_table.set_index('model').loc[m, 'bias_vs_WOA23']:+.2f} psu, "
                 f"RMSE {bias_table.set_index('model').loc[m, 'RMSE']:.2f}", color=colors[m])
    ax.set_xlabel("Longitude (°E)")
    ax.set_aspect("equal")
axes[0].set_ylabel("Latitude (°N)")
cb = fig1.colorbar(pm, ax=axes, orientation="vertical", extend="both", shrink=0.8)
cb.set_label("SSS bias (psu)")
fig1.suptitle("SSS bias, model 1971-2000 minus WOA23 (shared 1° ocean mask; colour scale clipped at 98th percentile of |bias|)")

# ---------------- Figure 2: trend maps
pooled2 = np.concatenate([np.abs(trend_maps[m].values).ravel() for m in models])
vmax2 = float(np.nanpercentile(pooled2, 98))
fig2, axes = plt.subplots(1, 3, figsize=(18, 4.5), constrained_layout=True, sharey=True)
for ax, m in zip(axes, models):
    t = trend_maps[m]
    pm = ax.pcolormesh(t["lon"], t["lat"], t.values, cmap="RdBu_r", vmin=-vmax2, vmax=vmax2, shading="auto")
    sg = sig_maps[m].values.astype(float)
    sg = np.where(np.isnan(t.values), 0.0, sg)
    ax.contourf(t["lon"], t["lat"], sg, levels=[0.5, 1.5], colors="none", hatches=["..."])
    ax.set_title(m, color=colors[m])
    ax.set_xlabel("Longitude (°E)")
    ax.set_aspect("equal")
axes[0].set_ylabel("Latitude (°N)")
cb = fig2.colorbar(pm, ax=axes, orientation="vertical", extend="both", shrink=0.8)
cb.set_label("SSS trend (psu/decade)")
fig2.suptitle("Annual-mean SSS trend 1950-2014 (stippled: p < 0.05; scale clipped at 98th percentile)")

# ---------------- Figure 3: global-mean series
fig3, ax = plt.subplots(figsize=(12, 5), constrained_layout=True)
for m in models:
    c = colors[m]
    yh = gm_hist[m]["year"].values
    ys = gm_ssp[m]["year"].values
    ax.plot(yh, gm_hist[m].values, color=c, lw=1.2, label=m)
    ax.plot(ys, gm_ssp[m].values, color=c, lw=1.2)
    ax.plot(yh, np.asarray(fits[m]["hist"]["fitted_data"]), color=c, lw=2.2, ls="--")
    ax.plot(ys, np.asarray(fits[m]["ssp"]["fitted_data"]), color=c, lw=2.2, ls="--")
ax.set_xlim(1948, 2102)
ax.set_xlabel("Year")
ax.set_ylabel("Global-mean SSS (psu)")
ax.set_title("Global-mean annual SSS: historical 1950-2014 and ssp245 2071-2100 (dashed: segment linear fits; 2015-2070 not shown)")
ax.legend(loc="best")
ax.grid(alpha=0.3)
