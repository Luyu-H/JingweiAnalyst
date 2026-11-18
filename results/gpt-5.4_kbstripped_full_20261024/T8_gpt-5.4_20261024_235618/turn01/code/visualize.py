# Visualization code only

# Figure 1: three global maps of 0-50 m DIC bias against GLODAP
pooled = np.concatenate([
    dic_bias_maps[m].values.ravel() for m in models
])
vmax = float(np.nanpercentile(np.abs(pooled), 98))

fig1, axes = plt.subplots(1, 3, figsize=(18, 4.8), constrained_layout=True)
for ax, model in zip(axes, models):
    da = dic_bias_maps[model]
    pcm = ax.pcolormesh(
        da["lon"].values,
        da["lat"].values,
        da.values,
        shading="auto",
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
    )
    ax.set_title(model)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

cbar = fig1.colorbar(pcm, ax=axes, orientation="vertical", shrink=0.9, extend="both")
cbar.set_label("DIC bias (µmol kg$^{-1}$)")
fig1.suptitle("Figure 1. 0-50 m DIC bias vs GLODAP on shared mask (shared scale, 98th-percentile clipped)")

# Figure 2: vertical profile of global mean DIC bias
fig2, ax = plt.subplots(figsize=(6.5, 5.5), constrained_layout=True)
for model in models:
    x = [depth_bias_profiles[model][d] for d in target_depths]
    y = target_depths
    ax.plot(x, y, marker="o", label=model)

ax.axvline(0, color="k", linewidth=1, linestyle="--")
ax.set_ylim(max(target_depths), min(target_depths))
ax.set_xlabel("DIC bias (µmol kg$^{-1}$)")
ax.set_ylabel("Depth (m)")
ax.set_title("Figure 2. Global mean DIC bias profile vs GLODAP")
ax.legend(loc="best")
