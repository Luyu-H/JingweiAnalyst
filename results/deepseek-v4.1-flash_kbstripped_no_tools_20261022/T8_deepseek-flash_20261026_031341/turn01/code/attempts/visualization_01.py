# ---------------------------------------------------------------------------
# Figure 1: 0-50 m DIC bias maps (model minus GLODAP) on the common mask
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(
    1, 3, figsize=(16, 4.5), constrained_layout=True
)

pooled = np.concatenate(
    [
        np.abs(map_metrics[m]["diff"].values.ravel())
        for m in models
        if map_metrics[m]["diff"] is not None
    ]
)
vmax = float(np.nanpercentile(pooled, 98))
if not np.isfinite(vmax) or vmax == 0:
    vmax = 1.0

pcm = None
for ax, m in zip(axes, models):
    diff = map_metrics[m]["diff"]
    pcm = ax.pcolormesh(
        glodap_lon,
        glodap_lat,
        diff.values,
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        shading="auto",
    )
    ax.set_title(m)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(glodap_lon.min(), glodap_lon.max())
    ax.set_ylim(glodap_lat.min(), glodap_lat.max())

cbar = fig.colorbar(
    pcm,
    ax=axes,
    orientation="horizontal",
    fraction=0.06,
    pad=0.06,
    extend="both",
)
cbar.set_label("DIC bias (µmol kg$^{-1}$)")
fig.suptitle(
    "0-50 m DIC bias vs GLODAP (2002 annual mean; shared scale, "
    "clipped at 98th percentile of |bias|)",
    fontsize=12,
)

# ---------------------------------------------------------------------------
# Figure 2: vertical profile of global mean DIC bias
# ---------------------------------------------------------------------------
fig2, ax2 = plt.subplots(figsize=(6, 5), constrained_layout=True)

for m in models:
    ax2.plot(
        model_results[m]["depth_biases"],
        depth_targets,
        marker="o",
        label=m,
    )

ax2.axvline(0.0, color="k", linestyle="--", linewidth=1)
ax2.invert_yaxis()
ax2.set_yticks(depth_targets)
ax2.set_ylabel("Depth (m)")
ax2.set_xlabel("DIC bias (µmol kg$^{-1}$)")
ax2.set_title("Global mean DIC bias vs GLODAP")
ax2.legend(loc="best")
ax2.grid(True, linestyle=":", alpha=0.5)
