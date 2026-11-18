# Figure 1: 0-50 m DIC bias maps, one panel per model, shared diverging colour scale.
fig1, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=True)

all_bias = np.concatenate(
    [bias_maps[m][np.isfinite(bias_maps[m])].ravel() for m in MODELS]
)
vmax = np.nanpercentile(np.abs(all_bias), 95)
if not np.isfinite(vmax) or vmax <= 0:
    vmax = 1.0

im = None
for ax, model in zip(axes, MODELS):
    im = ax.pcolormesh(
        target_lon,
        target_lat,
        bias_maps[model],
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        shading="auto",
    )
    ax.set_title(model)
    ax.set_xlabel("Longitude")
    ax.set_xlim([-180, 180])
    ax.set_ylim([-90, 90])

axes[0].set_ylabel("Latitude")
fig1.colorbar(
    im,
    ax=axes,
    orientation="horizontal",
    fraction=0.05,
    pad=0.15,
    label="DIC bias (µmol kg$^{-1}$)",
)
fig1.suptitle("0–50 m mean DIC bias vs GLODAP", y=1.02)

# Figure 2: vertical profile of global mean DIC bias.
fig2, ax2 = plt.subplots(figsize=(5.5, 6))
for model in MODELS:
    vals = [profile_bias[model][d] for d in TARGET_DEPTHS]
    ax2.plot(vals, TARGET_DEPTHS, marker="o", label=model)

ax2.axvline(0.0, color="k", lw=0.8)
for d in TARGET_DEPTHS:
    ax2.axhline(d, color="gray", lw=0.3, alpha=0.6)

ax2.invert_yaxis()
ax2.set_xlabel("DIC bias (µmol kg$^{-1}$)")
ax2.set_ylabel("Depth (m)")
ax2.set_title("Global mean DIC bias vs GLODAP")
ax2.legend()
