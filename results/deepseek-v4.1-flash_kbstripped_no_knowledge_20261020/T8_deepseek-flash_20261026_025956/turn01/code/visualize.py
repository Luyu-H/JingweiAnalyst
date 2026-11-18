# ----------------------------------------------------------------------
# Figure 1: 0-50 m DIC bias maps (model minus GLODAP)
# ----------------------------------------------------------------------
fig1, axes = plt.subplots(
    1, 3, figsize=(15, 4), constrained_layout=True
)

all_bias_values = np.concatenate(
    [
        model_results[model]["bias_map"].values.flatten()
        for model in models
    ]
)
vmax = float(np.nanpercentile(np.abs(all_bias_values), 98))

for ax, model in zip(axes, models):
    bias = model_results[model]["bias_map"]
    lon = bias.lon
    lat = bias.lat
    im = ax.pcolormesh(
        lon,
        lat,
        bias,
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        shading="auto",
    )
    ax.set_title(
        f"{model}\n0-50 m DIC bias (clipped at 98th pct)",
        fontsize=10,
    )
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.axhline(0, color="gray", lw=0.5, ls="--")

fig1.colorbar(
    im,
    ax=axes,
    orientation="horizontal",
    fraction=0.05,
    pad=0.05,
    extend="both",
    label="DIC bias (µmol kg$^{-1}$)",
)
fig1.suptitle(
    "Figure 1: 0-50 m DIC bias against GLODAP (model - GLODAP)",
    fontsize=12,
)

# ----------------------------------------------------------------------
# Figure 2: Vertical profile of global mean DIC bias
# ----------------------------------------------------------------------
fig2, ax2 = plt.subplots(figsize=(6, 5))

for model in models:
    biases = [
        model_results[model]["vertical_bias"][d] for d in target_depths
    ]
    ax2.plot(biases, target_depths, marker="o", label=model)

ax2.axvline(0, color="k", lw=0.8, ls="--")
ax2.set_yticks(target_depths)
ax2.invert_yaxis()
ax2.set_xlabel("DIC bias (µmol kg$^{-1}$)")
ax2.set_ylabel("Depth (m)")
ax2.set_title("Figure 2: Vertical profile of global mean DIC bias")
ax2.legend(loc="best")
ax2.grid(True, alpha=0.3)
