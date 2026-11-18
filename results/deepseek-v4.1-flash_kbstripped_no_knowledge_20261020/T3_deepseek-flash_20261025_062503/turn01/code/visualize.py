# ---------------------------------------------------------------------
# Figure 1: 2 x 2 panel of SST bias maps + WOA23 reference SST
# ---------------------------------------------------------------------
pooled = np.concatenate([np.abs(np.asarray(bias_maps[m].values)).ravel() for m in models])
vmax = float(np.nanpercentile(pooled, 98))

fig, axes = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)

panel_specs = [
    (axes[0, 0], models[0], bias_maps[models[0]]),
    (axes[0, 1], models[1], bias_maps[models[1]]),
    (axes[1, 0], models[2], bias_maps[models[2]]),
]

for ax, name, field in panel_specs:
    pc = ax.pcolormesh(
        lon1d, lat1d, np.asarray(field.values),
        cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto",
    )
    ax.set_title(f"{name} - WOA23 SST bias", fontsize=11)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

# shared diverging colourbar for the three bias panels
cb = fig.colorbar(pc, ax=[axes[0, 0], axes[0, 1], axes[1, 0]],
                  orientation="vertical", shrink=0.85, extend="both")
cb.set_label("SST bias (\u00b0C)")

# WOA23 reference panel (sequential)
axw = axes[1, 1]
pcw = axw.pcolormesh(
    lon1d, lat1d, np.asarray(woa_map.values),
    cmap="viridis", shading="auto",
)
axw.set_title("WOA23 SST, 0 m (1971-2000)", fontsize=11)
axw.set_xlabel("Longitude")
axw.set_ylabel("Latitude")
axw.set_xlim(-180, 180)
axw.set_ylim(-90, 90)
cbw = fig.colorbar(pcw, ax=axw, orientation="vertical", shrink=0.85)
cbw.set_label("SST (\u00b0C)")

fig.suptitle(
    "SST climatology 1971-2000 on the common WOA23 comparison cells\n"
    f"(bias colour scale shared, clipped at the 98th percentile of |bias| = {vmax:.2f} \u00b0C)",
    fontsize=12,
)

# ---------------------------------------------------------------------
# Figure 2: global mean SST per product, each over its own native grid
# ---------------------------------------------------------------------
names = models + ["WOA23"]
values = [native_global_mean[m] for m in models] + [woa_global_mean]
colors = ["#4C72B0", "#DD8452", "#55A868", "#8172B3"]

fig2, ax2 = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
bars = ax2.bar(names, values, color=colors, edgecolor="black", linewidth=0.6)

for b, v in zip(bars, values):
    ax2.text(b.get_x() + b.get_width() / 2.0, b.get_height() + 0.05,
             f"{v:.2f}", ha="center", va="bottom", fontsize=10)

ax2.set_ylabel("Global mean SST (\u00b0C)")
ax2.set_xlabel("Product")
ax2.set_title("Global mean SST climatology 1971-2000, each on its own native ocean grid")
ax2.set_ylim(min(values) - 0.8, max(values) + 0.8)
ax2.grid(axis="y", linestyle=":", alpha=0.6)
plt.show()
