# Visualization code only
# Figure 1: 3x2 maps of operator MLD minus mlotst on a common 1-degree grid.
# Figure 2: 3x2 scatter plots of operator MLD versus mlotst on native shared-mask cells.

fig1, axes1 = plt.subplots(
    3, 2, figsize=(13, 10), constrained_layout=True, sharex=True, sharey=True
)

cmap = plt.get_cmap("RdBu_r")
norm = plt.Normalize(-vmax, vmax)

for i, model_name in enumerate(model_order):
    for j, month_name in enumerate(months):
        ax = axes1[i, j]
        field = results[model_name]["diff_maps"][month_name]
        ax.pcolormesh(
            target_lon_edges,
            target_lat_edges,
            field,
            shading="flat",
            cmap=cmap,
            norm=norm,
        )
        ax.set_xlim(0.0, 360.0)
        ax.set_ylim(-90.0, 90.0)
        ax.set_xticks(np.arange(0.0, 361.0, 60.0))
        ax.set_yticks(np.arange(-90.0, 91.0, 30.0))
        ax.set_aspect("equal")
        ax.set_title(f"{model_name}  {month_name}", fontsize=10)
        if i == 2:
            ax.set_xlabel("Longitude (°E)")
        if j == 0:
            ax.set_ylabel("Latitude (°N)")
        ax.grid(True, linewidth=0.3, alpha=0.5)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig1.colorbar(
    sm,
    ax=axes1,
    orientation="vertical",
    fraction=0.025,
    pad=0.02,
    extend="both",
)
cbar.set_label("MLD difference (m)")

fig1.suptitle(
    "Operator MLD minus mlotst, 1971-2000\n"
    "conservative 1° area-weighted cell-centre regrid for display; "
    f"shared scale ±{vmax:.0f} m (clipped at pooled 98th percentile)",
    fontsize=11,
)

fig2, axes2 = plt.subplots(
    3, 2, figsize=(11, 11), constrained_layout=True, sharex=True, sharey=True
)

for i, model_name in enumerate(model_order):
    for j, month_name in enumerate(months):
        ax = axes2[i, j]
        x, y = results[model_name]["scatter"][month_name]
        ax.scatter(x, y, s=1.2, alpha=0.30, rasterized=True, color="tab:blue")
        ax.plot(
            [scatter_lim_min, scatter_lim_max],
            [scatter_lim_min, scatter_lim_max],
            "k--",
            linewidth=1.0,
        )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(scatter_lim_min, scatter_lim_max)
        ax.set_ylim(scatter_lim_min, scatter_lim_max)
        ax.set_xticks(scatter_ticks)
        ax.set_yticks(scatter_ticks)
        ax.set_xticklabels([f"{t:g}" for t in scatter_ticks])
        ax.set_yticklabels([f"{t:g}" for t in scatter_ticks])
        ax.grid(True, which="both", linestyle=":", linewidth=0.4, alpha=0.5)
        ax.set_aspect("equal", adjustable="box")

        met = results[model_name]["metrics"][month_name]
        ax.set_title(
            f"{model_name}  {month_name}\n"
            f"bias = {met['bias']:.1f} m, "
            f"RMS = {met['rms']:.1f} m, "
            f"r = {met['corr']:.3f}",
            fontsize=9,
        )

        if i == 2:
            ax.set_xlabel("Operator MLD (m)")
        if j == 0:
            ax.set_ylabel("mlotst (m)")

fig2.suptitle(
    "Operator MLD versus published mlotst, native shared-mask ocean cells",
    fontsize=11,
)
