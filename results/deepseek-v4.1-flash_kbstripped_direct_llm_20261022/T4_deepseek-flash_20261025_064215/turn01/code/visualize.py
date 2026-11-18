# Figure 1: operator MLD minus published mlotst, on a common regular map grid.
fig1, axes1 = plt.subplots(3, 2, figsize=(12, 11), sharex=True, sharey=True)

for i, model in enumerate(MODELS):
    for j, month in enumerate(MONTHS):
        ax = axes1[i, j]
        data = fig1_src[(model, month)]
        im = ax.pcolormesh(
            TARGET_LON,
            TARGET_LAT,
            data,
            cmap="RdBu_r",
            vmin=-fig1_vmax,
            vmax=fig1_vmax,
            shading="auto",
        )
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        ax.set_xticks([-180, -120, -60, 0, 60, 120, 180])
        ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
        ax.grid(True, linewidth=0.2, alpha=0.4)
        ax.set_title(f"{model} — {MONTH_NAMES[month]}")
        if i == 2:
            ax.set_xlabel("Longitude")
        if j == 0:
            ax.set_ylabel("Latitude")

fig1.suptitle("Operator MLD minus published mlotst (m)", fontsize=14)
cbar = fig1.colorbar(
    im,
    ax=axes1.ravel().tolist(),
    orientation="horizontal",
    fraction=0.04,
    pad=0.06,
)
cbar.set_label("Operator MLD - mlotst (m)")

# Figure 2: operator MLD vs published mlotst scatter, common logarithmic axes.
fig2, axes2 = plt.subplots(3, 2, figsize=(11, 12), sharex=True, sharey=True)

for i, model in enumerate(MODELS):
    for j, month in enumerate(MONTHS):
        ax = axes2[i, j]
        x, y, bias, rmsd, corr = fig2_data[(model, month)]

        ax.scatter(x, y, s=1, alpha=0.15, edgecolors="none", rasterized=True)
        ax.plot(fig2_lims, fig2_lims, "k--", linewidth=0.9)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(fig2_lims)
        ax.set_ylim(fig2_lims)
        ax.set_xticks(fig2_ticks)
        ax.set_yticks(fig2_ticks)
        ax.set_xticklabels([f"{t:g}" for t in fig2_ticks])
        ax.set_yticklabels([f"{t:g}" for t in fig2_ticks])
        ax.minorticks_off()
        ax.grid(True, which="major", linewidth=0.2, alpha=0.4)

        ax.set_title(f"{model} — {MONTH_NAMES[month]}")
        ax.text(
            0.04,
            0.96,
            f"bias = {bias:.1f} m\nRMSD = {rmsd:.1f} m\nr = {corr:.3f}",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=8,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8, edgecolor="none"),
        )

        if i == 2:
            ax.set_xlabel("Operator MLD (m)")
        if j == 0:
            ax.set_ylabel("Published mlotst (m)")

fig2.suptitle("Operator MLD vs published mlotst", fontsize=14)
