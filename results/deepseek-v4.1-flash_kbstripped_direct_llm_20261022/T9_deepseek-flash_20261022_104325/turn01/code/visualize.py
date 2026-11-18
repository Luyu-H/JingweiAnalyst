# -----------------------------------------------------------------------------
# Figure 1: PDO regression patterns
# -----------------------------------------------------------------------------
panel_order = ["ERSSTv5", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

pat_abs_max = max(
    np.nanmax(np.abs(all_results[name]["pattern"])) for name in panel_order
)
vmax = np.ceil(pat_abs_max * 10.0) / 10.0

fig1, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True, sharey=True)
axes = axes.ravel()

for ax, name in zip(axes, panel_order):
    res = all_results[name]
    mesh = ax.pcolormesh(
        target_lon,
        target_lat,
        res["pattern"],
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        shading="auto",
    )

    ax.set_title(name, fontsize=12)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xticks([110, 140, 170, 200, 230, 260])
    ax.set_xticklabels(["110E", "140E", "170E", "160W", "130W", "100W"])
    ax.set_yticks([20, 30, 40, 50, 60, 70])
    ax.grid(alpha=0.25)

    text = f"Var explained: {res['var_explained']:.1f}%"
    if name != "ERSSTv5":
        text += f"\nPattern corr: {res['pattern_corr_vs_ERSSTv5']:.3f}"
    ax.text(
        0.02,
        0.98,
        text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox=dict(facecolor="white", alpha=0.75, edgecolor="none"),
    )

cbar = fig1.colorbar(mesh, ax=axes.tolist(), shrink=0.9, pad=0.02)
cbar.set_label("°C per standard deviation")
fig1.suptitle("PDO regression patterns, 1950–2014", fontsize=14)
fig1.tight_layout()

# -----------------------------------------------------------------------------
# Figure 2: Standardized PDO indices
# -----------------------------------------------------------------------------
fig2, ax2 = plt.subplots(figsize=(13, 5))

colors = {
    "ERSSTv5": "C0",
    "CESM2": "C1",
    "GFDL-ESM4": "C2",
    "MPI-ESM1-2-LR": "C3",
}

for name in panel_order:
    res = all_results[name]
    ax2.plot(
        res["time_frac"],
        res["index"],
        label=name,
        color=colors[name],
        lw=1.0,
        alpha=0.85,
    )

ax2.plot(
    pub_time_frac,
    pub_std.values,
    label="Published PDO (ERSSTv6)",
    color="black",
    lw=1.5,
    ls="--",
)

ax2.plot(
    smooth_time_frac,
    smooth_ersst.values,
    label="ERSSTv5 9-yr running mean",
    color="C0",
    lw=3.0,
)

for yr in transition_years:
    ax2.axvline(yr, color="grey", ls=":", lw=1.0)
    ax2.text(
        yr,
        ax2.get_ylim()[1] * 0.92,
        str(yr),
        rotation=90,
        va="top",
        ha="right",
        fontsize=8,
        color="dimgrey",
    )

ax2.set_xlim(1950, 2015)
ax2.set_xlabel("Year")
ax2.set_ylabel("Standardized index")
ax2.set_title("PDO indices, 1950–2014")
ax2.grid(alpha=0.3)
ax2.legend(ncol=2, fontsize=8, loc="lower left")
fig2.tight_layout()
