# ---------------------------------------------------------------------------
# Figure 1: regression patterns over the PDO domain
# ---------------------------------------------------------------------------
order = ["ERSSTv5"] + models
fig1, axes = plt.subplots(2, 2, figsize=(13, 8))
axes = axes.ravel()

all_vals = np.concatenate(
    [amp_maps[s][np.isfinite(amp_maps[s])].ravel() for s in order]
)
vmax = np.nanmax(np.abs(all_vals))
levels = np.linspace(-vmax, vmax, 21)

cf = None
for ax, src in zip(axes, order):
    cf = ax.contourf(
        target_lon2d,
        target_lat2d,
        amp_maps[src],
        levels=levels,
        cmap="RdBu_r",
        extend="both",
    )
    ax.set_title(src)
    ax.set_xlabel("Longitude (0-360)")
    ax.set_ylabel("Latitude")
    ax.set_xlim(target_lon.min(), target_lon.max())
    ax.set_ylim(target_lat.min(), target_lat.max())

    text = f"Var explained = {100.0 * var_expl[src]:.1f}%"
    if src in pattern_corr:
        text += f"\nPattern corr = {pattern_corr[src]:.3f}"
    ax.text(
        0.02,
        0.97,
        text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"),
    )

fig1.colorbar(
    cf,
    ax=axes,
    orientation="horizontal",
    fraction=0.05,
    pad=0.08,
    label="°C per standard deviation",
)
fig1.suptitle("PDO regression patterns, 1950-2014", y=0.98)

# ---------------------------------------------------------------------------
# Figure 2: standardized PDO indices
# ---------------------------------------------------------------------------
fig2, ax = plt.subplots(figsize=(13, 5))

colors = {
    "ERSSTv5": "black",
    "CESM2": "tab:blue",
    "GFDL-ESM4": "tab:orange",
    "MPI-ESM1-2-LR": "tab:green",
}

for src in order:
    lw = 2.0 if src == "ERSSTv5" else 1.0
    alpha = 1.0 if src == "ERSSTv5" else 0.8
    ax.plot(
        dates,
        idx_all[src],
        label=src,
        color=colors[src],
        linewidth=lw,
        alpha=alpha,
    )

ax.plot(
    dates,
    pub_idx,
    label="Published PDO (NCEI, ERSSTv6)",
    color="tab:red",
    linestyle="--",
    linewidth=1.5,
)

ax.plot(
    smooth_obs.index,
    smooth_obs.values,
    label="ERSSTv5 9-yr centred mean",
    color="black",
    linewidth=2.8,
    alpha=0.8,
)

for yr in transition_years:
    x = pd.Timestamp(f"{yr}-01-01")
    ax.axvline(x, color="grey", linestyle=":", linewidth=1.0, alpha=0.8)
    ax.text(
        x,
        ax.get_ylim()[1] * 0.95,
        str(yr),
        rotation=90,
        va="top",
        ha="right",
        fontsize=8,
        color="dimgrey",
    )

ax.axhline(0.0, color="grey", linewidth=0.8, alpha=0.6)
ax.set_xlabel("Time")
ax.set_ylabel("Standardized index")
ax.set_title("PDO indices, 1950-2014")
ax.legend(loc="upper left", ncol=2, fontsize=8)
ax.grid(alpha=0.25)

fig2.tight_layout()
