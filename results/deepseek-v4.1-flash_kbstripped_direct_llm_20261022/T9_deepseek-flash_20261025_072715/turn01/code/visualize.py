# ---------------------------------------------------------------------
# Figure 1: PDO regression patterns
# ---------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True, sharey=True)

vmax = np.nanmax([
    np.nanmax(np.abs(results[name]["reg_grid"]))
    for name in source_order
])
if not np.isfinite(vmax) or vmax == 0:
    vmax = 1.0

for ax, name in zip(axes.flat, source_order):
    reg = results[name]["reg_grid"]
    im = ax.pcolormesh(
        lon_t, lat_t, reg,
        cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto"
    )
    ax.set_title(name)
    ax.set_xlabel("Longitude (0-360)")
    ax.set_ylabel("Latitude")
    ax.set_xlim(lon_t.min(), lon_t.max())
    ax.set_ylim(lat_t.min(), lat_t.max())

    text = f"Var. expl.: {results[name]['var_expl']:.1f}%"
    if name != "ERSSTv5":
        text += f"\nPattern corr.: {results[name]['pattern_corr']:.3f}"
    ax.text(
        0.02, 0.98, text,
        transform=ax.transAxes, va="top", ha="left", fontsize=9,
        bbox=dict(facecolor="white", alpha=0.75, edgecolor="none")
    )

cbar = fig.colorbar(
    im, ax=axes, orientation="horizontal",
    fraction=0.05, pad=0.06
)
cbar.set_label("°C per standard deviation")
fig.suptitle("PDO regression patterns, 1950-2014", y=0.98)

# ---------------------------------------------------------------------
# Figure 2: standardized indices, smoothed ERSSTv5, transitions
# ---------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(14, 4.5))

colors = {
    "ERSSTv5": "C0",
    "CESM2": "C1",
    "GFDL-ESM4": "C2",
    "MPI-ESM1-2-LR": "C3",
}

for name in source_order:
    ax.plot(
        time_common, results[name]["pc_std"],
        label=name, color=colors[name], lw=1.2
    )

ax.plot(
    time_common, pdo_pub_std,
    label="Published PDO", color="black", ls="--", lw=1.5
)

ax.plot(
    time_common, smoothed_ersst,
    label="ERSSTv5 9-yr running mean", color="red", lw=2.5
)

for yr in transitions:
    x = pd.Timestamp(f"{yr}-01-01")
    ax.axvline(x, color="0.5", ls=":", lw=1)
    ax.text(
        x, 1.02, str(yr),
        transform=ax.get_xaxis_transform(),
        ha="center", va="bottom", fontsize=8, color="0.35"
    )

ax.axhline(0, color="0.7", lw=0.8)
ax.set_ylabel("Standardized index")
ax.set_xlabel("Year")
ax.set_title("PDO indices, 1950-2014")
ax.legend(ncol=3, fontsize=8, loc="lower left")
