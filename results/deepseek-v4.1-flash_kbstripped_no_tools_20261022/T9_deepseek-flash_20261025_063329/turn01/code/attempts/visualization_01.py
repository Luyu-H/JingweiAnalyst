# ============================================================
# Figure 1: 2x2 regression patterns on the common 1-degree grid
# ============================================================
source_order = ["ERSSTv5", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

pooled = np.concatenate([patterns[name][shared_mask] for name in source_order])
vmax = float(np.nanpercentile(np.abs(pooled), 98))
if not np.isfinite(vmax) or vmax == 0.0:
    vmax = 1.0

fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)

for ax, name in zip(axes.ravel(), source_order):
    im = ax.pcolormesh(
        target_lon, target_lat, patterns[name],
        shading="auto", cmap="RdBu_r",
        vmin=-vmax, vmax=vmax
    )
    ax.set_xlim(110, 260)
    ax.set_ylim(20, 70)
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")
    ax.grid(True, color="gray", linewidth=0.3, alpha=0.5)

    ve = metrics[name]["variance_explained_pct"]
    if name == "ERSSTv5":
        title = f"ERSSTv5\nVar = {ve:.1f}%"
    else:
        title = f"{name}\nVar = {ve:.1f}%, r = {pattern_corr[name]:.3f}"
    ax.set_title(title, fontsize=10)

fig.colorbar(
    im, ax=axes.ravel().tolist(), orientation="horizontal",
    fraction=0.04, pad=0.06, extend="both",
    label="SST regression (°C per SD)"
)
fig.suptitle(
    "PDO regression patterns, 1950–2014 "
    "(colour clipped at pooled 98th percentile)",
    fontsize=12
)

# ============================================================
# Figure 2: standardized indices, smoothed ERSSTv5 index,
# and observed phase-transition years.
# ============================================================
fig2, ax2 = plt.subplots(figsize=(13, 5), constrained_layout=True)

for name in source_order:
    ax2.plot(common_time, indices[name], lw=0.8, alpha=0.65, label=name)

ax2.plot(
    pdo_std_series.index, pdo_std_series.values,
    lw=1.2, ls="--", color="black", alpha=0.9,
    label="Published PDO (std)"
)
ax2.plot(
    common_time, ersst_smoothed.values,
    lw=2.8, color="black",
    label="ERSSTv5 9-yr running mean"
)

ymin, ymax = ax2.get_ylim()
for year in pdo_transition_years:
    t = pd.Timestamp(f"{year}-01-01")
    ax2.axvline(t, color="grey", ls=":", lw=1.0, alpha=0.8)
    ax2.text(
        t, ymax * 0.95, str(year),
        rotation=90, va="top", ha="right",
        fontsize=8, color="grey"
    )

ax2.set_xlabel("Year")
ax2.set_ylabel("Standardized index")
ax2.set_title("PDO indices and ERSSTv5 phase transitions (1950–2014)")
ax2.grid(True, alpha=0.3)
ax2.legend(
    loc="upper left", bbox_to_anchor=(1.01, 1.0),
    borderaxespad=0.0, fontsize=8
)

plt.show()
