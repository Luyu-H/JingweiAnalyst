# Figure 1: PDO regression patterns on the common 1-degree grid.
fig1, axes1 = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
panel_order = ["ERSSTv5", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

pooled = np.concatenate([pdo_results[s]["pattern"].ravel() for s in panel_order])
vmax = float(np.nanpercentile(np.abs(pooled), 98))

im = None
for ax, s in zip(axes1.ravel(), panel_order):
    res = pdo_results[s]
    im = ax.pcolormesh(
        res["lon"], res["lat"], res["pattern"],
        cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto"
    )
    ax.set_title(s)
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")
    ax.set_xlim(110, 260)
    ax.set_ylim(20, 70)

    txt = f"Var exp: {res['var_explained']:.1f}%"
    if s in pattern_corrs:
        txt += f"\nPattern corr: {pattern_corrs[s]['corr']:.3f} (n={pattern_corrs[s]['n_shared']})"
    ax.text(
        0.02, 0.95, txt,
        transform=ax.transAxes, va="top", ha="left",
        fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8, edgecolor="none")
    )

fig1.colorbar(
    im, ax=axes1.ravel().tolist(), orientation="vertical", extend="both",
    label="SST regression (°C per std)", shrink=0.9
)
fig1.suptitle(
    "PDO regression patterns on the 1° bilinear grid, 1950–2014 "
    "(shared scale; clipped at 98th percentile)"
)

# Figure 2: standardized indices, published PDO, smoothed ERSSTv5, transitions.
fig2, ax2 = plt.subplots(figsize=(12, 6), constrained_layout=True)

for s in panel_order:
    ax2.plot(source_times[s], pdo_results[s]["pc"], label=s, lw=1.2)

ax2.plot(
    published_pdo_on_ersst_time.index,
    published_pdo_on_ersst_time.values,
    label="Published PDO", color="black", ls="--", lw=1.2
)
ax2.plot(
    smoothed_ersst.index,
    smoothed_ersst.values,
    label="ERSSTv5 9-yr running mean", color="black", lw=3.0
)

for yr in transition_years:
    ax2.axvline(pd.Timestamp(f"{yr}-01-01"), color="gray", ls=":", lw=1.0)
    ax2.text(
        pd.Timestamp(f"{yr}-01-01"), 0.95, str(yr),
        transform=ax2.get_xaxis_transform(),
        rotation=90, ha="right", va="top", fontsize=8, color="gray"
    )

ax2.axhline(0, color="gray", lw=0.8)
ax2.set_xlabel("Year")
ax2.set_ylabel("Standardized index")
ax2.set_title(
    "Standardized PDO indices and ERSSTv5 9-year smoothed curve, 1950–2014"
)
ax2.set_xlim(pd.Timestamp("1950-01-01"), pd.Timestamp("2014-12-31"))
ax2.legend(
    loc="upper left", bbox_to_anchor=(1.01, 1.0),
    borderaxespad=0, fontsize=9
)
