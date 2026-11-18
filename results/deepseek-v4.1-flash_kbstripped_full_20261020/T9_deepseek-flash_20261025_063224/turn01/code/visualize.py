# Visualization

# Figure 1: PDO regression patterns
def _to_plot_domain(da):
    p = da.copy()
    p = p.assign_coords(lon=(p.lon % 360)).sortby("lon")
    p = p.sortby("lat")
    p = p.sel(lat=slice(20, 70), lon=slice(110, 260))
    return p


panel_names = ["ERSSTv5", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
plot_patterns = [
    _to_plot_domain(results[name]["reg_pattern"]) for name in panel_names
]

pooled = np.concatenate([p.values.ravel() for p in plot_patterns])
vmax = max(float(np.nanpercentile(np.abs(pooled), 98)), 1e-6)

fig, axes = plt.subplots(
    2, 2, figsize=(12, 9), constrained_layout=True
)

for ax, name, p in zip(axes.ravel(), panel_names, plot_patterns):
    im = ax.pcolormesh(
        p.lon.values,
        p.lat.values,
        p.values,
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        shading="auto",
    )
    ax.set_xlim(110, 260)
    ax.set_ylim(20, 70)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.grid(True, color="0.8", lw=0.5)

    title = (
        f"{name}\n"
        f"Explained variance: {results[name]['explained_variance']:.1f}%"
    )
    if name != "ERSSTv5":
        title += f"\nPattern corr.: {pattern_corr[name]:.3f}"
    ax.set_title(title, fontsize=10)

cbar = fig.colorbar(
    im,
    ax=axes,
    orientation="horizontal",
    shrink=0.7,
    pad=0.06,
    extend="both",
    label="SST regression (degC per std)",
)
fig.suptitle(
    "PDO regression patterns on the 1-degree grid; "
    "shared scale clipped at the 98th percentile of |values|",
    fontsize=11,
)

# Figure 2: standardized indices and ERSSTv5 phase transitions
fig2, ax2 = plt.subplots(figsize=(14, 5), constrained_layout=True)

for name in panel_names:
    idx = results[name]["index"].squeeze(drop=True)
    ax2.plot(
        idx.time.values,
        idx.values,
        lw=1.0,
        alpha=0.8,
        label=name,
    )

ax2.plot(
    pub_pdo.time.values,
    pub_pdo.values,
    color="black",
    lw=1.3,
    ls="--",
    label="Published PDO (as published)",
)

ax2.plot(
    ersst_smoothed.time.values,
    ersst_smoothed.values,
    color="navy",
    lw=2.5,
    label="ERSSTv5 9-year running mean",
)

ax2.axhline(0.0, color="0.7", lw=0.6)

for y in transition_years:
    x = pd.Timestamp(f"{y}-01-01")
    ax2.axvline(x, color="0.5", ls=":", lw=1.0, zorder=0)
    ax2.text(
        x,
        1.0,
        str(y),
        transform=ax2.get_xaxis_transform(),
        ha="center",
        va="bottom",
        fontsize=7,
        rotation=90,
        color="0.3",
        clip_on=False,
    )

ax2.set_xlim(pd.Timestamp("1950-01-01"), pd.Timestamp("2015-01-01"))
ax2.set_xlabel("Year")
ax2.set_ylabel("Standardized index")
ax2.set_title(
    "PDO indices and ERSSTv5 phase transitions (dotted lines)"
)
ax2.legend(
    loc="upper center",
    bbox_to_anchor=(0.5, -0.15),
    ncol=3,
    fontsize=8,
    frameon=False,
)
