# Figure 1: shared-scale PDO regression patterns.
sources_order = ["ERSSTv5", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

pooled = []
for name in sources_order:
    pooled.append(np.asarray(sources[name]["regpat"].values).ravel())
pooled = np.concatenate(pooled)
vmax = float(np.nanpercentile(np.abs(pooled), 98))
if not np.isfinite(vmax) or vmax <= 0:
    vmax = 1.0

fig1, axes1 = plt.subplots(
    2, 2, figsize=(12, 8), constrained_layout=True
)
axes1 = axes1.ravel()

im = None
for ax, name in zip(axes1, sources_order):
    pat = sources[name]["regpat"]
    p = (
        pat.assign_coords(lon=(pat.lon % 360))
        .sortby("lon")
        .sortby("lat")
    )
    lon = p.lon.values
    lat = p.lat.values

    im = ax.pcolormesh(
        lon,
        lat,
        p.values,
        shading="auto",
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
    )

    ax.set_xlim(110, 260)
    ax.set_ylim(20, 70)
    ax.set_xticks([110, 140, 170, 200, 230, 260])
    ax.set_yticks([20, 30, 40, 50, 60, 70])
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")

    ev = sources[name]["var_exp_pct"]
    if name == "ERSSTv5":
        title = f"ERSSTv5\nVar expl. {ev:.1f}%"
    else:
        title = (
            f"{name}\n"
            f"Var expl. {ev:.1f}%, r = {sources[name]['pattern_r']:.3f}"
        )
    ax.set_title(title, fontsize=10)

fig1.suptitle(
    "PDO SST regression patterns on the common 1° grid "
    "(shared scale; clipped at 98th percentile of |pattern|)",
    fontsize=12,
)
cbar1 = fig1.colorbar(im, ax=axes1, shrink=0.9, extend="both")
cbar1.set_label("SST regression (°C per SD)")

# Figure 2: standardized time series, smoothed ERSSTv5, and transitions.
fig2, ax2 = plt.subplots(
    figsize=(13, 5.5), constrained_layout=True
)

colors = {
    "ERSSTv5": "C0",
    "CESM2": "C1",
    "GFDL-ESM4": "C2",
    "MPI-ESM1-2-LR": "C3",
}

for name in sources_order:
    da = sources[name]["index"]
    ax2.plot(
        da.time.values,
        da.values,
        lw=1.0,
        color=colors[name],
        label=name,
        alpha=0.8,
    )

ax2.plot(
    pdo_pub_plot.time.values,
    pdo_pub_plot.values,
    color="k",
    lw=1.8,
    ls="--",
    label="Published PDO",
)
ax2.plot(
    ersst_smooth_c.time.values,
    ersst_smooth_c.values,
    color="k",
    lw=2.5,
    label="ERSSTv5 9-yr mean",
)

for i, year in enumerate(transition_years):
    t = pd.Timestamp(f"{year}-07-01")
    ax2.axvline(
        t,
        color="0.3",
        ls=":",
        lw=0.9,
        label="9-yr sign change" if i == 0 else "_nolegend_",
    )
    ax2.text(
        t,
        1.02,
        str(year),
        transform=ax2.get_xaxis_transform(),
        rotation=90,
        va="bottom",
        ha="center",
        fontsize=7,
        color="0.3",
    )

ax2.axhline(0, color="grey", lw=0.5)
ax2.set_xlim(pd.Timestamp("1950-01-01"), pd.Timestamp("2015-01-01"))
ax2.set_xlabel("Year")
ax2.set_ylabel("Standardized index")
ax2.set_title(
    "Standardized PDO indices, 1950–2014, with ERSSTv5 "
    "9-year running mean and sign-change years"
)
ax2.legend(
    loc="upper center",
    bbox_to_anchor=(0.5, -0.18),
    ncol=4,
    frameon=False,
    fontsize=8,
)
