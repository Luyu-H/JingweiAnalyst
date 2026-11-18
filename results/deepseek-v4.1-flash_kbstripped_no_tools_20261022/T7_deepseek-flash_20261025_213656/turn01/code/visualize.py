# -----------------------------------------------------------------------------
# Figure 1: 4-panel NO3 vs PO4 scatter with OLS and Redfield lines
# -----------------------------------------------------------------------------
fig1, axes1 = plt.subplots(2, 2, figsize=(11, 9), constrained_layout=True)
axes1 = axes1.ravel()
rng = np.random.default_rng(42)

for ax, name in zip(axes1, SOURCE_ORDER):
    r = source_results[name]
    x = r["po4_vals"]
    y = r["no3_vals"]

    n_plot = 20000
    if x.size > n_plot:
        idx = rng.choice(x.size, size=n_plot, replace=False)
        xs, ys = x[idx], y[idx]
    else:
        xs, ys = x, y

    ax.scatter(
        xs,
        ys,
        s=3,
        alpha=0.25,
        color="tab:blue",
        edgecolors="none",
        rasterized=True,
    )

    xfit = np.linspace(0.0, fig1_xlim[1], 100)
    ax.plot(
        xfit,
        r["slope_a"] * xfit + r["intercept_b"],
        color="tab:red",
        lw=1.5,
        label="OLS fit",
    )
    ax.plot(
        xfit,
        16.0 * xfit,
        color="black",
        ls="--",
        lw=1.2,
        label="Redfield 16:1",
    )

    ax.set_xlim(fig1_xlim)
    ax.set_ylim(fig1_ylim)
    ax.set_xticks(np.linspace(0.0, fig1_xlim[1], 5))
    ax.set_yticks(np.linspace(0.0, fig1_ylim[1], 5))

    ax.set_xlabel("Phosphate (µmol kg$^{-1}$)")
    ax.set_ylabel("Nitrate (µmol kg$^{-1}$)")
    ax.set_title(name)
    ax.grid(True, alpha=0.3)

    txt = (
        f"a={r['slope_a']:.2f}\n"
        f"b={r['intercept_b']:.2f}\n"
        f"R²={r['r2']:.3f}"
    )
    ax.text(
        0.04,
        0.96,
        txt,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="none"),
    )
    ax.legend(loc="lower right", fontsize=8, framealpha=0.85)

fig1.suptitle(
    "200–1000 m nitrate vs phosphate: WOA23 and CMIP6 historical models",
    fontsize=13,
)

# -----------------------------------------------------------------------------
# Figure 2: 2x2 global maps of 200–1000 m mean N:P minus 16
# -----------------------------------------------------------------------------
fig2, axes2 = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
axes2 = axes2.ravel()
cmap = plt.get_cmap("RdBu_r")
mesh = None

for ax, name in zip(axes2, SOURCE_ORDER):
    data = source_results[name]["map_common_anom"]
    mesh = ax.pcolormesh(
        common_lon,
        common_lat,
        data,
        cmap=cmap,
        vmin=-vmax_map,
        vmax=vmax_map,
        shading="auto",
    )
    ax.set_xlim(0.0, 360.0)
    ax.set_ylim(-90.0, 90.0)
    ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(name)
    ax.grid(True, alpha=0.3)

cbar = fig2.colorbar(
    mesh,
    ax=axes2,
    orientation="horizontal",
    fraction=0.05,
    pad=0.05,
    extend="both",
)
cbar.set_label("N:P − 16")

fig2.suptitle(
    "200–1000 m mean N:P − 16; conservative area-weighted regrid to 1° for display "
    "(shared scale, clipped at 98th percentile of |anomaly|)",
    fontsize=11,
)
