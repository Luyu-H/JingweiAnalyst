# Figure 1: four-panel NO3 vs PO4 scatter with fitted and Redfield lines
fig1, axes1 = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
axes1 = axes1.ravel()

for ax, name in zip(axes1, scatter_order):
    x = scatter_samples[name]["x"]
    y = scatter_samples[name]["y"]

    # Thin very large scatters for legibility while preserving full regression stats already computed
    n = x.size
    if n > 80000:
        idx = np.linspace(0, n - 1, 80000).astype(int)
        xs = x[idx]
        ys = y[idx]
    else:
        xs = x
        ys = y

    ax.scatter(xs, ys, s=2, alpha=0.18, color="tab:blue", edgecolors="none", rasterized=False)

    xx = np.linspace(xlim[0], xlim[1], 200)
    reg = regression_stats[name]
    ax.plot(xx, reg["slope"] * xx + reg["intercept"], color="tab:red", lw=2, label="OLS fit")
    ax.plot(xx, 16.0 * xx, color="k", lw=1.5, ls="--", label="Redfield 16:1")

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_title(name)
    ax.set_xlabel("Phosphate (µmol kg$^{-1}$)")
    ax.set_ylabel("Nitrate (µmol kg$^{-1}$)")
    ax.text(
        0.03, 0.97,
        f"a = {reg['slope']:.2f}\n"
        f"b = {reg['intercept']:.2f} µmol kg$^{{-1}}$\n"
        f"$R^2$ = {reg['r2']:.3f}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=10,
        bbox=dict(facecolor="white", alpha=0.8, edgecolor="0.7")
    )
    ax.legend(loc="best", frameon=True)

fig1.suptitle("Figure 1. Intermediate-water nitrate vs phosphate, 200-1000 m", fontsize=14)

# Figure 2: 2x2 global maps of mean N:P minus 16 on a shared scale
fig2, axes2 = plt.subplots(2, 2, figsize=(14, 8), constrained_layout=True)
axes2 = axes2.ravel()

mappable = None
for ax, name in zip(axes2, scatter_order):
    da = map_fields_display[name]
    lon = da["lon"].values
    lat = da["lat"].values
    mappable = ax.pcolormesh(
        lon, lat, da.values,
        cmap="RdBu_r",
        vmin=-map_vmax,
        vmax=map_vmax,
        shading="auto",
    )
    ax.set_title(name)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(lon.min(), lon.max())
    ax.set_ylim(lat.min(), lat.max())

cbar = fig2.colorbar(mappable, ax=axes2, orientation="vertical", shrink=0.9, extend="both")
cbar.set_label("Mean N:P - 16 (1)")
fig2.suptitle("Figure 2. 200-1000 m mean N:P minus 16 on shared 98th-percentile-clipped scale", fontsize=14)

plt.show()
