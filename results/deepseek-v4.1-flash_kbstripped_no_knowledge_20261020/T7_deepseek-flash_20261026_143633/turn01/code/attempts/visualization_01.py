# ----------------------------------------------------------------------
# Figure 1: nitrate vs phosphate scatter and OLS regressions
# ----------------------------------------------------------------------
fig1, axes1 = plt.subplots(2, 2, figsize=(10, 10), constrained_layout=True)
axes1 = axes1.ravel()

# Common axis limits from all sources (subsampled points)
all_x = []
all_y = []
for src in sources_order:
    if src in results and "scatter_x" in results[src]:
        all_x.append(results[src]["scatter_x"])
        all_y.append(results[src]["scatter_y"])

if all_x:
    x_all = np.concatenate(all_x)
    y_all = np.concatenate(all_y)
    x_max = np.nanpercentile(x_all, 99.5)
    y_max = np.nanpercentile(y_all, 99.5)
    x_max = max(float(x_max) * 1.05, 1.0)
    y_max = max(float(y_max) * 1.05, 1.0)
else:
    x_max, y_max = 3.0, 45.0

x_min, y_min = 0.0, 0.0
x_ticks = np.linspace(x_min, x_max, 5)
y_ticks = np.linspace(y_min, y_max, 6)

for ax, src in zip(axes1, sources_order):
    if src not in results or "missing" in results[src] or "scatter_x" not in results[src]:
        ax.text(
            0.5, 0.5, f"{src}\nnot available",
            transform=ax.transAxes, ha="center", va="center", fontsize=12,
        )
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_xticks(x_ticks)
        ax.set_yticks(y_ticks)
        ax.set_title(src)
        ax.set_xlabel("Phosphate (µmol kg$^{-1}$)")
        ax.set_ylabel("Nitrate (µmol kg$^{-1}$)")
        continue

    r = results[src]
    ax.scatter(
        r["scatter_x"], r["scatter_y"],
        s=1, alpha=0.25, c="steelblue", edgecolors="none", rasterized=True,
    )

    xx = np.linspace(x_min, x_max, 200)
    yy_ols = r["slope"] * xx + r["intercept"]
    ax.plot(xx, yy_ols, color="crimson", lw=1.8, label="OLS fit")
    ax.plot(xx, 16.0 * xx, color="black", lw=1.2, ls="--", label="Redfield 16:1")

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)
    ax.set_title(src)
    ax.set_xlabel("Phosphate (µmol kg$^{-1}$)")
    ax.set_ylabel("Nitrate (µmol kg$^{-1}$)")
    ax.grid(True, alpha=0.3)

    ax.text(
        0.05, 0.95,
        f"a = {r['slope']:.2f}\n"
        f"b = {r['intercept']:.2f}\n"
        f"R$^2$ = {r['r2']:.3f}",
        transform=ax.transAxes, va="top", ha="left", fontsize=9,
        bbox=dict(boxstyle="round", fc="white", alpha=0.85),
    )
    ax.legend(loc="best", fontsize=7, framealpha=0.8)

fig1.suptitle(
    "Nitrate vs phosphate, 200–1000 m (points subsampled for display)",
    fontsize=14,
)

# ----------------------------------------------------------------------
# Figure 2: maps of 200-1000 m mean N:P minus 16
# ----------------------------------------------------------------------
fig2, axes2 = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
axes2 = axes2.ravel()

map_arrays = []
for src in sources_order:
    if src in results and "map_display" in results[src]:
        map_arrays.append(np.asarray(results[src]["map_display"].values).ravel())

if map_arrays:
    pooled = np.concatenate(map_arrays)
    vmax = float(np.nanpercentile(np.abs(pooled), 98))
    vmax = max(vmax, 0.5)
else:
    vmax = 1.0

cmap = plt.get_cmap("RdBu_r")
pcm = None

for ax, src in zip(axes2, sources_order):
    if src not in results or "map_display" not in results[src] or "missing" in results[src]:
        ax.text(
            0.5, 0.5, f"{src}\nnot available",
            transform=ax.transAxes, ha="center", va="center", fontsize=12,
        )
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        ax.set_xticks(np.arange(-180, 181, 60))
        ax.set_yticks(np.arange(-90, 91, 30))
        ax.set_title(src)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        continue

    da = results[src]["map_display"]
    lon = da.lon.values
    lat = da.lat.values

    pcm = ax.pcolormesh(
        lon, lat, da.values,
        cmap=cmap, vmin=-vmax, vmax=vmax, shading="auto",
    )
    ax.set_title(src)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xticks(np.arange(-180, 181, 60))
    ax.set_yticks(np.arange(-90, 91, 30))

if pcm is not None:
    cbar = fig2.colorbar(
        pcm, ax=axes2, orientation="horizontal",
        fraction=0.05, pad=0.05, extend="both",
    )
    cbar.set_label("N:P − 16 (dimensionless)")

fig2.suptitle(
    f"200–1000 m mean N:P − 16; shared scale clipped at 98th percentile "
    f"(|Δ| ≤ {vmax:.2f})",
    fontsize=14,
)
