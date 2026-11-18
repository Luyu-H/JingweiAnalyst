# -----------------------------------------------------------------------------
# Figure 1: nitrate versus phosphate, 200-1000 m, all valid cells equally weighted
# -----------------------------------------------------------------------------
all_x = np.concatenate([scatter_data[k]["po4"] for k in ORDER])
all_y = np.concatenate([scatter_data[k]["no3"] for k in ORDER])

xmax = np.nanpercentile(all_x, 99.5)
ymax = np.nanpercentile(all_y, 99.5)
xmax = max(float(xmax), 1.0) * 1.05
ymax = max(float(ymax), 1.0) * 1.05
xlim = (0.0, xmax)
ylim = (0.0, ymax)

fig1, axes1 = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
for ax, name in zip(axes1.ravel(), ORDER):
    d = scatter_data[name]
    res = results[name]

    ax.scatter(
        d["po4"], d["no3"],
        s=1.5, alpha=0.25, color="steelblue",
        rasterized=True, edgecolors="none"
    )

    xfit = np.linspace(xlim[0], xlim[1], 300)
    ax.plot(
        xfit, res["a"] * xfit + res["b"],
        color="crimson", lw=1.8, label="OLS fit"
    )
    ax.plot(
        xfit, 16.0 * xfit,
        color="black", lw=1.5, ls="--", label="Redfield 16:1"
    )

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel(r"Phosphate ($\mu$mol kg$^{-1}$)")
    ax.set_ylabel(r"Nitrate ($\mu$mol kg$^{-1}$)")
    ax.set_title(
        f"{name}\n"
        f"a={res['a']:.2f}, b={res['b']:.2f} $\\mu$mol kg$^{{-1}}$, "
        f"$R^2$={res['r2']:.3f}"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8)

fig1.suptitle(
    "Nitrate vs phosphate, 200-1000 m; all valid cells weighted equally"
)

# -----------------------------------------------------------------------------
# Figure 2: 200-1000 m mean N:P minus 16, display-only common regular grid
# -----------------------------------------------------------------------------
pooled = np.concatenate([map_fields[k]["field"].ravel() for k in ORDER])
vmax = float(np.nanpercentile(np.abs(pooled), 98.0))
if not np.isfinite(vmax) or vmax <= 0.0:
    vmax = 1.0

fig2, axes2 = plt.subplots(2, 2, figsize=(12, 7), constrained_layout=True)
pm = None
for ax, name in zip(axes2.ravel(), ORDER):
    mf = map_fields[name]
    pm = ax.pcolormesh(
        mf["lon"], mf["lat"], mf["field"],
        cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto"
    )
    ax.set_xlim(-180.0, 180.0)
    ax.set_ylim(-90.0, 90.0)
    ax.set_xticks(np.arange(-180, 181, 60))
    ax.set_yticks(np.arange(-90, 91, 30))
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(name)
    ax.set_facecolor("0.85")
    ax.grid(True, alpha=0.3)

fig2.colorbar(
    pm,
    ax=axes2,
    orientation="horizontal",
    fraction=0.05,
    pad=0.05,
    extend="both",
    label="N:P - 16 (unitless)",
)
fig2.suptitle(
    "200-1000 m mean N:P minus 16; display regridded, shared scale "
    "(clipped at pooled 98th percentile)"
)
