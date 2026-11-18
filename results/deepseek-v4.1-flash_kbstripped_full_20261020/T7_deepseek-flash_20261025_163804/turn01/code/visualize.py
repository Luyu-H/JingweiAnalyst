# ------------------------------------------------------------------
# Figure 1: NO3 vs PO4 scatter over 200-1000 m, common-mask sample
# ------------------------------------------------------------------
fit_redfield = lambda x: 16.0 * x

pooled_x = np.concatenate([scatter_data[s][0] for s in src_names])
pooled_y = np.concatenate([scatter_data[s][1] for s in src_names])
xmax = float(np.nanpercentile(pooled_x, 99.5))
ymax = float(np.nanpercentile(pooled_y, 99.5))

fig, axes = plt.subplots(2, 2, figsize=(11, 9), constrained_layout=True)
axes = axes.ravel()
xline = np.linspace(0, xmax, 200)

for ax, s in zip(axes, src_names):
    pv, nv = scatter_data[s]
    ax.scatter(pv, nv, s=1, alpha=0.15, color="#3b6ea5", rasterized=True,
               label="cells")
    O = ols[s]
    ax.plot(xline, O["a"] * xline + O["b"], color="#c0392b", lw=1.8,
            label=f"OLS fit (a={O['a']:.2f}, b={O['b']:.2f})")
    ax.plot(xline, fit_redfield(xline), color="k", lw=1.3, ls="--",
            label="Redfield 16:1")

    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    ax.set_xlabel("Phosphate (µmol kg$^{-1}$)")
    ax.set_ylabel("Nitrate (µmol kg$^{-1}$)")
    ax.set_title(f"{s}   (n = {O['n']} cells)")
    ax.grid(alpha=0.25)

    txt = (f"NO$_3$ = {O['a']:.2f}·PO$_4$ {O['b']:+.2f}\n"
           f"$R^2$ = {O['r2']:.3f}")
    ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", ha="left",
            fontsize=9, bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.85))
    ax.legend(loc="lower right", fontsize=7, framealpha=0.85)

fig.suptitle("Figure 1 — Nitrate vs phosphate, 200-1000 m "
             "(common WOA23-grid sample; equal weight per cell)")
plt.show()


# ------------------------------------------------------------------
# Figure 2: 200-1000 m vertical volume-weighted mean N:P minus 16
# ------------------------------------------------------------------
vals = np.concatenate([np.abs(np.asarray(np_map[s].values).ravel()) for s in src_names])
vals = vals[np.isfinite(vals)]
vmax = float(np.nanpercentile(vals, 98))

fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
axes = axes.ravel()
pc = None
for ax, s in zip(axes, src_names):
    f = np_map[s]
    lon = np.asarray(f["lon"].values)
    lat = np.asarray(f["lat"].values)
    pc = ax.pcolormesh(lon, lat, np.asarray(f.values), cmap="RdBu_r",
                       vmin=-vmax, vmax=vmax, shading="auto")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-80, 90)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(f"{s}")
    ax.grid(alpha=0.2)

cbar = fig.colorbar(pc, ax=axes, extend="both", shrink=0.85)
cbar.set_label("N:P − 16")

fig.suptitle("Figure 2 — 200-1000 m volume-weighted mean N:P minus 16 "
             f"(models bilinearly regridded for display; shared scale "
             f"±{vmax:.2f}, 98th pct)")
plt.show()
