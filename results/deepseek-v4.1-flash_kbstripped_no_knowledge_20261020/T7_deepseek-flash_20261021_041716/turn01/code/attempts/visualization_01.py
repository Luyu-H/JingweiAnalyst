# ------------------------------------------------------------------
# Figure 1 - NO3 vs PO4 scatter (200-1000 m), one panel per source
# ------------------------------------------------------------------
# common axis limits, taken from the pooled 99.5th percentiles
xhi = max(float(np.nanpercentile(stats[n]["po4"], 99.5)) for n in order)
yhi = max(float(np.nanpercentile(stats[n]["no3"], 99.5)) for n in order)
xhi = np.ceil(xhi * 2) / 2.0
yhi = np.ceil(yhi / 5.0) * 5.0

fig, axes = plt.subplots(2, 2, figsize=(11, 9.5), constrained_layout=True)
for ax, name in zip(axes.ravel(), order):
    r = results[name]
    ax.scatter(stats[name]["po4"], stats[name]["no3"], s=4, alpha=0.25,
               color="steelblue", linewidths=0, rasterized=True, label="cells (subsample)")

    xs = np.linspace(0.0, xhi, 100)
    ax.plot(xs, r["slope"] * xs + r["intercept"], color="crimson", lw=2.0,
            label=f"fit  NO3 = {r['slope']:.2f}*PO4 + {r['intercept']:.2f}")
    ax.plot(xs, 16.0 * xs, color="black", lw=1.6, ls="--", label="Redfield 16:1")

    ax.set_xlim(0, xhi)
    ax.set_ylim(0, yhi)
    ax.set_xticks(np.linspace(0, xhi, 6))
    ax.set_yticks(np.linspace(0, yhi, 6))
    ax.set_xlabel("Phosphate (umol kg$^{-1}$)")
    ax.set_ylabel("Nitrate (umol kg$^{-1}$)")
    ax.set_title(name)
    ax.text(0.97, 0.05,
            f"a = {r['slope']:.2f}\nb = {r['intercept']:.2f} umol kg$^{{-1}}$\n$R^2$ = {r['r2']:.3f}\n"
            f"mean N:P = {r['mean_np_ratio']:.2f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.6", alpha=0.85))
    ax.legend(loc="upper left", fontsize=7, framealpha=0.85)

fig.suptitle("Nitrate vs phosphate, 200-1000 m: fitted regression and Redfield 16:1", fontsize=13)

# ------------------------------------------------------------------
# Figure 2 - 200-1000 m mean N:P minus 16, common diverging scale
# ------------------------------------------------------------------
pooled = np.concatenate([maps[n].values.ravel()[np.isfinite(maps[n].values.ravel())]
                         for n in order])
vmax = float(np.nanpercentile(np.abs(pooled), 98))
vmax = max(vmax, 0.5)

fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), constrained_layout=True)
pcm = None
for ax, name in zip(axes.ravel(), order):
    f = maps[name]
    lat = f.coords[f.attrs.get("lat_coord", "lat")].values
    lon = f.coords[f.attrs.get("lon_coord", "lon")].values
    pcm = ax.pcolormesh(lon, lat, f.values, cmap="RdBu_r",
                        vmin=-vmax, vmax=vmax, shading="auto")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-70, 85)
    ax.set_xticks(np.arange(-180, 181, 60))
    ax.set_yticks(np.arange(-60, 90, 30))
    ax.set_title(name)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

fig.colorbar(pcm, ax=axes, label="N:P - 16 (unitless)", shrink=0.9, extend="both")
fig.suptitle(f"200-1000 m mean N:P minus 16 (shared scale, clipped at 98th pct = +/-{vmax:.1f})",
             fontsize=12)
