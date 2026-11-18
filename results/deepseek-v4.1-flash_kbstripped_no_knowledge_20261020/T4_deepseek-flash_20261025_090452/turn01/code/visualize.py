import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
months = [3, 9]
MONTH_NAME = {3: "March", 9: "September"}


def _coord(da, names, fallback):
    for n in names:
        if n in da.coords:
            return da.coords[n].values
    return da.coords[fallback].values


# ------------------------------------------------------------------ Figure 1
pooled = np.concatenate([np.abs(np.asarray(diff_map[k].values).ravel())
                         for k in diff_map])
vmax = float(np.nanpercentile(pooled, 98))
if not np.isfinite(vmax) or vmax <= 0:
    vmax = 1.0
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

fig, axes = plt.subplots(3, 2, figsize=(15, 10), constrained_layout=True)
for i, model in enumerate(models):
    for j, mm in enumerate(months):
        ax = axes[i, j]
        da = diff_map[(model, mm)]
        lon = _coord(da, ("lon", "longitude"), "lon")
        lat = _coord(da, ("lat", "latitude"), "lat")
        z = np.clip(np.asarray(da.values), -vmax, vmax)
        pm = ax.pcolormesh(lon, lat, z, cmap="RdBu_r", norm=norm,
                           shading="auto", rasterized=True)
        ax.set_xlim(-180, 180)
        ax.set_ylim(-80, 90)
        ax.set_title(f"{model} — {MONTH_NAME[mm]}", fontsize=11)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        r = results[(model, mm)]
        ax.text(0.02, 0.03,
                f"bias {r['bias']:.1f} m   RMSD {r['rmse']:.1f} m   r {r['correlation']:.3f}",
                transform=ax.transAxes, fontsize=8, va="bottom",
                bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.8))

cb = fig.colorbar(pm, ax=axes, orientation="vertical", shrink=0.9, extend="both")
cb.set_label("MLD difference (m)")
fig.suptitle("Operator MLD minus published mlotst, 1971–2000 March/September climatology\n"
             f"(display grid 1°; colour scale shared, clipped at the pooled 98th percentile "
             f"of |difference| = {vmax:.1f} m)", fontsize=12)
plt.show()


# ------------------------------------------------------------------ Figure 2
pairs = {}
allvals = []
for model in models:
    for mm in months:
        a = np.asarray(op_clim[(model, mm)].values).ravel()
        b = np.asarray(ml_clim[(model, mm)].values).ravel()
        ok = np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
        pairs[(model, mm)] = (a[ok], b[ok])
        if ok.any():
            allvals.append(np.concatenate([a[ok], b[ok]]))
allvals = np.concatenate(allvals)
lo = max(0.5, float(np.nanpercentile(allvals, 0.2)))
hi = float(np.nanmax(allvals)) * 1.1
lo = 10 ** np.floor(np.log10(lo))
hi = 10 ** np.ceil(np.log10(hi))
tick_lo = int(np.floor(np.log10(lo)))
tick_hi = int(np.ceil(np.log10(hi)))
ticks = [10.0 ** e for e in range(tick_lo, tick_hi + 1)]
ticks = [t for t in ticks if lo <= t <= hi]

fig2, axes2 = plt.subplots(3, 2, figsize=(13, 12), constrained_layout=True)
for i, model in enumerate(models):
    for j, mm in enumerate(months):
        ax = axes2[i, j]
        a, b = pairs[(model, mm)]
        ax.scatter(a, b, s=1, alpha=0.25, color="#1f4e79", edgecolors="none",
                   rasterized=True)
        ax.plot([lo, hi], [lo, hi], "k--", lw=1.0, label="y = x")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels([f"{t:g}" for t in ticks])
        ax.set_yticklabels([f"{t:g}" for t in ticks])
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("Operator MLD (m)")
        ax.set_ylabel("Published mlotst (m)")
        ax.set_title(f"{model} — {MONTH_NAME[mm]}", fontsize=11)
        r = results[(model, mm)]
        ax.text(0.04, 0.96,
                f"bias {r['bias']:.1f} m\nRMSD {r['rmse']:.1f} m\nr = {r['correlation']:.3f}",
                transform=ax.transAxes, fontsize=9, va="top", ha="left",
                bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
        ax.grid(True, which="both", lw=0.3, alpha=0.4)
        if i == 0 and j == 0:
            ax.legend(loc="lower right", fontsize=8)

fig2.suptitle("Operator MLD vs published mlotst, 1971–2000 March/September climatology\n"
              "(one point per common-mask ocean cell, identical logarithmic axes)",
              fontsize=12)
plt.show()
