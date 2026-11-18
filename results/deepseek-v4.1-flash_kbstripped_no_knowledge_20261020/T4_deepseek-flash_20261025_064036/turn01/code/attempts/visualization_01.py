import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

mon_list = ["March", "September"]

# ---------- Figure 1: difference maps ----------
pooled = []
for da in fig1_maps.values():
    v = da.values.astype(float).ravel()
    v = v[np.isfinite(v)]
    if v.size:
        pooled.append(v)
pooled = np.concatenate(pooled)
vmax = float(np.nanpercentile(np.abs(pooled), 98))
if not np.isfinite(vmax) or vmax <= 0:
    vmax = 100.0

fig1, axes1 = plt.subplots(len(models), 2, figsize=(13, 10.5), constrained_layout=True)

for i, m in enumerate(models):
    for j, mon_name in enumerate(mon_list):
        ax = axes1[i, j]
        da = fig1_maps[(m, mon_name)]
        vals = da.values
        mesh = ax.pcolormesh(da.lon.values, da.lat.values, vals,
                             cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                             shading="auto", extend="both")
        ax.set_title(f"{m} – {mon_name}", fontsize=10)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_xlim(-180, 180)
        ax.set_ylim(-80, 90)

cb = fig1.colorbar(mesh, ax=axes1, orientation="horizontal",
                   fraction=0.04, pad=0.04, extend="both")
cb.set_label("Operator MLD − mlotst (m)")

fig1.suptitle(
    f"Figure 1. Offline operator MLD minus published mlotst, 1971–2000 climatologies.\n"
    f"Shared diverging scale clipped at the 98th percentile of |difference| = {vmax:.1f} m. "
    f"Each panel regridded to a 1° longitude–latitude frame for display only.",
    fontsize=11)

# ---------- Figure 2: scatter operator MLD vs mlotst ----------
allx, ally = [], []
for m in models:
    for mon_name in mon_list:
        x, y = scatter[(m, mon_name)]
        good = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
        allx.append(x[good])
        ally.append(y[good])
allx = np.concatenate(allx)
ally = np.concatenate(ally)
lo = max(1.0, float(np.nanpercentile(np.concatenate([allx, ally]), 0.5)))
hi = float(np.nanpercentile(np.concatenate([allx, ally]), 99.5))
lo = 10 ** np.floor(np.log10(lo))
hi = 10 ** np.ceil(np.log10(hi))
ticks = np.logspace(np.log10(lo), np.log10(hi), 5)

fig2, axes2 = plt.subplots(len(models), 2, figsize=(12, 11), constrained_layout=True)

for i, m in enumerate(models):
    for j, mon_name in enumerate(mon_list):
        ax = axes2[i, j]
        x, y = scatter[(m, mon_name)]
        good = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
        ax.scatter(x[good], y[good], s=3, alpha=0.25, edgecolors="none",
                   color="steelblue", rasterized=True)
        ax.plot([lo, hi], [lo, hi], "k--", lw=1.0, label="y = x")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels([f"{t:g}" for t in ticks], fontsize=8, rotation=30)
        ax.set_yticklabels([f"{t:g}" for t in ticks], fontsize=8)
        ax.set_xlabel("Operator MLD (m)")
        ax.set_ylabel("mlotst (m)")
        ax.set_title(f"{m} – {mon_name}", fontsize=10)
        ax.minorticks_off()
        r = results[(m, mon_name)]
        ax.text(0.03, 0.97,
                f"bias = {r['bias']:.1f} m\nRMSE = {r['rmse']:.1f} m\nr = {r['corr']:.3f}",
                transform=ax.transAxes, va="top", ha="left", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                          edgecolor="0.7", alpha=0.85))
        ax.grid(True, which="major", ls=":", lw=0.4, alpha=0.5)

fig2.suptitle(
    "Figure 2. Operator MLD (x) against published mlotst (y), one point per common-valid ocean cell.\n"
    "Identical logarithmic axes in every panel; dashed line is y = x.",
    fontsize=11)
