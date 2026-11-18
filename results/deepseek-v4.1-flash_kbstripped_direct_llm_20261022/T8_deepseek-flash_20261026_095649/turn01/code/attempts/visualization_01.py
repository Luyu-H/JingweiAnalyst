import numpy as np
import matplotlib.pyplot as plt

# ---------------- Figure 1: 0-50 m DIC bias maps (shared diverging scale) ----------
def _disp(arr):
    """Roll the 20.5-379.5 grid so the map runs from 180W eastwards."""
    shift = int(np.searchsorted(np.mod(GLON, 360.0), 180.0))
    a = np.roll(arr, -shift, axis=1)
    lo = np.roll(GLON, -shift)
    lo = np.where(lo > 180.0, lo - 360.0, lo)
    return lo, a

allb = np.concatenate([BIAS_MAPS[m][np.isfinite(BIAS_MAPS[m])].ravel() for m in MODELS])
vmax = float(np.nanpercentile(np.abs(allb), 98))

fig, axes = plt.subplots(1, 3, figsize=(19, 4.4), sharey=True)
for ax, m in zip(axes, MODELS):
    lo, a = _disp(BIAS_MAPS[m])
    pm = ax.pcolormesh(lo, GLAT, a, cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto")
    ax.set_title("%s - GLODAP, 0-50 m DIC" % m, fontsize=10)
    ax.set_xlabel("longitude")
    ax.set_xlim(-180, 180)
    ax.set_xticks([-180, -120, -60, 0, 60, 120, 180])
axes[0].set_ylabel("latitude")
cb = fig.colorbar(pm, ax=axes, orientation="vertical", fraction=0.02, pad=0.01)
cb.set_label("DIC bias (µmol kg$^{-1}$)")
fig.suptitle("Figure 1 - Annual-mean 2002 surface (0-50 m) DIC bias against GLODAPv2.2016b", y=1.02)

# ---------------- Figure 2: vertical profile of the global mean DIC bias -----------
fig, ax = plt.subplots(figsize=(5.6, 6.0))
for m, c in zip(MODELS, ["#c0392b", "#1f77b4", "#2ca02c"]):
    prof = np.array([PROFILE[m][z] for z in DEPTHS])
    ax.plot(prof, DEPTHS, "-o", color=c, label=m, lw=1.8, ms=5)
    ax.plot(prof, DEPTHS, "o", color=c, ms=7, mfc="none", mew=1.2)   # sampled depths
ax.axvline(0.0, color="k", lw=1.0, ls="--")
ax.set_ylim(max(DEPTHS) * 1.05, 0)
ax.set_xlabel("DIC bias, model - GLODAP (µmol kg$^{-1}$)")
ax.set_ylabel("depth (m)")
ax.set_title("Figure 2 - Global mean DIC bias profile")
ax.grid(alpha=0.3)
ax.legend(frameon=False)
plt.tight_layout()
