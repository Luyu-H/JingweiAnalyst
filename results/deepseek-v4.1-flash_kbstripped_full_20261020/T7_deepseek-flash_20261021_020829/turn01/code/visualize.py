import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

SOURCES = ["WOA23", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
COLORS = {"WOA23": "black", "CESM2": "tab:blue",
          "GFDL-ESM4": "tab:green", "MPI-ESM1-2-LR": "tab:red"}

# ----------------------------------------------------------------------------
# Figure 1 - 4-panel NO3 vs PO4 scatter with fitted and Redfield lines
# ----------------------------------------------------------------------------
pool_x = np.concatenate([FIELDS[n]["scatter_x"] for n in SOURCES])
pool_y = np.concatenate([FIELDS[n]["scatter_y"] for n in SOURCES])
xmax = float(np.nanpercentile(pool_x, 99.5)) * 1.05
ymax = float(np.nanpercentile(pool_y, 99.5)) * 1.05
xmax = max(xmax, 1.0)
ymax = max(ymax, 16.0 * xmax * 0.6)

fig1, axes1 = plt.subplots(2, 2, figsize=(11, 10), constrained_layout=True)
axes1 = axes1.ravel()

for ax, name in zip(axes1, SOURCES):
    d = FIELDS[name]
    ax.scatter(d["scatter_x"], d["scatter_y"], s=2, alpha=0.20,
               color=COLORS[name], linewidths=0, rasterized=True,
               label="cells")
    xx = np.linspace(0, xmax, 200)
    ax.plot(xx, d["slope"] * xx + d["intercept"], color="tab:orange", lw=2,
            label=f"fit: a={d['slope']:.2f}, b={d['intercept']:.2f}")
    ax.plot(xx, 16.0 * xx, color="grey", ls="--", lw=1.6, label="Redfield 16:1")
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    ax.set_xlabel("Phosphate (µmol kg$^{-1}$)")
    ax.set_ylabel("Nitrate (µmol kg$^{-1}$)")
    ax.set_title(f"{name}   R² = {d['r2']:.3f}")
    ax.text(0.97, 0.05,
            f"a = {d['slope']:.2f}\nb = {d['intercept']:.2f} µmol kg$^{{-1}}$\n"
            f"R² = {d['r2']:.3f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
            bbox=dict(boxstyle="round", fc="white", alpha=0.85, ec="0.6"))
    ax.legend(loc="upper left", fontsize=7, framealpha=0.85)
    ax.grid(alpha=0.25)

fig1.suptitle("Nitrate vs phosphate, 200-1000 m (each source on its own grid)",
              fontsize=12)
fig1.supxlabel("Model concentrations converted to µmol kg$^{-1}$ at ρ = 1025 kg m$^{-3}$",
               fontsize=9)

# ----------------------------------------------------------------------------
# Figure 2 - global maps of 200-1000 m mean N:P minus 16
# ----------------------------------------------------------------------------
ref = DISPLAY_MAPS["WOA23"]
lon = np.asarray(ref["lon"].values, dtype=float)
lat = np.asarray(ref["lat"].values, dtype=float)
if lon[0] > lon[-1]:
    order = np.argsort(lon)
    lon = lon[order]
else:
    order = slice(None)

anom = {}
for name in SOURCES:
    f = DISPLAY_MAPS[name].transpose("lat", "lon")
    arr = np.asarray(f.compute().values, dtype=float)
    if not isinstance(order, slice):
        arr = arr[:, order]
    anom[name] = arr - 16.0

pooled = np.concatenate([np.abs(a[np.isfinite(a)]).ravel() for a in anom.values()])
vmax = float(np.nanpercentile(pooled, 98))
vmax = max(vmax, 1.0)

norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
cmap = plt.get_cmap("RdBu_r").copy()
cmap.set_bad("0.85")

fig2, axes2 = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
axes2 = axes2.ravel()

pm = None
for ax, name in zip(axes2, SOURCES):
    pm = ax.pcolormesh(lon, lat, np.ma.masked_invalid(anom[name]),
                       cmap=cmap, norm=norm, shading="auto")
    ax.set_title(f"{name}   (mean = {FIELDS[name]['headline']:.2f})")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-80, 90)
    ax.axhline(0, color="0.3", lw=0.5, alpha=0.5)

cbar = fig2.colorbar(pm, ax=axes2, orientation="vertical", extend="both",
                     shrink=0.9, pad=0.02)
cbar.set_label("N:P anomaly (unitless)")

fig2.suptitle("Mean N:P minus 16, 200-1000 m (shared scale, clipped at 98th "
              "percentile of |anomaly|; models regridded to the WOA23 1° grid, "
              "conservative)", fontsize=11)

plt.show()
