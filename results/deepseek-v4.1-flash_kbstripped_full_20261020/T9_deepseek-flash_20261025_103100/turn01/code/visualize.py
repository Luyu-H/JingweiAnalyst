import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def prep_pdo_pattern(da):
    """Put a regular 1-degree PDO-domain field on a continuous 110-260E axis."""
    p = da.assign_coords(lon=da.lon % 360).sortby("lon")
    if p.lat[0] > p.lat[-1]:
        p = p.sortby("lat")
    p = p.sel(lat=slice(20, 70), lon=slice(110, 260))
    return p


# ------------------------------------------------------------------
# Figure 1: regression patterns.
# ------------------------------------------------------------------
order = ["ERSSTv5", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
pooled = []
for name in order:
    p = prep_pdo_pattern(results[name]["reg_pattern"])
    pooled.append(np.asarray(p.values).ravel())
pooled = np.concatenate(pooled)
vmax = float(np.nanpercentile(np.abs(pooled), 98))
if not np.isfinite(vmax) or vmax <= 0:
    vmax = float(np.nanmax(np.abs(pooled)))
levels = np.linspace(-vmax, vmax, 21)

fig1, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
axes = axes.ravel()

cf = None
for ax, name in zip(axes, order):
    p = prep_pdo_pattern(results[name]["reg_pattern"])
    cf = ax.contourf(
        p.lon,
        p.lat,
        p.values,
        levels=levels,
        cmap="RdBu_r",
        extend="both",
    )
    ax.set_title(name)
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")
    ax.set_xlim(110, 260)
    ax.set_ylim(20, 70)

    text = f"Var. expl. = {results[name]['explained_variance']:.1f}%"
    if name != "ERSSTv5":
        text += f"\nPattern corr. = {results[name]['pattern_corr_vs_ersst']:.3f}"
    ax.text(
        0.02,
        0.98,
        text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.75, edgecolor="none"),
    )

cbar = fig1.colorbar(cf, ax=axes, orientation="vertical", shrink=0.9, pad=0.02)
cbar.set_label("°C per standard deviation")
fig1.suptitle(
    "PDO regression patterns (bilinear 1° grid; shared diverging scale, 98th-percentile clipping)",
    fontsize=12,
)

# ------------------------------------------------------------------
# Figure 2: standardized indices and smoothed ERSSTv5 index.
# ------------------------------------------------------------------
fig2, ax2 = plt.subplots(figsize=(12, 5.5), constrained_layout=True)

colors = {
    "ERSSTv5": "black",
    "CESM2": "tab:blue",
    "GFDL-ESM4": "tab:orange",
    "MPI-ESM1-2-LR": "tab:green",
    "Published PDO": "tab:red",
}

for name in order:
    idx = results[name]["index"]
    ax2.plot(
        idx.time,
        idx.values,
        label=name,
        color=colors[name],
        lw=1.0,
        alpha=0.75,
    )

ax2.plot(
    pdo_pub_std.time,
    pdo_pub_std.values,
    label="Published PDO",
    color=colors["Published PDO"],
    lw=1.2,
    alpha=0.8,
    linestyle="--",
)

sm = ersst_smoothed.dropna(dim="time")
ax2.plot(
    sm.time,
    sm.values,
    label="ERSSTv5 9-yr running mean",
    color="black",
    lw=2.6,
)

ymin, ymax = ax2.get_ylim()
ax2.set_ylim(ymin, max(ymax, 2.7))
for yr in transition_years:
    x = pd.Timestamp(f"{yr}-07-01")
    ax2.axvline(x, color="gray", ls=":", lw=0.9, alpha=0.75)
    ax2.text(
        x,
        ax2.get_ylim()[1] * 0.96,
        str(yr),
        rotation=90,
        fontsize=8,
        ha="center",
        va="top",
        color="gray",
    )

ax2.set_xlabel("Year")
ax2.set_ylabel("Standardized index")
ax2.set_title(
    "Standardized PDO indices and 9-year running mean of ERSSTv5 "
    "(detected transitions marked)"
)
ax2.set_xlim(pd.Timestamp("1950-01-01"), pd.Timestamp("2014-12-31"))
ax2.grid(True, alpha=0.3)
ax2.legend(loc="upper left", fontsize=8, ncol=2, framealpha=0.9)
