# ============================================================================
# Visualize: Figure 1 (NO3 vs PO4 scatter) and Figure 2 (N:P - 16 maps)
# ============================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# ---------------------------------------------------------------------------
# Figure 1: 4-panel NO3 vs PO4 scatter with OLS and Redfield lines
# ---------------------------------------------------------------------------
fig1, axes1 = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
axes1 = axes1.ravel()

# Determine common axis limits from all four sources
all_po4 = []
all_no3 = []
for src in source_order:
    s = sources.get(src)
    if s is None:
        continue
    valid = (
        np.isfinite(s["np_ratio"])
        & np.isfinite(s["nitrate"])
        & np.isfinite(s["phosphate"])
    )
    all_po4.append(s["phosphate"][valid])
    all_no3.append(s["nitrate"][valid])

all_po4 = np.concatenate(all_po4)
all_no3 = np.concatenate(all_no3)

x_lo, x_hi = np.nanpercentile(all_po4, [0.5, 99.5])
y_lo, y_hi = np.nanpercentile(all_no3, [0.5, 99.5])
x_lo = max(0.0, 0.95 * x_lo)
x_hi = 1.05 * x_hi
y_lo = max(0.0, 0.95 * y_lo)
y_hi = 1.05 * y_hi

xticks = np.linspace(x_lo, x_hi, 5)
yticks = np.linspace(y_lo, y_hi, 5)

for ax, src in zip(axes1, source_order):
    s = sources.get(src)
    if s is None:
        ax.set_visible(False)
        continue

    valid = (
        np.isfinite(s["np_ratio"])
        & np.isfinite(s["nitrate"])
        & np.isfinite(s["phosphate"])
    )
    ax.scatter(
        s["phosphate"][valid],
        s["nitrate"][valid],
        s=1.0,
        alpha=0.30,
        rasterized=True,
    )

    st = s["stats"]
    a = st["ols_slope"]
    b = st["ols_intercept"]
    r2 = st["r2"]

    xx = np.linspace(x_lo, x_hi, 200)
    ax.plot(xx, a * xx + b, color="C1", linewidth=1.8, label="OLS fit")
    ax.plot(
        xx,
        16.0 * xx,
        color="C2",
        linewidth=1.5,
        linestyle="--",
        label="Redfield 16:1",
    )

    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xlabel(r"PO$_4$ ($\mu$mol kg$^{-1}$)")
    ax.set_ylabel(r"NO$_3$ ($\mu$mol kg$^{-1}$)")
    ax.set_title(src)
    ax.text(
        0.03,
        0.97,
        f"a = {a:.2f}\nb = {b:.2f}\nR$^2$ = {r2:.3f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox=dict(facecolor="white", alpha=0.85, edgecolor="none"),
    )
    ax.legend(loc="lower right", fontsize=8)

fig1.suptitle("200-1000 m NO3 vs PO4: WOA23 and CMIP6 historical models")
plt.show()

# ---------------------------------------------------------------------------
# Figure 2: 2x2 maps of 200-1000 m mean N:P - 16
# ---------------------------------------------------------------------------
lon_reg = np.arange(-179.5, 180.0, 1.0)
lat_reg = np.arange(-89.5, 90.0, 1.0)
lon2d_reg, lat2d_reg = np.meshgrid(lon_reg, lat_reg)


def regrid_to_regular(lat2d, lon2d, field):
    """Nearest/linear interpolation of a curvilinear field to 1 deg lon-lat."""
    lat = np.asarray(lat2d).ravel()
    lon = np.asarray(lon2d).ravel()
    vals = np.asarray(field).ravel()

    lon = ((lon + 180.0) % 360.0) - 180.0
    valid = np.isfinite(vals) & np.isfinite(lat) & np.isfinite(lon)

    pts = np.column_stack([lon[valid], lat[valid]])
    vals_v = vals[valid]

    # Duplicate points across the dateline for interpolation near +/-180
    pts_all = np.vstack(
        [pts, pts + np.array([360.0, 0.0]), pts - np.array([360.0, 0.0])]
    )
    vals_all = np.concatenate([vals_v, vals_v, vals_v])

    reg = griddata(
        pts_all,
        vals_all,
        (lon2d_reg, lat2d_reg),
        method="linear",
        fill_value=np.nan,
    )
    if np.any(np.isnan(reg)):
        reg_nn = griddata(
            pts_all,
            vals_all,
            (lon2d_reg, lat2d_reg),
            method="nearest",
        )
        reg = np.where(np.isnan(reg), reg_nn, reg)
    return reg


# Shared diverging scale from pooled native-grid departures
all_dev = []
for src in source_order:
    s = sources.get(src)
    if s is None:
        continue
    v = s["np_ratio"][np.isfinite(s["np_ratio"])] - 16.0
    all_dev.append(v)
all_dev = np.concatenate(all_dev)
vmax = float(np.nanpercentile(np.abs(all_dev), 98))
if not np.isfinite(vmax) or vmax == 0.0:
    vmax = 1.0

fig2, axes2 = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
axes2 = axes2.ravel()

for ax, src in zip(axes2, source_order):
    s = sources.get(src)
    if s is None:
        ax.set_visible(False)
        continue

    if src == "WOA23":
        field = s["np_ratio"] - 16.0
        lon_plot = s["lon"][0, :]
        lat_plot = s["lat"][:, 0]
    else:
        field = regrid_to_regular(s["lat"], s["lon"], s["np_ratio"]) - 16.0
        lon_plot = lon_reg
        lat_plot = lat_reg

    im = ax.pcolormesh(
        lon_plot,
        lat_plot,
        field,
        shading="auto",
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
    )
    ax.set_xlim(-180.0, 180.0)
    ax.set_ylim(-90.0, 90.0)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(src)

fig2.colorbar(
    im,
    ax=axes2,
    orientation="vertical",
    fraction=0.025,
    pad=0.02,
    extend="both",
    label="N:P - 16",
)
fig2.suptitle(
    "200-1000 m mean N:P minus 16 (shared scale, clipped at 98th percentile)"
)
plt.show()
