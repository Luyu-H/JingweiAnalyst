# -----------------------------------------------------------------------------
# Figure 1: nitrate vs phosphate scatter, 200-1000 m, four sources
# -----------------------------------------------------------------------------
source_order = ["WOA23", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

scatter_data = {}
for name in source_order:
    s = sources[name]
    no3 = np.asarray(s["no3"], dtype=float)
    po4 = np.asarray(s["po4"], dtype=float)
    vol = np.asarray(s["vol"], dtype=float)
    m = (
        np.isfinite(no3)
        & np.isfinite(po4)
        & np.isfinite(vol)
        & (vol > 0)
        & (po4 > 0)
    )
    scatter_data[name] = (po4[m], no3[m])

xmax = max(np.nanmax(v[0]) for v in scatter_data.values() if v[0].size)
ymax = max(np.nanmax(v[1]) for v in scatter_data.values() if v[1].size)
xmax = 1.05 * xmax if np.isfinite(xmax) and xmax > 0 else 1.0
ymax = 1.05 * ymax if np.isfinite(ymax) and ymax > 0 else 1.0

fig1, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)

from matplotlib.ticker import MaxNLocator

for ax, name in zip(axes.ravel(), source_order):
    x, y = scatter_data[name]
    ax.scatter(
        x,
        y,
        s=2,
        alpha=0.30,
        c="steelblue",
        edgecolors="none",
        rasterized=True,
    )

    st = stats[name]
    xx = np.linspace(0.0, xmax, 200)
    ax.plot(
        xx,
        st["slope"] * xx + st["intercept"],
        color="crimson",
        lw=1.6,
        label="OLS fit",
    )
    ax.plot(xx, 16.0 * xx, color="black", ls="--", lw=1.3, label="Redfield 16:1")

    ax.set_xlim(0.0, xmax)
    ax.set_ylim(0.0, ymax)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.set_title(name)
    ax.set_xlabel(r"PO$_4$ ($\mu$mol kg$^{-1}$)")
    ax.set_ylabel(r"NO$_3$ ($\mu$mol kg$^{-1}$)")

    ax.text(
        0.03,
        0.97,
        f"a = {st['slope']:.2f}\nb = {st['intercept']:.2f}\nR$^2$ = {st['R2']:.3f}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8,
        bbox=dict(boxstyle="round", fc="white", alpha=0.85, ec="0.6"),
    )

    if name == source_order[0]:
        ax.legend(loc="lower right", fontsize=7, framealpha=0.9)

fig1.suptitle("200-1000 m NO$_3$ vs PO$_4$", fontsize=12)

# -----------------------------------------------------------------------------
# Figure 2: global maps of 200-1000 m mean N:P minus 16
# -----------------------------------------------------------------------------
lon_edges = np.arange(-180.0, 181.0, 1.0)
lat_edges = np.arange(-90.0, 91.0, 1.0)


def regrid_to_regular(lat2d, lon2d, values):
    lat = np.asarray(lat2d, dtype=float)
    lon = np.asarray(lon2d, dtype=float)
    val = np.asarray(values, dtype=float)

    lon_norm = np.mod(lon + 180.0, 360.0) - 180.0
    m = np.isfinite(val) & np.isfinite(lat) & np.isfinite(lon_norm)

    xi = np.floor((lat[m] - lat_edges[0]) / (lat_edges[1] - lat_edges[0])).astype(int)
    yi = np.floor((lon_norm[m] - lon_edges[0]) / (lon_edges[1] - lon_edges[0])).astype(int)
    xi = np.clip(xi, 0, lat_edges.size - 2)
    yi = np.clip(yi, 0, lon_edges.size - 2)

    sumv = np.zeros((lat_edges.size - 1, lon_edges.size - 1), dtype=float)
    cnt = np.zeros_like(sumv)

    np.add.at(sumv, (xi, yi), val[m])
    np.add.at(cnt, (xi, yi), 1.0)

    out = np.full_like(sumv, np.nan)
    nz = cnt > 0
    out[nz] = sumv[nz] / cnt[nz]
    return out


maps = {}
for name in source_order:
    s = sources[name]
    no3 = np.asarray(s["no3"], dtype=float)
    po4 = np.asarray(s["po4"], dtype=float)
    vol = np.asarray(s["vol"], dtype=float)

    ratio = np.full(no3.shape, np.nan, dtype=float)
    m = (
        np.isfinite(no3)
        & np.isfinite(po4)
        & (po4 > 0)
        & np.isfinite(vol)
        & (vol > 0)
    )
    ratio[m] = no3[m] / po4[m]

    maps[name] = regrid_to_regular(s["lat2d"], s["lon2d"], ratio)

all_anom = np.concatenate([(maps[name] - 16.0).ravel() for name in source_order])
vmax = np.nanmax(np.abs(all_anom))
if not np.isfinite(vmax) or vmax == 0:
    vmax = 1.0

from matplotlib.colors import Normalize

norm = Normalize(vmin=-vmax, vmax=vmax)

fig2, axes2 = plt.subplots(2, 2, figsize=(12, 7), constrained_layout=True)

for ax, name in zip(axes2.ravel(), source_order):
    field = maps[name] - 16.0
    im = ax.pcolormesh(
        lon_edges,
        lat_edges,
        field,
        shading="auto",
        cmap="RdBu_r",
        norm=norm,
    )
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_aspect("equal")
    ax.set_facecolor("lightgray")
    ax.set_title(name)
    ax.set_xticks([-180, -120, -60, 0, 60, 120, 180])
    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

fig2.colorbar(
    im,
    ax=axes2.ravel().tolist(),
    shrink=0.8,
    label="200-1000 m mean N:P - 16",
)
fig2.suptitle("Global 200-1000 m N:P anomaly relative to Redfield 16:1", fontsize=12)
