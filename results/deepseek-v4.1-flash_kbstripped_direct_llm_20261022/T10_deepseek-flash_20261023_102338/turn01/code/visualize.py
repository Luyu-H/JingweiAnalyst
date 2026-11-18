# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------
def _plot_field(ax, d, field, **kwargs):
    lon = np.asarray(d["lon"], dtype=float)
    lat = np.asarray(d["lat"], dtype=float)
    data = np.asarray(d[field], dtype=float)

    if lon.ndim == 2 and lat.ndim == 2:
        # For a regular lon/lat grid, sort longitude so the wrap-around seam is less ugly.
        if np.allclose(lon[0, :], lon[1, :]) if lon.shape[0] > 1 else True:
            order = np.argsort(((lon[0, :] + 180.0) % 360.0) - 180.0)
            lonp = ((lon[0, order] + 180.0) % 360.0) - 180.0
            latp = lat[:, 0]
            datap = data[:, order]
            ax.pcolormesh(lonp, latp, datap, shading="auto", **kwargs)
            return
        lonp = ((lon + 180.0) % 360.0) - 180.0
        ax.pcolormesh(lonp, lat, data, shading="auto", **kwargs)
    else:
        ax.pcolormesh(lon, lat, data, shading="auto", **kwargs)

    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

# ---------------------------------------------------------------------------
# Figure 1: historical OMZ thickness at O2 < 60 umol/kg
# ---------------------------------------------------------------------------
hist_keys = [
    ("GFDL-ESM4", "GFDL-ESM4 historical"),
    ("MPI-ESM1-2-LR", "MPI-ESM1-2-LR historical"),
    ("WOA23", "WOA23 1971-2000"),
]
hist_available = [(k, t) for k, t in hist_keys if k in hist_maps]
vmax_hist = max(
    np.nanmax(hist_maps[k]["thickness60"]) for k, _ in hist_available
)

fig1, axes1 = plt.subplots(
    1, len(hist_available), figsize=(6.0 * len(hist_available), 4.6),
    constrained_layout=True, squeeze=False
)
axes1 = axes1.ravel()
p1 = None
for ax, (key, title) in zip(axes1, hist_available):
    d = hist_maps[key]
    p1 = _plot_field(
        ax, d, "thickness60",
        cmap="viridis", vmin=0.0, vmax=vmax_hist
    )
    ax.set_title(title)
if p1 is not None:
    fig1.colorbar(
        p1, ax=list(axes1), shrink=0.82,
        label="OMZ thickness for O2 < 60 $\\mu$mol kg$^{-1}$ (m)"
    )

# ---------------------------------------------------------------------------
# Figure 2: ssp245 - historical change, plus CESM2 ssp245 alone
# ---------------------------------------------------------------------------
change_keys = [k for k in ["GFDL-ESM4", "MPI-ESM1-2-LR"] if k in change_maps]
vmax_change = max(
    np.nanmax(np.abs(change_maps[k]["thickness60_change"]))
    for k in change_keys
) if change_keys else 1.0

fig2, axes2 = plt.subplots(
    1, 3, figsize=(18.0, 4.6), constrained_layout=True, squeeze=False
)
axes2 = axes2.ravel()

p_change = None
for ax, model in zip(axes2[:2], ["GFDL-ESM4", "MPI-ESM1-2-LR"]):
    if model in change_maps:
        d = change_maps[model]
        p_change = _plot_field(
            ax, d, "thickness60_change",
            cmap="RdBu_r", vmin=-vmax_change, vmax=vmax_change
        )
        ax.set_title(f"{model} ssp245 2071-2100 minus historical 1971-2000")
    else:
        ax.set_title(f"{model}: change unavailable")
        ax.axis("off")

# Third panel: CESM2 ssp245 thickness, no historical comparison possible.
ax = axes2[2]
if "CESM2" in ssp_maps:
    d = ssp_maps["CESM2"]
    p_cesm = _plot_field(
        ax, d, "thickness60",
        cmap="viridis", vmin=0.0,
        vmax=float(np.nanmax(d["thickness60"]))
    )
    ax.set_title(
        "CESM2 ssp245 2071-2100 OMZ thickness\n"
        "no historical o2; change cannot be formed"
    )
    fig2.colorbar(
        p_cesm, ax=ax, shrink=0.82,
        label="OMZ thickness for O2 < 60 $\\mu$mol kg$^{-1}$ (m)"
    )
else:
    ax.set_title("CESM2 ssp245 unavailable")
    ax.axis("off")

if p_change is not None:
    fig2.colorbar(
        p_change, ax=axes2[:2], shrink=0.82,
        label="OMZ thickness change (m)"
    )
