# ---------------------------------------------------------------------------
# Figure 1: global nitracline-depth maps (WOA23 and three models)
# ---------------------------------------------------------------------------
fig1, axes1 = plt.subplots(2, 2, figsize=(13, 8), sharex=True, sharey=True)
map_sources = ["WOA23"] + models
map_data = [woa_nit] + [model_nit_on_woa[m] for m in models]

all_map_vals = np.concatenate(
    [d[np.isfinite(d)].ravel() for d in map_data]
)
vmin_map = np.nanmin(all_map_vals)
vmax_map = np.nanmax(all_map_vals)

for ax, name, data in zip(axes1.ravel(), map_sources, map_data):
    im1 = ax.pcolormesh(
        woa_lon, woa_lat, data,
        cmap="YlGnBu", vmin=vmin_map, vmax=vmax_map, shading="auto"
    )
    ax.set_title(name)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

fig1.colorbar(
    im1, ax=axes1, label="Nitracline depth (m)", shrink=0.85
)
fig1.suptitle("Nitracline depth: WOA23 vs CMIP6 historical models")

# ---------------------------------------------------------------------------
# Figure 2: equatorial Pacific (2S-2N) nitrate longitude-depth sections
# ---------------------------------------------------------------------------
fig2, axes2 = plt.subplots(4, 1, figsize=(13, 11), sharex=True, sharey=True)
section_sources = ["WOA23"] + models
section_arrays = [section_data[s] for s in section_sources]

all_sec_vals = np.concatenate(
    [a[np.isfinite(a)].ravel() for a in section_arrays]
)
vmin_sec = 0.0
vmax_sec = np.nanpercentile(all_sec_vals, 99.0)

for ax, name, sec in zip(axes2, section_sources, section_arrays):
    im2 = ax.pcolormesh(
        section_data["lon"], section_data["depth"], sec.T,
        cmap="viridis", vmin=vmin_sec, vmax=vmax_sec, shading="auto"
    )
    ax.plot(
        section_data["lon"], section_data["nit_line"][name],
        color="k", linewidth=1.5, label="nitracline depth"
    )
    ax.set_ylim(1000, 0)
    ax.set_xlim(120, 280)
    ax.set_ylabel("Depth (m)")
    ax.set_title(name)
    ax.legend(loc="lower right")

axes2[-1].set_xlabel("Longitude (degrees east)")

fig2.colorbar(
    im2, ax=axes2,
    label=r"Nitrate ($\mu$mol kg$^{-1}$)", shrink=0.85
)
fig2.suptitle("Equatorial Pacific nitrate section (2°S–2°N)")
