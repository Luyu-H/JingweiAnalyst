# Figure 1: 2x2 global nitracline depth maps on shared sequential colour scale
map_fields = {
    "WOA23": woa_nitracline,
    "CESM2": model_nitracline_regridded["CESM2"],
    "GFDL-ESM4": model_nitracline_regridded["GFDL-ESM4"],
    "MPI-ESM1-2-LR": model_nitracline_regridded["MPI-ESM1-2-LR"],
}

map_vals = np.concatenate([np.ravel(f.values[np.isfinite(f.values)]) for f in map_fields.values()])
vmin_map = 0.0
vmax_map = float(np.nanpercentile(map_vals, 98))

fig1, axes1 = plt.subplots(2, 2, figsize=(16, 8), constrained_layout=True)
axes1 = axes1.ravel()

for ax, name in zip(axes1, ["WOA23", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]):
    da = map_fields[name]
    im = ax.pcolormesh(da.lon, da.lat, da, cmap="viridis", vmin=vmin_map, vmax=vmax_map, shading="auto")
    ax.set_title(name)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(float(da.lon.min()), float(da.lon.max()))
    ax.set_ylim(float(da.lat.min()), float(da.lat.max()))

cbar1 = fig1.colorbar(im, ax=axes1, orientation="vertical", shrink=0.9, extend="max")
cbar1.set_label("Nitracline depth (m)")
fig1.suptitle("Figure 1. Global nitracline depth maps (shared scale, clipped at 98th percentile)")

# Figure 2: 4-panel equatorial Pacific longitude-depth nitrate sections with nitracline overlay
# Shared colour scale across sections
section_arrays = []
for name in section_order:
    sec = section_fields[name]
    vals = sec.values
    section_arrays.append(vals[np.isfinite(vals)])
section_vals = np.concatenate(section_arrays)
vmin_sec = float(np.nanpercentile(section_vals, 2))
vmax_sec = float(np.nanpercentile(section_vals, 98))

fig2, axes2 = plt.subplots(4, 1, figsize=(14, 14), constrained_layout=True, sharex=False)
for ax, name in zip(axes2, section_order):
    sec = section_fields[name]

    # Plot section on a continuous 0-360 longitude axis for the dateline-crossing Pacific
    if "lon" in sec.coords:
        sec_plot = sec.assign_coords(lon=sec.lon % 360).sortby("lon")
        lon = sec_plot.lon
        depth = sec_plot.depth
        pcm = ax.pcolormesh(lon, depth, sec_plot, cmap="viridis", vmin=vmin_sec, vmax=vmax_sec, shading="auto")
    else:
        lon_name = sec.attrs.get("lon_coord", "lon")
        lon = sec[lon_name]
        depth = sec.depth
        pcm = ax.pcolormesh(lon, depth, sec, cmap="viridis", vmin=vmin_sec, vmax=vmax_sec, shading="auto")

    ncl = section_ncl_lines[name]
    if "lon" in ncl.coords:
        ncl_plot = ncl.assign_coords(lon=ncl.lon % 360).sortby("lon")
        ax.plot(ncl_plot.lon, ncl_plot, color="k", linewidth=2, label="Nitracline")
    else:
        lon_name_ncl = ncl.attrs.get("lon_coord", "lon")
        ax.plot(ncl[lon_name_ncl], ncl, color="k", linewidth=2, label="Nitracline")

    ax.set_ylim(1000, 0)
    ax.set_xlim(120, 280)
    ax.set_ylabel("Depth (m)")
    ax.set_title(name)
    ax.legend(loc="best")

axes2[-1].set_xlabel("Longitude (°E, 0-360)")
cbar2 = fig2.colorbar(pcm, ax=axes2, orientation="vertical", shrink=0.95, extend="both")
cbar2.set_label("Nitrate (umol kg$^{-1}$)")
fig2.suptitle("Figure 2. Equatorial Pacific nitrate sections, 2°S-2°N mean (shared scale, clipped at 2nd-98th percentiles)")
