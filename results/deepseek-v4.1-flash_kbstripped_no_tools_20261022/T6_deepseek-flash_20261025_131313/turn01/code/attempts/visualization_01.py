# -----------------------------------------------------------------------------
# Figure 1: global maps of nitracline depth on the common WOA 1-degree grid
# -----------------------------------------------------------------------------
sources = ["WOA23"] + models

fields_for_scale = [woa_nitra_shared] + [model_nitra_shared[m] for m in models]
pooled_nitra = np.concatenate([
    f[np.isfinite(f)] for f in fields_for_scale
])
vmin_nitra = 0.0
vmax_nitra = float(np.nanpercentile(pooled_nitra, 98))

fig1, axes1 = plt.subplots(
    2, 2, figsize=(12, 8), constrained_layout=True
)

plot_fields = {
    "WOA23": woa_nitra_shared,
    "CESM2": model_nitra_shared["CESM2"],
    "GFDL-ESM4": model_nitra_shared["GFDL-ESM4"],
    "MPI-ESM1-2-LR": model_nitra_shared["MPI-ESM1-2-LR"],
}

for ax, source in zip(axes1.ravel(), sources):
    pcm1 = ax.pcolormesh(
        woa_lon, woa_lat, plot_fields[source],
        shading="auto", cmap="Blues",
        vmin=vmin_nitra, vmax=vmax_nitra,
    )
    ax.set_title(source)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_aspect("auto")

cbar1 = fig1.colorbar(
    pcm1, ax=axes1, label="Nitracline depth (m)", extend="max"
)
fig1.suptitle(
    "Nitracline depth (m), shared 0–98th percentile scale; land masked"
)

# -----------------------------------------------------------------------------
# Figure 2: equatorial Pacific nitrate sections, 2S-2N, 120E-80W, top 1000 m
# -----------------------------------------------------------------------------
section_sources = ["WOA23", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

pooled_section = np.concatenate([
    section_display[s].values[np.isfinite(section_display[s].values)]
    for s in section_sources
])
vmin_sec = float(np.nanpercentile(pooled_section, 2))
vmax_sec = float(np.nanpercentile(pooled_section, 98))

fig2, axes2 = plt.subplots(
    2, 2, figsize=(13, 8), constrained_layout=True
)

for ax, source in zip(axes2.ravel(), section_sources):
    sec = section_display[source]
    pcm2 = ax.pcolormesh(
        sec.lon.values, sec.depth.values, sec.values,
        shading="auto", cmap="viridis",
        vmin=vmin_sec, vmax=vmax_sec,
    )
    ax.plot(
        target_lon, nitra_line_section[source],
        color="red", linewidth=1.3,
    )
    ax.set_ylim(1000, 0)
    ax.set_xlim(120, 280)
    ax.set_title(f"{source}  (red = nitracline)")
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Depth (m)")

cbar2 = fig2.colorbar(
    pcm2, ax=axes2, label="Nitrate (µmol kg$^{-1}$)", extend="both"
)
fig2.suptitle(
    "Equatorial Pacific (2°S–2°N) nitrate; shared 2–98th percentile scale"
)
