# Figure 1: global nitracline-depth maps, shared sequential scale with deeper darker
plot_sources = ["WOA23"] + models
fields_plot = {
    src: nitri_fields[src].where(mask).compute() for src in plot_sources
}

all_vals = np.concatenate([fields_plot[src].values.ravel() for src in plot_sources])
vmax1 = float(np.nanpercentile(all_vals, 98))
vmin1 = 0.0

fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
axes = axes.ravel()
cmap1 = plt.cm.YlGnBu

for ax, src in zip(axes, plot_sources):
    f = fields_plot[src]
    pcm = ax.pcolormesh(
        f.lon,
        f.lat,
        f.values,
        shading="auto",
        cmap=cmap1,
        vmin=vmin1,
        vmax=vmax1,
    )
    ax.set_title(src)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

fig.colorbar(
    pcm,
    ax=axes.tolist(),
    label="Nitracline depth (m)",
    extend="max",
)
fig.suptitle(
    "Global nitracline depth (98th percentile colour clipping)",
    fontsize=12,
)

# Figure 2: equatorial Pacific nitrate sections, shared scale and nitracline overlays
plot_sources2 = ["WOA23"] + models
sec_plot = {
    src: section_fields[src].sel(depth=slice(0, 1000)) for src in plot_sources2
}

all_vals2 = np.concatenate([sec_plot[src].values.ravel() for src in plot_sources2])
vmax2 = float(np.nanpercentile(all_vals2, 98))
vmin2 = 0.0

fig2, axes2 = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
axes2 = axes2.ravel()
cmap2 = plt.cm.YlGnBu

for ax, src in zip(axes2, plot_sources2):
    sec = sec_plot[src]
    pcm2 = ax.pcolormesh(
        sec.lon,
        sec.depth,
        sec.values,
        shading="auto",
        cmap=cmap2,
        vmin=vmin2,
        vmax=vmax2,
    )

    line = section_nitri[src]
    ax.plot(line.lon, line.values, color="red", lw=1.6)

    ax.set_xlim(120, 280)
    ax.set_ylim(1000, 0)  # depth increases downward
    ax.set_title(src)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Depth (m)")
    ax.set_xticks([120, 180, 240, 280])
    ax.set_xticklabels(["120°E", "180°", "120°W", "80°W"])

fig2.colorbar(
    pcm2,
    ax=axes2.tolist(),
    label="Nitrate (mol m$^{-3}$)",
    extend="max",
)
fig2.suptitle(
    "Equatorial Pacific (2°S–2°N) nitrate; red line = nitracline "
    "(98th percentile colour clipping)",
    fontsize=12,
)
