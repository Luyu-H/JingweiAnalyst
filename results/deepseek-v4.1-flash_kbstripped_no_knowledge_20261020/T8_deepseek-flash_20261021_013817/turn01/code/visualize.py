# Visualization code only
# Figure 1: 0-50 m DIC bias maps, one panel per model, shared diverging scale
lats = dic_map_0_50['GLODAP'].coords['lat'].values
lons = dic_map_0_50['GLODAP'].coords['lon'].values

pooled = np.concatenate([
    np.abs(dic_bias_map[m].values).ravel() for m in models
])
vmax = float(np.nanpercentile(pooled, 98))
if not np.isfinite(vmax) or vmax == 0:
    vmax = float(np.nanmax(pooled)) if np.isfinite(np.nanmax(pooled)) and np.nanmax(pooled) > 0 else 1.0

fig, axes = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)
for ax, m in zip(axes, models):
    data = dic_bias_map[m].values
    pcm = ax.pcolormesh(
        lons, lats, data, shading='auto',
        cmap='RdBu_r', vmin=-vmax, vmax=vmax
    )
    ax.set_title(m)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_xlim([np.nanmin(lons), np.nanmax(lons)])
    ax.set_ylim([np.nanmin(lats), np.nanmax(lats)])

cbar = fig.colorbar(
    pcm, ax=axes, orientation='horizontal',
    fraction=0.05, pad=0.05, extend='both'
)
cbar.set_label('DIC bias (µmol kg$^{-1}$)')
fig.suptitle(
    f'0–50 m DIC bias vs GLODAP (model − GLODAP); '
    f'shared scale clipped at 98th percentile = {vmax:.0f} µmol kg$^{{-1}}$',
    fontsize=12
)

# Figure 2: vertical profile of global mean DIC bias
fig2, ax2 = plt.subplots(figsize=(6, 4.5), constrained_layout=True)
for m in models:
    xvals = [vertical_bias[m][d] for d in depths]
    ax2.plot(xvals, depths, marker='o', label=m)

ax2.axvline(0, color='k', linestyle='--', linewidth=1)
ax2.set_yticks(depths)
ax2.invert_yaxis()
ax2.set_xlabel('DIC bias (µmol kg$^{-1}$)')
ax2.set_ylabel('Depth (m)')
ax2.set_title('Vertical profile of global mean DIC bias (model − GLODAP)')
ax2.legend(loc='best')
ax2.grid(True, alpha=0.3)
