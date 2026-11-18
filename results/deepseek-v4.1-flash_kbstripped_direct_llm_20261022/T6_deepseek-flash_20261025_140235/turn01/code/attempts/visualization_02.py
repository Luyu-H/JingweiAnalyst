# ---------------------------------------------------------------------------
# Figure 1: global maps of nitracline depth
# ---------------------------------------------------------------------------
fig1, axes1 = plt.subplots(2, 2, figsize=(14, 10))
map_fields = {
    'WOA23': nitracline_woa,
    'CESM2': nitracline_models['CESM2'],
    'GFDL-ESM4': nitracline_models['GFDL-ESM4'],
    'MPI-ESM1-2-LR': nitracline_models['MPI-ESM1-2-LR'],
}
map_titles = ['WOA23', 'CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']

all_map_vals = []
for f in map_fields.values():
    v = f[np.isfinite(f)]
    if v.size:
        all_map_vals.append(v)
all_map_vals = np.concatenate(all_map_vals)
vmin_map = np.nanpercentile(all_map_vals, 2)
vmax_map = np.nanpercentile(all_map_vals, 98)
cmap_map = 'YlGnBu'   # deeper -> darker/heavier

for ax, key, title in zip(axes1.ravel(), map_fields.keys(), map_titles):
    field = np.where(ocean_mask, map_fields[key], np.nan)
    pcm = ax.pcolormesh(
        woa_lon, woa_lat, field,
        cmap=cmap_map, vmin=vmin_map, vmax=vmax_map, shading='auto'
    )
    ax.set_title(title)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_facecolor('lightgray')

fig1.colorbar(
    pcm, ax=axes1, orientation='horizontal', fraction=0.05, pad=0.05,
    label='Nitracline depth (m)'
)

# ---------------------------------------------------------------------------
# Figure 2: equatorial Pacific longitude-depth nitrate sections
# ---------------------------------------------------------------------------
fig2, axes2 = plt.subplots(4, 1, figsize=(14, 16), sharex=True)
sec_keys = ['WOA23', 'CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']
sec_titles = ['WOA23', 'CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']

all_sec_vals = []
for key in sec_keys:
    depth = woa_depth if key == 'WOA23' else model_depths[key]
    sec = sections_umol[key]
    mask = depth <= 1000
    vals = sec[mask, :]
    vals = vals[np.isfinite(vals)]
    if vals.size:
        all_sec_vals.append(vals)
all_sec_vals = np.concatenate(all_sec_vals)
vmin_sec = 0.0
vmax_sec = float(np.nanpercentile(all_sec_vals, 99))
cmap_sec = 'viridis'

for ax, key, title in zip(axes2, sec_keys, sec_titles):
    depth = woa_depth if key == 'WOA23' else model_depths[key]
    mask = depth <= 1000
    sec_plot = sections_umol[key][mask, :]
    depth_plot = depth[mask]

    pcm2 = ax.pcolormesh(
        lon_centers, depth_plot, sec_plot,
        cmap=cmap_sec, vmin=vmin_sec, vmax=vmax_sec, shading='auto'
    )
    ax.plot(lon_centers, nitracline_lines[key], color='red', linewidth=1.5, label='Nitracline')
    ax.set_title(title)
    ax.set_ylabel('Depth (m)')
    ax.invert_yaxis()
    ax.set_xlim(120, 280)
    ax.set_ylim(1000, 0)
    ax.legend(loc='lower right')

axes2[-1].set_xlabel('Longitude')
axes2[-1].set_xticks([120, 160, 180, 200, 240, 280])
axes2[-1].set_xticklabels(['120E', '160E', '180', '160W', '120W', '80W'])

fig2.colorbar(
    pcm2, ax=axes2, orientation='horizontal', fraction=0.05, pad=0.05,
    label='Nitrate (µmol kg⁻¹)'
)
