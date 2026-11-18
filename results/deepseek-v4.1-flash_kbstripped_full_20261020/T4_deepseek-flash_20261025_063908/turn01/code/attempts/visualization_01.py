# Figure 1: Maps of operator MLD minus mlotst
all_diffs = []
for model in models:
    for month in months:
        d = map_data[(model, month)]['diff']
        all_diffs.append(d[np.isfinite(d)])
all_diffs = np.concatenate(all_diffs)
vmax = float(np.nanpercentile(np.abs(all_diffs), 98))
if vmax == 0:
    vmax = 1.0

fig1, axes1 = plt.subplots(3, 2, figsize=(12, 10), constrained_layout=True)
for i, model in enumerate(models):
    for j, month in enumerate(months):
        ax = axes1[i, j]
        data = map_data[(model, month)]
        lon = data['lon']
        lat = data['lat']
        diff = data['diff']
        im = ax.pcolormesh(lon, lat, diff, cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='auto')
        ax.set_title(f'{model} - {"March" if month==3 else "September"}')
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        ax.set_aspect('equal')
        ax.grid(True, ls=':', lw=0.3, alpha=0.5)
fig1.colorbar(im, ax=axes1, orientation='horizontal', fraction=0.05, pad=0.05, label='MLD difference (m)')
fig1.suptitle(f'Operator MLD minus mlotst (shared scale, 98th percentile clip ±{vmax:.1f} m)', fontsize=14)

# Figure 2: Scatter plots
all_x = np.concatenate([scatter_data[k][0] for k in scatter_data])
all_y = np.concatenate([scatter_data[k][1] for k in scatter_data])
all_vals = np.concatenate([all_x, all_y])
min_val = np.nanmin(all_vals)
max_val = np.nanmax(all_vals)
xlim = (min_val * 0.8, max_val * 1.2)
ylim = (min_val * 0.8, max_val * 1.2)

# Generate nice log ticks
ticks = []
for exp in range(int(np.floor(np.log10(xlim[0]))), int(np.ceil(np.log10(xlim[1])))+1):
    for mult in [1, 2, 5]:
        val = mult * 10**exp
        if xlim[0] <= val <= xlim[1]:
            ticks.append(val)
ticks = np.array(ticks)

fig2, axes2 = plt.subplots(3, 2, figsize=(12, 10), constrained_layout=True)
for i, model in enumerate(models):
    for j, month in enumerate(months):
        ax = axes2[i, j]
        x, y = scatter_data[(model, month)]
        ax.scatter(x, y, s=1, alpha=0.3, edgecolors='none', rasterized=True)
        ax.plot([xlim[0], xlim[1]], [ylim[0], ylim[1]], 'k--', lw=1)  # y=x
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels([f'{t:g}' for t in ticks], rotation=45, fontsize=8)
        ax.set_yticklabels([f'{t:g}' for t in ticks], fontsize=8)
        ax.set_xlabel('Operator MLD (m)')
        ax.set_ylabel('mlotst (m)')
        ax.set_title(f'{model} - {"March" if month==3 else "September"}')
        r = results[(model, month)]
        text_str = f'Bias: {r["bias"]:.1f} m\nRMSE: {r["rmse"]:.1f} m\nCorr: {r["corr"]:.3f}'
        ax.text(0.05, 0.95, text_str, transform=ax.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), fontsize=9)
        ax.grid(True, which='both', ls=':', lw=0.5, alpha=0.5)
