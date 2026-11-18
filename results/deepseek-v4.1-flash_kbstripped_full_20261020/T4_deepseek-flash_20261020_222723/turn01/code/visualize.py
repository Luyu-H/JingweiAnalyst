# Figure 1: 3 x 2 difference maps on a common longitude-latitude display grid.
fig1, axes1 = plt.subplots(3, 2, figsize=(12, 10), constrained_layout=True)
for i, model in enumerate(models):
    for j, month in enumerate(months):
        ax = axes1[i, j]
        d = diff_regridded[(model, month)]
        lon = d.lon.values
        lat = d.lat.values

        pcm = ax.pcolormesh(
            lon, lat, d.values,
            cmap='RdBu_r', vmin=-fig1_vmax, vmax=fig1_vmax,
            shading='auto'
        )
        ax.set_title(
            f"{model} ({grid_labels[model]}) - {month_names[month]}\n"
            f"bias={store[(model, month)]['bias']:.1f} m, "
            f"RMSD={store[(model, month)]['rmse']:.1f} m, "
            f"r={store[(model, month)]['corr']:.3f}",
            fontsize=9
        )
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.set_xlim(float(np.nanmin(lon)), float(np.nanmax(lon)))
        ax.set_ylim(float(np.nanmin(lat)), float(np.nanmax(lat)))
        ax.set_aspect('auto')

fig1.colorbar(
    pcm, ax=axes1, orientation='horizontal', fraction=0.04, pad=0.06,
    label='Operator MLD - mlotst (m)'
)
fig1.suptitle(
    f"Figure 1: Operator MLD minus published mlotst, 1971-2000 climatologies "
    f"(shared scale; 98th pct |diff| = {fig1_vmax:.1f} m)",
    fontsize=12
)

# Figure 2: 3 x 2 cell-wise scatter plots with common logarithmic axes.
all_positive = []
for model in models:
    for month in months:
        x = store[(model, month)]['op'].values.ravel()
        y = store[(model, month)]['pub'].values.ravel()
        v = np.concatenate([x, y])
        v = v[np.isfinite(v) & (v > 0)]
        if v.size:
            all_positive.append(v)

if all_positive:
    all_positive = np.concatenate(all_positive)
    vmin = float(np.nanpercentile(all_positive, 1))
    vmax = float(np.nanpercentile(all_positive, 99.5))
    vmin = max(vmin, 0.1)
    vmax = max(vmax, vmin * 10)
else:
    vmin, vmax = 0.1, 1000.0

tick_candidates = np.array([0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000])
ticks = [float(t) for t in tick_candidates if vmin <= t <= vmax]
if not ticks:
    ticks = [vmin, vmax]

fig2, axes2 = plt.subplots(3, 2, figsize=(11, 11), constrained_layout=True)
for i, model in enumerate(models):
    for j, month in enumerate(months):
        ax = axes2[i, j]
        op = store[(model, month)]['op']
        pub = store[(model, month)]['pub']

        x = op.values.ravel()
        y = pub.values.ravel()
        valid = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)

        ax.scatter(
            x[valid], y[valid],
            s=1, alpha=0.25, color='C0', rasterized=True
        )
        ax.plot([vmin, vmax], [vmin, vmax], 'k--', lw=1.0)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(vmin, vmax)
        ax.set_ylim(vmin, vmax)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels([f'{t:g}' for t in ticks])
        ax.set_yticklabels([f'{t:g}' for t in ticks])

        ax.set_xlabel('Operator MLD (m)')
        ax.set_ylabel('mlotst (m)')
        ax.set_title(
            f"{model} ({grid_labels[model]}) - {month_names[month]}",
            fontsize=9
        )

        txt = (
            f"bias={store[(model, month)]['bias']:.1f} m\n"
            f"RMSD={store[(model, month)]['rmse']:.1f} m\n"
            f"r={store[(model, month)]['corr']:.3f}"
        )
        ax.text(
            0.03, 0.97, txt,
            transform=ax.transAxes, va='top', ha='left', fontsize=8,
            bbox=dict(boxstyle='round', fc='white', ec='0.8', alpha=0.9)
        )

fig2.suptitle(
    "Figure 2: Operator MLD vs published mlotst, 1971-2000 March and September climatologies",
    fontsize=12
)
