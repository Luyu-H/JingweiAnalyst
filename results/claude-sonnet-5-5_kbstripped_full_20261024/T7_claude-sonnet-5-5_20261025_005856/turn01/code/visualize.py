fig, axes = plt.subplots(2, 2, figsize=(11, 9), sharex=True, sharey=True, constrained_layout=True)
xmax = float(np.ceil(max(scatter[s]['x_p995'] for s in sources) * 1.1 * 2) / 2)
ymax = float(np.ceil(max(scatter[s]['y_p995'] for s in sources) * 1.1 / 5) * 5)
xx = np.linspace(0, xmax, 50)
for ax, s in zip(axes.ravel(), sources):
    sc = scatter[s]
    f = fits[s]
    ax.scatter(sc['x'], sc['y'], s=3, alpha=0.2, color='tab:blue', rasterized=True)
    ax.plot(xx, f['a'] * xx + f['b_umol_kg'], color='tab:red', lw=2, ls='-', label='OLS fit')
    ax.plot(xx, 16 * xx, color='k', lw=1.8, ls='--', label='Redfield 16:1')
    ax.text(0.03, 0.97, f"a = {f['a']:.2f}\nb = {f['b_umol_kg']:.2f} µmol/kg\nR² = {f['r2']:.3f}",
            transform=ax.transAxes, va='top', ha='left', fontsize=10,
            bbox=dict(boxstyle='round', fc='white', ec='gray', alpha=0.85))
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    ax.set_title(s)
    ax.legend(loc='lower right', fontsize=9)
    ax.grid(alpha=0.3)
for ax in axes[1, :]:
    ax.set_xlabel('PO4 (µmol/kg)')
for ax in axes[:, 0]:
    ax.set_ylabel('NO3 (µmol/kg)')
fig.suptitle('NO3 vs PO4, 200-1000 m, shared 1° mask (models at ρ = 1025 kg/m³)')

pooled = np.concatenate([map_data[s]['field'].ravel() for s in sources])
pooled = pooled[np.isfinite(pooled)]
vmax = float(np.nanpercentile(np.abs(pooled), 98))
fig2, axes2 = plt.subplots(2, 2, figsize=(14, 8), constrained_layout=True)
for ax, s in zip(axes2.ravel(), sources):
    md = map_data[s]
    pm = ax.pcolormesh(md['lon'], md['lat'], md['field'], cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='auto')
    ax.set_facecolor('lightgray')
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_aspect('equal')
    ax.set_title(f"{s}  (headline N:P = {stats[s]['np_headline_cell_by_cell']:.2f})")
    ax.set_xlabel('Longitude (°E)')
    ax.set_ylabel('Latitude (°N)')
    ax.grid(alpha=0.3, lw=0.5)
cb = fig2.colorbar(pm, ax=axes2, orientation='vertical', shrink=0.8, extend='both')
cb.set_label('N:P − 16')
fig2.suptitle(f'200-1000 m mean N:P minus 16, regridded to 1° for display; scale clipped at 98th percentile (±{vmax:.1f})')
plt.show()
