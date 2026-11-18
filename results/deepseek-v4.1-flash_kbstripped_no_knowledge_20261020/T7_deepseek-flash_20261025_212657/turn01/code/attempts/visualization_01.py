# -----------------------------------------------------------------------------
# Figure 1: Nitrate vs phosphate scatter with OLS and Redfield lines
# -----------------------------------------------------------------------------
fig1, axes1 = plt.subplots(2, 2, figsize=(12, 10))
axes1 = axes1.flatten()

for i, src in enumerate(sources):
    ax = axes1[i]
    x, y = figure1_data[src]
    n = len(x)
    if n > 5000:
        idx = np.linspace(0, n - 1, 5000, dtype=int)
        xp, yp = x[idx], y[idx]
    else:
        xp, yp = x, y

    ax.scatter(xp, yp, s=2, alpha=0.30, color='gray', label='data')

    xx = np.linspace(0, fig1_xmax, 200)
    a = figure1_reg[src]['a']
    b = figure1_reg[src]['b']
    r2 = figure1_reg[src]['r2']
    ax.plot(xx, a * xx + b, color='C0', lw=2.2, label='OLS fit')
    ax.plot(xx, 16.0 * xx, color='C1', lw=2.2, ls='--', label='Redfield 16:1')

    ax.set_xlim(0, fig1_xmax)
    ax.set_ylim(0, fig1_ymax)
    ax.set_xticks(np.linspace(0, fig1_xmax, 6))
    ax.set_yticks(np.linspace(0, fig1_ymax, 6))
    ax.set_xlabel('Phosphate (µmol kg$^{-1}$)')
    ax.set_ylabel('Nitrate (µmol kg$^{-1}$)')
    ax.set_title(src)
    ax.grid(True, alpha=0.30)

    ax.text(0.03, 0.97,
            f'a = {a:.2f}\nb = {b:.2f}\nR² = {r2:.3f}',
            transform=ax.transAxes, va='top', ha='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))
    ax.legend(loc='lower right', fontsize=8)

fig1.suptitle('Nitrate vs phosphate in 200–1000 m', fontsize=14)
fig1.tight_layout()

# -----------------------------------------------------------------------------
# Figure 2: Global maps of mean N:P minus 16
# -----------------------------------------------------------------------------
fig2, axes2 = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
axes2 = axes2.flatten()

for i, src in enumerate(sources):
    ax = axes2[i]
    field = figure2_data[src] - 16.0
    if field.dims != ('lat', 'lon'):
        field = field.transpose('lat', 'lon')

    lon = field.lon.values
    lat = field.lat.values
    data = field.values

    im = ax.pcolormesh(lon, lat, data, cmap='RdBu_r',
                       vmin=-fig2_vmax, vmax=fig2_vmax, shading='auto')
    ax.set_title(src)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xticks(np.arange(-180, 181, 60))
    ax.set_yticks(np.arange(-90, 91, 30))
    ax.grid(True, alpha=0.20)

cbar = fig2.colorbar(im, ax=axes2, orientation='horizontal',
                     fraction=0.05, pad=0.05, extend='both')
cbar.set_label('N:P − 16')
fig2.suptitle('Mean N:P − 16 in 200–1000 m (shared scale, clipped at 98th percentile)',
              fontsize=14)
