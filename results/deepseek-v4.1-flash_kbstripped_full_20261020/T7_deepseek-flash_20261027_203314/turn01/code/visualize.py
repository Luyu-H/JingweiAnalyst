import numpy as np
import matplotlib.pyplot as plt

SOURCES_ORDER = ['WOA23', 'CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']
COLORS = {'WOA23': 'k', 'CESM2': 'tab:blue', 'GFDL-ESM4': 'tab:green', 'MPI-ESM1-2-LR': 'tab:red'}

# ------------------------------------------------------------------
# Figure 1: NO3 vs PO4 scatter, 200-1000 m, one panel per source
# ------------------------------------------------------------------
fig1, axes1 = plt.subplots(2, 2, figsize=(11, 9.5), constrained_layout=True)

xall = np.concatenate([scatter_data[s][0] for s in SOURCES_ORDER])
yall = np.concatenate([scatter_data[s][1] for s in SOURCES_ORDER])
xmin, xmax = 0.0, float(np.nanpercentile(xall, 99.7))
ymin, ymax = 0.0, float(np.nanpercentile(yall, 99.7))
xticks = np.linspace(0, np.ceil(xmax), 6)
yticks = np.linspace(0, np.ceil(ymax), 6)

xline = np.linspace(0.0, np.ceil(xmax), 100)

for ax, s in zip(axes1.ravel(), SOURCES_ORDER):
    xs, ys = scatter_data[s]
    ax.scatter(xs, ys, s=2, alpha=0.18, color=COLORS[s], edgecolors='none',
               label='cells (200-1000 m)')
    R = results[s]
    ax.plot(xline, R['a'] * xline + R['b'], color='tab:orange', lw=2,
            label=f"fit: NO3 = {R['a']:.2f}*PO4 + {R['b']:.2f}")
    ax.plot(xline, 16.0 * xline, color='m', lw=1.4, ls='--', label='Redfield 16:1')
    ax.set_xlim(0, np.ceil(xmax))
    ax.set_ylim(0, np.ceil(ymax))
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xlabel('Phosphate (umol kg$^{-1}$)')
    ax.set_ylabel('Nitrate (umol kg$^{-1}$)')
    ax.set_title(f"{s}\na = {R['a']:.2f}, b = {R['b']:.2f}, R$^2$ = {R['r2']:.3f}")
    ax.legend(loc='upper left', fontsize=8, framealpha=0.85)
    ax.grid(alpha=0.25)

fig1.suptitle('Nitrate vs phosphate, 200-1000 m (each source, equal-cell-weight OLS)', fontsize=13)

# ------------------------------------------------------------------
# Figure 2: 200-1000 m mean N:P minus 16, common lon-lat frame
# ------------------------------------------------------------------
fig2, axes2 = plt.subplots(2, 2, figsize=(13, 7.5), constrained_layout=True)

diff_maps = {s: (maps_reg[s] - 16.0) for s in SOURCES_ORDER}
pooled = np.concatenate([np.asarray(diff_maps[s].values, dtype=float).ravel()
                         for s in SOURCES_ORDER])
vmax = float(np.nanpercentile(np.abs(pooled), 98))
if not np.isfinite(vmax) or vmax == 0:
    vmax = 1.0

pcm = None
for ax, s in zip(axes2.ravel(), SOURCES_ORDER):
    f = diff_maps[s]
    lon = np.asarray(f['lon'].values)
    lat = np.asarray(f['lat'].values)
    pcm = ax.pcolormesh(lon, lat, np.asarray(f.values), cmap='RdBu_r',
                        vmin=-vmax, vmax=vmax, shading='auto')
    ax.set_xlim(-180, 180)
    ax.set_ylim(-80, 85)
    ax.set_xticks([-180, -120, -60, 0, 60, 120, 180])
    ax.set_yticks([-80, -40, 0, 40, 80])
    ax.set_title(f"{s}", fontsize=11)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

cb = fig2.colorbar(pcm, ax=axes2, orientation='vertical', extend='both', shrink=0.9)
cb.set_label('N:P - 16')
fig2.suptitle(f'Mean N:P minus 16 over 200-1000 m (shared scale, clipped at '
              f'{vmax:.1f}; models regridded to the 1$^\\circ$ grid for display only)',
              fontsize=12)

print("figures built")
