import matplotlib.pyplot as plt
import numpy as np

models = ['CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']
months = [3, 9]
month_names = {3: 'March', 9: 'September'}

# ----------------------------------------------------------------------
# Figure 1: maps of operator MLD minus mlotst
# ----------------------------------------------------------------------
all_abs = []
for m in models:
    for mo in months:
        res = results[(m, mo)]
        mask = res['mask'] & np.isfinite(res['diff'])
        all_abs.append(np.abs(res['diff'][mask]))
all_abs = np.concatenate(all_abs)
vmax = np.nanpercentile(all_abs, 98)
vmax = max(vmax, 1.0)

fig, axes = plt.subplots(3, 2, figsize=(14, 12),
                         sharex=True, sharey=True)
for i, m in enumerate(models):
    for j, mo in enumerate(months):
        ax = axes[i, j]
        res = results[(m, mo)]
        lat2d, lon2d = get_latlon_2d(res['lat'], res['lon'])
        glon, glat, gdiff = regrid_to_regular(
            lat2d, lon2d, res['diff'], dlat=2, dlon=2
        )
        pc = ax.pcolormesh(glon, glat, gdiff, cmap='RdBu_r',
                           vmin=-vmax, vmax=vmax, shading='auto')
        ax.set_title(f'{m} {month_names[mo]}')
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        ax.set_xticks(np.arange(-180, 181, 60))
        ax.set_yticks(np.arange(-90, 91, 30))
        ax.grid(True, alpha=0.3)
        if j == 0:
            ax.set_ylabel('Latitude')
        if i == 2:
            ax.set_xlabel('Longitude')

fig.colorbar(pc, ax=axes, orientation='horizontal',
             fraction=0.05, pad=0.05,
             label='Operator MLD - mlotst (m)')
fig.suptitle('Operator MLD minus published mlotst', y=0.92)

# ----------------------------------------------------------------------
# Figure 2: scatter plots, operator MLD vs mlotst
# ----------------------------------------------------------------------
all_vals = []
for m in models:
    for mo in months:
        res = results[(m, mo)]
        mask = (res['mask']
                & (res['op_mean'] > 0)
                & (res['ml_mean'] > 0))
        all_vals.append(res['op_mean'][mask])
        all_vals.append(res['ml_mean'][mask])
all_vals = np.concatenate(all_vals)
all_vals = all_vals[np.isfinite(all_vals) & (all_vals > 0)]

vmin = np.nanmin(all_vals)
vmax2 = np.nanmax(all_vals)
vmin = 10 ** np.floor(np.log10(vmin))
vmax2 = 10 ** np.ceil(np.log10(vmax2))

ticks = []
for exp in range(int(np.log10(vmin)), int(np.log10(vmax2)) + 1):
    for mult in [1, 2, 5]:
        t = mult * 10**exp
        if vmin <= t <= vmax2:
            ticks.append(t)
ticks = sorted(set(ticks))

fig2, axes2 = plt.subplots(3, 2, figsize=(12, 12),
                           sharex=True, sharey=True)
for i, m in enumerate(models):
    for j, mo in enumerate(months):
        ax = axes2[i, j]
        res = results[(m, mo)]
        mask = (res['mask']
                & (res['op_mean'] > 0)
                & (res['ml_mean'] > 0))
        x = res['op_mean'][mask]
        y = res['ml_mean'][mask]

        ax.scatter(x, y, s=1, alpha=0.3, c='k', edgecolors='none')
        ax.plot([vmin, vmax2], [vmin, vmax2], 'k--', lw=0.8)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(vmin, vmax2)
        ax.set_ylim(vmin, vmax2)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels([f'{t:g}' for t in ticks])
        ax.set_yticklabels([f'{t:g}' for t in ticks])

        ax.set_title(f'{m} {month_names[mo]}')
        ax.text(0.05, 0.95,
                f'bias = {res["bias"]:.1f} m\n'
                f'RMSE = {res["rmse"]:.1f} m\n'
                f'r = {res["corr"]:.3f}',
                transform=ax.transAxes, va='top', fontsize=9)
        ax.grid(True, which='both', alpha=0.3)

        if j == 0:
            ax.set_ylabel('mlotst (m)')
        if i == 2:
            ax.set_xlabel('Operator MLD (m)')

fig2.suptitle('Operator MLD vs published mlotst', y=0.92)
