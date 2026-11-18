# ---------------------------------------------------------------------------
# Figure 1: maps of operator MLD minus published mlotst
# ---------------------------------------------------------------------------
fig1, axes1 = plt.subplots(3, 2, figsize=(14, 10), constrained_layout=True)

pooled = np.concatenate([np.abs(map_fields[(m, mn)].values).ravel()
                         for m in models for mn in month_names])
vmax = float(np.nanpercentile(pooled[np.isfinite(pooled)], 98))

mesh = None
for i, model in enumerate(models):
    for j, mname in enumerate(month_names):
        ax = axes1[i, j]
        d = map_fields[(model, mname)]
        mesh = d.plot(ax=ax, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                      add_colorbar=False, extend='both')
        ax.set_title(f'{model} — {mname}', fontsize=10)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xlim(float(d['lon'].min()), float(d['lon'].max()))
        ax.set_ylim(float(d['lat'].min()), float(d['lat'].max()))
        ax.tick_params(labelsize=8)

for ax in axes1[:, 0]:
    ax.set_ylabel('Latitude', fontsize=9)
for ax in axes1[-1, :]:
    ax.set_xlabel('Longitude', fontsize=9)

cb1 = fig1.colorbar(mesh, ax=axes1, shrink=0.6, pad=0.02,
                    label='Operator MLD - mlotst (m)')
fig1.suptitle('Figure 1. Operator MLD minus published mlotst, 1971-2000 climatologies '
              '(shared diverging scale, clipped at pooled 98th percentile '
              f'|$\\Delta$| = {vmax:.0f} m)', fontsize=12)

# ---------------------------------------------------------------------------
# Figure 2: native-grid scatter, operator MLD vs published mlotst
# ---------------------------------------------------------------------------
pos_all = []
for key, (opv, pubv) in scatter_vals.items():
    good = np.isfinite(opv) & np.isfinite(pubv) & (opv > 0) & (pubv > 0)
    pos_all.append(opv[good])
    pos_all.append(pubv[good])
pos_all = np.concatenate(pos_all)

lo = max(0.5, 10.0 ** np.floor(np.log10(np.nanmin(pos_all))))
hi = 10.0 ** np.ceil(np.log10(np.nanmax(pos_all)))
ticks = [t for t in [0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
         if lo <= t <= hi]

fig2, axes2 = plt.subplots(3, 2, figsize=(12, 12), constrained_layout=True)
for i, model in enumerate(models):
    for j, mname in enumerate(month_names):
        ax = axes2[i, j]
        opv, pubv = scatter_vals[(model, mname)]
        good = np.isfinite(opv) & np.isfinite(pubv) & (opv > 0) & (pubv > 0)
        ax.scatter(opv[good], pubv[good], s=0.6, alpha=0.15, c='k', linewidths=0)
        ax.plot([lo, hi], [lo, hi], 'r--', lw=1.0, label='y = x')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels([f'{t:g}' for t in ticks], fontsize=7, rotation=45)
        ax.set_yticklabels([f'{t:g}' for t in ticks], fontsize=7)
        ax.grid(True, which='both', alpha=0.2, lw=0.4)
        ax.set_title(f'{model} — {mname}', fontsize=10)

        d = results[(model, mname)]
        txt = (f'bias = {d["bias"]:+.1f} m\n'
               f'RMSE = {d["rmse"]:.1f} m\n'
               f'r = {d["corr"]:+.3f}')
        ax.text(0.03, 0.97, txt, transform=ax.transAxes, va='top', ha='left',
                fontsize=7.5, bbox=dict(facecolor='white', alpha=0.85, edgecolor='none'))

for ax in axes2[:, 0]:
    ax.set_ylabel('Published mlotst (m)', fontsize=9)
for ax in axes2[-1, :]:
    ax.set_xlabel('Operator MLD (m)', fontsize=9)

fig2.suptitle('Figure 2. Operator MLD vs published mlotst on the common valid native-grid '
              'cells (log-log, identical range)', fontsize=12)

plt.show()
