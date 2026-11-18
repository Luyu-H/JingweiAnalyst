import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------------ #
# Figure 1 - historical OMZ thickness at the 60 threshold
# ------------------------------------------------------------------ #
panels1 = [
    ('GFDL-ESM4  (1971-2000)',    th_gfdl_h_1d),
    ('MPI-ESM1-2-LR  (1971-2000)', th_mpi_h_1d),
    ('WOA23  (1971-2000)',        th_woa_1d),
]

pooled = np.concatenate([np.asarray(p[1].values).ravel() for p in panels1])
vmax1 = float(np.nanpercentile(pooled, 98))

fig1, axes1 = plt.subplots(1, 3, figsize=(16, 4.2), constrained_layout=True)
im1 = None
for ax, (title, fld) in zip(axes1, panels1):
    im1 = ax.pcolormesh(fld['lon'].values, fld['lat'].values,
                        np.asarray(fld.values), cmap='viridis',
                        vmin=0.0, vmax=vmax1, shading='auto')
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_yticks([-90, -45, 0, 45, 90])

cb1 = fig1.colorbar(im1, ax=axes1.tolist(), orientation='horizontal',
                    fraction=0.05, pad=0.10, extend='max')
cb1.set_label('OMZ thickness at O$_2$ < 60 $\\mu$mol kg$^{-1}$ (m)')
fig1.suptitle('Historical oxygen-minimum-zone thickness (shared scale, '
              f'clipped at 98th percentile = {vmax1:.0f} m)', fontsize=11)

# ------------------------------------------------------------------ #
# Figure 2 - ssp245 minus historical (GFDL, MPI) and CESM2 ssp245
# ------------------------------------------------------------------ #
ch_vals = np.concatenate([np.abs(np.asarray(ch_gfdl.values)).ravel(),
                          np.abs(np.asarray(ch_mpi.values)).ravel()])
vmax2 = float(np.nanpercentile(ch_vals, 98))
if vmax2 == 0 or not np.isfinite(vmax2):
    vmax2 = 1.0

cesm_vals = np.asarray(th_cesm_s_1d.values).ravel()
vmax3 = float(np.nanpercentile(cesm_vals, 98))

fig2, axes2 = plt.subplots(1, 3, figsize=(16, 4.2), constrained_layout=True)

im2a = axes2[0].pcolormesh(ch_gfdl['lon'].values, ch_gfdl['lat'].values,
                           np.asarray(ch_gfdl.values), cmap='RdBu_r',
                           vmin=-vmax2, vmax=vmax2, shading='auto')
axes2[0].set_title('GFDL-ESM4: ssp245 (2071-2100) $-$ historical (1971-2000)',
                   fontsize=10)

im2b = axes2[1].pcolormesh(ch_mpi['lon'].values, ch_mpi['lat'].values,
                           np.asarray(ch_mpi.values), cmap='RdBu_r',
                           vmin=-vmax2, vmax=vmax2, shading='auto')
axes2[1].set_title('MPI-ESM1-2-LR: ssp245 (2071-2100) $-$ historical (1971-2000)',
                   fontsize=10)

im2c = axes2[2].pcolormesh(th_cesm_s_1d['lon'].values, th_cesm_s_1d['lat'].values,
                           np.asarray(th_cesm_s_1d.values), cmap='viridis',
                           vmin=0.0, vmax=vmax3, shading='auto')
axes2[2].set_title('CESM2: ssp245 (2071-2100) only\n'
                   '(no historical O$_2$ - no change can be formed)', fontsize=10)

for ax in axes2:
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_yticks([-90, -45, 0, 45, 90])

cb2 = fig2.colorbar(im2a, ax=axes2[:2].tolist(), orientation='horizontal',
                    fraction=0.05, pad=0.10, extend='both')
cb2.set_label('$\\Delta$ OMZ thickness (m)')

cb3 = fig2.colorbar(im2c, ax=axes2[2], orientation='horizontal',
                    fraction=0.16, pad=0.33, extend='max')
cb3.set_label('OMZ thickness (m)')

fig2.suptitle('Projected change in oxygen-minimum-zone thickness at O$_2$ < 60 '
              f'$\\mu$mol kg$^{{-1}}$ (change clipped at 98th pct = {vmax2:.0f} m)',
              fontsize=11)

plt.show()
