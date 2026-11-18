# -----------------------------------------------------------------------------
# Figure 1: nitrate versus phosphate scatter for 200-1000 m
# -----------------------------------------------------------------------------
fig1, axes1 = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
axes1 = axes1.ravel()

all_x = np.concatenate([scatter_data[s][0] for s in SOURCES])
all_y = np.concatenate([scatter_data[s][1] for s in SOURCES])

xmax = np.nanpercentile(all_x, 99.8) * 1.05
ymax = np.nanpercentile(all_y, 99.8) * 1.05
xmax = max(float(xmax), 3.0)
ymax = max(float(ymax), 16.0 * xmax * 1.02)

xticks = np.linspace(0.0, xmax, 6)
yticks = np.linspace(0.0, ymax, 7)

for ax, name in zip(axes1, SOURCES):
    x, y = scatter_data[name]
    row = metrics_df.loc[metrics_df['source'] == name].iloc[0]

    ax.scatter(x, y, s=1, alpha=0.25, color='0.35', rasterized=True)

    xx = np.array([0.0, xmax])
    ax.plot(xx, row['slope'] * xx + row['intercept'],
            color='tab:red', lw=2.0, label='OLS fit')
    ax.plot(xx, 16.0 * xx, color='k', ls='--', lw=1.5,
            label='Redfield 16:1')

    ax.set_xlim(0.0, xmax)
    ax.set_ylim(0.0, ymax)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_title(name)
    ax.set_xlabel(r'PO$_4$ ($\mu$mol kg$^{-1}$)')
    ax.set_ylabel(r'NO$_3$ ($\mu$mol kg$^{-1}$)')

    ax.text(
        0.03, 0.97,
        f"a = {row['slope']:.2f}\nb = {row['intercept']:.2f}\n"
        f"$R^2$ = {row['r2']:.3f}",
        transform=ax.transAxes, va='top', fontsize=9
    )

    if name == SOURCES[0]:
        ax.legend(loc='lower right', fontsize=8)

# -----------------------------------------------------------------------------
# Figure 2: maps of 200-1000 m mean N:P minus 16
# -----------------------------------------------------------------------------
def regrid_np_minus_16(fld):
    lat = np.asarray(fld['lat'], dtype=float)
    lon = np.asarray(fld['lon'], dtype=float)
    ratio = np.asarray(fld['ratio'], dtype=float)
    vol = np.asarray(fld['vol'], dtype=float)

    lon_w = ((lon + 180.0) % 360.0) - 180.0
    i = np.floor((lat + 90.0) / 1.0).astype(int)
    j = np.floor((lon_w + 180.0) / 1.0).astype(int)

    valid = (
        np.isfinite(ratio) & np.isfinite(vol) & (vol > 0.0) &
        (i >= 0) & (i < 180) & (j >= 0) & (j < 360)
    )

    idx = i[valid] * 360 + j[valid]
    w = vol[valid]
    r = ratio[valid]

    sum_w = np.bincount(idx, weights=w, minlength=180 * 360)
    sum_r = np.bincount(idx, weights=r * w, minlength=180 * 360)

    grid = np.divide(
        sum_r, sum_w,
        out=np.full(180 * 360, np.nan),
        where=sum_w > 0
    ).reshape(180, 360)

    return grid - 16.0

grids = {name: regrid_np_minus_16(fields[name]) for name in SOURCES}

all_vals = np.concatenate([g[np.isfinite(g)].ravel() for g in grids.values()])
vmax = np.nanpercentile(np.abs(all_vals), 99.0) if len(all_vals) > 0 else 1.0
vmax = max(float(vmax), 1.0)

fig2, axes2 = plt.subplots(2, 2, figsize=(12, 7), constrained_layout=True)
axes2 = axes2.ravel()

lon_edges = np.arange(-180.0, 181.0, 1.0)
lat_edges = np.arange(-90.0, 91.0, 1.0)

im = None
for ax, name in zip(axes2, SOURCES):
    im = ax.pcolormesh(
        lon_edges, lat_edges, grids[name],
        cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='auto'
    )
    ax.set_xlim(-180.0, 180.0)
    ax.set_ylim(-90.0, 90.0)
    ax.set_aspect('equal')
    ax.set_xticks([-180, -120, -60, 0, 60, 120, 180])
    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
    ax.set_title(name)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

fig2.colorbar(
    im, ax=axes2.tolist(), orientation='horizontal',
    fraction=0.04, pad=0.06, label='Mean N:P - 16 (200-1000 m)'
)
