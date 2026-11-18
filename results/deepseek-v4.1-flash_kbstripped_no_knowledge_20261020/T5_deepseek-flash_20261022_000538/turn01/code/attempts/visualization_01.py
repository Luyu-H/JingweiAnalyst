# ---------------------------------------------------------------------------
# Shared model colours
# ---------------------------------------------------------------------------
colors = {
    'CESM2': 'C0',
    'GFDL-ESM4': 'C1',
    'MPI-ESM1-2-LR': 'C2',
}


def _x(da):
    """Return the x coordinate for a 1-D time series DataArray."""
    dim = da.dims[0]
    if dim in da.coords:
        return da[dim].values
    if 'year' in da.coords:
        return da['year'].values
    if 'time' in da.coords:
        return da['time'].dt.year.values
    return np.arange(da.sizes[dim])


# ---------------------------------------------------------------------------
# Figure 1: SSS bias (model 1971-2000 climatology minus WOA23)
# ---------------------------------------------------------------------------
bias_arrays = [bias_fields[m].values.ravel() for m in models]
vmax_bias = float(np.nanpercentile(np.abs(np.concatenate(bias_arrays)), 98))
if not np.isfinite(vmax_bias) or vmax_bias == 0:
    vmax_bias = 0.1

fig, axes = plt.subplots(1, 3, figsize=(16, 4.5), constrained_layout=True)
mappable_bias = None
for ax, model in zip(axes, models):
    b = bias_fields[model]
    lon = b['lon'].values
    lat = b['lat'].values
    mappable_bias = ax.pcolormesh(
        lon, lat, b.values,
        cmap='RdBu_r',
        vmin=-vmax_bias,
        vmax=vmax_bias,
        shading='auto',
    )
    ax.set_title(model)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

fig.colorbar(
    mappable_bias,
    ax=axes,
    shrink=0.8,
    label='SSS bias (psu)',
    extend='both',
)
fig.suptitle(
    f'SSS bias, 1971-2000 (model - WOA23); shared scale clipped at '
    f'{vmax_bias:.2f} psu (98th percentile)',
    fontsize=12,
)

# ---------------------------------------------------------------------------
# Figure 2: cell-by-cell linear trend of annual-mean SSS, 1950-2014
# ---------------------------------------------------------------------------
trend_arrays = [trend_fields[m].values.ravel() for m in models]
vmax_trend = float(np.nanpercentile(np.abs(np.concatenate(trend_arrays)), 98))
if not np.isfinite(vmax_trend) or vmax_trend == 0:
    vmax_trend = 0.01

fig, axes = plt.subplots(1, 3, figsize=(16, 4.5), constrained_layout=True)
mappable_trend = None
for ax, model in zip(axes, models):
    t = trend_fields[model]
    s = sig_fields[model]
    lon = t['lon'].values
    lat = t['lat'].values

    mappable_trend = ax.pcolormesh(
        lon, lat, t.values,
        cmap='RdBu_r',
        vmin=-vmax_trend,
        vmax=vmax_trend,
        shading='auto',
    )
    # Stipple cells significant at p < 0.05
    ax.contourf(
        lon, lat, s.values.astype(float),
        levels=[0.5, 1.5],
        colors='none',
        hatches=['..'],
        add_colorbar=False,
    )
    ax.set_title(model)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

fig.colorbar(
    mappable_trend,
    ax=axes,
    shrink=0.8,
    label='SSS trend (psu/decade)',
    extend='both',
)
fig.suptitle(
    f'Annual-mean SSS trend, 1950-2014; stippled where p < 0.05; '
    f'shared scale clipped at {vmax_trend:.3f} psu/decade',
    fontsize=12,
)

# ---------------------------------------------------------------------------
# Figure 3: global mean SSS, historical and ssp245 segments
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(12, 5.5), constrained_layout=True)

for model in models:
    h = hist_series[model]
    s = ssp_series[model]

    ax.plot(
        _x(h), h.values,
        color=colors[model],
        linestyle='-',
        linewidth=1.6,
        label=f'{model} historical',
    )
    ax.plot(
        _x(s), s.values,
        color=colors[model],
        linestyle='-',
        linewidth=1.6,
        label=f'{model} ssp245',
    )

    hf = hist_fits[model]['fitted_data']
    sf = ssp_fits[model]['fitted_data']
    ax.plot(
        _x(hf), hf.values,
        color=colors[model],
        linestyle='--',
        linewidth=1.0,
    )
    ax.plot(
        _x(sf), sf.values,
        color=colors[model],
        linestyle='--',
        linewidth=1.0,
    )

ax.set_xlabel('Year')
ax.set_ylabel('Global mean SSS (psu)')
ax.set_xlim(1948, 2102)
ax.grid(True, alpha=0.3)
ax.legend(
    loc='upper center',
    bbox_to_anchor=(0.5, -0.12),
    ncol=3,
    fontsize=8,
    frameon=False,
)
ax.set_title('Global mean SSS: historical 1950-2014 and ssp245 2071-2100')
