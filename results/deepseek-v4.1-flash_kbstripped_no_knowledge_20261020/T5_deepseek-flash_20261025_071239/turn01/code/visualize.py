# ---------------------------------------------------------------------------
# Consistent model colours
# ---------------------------------------------------------------------------
model_colors = {
    'CESM2': 'C0',
    'GFDL-ESM4': 'C1',
    'MPI-ESM1-2-LR': 'C2',
}

# ---------------------------------------------------------------------------
# Figure 1: 1971-2000 SSS bias maps (model minus WOA23)
# ---------------------------------------------------------------------------
fig1, axes1 = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)
bias_pool = np.concatenate(
    [bias_maps[model].values.ravel() for model in models]
)
vmax_bias = float(np.nanpercentile(np.abs(bias_pool), 98))
if not np.isfinite(vmax_bias) or vmax_bias == 0:
    vmax_bias = 0.1

mappable1 = None
for ax, model in zip(axes1, models):
    field = bias_maps[model]
    lon = field.lon.values
    lat = field.lat.values
    mappable1 = ax.pcolormesh(
        lon, lat, field.values,
        cmap='RdBu_r', vmin=-vmax_bias, vmax=vmax_bias,
        shading='auto'
    )
    ax.set_title(model, color=model_colors[model], fontsize=11)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

cbar1 = fig1.colorbar(
    mappable1, ax=axes1.ravel().tolist(),
    label='SSS bias (psu)', extend='both'
)
fig1.suptitle(
    '1971-2000 SSS bias (model minus WOA23); shared scale clipped at 98th percentile',
    fontsize=12
)

# ---------------------------------------------------------------------------
# Figure 2: 1950-2014 cell-by-cell linear SSS trend maps
# ---------------------------------------------------------------------------
fig2, axes2 = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)
trend_pool = np.concatenate(
    [trend_slope_reg[model].values.ravel() for model in models]
)
vmax_trend = float(np.nanpercentile(np.abs(trend_pool), 98))
if not np.isfinite(vmax_trend) or vmax_trend == 0:
    vmax_trend = 0.1

mappable2 = None
for ax, model in zip(axes2, models):
    field = trend_slope_reg[model]
    sig = trend_sig_reg[model]
    lon = field.lon.values
    lat = field.lat.values

    mappable2 = ax.pcolormesh(
        lon, lat, field.values,
        cmap='RdBu_r', vmin=-vmax_trend, vmax=vmax_trend,
        shading='auto'
    )
    # p < 0.05 stippling
    ax.contourf(
        lon, lat, sig.astype(float).values,
        levels=[0.5, 1.5], colors='none', hatches=['...']
    )
    ax.set_title(model, color=model_colors[model], fontsize=11)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

cbar2 = fig2.colorbar(
    mappable2, ax=axes2.ravel().tolist(),
    label='SSS trend (psu/decade)', extend='both'
)
fig2.suptitle(
    '1950-2014 annual-mean SSS linear trend; p < 0.05 stippled; '
    'shared scale clipped at 98th percentile',
    fontsize=12
)

# ---------------------------------------------------------------------------
# Figure 3: global-mean SSS annual series, historical and ssp245 segments
# ---------------------------------------------------------------------------
def series_xy(da):
    if 'year' in da.coords:
        x = da['year'].values
    elif 'time' in da.coords:
        x = da['time'].dt.year.values
    else:
        x = np.arange(da.sizes[list(da.dims)[0]])
    return x, da.values


fig3, ax3 = plt.subplots(figsize=(9, 5), constrained_layout=True)

for model in models:
    color = model_colors[model]

    # Historical segment
    xh, yh = series_xy(hist_annual_global[model])
    ax3.plot(xh, yh, color=color, lw=1.8,
             label=f'{model} historical')
    xhf, yhf = series_xy(hist_fit[model]['fitted_data'])
    ax3.plot(xhf, yhf, color=color, lw=1.0, ls='--')

    # SSP2-4.5 segment
    xs, ys = series_xy(ssp_annual_global[model])
    ax3.plot(xs, ys, color=color, lw=1.8, ls='-',
             label=f'{model} ssp245')
    xsf, ysf = series_xy(ssp_fit[model]['fitted_data'])
    ax3.plot(xsf, ysf, color=color, lw=1.0, ls=':')

# Leave the 2015-2070 gap visibly empty
ax3.axvspan(2015, 2070, color='0.9', alpha=0.25, zorder=0)

ax3.set_xlabel('Year')
ax3.set_ylabel('Global mean SSS (psu)')
ax3.set_title('Global mean sea surface salinity: historical and ssp245')
ax3.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)
ax3.grid(True, alpha=0.3)
ax3.set_xlim(1948, 2102)
