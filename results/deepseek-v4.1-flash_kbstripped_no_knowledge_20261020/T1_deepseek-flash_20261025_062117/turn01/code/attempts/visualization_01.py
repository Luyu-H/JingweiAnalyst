# Visualization code only

colors = {
    'ERSSTv5': 'black',
    'CESM2': 'tab:blue',
    'GFDL-ESM4': 'tab:orange',
    'MPI-ESM1-2-LR': 'tab:green',
}
styles = {
    'ERSSTv5': '-',
    'CESM2': '--',
    'GFDL-ESM4': '-.',
    'MPI-ESM1-2-LR': ':',
}

# Figure 1: seasonal cycle of Nino3.4 amplitude
fig1, ax1 = plt.subplots(figsize=(10, 6), constrained_layout=True)

for src in sources_order:
    ax1.plot(
        months,
        seasonal_std[src].values,
        label=src,
        color=colors[src],
        linestyle=styles[src],
        linewidth=3.0 if src == 'ERSSTv5' else 1.8,
        marker='o' if src == 'ERSSTv5' else None,
        markersize=4,
    )
    ax1.scatter(
        peak_month[src],
        peak_value[src],
        s=90,
        facecolors='none',
        edgecolors=colors[src],
        linewidths=2.0,
        zorder=5,
    )

# Dummy handle to explain the open-circle markers
ax1.scatter([], [], s=90, facecolors='none', edgecolors='k', linewidths=2.0, label='peak month')

ax1.set_xticks(months)
ax1.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                     'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
ax1.set_xlabel('Calendar month')
ax1.set_ylabel('Nino3.4 anomaly standard deviation (°C)')
ax1.set_title('Seasonal cycle of Nino3.4 amplitude, 1950–2014')
ax1.grid(True, alpha=0.3)
ax1.legend(title='Source', loc='best')

# Figure 2: stacked Nino3.4 anomaly time series
fig2, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=True, constrained_layout=True)

all_abs = np.concatenate([np.abs(nino_series[src].values) for src in sources_order])
ymax = float(np.nanmax(all_abs)) + 0.2
ymax = max(ymax, 1.2)

for ax, src in zip(axes, sources_order):
    da = nino_series[src]
    ax.plot(da.time.values, da.values, color=colors[src], linewidth=1.2)
    ax.axhline(1.0, color='red', linestyle='--', linewidth=0.8)
    ax.axhline(-1.0, color='red', linestyle='--', linewidth=0.8)
    ax.axhline(0.0, color='gray', linestyle=':', linewidth=0.6)
    ax.set_ylim(-ymax, ymax)
    ax.set_ylabel('Anomaly (°C)')
    ax.grid(True, alpha=0.3)
    ax.set_title(
        f"{src}: σ = {std_all[src]:.2f} °C; "
        f"|anom| > 1 °C: {exceed_n[src]}/780 months ({exceed_pct[src]:.1f}%)",
        loc='left',
        fontsize=10,
    )

axes[-1].set_xlabel('Year')
axes[-1].set_xlim(pd.Timestamp('1950-01-01'), pd.Timestamp('2015-01-01'))
fig2.suptitle('Nino3.4 SST anomaly time series, 1950–2014', fontsize=14)
