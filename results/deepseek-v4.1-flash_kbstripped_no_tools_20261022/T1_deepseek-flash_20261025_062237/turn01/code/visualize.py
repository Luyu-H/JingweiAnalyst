colors = {
    'ERSSTv5': 'black',
    'CESM2': 'tab:blue',
    'GFDL-ESM4': 'tab:orange',
    'MPI-ESM1-2-LR': 'tab:green',
}
linestyles = {
    'ERSSTv5': '-',
    'CESM2': '--',
    'GFDL-ESM4': '--',
    'MPI-ESM1-2-LR': '--',
}
markers = {
    'ERSSTv5': 'o',
    'CESM2': 's',
    'GFDL-ESM4': '^',
    'MPI-ESM1-2-LR': 'D',
}

# Figure 1: seasonal cycle of Nino 3.4 amplitude
fig1, ax1 = plt.subplots(figsize=(10, 6), constrained_layout=True)
xmonths = np.arange(1, 13)

for name in sources_order:
    r = results[name]
    lw = 3.0 if name == 'ERSSTv5' else 1.8
    ax1.plot(xmonths, r['monthly_sd'],
             linestyle=linestyles[name], marker=markers[name],
             markersize=5, linewidth=lw, color=colors[name],
             label=name, zorder=3)
    ax1.scatter(r['peak_month'], r['peak_sd'],
                marker='*', s=170, color=colors[name],
                edgecolor='black', linewidth=0.6, zorder=5)
    ax1.annotate(month_names_abbr[r['peak_month'] - 1],
                 (r['peak_month'], r['peak_sd']),
                 textcoords='offset points', xytext=(0, 8),
                 ha='center', fontsize=8, color=colors[name], zorder=6)

ax1.set_xticks(xmonths)
ax1.set_xticklabels(month_names_full, rotation=45, ha='right')
ax1.set_xlabel('Calendar month')
ax1.set_ylabel('Nino 3.4 anomaly standard deviation (degC)')
ax1.set_title('Seasonal cycle of Nino 3.4 amplitude, 1950-2014\n'
              '(monthly anomaly SD by calendar month; stars mark each peak)')
ax1.grid(True, alpha=0.3)
ax1.legend(loc='best')
ax1.set_ylim(bottom=0)

# Figure 2: stacked anomaly time series
max_abs = max(np.max(np.abs(results[name]['anom'])) for name in sources_order)
ylim = max_abs * 1.08

fig2, axes = plt.subplots(4, 1, figsize=(14, 10), sharex=True,
                          constrained_layout=True)

for ax, name in zip(axes, sources_order):
    r = results[name]
    ax.plot(r['time_years'], r['anom'],
            color=colors[name], linewidth=1.0)
    ax.axhline(1.0, color='red', linestyle='--', linewidth=1.0, alpha=0.8)
    ax.axhline(-1.0, color='red', linestyle='--', linewidth=1.0, alpha=0.8)
    ax.axhline(0.0, color='gray', linestyle=':', linewidth=0.8, alpha=0.6)
    ax.set_ylim(-ylim, ylim)
    ax.set_xlim(1950, 2015)
    ax.set_ylabel('Anomaly (degC)')
    ax.grid(True, alpha=0.25)
    ax.text(0.01, 0.95,
            f"{name}\nSD = {r['sd']:.2f} degC\n"
            f"|anom| > 1.0: {r['exceed_count']} months ({r['exceed_pct']:.1f}%)",
            transform=ax.transAxes, va='top', ha='left', fontsize=9,
            bbox=dict(facecolor='white', alpha=0.85, edgecolor='none', pad=2.0))

axes[-1].set_xlabel('Year')
axes[0].set_title('Nino 3.4 monthly anomaly series, 1950-2014 '
                  '(each source relative to its own 1950-2014 calendar-month climatology)')
