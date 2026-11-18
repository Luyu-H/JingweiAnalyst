months = np.arange(1, 13)
month_labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

colors = {
    'ERSSTv5': 'black',
    'CESM2': 'tab:blue',
    'GFDL-ESM4': 'tab:orange',
    'MPI-ESM1-2-LR': 'tab:green',
}

# Figure 1: seasonal cycle of Nino 3.4 amplitude
fig1, ax1 = plt.subplots(figsize=(9, 5))

for src in SOURCES:
    vals = nino34_stats[src]['seasonal_std']

    if src == 'ERSSTv5':
        ax1.plot(
            months, vals, color=colors[src], lw=3.0, marker='o', ms=5,
            label=src, zorder=10
        )
    else:
        ax1.plot(
            months, vals, color=colors[src], lw=1.6, marker='o', ms=4,
            label=src
        )

    pm = nino34_stats[src]['peak_month']
    pv = nino34_stats[src]['peak_std']
    ax1.scatter(
        pm, pv, color=colors[src],
        s=150 if src == 'ERSSTv5' else 80,
        marker='*',
        edgecolor='black' if src == 'ERSSTv5' else 'none',
        zorder=11
    )

ax1.set_xticks(months)
ax1.set_xticklabels(month_labels)
ax1.set_xlabel('Calendar month')
ax1.set_ylabel('Nino 3.4 anomaly standard deviation (degC)')
ax1.set_title('Seasonal cycle of Nino 3.4 amplitude, 1950-2014')
ax1.grid(True, alpha=0.3)
ax1.legend(frameon=False)
fig1.tight_layout()

# Figure 2: four stacked Nino 3.4 anomaly panels
fig2, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=True)

all_anoms = np.concatenate([nino34_stats[s]['anom_values'] for s in SOURCES])
ymax = np.nanmax(np.abs(all_anoms)) * 1.05

for ax, src in zip(axes, SOURCES):
    vals = nino34_stats[src]['anom_values']

    ax.plot(TIME_AXIS, vals, color=colors[src], lw=1.0)
    ax.axhline(1.0, color='gray', linestyle='--', linewidth=0.8)
    ax.axhline(-1.0, color='gray', linestyle='--', linewidth=0.8)
    ax.axhline(0.0, color='gray', linestyle=':', linewidth=0.5)
    ax.set_ylim(-ymax, ymax)
    ax.set_ylabel('Anomaly (degC)')

    ax.text(
        0.01, 0.95,
        f"{src}\n"
        f"SD = {nino34_stats[src]['std']:.2f} degC\n"
        f"|anom| > 1 degC: {nino34_stats[src]['exceed_count']} months "
        f"({nino34_stats[src]['exceed_pct']:.1f}%)",
        transform=ax.transAxes,
        va='top',
        ha='left',
        fontsize=9,
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
    )

axes[-1].set_xlabel('Year')
fig2.suptitle('Nino 3.4 SST anomalies, 1950-2014')
fig2.tight_layout(rect=[0, 0, 1, 0.97])
