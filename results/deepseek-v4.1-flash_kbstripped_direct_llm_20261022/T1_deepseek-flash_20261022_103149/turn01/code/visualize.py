# -----------------------------------------------------------------------------
# Figure 1: seasonal cycle of Nino 3.4 amplitude
# -----------------------------------------------------------------------------
colors = {
    'ERSSTv5': 'black',
    'CESM2': '#1f77b4',
    'GFDL-ESM4': '#d62728',
    'MPI-ESM1-2-LR': '#2ca02c',
}
linewidths = {
    'ERSSTv5': 3.0,
    'CESM2': 2.0,
    'GFDL-ESM4': 2.0,
    'MPI-ESM1-2-LR': 2.0,
}

fig1, ax1 = plt.subplots(figsize=(9, 5.5))

for name in source_order:
    ax1.plot(
        months,
        monthly_std[name],
        marker='o',
        ms=4,
        lw=linewidths[name],
        color=colors[name],
        label=name,
    )
    ax1.scatter(
        [peak_months[name]],
        [peak_values[name]],
        s=90,
        facecolor=colors[name],
        edgecolor='white',
        zorder=5,
    )
    label = 'OBS' if name == 'ERSSTv5' else name
    ax1.annotate(
        f'{label} peak: month {peak_months[name]}',
        xy=(peak_months[name], peak_values[name]),
        xytext=(0, 8),
        textcoords='offset points',
        ha='center',
        fontsize=8,
        color=colors[name],
    )

ax1.set_xticks(months)
ax1.set_xticklabels(
    ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
     'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
)
ax1.set_xlabel('Calendar month')
ax1.set_ylabel('Nino 3.4 SST anomaly std. dev. (°C)')
ax1.set_title('Seasonal cycle of Nino 3.4 amplitude, 1950–2014')
ax1.grid(True, alpha=0.3)
ax1.legend(frameon=False, ncol=4, loc='upper center')
fig1.tight_layout()

# -----------------------------------------------------------------------------
# Figure 2: Nino 3.4 anomaly time series in four stacked panels
# -----------------------------------------------------------------------------
fig2, axes = plt.subplots(4, 1, figsize=(13, 9), sharex=True, sharey=True)

all_min = min(float(a.min()) for a in anomalies.values())
all_max = max(float(a.max()) for a in anomalies.values())
span = max(all_max - all_min, 1.0)
ylim = (all_min - 0.08 * span, all_max + 0.08 * span)

for ax, name in zip(axes, source_order):
    a = anomalies[name]
    x = a.index.to_timestamp()

    ax.plot(x, a.values, color=colors[name], lw=1.0)
    ax.axhline(0.0, color='0.4', lw=0.7)
    ax.axhline(1.0, color='0.2', ls='--', lw=0.9)
    ax.axhline(-1.0, color='0.2', ls='--', lw=0.9)

    ax.set_ylim(ylim)
    ax.grid(True, alpha=0.25)
    ax.set_ylabel('Anomaly (°C)')
    ax.set_title(name, loc='left', fontsize=11, fontweight='bold')

    txt = (
        f"σ = {std_values[name]:.2f} °C\n"
        f"|anom| > 1.0 °C: {exceed_counts[name]} months "
        f"({exceed_pct_values[name]:.1f}%)"
    )
    ax.text(
        0.01,
        0.95,
        txt,
        transform=ax.transAxes,
        va='top',
        ha='left',
        fontsize=9,
        bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='0.7', alpha=0.85),
    )

axes[-1].set_xlim(pd.Timestamp('1950-01-01'), pd.Timestamp('2015-01-01'))
axes[-1].set_xlabel('Year')
fig2.suptitle('Nino 3.4 SST anomalies, 1950–2014', y=0.995)
fig2.tight_layout()
