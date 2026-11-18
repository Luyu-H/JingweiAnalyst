import matplotlib.pyplot as plt
import numpy as np

months = np.arange(1, 13)
months_abbr = [
    'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'
]

colors = {
    'ERSSTv5': 'black',
    'CESM2': 'tab:blue',
    'GFDL-ESM4': 'tab:orange',
    'MPI-ESM1-2-LR': 'tab:green',
}
linestyles = {
    'ERSSTv5': '-',
    'CESM2': '--',
    'GFDL-ESM4': '-.',
    'MPI-ESM1-2-LR': ':',
}
linewidths = {
    'ERSSTv5': 3.0,
    'CESM2': 1.8,
    'GFDL-ESM4': 1.8,
    'MPI-ESM1-2-LR': 1.8,
}

# Figure 1: seasonal amplitude cycle.
fig1, ax1 = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
for src in source_order:
    vals = month_std_dict[src]
    ax1.plot(
        months, vals,
        color=colors[src],
        linestyle=linestyles[src],
        linewidth=linewidths[src],
        marker='o',
        markersize=4,
        label=src,
    )
    peak_i = int(np.argmax(vals))
    ax1.scatter(
        months[peak_i], vals[peak_i],
        color=colors[src],
        s=110,
        marker='o',
        edgecolor='white',
        linewidth=1.0,
        zorder=5,
    )

ax1.set_xticks(months)
ax1.set_xticklabels(months_abbr)
ax1.set_xlabel('Calendar month')
ax1.set_ylabel('Std dev of Nino 3.4 anomaly (°C)')
ax1.set_title('Seasonal cycle of Nino 3.4 amplitude (1950–2014)')
ax1.grid(True, alpha=0.3)
ax1.legend(loc='best', frameon=True)

# Figure 2: four stacked anomaly time-series panels.
fig2, axes = plt.subplots(
    4, 1, figsize=(12, 10), sharex=True, constrained_layout=True
)

ylim_shared = (
    min(-1.05, ymin - 0.10),
    max(1.05, ymax + 0.10),
)

for ax, src in zip(axes, source_order):
    da = idx_dict[src]
    ax.plot(
        da.time, da.values,
        color=colors[src],
        linewidth=1.0,
    )
    ax.axhline(1.0, color='gray', linestyle='--', linewidth=0.8)
    ax.axhline(-1.0, color='gray', linestyle='--', linewidth=0.8)
    ax.axhline(0.0, color='gray', linestyle='-', linewidth=0.5, alpha=0.5)
    ax.set_ylim(ylim_shared)
    ax.set_ylabel('Anomaly (°C)')
    ax.grid(True, alpha=0.25)

    s = stats[src]
    ax.set_title(
        f"{src}   std = {s['std_all']:.2f} °C, "
        f"|anom| > 1: {s['exceed_pct']:.1f}% "
        f"({s['exceed_count']}/{s['n_months']})",
        loc='left',
        fontsize=10,
    )

axes[-1].set_xlabel('Year')
fig2.suptitle('Nino 3.4 monthly anomaly time series (1950–2014)', fontsize=13)
