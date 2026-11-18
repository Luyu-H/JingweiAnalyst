import numpy as np
import matplotlib.pyplot as plt

colors = {
    "ERSSTv5": "black",
    "CESM2": "tab:blue",
    "GFDL-ESM4": "tab:orange",
    "MPI-ESM1-2-LR": "tab:green",
}

# Figure 1: seasonal cycle of Nino 3.4 amplitude.
months = np.arange(1, 13)
fig1, ax1 = plt.subplots(figsize=(9.5, 5.2), constrained_layout=True)

for name in SOURCE_ORDER:
    s = results[name]
    vals = s["seasonal_std"]
    lw = 3.0 if name == "ERSSTv5" else 1.8
    ax1.plot(
        months, vals,
        marker="o", lw=lw, ms=4.5,
        color=colors[name], label=name,
    )
    peak_idx = s["peak_month"] - 1
    ax1.plot(
        months[peak_idx], vals[peak_idx],
        marker="*", ms=15, color=colors[name],
        markeredgecolor="black", markeredgewidth=0.6, zorder=5,
    )
    ax1.annotate(
        MONTH_ABBR[peak_idx],
        (months[peak_idx], vals[peak_idx]),
        xytext=(0, 8), textcoords="offset points",
        ha="center", fontsize=8,
    )

ax1.set_xticks(months)
ax1.set_xticklabels(MONTH_ABBR)
ax1.set_xlabel("Calendar month")
ax1.set_ylabel("Standard deviation of monthly Nino 3.4 anomaly (degC)")
ax1.set_title("Seasonal amplitude of Nino 3.4, 1950-2014")
ax1.grid(alpha=0.25)
ax1.legend(
    loc="upper left", bbox_to_anchor=(1.02, 1.0),
    frameon=False, borderaxespad=0.0,
)

# Figure 2: monthly Nino 3.4 anomaly series, stacked.
anomalies = [results[name]["anomaly"] for name in SOURCE_ORDER]
ymax = max(float(np.nanmax(np.abs(a))) for a in anomalies)
ylim = (-1.05 * ymax, 1.05 * ymax)

fig2, axes = plt.subplots(
    4, 1, figsize=(12.5, 8.0), sharex=True, constrained_layout=True
)

for ax, name in zip(axes, SOURCE_ORDER):
    s = results[name]
    ax.plot(target_time, s["anomaly"], lw=0.8, color=colors[name])
    ax.axhline(1.0, color="black", lw=0.7, ls="--")
    ax.axhline(-1.0, color="black", lw=0.7, ls="--")
    ax.axhline(0.0, color="0.6", lw=0.5)
    ax.set_ylim(ylim)
    ax.set_ylabel("Anomaly (degC)")
    ax.grid(alpha=0.20)
    ax.set_title(
        f"{name}: SD = {s['overall_std_2dp']:.2f} degC; "
        f"|anom|>1.0 = {s['exceed_pct_1dp']:.1f}%",
        loc="left", fontsize=10,
    )

axes[-1].set_xlabel("Year")
axes[-1].set_xlim(target_time[0], target_time[-1])
fig2.suptitle(
    "Monthly Nino 3.4 SST anomalies (own 1950-2014 monthly climatology)",
    fontsize=12,
)
