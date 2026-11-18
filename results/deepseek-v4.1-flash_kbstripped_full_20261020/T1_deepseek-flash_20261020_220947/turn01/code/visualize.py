# Figure 1: seasonal cycle of Nino 3.4 amplitude (standard deviation by calendar month).
month_labels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
month_nums = np.arange(1, 13)

colors = {
    "ERSSTv5": "black",
    "CESM2": "tab:blue",
    "GFDL-ESM4": "tab:orange",
    "MPI-ESM1-2-LR": "tab:green",
}
linestyles = {
    "ERSSTv5": "-",
    "CESM2": "--",
    "GFDL-ESM4": "-.",
    "MPI-ESM1-2-LR": ":",
}

fig1, ax1 = plt.subplots(figsize=(8.5, 5.2), constrained_layout=True)

for src in sources:
    vals = monthly_std[src]
    pm = peak_month[src]
    pv = peak_value[src]
    lw = 3.0 if src == "ERSSTv5" else 1.9
    ax1.plot(
        month_nums,
        vals,
        color=colors[src],
        linestyle=linestyles[src],
        linewidth=lw,
        marker="o",
        markersize=5,
        label=f"{src} (peak {month_labels[pm - 1]}, {pv:.2f} degC)",
    )
    ax1.scatter(
        [pm],
        [pv],
        s=110 if src == "ERSSTv5" else 65,
        color=colors[src],
        edgecolor="black",
        zorder=5,
    )

ax1.set_xticks(month_nums)
ax1.set_xticklabels(month_labels)
ax1.set_xlabel("Calendar month")
ax1.set_ylabel("Nino 3.4 anomaly standard deviation (degC)")
ax1.set_title("Seasonal cycle of ENSO amplitude, 1950-2014 (own monthly climatology)")
ax1.grid(True, alpha=0.3)
ax1.legend(loc="best", fontsize=8, frameon=True)

# Figure 2: four stacked Nino 3.4 anomaly time series on a shared y-axis range.
all_vals = np.concatenate([index_series[src].values.ravel() for src in sources])
ymin = float(np.nanmin(all_vals))
ymax = float(np.nanmax(all_vals))
pad = 0.06 * max(ymax - ymin, 0.5)
shared_ylim = (min(ymin - pad, -1.1), max(ymax + pad, 1.1))

fig2, axes = plt.subplots(
    4,
    1,
    figsize=(11, 9),
    sharex=True,
    constrained_layout=True,
)

for ax, src in zip(axes, sources):
    da = index_series[src]
    ax.plot(
        da["time"].values,
        da.values,
        color=colors[src],
        linewidth=1.0,
    )
    ax.axhline(0.0, color="gray", linewidth=0.5)
    ax.axhline(1.0, color="red", linestyle="--", linewidth=0.8)
    ax.axhline(-1.0, color="red", linestyle="--", linewidth=0.8)
    ax.set_ylim(shared_ylim)
    ax.set_ylabel("Anomaly (degC)")
    ax.grid(True, alpha=0.25)
    ax.set_title(
        f"{src}: overall std={overall_std[src]:.2f} degC; "
        f"|anomaly| > 1.0 degC in {exceed_pct[src]:.1f}% of months",
        loc="left",
        fontsize=10,
    )

axes[-1].set_xlabel("Year")
fig2.suptitle(
    "Nino 3.4 SST anomaly, 1950-2014 (dashed red lines: +/-1.0 degC)",
    fontsize=12,
)
