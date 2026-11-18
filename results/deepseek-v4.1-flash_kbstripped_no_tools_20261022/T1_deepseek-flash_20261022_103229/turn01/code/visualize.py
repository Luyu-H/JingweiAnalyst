month_labels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

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
markers = {
    "ERSSTv5": "o",
    "CESM2": "s",
    "GFDL-ESM4": "^",
    "MPI-ESM1-2-LR": "D",
}

# Figure 1: seasonal cycle of Nino 3.4 anomaly amplitude
fig1, ax1 = plt.subplots(figsize=(10, 6), constrained_layout=True)
for src in sources:
    s = stats[src]
    ax1.plot(
        np.arange(1, 13),
        s["std_month"],
        color=colors[src],
        linestyle=linestyles[src],
        marker=markers[src],
        linewidth=3.0 if src == "ERSSTv5" else 1.8,
        markersize=5,
        label=src,
    )
    pk = s["peak_month"]
    pv = s["peak_value"]
    ax1.scatter(
        [pk],
        [pv],
        s=90,
        color=colors[src],
        marker=markers[src],
        edgecolor="white",
        linewidth=1.0,
        zorder=5,
    )
    ax1.annotate(
        month_labels[pk - 1],
        xy=(pk, pv),
        xytext=(5, 5),
        textcoords="offset points",
        fontsize=8,
        color=colors[src],
    )

ax1.set_xticks(np.arange(1, 13))
ax1.set_xticklabels(month_labels)
ax1.set_xlabel("Calendar month")
ax1.set_ylabel("Nino 3.4 anomaly standard deviation (degC)")
ax1.set_title(
    "Seasonal cycle of Nino 3.4 amplitude, 1950-2014 "
    "(markers mark each curve's peak month)"
)
ax1.grid(True, alpha=0.3)
ax1.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=True)

# Figure 2: monthly Nino 3.4 anomaly time series, four stacked panels
all_anom = np.concatenate([stats[src]["anomaly"] for src in sources])
ymin = float(np.nanmin(all_anom))
ymax = float(np.nanmax(all_anom))
pad = 0.08 * (ymax - ymin)
ylim = (ymin - pad, ymax + pad)

fig2, axes = plt.subplots(
    4,
    1,
    figsize=(12, 10),
    sharex=True,
    sharey=True,
    constrained_layout=True,
)

for ax, src in zip(axes, sources):
    s = stats[src]
    ax.plot(s["decimal_year"], s["anomaly"], color=colors[src], linewidth=1.0)
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.8)
    ax.axhline(-1.0, color="gray", linestyle="--", linewidth=0.8)
    ax.axhline(0.0, color="gray", linestyle=":", linewidth=0.6)
    ax.set_ylim(ylim)
    ax.set_ylabel("Anomaly (degC)")
    ax.set_title(
        f"{src}: std = {s['std_all']:.2f} degC; "
        f"|anomaly| > 1.0: {s['n_exceed']} months ({s['pct']:.1f}%)"
    )
    ax.grid(True, alpha=0.3)

axes[-1].set_xlabel("Year")
axes[-1].set_xlim(1950, 2015)
axes[-1].set_xticks(np.arange(1950, 2016, 10))
fig2.suptitle(
    "Nino 3.4 SST anomaly, 1950-2014; each series relative to its own "
    "1950-2014 monthly climatology",
    fontsize=12,
)
