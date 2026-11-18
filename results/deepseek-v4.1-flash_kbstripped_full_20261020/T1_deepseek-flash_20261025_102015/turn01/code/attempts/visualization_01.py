# Figure 1: seasonal cycle of Nino 3.4 amplitude
colors = {
    "ERSSTv5": "black",
    "CESM2": "tab:red",
    "GFDL-ESM4": "tab:green",
    "MPI-ESM1-2-LR": "tab:blue",
}

fig1, ax1 = plt.subplots(figsize=(8, 5), constrained_layout=True)
for name in ["ERSSTv5"] + models:
    vals = monthly_std[name]
    ax1.plot(
        np.arange(1, 13),
        vals,
        marker="o",
        lw=2.2 if name == "ERSSTv5" else 1.4,
        color=colors[name],
        label=name,
    )
    pm = peak_month[name]
    ax1.scatter(
        pm,
        vals[pm - 1],
        s=90,
        facecolors="none",
        edgecolors=colors[name],
        linewidth=2,
        zorder=5,
    )
    ax1.annotate(
        month_labels[pm - 1],
        (pm, vals[pm - 1]),
        xytext=(0, 8),
        textcoords="offset points",
        ha="center",
        fontsize=9,
        color=colors[name],
    )

ax1.set_xticks(np.arange(1, 13))
ax1.set_xticklabels(month_labels)
ax1.set_xlabel("Calendar month")
ax1.set_ylabel("Standard deviation of Niño 3.4 anomaly (°C)")
ax1.set_title("Seasonal phase locking of ENSO amplitude, 1950–2014")
ax1.legend(loc="best")
ax1.grid(alpha=0.3)

# Figure 2: stacked monthly Nino 3.4 anomaly series
order = ["ERSSTv5"] + models
all_vals = np.concatenate([indices[name].values for name in order])
ymin = float(np.nanmin(all_vals))
ymax = float(np.nanmax(all_vals))
pad = 0.08 * (ymax - ymin) if ymax > ymin else 0.5
ylim = (ymin - pad, ymax + pad)

fig2, axes = plt.subplots(
    4,
    1,
    figsize=(12, 10),
    sharex=True,
    constrained_layout=True,
)

for ax, name in zip(axes, order):
    idx = indices[name]
    ax.plot(
        idx["time"].values,
        idx.values,
        lw=0.8,
        color=colors[name],
    )
    ax.axhline(1.0, color="0.4", ls="--", lw=0.9)
    ax.axhline(-1.0, color="0.4", ls="--", lw=0.9)
    ax.axhline(0.0, color="0.6", ls=":", lw=0.7)
    ax.set_ylim(ylim)
    ax.set_ylabel("Niño 3.4 anomaly (°C)")
    ax.text(
        0.01,
        0.96,
        f"{name}\n"
        f"SD = {std_all[name]:.2f} °C\n"
        f"|anom| > 1°C = {exceed_pct[name]:.1f}%",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox=dict(facecolor="white", alpha=0.85, edgecolor="none"),
    )
    ax.grid(alpha=0.25)

axes[0].text(
    0.99,
    0.95,
    "dashed lines: +1°C and −1°C",
    transform=axes[0].transAxes,
    ha="right",
    va="top",
    fontsize=8,
    color="0.3",
)
axes[-1].set_xlabel("Year")
axes[-1].set_xlim(pd.Timestamp("1950-01-01"), pd.Timestamp("2015-01-01"))
fig2.suptitle("Monthly Niño 3.4 anomaly, 1950–2014", fontsize=14)
