# Visualization code only

plot_order = ["ERSSTv5", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
plot_colors = {
    "ERSSTv5": "black",
    "CESM2": "tab:blue",
    "GFDL-ESM4": "tab:orange",
    "MPI-ESM1-2-LR": "tab:green",
}
plot_linestyles = {
    "ERSSTv5": "-",
    "CESM2": "-",
    "GFDL-ESM4": "--",
    "MPI-ESM1-2-LR": "-.",
}

# Figure 1: seasonal cycle of Niño 3.4 anomaly standard deviation
fig1, ax1 = plt.subplots(figsize=(10, 5), constrained_layout=True)

x = np.arange(1, 13)
xlabels = [month_names[m] for m in x]

for src in plot_order:
    y = month_std_values[src]
    lw = 2.8 if src == "ERSSTv5" else 2.0
    ax1.plot(
        x, y,
        color=plot_colors[src],
        linestyle=plot_linestyles[src],
        linewidth=lw,
        label=src,
    )
    ax1.scatter(
        [peak_month[src]],
        [peak_value[src]],
        color=plot_colors[src],
        s=50 if src == "ERSSTv5" else 35,
        zorder=4,
    )

ax1.set_xticks(x)
ax1.set_xticklabels(xlabels)
ax1.set_xlabel("Calendar month")
ax1.set_ylabel("Std. dev. of Niño 3.4 anomaly (°C)")
ax1.set_title("Figure 1. Seasonal cycle of Niño 3.4 amplitude, 1950–2014")
ax1.grid(True, alpha=0.3)
ax1.legend(loc="best")

# Figure 2: stacked Niño 3.4 anomaly time series
all_vals = np.concatenate([np.asarray(nino34_anom[src].values).ravel() for src in plot_order])
ymax = float(np.nanmax(np.abs(all_vals)))
ymax = max(1.2, np.ceil(ymax * 10) / 10)

fig2, axes = plt.subplots(4, 1, figsize=(12, 9), sharex=True, sharey=True, constrained_layout=True)

for ax, src in zip(axes, plot_order):
    da = nino34_anom[src]
    time_vals = pd.to_datetime(da["time"].values)
    y_vals = np.asarray(da.values)

    ax.plot(time_vals, y_vals, color=plot_colors[src], linewidth=1.0)
    ax.axhline(0.0, color="0.4", linewidth=0.8)
    ax.axhline(1.0, color="red", linestyle="--", linewidth=0.9)
    ax.axhline(-1.0, color="red", linestyle="--", linewidth=0.9)
    ax.set_ylim(-ymax, ymax)
    ax.set_ylabel("°C")
    ax.grid(True, alpha=0.25)
    ax.set_title(src, loc="left")
    ax.text(
        0.995, 0.92,
        f"σ = {std_full[src]:.2f} °C | |anom| > 1°C: {exceed_percent[src]:.1f}%",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=2.5),
    )

axes[-1].set_xlabel("Year")
fig2.suptitle("Figure 2. Niño 3.4 SST anomaly time series, 1950–2014", y=1.02)
