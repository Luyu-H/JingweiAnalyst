# ---------------------------------------------------------------------------
# Figure 1 - seasonal cycle of Nino 3.4 amplitude
# ---------------------------------------------------------------------------
months = np.arange(1, 13)
month_labels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
model_colors = {"CESM2": "tab:red",
                "GFDL-ESM4": "tab:green",
                "MPI-ESM1-2-LR": "tab:orange"}

fig, ax = plt.subplots(figsize=(9, 5.5))

ax.plot(months, calendar_std["ERSSTv5"], color="black", lw=2.8,
        marker="o", ms=7, zorder=5, label="ERSSTv5 (observations)")

for m in model_names:
    ax.plot(months, calendar_std[m], color=model_colors[m], lw=1.7,
            ls="--", marker="s", ms=4, alpha=0.9, label=m)

# mark each curve's peak month
ax.scatter([peak_month["ERSSTv5"]], [peak_value["ERSSTv5"]],
           marker="*", s=260, color="black", edgecolor="gold",
           linewidth=1.0, zorder=6)
for m in model_names:
    ax.scatter([peak_month[m]], [peak_value[m]],
               marker="*", s=200, color=model_colors[m],
               edgecolor="black", linewidth=0.6, zorder=6)

ax.set_xticks(months)
ax.set_xticklabels(month_labels)
ax.set_xlabel("Calendar month")
ax.set_ylabel("Nino 3.4 SST anomaly std (degC)")
ax.set_title("Figure 1 - Seasonal cycle of Nino 3.4 amplitude\n"
             "1950-2014, self monthly climatology; stars mark each curve's peak month")
ax.grid(alpha=0.3)
ax.legend(loc="upper left", framealpha=0.9)
fig.tight_layout()

# ---------------------------------------------------------------------------
# Figure 2 - the four monthly anomaly series, stacked, shared axes
# ---------------------------------------------------------------------------
all_vals = np.concatenate([np.asarray(series[s].values, dtype=float) for s in sources])
ylim = float(np.nanmax(np.abs(all_vals))) * 1.10

fig2, axes = plt.subplots(4, 1, figsize=(13, 10), sharex=True, sharey=True,
                          constrained_layout=True)

for ax, src in zip(axes, sources):
    s = series[src]
    color = "black" if src == "ERSSTv5" else model_colors[src]
    ax.axhline(0.0, color="0.5", lw=0.8)
    ax.axhline(1.0, color="tab:blue", lw=1.0, ls="--")
    ax.axhline(-1.0, color="tab:blue", lw=1.0, ls="--")
    ax.plot(s["time"].values, np.asarray(s.values, dtype=float),
            color=color, lw=0.9)
    ax.fill_between(s["time"].values, 1.0,
                    np.asarray(s.values, dtype=float),
                    where=np.asarray(s.values, dtype=float) > 1.0,
                    color="tab:red", alpha=0.35, interpolate=True)
    ax.fill_between(s["time"].values, -1.0,
                    np.asarray(s.values, dtype=float),
                    where=np.asarray(s.values, dtype=float) < -1.0,
                    color="tab:blue", alpha=0.35, interpolate=True)
    ax.set_ylim(-ylim, ylim)
    ax.set_ylabel("Anomaly (degC)")
    ax.grid(alpha=0.25)
    ax.set_title(
        f"{src}   |   std = {overall_std[src]:.2f} degC   |   "
        f"|anomaly| > 1.0 degC: {exceed_count[src]} months "
        f"({exceed_pct[src]:.1f}%)",
        loc="left", fontsize=10)

axes[0].text(0.995, 0.05, "+1.0 / -1.0 degC thresholds dashed",
             transform=axes[0].transAxes, ha="right", va="bottom", fontsize=8,
             bbox=dict(boxstyle="round", fc="white", alpha=0.8))
axes[-1].set_xlabel("Time")
fig2.suptitle("Figure 2 - Monthly Nino 3.4 SST anomaly, 1950-2014 "
              "(each source vs its own 1950-2014 monthly climatology)")
