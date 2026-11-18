from matplotlib.ticker import FuncFormatter

# ---------------------------------------------------------------------------
# Figure 1: maps of operator MLD minus mlotst
# ---------------------------------------------------------------------------
fig1, axes1 = plt.subplots(
    3, 2, figsize=(13, 10), constrained_layout=True
)

all_diff_vals = []
for model in MODELS:
    for mname, _ in MONTHS:
        d = results[model]["months"][mname]["diff_reg"]
        all_diff_vals.append(d[np.isfinite(d)])
all_diff_vals = (
    np.concatenate(all_diff_vals) if all_diff_vals else np.array([0.0])
)
vmax = np.nanpercentile(np.abs(all_diff_vals), 98)
if not np.isfinite(vmax) or vmax <= 0:
    vmax = 1.0

cmap = plt.get_cmap("RdBu_r").copy()
cmap.set_bad("lightgray")

im = None
for i, model in enumerate(MODELS):
    for j, (mname, _) in enumerate(MONTHS):
        ax = axes1[i, j]
        d = results[model]["months"][mname]["diff_reg"]
        im = ax.imshow(
            d,
            extent=[-180, 180, -90, 90],
            origin="lower",
            cmap=cmap,
            vmin=-vmax,
            vmax=vmax,
            aspect="auto",
            interpolation="nearest",
        )
        ax.set_facecolor("lightgray")
        ax.set_title(f"{model} — {mname}")
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        ax.set_xticks([-180, -120, -60, 0, 60, 120, 180])
        ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

cbar = fig1.colorbar(
    im, ax=axes1.ravel().tolist(), shrink=0.7, pad=0.02
)
cbar.set_label("Operator MLD − mlotst (m)")
fig1.suptitle(
    "Operator MLD minus published mlotst, 1971–2000 climatology",
    y=1.02,
)

# ---------------------------------------------------------------------------
# Figure 2: operator MLD vs mlotst scatter, common log axes
# ---------------------------------------------------------------------------
all_pos = []
for model in MODELS:
    for mname, _ in MONTHS:
        r = results[model]["months"][mname]
        v = r["valid"]
        x = r["op_clim"][v]
        y = r["ml_clim"][v]
        all_pos.append(x[x > 0])
        all_pos.append(y[y > 0])
all_pos = (
    np.concatenate([a for a in all_pos if a.size])
    if all_pos else np.array([1.0])
)

vmin = np.nanpercentile(all_pos, 0.5)
vmax = np.nanpercentile(all_pos, 99.5)
vmin = max(float(vmin), 0.05)
vmax = max(float(vmax), vmin * 10.0)
vmin = 10.0 ** np.floor(np.log10(vmin))
vmax = 10.0 ** np.ceil(np.log10(vmax))

tick_candidates = np.array(
    [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000]
)
ticks = tick_candidates[
    (tick_candidates >= vmin) & (tick_candidates <= vmax)
]

fig2, axes2 = plt.subplots(
    3, 2, figsize=(12, 11), constrained_layout=True
)

for i, model in enumerate(MODELS):
    for j, (mname, _) in enumerate(MONTHS):
        ax = axes2[i, j]
        r = results[model]["months"][mname]
        v = r["valid"]
        x = r["op_clim"][v]
        y = r["ml_clim"][v]
        pos = (x > 0) & (y > 0)

        ax.scatter(
            x[pos], y[pos],
            s=2, alpha=0.25, c="k", edgecolors="none", rasterized=True,
        )
        ax.plot([vmin, vmax], [vmin, vmax], "k--", lw=1.0)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(vmin, vmax)
        ax.set_ylim(vmin, vmax)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.minorticks_off()
        ax.xaxis.set_major_formatter(
            FuncFormatter(lambda val, pos: f"{val:g}")
        )
        ax.yaxis.set_major_formatter(
            FuncFormatter(lambda val, pos: f"{val:g}")
        )
        ax.set_xlabel("Operator MLD (m)")
        ax.set_ylabel("mlotst (m)")
        ax.set_title(f"{model} — {mname}")
        ax.grid(True, which="major", alpha=0.3)

        txt = (
            f"bias = {r['bias_w']:.1f} m\n"
            f"RMSE = {r['rmse_w']:.1f} m\n"
            f"r = {r['corr_w']:.3f}"
        )
        ax.text(
            0.04, 0.96, txt,
            transform=ax.transAxes,
            va="top", ha="left",
            fontsize=8,
            bbox=dict(boxstyle="round", fc="white", ec="0.5", alpha=0.8),
        )

fig2.suptitle(
    "Operator MLD vs published mlotst, 1971–2000 monthly climatologies",
    y=1.02,
)
