# Visualization code only
import matplotlib.ticker as mticker

# ---------------------------------------------------------------------------
# Figure 1: maps of operator MLD minus mlotst on a common lon-lat frame.
# ---------------------------------------------------------------------------
reg_fields = {}
pooled_abs = []

for model in MODELS:
    for month in MONTHS:
        if (model, month) not in plot_data:
            continue

        d = plot_data[(model, month)]
        LON, LAT, Z, M = interp_native_to_regular(
            d["lon2d"],
            d["lat2d"],
            d["diff"],
            d["common"],
            d["ocean_mask"],
            reg_lon,
            reg_lat,
        )

        reg_fields[(model, month)] = (LON, LAT, Z, M)

        finite = np.isfinite(Z)
        if np.any(finite):
            pooled_abs.append(np.abs(Z[finite]))

if pooled_abs:
    pooled_abs = np.concatenate(pooled_abs)
    vmax = float(np.nanpercentile(pooled_abs, 98.0))
else:
    vmax = 1.0

if not np.isfinite(vmax) or vmax <= 0.0:
    vmax = 1.0

fig, axes = plt.subplots(
    len(MODELS),
    len(MONTHS),
    figsize=(13.5, 10.5),
    constrained_layout=True,
)

for i, model in enumerate(MODELS):
    for j, month in enumerate(MONTHS):
        ax = axes[i, j]

        if (model, month) not in reg_fields:
            ax.axis("off")
            continue

        LON, LAT, Z, M = reg_fields[(model, month)]

        ax.pcolormesh(
            LON,
            LAT,
            Z,
            cmap="RdBu_r",
            vmin=-vmax,
            vmax=vmax,
            shading="auto",
        )

        if np.any(np.isfinite(M)):
            ax.contour(
                LON,
                LAT,
                M,
                levels=[0.5],
                colors="k",
                linewidths=0.35,
            )

        ax.set_xlim(0.0, 360.0)
        ax.set_ylim(-90.0, 90.0)
        ax.set_xticks(np.arange(0.0, 361.0, 60.0))
        ax.set_yticks(np.arange(-90.0, 91.0, 30.0))
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_facecolor("white")
        ax.set_title(f"{model} {MONTH_NAMES[month]}", fontsize=9)

mappable = plt.cm.ScalarMappable(
    norm=plt.Normalize(vmin=-vmax, vmax=vmax),
    cmap="RdBu_r",
)
mappable.set_array([])

cbar = fig.colorbar(
    mappable,
    ax=axes.ravel().tolist(),
    orientation="horizontal",
    fraction=0.04,
    pad=0.04,
    extend="both",
)
cbar.set_label("MLD_op - mlotst (m)")

fig.suptitle(
    "Operator MLD minus published mlotst (m)\n"
    f"Shared diverging scale, clipped at pooled 98th percentile of |difference| = {vmax:.1f} m",
    fontsize=11,
)

# ---------------------------------------------------------------------------
# Figure 2: scatter plots, operator MLD vs mlotst, log axes.
# ---------------------------------------------------------------------------
positive_values = []

for model in MODELS:
    for month in MONTHS:
        if (model, month) not in scatter_data:
            continue

        d = scatter_data[(model, month)]
        m = d["common"] & (d["op"] > 0.0) & (d["ml"] > 0.0)

        if np.any(m):
            positive_values.append(d["op"][m])
            positive_values.append(d["ml"][m])

if positive_values:
    positive_values = np.concatenate(positive_values)
    lo = float(np.nanmin(positive_values))
    hi = float(np.nanmax(positive_values))
    vmin = 10.0 ** np.floor(np.log10(lo))
    vmax_sc = 10.0 ** np.ceil(np.log10(hi))
else:
    vmin = 0.1
    vmax_sc = 1000.0

if vmax_sc <= vmin:
    vmax_sc = vmin * 10.0

ticks = []
for decade in range(
    int(np.floor(np.log10(vmin))), int(np.ceil(np.log10(vmax_sc))) + 1
):
    for mult in [1.0, 2.0, 5.0]:
        t = mult * (10.0 ** decade)
        if vmin <= t <= vmax_sc:
            ticks.append(t)
ticks = np.asarray(ticks)

fig2, axes2 = plt.subplots(
    len(MODELS),
    len(MONTHS),
    figsize=(11.5, 10.5),
    constrained_layout=True,
)

for i, model in enumerate(MODELS):
    for j, month in enumerate(MONTHS):
        ax = axes2[i, j]

        if (model, month) not in scatter_data:
            ax.axis("off")
            continue

        d = scatter_data[(model, month)]
        common = d["common"]
        op = d["op"]
        ml = d["ml"]

        valid = common & (op > 0.0) & (ml > 0.0)
        omitted = int(np.sum(common & ((op <= 0.0) | (ml <= 0.0))))

        ax.scatter(
            op[valid],
            ml[valid],
            s=1.0,
            alpha=0.25,
            color="tab:blue",
            edgecolors="none",
        )

        ax.plot([vmin, vmax_sc], [vmin, vmax_sc], "k--", lw=0.8)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(vmin, vmax_sc)
        ax.set_ylim(vmin, vmax_sc)

        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda v, pos: f"{v:g}")
        )
        ax.yaxis.set_major_formatter(
            mticker.FuncFormatter(lambda v, pos: f"{v:g}")
        )

        ax.grid(True, which="major", alpha=0.2)

        row = results_df[
            (results_df["model"] == model) & (results_df["month"] == month)
        ].iloc[0]

        ax.set_title(
            f"{model} {MONTH_NAMES[month]}\n"
            f"bias={row['bias']:.1f} m, RMSD={row['rmsd']:.1f} m, r={row['corr']:.3f}\n"
            f"omitted non-positive: {omitted}",
            fontsize=8,
        )

        ax.set_xlabel("Operator MLD (m)")
        ax.set_ylabel("mlotst (m)")

fig2.suptitle(
    "Operator MLD vs published mlotst: one point per common-mask ocean cell\n"
    "Identical logarithmic axes and ticks in all panels; dashed line is y = x",
    fontsize=11,
)
