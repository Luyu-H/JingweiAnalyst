import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

model_colors = {
    "CESM2": "tab:blue",
    "GFDL-ESM4": "tab:orange",
    "MPI-ESM1-2-LR": "tab:green",
}

# ---------------------------------------------------------------------------
# Figure 1: SSS bias maps (model 1971-2000 climatology minus WOA23)
# ---------------------------------------------------------------------------
fig1, axes1 = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)

absdiff = np.concatenate([np.abs(bias_map[m].values).ravel() for m in models])
absdiff = absdiff[np.isfinite(absdiff)]
vmax_bias = float(np.nanpercentile(absdiff, 98)) if absdiff.size else 0.1
if vmax_bias <= 0:
    vmax_bias = 0.1

for ax, model in zip(axes1, models):
    da = bias_map[model]
    im1 = ax.pcolormesh(
        da.lon.values,
        da.lat.values,
        da.values,
        shading="auto",
        cmap="RdBu_r",
        vmin=-vmax_bias,
        vmax=vmax_bias,
    )
    ax.set_title(model)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_aspect("auto")

cbar1 = fig1.colorbar(
    im1,
    ax=axes1,
    orientation="horizontal",
    fraction=0.05,
    pad=0.06,
    extend="both",
)
cbar1.set_label("SSS bias (psu)")
fig1.suptitle(
    "SSS bias vs WOA23 1971-2000 (psu); shared scale clipped at 98th pct "
    f"|bias|={vmax_bias:.2f} psu"
)

# ---------------------------------------------------------------------------
# Figure 2: historical 1950-2014 annual-mean SSS trend maps with significance
# ---------------------------------------------------------------------------
fig2, axes2 = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)

abs_slope = np.concatenate([np.abs(slope_regrid_hist[m].values).ravel() for m in models])
abs_slope = abs_slope[np.isfinite(abs_slope)]
vmax_slope = float(np.nanpercentile(abs_slope, 98)) if abs_slope.size else 0.1
if vmax_slope <= 0:
    vmax_slope = 0.1

for ax, model in zip(axes2, models):
    da = slope_regrid_hist[model]
    im2 = ax.pcolormesh(
        da.lon.values,
        da.lat.values,
        da.values,
        shading="auto",
        cmap="RdBu_r",
        vmin=-vmax_slope,
        vmax=vmax_slope,
    )

    sig = sig_regrid_hist[model]
    lon2d, lat2d = np.meshgrid(sig.lon.values, sig.lat.values)
    sv = sig.values.astype(bool)
    ax.scatter(
        lon2d[sv],
        lat2d[sv],
        s=0.05,
        c="k",
        marker=".",
        linewidths=0,
        alpha=0.5,
    )

    ax.set_title(model)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_aspect("auto")

cbar2 = fig2.colorbar(
    im2,
    ax=axes2,
    orientation="horizontal",
    fraction=0.05,
    pad=0.06,
    extend="both",
)
cbar2.set_label("SSS trend (psu per decade)")
fig2.suptitle(
    "Annual-mean SSS trend 1950-2014; stippling p<0.05; "
    f"shared 98th pct |slope|={vmax_slope:.3g} psu/dec"
)

# ---------------------------------------------------------------------------
# Figure 3: global mean SSS time series, historical and ssp245
# ---------------------------------------------------------------------------
fig3, ax3 = plt.subplots(figsize=(10, 5), constrained_layout=True)

for model in models:
    color = model_colors[model]

    sh = hist_global_annual[model]
    ss = ssp_global_annual[model]

    ax3.plot(
        sh.year.values,
        sh.values,
        color=color,
        lw=1.2,
        marker=".",
        ms=2,
    )
    ax3.plot(
        ss.year.values,
        ss.values,
        color=color,
        lw=1.2,
        marker=".",
        ms=2,
    )

    fh = trend_hist[model]["fitted_data"]
    fs = trend_ssp[model]["fitted_data"]

    ax3.plot(fh.year.values, fh.values, color=color, lw=2.0, ls="--")
    ax3.plot(fs.year.values, fs.values, color=color, lw=2.0, ls="--")

ax3.set_xlabel("Year")
ax3.set_ylabel("Global mean SSS (psu)")
ax3.set_xlim(1945, 2105)
ax3.grid(True, alpha=0.3)

handles = [
    Line2D([0], [0], color=model_colors[m], lw=2.0)
    for m in models
]
ax3.legend(
    handles,
    models,
    loc="upper left",
    bbox_to_anchor=(1.02, 1.0),
    frameon=False,
)
fig3.suptitle(
    "Global mean SSS: historical 1950-2014 and ssp245 2071-2100; "
    "dashed lines are linear fits"
)
