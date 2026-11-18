# Figure 1: SSS bias maps
fig1, axes1 = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)
pooled_bias = np.concatenate([bias_fields[m].values.ravel() for m in MODELS])
vmax_bias = float(np.nanpercentile(np.abs(pooled_bias), 98))
if not np.isfinite(vmax_bias) or vmax_bias == 0:
    vmax_bias = 0.1

for ax, model in zip(axes1, MODELS):
    da = bias_fields[model]
    im = ax.pcolormesh(
        woa_lon, woa_lat, da.values,
        cmap="RdBu_r", vmin=-vmax_bias, vmax=vmax_bias, shading="auto",
    )
    ax.set_title(model, color=COLORS[model], fontweight="bold")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

cbar1 = fig1.colorbar(
    im, ax=axes1, orientation="horizontal", fraction=0.05, pad=0.06, extend="both"
)
cbar1.set_label("SSS bias (psu)")
fig1.suptitle(
    "1971-2000 SSS bias (model - WOA23); shared scale clipped at pooled 98th percentile"
)

# Figure 2: 1950-2014 annual-mean SSS trend maps with p<0.05 stippling
fig2, axes2 = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)
pooled_trend = np.concatenate([trend_fields[m].values.ravel() for m in MODELS])
vmax_trend = float(np.nanpercentile(np.abs(pooled_trend), 98))
if not np.isfinite(vmax_trend) or vmax_trend == 0:
    vmax_trend = 0.01

for ax, model in zip(axes2, MODELS):
    da = trend_fields[model]
    im2 = ax.pcolormesh(
        woa_lon, woa_lat, da.values,
        cmap="RdBu_r", vmin=-vmax_trend, vmax=vmax_trend, shading="auto",
    )
    sig = trend_sig_fields[model].values.astype(float)
    cs = ax.contourf(
        woa_lon, woa_lat, sig,
        levels=[0.5, 1.5], colors="none", hatches=["...."], add_colorbar=False,
    )
    for coll in cs.collections:
        coll.set_edgecolor("k")
        coll.set_linewidth(0.0)
    ax.set_title(model, color=COLORS[model], fontweight="bold")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

cbar2 = fig2.colorbar(
    im2, ax=axes2, orientation="horizontal", fraction=0.05, pad=0.06, extend="both"
)
cbar2.set_label("SSS trend (psu/decade)")
fig2.suptitle(
    "1950-2014 annual-mean SSS trend; stippling marks p<0.05 (ESS); shared scale clipped at 98th percentile"
)

# Figure 3: global mean SSS time series with separate segment fits
fig3, ax3 = plt.subplots(figsize=(9, 5.2), constrained_layout=True)

for model in MODELS:
    c = COLORS[model]
    sh = series_hist[model]
    ss = series_ssp[model]

    xh = np.asarray(sh["year"].values, dtype=float)
    yh = np.asarray(sh.values, dtype=float)
    xs = np.asarray(ss["year"].values, dtype=float)
    ys = np.asarray(ss.values, dtype=float)

    ax3.plot(xh, yh, color=c, linewidth=1.8, label=f"{model} historical")
    ax3.plot(xs, ys, color=c, linewidth=1.8, label=f"{model} ssp245")

    fh = fit_params_hist[model]
    fs = fit_params_ssp[model]
    ax3.plot(
        xh, fh["intercept"] + fh["slope"] * xh,
        color=c, linestyle="--", linewidth=1.0,
    )
    ax3.plot(
        xs, fs["intercept"] + fs["slope"] * xs,
        color=c, linestyle="--", linewidth=1.0,
    )

# Leave 2015-2070 empty: the two segments are plotted separately with no connecting line.
ax3.set_xlim(1945, 2105)
ax3.set_xlabel("Year")
ax3.set_ylabel("Global mean SSS (psu)")
ax3.set_title("Global mean 0-10 m SSS; dashed lines are linear fits with ESS")
ax3.legend(loc="best", fontsize=8, ncol=2)
ax3.grid(alpha=0.25)
