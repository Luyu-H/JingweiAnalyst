import numpy as _np
import matplotlib as _mpl

# ---------------------------------------------------------------------------
# Figure 1 - global maps of operator MLD minus mlotst (display grid only)
# ---------------------------------------------------------------------------
pooled = []
for model in models:
    for mname, _ in month_defs:
        v = panels[(model, mname)]["diff_map"].values
        v = v[_np.isfinite(v)]
        if v.size:
            pooled.append(_np.abs(v))
pooled = _np.concatenate(pooled)
vmax = float(_np.nanpercentile(pooled, 98))

fig1, axes1 = plt.subplots(3, 2, figsize=(14, 11), constrained_layout=True)
cmap = plt.get_cmap("RdBu_r")
pm = None
for i, model in enumerate(models):
    for j, (mname, _) in enumerate(month_defs):
        ax = axes1[i, j]
        da = panels[(model, mname)]["diff_map"]
        lon = _np.asarray(da["lon"].values)
        lat = _np.asarray(da["lat"].values)
        pm = ax.pcolormesh(lon, lat, da.values, cmap=cmap, vmin=-vmax, vmax=vmax,
                           shading="auto", rasterized=True)
        ax.set_facecolor("0.85")
        ax.set_xlim(-180, 180); ax.set_ylim(-80, 90)
        d = results[(model, mname)]
        ax.set_title("%s  %s\nbias %+.1f m, RMSD %.1f m, r %.3f"
                     % (model, mname, d["bias"], d["rmse"], d["corr"]), fontsize=10)
        ax.set_xticks([-180, -120, -60, 0, 60, 120, 180])
        ax.set_yticks([-60, -30, 0, 30, 60])
        ax.tick_params(labelsize=8)
        if j == 0:
            ax.set_ylabel("latitude")
        if i == 2:
            ax.set_xlabel("longitude")

cb1 = fig1.colorbar(pm, ax=axes1, orientation="horizontal", fraction=0.035, pad=0.03,
                    extend="both")
cb1.set_label("operator MLD - mlotst (m)")
fig1.suptitle("Figure 1  Offline MLD minus published mlotst, March (left) and September (right)\n"
              "native-grid fields, remapped to a 1-deg grid for display only; shared scale clipped at the pooled 98th percentile (%.0f m)" % vmax,
              fontsize=12)

# ---------------------------------------------------------------------------
# Figure 2 - scatter of operator MLD vs mlotst on native common-mask cells
# ---------------------------------------------------------------------------
pos_all = []
for model in models:
    for mname, _ in month_defs:
        x = panels[(model, mname)]["x"]; y = panels[(model, mname)]["y"]
        sel = _np.isfinite(x) & _np.isfinite(y) & (x > 0) & (y > 0)
        if sel.any():
            pos_all.append(x[sel]); pos_all.append(y[sel])
pos_all = _np.concatenate(pos_all)
lo = float(_np.nanpercentile(pos_all, 0.5))
hi = float(_np.nanpercentile(pos_all, 99.5))
lo = max(lo, 0.5); hi = max(hi, lo * 10.0)

log_ticks = _np.array([1, 3, 10, 30, 100, 300, 1000, 3000, 10000])
log_ticks = log_ticks[(log_ticks >= lo) & (log_ticks <= hi)]

fig2, axes2 = plt.subplots(3, 2, figsize=(12, 12), constrained_layout=True)
for i, model in enumerate(models):
    for j, (mname, _) in enumerate(month_defs):
        ax = axes2[i, j]
        x = panels[(model, mname)]["x"]; y = panels[(model, mname)]["y"]
        sel = _np.isfinite(x) & _np.isfinite(y) & (x > 0) & (y > 0)
        ax.scatter(x[sel], y[sel], s=1.5, alpha=0.25, c="tab:blue", edgecolors="none",
                   rasterized=True)
        ax.plot([lo, hi], [lo, hi], "k--", lw=1.0, label="1:1")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_xticks(log_ticks); ax.set_yticks(log_ticks)
        ax.set_xticklabels([("%g" % t) for t in log_ticks], fontsize=8)
        ax.set_yticklabels([("%g" % t) for t in log_ticks], fontsize=8)
        ax.grid(True, which="both", lw=0.3, alpha=0.4)
        d = results[(model, mname)]
        ax.set_title("%s  %s" % (model, mname), fontsize=10)
        ax.text(0.03, 0.97,
                "bias %+.1f m\nRMSD %.1f m\nr %.3f\n>20 m: %.1f %%"
                % (d["bias"], d["rmse"], d["corr"], d["exceed_pct"]),
                transform=ax.transAxes, va="top", ha="left", fontsize=9,
                bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=0.9))
        if j == 0:
            ax.set_ylabel("mlotst (m)")
        if i == 2:
            ax.set_xlabel("operator MLD (m)")
        if i == 0 and j == 0:
            ax.legend(loc="lower right", fontsize=8)

fig2.suptitle("Figure 2  Offline MLD vs published mlotst on native common-mask cells\n"
              "identical logarithmic axes in all panels; rows CESM2, GFDL-ESM4, MPI-ESM1-2-LR",
              fontsize=12)
