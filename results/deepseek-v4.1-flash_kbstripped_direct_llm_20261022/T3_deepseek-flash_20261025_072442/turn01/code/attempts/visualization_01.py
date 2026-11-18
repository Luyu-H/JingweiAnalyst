import numpy as np
import matplotlib.pyplot as plt

X, Y = np.meshgrid(woa_lon, woa_lat)

bias_fields = {m: np.where(common, tos_on_woa[m] - woa_arr, np.nan) for m in MODELS}
woa_field   = np.where(common, woa_arr, np.nan)

allb = np.concatenate([bias_fields[m][np.isfinite(bias_fields[m])].ravel() for m in MODELS])
vmax = np.nanpercentile(np.abs(allb), 99)
vmax = np.ceil(vmax * 10) / 10.0

# ---------------- Figure 1 ----------------
fig, axes = plt.subplots(2, 2, figsize=(15, 9), constrained_layout=True)

for ax, m in zip([axes[0, 0], axes[0, 1], axes[1, 0]], MODELS):
    pm = ax.pcolormesh(X, Y, bias_fields[m], cmap="RdBu_r",
                       vmin=-vmax, vmax=vmax, shading="auto")
    cb = fig.colorbar(pm, ax=ax, orientation="horizontal", pad=0.06, fraction=0.05)
    cb.set_label("SST bias (model $-$ WOA23) [degrees Celsius]")
    ax.set_title(f"{m}: SST bias vs WOA23 (1971-2000)")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

ax = axes[1, 1]
pm = ax.pcolormesh(X, Y, woa_field, cmap="viridis", shading="auto")
cb = fig.colorbar(pm, ax=ax, orientation="horizontal", pad=0.06, fraction=0.05)
cb.set_label("WOA23 SST [degrees Celsius]")
ax.set_title("WOA23 0 m annual mean SST, 1971-2000 (reference)")
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
ax.set_xlim(-180, 180)
ax.set_ylim(-90, 90)

fig.suptitle("CMIP6 historical SST climatology (1971-2000) against WOA23", fontsize=14)

# ---------------- Figure 2 ----------------
labels = MODELS + ["WOA23"]
vals   = [native_global_mean[m] for m in MODELS] + [woa_global_mean]

fig2, ax2 = plt.subplots(figsize=(8, 5.5), constrained_layout=True)
colors = ["#4C72B0", "#55A868", "#C44E52", "#8172B2"]
bars = ax2.bar(labels, vals, color=colors, edgecolor="black", linewidth=0.8)
for b, v in zip(bars, vals):
    ax2.annotate(f"{v:.2f}", (b.get_x() + b.get_width() / 2, v),
                 xytext=(0, 3), textcoords="offset points",
                 ha="center", va="bottom", fontsize=11)
ax2.set_ylabel("Global mean SST [degrees Celsius]")
ax2.set_title("Global mean SST, 1971-2000 climatology\n(each over its own native ocean grid)")
ax2.set_ylim(0, max(vals) * 1.18)
ax2.grid(axis="y", alpha=0.3)
