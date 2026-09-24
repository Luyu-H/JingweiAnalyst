"""
Reference answer for case T1 - Multi-model bias of the global SST climatology.

Run from the repository root with the project's conda environment:

    source /home/lorraine/miniconda3/bin/activate base && conda activate jingwei_QA
    PYTHONPATH=. python cases/reference_ans/T1/reference_T1.py

Writes answer_T1.json, answer_T1.md, fig1_bias_maps.png and fig2_global_mean_bars.png
next to this file.

Every physical quantity is produced by a tool from src/tools; nothing here
re-implements an average, a regridding or an error metric by hand.
"""

import json
from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from src.tools.registry import auto_import_tools
auto_import_tools()
from src.tools.data_loader import load_regular_data, load_climatology_data, load_landsea_mask
from src.tools.preprocessing import regrid_horizontal, build_comparison_mask
from src.tools.basic_calculator import (compute_basic_stat, temporal_aggregate,
                                        compute_difference_metrics, analyze_correlation)

OUT = Path(__file__).resolve().parent
MODELS = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
TIME_RANGE = ["1971-01", "2000-12"]

# The loaders attach provenance attrs that tell the weighting tools which grid a
# field is on; without keep_attrs a plain xarray reduction would drop them.
xr.set_options(keep_attrs=True)


# ---------------------------------------------------------------------------
# Step 1-2. Model climatology on the native grid, and its native global mean
# ---------------------------------------------------------------------------
model_clim, native_mean = {}, {}
for model in MODELS:
    tos = load_regular_data("cmip6_historical", "tos", model=model, time_range=TIME_RANGE)
    assert tos.sizes["time"] == 360, f"{model}: {tos.sizes['time']} months, expected 360"

    # One field for the whole window, each month weighted by its length in the
    # model's own calendar, so a 31-day and a 28-day month do not count alike.
    clim = temporal_aggregate(tos, freq="mean", op="mean")
    incomplete = clim.attrs.get("incomplete_periods")
    assert not incomplete, f"{model}: incomplete years in the window: {incomplete}"
    model_clim[model] = clim

    # areacello weighting is picked up from the provenance attrs, which is the
    # only correct area weighting on CESM2's and MPI's curvilinear grids.
    native_mean[model] = float(
        compute_basic_stat(clim, op="area_weighted_mean", dim=["lat", "lon"]).compute()
    )
    print(f"[1-2] {model:14s} native global mean SST = {native_mean[model]:.4f} degC")


# ---------------------------------------------------------------------------
# Step 3. WOA23 reference field and its global mean
# ---------------------------------------------------------------------------
woa = load_climatology_data("woa23", "temperature", time_slice="annual",
                            period="1971-2000", depth_range=[0, 0]).squeeze("depth", drop=True)
woa_mean = float(compute_basic_stat(woa, op="area_weighted_mean", dim=["lat", "lon"]).compute())
print(f"[3]   WOA23 1971-2000 0 m global mean SST = {woa_mean:.4f} degC")

landsea = load_landsea_mask()          # 1 = land, ocean is (value > 1); used to shade land in fig. 1


# ---------------------------------------------------------------------------
# Step 4. Conservative regridding onto the WOA23 1-degree grid, common mask
# ---------------------------------------------------------------------------
regridded = {}
for model in MODELS:
    regridded[model] = regrid_horizontal(model_clim[model], target="woa23", method="conservative")
    print(f"[4]   {model:14s} regridded, {int(regridded[model].notnull().sum())} valid cells")

# One mask for all three models, so the three scores rest on an identical sample.
# build_comparison_mask also applies the WOA23 land-sea mask with its own
# semantics (1 = land, ocean is value > 1).
common = build_comparison_mask([regridded[m] for m in MODELS] + [woa], ocean_only=True)
n_common = int(common.sum())
print(f"[4]   common comparison mask: {n_common} cells ({common.attrs['ocean_mask']})")

woa_common = woa.where(common)
model_common = {m: regridded[m].where(common) for m in MODELS}


# ---------------------------------------------------------------------------
# Step 5. Bias, RMSE and pattern correlation on the common mask
# ---------------------------------------------------------------------------
metrics, bias_field = {}, {}
woa_common_mean = None
for model in MODELS:
    diff = compute_difference_metrics(model_common[model], woa_common,
                                      dim=["lat", "lon"], weighting="area")
    corr = analyze_correlation(model_common[model], woa_common, dim=["lat", "lon"],
                               method="pearson", weighting="area")
    metrics[model] = {
        "delta_sst_degC": float(diff["bias"]),
        "rmse_degC": float(diff["rmse"]),
        "centered_rmse_degC": float(diff["centered_rmse"]),
        "pattern_correlation": float(corr["correlation"]),
        "common_mask_mean_model_degC": float(diff["mean_data"]),
        "common_mask_mean_woa23_degC": float(diff["mean_reference"]),
        "n_valid_cells": int(diff["n_valid"]),
        "weighting": str(diff["weighting"]),
    }
    bias_field[model] = diff["difference"]
    woa_common_mean = float(diff["mean_reference"])
    print(f"[5]   {model:14s} dSST={metrics[model]['delta_sst_degC']:+.4f}  "
          f"RMSE={metrics[model]['rmse_degC']:.4f}  r={metrics[model]['pattern_correlation']:.5f}")

# All three share one mask, so all three see the same WOA23 mean.
assert len({round(metrics[m]["common_mask_mean_woa23_degC"], 10) for m in MODELS}) == 1
assert len({metrics[m]["n_valid_cells"] for m in MODELS}) == 1


# ---------------------------------------------------------------------------
# Step 6. Ranking by the magnitude of the bias
# ---------------------------------------------------------------------------
ranking = sorted(MODELS, key=lambda m: abs(metrics[m]["delta_sst_degC"]))
print(f"[6]   ranking by |dSST|, smallest first: {' < '.join(ranking)}")


# ---------------------------------------------------------------------------
# Figure 1. 2x2 global maps - three bias fields plus the WOA23 reference
# ---------------------------------------------------------------------------
LAND = "#d9d9d6"        # neutral land fill, recessive against both colour scales
INK = "#2f3437"

lon = woa["lon"].values
lat = woa["lat"].values
land_layer = xr.where(landsea == 1, 1.0, np.nan)

# One symmetric scale shared by the three bias panels, so the panels are
# comparable; rounded up from the 99th percentile of |bias| over all models.
all_bias = np.concatenate([np.abs(bias_field[m].values[np.isfinite(bias_field[m].values)])
                           for m in MODELS])
vmax = float(np.ceil(np.percentile(all_bias, 99) * 2) / 2)
print(f"[fig1] shared diverging scale +/-{vmax:.1f} degC "
      f"(99th percentile of |bias| = {np.percentile(all_bias, 99):.2f})")

fig, axes = plt.subplots(2, 2, figsize=(15.0, 8.8))
# Manual placement rather than constrained_layout: the two colour bars below
# cover different panel sets, and an automatic layout runs them into each other.
fig.subplots_adjust(left=0.055, right=0.985, bottom=0.145, top=0.855, wspace=0.14, hspace=0.42)
panels = [(MODELS[0], axes[0, 0]), (MODELS[1], axes[0, 1]), (MODELS[2], axes[1, 0])]
bias_norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

for model, ax in panels:
    ax.pcolormesh(lon, lat, land_layer, cmap=matplotlib.colors.ListedColormap([LAND]),
                  shading="auto", zorder=0)
    bias_mesh = ax.pcolormesh(lon, lat, bias_field[model].values, cmap="RdBu_r",
                              norm=bias_norm, shading="auto", zorder=1)
    ax.set_title(f"{model} \u2212 WOA23\n"
                 f"\u0394SST = {metrics[model]['delta_sst_degC']:+.2f} \u00b0C,   "
                 f"RMSE = {metrics[model]['rmse_degC']:.3f} \u00b0C,   "
                 f"r = {metrics[model]['pattern_correlation']:.3f}",
                 fontsize=11, color=INK)

ax = axes[1, 1]
ax.pcolormesh(lon, lat, land_layer, cmap=matplotlib.colors.ListedColormap([LAND]),
              shading="auto", zorder=0)
ref_mesh = ax.pcolormesh(lon, lat, woa_common.values, cmap="viridis", shading="auto", zorder=1)
ax.set_title("WOA23 reference\n"
             f"1971\u20132000 annual mean at 0 m   "
             f"(mask mean = {woa_common_mean:.2f} \u00b0C)", fontsize=11, color=INK)

for ax in axes.flat:
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xticks([-180, -120, -60, 0, 60, 120, 180])
    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
    ax.set_xlabel("Longitude (\u00b0E)", color=INK)
    ax.set_ylabel("Latitude (\u00b0N)", color=INK)
    ax.tick_params(colors=INK, labelsize=9)
    for spine in ax.spines.values():
        spine.set_color("#c8ccce")

# One colour bar for the three bias panels, because they share one scale, and a
# separate one for the reference panel, which does not.
cax_bias = fig.add_axes([0.055, 0.065, 0.50, 0.020])
cb = fig.colorbar(bias_mesh, cax=cax_bias, orientation="horizontal", extend="both")
cb.set_label("SST bias, model \u2212 WOA23 (\u00b0C) \u2014 one shared scale for the three model panels",
             color=INK, fontsize=10)
cb.ax.tick_params(colors=INK, labelsize=9)

cax_ref = fig.add_axes([0.70, 0.065, 0.285, 0.020])
cb_ref = fig.colorbar(ref_mesh, cax=cax_ref, orientation="horizontal")
cb_ref.set_label("WOA23 SST (\u00b0C)", color=INK, fontsize=10)
cb_ref.ax.tick_params(colors=INK, labelsize=9)

fig.suptitle("CMIP6 historical SST climatology bias against WOA23, 1971\u20132000\n"
             "conservatively regridded to the WOA23 1\u00b0 grid, "
             f"common ocean mask of {n_common} cells",
             fontsize=13, color=INK, y=0.975)
fig.savefig(OUT / "fig1_bias_maps.png", dpi=160, facecolor="white")
plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2. Global mean SST of steps 2 and 3
# ---------------------------------------------------------------------------
labels = MODELS + ["WOA23"]
values = [native_mean[m] for m in MODELS] + [woa_mean]
# Models share one hue (identity: "CMIP6 model"); the observational reference is
# the one entity that differs, so it carries the second hue.
colors = ["#4b6bb5", "#4b6bb5", "#4b6bb5", "#b5643c"]

fig, ax = plt.subplots(figsize=(8.0, 5.2), constrained_layout=True)
bars = ax.bar(labels, values, color=colors, width=0.62, zorder=3)
for bar, value in zip(bars, values):
    # A white pad lets the reference line below break around the label instead
    # of running through it; the label sits above the line in the stacking order.
    ax.text(bar.get_x() + bar.get_width() / 2, value + 0.10, f"{value:.2f}",
            ha="center", va="bottom", fontsize=12, color=INK, zorder=5,
            bbox=dict(facecolor="white", edgecolor="none", boxstyle="square,pad=0.15"))

# The four bars span barely 1 degC on a zero-based axis, so the observed value is
# drawn across them: it is what makes the sign of each model's offset legible
# without truncating the axis, which would exaggerate the differences.
reference_line = ax.axhline(
    woa_mean, color="#b5643c", linestyle="--", linewidth=1.4, zorder=4,
    label=f"WOA23 observed reference, {woa_mean:.2f} \u00b0C")

ax.set_ylabel("Global mean SST (\u00b0C)", color=INK)
ax.set_ylim(0, max(values) + 3.4)
ax.set_title("Global mean SST climatology, 1971\u20132000\n"
             "each product area-weighted on its own native ocean grid",
             fontsize=12, color=INK)
ax.legend(handles=[reference_line], loc="upper left", frameon=False,
          fontsize=10, labelcolor=INK)
ax.grid(axis="y", color="#e4e7e9", linewidth=0.8, zorder=0)
ax.set_axisbelow(True)
ax.tick_params(colors=INK)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
for side in ("left", "bottom"):
    ax.spines[side].set_color("#c8ccce")

fig.savefig(OUT / "fig2_global_mean_bars.png", dpi=160, facecolor="white")
plt.close(fig)


# ---------------------------------------------------------------------------
# Numeric answer
# ---------------------------------------------------------------------------
answer = {
    "case_id": "T1",
    "title": "Multi-model bias of the global sea surface temperature climatology",
    "period": "1971-01 to 2000-12",
    "models": MODELS,
    "reference": "WOA23 temperature, 1971-2000 climate normal, annual, objectively analysed, 0 m",

    "native_global_mean_sst_degC": {m: round(native_mean[m], 2) for m in MODELS},
    "woa23_global_mean_sst_degC": round(woa_mean, 2),

    "common_mask": {
        "definition": "WOA23 landsea_01 ocean (value > 1) AND valid WOA23 0 m value AND "
                      "valid conservatively regridded climatology in all three models",
        "grid": "WOA23 1-degree",
        "n_cells": n_common,
    },
    "common_mask_mean_sst_degC": {
        **{m: round(metrics[m]["common_mask_mean_model_degC"], 2) for m in MODELS},
        "WOA23": round(woa_common_mean, 2),
    },
    "delta_sst_degC": {m: round(metrics[m]["delta_sst_degC"], 2) for m in MODELS},
    "rmse_degC": {m: round(metrics[m]["rmse_degC"], 3) for m in MODELS},
    "pattern_correlation": {m: round(metrics[m]["pattern_correlation"], 3) for m in MODELS},
    "ranking_by_abs_delta_sst": ranking,

    "full_precision": {m: metrics[m] for m in MODELS},
    "full_precision_native_global_mean_sst_degC": native_mean,
    "full_precision_woa23_global_mean_sst_degC": woa_mean,

    "figures": ["fig1_bias_maps.png", "fig2_global_mean_bars.png"],
    "method_notes": [
        "Model climatology: the day-weighted mean of all 360 monthly fields over the window, "
        "on the model's native grid (temporal_aggregate freq='mean').",
        "Native global means are area-weighted with each model's published areacello and the "
        "WOA23 value with the 1-degree cell areas, so the four numbers cover each product's "
        "own ocean mask and their differences are NOT delta_SST.",
        "delta_SST, RMSE and the pattern correlation are all area-weighted on the common mask, "
        "so delta_SST is exactly the difference of the two common-mask means.",
        "Regridding is conservative onto the WOA23 1-degree grid; RMSE and the pattern "
        "correlation depend on that choice.",
    ],
    "tools_used": [
        "load_regular_data", "load_climatology_data", "load_landsea_mask",
        "temporal_aggregate", "compute_basic_stat", "regrid_horizontal",
        "build_comparison_mask",
        "compute_difference_metrics", "analyze_correlation",
    ],
}

(OUT / "answer_T1.json").write_text(json.dumps(answer, indent=2) + "\n")


md = [
    "# Reference answer - T1: Multi-model bias of the global SST climatology",
    "",
    f"Period 1971-01 to 2000-12. Reference: {answer['reference']}.",
    f"Common comparison mask: {n_common} cells on the WOA23 1-degree grid.",
    "",
    "## Step 2-3. Global mean SST, each product on its own native ocean grid",
    "",
    "| Product | global mean SST (degC) |",
    "|---|---|",
]
for label, value in zip(labels, values):
    md.append(f"| {label} | {value:.2f} |")

md += [
    "",
    "These four numbers each cover that product's own ocean mask, so their differences are NOT",
    "delta_SST. delta_SST is defined on the common mask below.",
    "",
    "## Step 5-6. On the common comparison mask",
    "",
    "The mean SST column is given to 3 decimals so that delta_SST can be checked against it. The",
    "graded value of delta_SST is the 2-decimal one, and of RMSE and the correlation the 3-decimal",
    "one; rounding the means to 2 decimals first would not reproduce them.",
    "",
    "| Model | mean SST (degC) | delta_SST (degC) | RMSE (degC) | pattern correlation |",
    "|---|---|---|---|---|",
]
for m in MODELS:
    md.append(f"| {m} | {metrics[m]['common_mask_mean_model_degC']:.3f} | "
              f"{metrics[m]['delta_sst_degC']:+.2f} | {metrics[m]['rmse_degC']:.3f} | "
              f"{metrics[m]['pattern_correlation']:.3f} |")
md.append(f"| WOA23 (reference) | {woa_common_mean:.3f} | - | - | - |")
md += [
    "",
    f"Ranking by |delta_SST|, smallest first: **{' < '.join(ranking)}**.",
    "",
    "delta_SST at full precision: "
    + ", ".join(f"{m} {metrics[m]['delta_sst_degC']:+.4f}" for m in MODELS)
    + ". GFDL-ESM4's value sits close to the -0.085 rounding boundary, so grade it against the",
    "full-precision number in `answer_T1.json` rather than the rounded string.",
    "",
    "## Figures",
    "",
    "- `fig1_bias_maps.png` - 2x2 global maps: three bias fields on one shared "
    f"diverging scale (+/-{vmax:.1f} degC) plus the WOA23 reference field.",
    "- `fig2_global_mean_bars.png` - bar chart of the four native-grid global means, "
    "each annotated with its value.",
]
(OUT / "answer_T1.md").write_text("\n".join(md) + "\n")

print("\nWrote answer_T1.json, answer_T1.md, fig1_bias_maps.png, fig2_global_mean_bars.png")
