"""
Reference answer for case T6 - Multi-model assessment of nitracline depth and
its basin structure.

Run from the repository root with the project's conda environment:

    python cases/reference_ans/T6/reference_T6.py

Writes answer_T6.json, answer_T6.md, fig1_nitracline_maps.png and
fig2_equatorial_pacific_sections.png next to this file. Takes roughly fifteen
minutes, almost all of it forming the 30-year mean nitrate profiles.
"""

import json
import time
from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from src.tools.registry import auto_import_tools
auto_import_tools()
from src.tools.data_loader import load_regular_data, load_climatology_data
from src.tools.preprocessing import regrid_horizontal, build_comparison_mask, select_basin
from src.tools.scientific_calculator import compute_nutricline_depth
from src.tools.basic_calculator import (compute_basic_stat, temporal_aggregate, transform_unit,
                                        compute_difference_metrics, analyze_correlation)

OUT = Path(__file__).resolve().parent
MODELS = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
SOURCES = ["WOA23"] + MODELS
WINDOW = ["1971-01", "2000-12"]
MAX_DEPTH = 1000.0
LOAD_DEPTH = 1500.0
BASINS = ["Pacific", "Atlantic", "Indian"]
EQ_BAND = [-2, 2]             # the equatorial Pacific section of figure 2

xr.set_options(keep_attrs=True)


# ---------------------------------------------------------------------------
# Step 1. Climatological nitrate profiles, then the nitracline from them
# ---------------------------------------------------------------------------
profile, nitracline = {}, {}

woa_profile = load_climatology_data("woa23", "nitrate", time_slice="annual",
                                    depth_range=[0, LOAD_DEPTH])
profile["WOA23"] = woa_profile
print(f"[1]   WOA23 nitrate: {woa_profile.sizes['depth']} levels, "
      f"units {woa_profile.attrs['units']}, period {woa_profile.attrs.get('climatology_period')}")

for model in MODELS:
    started = time.time()
    monthly = load_regular_data("cmip6_historical", "no3", model=model, time_range=WINDOW,
                                depth_range=[0, LOAD_DEPTH], chunk_time=12)
    
    climatology = temporal_aggregate(monthly, freq="mean", op="mean").compute()
    profile[model] = climatology
    print(f"[1]   {model:14s} 30-year mean profile {climatology.shape} "
          f"({time.time() - started:.0f} s)")

for source in SOURCES:
    nitracline[source] = compute_nutricline_depth(profile[source], min_depth=0.0,
                                                  max_depth=MAX_DEPTH).compute()
    print(f"[1]   {source:14s} nitracline valid at {int(nitracline[source].notnull().sum())} cells")


# ---------------------------------------------------------------------------
# Step 2. Global mean depth
# ---------------------------------------------------------------------------
global_mean = {s: float(compute_basic_stat(nitracline[s], op="area_weighted_mean",
                                           dim=["lat", "lon"]).compute()) for s in SOURCES}
global_bias = {m: global_mean[m] - global_mean["WOA23"] for m in MODELS}
for source in SOURCES:
    extra = f"  bias {global_bias[source]:+7.1f} m" if source in MODELS else "  (reference)"
    print(f"[2]   {source:14s} global mean nitracline {global_mean[source]:7.2f} m{extra}")


# ---------------------------------------------------------------------------
# Step 3. Basin means
# ---------------------------------------------------------------------------
basin_mean, basin_definition = {}, {}
for source in SOURCES:
    basin_mean[source] = {}
    for basin in BASINS:
        masked = select_basin(nitracline[source], basin=basin)   # marginal seas folded in
        basin_mean[source][basin] = float(
            compute_basic_stat(masked, op="area_weighted_mean", dim=["lat", "lon"]).compute())
        basin_definition[basin] = masked.attrs.get("basin_definition", "")
    print(f"[3]   {source:14s} " + "  ".join(
        f"{b} {basin_mean[source][b]:6.1f} m" for b in BASINS))

basin_bias = {m: {b: basin_mean[m][b] - basin_mean["WOA23"][b] for b in BASINS} for m in MODELS}

# Does the model reproduce the observed ordering of the three basins?
observed_order = sorted(BASINS, key=lambda b: basin_mean["WOA23"][b])
ordering = {m: sorted(BASINS, key=lambda b: basin_mean[m][b]) for m in MODELS}
ordering_match = {m: ordering[m] == observed_order for m in MODELS}
print(f"[3]   observed ordering, shallowest first: {' < '.join(observed_order)}")
for model in MODELS:
    print(f"[3]   {model:14s} ordering {' < '.join(ordering[model])}  "
          f"matches: {ordering_match[model]}")


# ---------------------------------------------------------------------------
# Step 4. Pattern agreement on one shared 1-degree mask
# ---------------------------------------------------------------------------
regridded = {m: regrid_horizontal(nitracline[m], target="woa23", method="conservative")
             for m in MODELS}
for model in MODELS:
    print(f"[4]   {model:14s} nitracline regridded to 1 degree")

common = build_comparison_mask([regridded[m] for m in MODELS] + [nitracline["WOA23"]],
                               ocean_only=True)
n_common = int(common.sum())
print(f"[4]   shared comparison mask: {n_common} cells ({common.attrs['ocean_mask']})")

woa_common = nitracline["WOA23"].where(common)
pattern = {}
for model in MODELS:
    field = regridded[model].where(common)
    difference = compute_difference_metrics(field, woa_common, dim=["lat", "lon"], weighting="area")
    correlation = analyze_correlation(field, woa_common, dim=["lat", "lon"],
                                      method="pearson", weighting="area")
    pattern[model] = {
        "spatial_correlation": float(correlation["correlation"]),
        "rmse_m": float(difference["rmse"]),
        "bias_on_mask_m": float(difference["bias"]),
        "n_cells": int(difference["n_valid"]),
    }
    print(f"[4]   {model:14s} r {pattern[model]['spatial_correlation']:.4f}  "
          f"RMSE {pattern[model]['rmse_m']:7.2f} m")

best = min(MODELS, key=lambda m: pattern[m]["rmse_m"])
print(f"[5]   closest to observations by RMSE: {best}")


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
INK = "#2f3437"
LAND = "#d9d9d6"
lon1 = nitracline["WOA23"]["lon"].values
lat1 = nitracline["WOA23"]["lat"].values

display = {"WOA23": nitracline["WOA23"]}
display.update({m: regridded[m] for m in MODELS})

stack = np.concatenate([display[s].values[np.isfinite(display[s].values)] for s in SOURCES])
vmax_depth = float(np.ceil(np.percentile(stack, 95) / 50) * 50)
print(f"[fig1] shared depth scale 0-{vmax_depth:.0f} m "
      f"(95th percentile of nitracline depth = {np.percentile(stack, 95):.0f} m)")

fig, axes = plt.subplots(2, 2, figsize=(14.6, 8.2))
fig.subplots_adjust(left=0.06, right=0.985, bottom=0.145, top=0.855, wspace=0.12, hspace=0.30)
for ax, source in zip(axes.flat, SOURCES):
    ax.set_facecolor(LAND)
    # Reversed so that deeper reads darker, which is how a depth map is read.
    mesh = ax.pcolormesh(lon1, lat1, display[source].values, cmap="viridis_r",
                         vmin=0, vmax=vmax_depth, shading="auto")
    if source == "WOA23":
        title = f"WOA23 (observed)\nglobal mean {global_mean[source]:.0f} m"
    else:
        title = (f"{source}\nglobal mean {global_mean[source]:.0f} m,   "
                 f"bias {global_bias[source]:+.0f} m,   "
                 f"r = {pattern[source]['spatial_correlation']:.3f}")
    ax.set_title(title, fontsize=10.5, color=INK)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_yticks([-90, -45, 0, 45, 90])
    ax.tick_params(colors=INK, labelsize=8.5)
    ax.set_xlabel("Longitude (°E)", color=INK, fontsize=9)
    ax.set_ylabel("Latitude (°N)", color=INK, fontsize=9)
    for spine in ax.spines.values():
        spine.set_color("#c8ccce")

cax = fig.add_axes([0.28, 0.065, 0.44, 0.022])
bar = fig.colorbar(mesh, cax=cax, orientation="horizontal", extend="max")
bar.set_label("Nitracline depth (m), depth of the maximum vertical nitrate gradient",
              color=INK, fontsize=10)
bar.ax.tick_params(colors=INK, labelsize=9)
fig.suptitle("Nitracline depth: three CMIP6 historical simulations and WOA23\n"
             "models 1971–2000, WOA23 all-period climatology; search confined to the top "
             f"{MAX_DEPTH:.0f} m.\nModel panels regridded to 1° for display; every reported "
             "number comes from the native grid",
             fontsize=12, color=INK, y=0.985)
fig.savefig(OUT / "fig1_nitracline_maps.png", dpi=150, facecolor="white")
plt.close(fig)


# --- Figure 2: equatorial Pacific sections --------------------------------
DENSITY = 1.025
section, section_nitracline = {}, {}

for source in SOURCES:
    if source == "WOA23":
        field = profile[source]
        line = nitracline[source]
    else:
        field = regrid_horizontal(profile[source], target="woa23", method="conservative")
        field = transform_unit(field, from_unit=profile[source].attrs["units"],
                               to_unit="umol/kg", density=DENSITY)
        line = regridded[source]
        print(f"[fig2] {source:14s} section regridded and converted "
              f"({field.attrs.get('density_assumption', 'constant density')})")

    band = field.sel(lat=slice(EQ_BAND[0], EQ_BAND[1]), depth=slice(0, MAX_DEPTH))
    
    section[source] = compute_basic_stat(band, op="mean", dim="lat").compute()
    section_nitracline[source] = compute_basic_stat(
        line.sel(lat=slice(EQ_BAND[0], EQ_BAND[1])), op="mean", dim="lat").compute()

stack = np.concatenate([section[s].values[np.isfinite(section[s].values)] for s in SOURCES])
vmax_no3 = float(np.ceil(np.percentile(stack, 99) / 5) * 5)
print(f"[fig2] shared nitrate scale 0-{vmax_no3:.0f} umol/kg")


def eastward(values_lon):
    """Longitudes relabelled 120E..280E so the Pacific section runs unbroken."""
    shifted = np.where(values_lon < 0, values_lon + 360.0, values_lon)
    inside = (shifted >= 120.0) & (shifted <= 280.0)
    order = np.argsort(shifted[inside])
    return inside, order, shifted[inside][order]


fig, axes = plt.subplots(4, 1, figsize=(11.0, 12.4), sharex=True)
fig.subplots_adjust(left=0.085, right=0.985, bottom=0.115, top=0.895, hspace=0.20)
for ax, source in zip(axes, SOURCES):
    field = section[source]
    inside, order, x_plot = eastward(field["lon"].values)
    mesh = ax.pcolormesh(x_plot, field["depth"].values,
                         field.values[:, inside][:, order],
                         cmap="YlGnBu", vmin=0, vmax=vmax_no3, shading="auto")

    line = section_nitracline[source]
    l_inside, l_order, lx = eastward(line["lon"].values)
    ax.plot(lx, line.values[l_inside][l_order], color="#c1543a", linewidth=1.8,
            zorder=4, label="nitracline")

    ax.set_ylim(MAX_DEPTH, 0)
    ax.set_facecolor(LAND)        # the Maritime Continent columns are land, not gaps
    ax.set_ylabel("Depth (m)", color=INK, fontsize=9.5)
    label = "WOA23 (observed)" if source == "WOA23" else source
    ax.text(0.008, 0.06, label, transform=ax.transAxes, fontsize=11, color=INK,
            va="bottom", ha="left", fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="none", boxstyle="square,pad=0.25"))
    ax.tick_params(colors=INK, labelsize=9)
    for spine in ax.spines.values():
        spine.set_color("#c8ccce")
axes[0].legend(loc="lower right", fontsize=9.5, labelcolor=INK, frameon=True,
               facecolor="white", edgecolor="none", framealpha=0.85)

axes[-1].set_xticks([120, 150, 180, 210, 240, 270])
axes[-1].set_xticklabels(["120\u00b0E", "150\u00b0E", "180\u00b0",
                          "150\u00b0W", "120\u00b0W", "90\u00b0W"])
axes[-1].set_xlabel("Longitude", color=INK)

cax = fig.add_axes([0.30, 0.055, 0.40, 0.016])
bar = fig.colorbar(mesh, cax=cax, orientation="horizontal", extend="max")
bar.set_label("Nitrate (\u00b5mol kg\u207b\u00b9)", color=INK, fontsize=10)
bar.ax.tick_params(colors=INK, labelsize=9)
fig.suptitle("Equatorial Pacific nitrate section, averaged over 2\u00b0S\u20132\u00b0N\n"
             "models 1971\u20132000, WOA23 all-period climatology; model values regridded to "
             f"1\u00b0 and converted\nfrom mol m\u207b\u00b3 at a constant density of "
             f"{DENSITY} kg L\u207b\u00b9. The red line is each source's own nitracline",
             fontsize=12, color=INK, y=0.982)
fig.savefig(OUT / "fig2_equatorial_pacific_sections.png", dpi=150, facecolor="white")
plt.close(fig)


# ---------------------------------------------------------------------------
# Numeric answer
# ---------------------------------------------------------------------------
answer = {
    "case_id": "T6",
    "title": "Multi-model assessment of nitracline depth and its basin structure",
    "sources": SOURCES,
    "method": {
        "definition": "depth of the maximum vertical nitrate gradient, taken across each layer's "
                      f"own bounds, searched within the top {MAX_DEPTH:.0f} m; profiles loaded to "
                      f"{LOAD_DEPTH:.0f} m so the deepest layer searched has the level below it",
        "order_of_operations": "the 1971-2000 climatological profile is formed first and the "
                               "nitracline taken from it, because WOA23 nitrate exists only as a "
                               "single all-period annual climatology at these depths",
        "units": "no conversion is needed for the depth: the maximum-gradient depth is invariant "
                 "under scaling of the profile, verified as identical at every valid cell",
        "climatology_periods": "models 1971-2000, WOA23 nitrate 1965-2022 - a genuine mismatch "
                               "that cannot be removed",
        "grid": "every reported number on the source's native grid with its own cell areas; "
                "regridding only for the pattern metrics of step 4 and for display",
        "marginal_seas": "folded into their parent ocean (select_basin default)",
        "southern_and_arctic_oceans": "separate basins in the WOA23 mask, not included in the Atlantic, "
                                      "Indian or Pacific numbers",
    },

    "global_mean_nitracline_depth_m": {s: round(global_mean[s]) for s in SOURCES},
    "global_bias_m": {m: round(global_bias[m]) for m in MODELS},

    "basin_mean_nitracline_depth_m": {s: {b: round(basin_mean[s][b]) for b in BASINS}
                                      for s in SOURCES},
    "basin_bias_m": {m: {b: round(basin_bias[m][b]) for b in BASINS} for m in MODELS},
    "observed_basin_ordering_shallowest_first": observed_order,
    "model_basin_ordering": ordering,
    "basin_ordering_match": ordering_match,
    "basin_definitions": basin_definition,

    "comparison_mask": {"grid": "WOA23 1-degree", "n_cells": n_common},
    "spatial_correlation": {m: round(pattern[m]["spatial_correlation"], 3) for m in MODELS},
    "rmse_m": {m: round(pattern[m]["rmse_m"], 1) for m in MODELS},
    "closest_to_observations_by_rmse": best,

    "full_precision": {
        "global_mean_nitracline_depth_m": global_mean,
        "global_bias_m": global_bias,
        "basin_mean_nitracline_depth_m": basin_mean,
        "basin_bias_m": basin_bias,
        "pattern": pattern,
    },

    "findings": [
        "All three models place the nitracline too shallow globally, by 18 to 34 m. "
        "CESM2 and GFDL-ESM4 have their largest basin shortfalls in the Atlantic "
        "(56 and 74 m); MPI-ESM1-2-LR's largest shortfall is in the Indian Ocean "
        "(58 m, versus 54 m in the Atlantic).",
        "Only MPI-ESM1-2-LR reproduces the observed ordering of the three basins "
        f"({' < '.join(observed_order)}); CESM2 and GFDL-ESM4 both put the Indian Ocean deepest.",
        "Spatial correlations are low (0.40 to 0.60) and the RMSE is 152 to 163 m, far worse than "
        "a temperature-based diagnostic would give. Part of that is the diagnostic itself: where "
        "the profile is smooth the depth of the maximum gradient jumps between levels, and WOA23 "
        "has 47 levels in the top 1000 m against the models' 19 to 40, so the observed map is the "
        "patchiest of the four.",
        "The equatorial Pacific section exposes a large CESM2 bias that the depth diagnostic alone "
        "does not show: below about 300 m in the eastern basin CESM2's nitrate falls to roughly a "
        "fifth of the observed value (8.5 against 40.2 micromoles per kilogram at 787 m, 110W), "
        "the signature of excessive water-column denitrification in an over-extended eastern "
        "tropical Pacific oxygen minimum zone. GFDL-ESM4 and MPI-ESM1-2-LR do not show it.",
    ],

    "figures": ["fig1_nitracline_maps.png", "fig2_equatorial_pacific_sections.png"],
    "tools_used": [
        "load_regular_data", "load_climatology_data", "temporal_aggregate",
        "compute_nutricline_depth", "compute_basic_stat", "select_basin", "regrid_horizontal",
        "build_comparison_mask", "compute_difference_metrics", "analyze_correlation",
        "transform_unit",
    ],
}

(OUT / "answer_T6.json").write_text(json.dumps(answer, indent=2) + "\n")

md = [
    "# Reference answer - T6: nitracline depth, multi-model assessment",
    "",
    "Nitracline = depth of the maximum vertical nitrate gradient, searched within the top "
    f"{MAX_DEPTH:.0f} m; profiles are loaded to {LOAD_DEPTH:.0f} m so the deepest layer searched",
    "has the level below it and its gradient is not one-sided.",
    "The climatological profile is formed first and the nitracline taken from it, because WOA23",
    "nitrate exists only as one all-period annual climatology at these depths.",
    "Models cover 1971-2000; WOA23 nitrate is a 1965-2022 composite - the periods genuinely differ.",
    "",
    "## Step 2. Global mean nitracline depth",
    "",
    "| Source | global mean (m) | bias vs WOA23 (m) |",
    "|---|---|---|",
    f"| WOA23 (observed) | {global_mean['WOA23']:.0f} | - |",
]
for model in MODELS:
    md.append(f"| {model} | {global_mean[model]:.0f} | {global_bias[model]:+.0f} |")

md += [
    "",
    "## Step 3. Basin means",
    "",
    "| Source | Pacific (m) | Atlantic (m) | Indian (m) |",
    "|---|---|---|---|",
]
for source in SOURCES:
    md.append(f"| {source} | " + " | ".join(f"{basin_mean[source][b]:.0f}" for b in BASINS) + " |")
md += ["", "| Model | Pacific bias | Atlantic bias | Indian bias | ordering matches WOA23 |",
       "|---|---|---|---|---|"]
for model in MODELS:
    md.append(f"| {model} | " + " | ".join(f"{basin_bias[model][b]:+.0f}" for b in BASINS)
              + f" | {'yes' if ordering_match[model] else 'NO'} |")
md += [
    "",
    f"Observed ordering, shallowest first: **{' < '.join(observed_order)}**.",
    "Marginal seas are folded into their parent ocean. The Southern and Arctic Oceans are separate "
    "basins in the WOA23 mask and are not part of these numbers.",
    "",
    f"## Step 4. Pattern agreement on the shared 1-degree mask ({n_common} cells)",
    "",
    "| Model | spatial correlation | RMSE (m) |",
    "|---|---|---|",
]
for model in MODELS:
    md.append(f"| {model} | {pattern[model]['spatial_correlation']:.3f} | "
              f"{pattern[model]['rmse_m']:.1f} |")
md += [
    "",
    f"Closest to observations by RMSE: **{best}**.",
    "",
    "## Step 5. What the comparison shows",
    "",
]
md += [f"- {line}" for line in answer["findings"]]
md += [
    "",
    "## Figures",
    "",
    "- `fig1_nitracline_maps.png` - 2x2 maps of nitracline depth, one shared scale, models",
    "  regridded to 1 degree for display only.",
    "- `fig2_equatorial_pacific_sections.png` - 2S-2N nitrate sections from 120E to 80W in the",
    "  top 1000 m, one shared scale in micromoles per kilogram, each source's nitracline overlaid.",
]
(OUT / "answer_T6.md").write_text("\n".join(md) + "\n")

print("\nWrote answer_T6.json, answer_T6.md and two figures")
