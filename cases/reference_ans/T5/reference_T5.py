"""
Reference answer for case T5 - Sea surface salinity: mean-state bias,
historical trend and ssp245 scenario change.

Run from the repository root with the project's conda environment:

    source /home/lorraine/miniconda3/bin/activate base && conda activate jingwei_QA
    PYTHONPATH=. python cases/reference_ans/T5/reference_T5.py

Writes answer_T5.json, answer_T5.md and three figures next to this file.
Takes roughly half an hour, most of it reading the surface level out of the
monthly salinity archives.

Every physical quantity is produced by a tool from src/tools. The annual SSS
field of each model and experiment is computed once and held in memory, because
every later step is derived from it and the read is what costs the time.
"""

import json
import time
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
from src.tools.preprocessing import regrid_horizontal, build_comparison_mask, select_basin
from src.tools.scientific_calculator import compute_mixed_layer_depth  # noqa: F401  (tool registry)
from src.tools.basic_calculator import (compute_basic_stat, temporal_aggregate, fit_1D_trend,
                                        fit_trend_map, compute_difference_metrics,
                                        analyze_correlation)

OUT = Path(__file__).resolve().parent
MODELS = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
SURFACE = [0, 10]                      # metres, the 0-10 m layer that defines SSS
HIST_WINDOW = ["1950-01", "2014-12"]
SSP_WINDOW = ["2071-01", "2100-12"]
CLIM_YEARS = (1971, 2000)              # the WOA23 normal, for the mean-state comparison
BASE_YEARS = (1995, 2014)              # the historical baseline for the scenario change
BASINS = ["Atlantic", "Indian", "Pacific"]

xr.set_options(keep_attrs=True)


def annual_sss(dataset, model, window):
    """Annual 0-10 m volume-weighted mean salinity fields, on the native grid."""
    salinity = load_regular_data(dataset, "so", model=model, time_range=window,
                                 depth_range=SURFACE, chunk_time=12)
    # Volume weighting, so that unequal layer thicknesses do not count alike;
    # with a single level in the layer it reduces to that level, as it should.
    surface = compute_basic_stat(salinity, op="volume_weighted_mean", dim="depth")
    annual = temporal_aggregate(surface, freq="year", op="mean")
    incomplete = annual.attrs.get("incomplete_periods")
    assert not incomplete, f"{model} {dataset}: {incomplete}"
    return annual.compute()


# ---------------------------------------------------------------------------
# Step 1. The annual SSS fields, computed once
# ---------------------------------------------------------------------------
fields = {}
for model in MODELS:
    for dataset, window, label in (("cmip6_historical", HIST_WINDOW, "historical"),
                                   ("cmip6_ssp245", SSP_WINDOW, "ssp245")):
        started = time.time()
        fields[(model, label)] = annual_sss(dataset, model, window)
        print(f"[1]   {model:14s} {label:10s} {fields[(model, label)].sizes['year']:3d} annual fields "
              f"({time.time() - started:.0f} s)")


# ---------------------------------------------------------------------------
# Step 2. Mean-state bias against WOA23
# ---------------------------------------------------------------------------
woa_profile = load_climatology_data("woa23", "salinity", time_slice="annual",
                                    period="1971-2000", depth_range=SURFACE)
woa_sss = compute_basic_stat(woa_profile, op="volume_weighted_mean", dim="depth").compute()
woa_global = float(compute_basic_stat(woa_sss, op="area_weighted_mean", dim=["lat", "lon"]).compute())
print(f"[2]   WOA23 global mean SSS = {woa_global:.4f} psu")

model_clim, global_mean, mean_bias = {}, {}, {}
for model in MODELS:
    historical = fields[(model, "historical")]
    window = historical.sel(year=slice(*CLIM_YEARS))
    assert window.sizes["year"] == 30, f"{model}: {window.sizes['year']} years in the normal"
    clim = window.mean("year")
    clim.attrs = dict(historical.attrs)
    model_clim[model] = clim
    global_mean[model] = float(
        compute_basic_stat(clim, op="area_weighted_mean", dim=["lat", "lon"]).compute())
    mean_bias[model] = global_mean[model] - woa_global
    print(f"[2]   {model:14s} global mean SSS = {global_mean[model]:.4f} psu  "
          f"bias {mean_bias[model]:+.4f} psu")


# ---------------------------------------------------------------------------
# Step 3. Pattern agreement on one shared 1-degree mask
# ---------------------------------------------------------------------------
regridded = {m: regrid_horizontal(model_clim[m], target="woa23", method="conservative")
             for m in MODELS}
for model in MODELS:
    print(f"[3]   {model:14s} regridded to 1 degree")

common = build_comparison_mask([regridded[m] for m in MODELS] + [woa_sss], ocean_only=True)
n_common = int(common.sum())
print(f"[3]   shared comparison mask: {n_common} cells ({common.attrs['ocean_mask']})")

woa_common = woa_sss.where(common)
pattern, bias_field = {}, {}
for model in MODELS:
    field = regridded[model].where(common)
    difference = compute_difference_metrics(field, woa_common, dim=["lat", "lon"], weighting="area")
    correlation = analyze_correlation(field, woa_common, dim=["lat", "lon"],
                                      method="pearson", weighting="area")
    pattern[model] = {
        "rmse_psu": float(difference["rmse"]),
        "spatial_correlation": float(correlation["correlation"]),
        "bias_on_mask_psu": float(difference["bias"]),
        "n_cells": int(difference["n_valid"]),
    }
    bias_field[model] = difference["difference"]
    print(f"[3]   {model:14s} RMSE {pattern[model]['rmse_psu']:.4f} psu  "
          f"r {pattern[model]['spatial_correlation']:.4f}")


# ---------------------------------------------------------------------------
# Step 4. Global mean trends, historical and ssp245
# ---------------------------------------------------------------------------
series, trends = {}, {}
for model in MODELS:
    for label in ("historical", "ssp245"):
        annual = fields[(model, label)]
        global_series = compute_basic_stat(annual, op="area_weighted_mean",
                                           dim=["lat", "lon"]).compute()
        series[(model, label)] = global_series
        fit = fit_1D_trend(global_series, method="linear")
        trend = fit["trend"]
        # fit_1D_trend reports per x unit, which is one year on this axis.
        trends[(model, label)] = {
            "slope_psu_per_decade": trend["slope"] * 10.0,
            "ci95_low_psu_per_decade": trend["ci95_low"] * 10.0,
            "ci95_high_psu_per_decade": trend["ci95_high"] * 10.0,
            "r2": fit["metrics"]["r2"],
            "p_value": trend["p_value"],
            "ci_excludes_zero": bool(trend["significant"]),
            "n": int(trend["n"]),
            "n_effective": float(trend["n_effective"]),
            "fitted": fit["fitted_data"],
        }
        entry = trends[(model, label)]
        print(f"[4]   {model:14s} {label:10s} {entry['slope_psu_per_decade']:+.5f} psu/decade  "
              f"CI [{entry['ci95_low_psu_per_decade']:+.5f}, {entry['ci95_high_psu_per_decade']:+.5f}]  "
              f"R2 {entry['r2']:.3f}  CI excludes zero: {entry['ci_excludes_zero']}")

trend_ratio = {}
for model in MODELS:
    historical = trends[(model, "historical")]
    scenario = trends[(model, "ssp245")]
    trend_ratio[model] = {
        "ratio": scenario["slope_psu_per_decade"] / historical["slope_psu_per_decade"],
        "denominator_significant": historical["ci_excludes_zero"],
        "interpretable": historical["ci_excludes_zero"],
    }
    note = ("" if historical["ci_excludes_zero"]
            else "  <- denominator not distinguishable from zero, not a meaningful multiple")
    print(f"[4]   {model:14s} ssp245/historical slope ratio = "
          f"{trend_ratio[model]['ratio']:.1f}{note}")


# ---------------------------------------------------------------------------
# Step 5. Scenario change by basin
# ---------------------------------------------------------------------------
# Which ensemble member each experiment is, read from the files. A scenario run
# continues the historical member named in its parent_variant_label; where the
# archive's historical run is a different member, the difference below mixes the
# forced change with the internal variability between two members.
members = {(m, label): fields[(m, label)].attrs.get("variant_label", "unknown")
           for m in MODELS for label in ("historical", "ssp245")}
mismatched = [m for m in MODELS if members[(m, "historical")] != members[(m, "ssp245")]]
member_note = (
    "Ensemble members: " + "; ".join(
        f"{m} historical {members[(m, 'historical')]}, ssp245 {members[(m, 'ssp245')]}" for m in MODELS)
    + ". " + (
        " and ".join(mismatched) + (" uses" if len(mismatched) == 1 else " use") + " different members in "
        "the two experiments (the ssp245 run continues a historical member that is not in the archive), "
        "so its delta SSS contains the "
        "difference in internal variability between two members as well as the forced change."
        if mismatched else "Every model continues the same member into ssp245.")
)
print(f"[5]   {member_note}")

delta_sss, delta_field = {}, {}
for model in MODELS:
    scenario = fields[(model, "ssp245")].mean("year")
    baseline = fields[(model, "historical")].sel(year=slice(*BASE_YEARS)).mean("year")
    difference = scenario - baseline
    difference.attrs = dict(fields[(model, "historical")].attrs)
    delta_field[model] = difference

    entry = {"global": float(compute_basic_stat(difference, op="area_weighted_mean",
                                                dim=["lat", "lon"]).compute())}
    for basin in BASINS:
        masked = select_basin(difference, basin=basin)          # marginal seas folded in
        entry[basin] = float(compute_basic_stat(masked, op="area_weighted_mean",
                                                dim=["lat", "lon"]).compute())
        entry[f"{basin}_definition"] = masked.attrs.get("basin_definition", "")
    delta_sss[model] = entry
    print(f"[5]   {model:14s} dSSS global {entry['global']:+.4f}  "
          + "  ".join(f"{b} {entry[b]:+.4f}" for b in BASINS))


# ---------------------------------------------------------------------------
# Trend maps for figure 2, fitted on the native grid then regridded for display
# ---------------------------------------------------------------------------
trend_maps = {}
for model in MODELS:
    fitted = fit_trend_map(fields[(model, "historical")], dim="year", alpha=0.05)
    slope = fitted["slope"]
    # fit_trend_map's per_decade applies to a date axis; this axis is integer
    # years, so the slope is per year and is converted here.
    slope = slope * 10.0
    slope.attrs = dict(fields[(model, "historical")].attrs)
    significant = fitted["significant"].astype("float64")
    significant.attrs = dict(fields[(model, "historical")].attrs)
    trend_maps[model] = {
        "slope": regrid_horizontal(slope, target="woa23", method="conservative"),
        # A 0/1 flag is categorical, so it is carried across by nearest cell.
        "significant": regrid_horizontal(significant, target="woa23", method="nearest_s2d"),
    }
    print(f"[fig2] {model:14s} trend map fitted and regridded for display")


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
INK = "#2f3437"
LAND = "#d9d9d6"
COLOUR = {"CESM2": "#0072B2", "GFDL-ESM4": "#009E73", "MPI-ESM1-2-LR": "#E69F00"}
lon = woa_sss["lon"].values
lat = woa_sss["lat"].values


def map_row(fig, axes, values_of, vmax, cmap, label, stipple_of=None):
    for ax, model in zip(axes, MODELS):
        ax.set_facecolor(LAND)
        mesh = ax.pcolormesh(lon, lat, values_of(model), cmap=cmap,
                             norm=TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax),
                             shading="auto")
        if stipple_of is not None:
            flag = np.nan_to_num(stipple_of(model)) > 0.5
            # Sample the significant cells on a regular lattice, so the stipple
            # reads as an even texture instead of the diagonal streaks a
            # flattened every-nth-element thinning produces.
            lattice = np.zeros_like(flag, dtype=bool)
            lattice[::3, ::3] = True
            ys, xs = np.where(flag & lattice)
            ax.plot(lon[xs], lat[ys], ".", color="#33383b",
                    markersize=0.5, alpha=0.55, linestyle="none")
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        ax.set_xticks([-180, -90, 0, 90, 180])
        ax.set_yticks([-90, -45, 0, 45, 90])
        ax.tick_params(colors=INK, labelsize=8.5)
        ax.set_xlabel("Longitude (°E)", color=INK, fontsize=9)
        for spine in ax.spines.values():
            spine.set_color("#c8ccce")
    axes[0].set_ylabel("Latitude (°N)", color=INK, fontsize=9)
    return mesh


# --- Figure 1: mean-state bias ---------------------------------------------
stack = np.concatenate([np.abs(bias_field[m].values[np.isfinite(bias_field[m].values)])
                        for m in MODELS])
vmax_bias = float(np.ceil(np.percentile(stack, 98) * 2) / 2)

fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.6))
fig.subplots_adjust(left=0.045, right=0.99, bottom=0.235, top=0.80, wspace=0.13)
mesh = map_row(fig, axes, lambda m: bias_field[m].values, vmax_bias, "RdBu_r", "psu")
for ax, model in zip(axes, MODELS):
    ax.set_title(f"{model}\nbias {pattern[model]['bias_on_mask_psu']:+.2f} psu,   "
                 f"RMSE {pattern[model]['rmse_psu']:.2f} psu,   "
                 f"r = {pattern[model]['spatial_correlation']:.3f}",
                 fontsize=10.5, color=INK)
cax = fig.add_axes([0.30, 0.10, 0.40, 0.035])
bar = fig.colorbar(mesh, cax=cax, orientation="horizontal", extend="both")
bar.set_label("SSS bias, model − WOA23 (psu)", color=INK, fontsize=10)
bar.ax.tick_params(colors=INK, labelsize=9)
fig.suptitle("Sea surface salinity bias against WOA23, 1971–2000 climatology\n"
             f"0–10 m volume-weighted; conservatively regridded to 1°; "
             f"shared mask of {n_common} cells",
             fontsize=12.5, color=INK, y=0.975)
fig.savefig(OUT / "fig1_sss_bias_maps.png", dpi=150, facecolor="white")
plt.close(fig)

# --- Figure 2: historical trend, with significance stippled ----------------
stack = np.concatenate([np.abs(trend_maps[m]["slope"].values[np.isfinite(trend_maps[m]["slope"].values)])
                        for m in MODELS])
percentiles = {q: float(np.percentile(stack, q)) for q in (50, 75, 90, 95, 98, 99)}
print("[fig2] |trend| percentiles (psu/decade): "
      + "  ".join(f"p{q}={v:.3f}" for q, v in percentiles.items()))
# The 90th percentile, rounded up to a readable step: the colour bar extends at
# both ends, so the coastal extremes are still shown, just saturated.
vmax_trend = float(np.ceil(percentiles[90] * 20) / 20)
print(f"[fig2] shared diverging scale +/-{vmax_trend:.2f} psu/decade "
      f"(90th percentile of |trend| = {percentiles[90]:.3f})")

fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.6))
fig.subplots_adjust(left=0.045, right=0.99, bottom=0.235, top=0.80, wspace=0.13)
mesh = map_row(fig, axes, lambda m: trend_maps[m]["slope"].values, vmax_trend, "BrBG_r", "psu/decade",
               stipple_of=lambda m: trend_maps[m]["significant"].values)
for ax, model in zip(axes, MODELS):
    ax.set_title(f"{model}", fontsize=11, color=INK)
cax = fig.add_axes([0.30, 0.10, 0.40, 0.035])
bar = fig.colorbar(mesh, cax=cax, orientation="horizontal", extend="both")
bar.set_label("SSS trend, historical 1950–2014 (psu decade⁻¹)", color=INK, fontsize=10)
bar.ax.tick_params(colors=INK, labelsize=9)
fig.suptitle("Historical sea surface salinity trend, 1950–2014\n"
             "fitted on annual means at every cell of the native grid, regridded to 1° for "
             "display; stippling marks p < 0.05 after the serial-correlation correction",
             fontsize=12.5, color=INK, y=0.975)
fig.savefig(OUT / "fig2_sss_trend_maps.png", dpi=150, facecolor="white")
plt.close(fig)

# --- Figure 3: global mean series, two segments, gap left empty ------------
fig, ax = plt.subplots(figsize=(13.0, 6.0), constrained_layout=True)
for model in MODELS:
    for label, style in (("historical", "-"), ("ssp245", "-")):
        values = series[(model, label)]
        years = values["year"].values
        ax.plot(years, values.values, style, color=COLOUR[model], linewidth=1.3, alpha=0.55,
                label=model if label == "historical" else None)
        fitted = trends[(model, label)]["fitted"]
        ax.plot(years, np.asarray(fitted), "--", color=COLOUR[model], linewidth=2.2, zorder=4)

ax.axvspan(2015, 2070, color="#f4f5f6", zorder=0)
low, high = ax.get_ylim()
ax.set_ylim(low, high + 0.09 * (high - low))     # room for the gap label
ax.text(2042.5, ax.get_ylim()[1], "not computed\n(2015–2070)", ha="center", va="top",
        fontsize=9.5, color="#8c9295")
ax.set_xlabel("Year", color=INK)
ax.set_ylabel("Global area-weighted mean SSS (psu)", color=INK)
ax.set_title("Global mean sea surface salinity: historical 1950–2014 and ssp245 2071–2100\n"
             "thin lines the annual series, dashed the least-squares fit of each segment; "
             "the two segments are never joined",
             fontsize=12, color=INK)
ax.legend(loc="center left", frameon=False, fontsize=10, labelcolor=INK)
ax.grid(axis="y", color="#eceff0", linewidth=0.9)
ax.set_axisbelow(True)
ax.tick_params(colors=INK)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
for side in ("left", "bottom"):
    ax.spines[side].set_color("#c8ccce")
fig.savefig(OUT / "fig3_global_mean_timeseries.png", dpi=150, facecolor="white")
plt.close(fig)


# ---------------------------------------------------------------------------
# Numeric answer
# ---------------------------------------------------------------------------
def three_sig(x):
    return float(f"{x:.3g}")


answer = {
    "case_id": "T5",
    "title": "Sea surface salinity: mean-state bias, historical trend and ssp245 scenario change",
    "models": MODELS,
    "definition_of_sss": "0-10 m volume-weighted mean salinity, on each source's own grid; "
                         "one level for CESM2 (5 m) and MPI-ESM1-2-LR (6 m), two for GFDL-ESM4 "
                         "(2.5 and 10 m)",
    "woa23_global_mean_sss_psu": round(woa_global, 2),
    "global_mean_sss_psu": {m: round(global_mean[m], 2) for m in MODELS},
    "mean_state_bias_psu": {m: round(mean_bias[m], 2) for m in MODELS},

    "comparison_mask": {"grid": "WOA23 1-degree", "n_cells": n_common,
                        "definition": common.attrs["description"]},
    "rmse_vs_woa23_psu": {m: round(pattern[m]["rmse_psu"], 2) for m in MODELS},
    "spatial_correlation_vs_woa23": {m: round(pattern[m]["spatial_correlation"], 3) for m in MODELS},

    "trends_psu_per_decade": {
        m: {label: {
                "slope": three_sig(trends[(m, label)]["slope_psu_per_decade"]),
                "ci95": [three_sig(trends[(m, label)]["ci95_low_psu_per_decade"]),
                         three_sig(trends[(m, label)]["ci95_high_psu_per_decade"])],
                "r2": round(trends[(m, label)]["r2"], 3),
                "p_value": three_sig(trends[(m, label)]["p_value"]),
                "ci_excludes_zero": trends[(m, label)]["ci_excludes_zero"],
            } for label in ("historical", "ssp245")}
        for m in MODELS},
    "trend_ratio_ssp245_over_historical": {
        m: {"ratio": round(trend_ratio[m]["ratio"], 1),
            "denominator_ci_excludes_zero": trend_ratio[m]["denominator_significant"],
            "interpretable_as_a_multiple": trend_ratio[m]["interpretable"]}
        for m in MODELS},

    "delta_sss_psu": {m: {k: round(v, 2) for k, v in delta_sss[m].items()
                          if not k.endswith("_definition")} for m in MODELS},
    "basin_definitions": {b: delta_sss[MODELS[0]][f"{b}_definition"] for b in BASINS},
    "marginal_seas": "folded into their parent ocean (select_basin default)",
    "southern_and_arctic_oceans": "separate basins in the WOA23 mask, not included in the Atlantic, "
                                  "Indian or Pacific numbers",
    "ensemble_members": {f"{m}|{label}": members[(m, label)]
                         for m in MODELS for label in ("historical", "ssp245")},
    "member_note": member_note,

    "full_precision": {
        "global_mean_sss_psu": global_mean,
        "woa23_global_mean_sss_psu": woa_global,
        "mean_state_bias_psu": mean_bias,
        "pattern": pattern,
        "trends": {f"{m}|{label}": {k: v for k, v in trends[(m, label)].items() if k != "fitted"}
                   for m in MODELS for label in ("historical", "ssp245")},
        "trend_ratio": {m: trend_ratio[m]["ratio"] for m in MODELS},
        "delta_sss": {m: {k: v for k, v in delta_sss[m].items() if not k.endswith("_definition")}
                      for m in MODELS},
    },

    "figures": ["fig1_sss_bias_maps.png", "fig2_sss_trend_maps.png",
                "fig3_global_mean_timeseries.png"],
    "method_notes": [
        "SSS is the 0-10 m volume-weighted mean everywhere, models and WOA23 alike.",
        "Trends are fitted on ANNUAL means with a serial-correlation correction, and converted "
        "to psu per decade; fit_1D_trend returns per year on this axis, fit_trend_map per decade.",
        "Every global and basin mean is taken on the model's native grid with its own areacello; "
        "the only regridding is onto the 1-degree grid for the pattern metrics of step 3 and for "
        "the display maps.",
        "Basin masks are carried onto the model grid by nearest cell, with marginal seas folded "
        "into their parent ocean. The Southern and Arctic Oceans are separate basins in the WOA23 "
        "mask and are not part of the three ocean numbers.",
        member_note,
        "The global mean of 0-10 m salinity is not a conserved quantity, unlike the full-depth "
        "mean, so a non-zero global surface trend is physical rather than a numerical artefact.",
    ],
    "tools_used": [
        "load_regular_data", "load_climatology_data", "load_landsea_mask", "compute_basic_stat",
        "temporal_aggregate", "regrid_horizontal", "build_comparison_mask",
        "compute_difference_metrics", "analyze_correlation", "fit_1D_trend", "fit_trend_map",
        "select_basin",
    ],
}

(OUT / "answer_T5.json").write_text(json.dumps(answer, indent=2) + "\n")

md = [
    "# Reference answer - T5: sea surface salinity bias, trend and scenario change",
    "",
    "SSS is the 0-10 m volume-weighted mean salinity throughout, on each source's own grid.",
    "Models are CESM2, GFDL-ESM4 and MPI-ESM1-2-LR; observations are WOA23 1971-2000.",
    "",
    "## Step 2. Mean state, 1971-2000",
    "",
    "| Source | global mean SSS (psu) | bias vs WOA23 (psu) |",
    "|---|---|---|",
    f"| WOA23 (observed) | {woa_global:.2f} | - |",
]
for model in MODELS:
    md.append(f"| {model} | {global_mean[model]:.2f} | {mean_bias[model]:+.2f} |")

md += [
    "",
    f"## Step 3. Pattern agreement on the shared 1-degree mask ({n_common} cells)",
    "",
    "| Model | RMSE (psu) | spatial correlation |",
    "|---|---|---|",
]
for model in MODELS:
    md.append(f"| {model} | {pattern[model]['rmse_psu']:.2f} | "
              f"{pattern[model]['spatial_correlation']:.3f} |")

md += [
    "",
    "## Step 4. Global mean SSS trend",
    "",
    "| Model | Period | trend (psu/decade) | 95% CI | R2 | CI excludes zero |",
    "|---|---|---|---|---|---|",
]
for model in MODELS:
    for label, span in (("historical", "1950-2014"), ("ssp245", "2071-2100")):
        e = trends[(model, label)]
        md.append(f"| {model} | {label} {span} | {e['slope_psu_per_decade']:.3g} | "
                  f"[{e['ci95_low_psu_per_decade']:.3g}, {e['ci95_high_psu_per_decade']:.3g}] | "
                  f"{e['r2']:.3f} | {'yes' if e['ci_excludes_zero'] else 'NO'} |")

md += [
    "",
    "| Model | ssp245 / historical slope ratio | interpretable? |",
    "|---|---|---|",
]
for model in MODELS:
    r = trend_ratio[model]
    md.append(f"| {model} | {r['ratio']:.1f} | "
              f"{'yes' if r['interpretable'] else 'NO - denominator not distinguishable from zero'} |")

md += [
    "",
    "**Every historical trend has a confidence interval that includes zero**, so none of the three",
    "is distinguishable from no trend at all. The ssp245 trends are an order of magnitude larger",
    "and two of the three exclude zero. The ratio is therefore reported but is NOT a meaningful",
    "multiple for any model: its denominator is a number consistent with zero, and the ratio's",
    "magnitude is set by how close to zero that denominator happens to fall.",
    "",
    "## Step 5. Scenario change, ssp245 2071-2100 minus historical 1995-2014",
    "",
    "Marginal seas are folded into their parent ocean. The Southern and Arctic Oceans are separate",
    "basins in the WOA23 mask and are not part of the Atlantic, Indian or Pacific numbers.",
    "",
    member_note,
    "",
    "| Model | Global | Atlantic | Indian | Pacific |",
    "|---|---|---|---|---|",
]
for model in MODELS:
    e = delta_sss[model]
    md.append(f"| {model} | {e['global']:+.2f} | {e['Atlantic']:+.2f} | "
              f"{e['Indian']:+.2f} | {e['Pacific']:+.2f} |")

md += [
    "",
    "## Figures",
    "",
    "- `fig1_sss_bias_maps.png` - model minus WOA23 SSS bias, three panels on one shared scale.",
    "- `fig2_sss_trend_maps.png` - historical 1950-2014 SSS trend in psu/decade, p < 0.05 stippled.",
    "- `fig3_global_mean_timeseries.png` - global mean SSS, the two experiments drawn as separate",
    "  segments with the 2015-2070 gap left empty, each with its own fitted line.",
]
(OUT / "answer_T5.md").write_text("\n".join(md) + "\n")

print("\nWrote answer_T5.json, answer_T5.md and three figures")
