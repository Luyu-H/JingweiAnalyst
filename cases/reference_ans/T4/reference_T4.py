"""
Reference answer for case T4 - Consistency of a mixed layer depth operator with
each model's own mlotst output.

Run from the repository root with the project's conda environment:

    source /home/lorraine/miniconda3/bin/activate base && conda activate jingwei_QA
    PYTHONPATH=. python cases/reference_ans/T4/reference_T4.py

Writes answer_T4.json, answer_T4.md, fig1_difference_maps.png and
fig2_scatter_panels.png next to this file. Takes roughly twenty minutes, most
of it the 180 monthly MLD computations.

Every physical quantity is produced by a tool from src/tools. Nothing here
re-implements an MLD, an area weighting or an error metric by hand, and nothing
is regridded except the difference fields drawn in figure 1.
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
from src.tools.data_loader import load_regular_data
from src.tools.preprocessing import regrid_horizontal
from src.tools.scientific_calculator import compute_mixed_layer_depth
from src.tools.basic_calculator import (compute_difference_metrics, analyze_correlation,
                                        compute_exceedance_fraction, compute_basic_stat)

OUT = Path(__file__).resolve().parent
MODELS = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
MONTHS = {3: "March", 9: "September"}
WINDOW = ["1971-01", "2000-12"]
YEARS = range(1971, 2001)
DENSITY_THRESHOLD = 0.03      # kg/m3, the CMIP6 protocol's criterion for mlotst
LARGE_THRESHOLD = 0.125       # kg/m3, the larger criterion some studies use; sensitivity check only
SENSITIVITY_MONTH = "1971-03"
DIFF_THRESHOLD = 20.0         # m, the disagreement the case counts area for

xr.set_options(keep_attrs=True)

climatology = {}   # (model, month) -> {"mld": field, "mlotst": field}
metrics = {}       # (model, month) -> dict of numbers
ocean_area = {}    # model -> total ocean area in m2
sensitivity = {}   # model -> threshold -> bias and RMSD for SENSITIVITY_MONTH


for model in MODELS:
    started = time.time()
    # Opened once and kept lazy; each month is pulled out of these below, which
    # avoids reopening the multi-file archives sixty times over.
    thetao = load_regular_data("cmip6_historical", "thetao", model=model,
                               time_range=WINDOW, chunk_time=1)
    salinity = load_regular_data("cmip6_historical", "so", model=model,
                                 time_range=WINDOW, chunk_time=1)
    published = load_regular_data("cmip6_historical", "mlotst", model=model,
                                  time_range=WINDOW, chunk_time=1)

    # The models publish mlotst as sigma-T referenced to the ocean surface, so
    # the operator is referenced to the shallowest level rather than to the 10 m
    # of the observational criterion.
    reference_depth = float(thetao["depth"].values[0])

    areas = load_regular_data("cmip6_historical", "areacello", model=model)
    ocean_area[model] = float(compute_basic_stat(areas, op="sum").compute())

    for month, month_name in MONTHS.items():
        operator_years, published_years = [], []
        for year in YEARS:
            stamp = f"{year}-{month:02d}"
            temp_m = thetao.sel(time=stamp)
            salt_m = salinity.sel(time=stamp)
            # One MLD per month, averaged afterwards: the MLD of an averaged
            # profile is not the average of the MLDs.
            mld_m = compute_mixed_layer_depth(
                temp_m, salt_m, ref_depth=reference_depth, method="density",
                density_threshold=DENSITY_THRESHOLD).squeeze(drop=True)
            operator_years.append(mld_m.compute())
            published_years.append(published.sel(time=stamp).squeeze(drop=True).compute())

        # skipna=False keeps a cell only where all 30 years are valid, so both
        # climatologies rest on the same set of columns.
        operator = xr.concat(operator_years, dim="year").mean("year", skipna=False)
        publication = xr.concat(published_years, dim="year").mean("year", skipna=False)
        operator.attrs = dict(operator_years[0].attrs)
        publication.attrs = dict(published_years[0].attrs)
        climatology[(model, month)] = {"mld": operator, "mlotst": publication}

        difference = compute_difference_metrics(operator, publication,
                                                dim=["lat", "lon"], weighting="area")
        correlation = analyze_correlation(operator, publication, dim=["lat", "lon"],
                                          method="pearson", weighting="area")
        exceedance = compute_exceedance_fraction(difference["difference"],
                                                 threshold=DIFF_THRESHOLD,
                                                 comparison="abs_greater",
                                                 dim=["lat", "lon"], weighting="area")
        shared_area = float(exceedance["total_amount"])
        metrics[(model, month)] = {
            "mean_difference_m": float(difference["bias"]),
            "rmsd_m": float(difference["rmse"]),
            "spatial_correlation": float(correlation["correlation"]),
            "area_percent_abs_diff_gt_20m": float(exceedance["fraction_percent"]),
            "shared_area_m2": shared_area,
            "shared_area_percent_of_ocean": 100.0 * shared_area / ocean_area[model],
            "n_shared_cells": int(difference["n_valid"]),
            "mean_operator_mld_m": float(difference["mean_data"]),
            "mean_mlotst_m": float(difference["mean_reference"]),
            "weighting": str(difference["weighting"]),
        }
        entry = metrics[(model, month)]
        print(f"[{model:14s} {month_name:9s}] bias {entry['mean_difference_m']:+7.2f} m   "
              f"RMSD {entry['rmsd_m']:7.2f} m   r {entry['spatial_correlation']:.4f}   "
              f"|diff|>{DIFF_THRESHOLD:.0f} m over {entry['area_percent_abs_diff_gt_20m']:5.2f}% "
              f"of {entry['shared_area_percent_of_ocean']:.1f}% shared ocean")

    # Threshold sensitivity, one month: the same operator at the larger
    # 0.125 kg/m3 some studies use, against the model's own mlotst for March
    # 1971. This is what backs the statement that 0.03 kg/m3 - the protocol's
    # value (Griffies et al. 2016, H24) - is the one mlotst is consistent with.
    temp_s = thetao.sel(time=SENSITIVITY_MONTH).compute()
    salt_s = salinity.sel(time=SENSITIVITY_MONTH).compute()
    published_s = published.sel(time=SENSITIVITY_MONTH).squeeze(drop=True).compute()
    sensitivity[model] = {}
    for threshold in (DENSITY_THRESHOLD, LARGE_THRESHOLD):
        mld_s = compute_mixed_layer_depth(temp_s, salt_s, ref_depth=reference_depth, method="density",
                                          density_threshold=threshold).squeeze(drop=True)
        check = compute_difference_metrics(mld_s, published_s, dim=["lat", "lon"], weighting="area")
        sensitivity[model][threshold] = {"bias_m": float(check["bias"]), "rmsd_m": float(check["rmse"])}
    print(f"  {model} {SENSITIVITY_MONTH} RMSD at {DENSITY_THRESHOLD} kg/m3 "
          f"{sensitivity[model][DENSITY_THRESHOLD]['rmsd_m']:.1f} m, at {LARGE_THRESHOLD} kg/m3 "
          f"{sensitivity[model][LARGE_THRESHOLD]['rmsd_m']:.1f} m")

    print(f"  {model} done in {time.time() - started:.0f} s "
          f"(reference depth {reference_depth:.2f} m)")


# ---------------------------------------------------------------------------
# Step 6. Which model is most internally consistent
# ---------------------------------------------------------------------------
# Ranked on the RMSD averaged over the two months, which penalises disagreement
# of either sign rather than letting a warm and a cold region cancel.
mean_rmsd = {m: np.mean([metrics[(m, month)]["rmsd_m"] for month in MONTHS]) for m in MODELS}
negative_biases = sum(entry["mean_difference_m"] < 0 for entry in metrics.values())
ranking = sorted(MODELS, key=lambda m: mean_rmsd[m])
print("\nconsistency ranking (mean RMSD over the two months, best first): "
      + ", ".join(f"{m} {mean_rmsd[m]:.1f} m" for m in ranking))


# ---------------------------------------------------------------------------
# Figure 1. Difference maps, regridded to 1 degree FOR DISPLAY ONLY
# ---------------------------------------------------------------------------
INK = "#2f3437"
LAND = "#d9d9d6"

display = {}
for model in MODELS:
    for month in MONTHS:
        field = climatology[(model, month)]["mld"] - climatology[(model, month)]["mlotst"]
        field.attrs = dict(climatology[(model, month)]["mld"].attrs)
        display[(model, month)] = regrid_horizontal(field, target="woa23", method="conservative")
    print(f"[fig1] {model} difference fields regridded for display")

stack = np.concatenate([np.abs(display[k].values[np.isfinite(display[k].values)])
                        for k in display])
vmax = float(np.ceil(np.percentile(stack, 98) / 5) * 5)
print(f"[fig1] shared diverging scale +/-{vmax:.0f} m "
      f"(98th percentile of |difference| = {np.percentile(stack, 98):.1f} m)")

lon = display[(MODELS[0], 3)]["lon"].values
lat = display[(MODELS[0], 3)]["lat"].values

fig, axes = plt.subplots(3, 2, figsize=(14.0, 11.2))
fig.subplots_adjust(left=0.06, right=0.985, bottom=0.115, top=0.885, wspace=0.11, hspace=0.30)
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

for row, model in enumerate(MODELS):
    for column, month in enumerate(MONTHS):
        ax = axes[row, column]
        values = display[(model, month)].values
        ax.set_facecolor(LAND)
        mesh = ax.pcolormesh(lon, lat, values, cmap="RdBu_r", norm=norm, shading="auto")
        entry = metrics[(model, month)]
        ax.set_title(f"{model} — {MONTHS[month]}\n"
                     f"bias {entry['mean_difference_m']:+.1f} m,   "
                     f"RMSD {entry['rmsd_m']:.1f} m,   "
                     f"r = {entry['spatial_correlation']:.3f}",
                     fontsize=10.5, color=INK)
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        ax.set_xticks([-180, -90, 0, 90, 180])
        ax.set_yticks([-90, -45, 0, 45, 90])
        ax.tick_params(colors=INK, labelsize=8.5)
        if row == 2:
            ax.set_xlabel("Longitude (°E)", color=INK, fontsize=9)
        if column == 0:
            ax.set_ylabel("Latitude (°N)", color=INK, fontsize=9)
        for spine in ax.spines.values():
            spine.set_color("#c8ccce")

cax = fig.add_axes([0.25, 0.055, 0.50, 0.016])
bar = fig.colorbar(mesh, cax=cax, orientation="horizontal", extend="both")
bar.set_label(f"Operator MLD − published mlotst (m), one shared scale for all six panels",
              color=INK, fontsize=10)
bar.ax.tick_params(colors=INK, labelsize=9)

fig.suptitle("Mixed layer depth operator against each model's own mlotst, 1971–2000 climatology\n"
             f"potential density threshold {DENSITY_THRESHOLD} kg/m³ referenced to the shallowest "
             "level; every number computed on the native grid,\nthe maps regridded to 1° for "
             "display only",
             fontsize=12.5, color=INK, y=0.975)

fig.savefig(OUT / "fig1_difference_maps.png", dpi=150, facecolor="white")
plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2. Scatter of the two fields, on the native grid
# ---------------------------------------------------------------------------
lower, upper = 2.0, 4000.0
ticks = [10, 100, 1000]

fig, axes = plt.subplots(3, 2, figsize=(10.4, 14.0))
fig.subplots_adjust(left=0.085, right=0.975, bottom=0.06, top=0.905, wspace=0.16, hspace=0.26)

for row, model in enumerate(MODELS):
    for column, month in enumerate(MONTHS):
        ax = axes[row, column]
        operator = climatology[(model, month)]["mld"].values.ravel()
        publication = climatology[(model, month)]["mlotst"].values.ravel()
        both = np.isfinite(operator) & np.isfinite(publication)

        ax.plot([lower, upper], [lower, upper], color="#9aa0a3", linewidth=1.2,
                linestyle="--", zorder=3)
        ax.scatter(operator[both], publication[both], s=3.0, color="#1f4e79",
                   alpha=0.16, linewidths=0, zorder=2)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lower, upper)
        ax.set_ylim(lower, upper)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels([str(t) for t in ticks])
        ax.set_yticklabels([str(t) for t in ticks])
        ax.set_aspect("equal", adjustable="box")

        entry = metrics[(model, month)]
        ax.set_title(f"{model} — {MONTHS[month]}", fontsize=10.5, color=INK)
        ax.text(0.04, 0.96,
                f"bias {entry['mean_difference_m']:+.1f} m\n"
                f"RMSD {entry['rmsd_m']:.1f} m\n"
                f"r = {entry['spatial_correlation']:.3f}",
                transform=ax.transAxes, fontsize=9, color=INK, va="top", ha="left",
                bbox=dict(facecolor="white", edgecolor="none", boxstyle="square,pad=0.25"))
        if row == 2:
            ax.set_xlabel("Operator MLD (m)", color=INK, fontsize=9.5)
        if column == 0:
            ax.set_ylabel("Published mlotst (m)", color=INK, fontsize=9.5)
        ax.grid(color="#eceff0", linewidth=0.8, which="major")
        ax.set_axisbelow(True)
        ax.tick_params(colors=INK, labelsize=9)
        for spine in ax.spines.values():
            spine.set_color("#c8ccce")

fig.suptitle("Operator MLD against published mlotst, cell by cell on each model's native grid\n"
             "1971–2000 climatology; logarithmic axes, identical range in every panel, "
             "dashed line is y = x",
             fontsize=12.5, color=INK, y=0.972)

fig.savefig(OUT / "fig2_scatter_panels.png", dpi=150, facecolor="white")
plt.close(fig)


# ---------------------------------------------------------------------------
# Numeric answer
# ---------------------------------------------------------------------------
def rounded(model, month):
    entry = metrics[(model, month)]
    return {
        "mean_difference_m": round(entry["mean_difference_m"], 1),
        "rmsd_m": round(entry["rmsd_m"], 1),
        "spatial_correlation": round(entry["spatial_correlation"], 3),
        "area_percent_abs_diff_gt_20m": round(entry["area_percent_abs_diff_gt_20m"], 1),
        "shared_area_percent_of_ocean": round(entry["shared_area_percent_of_ocean"], 1),
        "n_shared_cells": entry["n_shared_cells"],
    }


answer = {
    "case_id": "T4",
    "title": "Consistency of a mixed layer depth operator with each model's own mlotst output",
    "window": "1971-01 to 2000-12, the March and September climatologies",
    "models": MODELS,
    "method": {
        "criterion": f"potential density threshold {DENSITY_THRESHOLD} kg/m3, "
                     "referenced to each model's shallowest level, crossing depth interpolated "
                     "between the bracketing levels",
        "order_of_operations": "MLD computed for each of the 30 months, then averaged; the "
                               "profiles are never averaged first",
        "climatology_rule": "a cell enters either climatology only if all 30 years are valid there",
        "grid": "every number on the model's native grid, weighted by its own areacello; "
                "no regridding except the display maps in figure 1",
    },

    "results": {model: {MONTHS[month]: rounded(model, month) for month in MONTHS}
                for model in MODELS},

    "mean_rmsd_over_both_months_m": {m: round(float(mean_rmsd[m]), 1) for m in MODELS},
    "consistency_ranking": ranking,
    "consistency_ranking_criterion": "mean RMSD over March and September, smallest first",

    "full_precision": {f"{model}|{MONTHS[month]}": metrics[(model, month)]
                       for model in MODELS for month in MONTHS},
    "total_ocean_area_m2": ocean_area,

    "threshold_sensitivity": {
        "month": SENSITIVITY_MONTH,
        "note": "operator against the same month's mlotst at the protocol's 0.03 kg/m3 and at the "
                "larger 0.125 kg/m3; not part of the graded answer",
        "results": {m: {f"{t} kg/m3": {k: round(v, 1) for k, v in sensitivity[m][t].items()}
                        for t in sensitivity[m]} for m in MODELS},
    },

    "figures": ["fig1_difference_maps.png", "fig2_scatter_panels.png"],
    "method_notes": [
        "The MLD is computed month by month and averaged afterwards. One MLD taken from the "
        "30-year mean profiles would be systematically shallower and is the main error to avoid.",
        "The reference depth is each model's shallowest level, because all three publish mlotst "
        "as sigma-T referenced to the ocean surface, not to 10 m.",
        f"The density threshold is {DENSITY_THRESHOLD} kg/m3, the value the CMIP6 protocol defines "
        "mlotst with: a buoyancy difference of 0.0003 m s-2, which Griffies et al. (2016, H24) equate "
        "to 0.03 kg/m3 and attribute to Levitus (1982). At the larger 0.125 kg/m3 some studies use, the "
        f"{SENSITIVITY_MONTH} RMSD against mlotst rises from "
        + ", ".join(f"{m} {sensitivity[m][DENSITY_THRESHOLD]['rmsd_m']:.0f} to "
                    f"{sensitivity[m][LARGE_THRESHOLD]['rmsd_m']:.0f} m" for m in MODELS) + ".",
        "Area percentages are fractions of ocean AREA, weighted by each model's areacello, not "
        "fractions of grid cells; on a curvilinear grid the two differ substantially.",
        "Cells where the threshold is never reached leave the operator MLD undefined; they are "
        "excluded from both climatologies and from the denominator, and the shared-area "
        "percentage records how much was dropped. Excluding them is necessary, not only convenient: "
        "the models' own mlotst encodes such columns differently (GFDL-ESM4 and MPI-ESM1-2-LR "
        "report roughly the full column depth, CESM2 mostly 0 m), so no single fill value is fair "
        "to all three.",
    ],
    "why_they_differ": [
        "Both fields come from the same simulation, so the residual is a property of the two "
        "diagnostics. The protocol computes mlotst online at every time step, from a surface parcel "
        "moved to each level's own pressure, and archives the monthly mean of those instantaneous "
        "depths (Griffies et al. 2016, H24). The operator takes a sigma0 difference on the "
        "monthly-mean profile. A monthly mean is itself an averaged profile, so the smoothing that "
        "makes a 30-year-mean profile too shallow also makes the monthly one shallower than the "
        "mean of the instantaneous depths: the bias is expected to be negative, and it is negative in "
        f"{negative_biases} of the {len(metrics)} model-and-month combinations.",
        "GFDL-ESM4 computes mlotst on its native 75 hybrid layers with a 2 m top layer, but its "
        "archived thetao and so were remapped to 35 z levels and regridded to 1 degree, so the "
        "operator sees a coarser profile than the model's own diagnostic did. Part of its larger "
        "RMSD is therefore the archive's resolution, and the ranking measures consistency between "
        "the archived fields, not the quality of the model's mixed layer.",
    ],
    "tools_used": [
        "load_regular_data", "compute_mixed_layer_depth", "compute_basic_stat",
        "compute_difference_metrics", "analyze_correlation", "compute_exceedance_fraction",
        "regrid_horizontal",
    ],
}

(OUT / "answer_T4.json").write_text(json.dumps(answer, indent=2) + "\n")

md = [
    "# Reference answer - T4: MLD operator against each model's own mlotst",
    "",
    "1971-2000 March and September climatologies. Potential density threshold "
    f"{DENSITY_THRESHOLD} kg/m3",
    "referenced to each model's shallowest level. MLD computed month by month and then averaged.",
    "Every number on the model's native grid, area-weighted by its own areacello.",
    "",
    "## Steps 3-5. Bias, error, correlation and area of disagreement",
    "",
    "| Model | Month | bias (m) | RMSD (m) | correlation | area with abs diff > 20 m | shared area |",
    "|---|---|---|---|---|---|---|",
]
for model in MODELS:
    for month in MONTHS:
        e = rounded(model, month)
        md.append(f"| {model} | {MONTHS[month]} | {e['mean_difference_m']:+.1f} | "
                  f"{e['rmsd_m']:.1f} | {e['spatial_correlation']:.3f} | "
                  f"{e['area_percent_abs_diff_gt_20m']:.1f}% | "
                  f"{e['shared_area_percent_of_ocean']:.1f}% of ocean |")

md += [
    "",
    "The bias is the area-weighted mean of (operator MLD - mlotst); a negative value means the",
    "operator puts the mixed layer shallower than the model's own diagnostic. The shared-area",
    "column is the ocean area where both fields are defined in all 30 years, as a percentage of",
    "the model's total ocean area.",
    "",
    "## Step 6. Consistency ranking",
    "",
    "| Model | mean RMSD over both months (m) |",
    "|---|---|",
]
for model in ranking:
    md.append(f"| {model} | {mean_rmsd[model]:.1f} |")
md += [
    "",
    f"Most internally consistent: **{ranking[0]}**. Least: **{ranking[-1]}**.",
    "",
    "## Why the two are not expected to agree exactly",
    "",
]
md += [f"- {line}" for line in answer["why_they_differ"]]
md += [
    f"- Threshold check, {SENSITIVITY_MONTH}: RMSD against mlotst at {DENSITY_THRESHOLD} and at "
    f"{LARGE_THRESHOLD} kg/m3 is "
    + "; ".join(f"{m} {sensitivity[m][DENSITY_THRESHOLD]['rmsd_m']:.1f} and "
                f"{sensitivity[m][LARGE_THRESHOLD]['rmsd_m']:.1f} m" for m in MODELS)
    + ". The protocol's 0.03 kg/m3 is the one mlotst is consistent with.",
    "",
    "## Figures",
    "",
    "- `fig1_difference_maps.png` - 3x2 maps of the difference field, one shared diverging scale",
    f"  (+/-{vmax:.0f} m), regridded to 1 degree for display only.",
    "- `fig2_scatter_panels.png` - 3x2 scatter of operator MLD against mlotst on the native grid,",
    "  logarithmic axes with identical ranges, dashed y = x.",
]
(OUT / "answer_T4.md").write_text("\n".join(md) + "\n")

print("\nWrote answer_T4.json, answer_T4.md, fig1_difference_maps.png, fig2_scatter_panels.png")
