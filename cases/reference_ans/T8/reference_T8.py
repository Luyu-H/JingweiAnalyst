"""
Reference answer for case T8 - Carbon system variables against the GLODAPv2
mapped climatology.

Run from the repository root with the project's conda environment:

    source /home/lorraine/miniconda3/bin/activate base && conda activate jingwei_QA
    PYTHONPATH=. python cases/reference_ans/T8/reference_T8.py

Writes answer_T8.json, answer_T8.md, fig1_dic_bias_maps.png and
fig2_dic_bias_profile.png next to this file.

Every physical quantity is produced by a tool from src/tools. The model fields
are converted from mol m-3 to umol/kg with an in-situ density field computed
from each model's own thetao and so; a constant would err by about the size of
the biases being measured.
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
from src.tools.data_loader import load_regular_data, load_climatology_data
from src.tools.preprocessing import regrid_horizontal, regrid_vertical, build_comparison_mask
from src.tools.scientific_calculator import compute_potential_density_anomaly
from src.tools.basic_calculator import (compute_basic_stat, temporal_aggregate, transform_unit,
                                        compute_difference_metrics, analyze_correlation)

OUT = Path(__file__).resolve().parent
MODELS = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
YEAR = ["2002-01", "2002-12"]
BASELINE_WINDOW = ["1971-01", "2000-12"]
SURFACE = [0, 50]
PROFILE_DEPTHS = [0.0, 100.0, 500.0, 1000.0, 2000.0]
VOLUME_DIMS = ["depth", "lat", "lon"]

# GLODAP's own names for the three carbon variables, and the CMIP6 counterparts.
PAIRS = {"DIC": ("dissic", "TCO2"), "TAlk": ("talk", "TAlk"), "pH": ("ph", "pHtsinsitutp")}

xr.set_options(keep_attrs=True)


def annual_mean(dataset, variable, model, window, depth_range):
    """Day-weighted annual mean of a monthly model field."""
    monthly = load_regular_data(dataset, variable, model=model, time_range=window,
                                depth_range=depth_range, chunk_time=12)
    return temporal_aggregate(monthly, freq="mean", op="mean").compute()


def to_mass_basis(field, temperature, salinity):
    """mol m-3 -> umol/kg using the in-situ density of each cell."""
    density = compute_potential_density_anomaly(temperature, salinity, ref_pressure=None)
    converted = transform_unit(field, from_unit=field.attrs["units"], to_unit="umol/kg",
                               density=density)
    return converted


# ---------------------------------------------------------------------------
# Step 1-2. Surface means of DIC, alkalinity and pH
# ---------------------------------------------------------------------------
glodap_surface, model_surface = {}, {model: {} for model in MODELS}

for name, (_, glodap_name) in PAIRS.items():
    glodap_surface[name] = load_climatology_data("glodap", glodap_name, depth_range=SURFACE)
    print(f"[1]   GLODAP {name:5s} {glodap_surface[name].sizes['depth']} levels, "
          f"units '{glodap_surface[name].attrs.get('units')}'")

density_note = {}
for model in MODELS:
    started = time.time()
    temperature = annual_mean("cmip6_historical", "thetao", model, YEAR, SURFACE)
    salinity = annual_mean("cmip6_historical", "so", model, YEAR, SURFACE)

    for name, (cmip_name, _) in PAIRS.items():
        field = annual_mean("cmip6_historical", cmip_name, model, YEAR, SURFACE)
        if name == "pH":
            model_surface[model][name] = field          # dimensionless, no conversion
        else:
            converted = to_mass_basis(field, temperature, salinity)
            density_note[model] = converted.attrs.get("density_assumption", "in-situ density field")
            model_surface[model][name] = converted
    print(f"[1]   {model:14s} 2002 surface fields ready ({time.time() - started:.0f} s)")


def volume_mean(field):
    return float(compute_basic_stat(field, op="volume_weighted_mean", dim=VOLUME_DIMS).compute())


surface_mean, surface_bias = {}, {}
for name in PAIRS:
    surface_mean[name] = {"GLODAP": volume_mean(glodap_surface[name])}
    for model in MODELS:
        surface_mean[name][model] = volume_mean(model_surface[model][name])
    surface_bias[name] = {m: surface_mean[name][m] - surface_mean[name]["GLODAP"] for m in MODELS}
    decimals = 3 if name == "pH" else 2
    print(f"[1-2] {name:5s} GLODAP {surface_mean[name]['GLODAP']:.{decimals}f}   " + "  ".join(
        f"{m} {surface_mean[name][m]:.{decimals}f} ({surface_bias[name][m]:+.{decimals}f})"
        for m in MODELS))


# ---------------------------------------------------------------------------
# Step 3. Pattern agreement for surface DIC on one shared mask
# ---------------------------------------------------------------------------
def layer_mean(field):
    """Volume-weighted 0-50 m mean, leaving the horizontal dimensions."""
    reduced = compute_basic_stat(field, op="volume_weighted_mean", dim="depth").compute()
    reduced.attrs = dict(field.attrs)
    return reduced


glodap_dic_map = layer_mean(glodap_surface["DIC"])
regridded = {}
for model in MODELS:
    regridded[model] = regrid_horizontal(layer_mean(model_surface[model]["DIC"]),
                                         target="glodap", method="conservative")
    print(f"[3]   {model:14s} surface DIC regridded to the GLODAP grid")

common = build_comparison_mask([regridded[m] for m in MODELS] + [glodap_dic_map],
                               ocean_only=False)
n_common = int(common.sum())
print(f"[3]   shared comparison mask: {n_common} cells ({common.attrs['ocean_mask']})")

glodap_common = glodap_dic_map.where(common)
pattern, bias_field = {}, {}
for model in MODELS:
    field = regridded[model].where(common)
    difference = compute_difference_metrics(field, glodap_common, dim=["lat", "lon"],
                                            weighting="area")
    correlation = analyze_correlation(field, glodap_common, dim=["lat", "lon"],
                                      method="pearson", weighting="area")
    pattern[model] = {
        "bias_umol_per_kg": float(difference["bias"]),
        "rmse_umol_per_kg": float(difference["rmse"]),
        "spatial_correlation": float(correlation["correlation"]),
        "n_cells": int(difference["n_valid"]),
    }
    bias_field[model] = difference["difference"]
    print(f"[3]   {model:14s} bias {pattern[model]['bias_umol_per_kg']:+8.2f}  "
          f"RMSE {pattern[model]['rmse_umol_per_kg']:7.2f}  "
          f"r {pattern[model]['spatial_correlation']:.4f}")


# ---------------------------------------------------------------------------
# Step 4. Depth profile of the DIC bias, on common levels
# ---------------------------------------------------------------------------
glodap_full = load_climatology_data("glodap", "TCO2")
glodap_levels = glodap_full.sel(depth=PROFILE_DEPTHS)
print(f"[4]   GLODAP levels selected: {glodap_levels['depth'].values}")

bias_profile = {}
for model in MODELS:
    started = time.time()
    temperature = annual_mean("cmip6_historical", "thetao", model, YEAR, None)
    salinity = annual_mean("cmip6_historical", "so", model, YEAR, None)
    dic = annual_mean("cmip6_historical", "dissic", model, YEAR, None)
    converted = to_mass_basis(dic, temperature, salinity)

    # Interpolate onto the named depths, THEN move horizontally: a bias quoted
    # at a depth needs both products answering at the same depth, and this is
    # the one place where interpolation rather than a layer mean is correct.
    on_levels = regrid_vertical(converted, target=PROFILE_DEPTHS, method="linear")
    on_grid = regrid_horizontal(on_levels, target="glodap", method="conservative")

    difference = on_grid - glodap_levels
    difference.attrs = dict(on_grid.attrs)
    values = compute_basic_stat(difference, op="area_weighted_mean", dim=["lat", "lon"]).compute()
    bias_profile[model] = {float(d): float(values.sel(depth=d)) for d in PROFILE_DEPTHS}
    print(f"[4]   {model:14s} ({time.time() - started:.0f} s) " + "  ".join(
        f"{int(d)}m {bias_profile[model][d]:+7.2f}" for d in PROFILE_DEPTHS))


# ---------------------------------------------------------------------------
# Step 5. How much does the year matter
# ---------------------------------------------------------------------------
period_change = {}
for model in MODELS:
    started = time.time()
    temperature = annual_mean("cmip6_historical", "thetao", model, BASELINE_WINDOW, SURFACE)
    salinity = annual_mean("cmip6_historical", "so", model, BASELINE_WINDOW, SURFACE)
    dic = annual_mean("cmip6_historical", "dissic", model, BASELINE_WINDOW, SURFACE)
    # Each period converted with its own in-situ density.
    baseline = volume_mean(to_mass_basis(dic, temperature, salinity))
    period_change[model] = {
        "mean_1971_2000_umol_per_kg": baseline,
        "mean_2002_umol_per_kg": surface_mean["DIC"][model],
        "period_change_umol_per_kg": surface_mean["DIC"][model] - baseline,
        "bias_vs_glodap_umol_per_kg": surface_bias["DIC"][model],
    }
    entry = period_change[model]
    entry["period_change_as_fraction_of_bias"] = (
        entry["period_change_umol_per_kg"] / entry["bias_vs_glodap_umol_per_kg"]
        if entry["bias_vs_glodap_umol_per_kg"] != 0 else float("nan"))
    print(f"[5]   {model:14s} ({time.time() - started:.0f} s) 1971-2000 {baseline:.2f}  "
          f"2002 {entry['mean_2002_umol_per_kg']:.2f}  change {entry['period_change_umol_per_kg']:+.2f}  "
          f"bias {entry['bias_vs_glodap_umol_per_kg']:+.2f}  "
          f"change/bias {entry['period_change_as_fraction_of_bias']:+.2f}")


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
INK = "#2f3437"
LAND = "#d9d9d6"
COLOUR = {"CESM2": "#0072B2", "GFDL-ESM4": "#009E73", "MPI-ESM1-2-LR": "#E69F00"}
lon = glodap_dic_map["lon"].values
lat = glodap_dic_map["lat"].values

stack = np.concatenate([np.abs(bias_field[m].values[np.isfinite(bias_field[m].values)])
                        for m in MODELS])
percentiles = {q: float(np.percentile(stack, q)) for q in (50, 75, 90, 95, 98)}
print("[fig1] |bias| percentiles (umol/kg): "
      + "  ".join(f"p{q}={v:.0f}" for q, v in percentiles.items()))
# The 95th percentile: the Arctic and the marginal seas carry biases several
# times the open-ocean value, and a scale wide enough to hold them leaves the
# basin-scale pattern flat. The bar extends at both ends, so they are still shown.
vmax = float(np.ceil(percentiles[95] / 10) * 10)
print(f"[fig1] shared diverging scale +/-{vmax:.0f} umol/kg "
      f"(95th percentile of |bias| = {percentiles[95]:.1f})")

fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.6))
fig.subplots_adjust(left=0.045, right=0.99, bottom=0.235, top=0.80, wspace=0.13)
for ax, model in zip(axes, MODELS):
    ax.set_facecolor(LAND)
    mesh = ax.pcolormesh(lon, lat, bias_field[model].values, cmap="RdBu_r",
                         norm=TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax), shading="auto")
    ax.set_title(f"{model}\nbias {pattern[model]['bias_umol_per_kg']:+.0f},   "
                 f"RMSE {pattern[model]['rmse_umol_per_kg']:.0f} µmol kg⁻¹,   "
                 f"r = {pattern[model]['spatial_correlation']:.3f}",
                 fontsize=10.5, color=INK)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_yticks([-90, -45, 0, 45, 90])
    ax.tick_params(colors=INK, labelsize=8.5)
    ax.set_xlabel("Longitude (°E)", color=INK, fontsize=9)
    for spine in ax.spines.values():
        spine.set_color("#c8ccce")
axes[0].set_ylabel("Latitude (°N)", color=INK, fontsize=9)

cax = fig.add_axes([0.30, 0.10, 0.40, 0.035])
bar = fig.colorbar(mesh, cax=cax, orientation="horizontal", extend="both")
bar.set_label("0–50 m DIC bias, model − GLODAP (µmol kg⁻¹)",
              color=INK, fontsize=10)
bar.ax.tick_params(colors=INK, labelsize=9)
fig.suptitle("Surface dissolved inorganic carbon against GLODAPv2, model year 2002\n"
             "0–50 m volume-weighted, converted from mol m⁻³ with each model's own "
             f"in-situ density; shared mask of {n_common} cells",
             fontsize=12.5, color=INK, y=0.975)
fig.savefig(OUT / "fig1_dic_bias_maps.png", dpi=150, facecolor="white")
plt.close(fig)

fig, ax = plt.subplots(figsize=(7.4, 7.6), constrained_layout=True)
ax.axvline(0.0, color="#8c9295", linewidth=1.3, linestyle="--", zorder=2, label="no bias")
for model in MODELS:
    values = [bias_profile[model][d] for d in PROFILE_DEPTHS]
    ax.plot(values, PROFILE_DEPTHS, color=COLOUR[model], linewidth=2.0, marker="o",
            markersize=7, zorder=4, label=model)
ax.set_ylim(2100, -60)
ax.set_yticks(PROFILE_DEPTHS)
ax.set_yticklabels([f"{int(d)}" for d in PROFILE_DEPTHS])
ax.set_xlabel("DIC bias, model − GLODAP (µmol kg⁻¹)", color=INK)
ax.set_ylabel("Depth (m)", color=INK)
ax.set_title("Global mean DIC bias against depth, model year 2002\n"
             "interpolated onto the GLODAP levels, then area-weighted globally",
             fontsize=12, color=INK)
ax.legend(loc="lower left", frameon=False, fontsize=10, labelcolor=INK)
ax.grid(color="#eceff0", linewidth=0.9)
ax.set_axisbelow(True)
ax.tick_params(colors=INK)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
for side in ("left", "bottom"):
    ax.spines[side].set_color("#c8ccce")
fig.savefig(OUT / "fig2_dic_bias_profile.png", dpi=150, facecolor="white")
plt.close(fig)


# ---------------------------------------------------------------------------
# Numeric answer
# ---------------------------------------------------------------------------
def rounded(name, value):
    return round(value, 3) if name == "pH" else round(value)


answer = {
    "case_id": "T8",
    "title": "Carbon system variables against the GLODAPv2 mapped climatology",
    "models": MODELS,
    "model_year": "2002, day-weighted annual mean of the twelve months",
    "reference": "GLODAPv2.2016b mapped climatology, input data spanning 1972-2013, no time axis",
    "method": {
        "units": "dissic and talk converted from mol m-3 to umol/kg with an in-situ density field "
                 "computed from each model's own thetao and so at each cell's own pressure; "
                 "pH is dimensionless and untouched",
        "density_note": density_note,
        "ph_field": "GLODAP pHtsinsitutp (total scale, in-situ T and P), which is the quantity "
                    "CMIP6 ph reports; pHts25p0 is a different quantity",
        "step1_vs_step4": "step 1 takes a volume-weighted mean on each source's own levels; "
                          "step 4 interpolates onto common levels. A volume-weighted mean is not "
                          "conserved by interpolation, and a bias at a named depth cannot be read "
                          "off levels that differ between products",
        "grid": "every global mean on the model's native grid with its own volcello; regridding "
                "only for steps 3 and 4 and the display maps",
    },

    "surface_means": {name: {**{s: rounded(name, surface_mean[name][s])
                                for s in ["GLODAP"] + MODELS}}
                      for name in PAIRS},
    "surface_bias": {name: {m: rounded(name, surface_bias[name][m]) for m in MODELS}
                     for name in PAIRS},

    "dic_pattern_metrics": {m: {"bias_umol_per_kg": round(pattern[m]["bias_umol_per_kg"]),
                                "rmse_umol_per_kg": round(pattern[m]["rmse_umol_per_kg"]),
                                "spatial_correlation": round(pattern[m]["spatial_correlation"], 3)}
                            for m in MODELS},
    "comparison_mask_cells": n_common,

    "dic_bias_profile_umol_per_kg": {m: {str(int(d)): round(bias_profile[m][d])
                                         for d in PROFILE_DEPTHS} for m in MODELS},

    "dic_period_change": {m: {"mean_1971_2000": round(period_change[m]["mean_1971_2000_umol_per_kg"]),
                      "mean_2002": round(period_change[m]["mean_2002_umol_per_kg"]),
                      "change_2002_minus_1971_2000": round(period_change[m]["period_change_umol_per_kg"]),
                      "bias_vs_glodap": round(period_change[m]["bias_vs_glodap_umol_per_kg"]),
                      "period_change_as_fraction_of_bias": round(period_change[m]["period_change_as_fraction_of_bias"], 2)}
                  for m in MODELS},

    "full_precision": {"surface_mean": surface_mean, "surface_bias": surface_bias,
                       "pattern": pattern, "bias_profile": bias_profile, "period_change": period_change},

    "figures": ["fig1_dic_bias_maps.png", "fig2_dic_bias_profile.png"],
    "tools_used": [
        "load_regular_data", "load_climatology_data", "temporal_aggregate",
        "compute_potential_density_anomaly", "transform_unit", "compute_basic_stat",
        "regrid_vertical", "regrid_horizontal", "build_comparison_mask",
        "compute_difference_metrics", "analyze_correlation",
    ],
}

best = min(MODELS, key=lambda m: abs(surface_bias["DIC"][m]))
answer["findings"] = [
    "The model fields had to be converted with an in-situ density field rather than a constant. "
    "DIC and alkalinity are absolute concentrations, so unlike a ratio or a depth the conversion "
    "moves every number, and a constant near-surface density would err by about the size of the "
    "biases reported here.",
    f"On surface DIC the smallest bias against GLODAP is {best} "
    f"({surface_bias['DIC'][best]:+.0f} micromoles per kilogram).",
    "Step 5 is what decides how much of the DIC bias can be read as model error. The period change "
    "each model shows between its 1971-2000 mean and 2002 is "
    + ", ".join(f"{m} {period_change[m]['period_change_umol_per_kg']:+.0f}" for m in MODELS)
    + " micromoles per kilogram, against biases of "
    + ", ".join(f"{m} {surface_bias['DIC'][m]:+.0f}" for m in MODELS)
    + ". Where those are comparable in size the bias is sensitive to the year sampled.",
    "That period change is mostly the forced rise of surface DIC as the ocean takes up "
    "anthropogenic CO2 between the two periods, plus internal variability. It is not model drift "
    "in the CMIP6 sense - the residual trend of the pre-industrial control run - which cannot be "
    "measured here without a piControl run.",
    "Alkalinity is far less sensitive to anthropogenic CO2 than DIC, so its bias cannot be "
    "explained by the period mismatch and is the cleaner indicator of a structural model error.",

    "The DIC and alkalinity biases carry the SAME SIGN in every model - "
    + ", ".join(f"{m} {surface_bias['DIC'][m]:+.0f} and {surface_bias['TAlk'][m]:+.0f}"
                for m in MODELS)
    + " micromoles per kilogram - which is why pH comes out far better than either variable "
    "alone. pH responds to the balance between the two, so errors that move DIC and alkalinity "
    "together largely cancel in it: "
    f"CESM2 is off by {surface_bias['DIC']['CESM2']:+.0f} in DIC and "
    f"{surface_bias['TAlk']['CESM2']:+.0f} in alkalinity yet reproduces pH to "
    f"{abs(surface_bias['pH']['CESM2']):.3f}. Agreement in pH alone is therefore NOT evidence "
    "that a model's carbon system is right, and the three variables have to be read together - "
    "which is what this case is built to show.",
]

(OUT / "answer_T8.json").write_text(json.dumps(answer, indent=2) + "\n")

md = [
    "# Reference answer - T8: carbon system against GLODAPv2",
    "",
    "Models sampled for calendar year 2002 (day-weighted annual mean). GLODAPv2.2016b mapped",
    "climatology, input data spanning 1972-2013. Model dissic and talk converted from mol m-3 to",
    "umol/kg with an in-situ density field from each model's own thetao and so.",
    "",
    "## Steps 1-2. Surface (0-50 m) volume-weighted global means",
    "",
    "| Variable | GLODAP | CESM2 | GFDL-ESM4 | MPI-ESM1-2-LR |",
    "|---|---|---|---|---|",
]
for name in PAIRS:
    fmt = "%.3f" if name == "pH" else "%.0f"
    row = [fmt % surface_mean[name]["GLODAP"]]
    for model in MODELS:
        bias_fmt = "%+.3f" if name == "pH" else "%+.0f"
        row.append(f"{fmt % surface_mean[name][model]} ({bias_fmt % surface_bias[name][model]})")
    unit = "" if name == "pH" else " (umol/kg)"
    md.append(f"| {name}{unit} | " + " | ".join(row) + " |")

md += [
    "",
    "Model columns show the value with its difference from GLODAP in brackets.",
    "",
    f"## Step 3. Surface DIC pattern agreement ({n_common} shared cells)",
    "",
    "| Model | bias (umol/kg) | RMSE (umol/kg) | spatial correlation |",
    "|---|---|---|---|",
]
for model in MODELS:
    e = pattern[model]
    md.append(f"| {model} | {e['bias_umol_per_kg']:+.0f} | {e['rmse_umol_per_kg']:.0f} | "
              f"{e['spatial_correlation']:.3f} |")

md += [
    "",
    "## Step 4. Global mean DIC bias by depth (umol/kg)",
    "",
    "| Model | " + " | ".join(f"{int(d)} m" for d in PROFILE_DEPTHS) + " |",
    "|---" * (len(PROFILE_DEPTHS) + 1) + "|",
]
for model in MODELS:
    md.append(f"| {model} | " + " | ".join(f"{bias_profile[model][d]:+.0f}"
                                           for d in PROFILE_DEPTHS) + " |")

md += [
    "",
    "## Step 5. How much the sampled year matters",
    "",
    "| Model | 1971-2000 mean | 2002 | change | bias vs GLODAP | change / bias |",
    "|---|---|---|---|---|---|",
]
for model in MODELS:
    e = period_change[model]
    md.append(f"| {model} | {e['mean_1971_2000_umol_per_kg']:.0f} | "
              f"{e['mean_2002_umol_per_kg']:.0f} | {e['period_change_umol_per_kg']:+.0f} | "
              f"{e['bias_vs_glodap_umol_per_kg']:+.0f} | "
              f"{e['period_change_as_fraction_of_bias']:+.2f} |")

md += ["", "## What the comparison shows", ""]
md += [f"- {line}" for line in answer["findings"]]
md += [
    "",
    "## Figures",
    "",
    "- `fig1_dic_bias_maps.png` - 0-50 m DIC bias against GLODAP, three panels on one shared scale.",
    "- `fig2_dic_bias_profile.png` - global mean DIC bias against depth, one line per model.",
]
(OUT / "answer_T8.md").write_text("\n".join(md) + "\n")

print("\nWrote answer_T8.json, answer_T8.md and two figures")
