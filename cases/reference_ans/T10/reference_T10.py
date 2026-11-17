"""
Reference answer for case T10 - Global oxygen minimum zones: model assessment
and the ssp245 change.

Run from the repository root with the project's conda environment:

    python cases/reference_ans/T10/reference_T10.py

Writes answer_T10.json, answer_T10.md, fig1_omz_thickness_historical.png and
fig2_omz_thickness_change.png next to this file. Takes roughly half an hour,
almost all of it reading five full-depth monthly archives of o2, thetao and so.
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
from src.tools.preprocessing import regrid_horizontal, select_basin
from src.tools.scientific_calculator import (compute_potential_density_anomaly,
                                             identify_oxygen_minimum_zone,
                                             convert_concentration_to_content)
from src.tools.basic_calculator import (compute_basic_stat, temporal_aggregate, transform_unit,
                                        compute_exceedance_fraction, locate_extreme)

OUT = Path(__file__).resolve().parent
THRESHOLDS = [20.0, 60.0]
MAIN_THRESHOLD = 60.0
BASINS = ["Pacific", "Indian", "Atlantic"]
VOLUME_DIMS = ["depth", "lat", "lon"]

# CESM2 has no historical o2; the run list says exactly which combinations exist.
RUNS = [
    ("GFDL-ESM4",     "historical", "cmip6_historical", ["1971-01", "2000-12"]),
    ("MPI-ESM1-2-LR", "historical", "cmip6_historical", ["1971-01", "2000-12"]),
    ("CESM2",         "ssp245",     "cmip6_ssp245",     ["2071-01", "2100-12"]),
    ("GFDL-ESM4",     "ssp245",     "cmip6_ssp245",     ["2071-01", "2100-12"]),
    ("MPI-ESM1-2-LR", "ssp245",     "cmip6_ssp245",     ["2071-01", "2100-12"]),
]
HISTORICAL_MODELS = ["GFDL-ESM4", "MPI-ESM1-2-LR"]

xr.set_options(keep_attrs=True)

oxygen = {}
oxygen_per_volume = {}
density_note = {}


# ---------------------------------------------------------------------------
# Step 0. The oxygen climatologies, observations and models
# ---------------------------------------------------------------------------
oxygen["WOA23"] = load_climatology_data("woa23", "oxygen", time_slice="annual",
                                        period="1971-2000")
print(f"[0]   WOA23 oxygen: {oxygen['WOA23'].sizes['depth']} levels, "
      f"units {oxygen['WOA23'].attrs['units']}, "
      f"period {oxygen['WOA23'].attrs.get('climatology_period')}")

for model, experiment, dataset, window in RUNS:
    label = f"{model}|{experiment}"
    started = time.time()
    raw = temporal_aggregate(
        load_regular_data(dataset, "o2", model=model, time_range=window, chunk_time=12),
        freq="mean", op="mean").compute()
    temperature = temporal_aggregate(
        load_regular_data(dataset, "thetao", model=model, time_range=window, chunk_time=12),
        freq="mean", op="mean").compute()
    salinity = temporal_aggregate(
        load_regular_data(dataset, "so", model=model, time_range=window, chunk_time=12),
        freq="mean", op="mean").compute()

    # In-situ density at each cell's own pressure.
    density = compute_potential_density_anomaly(temperature, salinity, ref_pressure=None)
    oxygen[label] = transform_unit(raw, from_unit=raw.attrs["units"], to_unit="umol/kg",
                                   density=density)

    oxygen_per_volume[label] = raw
    density_note[label] = oxygen[label].attrs.get("density_assumption", "in-situ density field")
    print(f"[0]   {label:28s} {raw.shape} -> umol/kg  ({time.time() - started:.0f} s)")

SOURCES = ["WOA23"] + [f"{m}|{e}" for m, e, _, _ in RUNS]


# ---------------------------------------------------------------------------
# Step 1. OMZ volume at both thresholds
# ---------------------------------------------------------------------------
volume = {}
for source in SOURCES:
    volume[source] = {}
    for threshold in THRESHOLDS:
        result = compute_exceedance_fraction(oxygen[source], threshold=threshold,
                                             comparison="less", dim=VOLUME_DIMS,
                                             weighting="volume")
        volume[source][threshold] = {
            "volume_m3": float(result["exceeding_amount"]),
            "percent_of_ocean": float(result["fraction_percent"]),
            "criterion": str(result["criterion"]),
        }
    print(f"[1]   {source:28s} " + "  ".join(
        f"<{t:.0f}: {volume[source][t]['volume_m3'] / 1e15:7.3f}e15 m3 "
        f"({volume[source][t]['percent_of_ocean']:5.2f}%)" for t in THRESHOLDS))


# ---------------------------------------------------------------------------
# Step 2. Thickest point at the 60 threshold
# ---------------------------------------------------------------------------
def geographic(coordinates):
    """Latitude and longitude out of locate_extreme, whatever the grid calls them.

    A curvilinear model names its 2-D coordinates after its own convention:
    MPI-ESM1-2-LR uses latitude/longitude where CESM2 and the regular grids use
    lat/lon, and the index pair beside them is not geographic at all.
    """
    def pick(*names):
        for name in names:
            if name in coordinates:
                return float(coordinates[name])
        raise KeyError(f"no geographic coordinate among {names} in {list(coordinates)}")
    return pick("lat", "latitude"), pick("lon", "longitude")


thickness, thickest = {}, {}
for source in SOURCES:
    field, _ = identify_oxygen_minimum_zone(oxygen[source], threshold=MAIN_THRESHOLD)
    thickness[source] = field.compute() if hasattr(field, "compute") else field
    extreme = locate_extreme(thickness[source], op="max")
    latitude, longitude = geographic(extreme["coordinates"])
    thickest[source] = {
        "thickness_m": float(extreme["value"]),
        "latitude": latitude,
        "longitude": longitude,
        "ties": int(extreme["ties"]),
    }
    entry = thickest[source]
    print(f"[2]   {source:28s} max thickness {entry['thickness_m']:7.1f} m at "
          f"{entry['latitude']:+6.2f} N, {entry['longitude']:+7.2f} E")


# ---------------------------------------------------------------------------
# Step 3. Scenario change, for the two models that have both experiments
# ---------------------------------------------------------------------------
scenario_change = {}
for model in HISTORICAL_MODELS:
    scenario_change[model] = {}
    for threshold in THRESHOLDS:
        historical = volume[f"{model}|historical"][threshold]["volume_m3"]
        scenario = volume[f"{model}|ssp245"][threshold]["volume_m3"]
        scenario_change[model][threshold] = {
            "historical_m3": historical,
            "ssp245_m3": scenario,
            "change_percent": 100.0 * (scenario - historical) / historical,
        }
    print(f"[3]   {model:14s} " + "  ".join(
        f"<{t:.0f}: {scenario_change[model][t]['change_percent']:+7.2f}%" for t in THRESHOLDS))


# ---------------------------------------------------------------------------
# Step 4. Oxygen inventory, models only
# ---------------------------------------------------------------------------
inventory = {}
for model, experiment, _, _ in RUNS:
    label = f"{model}|{experiment}"

    content = convert_concentration_to_content(oxygen_per_volume[label])
    inventory[label] = float(compute_basic_stat(content, op="sum").compute())
    print(f"[4]   {label:28s} inventory {inventory[label]:.4e} mol")


# ---------------------------------------------------------------------------
# Step 5. OMZ volume by basin, at the 60 threshold
# ---------------------------------------------------------------------------
basin_volume, basin_definition = {}, {}
for source in ["WOA23"] + [f"{m}|historical" for m in HISTORICAL_MODELS]:
    basin_volume[source] = {}
    for basin in BASINS:
        masked = select_basin(oxygen[source], basin=basin)      # marginal seas folded in
        result = compute_exceedance_fraction(masked, threshold=MAIN_THRESHOLD,
                                             comparison="less", dim=VOLUME_DIMS,
                                             weighting="volume")
        basin_volume[source][basin] = float(result["exceeding_amount"])
        basin_definition[basin] = masked.attrs.get("basin_definition", "")
    print(f"[5]   {source:28s} " + "  ".join(
        f"{b} {basin_volume[source][b] / 1e15:6.3f}e15" for b in BASINS))


# ---------------------------------------------------------------------------
# Step 6. What kind of OMZ, from the observations
# ---------------------------------------------------------------------------
_, observed_mask = identify_oxygen_minimum_zone(oxygen["WOA23"], threshold=MAIN_THRESHOLD)
aou = load_climatology_data("woa23", "aou", time_slice="annual", period="1971-2000")
saturation = load_climatology_data("woa23", "o2sat", time_slice="annual", period="1971-2000")
character = {
    "mean_aou_umol_per_kg": float(compute_basic_stat(aou.where(observed_mask),
                                                     op="volume_weighted_mean",
                                                     dim=VOLUME_DIMS).compute()),
    "mean_percent_saturation": float(compute_basic_stat(saturation.where(observed_mask),
                                                        op="volume_weighted_mean",
                                                        dim=VOLUME_DIMS).compute()),
}
print(f"[6]   inside the observed OMZ: AOU {character['mean_aou_umol_per_kg']:.2f} umol/kg, "
      f"saturation {character['mean_percent_saturation']:.2f}%")

better = min(HISTORICAL_MODELS,
             key=lambda m: abs(volume[f"{m}|historical"][MAIN_THRESHOLD]["volume_m3"]
                               - volume["WOA23"][MAIN_THRESHOLD]["volume_m3"]))
print(f"[7]   closer to the observed OMZ volume at {MAIN_THRESHOLD:.0f}: {better}")


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
INK = "#2f3437"
LAND = "#d9d9d6"

display = {}
for source in SOURCES:
    if source == "WOA23":
        display[source] = thickness[source]
    else:
        display[source] = regrid_horizontal(thickness[source], target="woa23",
                                            method="conservative")
        print(f"[fig] {source:28s} thickness regridded for display")

lon = display["WOA23"]["lon"].values
lat = display["WOA23"]["lat"].values


def draw_map(ax, values, cmap, norm=None, vmin=None, vmax=None):
    ax.set_facecolor(LAND)
    mesh = ax.pcolormesh(lon, lat, values, cmap=cmap, norm=norm, vmin=vmin, vmax=vmax,
                         shading="auto")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_yticks([-90, -45, 0, 45, 90])
    ax.tick_params(colors=INK, labelsize=8.5)
    ax.set_xlabel("Longitude (°E)", color=INK, fontsize=9)
    for spine in ax.spines.values():
        spine.set_color("#c8ccce")
    return mesh


# --- Figure 1: historical OMZ thickness ------------------------------------
historical_panels = [f"{m}|historical" for m in HISTORICAL_MODELS] + ["WOA23"]
stack = np.concatenate([display[s].values[np.isfinite(display[s].values)]
                        for s in historical_panels])
vmax_thickness = float(np.ceil(np.percentile(stack[stack > 0], 99) / 250) * 250)
print(f"[fig1] shared scale 0-{vmax_thickness:.0f} m "
      f"(99th percentile of non-zero thickness = {np.percentile(stack[stack > 0], 99):.0f})")

fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.6))
fig.subplots_adjust(left=0.045, right=0.99, bottom=0.235, top=0.80, wspace=0.13)
for ax, source in zip(axes, historical_panels):
    mesh = draw_map(ax, display[source].values, "magma_r", vmin=0, vmax=vmax_thickness)
    label = "WOA23 (observed)" if source == "WOA23" else source.split("|")[0]
    ax.set_title(f"{label}\nOMZ volume {volume[source][MAIN_THRESHOLD]['volume_m3'] / 1e15:.3g}"
                 f"×10¹⁵ m³,   max {thickest[source]['thickness_m']:.0f} m",
                 fontsize=10.5, color=INK)
axes[0].set_ylabel("Latitude (°N)", color=INK, fontsize=9)

cax = fig.add_axes([0.30, 0.10, 0.40, 0.035])
bar = fig.colorbar(mesh, cax=cax, orientation="horizontal", extend="max")
bar.set_label(f"OMZ thickness at O₂ < {MAIN_THRESHOLD:.0f} µmol kg⁻¹ (m)",
              color=INK, fontsize=10)
bar.ax.tick_params(colors=INK, labelsize=9)
fig.suptitle("Oxygen minimum zone thickness, historical 1971–2000 against WOA23\n"
             "computed on each source's own vertical levels; model panels regridded to 1° "
             "for display only",
             fontsize=12.5, color=INK, y=0.975)
fig.savefig(OUT / "fig1_omz_thickness_historical.png", dpi=150, facecolor="white")
plt.close(fig)

# --- Figure 2: the scenario picture ----------------------------------------
change_field = {m: display[f"{m}|ssp245"] - display[f"{m}|historical"]
                for m in HISTORICAL_MODELS}
stack = np.concatenate([np.abs(change_field[m].values[np.isfinite(change_field[m].values)])
                        for m in HISTORICAL_MODELS])
vmax_change = float(np.ceil(np.percentile(stack, 98) / 50) * 50)
print(f"[fig2] change scale +/-{vmax_change:.0f} m "
      f"(98th percentile of |change| = {np.percentile(stack, 98):.0f})")

fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.9))
fig.subplots_adjust(left=0.045, right=0.99, bottom=0.27, top=0.78, wspace=0.13)
for ax, model in zip(axes[:2], HISTORICAL_MODELS):
    mesh_change = draw_map(ax, change_field[model].values, "RdBu_r",
                           norm=TwoSlopeNorm(vmin=-vmax_change, vcenter=0.0, vmax=vmax_change))
    ax.set_title(f"{model}\nssp245 − historical,   volume "
                 f"{scenario_change[model][MAIN_THRESHOLD]['change_percent']:+.1f}%",
                 fontsize=10.5, color=INK)
axes[0].set_ylabel("Latitude (°N)", color=INK, fontsize=9)

mesh_cesm = draw_map(axes[2], display["CESM2|ssp245"].values, "magma_r",
                     vmin=0, vmax=vmax_thickness)
axes[2].set_title("CESM2 — ssp245 thickness only\n"
                  "no historical o₂ is published, so no change can be formed",
                  fontsize=10.5, color=INK)

cax = fig.add_axes([0.075, 0.115, 0.52, 0.035])
bar = fig.colorbar(mesh_change, cax=cax, orientation="horizontal", extend="both")
bar.set_label("Change in OMZ thickness, ssp245 2071–2100 − historical 1971–2000 (m)",
              color=INK, fontsize=10)
bar.ax.tick_params(colors=INK, labelsize=9)

cax2 = fig.add_axes([0.695, 0.115, 0.255, 0.035])
bar2 = fig.colorbar(mesh_cesm, cax=cax2, orientation="horizontal", extend="max")
bar2.set_label("CESM2 ssp245 OMZ thickness (m)", color=INK, fontsize=10)
bar2.ax.tick_params(colors=INK, labelsize=9)

fig.suptitle(f"Oxygen minimum zone under ssp245, threshold O₂ < {MAIN_THRESHOLD:.0f} "
             "µmol kg⁻¹\n"
             "the two models with both experiments show the change; CESM2 can only show its "
             "scenario state",
             fontsize=12.5, color=INK, y=0.975)
fig.savefig(OUT / "fig2_omz_thickness_change.png", dpi=150, facecolor="white")
plt.close(fig)


# ---------------------------------------------------------------------------
# Numeric answer
# ---------------------------------------------------------------------------
def sig3(x):
    return float(f"{x:.3g}")


answer = {
    "case_id": "T10",
    "title": "Global oxygen minimum zones: model assessment and the ssp245 change",
    "sources": SOURCES,
    "method": {
        "thresholds_umol_per_kg": THRESHOLDS,
        "criterion": "oxygen < threshold (strictly less than), for each listed threshold",
        "units": "model o2 converted from mol m-3 to umol/kg with an in-situ density field from "
                 "each model's own thetao and so; a constant moves the third significant figure",
        "density_note": density_note,
        "observations": "WOA23 1971-2000 climate normal, NOT the all-period composite; the two "
                        "differ by about 5 per cent in OMZ volume",
        "grid": "every volume, inventory and thickness on the source's own grid and own vertical "
                "levels; regridding is horizontal and for the figures only",
        "thickness": "the OMZ thickness sums the thickness of each sub-threshold cell, taken for the "
                     "models as volcello / areacello so that partial cells at the sea floor count only "
                     "the water they hold; WOA23 has no partial cells and uses its layer bounds",
        "inventory": "integrated from the model's native mol m-3 field, so no density enters "
                     "the volume integral at all; the umol/kg field is not used here because "
                     "converting it back would apply a constant density",
        "marginal_seas": "folded into their parent ocean (select_basin default)",
        "southern_and_arctic_oceans": "separate basins in the WOA23 mask, not included in the Atlantic, "
                                      "Indian or Pacific numbers",
        "cesm2": "publishes no historical o2, so it appears only in the ssp245 results",
    },

    "omz_volume_1e15_m3": {s: {str(int(t)): sig3(volume[s][t]["volume_m3"] / 1e15)
                               for t in THRESHOLDS} for s in SOURCES},
    "omz_percent_of_ocean_volume": {s: {str(int(t)): round(volume[s][t]["percent_of_ocean"], 2)
                                        for t in THRESHOLDS} for s in SOURCES},
    "max_thickness": {s: {"thickness_m": round(thickest[s]["thickness_m"], 1),
                          "latitude": round(thickest[s]["latitude"], 2),
                          "longitude": round(thickest[s]["longitude"], 2)} for s in SOURCES},
    "scenario_volume_change_percent": {
        m: {str(int(t)): round(scenario_change[m][t]["change_percent"], 2) for t in THRESHOLDS}
        for m in HISTORICAL_MODELS},
    "oxygen_inventory_mol": {k: sig3(v) for k, v in inventory.items()},
    "omz_volume_by_basin_1e15_m3": {s: {b: sig3(basin_volume[s][b] / 1e15) for b in BASINS}
                                    for s in basin_volume},
    "observed_omz_character": {"mean_aou_umol_per_kg": round(character["mean_aou_umol_per_kg"], 1),
                               "mean_percent_saturation": round(character["mean_percent_saturation"], 1)},
    "closer_to_observed_volume": better,
    "basin_definitions": basin_definition,

    "full_precision": {"volume": volume, "max_thickness": thickest,
                       "scenario_change": scenario_change, "inventory": inventory,
                       "basin_volume": basin_volume, "character": character},

    "figures": ["fig1_omz_thickness_historical.png", "fig2_omz_thickness_change.png"],
    "tools_used": [
        "load_regular_data", "load_climatology_data", "temporal_aggregate",
        "compute_potential_density_anomaly", "transform_unit", "identify_oxygen_minimum_zone",
        "compute_exceedance_fraction", "convert_concentration_to_content", "compute_basic_stat",
        "locate_extreme", "select_basin", "regrid_horizontal",
    ],
}

observed60 = volume["WOA23"][MAIN_THRESHOLD]["volume_m3"] / 1e15
answer["findings"] = [
    f"Both historical models overestimate the OMZ. At the {MAIN_THRESHOLD:.0f} threshold the "
    f"observed volume is {observed60:.3g}e15 m3 against "
    + ", ".join(f"{m} {volume[f'{m}|historical'][MAIN_THRESHOLD]['volume_m3'] / 1e15:.3g}e15"
                for m in HISTORICAL_MODELS)
    + f"; {better} is the closer of the two. Overestimating the OMZ is a documented common bias.",

    "Under ssp245 the two models do not agree on the sign of the change: "
    + "; ".join(f"{m} {scenario_change[m][MAIN_THRESHOLD]['change_percent']:+.1f}% at "
                f"{MAIN_THRESHOLD:.0f} and "
                f"{scenario_change[m][THRESHOLDS[0]]['change_percent']:+.1f}% at "
                f"{THRESHOLDS[0]:.0f}" for m in HISTORICAL_MODELS)
    + ". The threshold matters as much as the model: the suboxic core and the hypoxic envelope "
    "need not move together, which is why both thresholds are asked for. Each change compares two "
    "periods of one continuous run (r1i1p1f1 in both experiments), but it also contains whatever "
    "drift the model's pre-industrial control still carries, which cannot be removed without a "
    "piControl run.",

    f"Inside the observed OMZ the mean apparent oxygen utilization is "
    f"{character['mean_aou_umol_per_kg']:.1f} micromoles per kilogram and the mean saturation is "
    f"{character['mean_percent_saturation']:.1f} per cent. A high AOU with a low saturation says "
    "the zone is set by the respiration of sinking organic matter in water that has been out of "
    "contact with the atmosphere for a long time, rather than by cold, poorly ventilated water "
    "that was already low in oxygen when it left the surface.",

    "CESM2 could not be assessed against observations or given a scenario change at all, because "
    "the historical experiment publishes no o2 for it. Its ssp245 OMZ is reported and drawn, but "
    "it cannot be told whether that state is a change or a pre-existing bias.",
]

(OUT / "answer_T10.json").write_text(json.dumps(answer, indent=2) + "\n")

def shown(label):
    """A source label safe to put in a markdown cell.

    The internal labels are 'MODEL|experiment'; a bare pipe inside a cell would
    be read as a column separator and break the table.
    """
    return label.replace("|", " / ")


md = [
    "# Reference answer - T10: global oxygen minimum zones and the ssp245 change",
    "",
    "Thresholds 20 and 60 umol/kg, strictly less than. Model oxygen converted from mol m-3 with an",
    "in-situ density field from each model's own thetao and so. Observations are the WOA23",
    "1971-2000 climate normal. Every volume and inventory on the source's own grid and levels.",
    "",
    "CESM2 publishes no historical o2, so it appears only under ssp245.",
    "",
    "## Step 1. OMZ volume (1e15 m3)",
    "",
    "| Source | < 20 umol/kg | < 60 umol/kg |",
    "|---|---|---|",
]
for source in SOURCES:
    md.append(f"| {shown(source)} | {volume[source][20.0]['volume_m3'] / 1e15:.3g} | "
              f"{volume[source][60.0]['volume_m3'] / 1e15:.3g} |")

md += [
    "",
    "## Step 2. Thickest point at the 60 threshold",
    "",
    "| Source | max thickness (m) | latitude | longitude |",
    "|---|---|---|---|",
]
for source in SOURCES:
    e = thickest[source]
    md.append(f"| {shown(source)} | {e['thickness_m']:.1f} | {e['latitude']:+.2f} | "
              f"{e['longitude']:+.2f} |")

md += [
    "",
    "## Step 3. Scenario change in OMZ volume, ssp245 minus historical",
    "",
    "| Model | < 20 umol/kg | < 60 umol/kg |",
    "|---|---|---|",
]
for model in HISTORICAL_MODELS:
    md.append(f"| {model} | {scenario_change[model][20.0]['change_percent']:+.2f}% | "
              f"{scenario_change[model][60.0]['change_percent']:+.2f}% |")

md += [
    "",
    "## Step 4. Global dissolved oxygen inventory",
    "",
    "| Model and experiment | inventory (mol) |",
    "|---|---|",
]
for label, value in inventory.items():
    md.append(f"| {shown(label)} | {value:.3g} |")

md += [
    "",
    "## Step 5. OMZ volume by basin at the 60 threshold (1e15 m3)",
    "",
    "| Source | Pacific | Indian | Atlantic |",
    "|---|---|---|---|",
]
for source in basin_volume:
    md.append(f"| {shown(source)} | " + " | ".join(f"{basin_volume[source][b] / 1e15:.3g}"
                                            for b in BASINS) + " |")

md += [
    "",
    "Marginal seas are folded into their parent ocean. The Southern and Arctic Oceans are separate "
    "basins in the WOA23 mask and are not part of these numbers.",
    "",
    "## Step 6. What kind of OMZ, from the observations",
    "",
    f"Inside the observed OMZ: mean AOU {character['mean_aou_umol_per_kg']:.1f} umol/kg, "
    f"mean saturation {character['mean_percent_saturation']:.1f}%.",
    "",
    "## What the comparison shows",
    "",
]
md += [f"- {line}" for line in answer["findings"]]
md += [
    "",
    "## Figures",
    "",
    "- `fig1_omz_thickness_historical.png` - OMZ thickness for the two historical models and",
    "  WOA23 on one shared scale.",
    "- `fig2_omz_thickness_change.png` - the ssp245 minus historical change for the two models",
    "  that have both experiments, plus CESM2's ssp245 state on its own scale.",
]
(OUT / "answer_T10.md").write_text("\n".join(md) + "\n")

print("\nWrote answer_T10.json, answer_T10.md and two figures")
