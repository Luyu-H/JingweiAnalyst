"""
Reference answer for case T3 - Multi-model assessment of simulated ENSO
amplitude in CMIP6.

Run from the repository root with the project's conda environment:

    source /home/lorraine/miniconda3/bin/activate base && conda activate jingwei_QA
    PYTHONPATH=. python cases/reference_ans/T3/reference_T3.py

Writes answer_T3.json, answer_T3.md, fig1_seasonal_amplitude.png and
fig2_nino34_panels.png next to this file.

Every physical quantity is produced by a tool from src/tools; nothing here
re-implements an index, an anomaly, a standard deviation or an exceedance
fraction by hand.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from src.tools.registry import auto_import_tools
auto_import_tools()
from src.tools.data_loader import load_regular_data
from src.tools.scientific_calculator import compute_enso_nino_index
from src.tools.basic_calculator import (compute_basic_stat, temporal_aggregate,
                                        compute_exceedance_fraction)

OUT = Path(__file__).resolve().parent
WINDOW = ["1950-01", "2014-12"]
MODELS = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
SOURCES = ["ERSSTv5"] + MODELS          # observations first, then the models
N_MONTHS = 780
THRESHOLD = 1.0

xr.set_options(keep_attrs=True)


# ---------------------------------------------------------------------------
# Step 1. Nino 3.4 index for each source, on its own grid and own climatology
# ---------------------------------------------------------------------------
index = {}

sst = load_regular_data("ersst", "sst", time_range=WINDOW)
index["ERSSTv5"] = compute_enso_nino_index(sst, nino_region="nino3.4",
                                           base_period=WINDOW).compute()

for model in MODELS:
    # The region is cut out with a mask and area-weighted with the model's own
    # areacello, so the curvilinear grids need no regridding and are never
    # sliced by coordinate value.
    tos = load_regular_data("cmip6_historical", "tos", model=model, time_range=WINDOW)
    index[model] = compute_enso_nino_index(tos, nino_region="nino3.4",
                                           base_period=WINDOW).compute()

for source in SOURCES:
    series = index[source]
    assert series.sizes["time"] == N_MONTHS, f"{source}: {series.sizes['time']} months"
    print(f"[1]   {source:14s} n={series.sizes['time']}  "
          f"baseline={series.attrs['climatology_baseline']}")


# ---------------------------------------------------------------------------
# Step 2. Overall amplitude, bias against observations and ranking
# ---------------------------------------------------------------------------
std = {s: float(compute_basic_stat(index[s], op="std", dim="time").compute()) for s in SOURCES}
obs_std = std["ERSSTv5"]
bias = {m: std[m] - obs_std for m in MODELS}
ratio = {m: std[m] / obs_std for m in MODELS}
ranking = sorted(MODELS, key=lambda m: abs(bias[m]))

for source in SOURCES:
    extra = (f"  bias {bias[source]:+.4f}  ratio {ratio[source]:.3f}"
             if source in MODELS else "  (reference)")
    print(f"[2]   {source:14s} std = {std[source]:.4f} degC{extra}")
print(f"[2]   ranking by |amplitude bias|, closest first: {' < '.join(ranking)}")


# ---------------------------------------------------------------------------
# Step 3. Seasonal phase locking: one standard deviation per calendar month
# ---------------------------------------------------------------------------
MONTH_NAMES = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

monthly_std, phase_locking = {}, {}
for source in SOURCES:
    # freq='climatology' groups by calendar month; op='std' is deliberately
    # unweighted, since this is a spread across years, not a time average.
    series = temporal_aggregate(index[source], freq="climatology", op="std").compute()
    values = series.values.astype(float)
    assert values.shape == (12,)
    monthly_std[source] = values

    peak = int(np.argmax(values))
    trough = int(np.argmin(values))
    phase_locking[source] = {
        "peak_month": peak + 1,
        "peak_month_name": MONTH_NAMES[peak],
        "peak_std_degC": float(values[peak]),
        "trough_month": trough + 1,
        "trough_month_name": MONTH_NAMES[trough],
        "trough_std_degC": float(values[trough]),
        "locking_ratio": float(values[peak] / values[trough]),
    }
    print(f"[3]   {source:14s} peak {MONTH_NAMES[peak]} {values[peak]:.3f}  "
          f"trough {MONTH_NAMES[trough]} {values[trough]:.3f}  "
          f"ratio {values[peak] / values[trough]:.2f}")


# ---------------------------------------------------------------------------
# Step 4. Frequency of large events
# ---------------------------------------------------------------------------
exceedance = {}
for source in SOURCES:
    result = compute_exceedance_fraction(index[source], threshold=THRESHOLD,
                                         comparison="abs_greater", dim="time",
                                         weighting="none")
    exceedance[source] = {
        "percent": float(result["fraction_percent"]),
        "n_exceeding": int(result["n_exceeding"]),
        "n_valid": int(result["n_valid"]),
        "criterion": str(result["criterion"]),
    }
    print(f"[4]   {source:14s} |anom| > {THRESHOLD} degC in "
          f"{exceedance[source]['n_exceeding']:3d}/{exceedance[source]['n_valid']} months "
          f"= {exceedance[source]['percent']:.2f}%")


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
# Okabe-Ito hues for the three models, validated for colour-vision deficiency;
# the observations carry neutral ink and a heavier line, so the reference is
# never identified by hue alone.
COLOUR = {"ERSSTv5": "#2f3437", "CESM2": "#0072B2",
          "GFDL-ESM4": "#009E73", "MPI-ESM1-2-LR": "#E69F00"}
INK = "#2f3437"
GRID = "#eceff0"
months = np.arange(1, 13)

# --- Figure 1: the seasonal cycle of amplitude -----------------------------
fig, ax = plt.subplots(figsize=(9.2, 5.8), constrained_layout=True)
for source in SOURCES:
    observed = source == "ERSSTv5"
    ax.plot(months, monthly_std[source], color=COLOUR[source],
            linewidth=3.0 if observed else 1.9,
            marker="o", markersize=7 if observed else 5.5,
            zorder=5 if observed else 4,
            label=f"{source} (observed)" if observed else source)
    peak = phase_locking[source]["peak_month"]
    ax.plot(peak, monthly_std[source][peak - 1], marker="o",
            markersize=12 if observed else 10, markerfacecolor="none",
            markeredgecolor=COLOUR[source], markeredgewidth=1.8, zorder=6)

ax.set_xticks(months)
ax.set_xticklabels(MONTH_NAMES)
ax.set_xlim(0.6, 12.4)
ax.set_ylim(0, max(v.max() for v in monthly_std.values()) * 1.18)
ax.set_xlabel("Calendar month", color=INK)
ax.set_ylabel("Nino 3.4 anomaly standard deviation (°C)", color=INK)
ax.set_title("Seasonal cycle of ENSO amplitude, 1950–2014\n"
             "standard deviation within each calendar month; "
             "a ring marks each source's peak month",
             fontsize=12, color=INK)
ax.legend(loc="upper center", frameon=False, fontsize=10, ncol=2, labelcolor=INK)
ax.grid(axis="y", color=GRID, linewidth=0.9)
ax.set_axisbelow(True)
ax.tick_params(colors=INK)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
for side in ("left", "bottom"):
    ax.spines[side].set_color("#c8ccce")

fig.savefig(OUT / "fig1_seasonal_amplitude.png", dpi=160, facecolor="white")
plt.close(fig)

# --- Figure 2: the four series on one shared vertical scale ----------------
extreme = max(float(np.nanmax(np.abs(index[s].values))) for s in SOURCES)
# Round up to the next half degree with a little headroom, so the shared
# scale contains every series without leaving the panels mostly empty.
limit = np.ceil((extreme + 0.15) * 2) / 2
time = pd.to_datetime(index["ERSSTv5"]["time"].values)

fig, axes = plt.subplots(4, 1, figsize=(13.5, 10.0), sharex=True, sharey=True)
fig.subplots_adjust(left=0.075, right=0.985, bottom=0.075, top=0.90, hspace=0.18)

for ax, source in zip(axes, SOURCES):
    values = index[source].values
    ax.axhline(0.0, color="#b9bec1", linewidth=0.8, zorder=1)
    for level in (THRESHOLD, -THRESHOLD):
        ax.axhline(level, color="#8c9295", linewidth=0.9, linestyle="--", zorder=1)
    ax.fill_between(time, 0.0, values, color=COLOUR[source], alpha=0.30,
                    linewidth=0, zorder=2)
    ax.plot(time, values, color=COLOUR[source], linewidth=0.9, zorder=3)

    # A white pad keeps both labels readable where a large event reaches up
    # into the top of the panel, as CESM2's do.
    pad = dict(facecolor="white", edgecolor="none", boxstyle="square,pad=0.25")
    label = f"{source} (observed)" if source == "ERSSTv5" else source
    ax.text(0.008, 0.94, label, transform=ax.transAxes, fontsize=11.5,
            color=COLOUR[source], va="top", ha="left", fontweight="bold",
            bbox=pad, zorder=6)
    ax.text(0.992, 0.94,
            f"σ = {std[source]:.2f} °C     "
            f"|anomaly| > {THRESHOLD:.1f} °C in {exceedance[source]['percent']:.1f}% of months",
            transform=ax.transAxes, fontsize=10, color=INK, va="top", ha="right",
            bbox=pad, zorder=6)

    ax.set_ylim(-limit, limit)
    ax.set_ylabel("Anomaly (°C)", color=INK, fontsize=10)
    ax.grid(axis="y", color=GRID, linewidth=0.9)
    ax.set_axisbelow(True)
    ax.tick_params(colors=INK)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#c8ccce")

axes[-1].set_xlabel("Year", color=INK)
fig.suptitle("Nino 3.4 anomaly, 1950–2014: ERSSTv5 and three CMIP6 historical simulations\n"
             f"every panel on the identical vertical scale ±{limit:.2f} °C, "
             "each source on its own grid and its own 1950–2014 monthly climatology",
             fontsize=13, color=INK, y=0.975)

fig.savefig(OUT / "fig2_nino34_panels.png", dpi=160, facecolor="white")
plt.close(fig)


# ---------------------------------------------------------------------------
# Numeric answer
# ---------------------------------------------------------------------------
answer = {
    "case_id": "T3",
    "title": "Multi-model assessment of simulated ENSO amplitude in CMIP6",
    "window": "1950-01 to 2014-12",
    "n_months": N_MONTHS,
    "sources": SOURCES,
    "method": {
        "region": index["ERSSTv5"].attrs["region_bounds"],
        "weighting": "area-weighted on each source's own native grid and own cell areas; no regridding",
        "baseline": "each source's own monthly climatology over 1950-01 to 2014-12",
        "std_convention": "population standard deviation, ddof = 0",
        "exceedance_criterion": exceedance["ERSSTv5"]["criterion"],
    },

    "nino34_std_degC": {s: round(std[s], 2) for s in SOURCES},
    "amplitude_bias_degC": {m: round(bias[m], 2) for m in MODELS},
    "amplitude_ratio": {m: round(ratio[m], 2) for m in MODELS},
    "ranking_by_abs_amplitude_bias": ranking,

    "monthly_std_degC": {s: {MONTH_NAMES[i]: round(float(monthly_std[s][i]), 2) for i in range(12)}
                         for s in SOURCES},
    "phase_locking": {s: {**phase_locking[s],
                          "peak_std_degC": round(phase_locking[s]["peak_std_degC"], 2),
                          "trough_std_degC": round(phase_locking[s]["trough_std_degC"], 2),
                          "locking_ratio": round(phase_locking[s]["locking_ratio"], 2)}
                      for s in SOURCES},

    "exceedance_percent_abs_gt_1": {s: round(exceedance[s]["percent"], 1) for s in SOURCES},
    "exceedance_count_abs_gt_1": {s: exceedance[s]["n_exceeding"] for s in SOURCES},

    "full_precision": {
        "nino34_std_degC": std,
        "amplitude_bias_degC": bias,
        "amplitude_ratio": ratio,
        "monthly_std_degC": {s: [float(v) for v in monthly_std[s]] for s in SOURCES},
        "phase_locking": phase_locking,
        "exceedance_percent": {s: exceedance[s]["percent"] for s in SOURCES},
    },

    "figures": ["fig1_seasonal_amplitude.png", "fig2_nino34_panels.png"],
    "method_notes": [
        "Each source's index is area-weighted over 5S-5N / 170W-120W on its own native grid with "
        "its own cell areas; CESM2 and MPI-ESM1-2-LR are curvilinear and are masked, never sliced "
        "by coordinate value and never regridded.",
        "Each source is referenced to its own monthly climatology over the same 1950-2014 window, "
        "so a difference in mean state cannot masquerade as an amplitude error.",
        "Standard deviations use ddof = 0, both overall and within each calendar month.",
        "The exceedance test is strictly greater than 1.0 degC in absolute value.",
    ],
    "tools_used": [
        "load_regular_data", "compute_enso_nino_index", "compute_basic_stat",
        "temporal_aggregate", "compute_exceedance_fraction",
    ],
}

(OUT / "answer_T3.json").write_text(json.dumps(answer, indent=2) + "\n")


md = [
    "# Reference answer - T3: multi-model assessment of simulated ENSO amplitude",
    "",
    f"Window 1950-01 to 2014-12, {N_MONTHS} months. Nino 3.4 index area-weighted over "
    "5S-5N / 170W-120W",
    "on each source's own native grid, anomalies against that source's own 1950-2014 monthly",
    "climatology. Standard deviations use ddof = 0.",
    "",
    "## Step 2. Overall amplitude",
    "",
    "| Source | std (degC) | bias vs obs (degC) | ratio to obs |",
    "|---|---|---|---|",
    f"| ERSSTv5 (observed) | {std['ERSSTv5']:.2f} | - | - |",
]
for m in MODELS:
    md.append(f"| {m} | {std[m]:.2f} | {bias[m]:+.2f} | {ratio[m]:.2f} |")
md += [
    "",
    f"Ranking by |amplitude bias|, closest first: **{' < '.join(ranking)}**.",
    "",
    "## Step 3. Seasonal phase locking",
    "",
    "Standard deviation within each calendar month (degC):",
    "",
    "| Source | " + " | ".join(MONTH_NAMES) + " |",
    "|---" * 13 + "|",
]
for s in SOURCES:
    md.append(f"| {s} | " + " | ".join(f"{v:.2f}" for v in monthly_std[s]) + " |")
md += [
    "",
    "| Source | peak month | peak (degC) | trough month | trough (degC) | peak/trough |",
    "|---|---|---|---|---|---|",
]
for s in SOURCES:
    p = phase_locking[s]
    md.append(f"| {s} | {p['peak_month_name']} | {p['peak_std_degC']:.2f} | "
              f"{p['trough_month_name']} | {p['trough_std_degC']:.2f} | {p['locking_ratio']:.2f} |")
md += [
    "",
    "## Step 4. Months with |anomaly| > 1.0 degC",
    "",
    "| Source | months | percentage |",
    "|---|---|---|",
]
for s in SOURCES:
    md.append(f"| {s} | {exceedance[s]['n_exceeding']} / {exceedance[s]['n_valid']} | "
              f"{exceedance[s]['percent']:.1f}% |")
md += [
    "",
    "## Figures",
    "",
    "- `fig1_seasonal_amplitude.png` - the 12-month amplitude cycle, four curves on one axis,",
    "  each source's peak month ringed.",
    f"- `fig2_nino34_panels.png` - the four series stacked on one identical vertical scale "
    f"(+/-{limit:.2f} degC),",
    "  with the +/-1.0 degC levels marked and each panel annotated with its std and exceedance.",
]
(OUT / "answer_T3.md").write_text("\n".join(md) + "\n")

print("\nWrote answer_T3.json, answer_T3.md, fig1_seasonal_amplitude.png, fig2_nino34_panels.png")
