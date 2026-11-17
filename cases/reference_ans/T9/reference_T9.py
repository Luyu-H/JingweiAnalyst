"""
Reference answer for case T9 - Multi-model assessment of the PDO pattern and its decadal statistics.

Run from the repository root with the project's conda environment:

    python cases/reference_ans/T9/reference_T9.py

Writes answer_T9.json, answer_T9.md, fig1_pdo_patterns.png and fig2_pdo_timeseries.png next to this file.
"""

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from src.tools.registry import auto_import_tools
auto_import_tools()
from src.tools.data_loader import load_regular_data
from src.tools.preprocessing import regrid_horizontal
from src.tools.scientific_calculator import compute_pdo_index
from src.tools.basic_calculator import compute_basic_stat, temporal_aggregate, analyze_correlation

OUT = Path(__file__).resolve().parent
MODELS = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
SOURCES = ["ERSSTv5"] + MODELS
WINDOW = ["1950-01", "2014-12"]
SMOOTH_MONTHS = 108
MERGE_YEARS = 3

xr.set_options(keep_attrs=True)


# ---------------------------------------------------------------------------
# Step 1. The PDO of each source, all on one 1-degree grid
# ---------------------------------------------------------------------------
index, pattern, variance = {}, {}, {}

for source in SOURCES:
    started = time.time()
    if source == "ERSSTv5":
        field = load_regular_data("ersst", "sst", time_range=WINDOW)
    else:
        field = load_regular_data("cmip6_historical", "tos", model=source, time_range=WINDOW)
    native = field.attrs.get("grid_type", "regular")
    on_grid = regrid_horizontal(field, target="1deg", method="bilinear")

    pc, eigenvector, explained, regression = compute_pdo_index(on_grid, standardize=True)
    index[source] = pc.compute()
    pattern[source] = regression
    variance[source] = explained
    print(f"[1]   {source:14s} native {native:12s} -> 1 degree, "
          f"EOF1 {explained:.2f}%  ({time.time() - started:.0f} s)")


# ---------------------------------------------------------------------------
# Step 2-3. Explained variance, amplitude and persistence
# ---------------------------------------------------------------------------
def amplitude(field):
    """Area-weighted RMS of the regression pattern, in degC per standard deviation."""
    mean_square = compute_basic_stat(field ** 2, op="area_weighted_mean",
                                     dim=["lat", "lon"]).compute()
    return float(np.sqrt(float(mean_square)))


statistics = {}
for source in SOURCES:
    series = index[source]
    # Lag-1-year autocorrelation: the monthly series against itself 12 steps on.
    lagged = analyze_correlation(series, series, dim="time", method="pearson", lag=12)
    statistics[source] = {
        "explained_variance_percent": variance[source],
        "amplitude_degC_per_sd": amplitude(pattern[source]),
        "lag1yr_autocorrelation": float(lagged["correlation"]),
        "index_sd": float(series.std()),
    }
    entry = statistics[source]
    print(f"[2-3] {source:14s} var {entry['explained_variance_percent']:5.2f}%  "
          f"amplitude {entry['amplitude_degC_per_sd']:.4f} degC/sd  "
          f"lag-1yr r {entry['lag1yr_autocorrelation']:+.4f}  "
          f"(index sd {entry['index_sd']:.3f}, 1 by construction)")


# ---------------------------------------------------------------------------
# Step 4. The reconstruction against the published index
# ---------------------------------------------------------------------------
published = load_regular_data("cpc_indices", product="pdo", time_range=WINDOW)
recon_vs_published = float(analyze_correlation(index["ERSSTv5"], published, dim="time",
                                               method="pearson")["correlation"])
print(f"[4]   ERSSTv5 reconstruction vs published PDO: r = {recon_vs_published:.4f} "
      f"over {published.sizes['time']} months "
      f"(published is ERSST v6, the gridded field is v5)")


# ---------------------------------------------------------------------------
# Step 5. Pattern agreement with the observations
# ---------------------------------------------------------------------------
pattern_correlation = {}
for model in MODELS:
    pattern_correlation[model] = float(analyze_correlation(
        pattern[model], pattern["ERSSTv5"], dim=["lat", "lon"],
        method="pearson", weighting="area")["correlation"])
    print(f"[5]   {model:14s} pattern correlation with observations = "
          f"{pattern_correlation[model]:+.4f}")

best = max(MODELS, key=lambda m: pattern_correlation[m])


# ---------------------------------------------------------------------------
# Step 6. Observed phase transitions
# ---------------------------------------------------------------------------
smoothed = temporal_aggregate(index["ERSSTv5"], freq="rolling", op="mean",
                              window=SMOOTH_MONTHS, center=True)
values = smoothed.values
stamps = pd.to_datetime(smoothed["time"].values)
valid = np.isfinite(values)
values, stamps = values[valid], stamps[valid]

crossings = [stamps[i].year for i in range(1, len(values))
             if np.sign(values[i]) != np.sign(values[i - 1])]
transitions = []
for year in crossings:
    if not transitions or year - transitions[-1] >= MERGE_YEARS:
        transitions.append(year)
print(f"[6]   sign changes of the {SMOOTH_MONTHS // 12}-year smoothed index: {crossings}")
print(f"[6]   merged transitions (< {MERGE_YEARS} yr apart combined): {transitions}")


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
INK = "#2f3437"
LAND = "#d9d9d6"
COLOUR = {"ERSSTv5": "#2f3437", "CESM2": "#0072B2",
          "GFDL-ESM4": "#009E73", "MPI-ESM1-2-LR": "#E69F00"}
PUBLISHED = "#c1543a"

# --- Figure 1: the regression patterns -------------------------------------
stack = np.concatenate([np.abs(pattern[s].values[np.isfinite(pattern[s].values)])
                        for s in SOURCES])
vmax = float(np.ceil(np.percentile(stack, 98) * 10) / 10)
print(f"[fig1] shared diverging scale +/-{vmax:.1f} degC/sd "
      f"(98th percentile of |pattern| = {np.percentile(stack, 98):.2f})")

lon = pattern["ERSSTv5"]["lon"].values
lat = pattern["ERSSTv5"]["lat"].values

order = np.argsort(np.where(lon < 0, lon + 360.0, lon))
lon_plot = np.where(lon < 0, lon + 360.0, lon)[order]

fig, axes = plt.subplots(2, 2, figsize=(14.0, 7.4))
fig.subplots_adjust(left=0.06, right=0.985, bottom=0.15, top=0.845, wspace=0.13, hspace=0.42)
for ax, source in zip(axes.flat, SOURCES):
    ax.set_facecolor(LAND)
    mesh = ax.pcolormesh(lon_plot, lat, pattern[source].values[:, order], cmap="RdBu_r",
                         norm=TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax), shading="auto")
    if source == "ERSSTv5":
        title = f"ERSSTv5 (observed)\nEOF1 {variance[source]:.1f}% of variance"
    else:
        title = (f"{source}\nEOF1 {variance[source]:.1f}%,   "
                 f"pattern r = {pattern_correlation[source]:.3f}")
    ax.set_title(title, fontsize=10.5, color=INK)
    ax.set_xticks([120, 160, 200, 240])
    ax.set_xticklabels(["120°E", "160°E", "160°W", "120°W"])
    ax.set_yticks([20, 35, 50, 65])
    ax.tick_params(colors=INK, labelsize=8.5)
    ax.set_ylabel("Latitude (°N)", color=INK, fontsize=9)
    for spine in ax.spines.values():
        spine.set_color("#c8ccce")
for ax in axes[1]:                    # only the bottom row carries the x label
    ax.set_xlabel("Longitude", color=INK, fontsize=9)

cax = fig.add_axes([0.28, 0.065, 0.44, 0.025])
bar = fig.colorbar(mesh, cax=cax, orientation="horizontal", extend="both")
bar.set_label("SST anomaly regressed on the standardized PDO index "
              "(°C per standard deviation)", color=INK, fontsize=10)
bar.ax.tick_params(colors=INK, labelsize=9)
fig.suptitle("PDO spatial pattern, 1950–2014, 20°N–70°N / 110°E–100°W\n"
             "every source regridded to a common 1° grid before the EOF; positive phase is "
             "cool in the central and western North Pacific",
             fontsize=12.5, color=INK, y=0.975)
fig.savefig(OUT / "fig1_pdo_patterns.png", dpi=150, facecolor="white")
plt.close(fig)

# --- Figure 2: the index series -------------------------------------------
fig, ax = plt.subplots(figsize=(14.5, 6.2), constrained_layout=True)
ax.axhline(0.0, color="#b9bec1", linewidth=0.9, zorder=1)

for year in transitions:
    ax.axvline(pd.Timestamp(f"{year}-01-01"), color="#8c9295", linewidth=1.0,
               linestyle=":", zorder=2)
    ax.text(pd.Timestamp(f"{year}-01-01"), 3.55, f" {year}", fontsize=9,
            color="#6d7376", va="top", ha="left")

for source in SOURCES:
    observed = source == "ERSSTv5"
    ax.plot(pd.to_datetime(index[source]["time"].values), index[source].values,
            color=COLOUR[source], linewidth=1.5 if observed else 0.8,
            alpha=0.95 if observed else 0.55, zorder=4 if observed else 3,
            label=f"{source} (observed)" if observed else source)

ax.plot(pd.to_datetime(published["time"].values), published.values, color=PUBLISHED,
        linewidth=0.9, alpha=0.75, zorder=4, label="published PDO (NCEI, ERSST v6)")
ax.plot(stamps, values, color="#111111", linewidth=2.6, zorder=6,
        label=f"ERSSTv5, {SMOOTH_MONTHS // 12}-yr running mean")

ax.set_ylim(-3.8, 3.8)
ax.set_xlabel("Year", color=INK)
ax.set_ylabel("Standardized PDO index", color=INK)
ax.set_title("PDO index, 1950–2014: three CMIP6 models, the ERSSTv5 reconstruction and the "
             "published series\n"
             f"dotted lines mark the phase transitions of step 6; "
             f"reconstruction against published r = {recon_vs_published:.3f}",
             fontsize=12, color=INK)
ax.legend(loc="lower left", frameon=False, fontsize=9, ncol=3, labelcolor=INK)
ax.grid(axis="y", color="#eceff0", linewidth=0.9)
ax.set_axisbelow(True)
ax.tick_params(colors=INK)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
for side in ("left", "bottom"):
    ax.spines[side].set_color("#c8ccce")
fig.savefig(OUT / "fig2_pdo_timeseries.png", dpi=150, facecolor="white")
plt.close(fig)


# ---------------------------------------------------------------------------
# Numeric answer
# ---------------------------------------------------------------------------
answer = {
    "case_id": "T9",
    "title": "Multi-model assessment of the PDO pattern and its decadal statistics",
    "window": "1950-01 to 2014-12",
    "sources": SOURCES,
    "method": {
        "definition": "leading EOF of monthly North Pacific SST anomalies, 20N-70N / 110E-100W, "
                      "global mean anomaly removed, cells weighted by the square root of area",
        "common_grid": "every source regridded to 1 degree (bilinear) BEFORE the EOF, the "
                       "observations included; patterns on different grids cannot be correlated "
                       "and the explained variance is itself grid-dependent",
        "amplitude": "area-weighted RMS of the pattern regressed on the standardized index, in "
                     "degC per standard deviation; the index's own standard deviation is 1 by "
                     "construction and the unit-norm eigenvector carries no physical magnitude",
        "sign_convention": "positive phase is cool in the central and western North Pacific",
        "phase_transitions": f"sign changes of the {SMOOTH_MONTHS // 12}-year centred running "
                             f"mean, crossings under {MERGE_YEARS} years apart merged",
    },

    "explained_variance_percent": {s: round(variance[s], 1) for s in SOURCES},
    "amplitude_degC_per_sd": {s: round(statistics[s]["amplitude_degC_per_sd"], 3) for s in SOURCES},
    "lag1yr_autocorrelation": {s: round(statistics[s]["lag1yr_autocorrelation"], 3)
                               for s in SOURCES},
    "recon_vs_published_r": round(recon_vs_published, 3),
    "pattern_correlation_vs_obs": {m: round(pattern_correlation[m], 3) for m in MODELS},
    "phase_transition_years": transitions,
    "phase_transition_raw_crossings": crossings,
    "best_pattern_match": best,

    "full_precision": {"statistics": statistics, "pattern_correlation": pattern_correlation,
                       "recon_vs_published_r": recon_vs_published},

    "figures": ["fig1_pdo_patterns.png", "fig2_pdo_timeseries.png"],
    "tools_used": [
        "load_regular_data", "regrid_horizontal", "compute_pdo_index",
        "compute_basic_stat", "temporal_aggregate", "analyze_correlation",
    ],
}

answer["findings"] = [
    f"The ERSSTv5 reconstruction correlates with the published series at r = "
    f"{recon_vs_published:.3f}. That is high, but it is NOT a clean operator validation in the "
    "way case T2 is: the published index comes from ERSST v6 and the gridded field here is v5, "
    "so the residual mixes a product-version difference with any difference in method. It shows "
    "the reconstruction is sound; it cannot separate the two sources of disagreement.",
    f"The models place the mode's pattern well: correlations with observations are "
    + ", ".join(f"{m} {pattern_correlation[m]:.3f}" for m in MODELS)
    + f", with {best} closest.",
    "Reporting the standard deviation of the standardized index would have been vacuous - it is "
    "1 for every source by construction. The amplitude has to come from the regression pattern, "
    "which is the only output in physical units: observations show "
    f"{statistics['ERSSTv5']['amplitude_degC_per_sd']:.3f} degC per standard deviation against "
    + ", ".join(f"{m} {statistics[m]['amplitude_degC_per_sd']:.3f}" for m in MODELS) + ".",
    f"The observed phase transitions are {transitions}. The 1976 one is the documented Pacific "
    "climate shift; a 25-year window such as 1990-2014 would not contain it at all, which is why "
    "the full record is used.",

    "The models' index series are NOT expected to track the observed one in time, and figure 2 "
    "shows that they do not. These are free-running coupled simulations, so their internal "
    "variability has its own phasing and is not synchronised with the real world; only a run "
    "nudged to observations could be. What is comparable is the PATTERN and the STATISTICS - the "
    "variance share, the amplitude, the persistence - which is what steps 2, 3 and 5 measure. A "
    "low time correlation between a model index and the observed one is the expected result and "
    "not a model failure; the published series is the only line in figure 2 that should track "
    "the observed reconstruction, because it is built from nearly the same data.",
]

(OUT / "answer_T9.json").write_text(json.dumps(answer, indent=2) + "\n")

md = [
    "# Reference answer - T9: PDO pattern and decadal statistics",
    "",
    "Window 1950-01 to 2014-12. Leading EOF of North Pacific SST anomalies, 20N-70N / 110E-100W.",
    "Every source, the observations included, regridded to a common 1-degree grid before the EOF.",
    "",
    "## Steps 2-3. Variance, amplitude and persistence",
    "",
    "| Source | EOF1 variance | amplitude (degC per sd) | lag-1-year autocorrelation |",
    "|---|---|---|---|",
]
for source in SOURCES:
    e = statistics[source]
    md.append(f"| {source} | {e['explained_variance_percent']:.1f}% | "
              f"{e['amplitude_degC_per_sd']:.3f} | {e['lag1yr_autocorrelation']:+.3f} |")

md += [
    "",
    "The standard deviation of the standardized index is 1 for every source by construction and",
    "is deliberately not reported as a result. The amplitude comes from the regression pattern,",
    "the only output carrying physical units.",
    "",
    "## Step 4. Reconstruction against the published index",
    "",
    f"Pearson r = **{recon_vs_published:.3f}** over {published.sizes['time']} months.",
    "",
    "The published PDO is built from ERSST **v6** while the gridded field used here is **v5**, so",
    "this is not a pure operator validation: the residual contains the difference between two SST",
    "product versions as well as any difference in method.",
    "",
    "## Step 5. Pattern agreement",
    "",
    "| Model | area-weighted spatial correlation with ERSSTv5 |",
    "|---|---|",
]
for model in MODELS:
    md.append(f"| {model} | {pattern_correlation[model]:.3f} |")

md += [
    "",
    f"Closest to observations: **{best}**.",
    "",
    "## Step 6. Observed phase transitions",
    "",
    f"Sign changes of the {SMOOTH_MONTHS // 12}-year centred running mean: {crossings}.",
    f"After merging crossings less than {MERGE_YEARS} years apart: **{transitions}**.",
    "",
    "## What the comparison shows",
    "",
]
md += [f"- {line}" for line in answer["findings"]]
md += [
    "",
    "## Figures",
    "",
    "- `fig1_pdo_patterns.png` - the regression pattern for the four sources on one shared scale.",
    "- `fig2_pdo_timeseries.png` - the standardized index for four sources plus the published",
    "  series, with the smoothed observed curve and the transition years marked.",
]
(OUT / "answer_T9.md").write_text("\n".join(md) + "\n")

print("\nWrote answer_T9.json, answer_T9.md and two figures")
