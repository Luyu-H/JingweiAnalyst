"""
Reference answer for case T2 - Validating a reconstructed Nino 3.4 index
against the CPC published series.

Run from the repository root with the project's conda environment:

    source /home/lorraine/miniconda3/bin/activate base && conda activate jingwei_QA
    PYTHONPATH=. python cases/reference_ans/T2/reference_T2.py

Writes answer_T2.json, answer_T2.md, fig1_series_overlay.png and
fig2_scatter.png next to this file.

Every physical quantity is produced by a tool from src/tools; nothing here
re-implements an index, an anomaly, an error metric or an event rule by hand.
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
from src.tools.basic_calculator import (temporal_aggregate, detect_threshold_events,
                                        compute_difference_metrics, analyze_correlation)

OUT = Path(__file__).resolve().parent
WINDOW = ["1950-01", "2023-12"]          # the record the case is scored on
OFFSET_WINDOW = slice("1950-01", "2020-12")   # where the CPC 1981-2010 series ends
BASE_9120 = ["1991-01", "2020-12"]
BASE_8110 = ["1981-01", "2010-12"]

xr.set_options(keep_attrs=True)


# ---------------------------------------------------------------------------
# Step 1. Reconstruct the index from the raw ERSSTv5 field
# ---------------------------------------------------------------------------
sst = load_regular_data("ersst", "sst", time_range=WINDOW)
assert sst.sizes["time"] == 888, f"{sst.sizes['time']} months, expected 888"

# The Nino 3.4 box is cut out and area-weighted on ERSSTv5's own 2-degree grid;
# base_period forces the field's own monthly climatology over that window.
recon_9120 = compute_enso_nino_index(sst, nino_region="nino3.4", base_period=BASE_9120).compute()
recon_8110 = compute_enso_nino_index(sst, nino_region="nino3.4", base_period=BASE_8110).compute()
print(f"[1]   RECON_9120 baseline: {recon_9120.attrs['climatology_baseline']}")
print(f"[1]   region: {recon_9120.attrs['region_bounds']}, weighting: {recon_9120.attrs['weighting']}")


# ---------------------------------------------------------------------------
# Step 2. Agreement with the CPC 1991-2020 series
# ---------------------------------------------------------------------------
cpc_9120 = load_regular_data("cpc_indices", "NINO3.4", product="ersst5_91-20",
                             field="anomaly", time_range=WINDOW)
assert cpc_9120.sizes["time"] == recon_9120.sizes["time"]

diff = compute_difference_metrics(recon_9120, cpc_9120, dim="time", weighting="none")
corr = analyze_correlation(recon_9120, cpc_9120, dim="time", method="pearson")

pearson_r = float(corr["correlation"])
rmse = float(diff["rmse"])
max_abs = float(diff["max_abs_error"])
max_abs_at = diff["max_abs_error_at"]["time"]
print(f"[2]   r = {pearson_r:.8f}   RMSE = {rmse:.5f} degC   "
      f"max|dev| = {max_abs:.5f} degC at {max_abs_at}")


# ---------------------------------------------------------------------------
# Step 3-4. ONI criterion, on the reconstruction and on the official ONI
# ---------------------------------------------------------------------------
def el_nino_and_la_nina(series, label):
    """The operational ONI rule: >=5 consecutive overlapping seasons past +/-0.5."""
    warm = detect_threshold_events(series, threshold=0.5, direction="above", min_duration=5)
    cold = detect_threshold_events(series, threshold=-0.5, direction="below", min_duration=5)
    print(f"[3-4] {label:12s} El Nino {len(warm):2d} events, La Nina {len(cold):2d} events")
    return warm, cold


def top_three(events, by_abs=False):
    """The three strongest events, as plain records for the answer file."""
    key = events["peak_value"].abs() if by_abs else events["peak_value"]
    ordered = events.loc[key.sort_values(ascending=False).index].head(3)
    return [
        {"peak_time": pd.Timestamp(row.peak_time).strftime("%Y-%m"),
         "peak_value_degC": round(float(row.peak_value), 2),
         "start": pd.Timestamp(row.start).strftime("%Y-%m"),
         "end": pd.Timestamp(row.end).strftime("%Y-%m"),
         "duration_seasons": int(row.duration)}
        for row in ordered.itertuples()
    ]


# One value per overlapping 3-month season, which is what the ONI rule counts.
running_9120 = temporal_aggregate(recon_9120, freq="rolling", op="mean", window=3, center=True)
recon_warm, recon_cold = el_nino_and_la_nina(running_9120, "RECON_9120")

# The ONI is defined on Nino3.4 alone, so it takes no region.
oni = load_regular_data("cpc_indices", product="oni", field="anomaly", time_range=WINDOW)
oni_warm, oni_cold = el_nino_and_la_nina(oni, "CPC ONI")


# ---------------------------------------------------------------------------
# Step 5-6. Baseline offset, reconstructed and published
# ---------------------------------------------------------------------------
cpc_8110 = load_regular_data("cpc_indices", "NINO3.4", product="ersst5_81-10", field="anomaly")

offset_recon = compute_difference_metrics(recon_8110.sel(time=OFFSET_WINDOW),
                                          recon_9120.sel(time=OFFSET_WINDOW),
                                          dim="time", weighting="none")
offset_cpc = compute_difference_metrics(cpc_8110.sel(time=OFFSET_WINDOW),
                                        cpc_9120.sel(time=OFFSET_WINDOW),
                                        dim="time", weighting="none")
offset_recon_value = float(offset_recon["bias"])
offset_cpc_value = float(offset_cpc["bias"])
n_offset = int(offset_recon["n_valid"])
assert n_offset == int(offset_cpc["n_valid"]), "the two offsets must cover the same months"
print(f"[5-6] offset 8110-9120: reconstructed {offset_recon_value:+.5f} degC, "
      f"published {offset_cpc_value:+.5f} degC, over {n_offset} months "
      f"(agree to {abs(offset_recon_value - offset_cpc_value):.5f} degC)")


# ---------------------------------------------------------------------------
# Figure 1. Series overlay with event shading, plus the residual
# ---------------------------------------------------------------------------
INK = "#2f3437"
WARM, COOL = "#c1543a", "#3b6ea8"
RECON_C, CPC_C = "#1f4e79", "#d98b3a"

time = pd.to_datetime(recon_9120["time"].values)
recon_v = recon_9120.values
cpc_v = cpc_9120.values
run_v = running_9120.values
residual = recon_v - cpc_v

fig, (ax_top, ax_bot) = plt.subplots(
    2, 1, figsize=(15.0, 7.6), sharex=True,
    gridspec_kw={"height_ratios": [3.0, 1.0], "hspace": 0.12})

# Shade each detected event between the threshold and the running-mean curve.
for events, threshold, colour in ((recon_warm, 0.5, WARM), (recon_cold, -0.5, COOL)):
    for row in events.itertuples():
        inside = (time >= pd.Timestamp(row.start)) & (time <= pd.Timestamp(row.end))
        ax_top.fill_between(time[inside], threshold, run_v[inside], color=colour,
                            alpha=0.55, linewidth=0, zorder=2)

ax_top.axhline(0.0, color="#b9bec1", linewidth=0.8, zorder=1)
for threshold in (0.5, -0.5):
    ax_top.axhline(threshold, color="#8c9295", linewidth=0.9, linestyle="--", zorder=1)
# Against the right margin with a white pad: the shading reaches the thresholds
# almost everywhere, so a label sitting on the curve is unreadable.
label_pad = dict(facecolor="white", edgecolor="none", boxstyle="square,pad=0.15")
ax_top.text(time[-1], 0.5, " +0.5 °C", fontsize=9, color="#6d7376",
            va="center", ha="left", bbox=label_pad, zorder=6, clip_on=False)
ax_top.text(time[-1], -0.5, " −0.5 °C", fontsize=9, color="#6d7376",
            va="center", ha="left", bbox=label_pad, zorder=6, clip_on=False)

ax_top.plot(time, cpc_v, color=CPC_C, linewidth=1.9, zorder=3,
            label="CPC 1991–2020 monthly anomaly")
ax_top.plot(time, recon_v, color=RECON_C, linewidth=0.9, zorder=4,
            label="Reconstruction from ERSSTv5 (1991–2020 baseline)")
ax_top.plot(time, run_v, color="#111111", linewidth=1.3, zorder=5,
            label="3-month running mean of the reconstruction (ONI rule)")

ax_top.set_ylabel("Nino 3.4 anomaly (°C)", color=INK)
ax_top.set_ylim(-3.0, 3.2)
ax_top.legend(loc="upper left", frameon=False, fontsize=9.5, ncol=2, labelcolor=INK)
ax_top.set_title(
    "Nino 3.4 index rebuilt from ERSSTv5 against the CPC published series, 1950–2023\n"
    f"r = {pearson_r:.6f},   RMSE = {rmse:.4f} °C   |   shading: "
    f"{len(recon_warm)} El Nino and {len(recon_cold)} La Nina events under the ONI criterion",
    fontsize=12, color=INK)

ax_bot.axhline(0.0, color="#b9bec1", linewidth=0.8)
ax_bot.plot(time, residual, color="#5c5f61", linewidth=0.8)
worst = int(np.nanargmax(np.abs(residual)))
ax_bot.plot(time[worst], residual[worst], "o", color=WARM, markersize=5.5, zorder=3)
ax_bot.annotate(f"largest |deviation| {max_abs:.4f} °C, {max_abs_at}",
                xy=(time[worst], residual[worst]), xytext=(12, 10),
                textcoords="offset points", fontsize=9, color=WARM)
ax_bot.set_ylabel("Reconstruction\n− CPC (°C)", color=INK, fontsize=10)
ax_bot.set_xlabel("Year", color=INK)
limit = float(np.nanmax(np.abs(residual))) * 1.6
ax_bot.set_ylim(-limit, limit)

for ax in (ax_top, ax_bot):
    ax.grid(axis="y", color="#eceff0", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(colors=INK)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#c8ccce")

fig.savefig(OUT / "fig1_series_overlay.png", dpi=160, facecolor="white", bbox_inches="tight")
plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2. Scatter of the reconstruction against CPC, on identical axes
# ---------------------------------------------------------------------------
span = float(np.nanmax(np.abs(np.concatenate([recon_v, cpc_v]))))
edge = np.ceil(span * 2) / 2 + 0.25
ticks = np.arange(-np.floor(edge), np.floor(edge) + 0.5, 1.0)

fig, ax = plt.subplots(figsize=(6.6, 6.6), constrained_layout=True)
ax.plot([-edge, edge], [-edge, edge], color="#9aa0a3", linewidth=1.2, linestyle="--",
        zorder=2, label="1:1")
ax.scatter(cpc_v, recon_v, s=11, color=RECON_C, alpha=0.45, linewidths=0, zorder=3)

ax.set_xlim(-edge, edge)
ax.set_ylim(-edge, edge)
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_aspect("equal", adjustable="box")
ax.set_xlabel("CPC 1991–2020 Nino 3.4 anomaly (°C)", color=INK)
ax.set_ylabel("Reconstruction from ERSSTv5 (°C)", color=INK)
ax.set_title(f"Monthly values, 1950–2023 ({len(cpc_v)} months)", fontsize=12, color=INK)
ax.legend(loc="upper left", frameon=False, fontsize=10, labelcolor=INK)
ax.text(0.97, 0.05,
        f"r = {pearson_r:.6f}\nRMSE = {rmse:.4f} °C\nmax |dev| = {max_abs:.4f} °C",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=10.5, color=INK)
ax.grid(color="#eceff0", linewidth=0.8)
ax.set_axisbelow(True)
ax.tick_params(colors=INK)
for spine in ax.spines.values():
    spine.set_color("#c8ccce")

fig.savefig(OUT / "fig2_scatter.png", dpi=160, facecolor="white")
plt.close(fig)


# ---------------------------------------------------------------------------
# Numeric answer
# ---------------------------------------------------------------------------
answer = {
    "case_id": "T2",
    "title": "Validating a reconstructed Nino 3.4 index against the CPC published series",
    "window": "1950-01 to 2023-12",
    "reconstruction": {
        "source": "ERSSTv5 monthly SST on its native 2-degree grid",
        "region": recon_9120.attrs["region_bounds"],
        "weighting": recon_9120.attrs["weighting"],
        "baseline": recon_9120.attrs["climatology_baseline"],
    },

    "pearson_r_recon_vs_cpc": round(pearson_r, 6),
    "rmse_recon_vs_cpc_degC": round(rmse, 4),
    "max_abs_deviation_degC": round(max_abs, 4),
    "max_abs_deviation_at": max_abs_at,

    "n_el_nino_recon": int(len(recon_warm)),
    "top3_el_nino_recon": top_three(recon_warm),
    "n_la_nina_recon": int(len(recon_cold)),

    "n_el_nino_cpc_oni": int(len(oni_warm)),
    "top3_el_nino_cpc_oni": top_three(oni_warm),
    "n_la_nina_cpc_oni": int(len(oni_cold)),

    "baseline_offset_window": "1950-01 to 2020-12",
    "baseline_offset_n_months": n_offset,
    "baseline_offset_reconstructed_degC": round(offset_recon_value, 3),
    "baseline_offset_published_degC": round(offset_cpc_value, 3),
    "baseline_offset_agreement_degC": round(abs(offset_recon_value - offset_cpc_value), 4),

    "full_precision": {
        "pearson_r": pearson_r,
        "rmse_degC": rmse,
        "max_abs_deviation_degC": max_abs,
        "bias_recon_minus_cpc_degC": float(diff["bias"]),
        "std_ratio": float(diff["std_ratio"]),
        "baseline_offset_reconstructed_degC": offset_recon_value,
        "baseline_offset_published_degC": offset_cpc_value,
    },

    "figures": ["fig1_series_overlay.png", "fig2_scatter.png"],
    "method_notes": [
        "The index is the area-weighted SST anomaly over 5S-5N, 170W-120W on ERSSTv5's native "
        "2-degree grid with its own cell areas; the field is never regridded.",
        "The anomaly removes a monthly climatology built from the same field over the stated "
        "baseline window, so each calendar month is referenced to its own 30-year mean.",
        "The ONI criterion is applied to the 3-month centred running mean: a run of at least 5 "
        "consecutive overlapping seasons at or beyond +/-0.5 degC.",
        "The two baseline offsets are both taken over 1950-01 to 2020-12, where the CPC 1981-2010 "
        "series ends, so they cover the identical months and can be checked against each other.",
        "The El Nino counts from the reconstruction and from the CPC ONI differ by design: the ONI "
        "re-centres its baseline every five years and so carries no warming trend, while the "
        "reconstruction is referenced to one fixed 1991-2020 climatology.",
    ],
    "tools_used": [
        "load_regular_data", "compute_enso_nino_index", "temporal_aggregate",
        "detect_threshold_events", "compute_difference_metrics", "analyze_correlation",
    ],
}

(OUT / "answer_T2.json").write_text(json.dumps(answer, indent=2) + "\n")


def event_rows(records):
    return [f"| {r['peak_time']} | {r['peak_value_degC']:.2f} | {r['start']} – {r['end']} | "
            f"{r['duration_seasons']} |" for r in records]


md = [
    "# Reference answer - T2: Nino 3.4 index reconstruction validation",
    "",
    "Window 1950-01 to 2023-12, 888 months. Index rebuilt from ERSSTv5 monthly SST, area-weighted",
    "over 5S-5N / 170W-120W on the native 2-degree grid, anomalies against a 1991-2020 monthly",
    "climatology from the same field.",
    "",
    "## Step 2. Agreement with the CPC 1991-2020 series",
    "",
    "| Quantity | Value |",
    "|---|---|",
    f"| Pearson r | {pearson_r:.6f} |",
    f"| RMSE | {rmse:.4f} degC |",
    f"| Largest absolute deviation | {max_abs:.4f} degC, at {max_abs_at} |",
    f"| Mean bias (reconstruction - CPC) | {float(diff['bias']):+.4f} degC |",
    "",
    "CPC publishes the series to 2 decimal places, and the reconstruction sits well inside that",
    "last digit, so the correlation is quoted to 6 places and the errors to 4; at the 4 and 3",
    "places of the original draft these would read 1.0000 and 0.007 and say nothing.",
    "",
    "## Step 3. El Nino events in the reconstruction, ONI criterion",
    "",
    f"**{len(recon_warm)} El Nino events** (and {len(recon_cold)} La Nina events) over 1950-2023.",
    "",
    "| Peak | Peak value (degC) | Event span | Seasons |",
    "|---|---|---|---|",
]
md += event_rows(answer["top3_el_nino_recon"])
md += [
    "",
    "## Step 4. The same criterion on the CPC ONI",
    "",
    f"**{len(oni_warm)} El Nino events** (and {len(oni_cold)} La Nina events) over the same window.",
    "",
    "| Peak | Peak value (degC) | Event span | Seasons |",
    "|---|---|---|---|",
]
md += event_rows(answer["top3_el_nino_cpc_oni"])
md += [
    "",
    f"The counts differ ({len(recon_warm)} against {len(oni_warm)}) because the ONI re-centres its",
    "baseline every five years and therefore carries no long-term warming trend, whereas the",
    "reconstruction is referenced to a single fixed 1991-2020 climatology. Against that fixed",
    "baseline the 1950s-1980s read cold, and the early El Nino events fall short of +0.5 degC.",
    "The three strongest events are the same three in both, at the same peak months.",
    "",
    "## Steps 5-6. Baseline offset, 1981-2010 relative to 1991-2020",
    "",
    f"Both means taken over 1950-01 to 2020-12 ({n_offset} months), where the CPC 1981-2010 series ends.",
    "",
    "| Source | Mean offset (degC) |",
    "|---|---|",
    f"| Reconstructed, RECON_8110 - RECON_9120 | {offset_recon_value:+.3f} |",
    f"| Published, CPC 81-10 - CPC 91-20 | {offset_cpc_value:+.3f} |",
    "",
    f"The two agree to {abs(offset_recon_value - offset_cpc_value):.4f} degC. The reconstruction",
    "therefore reproduces CPC's baseline arithmetic, not merely the variability of the index.",
    "",
    "## Figures",
    "",
    "- `fig1_series_overlay.png` - the two monthly series overlaid with the 3-month running mean,",
    "  warm/cool shading over each detected event, and the reconstruction-minus-CPC residual below.",
    "- `fig2_scatter.png` - reconstruction against CPC on identical square axes with a 1:1 line.",
]
(OUT / "answer_T2.md").write_text("\n".join(md) + "\n")

print("\nWrote answer_T2.json, answer_T2.md, fig1_series_overlay.png, fig2_scatter.png")
