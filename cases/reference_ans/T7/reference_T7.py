"""
Reference answer for case T7 - Departure of intermediate-water N:P from the
Redfield ratio.

Run from the repository root with the project's conda environment:

    python cases/reference_ans/T7/reference_T7.py

Writes answer_T7.json, answer_T7.md, fig1_np_scatter.png and fig2_np_departure_maps.png next to this file. Takes roughly fifteen minutes, almost all of it forming the 30-year mean nutrient fields.
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
from src.tools.scientific_calculator import compute_ratio
from src.tools.basic_calculator import (compute_basic_stat, temporal_aggregate, transform_unit,
                                        compute_exceedance_fraction, fit_multivariate_relation)

OUT = Path(__file__).resolve().parent
MODELS = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
SOURCES = ["WOA23"] + MODELS
WINDOW = ["1971-01", "2000-12"]
DEPTHS = [200, 1000]
BASINS = ["Pacific", "Atlantic", "Indian"]
REDFIELD = 16.0
BAND = 2.0
DENITRIFICATION_NP = 12.0
DENSITY = 1.025
VOLUME_DIMS = ["depth", "lat", "lon"]

xr.set_options(keep_attrs=True)

nitrate, phosphate, np_ratio = {}, {}, {}


# ---------------------------------------------------------------------------
# Step 1. Climatological nutrient fields and the cell-by-cell N:P
# ---------------------------------------------------------------------------
nitrate["WOA23"] = load_climatology_data("woa23", "nitrate", time_slice="annual",
                                         depth_range=DEPTHS)
phosphate["WOA23"] = load_climatology_data("woa23", "phosphate", time_slice="annual",
                                           depth_range=DEPTHS)
print(f"[1]   WOA23 {nitrate['WOA23'].sizes['depth']} levels, "
      f"units {nitrate['WOA23'].attrs['units']}, "
      f"period {nitrate['WOA23'].attrs.get('climatology_period')}")

for model in MODELS:
    started = time.time()
    fields = {}
    for variable in ("no3", "po4"):
        monthly = load_regular_data("cmip6_historical", variable, model=model,
                                    time_range=WINDOW, depth_range=DEPTHS, chunk_time=12)
        fields[variable] = temporal_aggregate(monthly, freq="mean", op="mean").compute()
    nitrate[model], phosphate[model] = fields["no3"], fields["po4"]
    print(f"[1]   {model:14s} 30-year mean nutrients {fields['no3'].shape} "
          f"({time.time() - started:.0f} s)")

for source in SOURCES:
    np_ratio[source] = compute_ratio(nitrate[source], phosphate[source])
    print(f"[1]   {source:14s} N:P field, units '{np_ratio[source].attrs.get('units')}', "
          f"valid {int(np_ratio[source].notnull().sum())} cells")


# ---------------------------------------------------------------------------
# Step 2. Mean N:P, both averages
# ---------------------------------------------------------------------------
def volume_mean(field, dims=VOLUME_DIMS):
    return float(compute_basic_stat(field, op="volume_weighted_mean", dim=dims).compute())


mean_of_ratios, ratio_of_means = {}, {}
for source in SOURCES:
    mean_of_ratios[source] = volume_mean(np_ratio[source])
    ratio_of_means[source] = volume_mean(nitrate[source]) / volume_mean(phosphate[source])
    print(f"[2]   {source:14s} mean-of-ratios {mean_of_ratios[source]:6.3f}   "
          f"ratio-of-means {ratio_of_means[source]:6.3f}")


# ---------------------------------------------------------------------------
# Step 3. Stoichiometric regression, in umol/kg so the intercept is comparable
# ---------------------------------------------------------------------------
regression = {}
for source in SOURCES:
    if source == "WOA23":
        n_field, p_field = nitrate[source], phosphate[source]
        density_note = "already per unit mass"
    else:
        # Both nutrients scale by the same factor, so only the intercept moves.
        n_field = transform_unit(nitrate[source], from_unit=nitrate[source].attrs["units"],
                                 to_unit="umol/kg", density=DENSITY)
        p_field = transform_unit(phosphate[source], from_unit=phosphate[source].attrs["units"],
                                 to_unit="umol/kg", density=DENSITY)
        density_note = n_field.attrs.get("density_assumption", f"constant {DENSITY} kg/L")

    fit = fit_multivariate_relation(target=n_field, predictors=[p_field],
                                    method="linear", validation="blocked")
    coefficients = fit["coefficients"]
    slope = float(coefficients[[k for k in coefficients if k != "intercept"][0]])
    regression[source] = {
        "slope_a": slope,
        "intercept_b_umol_per_kg": float(coefficients["intercept"]),
        "r2": float(fit["metrics"]["r2"]),
        "equation": fit["equation"],
        "density_assumption": density_note,
    }
    print(f"[3]   {source:14s} NO3 = {slope:.4f} * PO4 + "
          f"{regression[source]['intercept_b_umol_per_kg']:+.4f}   "
          f"R2 {regression[source]['r2']:.4f}")


# ---------------------------------------------------------------------------
# Step 4. Volume inside the Redfield band
# ---------------------------------------------------------------------------
within_band = {}
for source in SOURCES:
    result = compute_exceedance_fraction(np_ratio[source], threshold=BAND,
                                         comparison="abs_less", reference=REDFIELD,
                                         dim=VOLUME_DIMS, weighting="volume")
    within_band[source] = {
        "percent": float(result["fraction_percent"]),
        "criterion": str(result["criterion"]),
        "total_volume_m3": float(result["total_amount"]),
    }
    print(f"[4]   {source:14s} {within_band[source]['percent']:5.2f}% of volume inside "
          f"{within_band[source]['criterion']}")


# ---------------------------------------------------------------------------
# Step 5. Basin structure and the denitrification signal
# ---------------------------------------------------------------------------
basin_mean, low_np, basin_definition = {}, {}, {}
for source in SOURCES:
    basin_mean[source] = {}
    for basin in BASINS:
        masked = select_basin(np_ratio[source], basin=basin)   # marginal seas folded in
        basin_mean[source][basin] = volume_mean(masked)
        basin_definition[basin] = masked.attrs.get("basin_definition", "")

    low_np[source] = {}
    for label, field in (("global", np_ratio[source]),
                         ("Pacific", select_basin(np_ratio[source], basin="Pacific")),
                         ("Indian", select_basin(np_ratio[source], basin="Indian"))):
        result = compute_exceedance_fraction(field, threshold=DENITRIFICATION_NP,
                                             comparison="less", dim=VOLUME_DIMS,
                                             weighting="volume")
        low_np[source][label] = float(result["fraction_percent"])

    print(f"[5]   {source:14s} " + "  ".join(f"{b} {basin_mean[source][b]:5.2f}" for b in BASINS)
          + f"   |   N:P<12 by volume: global {low_np[source]['global']:5.2f}%  "
            f"Pacific {low_np[source]['Pacific']:5.2f}%  Indian {low_np[source]['Indian']:5.2f}%")


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
INK = "#2f3437"
LAND = "#d9d9d6"
COLOUR = {"WOA23": "#2f3437", "CESM2": "#0072B2",
          "GFDL-ESM4": "#009E73", "MPI-ESM1-2-LR": "#E69F00"}

# --- Figure 1: nitrate against phosphate ----------------------------------
scatter_fields = {}
for source in SOURCES:
    if source == "WOA23":
        scatter_fields[source] = (nitrate[source], phosphate[source])
    else:
        scatter_fields[source] = (
            transform_unit(nitrate[source], from_unit=nitrate[source].attrs["units"],
                           to_unit="umol/kg", density=DENSITY),
            transform_unit(phosphate[source], from_unit=phosphate[source].attrs["units"],
                           to_unit="umol/kg", density=DENSITY))

p_max = max(float(np.nanpercentile(scatter_fields[s][1].values, 99.9)) for s in SOURCES)
n_max = max(float(np.nanpercentile(scatter_fields[s][0].values, 99.9)) for s in SOURCES)
p_edge = np.ceil(p_max * 2) / 2
n_edge = np.ceil(n_max / 5) * 5
print(f"[fig1] axis ranges: PO4 0-{p_edge:.1f}, NO3 0-{n_edge:.0f} umol/kg")

fig, axes = plt.subplots(2, 2, figsize=(11.4, 10.6), constrained_layout=True)
for ax, source in zip(axes.flat, SOURCES):
    n_values = scatter_fields[source][0].values.ravel()
    p_values = scatter_fields[source][1].values.ravel()
    both = np.isfinite(n_values) & np.isfinite(p_values)

    ax.scatter(p_values[both], n_values[both], s=2.0, color=COLOUR[source],
               alpha=0.10, linewidths=0, zorder=2)
    x = np.linspace(0, p_edge, 50)
    ax.plot(x, REDFIELD * x, color="#8c9295", linewidth=1.5, linestyle=":", zorder=4,
            label=f"Redfield {REDFIELD:.0f}:1")
    entry = regression[source]
    ax.plot(x, entry["slope_a"] * x + entry["intercept_b_umol_per_kg"],
            color="#c1543a", linewidth=1.8, linestyle="--", zorder=5, label="fitted")

    ax.set_xlim(0, p_edge)
    ax.set_ylim(0, n_edge)
    label = "WOA23 (observed)" if source == "WOA23" else source
    ax.set_title(f"{label}\n"
                 f"a = {entry['slope_a']:.2f},   b = {entry['intercept_b_umol_per_kg']:+.2f} "
                 f"µmol kg⁻¹,   R² = {entry['r2']:.3f}",
                 fontsize=10.5, color=INK)
    ax.set_xlabel("Phosphate (µmol kg⁻¹)", color=INK, fontsize=9.5)
    ax.set_ylabel("Nitrate (µmol kg⁻¹)", color=INK, fontsize=9.5)
    ax.legend(loc="upper left", frameon=False, fontsize=9, labelcolor=INK)
    ax.grid(color="#eceff0", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(colors=INK, labelsize=9)
    for spine in ax.spines.values():
        spine.set_color("#c8ccce")

fig.suptitle("Nitrate against phosphate, 200–1000 m, every valid cell\n"
             "models 1971–2000, WOA23 all-period climatology; model values converted from "
             f"mol m⁻³ at a constant {DENSITY} kg L⁻¹.\n"
             "Dotted the Redfield 16:1 line through the origin, dashed each source's own fit",
             fontsize=12, color=INK)
fig.savefig(OUT / "fig1_np_scatter.png", dpi=150, facecolor="white")
plt.close(fig)

# --- Figure 2: departure from Redfield, mapped ----------------------------
# The depth average is volume-weighted, so thick deep layers do not count the
# same as thin ones; the maps are regridded only for display.
column = {}
for source in SOURCES:
    depth_mean = compute_basic_stat(np_ratio[source], op="volume_weighted_mean",
                                    dim="depth").compute()
    depth_mean.attrs = dict(np_ratio[source].attrs)
    if source == "WOA23":
        column[source] = depth_mean
    else:
        column[source] = regrid_horizontal(depth_mean, target="woa23", method="conservative")
        print(f"[fig2] {source:14s} departure map regridded for display")

stack = np.concatenate([np.abs(column[s].values[np.isfinite(column[s].values)] - REDFIELD)
                        for s in SOURCES])
vmax = float(np.ceil(np.percentile(stack, 98)))
print(f"[fig2] shared diverging scale +/-{vmax:.0f} "
      f"(98th percentile of |N:P - 16| = {np.percentile(stack, 98):.2f})")

lon = column["WOA23"]["lon"].values
lat = column["WOA23"]["lat"].values

fig, axes = plt.subplots(2, 2, figsize=(14.6, 8.2))
fig.subplots_adjust(left=0.06, right=0.985, bottom=0.145, top=0.855, wspace=0.12, hspace=0.30)
for ax, source in zip(axes.flat, SOURCES):
    ax.set_facecolor(LAND)
    mesh = ax.pcolormesh(lon, lat, column[source].values - REDFIELD, cmap="PuOr_r",
                         norm=TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax), shading="auto")
    label = "WOA23 (observed)" if source == "WOA23" else source
    ax.set_title(f"{label}\nmean N:P {mean_of_ratios[source]:.2f},   "
                 f"{low_np[source]['global']:.1f}% of volume below N:P 12",
                 fontsize=10.5, color=INK)
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
bar = fig.colorbar(mesh, cax=cax, orientation="horizontal", extend="both")
bar.set_label("200–1000 m volume-weighted N:P minus 16 "
              "(negative = poorer in nitrate than Redfield)", color=INK, fontsize=10)
bar.ax.tick_params(colors=INK, labelsize=9)
fig.suptitle("Departure of intermediate-water N:P from the Redfield ratio\n"
             "models 1971–2000, WOA23 all-period climatology; model panels regridded to "
             "1° for display,\nevery reported number computed on the native grid",
             fontsize=12, color=INK, y=0.985)
fig.savefig(OUT / "fig2_np_departure_maps.png", dpi=150, facecolor="white")
plt.close(fig)


# ---------------------------------------------------------------------------
# Numeric answer
# ---------------------------------------------------------------------------
answer = {
    "case_id": "T7",
    "title": "Departure of intermediate-water N:P from the Redfield ratio",
    "sources": SOURCES,
    "depth_range_m": DEPTHS,
    "method": {
        "climatology_periods": "models 1971-2000; WOA23 nitrate and phosphate are the 1965-2022 "
                               "all-period composite, the only one published for the nutrients - "
                               "a genuine mismatch that cannot be removed",
        "order_of_operations": "the climatological nutrient fields are formed first and every "
                               "stoichiometric quantity computed from them",
        "units": "the N:P ratio needs no conversion, since both nutrients of a source share a "
                 "molar basis and the density factor cancels; the regression intercept does, and "
                 f"the model fields were converted at a constant {DENSITY} kg/L",
        "grid": "every reported number on the source's native grid with its own cell volumes; "
                "regridding only for the display maps",
        "marginal_seas": "folded into their parent ocean (select_basin default)",
        "southern_and_arctic_oceans": "separate basins in the WOA23 mask, not included in the Atlantic, "
                                      "Indian or Pacific numbers",
        "redfield_band": within_band["WOA23"]["criterion"],
    },

    "mean_np_ratio_of_cell_ratios": {s: round(mean_of_ratios[s], 2) for s in SOURCES},
    "mean_np_ratio_of_means": {s: round(ratio_of_means[s], 2) for s in SOURCES},

    "regression": {s: {"slope_a": round(regression[s]["slope_a"], 2),
                       "intercept_b_umol_per_kg": round(regression[s]["intercept_b_umol_per_kg"], 2),
                       "r2": round(regression[s]["r2"], 3)} for s in SOURCES},

    "volume_percent_within_redfield_band": {s: round(within_band[s]["percent"], 1) for s in SOURCES},

    "basin_mean_np": {s: {b: round(basin_mean[s][b], 2) for b in BASINS} for s in SOURCES},
    "volume_percent_np_below_12": {s: {k: round(v, 1) for k, v in low_np[s].items()}
                                   for s in SOURCES},
    "basin_definitions": basin_definition,

    "full_precision": {
        "mean_of_ratios": mean_of_ratios,
        "ratio_of_means": ratio_of_means,
        "regression": regression,
        "within_band": within_band,
        "basin_mean_np": basin_mean,
        "volume_percent_np_below_12": low_np,
    },

    "figures": ["fig1_np_scatter.png", "fig2_np_departure_maps.png"],
    "tools_used": [
        "load_regular_data", "load_climatology_data", "temporal_aggregate", "compute_ratio",
        "compute_basic_stat", "compute_exceedance_fraction", "fit_multivariate_relation",
        "select_basin", "transform_unit", "regrid_horizontal",
    ],
}

observed = mean_of_ratios["WOA23"]
answer["findings"] = [
    f"Observed intermediate-water N:P is {observed:.2f}, clearly BELOW the canonical Redfield 16, "
    f"and the observed regression slope is {regression['WOA23']['slope_a']:.2f} with a positive "
    f"intercept of {regression['WOA23']['intercept_b_umol_per_kg']:+.2f} micromoles per kilogram. "
    "The 16:1 ratio is an uptake ratio for organic matter, not the ratio of the dissolved pool in "
    "aged intermediate water, so this is the expected result rather than an error.",
    "The two averages of step 2 differ for every source, which is why both are reported: the "
    "volume-weighted mean of the cell ratios gives each parcel the same say, while the ratio of "
    "the volume-weighted means is pulled towards the nutrient-rich water.",
    "The regression of step 3 weights cells equally and the mean of step 2 weights them by "
    "volume, so the slope and the mean ratio answer different questions and do not coincide.",

    "All three models place the low-N:P water in the right places - the eastern tropical "
    "Pacific, the Arabian Sea and the Bay of Bengal - so the denitrification signal is "
    f"reproduced in pattern. All three overdo it in amount: {low_np['WOA23']['global']:.1f} per "
    "cent of the observed 200-1000 m volume sits below N:P 12, against "
    + ", ".join(f"{m} {low_np[m]['global']:.1f}" for m in MODELS)
    + " per cent in the models.",

    f"CESM2 is the outlier. Its nitrate-phosphate relationship has largely broken down: the "
    f"slope is {regression['CESM2']['slope_a']:.2f} against an observed "
    f"{regression['WOA23']['slope_a']:.2f}, and R-squared falls to "
    f"{regression['CESM2']['r2']:.3f} against {regression['WOA23']['r2']:.3f}. The scatter shows "
    "why: a large branch of its water runs to high phosphate with nitrate collapsing towards "
    f"zero, and {low_np['CESM2']['Pacific']:.1f} per cent of its Pacific intermediate water "
    "falls below N:P 12. This is the same excessive eastern-tropical-Pacific denitrification "
    "that case T6 finds as a deep nitrate deficit in the same model.",

    f"GFDL-ESM4 reproduces the observed mean N:P almost exactly "
    f"({mean_of_ratios['GFDL-ESM4']:.2f} against {mean_of_ratios['WOA23']:.2f}) and keeps the "
    f"tightest nitrate-phosphate relationship of the three (R-squared "
    f"{regression['GFDL-ESM4']['r2']:.3f}). It holds MORE of its volume inside the Redfield band "
    f"than the observations do ({within_band['GFDL-ESM4']['percent']:.1f} against "
    f"{within_band['WOA23']['percent']:.1f} per cent), so it is not merely closer to observations "
    "but closer to Redfield than the real ocean is.",

    f"MPI-ESM1-2-LR has almost no water in the Redfield band at all "
    f"({within_band['MPI-ESM1-2-LR']['percent']:.1f} per cent), despite a mean N:P of "
    f"{mean_of_ratios['MPI-ESM1-2-LR']:.2f}: its distribution is tightly clustered below 14 "
    "rather than spread across the band. A mean close to the observed one can hide a "
    "distribution that is not.",

    "Every source puts the Atlantic highest and the Pacific lowest, so the inter-basin ordering "
    "of N:P is reproduced by all three models even where the absolute values are not.",

    "Part of the spread is model structure rather than circulation. GFDL-ESM4's COBALTv2 uses fixed "
    "N:P ratios for each plankton type and MPI-ESM1-2-LR's HAMOCC6 a prescribed Redfield composition "
    "of organic matter, so in those two models dissolved N:P can depart from the uptake ratio only "
    "through the nitrogen cycle - N2 fixation, denitrification, deposition - and transport; CESM2's "
    "MARBL has variable C:P stoichiometry (Dunne et al. 2020; Mauritsen et al. 2019; Danabasoglu "
    "et al. 2020).",
]
(OUT / "answer_T7.json").write_text(json.dumps(answer, indent=2) + "\n")


md = [
    "# Reference answer - T7: intermediate-water N:P against Redfield",
    "",
    f"Depth range {DEPTHS[0]}-{DEPTHS[1]} m. Models 1971-2000; WOA23 nitrate and phosphate are the",
    "1965-2022 all-period composite, the only climatology published for the nutrients.",
    "All numbers on each source's native grid, weighted by its own cell volumes.",
    "",
    "## Step 2. Mean N:P, two different averages",
    "",
    "| Source | mean of cell ratios | ratio of means |",
    "|---|---|---|",
]
for source in SOURCES:
    md.append(f"| {source} | {mean_of_ratios[source]:.2f} | {ratio_of_means[source]:.2f} |")

md += [
    "",
    "The first gives every parcel of water the same say and is the headline value; the second is",
    "pulled towards the richest water. They are different quantities, not two estimates of one.",
    "",
    "## Step 3. NO3 = a x PO4 + b, unweighted OLS over all valid cells",
    "",
    "| Source | a | b (umol/kg) | R2 |",
    "|---|---|---|---|",
]
for source in SOURCES:
    e = regression[source]
    md.append(f"| {source} | {e['slope_a']:.2f} | {e['intercept_b_umol_per_kg']:+.2f} | "
              f"{e['r2']:.3f} |")

md += [
    "",
    "## Step 4. Volume inside the Redfield band (14 < N:P < 18)",
    "",
    "| Source | volume within band |",
    "|---|---|",
]
for source in SOURCES:
    md.append(f"| {source} | {within_band[source]['percent']:.1f}% |")

md += [
    "",
    "## Step 5. Basin structure and the denitrification signal",
    "",
    "| Source | Pacific | Atlantic | Indian | N:P<12 global | N:P<12 Pacific | N:P<12 Indian |",
    "|---|---|---|---|---|---|---|",
]
for source in SOURCES:
    md.append(f"| {source} | " + " | ".join(f"{basin_mean[source][b]:.2f}" for b in BASINS)
              + f" | {low_np[source]['global']:.1f}% | {low_np[source]['Pacific']:.1f}% "
                f"| {low_np[source]['Indian']:.1f}% |")
md += ["", "Marginal seas are folded into their parent ocean. The Southern and Arctic Oceans are separate "
    "basins in the WOA23 mask and are not part of these numbers.", "",
       "## What the comparison shows", ""]
md += [f"- {line}" for line in answer["findings"]]
md += [
    "",
    "## Figures",
    "",
    "- `fig1_np_scatter.png` - nitrate against phosphate for the four sources on identical axes,",
    "  each with its own fit and the Redfield 16:1 line.",
    "- `fig2_np_departure_maps.png` - 200-1000 m volume-weighted N:P minus 16, one shared",
    "  diverging scale, models regridded for display only.",
]
(OUT / "answer_T7.md").write_text("\n".join(md) + "\n")

print("\nWrote answer_T7.json, answer_T7.md and two figures")
