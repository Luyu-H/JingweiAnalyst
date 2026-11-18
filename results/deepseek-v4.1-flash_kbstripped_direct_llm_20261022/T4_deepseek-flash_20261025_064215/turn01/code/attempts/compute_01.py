import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
ROOT_HIST = "/data1/lorraine-20419/benchmark_cmip6/historical"
Y0, Y1 = 1971, 2000
MONTHS = [3, 9]
MONTH_NAMES = {3: "March", 9: "September"}
MODELS = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

META = {
    "CESM2": dict(y="nlat", x="nlon", lat="lat", lon="lon", depth_scale=0.01),
    "GFDL-ESM4": dict(y="lat", x="lon", lat="lat", lon="lon", depth_scale=1.0),
    "MPI-ESM1-2-LR": dict(y="j", x="i", lat="latitude", lon="longitude", depth_scale=1.0),
}

# mlotst in CMIP6 is defined with a sigma-theta threshold of 0.03 kg m-3
# relative to a 10 m reference depth.
REF_DEPTH = 10.0
DENSITY_THRESHOLD = 0.03

TARGET_LON = np.arange(-179.5, 180.0, 1.0)
TARGET_LAT = np.arange(-89.5, 90.0, 1.0)


# ---------------------------------------------------------------------------
# File discovery / opening helpers
# ---------------------------------------------------------------------------
def files_for_variable(model, var):
    """Return {variant: [files]} for a model/variable directory."""
    d = os.path.join(ROOT_HIST, model, var)
    fs = sorted(glob.glob(os.path.join(d, f"{var}_*_{model}_historical_*_*.nc")))
    by = {}
    for f in fs:
        parts = os.path.basename(f).split("_")
        if len(parts) >= 6:
            variant = parts[4]
            by.setdefault(variant, []).append(f)
    return by


def open_time_files(files):
    """Open one or several monthly CMIP6 files with cftime time decoding."""
    if len(files) == 1:
        return xr.open_dataset(files[0], use_cftime=True, chunks={"time": 1})
    return xr.open_mfdataset(
        files, combine="by_coords", use_cftime=True, chunks={"time": 1}
    )


def open_fixed_file(model, var):
    """Open an Ofx fixed field, preferring Ofx and r1i1p1f1 if present."""
    d = os.path.join(ROOT_HIST, model, var)
    fs = sorted(glob.glob(os.path.join(d, f"{var}_*_{model}_historical_*.nc")))
    fs_ofx = [f for f in fs if "_Ofx_" in os.path.basename(f)]
    if fs_ofx:
        fs = fs_ofx
    fs_r1 = [f for f in fs if "r1i1p1f1" in os.path.basename(f)]
    if fs_r1:
        fs = fs_r1
    if not fs:
        raise FileNotFoundError(f"No file found for {model}/{var}")
    return xr.open_dataset(fs[0], use_cftime=True)


def load_area(model):
    """Load the model's ocean-cell areas, masking CESM2 with sftof."""
    y, x = META[model]["y"], META[model]["x"]

    ds = open_fixed_file(model, "areacello")
    area = ds["areacello"].squeeze().transpose(y, x).values.astype(np.float64)
    area = np.where(np.isfinite(area), area, 0.0)

    if model == "CESM2":
        ds2 = open_fixed_file(model, "sftof")
        sft = ds2["sftof"].squeeze().transpose(y, x).values.astype(np.float64)
        area = np.where(sft > 0.0, area, 0.0)
    else:
        area = np.where(area > 0.0, area, 0.0)

    return area


# ---------------------------------------------------------------------------
# Density and mixed-layer-depth operator
# ---------------------------------------------------------------------------
def sigma0_from_TS(T, S):
    """
    EOS-80 potential density anomaly sigma0 (kg m-3), using potential
    temperature T (degC) and practical salinity S.
    """
    T = np.asarray(T, dtype=np.float64)
    S = np.asarray(S, dtype=np.float64)

    # CMIP6 salinity is usually stored as practical-salinity numbers (~35).
    # Some files may encode it as 0.035 or as 35000; normalise defensively.
    finite = np.isfinite(S)
    if np.any(finite):
        med = np.nanmedian(S)
        if med < 1.0:
            S = S * 1000.0
        elif med > 100.0:
            S = S / 1000.0

    S = np.where(S < 0.0, np.nan, S)

    rho_w = (
        999.842594
        + 6.793952e-2 * T
        - 9.095290e-3 * T**2
        + 1.001685e-4 * T**3
        - 1.120083e-6 * T**4
        + 6.536332e-9 * T**5
    )
    A = (
        0.824493
        - 4.0899e-3 * T
        + 7.6438e-5 * T**2
        - 8.2467e-7 * T**3
        + 5.3875e-9 * T**4
    )
    B = -5.72466e-3 + 1.0227e-4 * T - 1.6546e-6 * T**2
    C = 4.8314e-4

    rho = rho_w + A * S + B * S**1.5 + C * S**2
    return rho - 1000.0


def compute_mld(T, S, depth, ref_depth=REF_DEPTH, threshold=DENSITY_THRESHOLD):
    """
    Offline mixed layer depth from T/S profiles.

    T, S : shape (lev, ny, nx)
    depth: 1-D positive-down depth axis in metres

    Returns MLD shape (ny, nx).  The MLD is the first depth at which sigma0
    exceeds sigma0 at ref_depth by threshold, linearly interpolated between
    bracketing levels.  If no such crossing exists, the deepest valid level is
    used.
    """
    sig = sigma0_from_TS(T, S)
    L, ny, nx = sig.shape

    depth = np.asarray(depth, dtype=np.float64)
    if depth[0] > depth[-1]:
        depth = depth[::-1]
        sig = sig[::-1]

    valid = np.isfinite(sig)
    depth3 = depth[:, None, None]

    # sigma at the reference depth, using the nearest valid levels above/below.
    above_mask = valid & (depth3 <= ref_depth)
    below_mask = valid & (depth3 >= ref_depth)

    has_above = above_mask.any(axis=0)
    has_below = below_mask.any(axis=0)

    idx_above = L - 1 - np.argmax(above_mask[::-1, :, :], axis=0)
    idx_below = np.argmax(below_mask, axis=0)

    idx_above_safe = np.where(has_above, idx_above, 0)
    idx_below_safe = np.where(has_below, idx_below, 0)

    sig_above = np.take_along_axis(
        sig, idx_above_safe[None, :, :], axis=0
    )[0]
    sig_below = np.take_along_axis(
        sig, idx_below_safe[None, :, :], axis=0
    )[0]

    d_above = depth[idx_above_safe]
    d_below = depth[idx_below_safe]

    denom_ref = d_below - d_above
    frac_ref = np.where(
        denom_ref != 0.0, (ref_depth - d_above) / denom_ref, 0.0
    )
    sig_both = sig_above + frac_ref * (sig_below - sig_above)

    sig_ref = np.where(
        has_above & has_below,
        sig_both,
        np.where(has_above, sig_above, np.where(has_below, sig_below, np.nan)),
    )

    delta = sig - sig_ref[None, :, :]

    # First level at or below ref_depth where the threshold is reached.
    cond = valid & (depth3 >= ref_depth) & (delta >= threshold)
    has_cross = cond.any(axis=0)
    j = np.argmax(cond, axis=0)

    # Default for no crossing: deepest valid level.
    rev_valid = valid[::-1, :, :]
    has_any = valid.any(axis=0)
    last_valid = L - 1 - np.argmax(rev_valid, axis=0)
    bottom_depth = np.where(has_any, depth[last_valid], np.nan)
    mld = np.array(bottom_depth, dtype=np.float64)

    if np.any(has_cross):
        jj = np.where(has_cross, j, 0)
        depth_j = depth[jj]
        delta_j = np.take_along_axis(delta, jj[None, :, :], axis=0)[0]

        jm1 = np.maximum(jj - 1, 0)
        depth_jm1 = depth[jm1]
        delta_jm1 = np.take_along_axis(delta, jm1[None, :, :], axis=0)[0]

        # If the previous level is shallower than the reference depth, the
        # true lower bracket is the reference depth itself (delta = 0).
        use_ref = (depth_jm1 < ref_depth) & (depth_j >= ref_depth)
        depth_jm1 = np.where(use_ref, ref_depth, depth_jm1)
        delta_jm1 = np.where(use_ref, 0.0, delta_jm1)

        denom = delta_j - delta_jm1
        with np.errstate(divide="ignore", invalid="ignore"):
            frac = (threshold - delta_jm1) / denom
        cross = depth_jm1 + frac * (depth_j - depth_jm1)

        bad = ~np.isfinite(cross)
        cross = np.where(bad, depth_j, cross)
        cross = np.where(jj > 0, cross, depth[0])
        mld = np.where(has_cross, cross, mld)

    return mld


# ---------------------------------------------------------------------------
# Regridding for display only
# ---------------------------------------------------------------------------
def regrid_to_regular(lon2d, lat2d, field, tlon, tlat):
    """Nearest-neighbour regridding to a regular lon/lat display grid."""
    try:
        from scipy.interpolate import griddata

        pts = np.column_stack((lon2d.ravel(), lat2d.ravel()))
        vals = field.ravel()
        good = np.isfinite(vals) & np.isfinite(pts[:, 0]) & np.isfinite(pts[:, 1])
        if not np.any(good):
            return np.full((len(tlat), len(tlon)), np.nan)

        lon_grid, lat_grid = np.meshgrid(tlon, tlat)
        out = griddata(
            pts[good], vals[good], (lon_grid, lat_grid), method="nearest"
        )
        return out
    except Exception:
        # Fallback: simple area-of-influence binning.
        out_sum = np.zeros((len(tlat), len(tlon)), dtype=np.float64)
        out_cnt = np.zeros_like(out_sum)
        dlon = tlon[1] - tlon[0]
        dlat = tlat[1] - tlat[0]
        ix = np.floor((lon2d - (tlon[0] - dlon / 2.0)) / dlon).astype(int)
        iy = np.floor((lat2d - (tlat[0] - dlat / 2.0)) / dlat).astype(int)
        good = (
            np.isfinite(field)
            & (ix >= 0)
            & (ix < len(tlon))
            & (iy >= 0)
            & (iy < len(tlat))
        )
        np.add.at(out_sum, (iy[good], ix[good]), field[good])
        np.add.at(out_cnt, (iy[good], ix[good]), 1.0)
        return np.where(out_cnt > 0.0, out_sum / out_cnt, np.nan)


def coarse_disagreement(lon2d, lat2d, area, diff, mask, lon_step=10.0, lat_step=10.0):
    """Return the top coarse boxes by area fraction with |diff| > 20 m."""
    lon_edges = np.arange(-180.0, 180.0 + lon_step, lon_step)
    lat_edges = np.arange(-90.0, 90.0 + lat_step, lat_step)

    ix = np.clip(np.digitize(lon2d, lon_edges) - 1, 0, len(lon_edges) - 2)
    iy = np.clip(np.digitize(lat2d, lat_edges) - 1, 0, len(lat_edges) - 2)

    good = mask & np.isfinite(diff) & np.isfinite(area) & (area > 0.0)
    bad = good & (np.abs(diff) > 20.0)

    area_valid = np.zeros((len(lat_edges) - 1, len(lon_edges) - 1), dtype=np.float64)
    area_bad = np.zeros_like(area_valid)

    np.add.at(area_valid, (iy[good], ix[good]), area[good])
    np.add.at(area_bad, (iy[bad], ix[bad]), area[bad])

    frac = np.where(area_valid > 0.0, 100.0 * area_bad / area_valid, np.nan)
    idx = np.argwhere(np.isfinite(frac) & (area_valid > 0.0))
    if idx.size == 0:
        return pd.DataFrame()

    vals = frac[idx[:, 0], idx[:, 1]]
    order = np.argsort(vals)[::-1][:5]

    rows = []
    for o in order:
        i, j = idx[o]
        rows.append(
            dict(
                lat_min=float(lat_edges[i]),
                lat_max=float(lat_edges[i + 1]),
                lon_min=float(lon_edges[j]),
                lon_max=float(lon_edges[j + 1]),
                area_fraction_gt20=float(frac[i, j]),
                area_valid=float(area_valid[i, j]),
            )
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Main computation
# ---------------------------------------------------------------------------
results = {}
fig1_src = {}
fig2_data = {}
region_stats = {}
model_variant = {}

for model in MODELS:
    by_t = files_for_variable(model, "thetao")
    by_s = files_for_variable(model, "so")
    by_m = files_for_variable(model, "mlotst")

    candidates = sorted(set(by_t) & set(by_s) & set(by_m))
    if not candidates:
        candidates = sorted(set(by_t) | set(by_s) | set(by_m))
    if not candidates:
        raise FileNotFoundError(f"No historical thetao/so/mlotst files for {model}")

    variant = "r1i1p1f1" if "r1i1p1f1" in candidates else candidates[0]
    model_variant[model] = variant

    ds_t = open_time_files(by_t[variant])
    ds_s = open_time_files(by_s[variant])
    ds_ml = open_time_files(by_m[variant])

    y_dim = META[model]["y"]
    x_dim = META[model]["x"]

    ny = ds_t.sizes[y_dim]
    nx = ds_t.sizes[x_dim]

    area = load_area(model)
    total_ocean_area = float(np.nansum(area))

    # 2-D longitude/latitude for display and region diagnostics.
    lon_name = META[model]["lon"]
    lat_name = META[model]["lat"]
    lon_coord = ds_t[lon_name].values
    lat_coord = ds_t[lat_name].values
    if lon_coord.ndim == 1 and lat_coord.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon_coord, lat_coord)
    elif lon_coord.ndim == 2 and lat_coord.ndim == 2:
        lon2d, lat2d = lon_coord, lat_coord
    else:
        lon2d, lat2d = np.meshgrid(np.ravel(lon_coord), np.ravel(lat_coord))

    lon2d = ((lon2d + 180.0) % 360.0) - 180.0

    # Depth axis in metres.
    lev_da = ds_t["lev"]
    lev_units = str(lev_da.attrs.get("units", "")).lower()
    lev_vals = lev_da.values.astype(np.float64)
    if "cm" in lev_units:
        lev_m = lev_vals * 0.01
    elif lev_units.strip() in {"m", "meter", "metre", "meters", "metres"}:
        lev_m = lev_vals * 1.0
    else:
        lev_m = lev_vals * META[model]["depth_scale"]

    for month in MONTHS:
        cond_t = (
            (ds_t.time.dt.year >= Y0)
            & (ds_t.time.dt.year <= Y1)
            & (ds_t.time.dt.month == month)
        )
        cond_s = (
            (ds_s.time.dt.year >= Y0)
            & (ds_s.time.dt.year <= Y1)
            & (ds_s.time.dt.month == month)
        )
        cond_m = (
            (ds_ml.time.dt.year >= Y0)
            & (ds_ml.time.dt.year <= Y1)
            & (ds_ml.time.dt.month == month)
        )

        ds_t_m = ds_t.sel(time=cond_t)
        ds_s_m = ds_s.sel(time=cond_s)
        ds_ml_m = ds_ml.sel(time=cond_m)

        n = ds_t_m.sizes["time"]
        if ds_s_m.sizes["time"] != n:
            ds_s_m = ds_s_m.sel(time=ds_t_m.time)
        if ds_ml_m.sizes["time"] != n:
            ds_ml_m = ds_ml_m.sel(time=ds_t_m.time)

        sum_op = np.zeros((ny, nx), dtype=np.float64)
        sum_ml = np.zeros((ny, nx), dtype=np.float64)
        count = np.zeros((ny, nx), dtype=np.int16)

        for i in range(n):
            T = (
                ds_t_m["thetao"]
                .isel(time=i)
                .transpose("lev", y_dim, x_dim)
                .load()
                .values
            )
            S = (
                ds_s_m["so"]
                .isel(time=i)
                .transpose("lev", y_dim, x_dim)
                .load()
                .values
            )
            ml = (
                ds_ml_m["mlotst"]
                .isel(time=i)
                .transpose(y_dim, x_dim)
                .load()
                .values
            )

            T = np.where(T > 1.0e20, np.nan, T)
            S = np.where(S > 1.0e20, np.nan, S)
            ml = np.where(ml > 1.0e20, np.nan, ml)

            mld = compute_mld(T, S, lev_m)
            valid = np.isfinite(mld) & np.isfinite(ml)

            if np.any(valid):
                sum_op[valid] += mld[valid]
                sum_ml[valid] += ml[valid]
                count[valid] += 1

        valid_all = count == n

        op_clim = np.full((ny, nx), np.nan, dtype=np.float64)
        ml_clim = np.full((ny, nx), np.nan, dtype=np.float64)
        op_clim[valid_all] = sum_op[valid_all] / float(n)
        ml_clim[valid_all] = sum_ml[valid_all] / float(n)

        diff = op_clim - ml_clim
        mask = valid_all & np.isfinite(op_clim) & np.isfinite(ml_clim)

        x = op_clim[mask]
        y = ml_clim[mask]
        a = area[mask]

        if x.size > 0:
            bias = float(np.mean(x - y))
            rmsd = float(np.sqrt(np.mean((x - y) ** 2)))
            corr = float(np.corrcoef(x, y)[0, 1]) if x.size > 1 else np.nan

            area_valid = float(np.sum(a))
            frac_gt20 = (
                float(np.sum(a[np.abs(x - y) > 20.0]) / area_valid * 100.0)
                if area_valid > 0.0
                else np.nan
            )
            frac_valid_area = (
                float(area_valid / total_ocean_area * 100.0)
                if total_ocean_area > 0.0
                else np.nan
            )
        else:
            bias = np.nan
            rmsd = np.nan
            corr = np.nan
            area_valid = 0.0
            frac_gt20 = np.nan
            frac_valid_area = np.nan

        results[(model, month)] = dict(
            model=model,
            month=month,
            month_name=MONTH_NAMES[month],
            bias=bias,
            rmsd=rmsd,
            corr=corr,
            frac_gt20=frac_gt20,
            frac_valid_area=frac_valid_area,
            n_cells=int(x.size),
            area_valid=area_valid,
            total_ocean_area=total_ocean_area,
            diff=diff,
            op_clim=op_clim,
            ml_clim=ml_clim,
            x=x,
            y=y,
            lon2d=lon2d,
            lat2d=lat2d,
            area=area,
            mask=mask,
        )

        fig2_data[(model, month)] = (x, y, bias, rmsd, corr)
        fig1_src[(model, month)] = regrid_to_regular(
            lon2d, lat2d, diff, TARGET_LON, TARGET_LAT
        )
        region_stats[(model, month)] = coarse_disagreement(
            lon2d, lat2d, area, diff, mask
        )

# ---------------------------------------------------------------------------
# Tables / rankings / plotting ranges
# ---------------------------------------------------------------------------
rows = []
for key, r in results.items():
    rows.append(
        {
            "model": r["model"],
            "month": r["month"],
            "month_name": r["month_name"],
            "bias_m": r["bias"],
            "rmsd_m": r["rmsd"],
            "corr": r["corr"],
            "area_pct_absdiff_gt20": r["frac_gt20"],
            "valid_area_pct_of_ocean": r["frac_valid_area"],
            "n_cells": r["n_cells"],
            "area_valid_m2": r["area_valid"],
            "total_ocean_area_m2": r["total_ocean_area"],
        }
    )

results_df = pd.DataFrame(rows).sort_values(["model", "month"]).reset_index(drop=True)

model_rmsd = results_df.groupby("model")["rmsd_m"].mean().sort_values()
ranking = list(model_rmsd.index)
ranking_criterion = "mean of March and September RMSD (lower is better)"

all_fig1 = np.concatenate([fig1_src[k].ravel() for k in fig1_src])
all_fig1 = all_fig1[np.isfinite(all_fig1)]
if all_fig1.size:
    fig1_vmax = float(np.nanpercentile(np.abs(all_fig1), 98.0))
else:
    fig1_vmax = 100.0
if not np.isfinite(fig1_vmax) or fig1_vmax <= 0.0:
    fig1_vmax = 100.0

all_scatter = np.concatenate(
    [np.concatenate([fig2_data[k][0], fig2_data[k][1]]) for k in fig2_data]
)
all_scatter = all_scatter[np.isfinite(all_scatter) & (all_scatter > 0.0)]
if all_scatter.size:
    fig2_lim_min = 10.0 ** np.floor(np.log10(np.nanmin(all_scatter)))
    fig2_lim_max = 10.0 ** np.ceil(np.log10(np.nanmax(all_scatter)))
else:
    fig2_lim_min, fig2_lim_max = 0.1, 10000.0

fig2_lim_min = max(float(fig2_lim_min), 0.1)
fig2_lim_max = float(fig2_lim_max)
if fig2_lim_max <= fig2_lim_min:
    fig2_lim_max = fig2_lim_min * 10.0
fig2_lims = (fig2_lim_min, fig2_lim_max)

_tick_candidates = np.array(
    [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000]
)
fig2_ticks = [
    float(t) for t in _tick_candidates if fig2_lims[0] <= t <= fig2_lims[1]
]
if len(fig2_ticks) < 2:
    fig2_ticks = [fig2_lims[0], fig2_lims[1]]

print(results_df.to_string(index=False))
print("\nRanking (most to least internally consistent):")
for i, model in enumerate(ranking, 1):
    print(f"{i}. {model}  mean RMSD = {model_rmsd.loc[model]:.1f} m")
