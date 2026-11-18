import os
import numpy as np
import pandas as pd
import xarray as xr

# ---------------------------------------------------------------------
# 1. Open raw ERSSTv5 and build the Nino 3.4 regional index
# ---------------------------------------------------------------------
ersst_path = "/data1/lorraine-20419/ocean_data/ERSSTv5/sst.mnmean.nc"
cpc_root = "/data1/lorraine-20419/ocean_data/CPC"

ds = xr.open_dataset(ersst_path)
sst = ds["sst"].sel(time=slice("1950-01-01", "2023-12-31"))

lat_vals = ds["lat"].values
lon_vals = ds["lon"].values

# Nino 3.4: 5S-5N, 170W-120W.  ERSSTv5 is 0-360, descending latitude.
if lat_vals[0] > lat_vals[-1]:
    lat_slice = slice(5, -5)
else:
    lat_slice = slice(-5, 5)

if lon_vals.min() < 0:
    lon_start, lon_end = -170, -120
else:
    lon_start, lon_end = 190, 240

if lon_vals[0] > lon_vals[-1]:
    lon_slice = slice(lon_end, lon_start)
else:
    lon_slice = slice(lon_start, lon_end)

sst_reg = sst.sel(lat=lat_slice, lon=lon_slice)


def build_nino34_index(sst_da, clim_start, clim_end, name):
    """Monthly Nino 3.4 index from SST, using a fixed monthly climatology."""
    clim = (
        sst_da.sel(time=slice(clim_start, clim_end))
        .groupby("time.month")
        .mean(dim="time")
    )
    anom = sst_da.groupby("time.month") - clim

    # cos(latitude) area weights, restricted to cells valid in the data.
    w = np.cos(np.deg2rad(sst_da["lat"]))
    valid = np.isfinite(anom)
    num = (anom * w).where(valid).sum(dim=("lat", "lon"), skipna=True)
    den = w.where(valid).sum(dim=("lat", "lon"), skipna=True)
    idx = (num / den).rename(name)
    return idx.compute()


recon_9120_da = build_nino34_index(
    sst_reg, "1991-01-01", "2020-12-31", "RECON_9120"
)
recon_8110_da = build_nino34_index(
    sst_reg, "1981-01-01", "2010-12-31", "RECON_8110"
)

time_index = pd.date_range("1950-01-01", "2023-12-01", freq="MS")

recon_9120_series = pd.Series(
    recon_9120_da.values,
    index=pd.to_datetime(recon_9120_da["time"].values),
    name="RECON_9120",
).reindex(time_index)

recon_8110_series = pd.Series(
    recon_8110_da.values,
    index=pd.to_datetime(recon_8110_da["time"].values),
    name="RECON_8110",
).reindex(time_index)

# ---------------------------------------------------------------------
# 2. Read CPC monthly Nino indices
# ---------------------------------------------------------------------
def read_cpc_nino_monthly(path):
    rows = []
    with open(path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 10:
                continue
            try:
                yr = int(parts[0])
                mon = int(parts[1])
            except ValueError:
                continue
            vals = []
            ok = True
            for x in parts[2:10]:
                try:
                    vals.append(float(x))
                except ValueError:
                    ok = False
                    break
            if not ok:
                continue
            rows.append([yr, mon] + vals)

    cols = [
        "YR", "MON",
        "NINO1+2", "NINO1+2_ANOM",
        "NINO3", "NINO3_ANOM",
        "NINO4", "NINO4_ANOM",
        "NINO3.4", "NINO3.4_ANOM",
    ]
    df = pd.DataFrame(rows, columns=cols)

    for c in cols[2:]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
        # CPC missing-value sentinels are around -999; real anomalies are O(1).
        df.loc[df[c] <= -90, c] = np.nan

    df["time"] = pd.to_datetime(dict(year=df["YR"], month=df["MON"], day=1))
    return df


nino91 = read_cpc_nino_monthly(
    os.path.join(cpc_root, "ersst5.nino.mth.91-20.ascii")
)
nino81 = read_cpc_nino_monthly(
    os.path.join(cpc_root, "ersst5.nino.mth.81-10.ascii")
)


def monthly_nino34_series(df, name):
    s = (
        df.drop_duplicates("time")
        .set_index("time")["NINO3.4_ANOM"]
        .sort_index()
        .reindex(time_index)
    )
    s.name = name
    return s


cpc91_series = monthly_nino34_series(nino91, "CPC_9120")
cpc81_series = monthly_nino34_series(nino81, "CPC_8110")

# ---------------------------------------------------------------------
# 3. Read CPC ONI seasonal file and map seasons to centre months
# ---------------------------------------------------------------------
season_to_month = {
    "DJF": 1, "JFM": 2, "FMA": 3, "MAM": 4,
    "AMJ": 5, "MJJ": 6, "JJA": 7, "JAS": 8,
    "ASO": 9, "SON": 10, "OND": 11, "NDJ": 12,
}


def read_oni(path):
    rows = []
    with open(path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            seas = parts[0]
            if seas not in season_to_month:
                continue
            try:
                yr = int(parts[1])
                total = float(parts[2])
                anom = float(parts[3])
            except ValueError:
                continue
            if anom <= -90:
                anom = np.nan
            month = season_to_month[seas]
            rows.append(
                [seas, yr, total, anom,
                 pd.Timestamp(year=yr, month=month, day=1)]
            )
    return pd.DataFrame(
        rows, columns=["SEAS", "YR", "TOTAL", "ANOM", "time"]
    )


oni_df = read_oni(os.path.join(cpc_root, "oni.ascii.txt"))
oni_series = (
    oni_df.drop_duplicates("time")
    .set_index("time")["ANOM"]
    .sort_index()
    .reindex(time_index)
)
oni_series.name = "ONI"

# ---------------------------------------------------------------------
# 4. Running means and event detection
# ---------------------------------------------------------------------
recon_9120_run = recon_9120_series.rolling(
    window=3, center=True, min_periods=3
).mean()
cpc91_run = cpc91_series.rolling(
    window=3, center=True, min_periods=3
).mean()


def detect_threshold_runs(series, threshold, min_duration=5, direction="high"):
    """Find contiguous runs of at least min_duration months beyond threshold."""
    vals = series.values.astype(float)
    times = pd.to_datetime(series.index)

    if direction == "high":
        cond = np.isfinite(vals) & (vals >= threshold)
    else:
        cond = np.isfinite(vals) & (vals <= threshold)

    events = []
    i = 0
    n = len(vals)
    while i < n:
        if cond[i]:
            j = i
            while j + 1 < n and cond[j + 1]:
                j += 1
            length = j - i + 1
            if length >= min_duration:
                seg = vals[i:j + 1]
                if direction == "high":
                    k = i + int(np.nanargmax(seg))
                else:
                    k = i + int(np.nanargmin(seg))
                events.append(
                    {
                        "start": times[i],
                        "end": times[j],
                        "duration": length,
                        "peak": float(vals[k]),
                        "peak_time": times[k],
                    }
                )
            i = j + 1
        else:
            i += 1
    return events


recon_elnino_events = detect_threshold_runs(
    recon_9120_run, 0.5, min_duration=5, direction="high"
)
recon_lanina_events = detect_threshold_runs(
    recon_9120_run, -0.5, min_duration=5, direction="low"
)

# ONI is already the official 3-month running mean: threshold it directly.
oni_elnino_events = detect_threshold_runs(
    oni_series, 0.5, min_duration=5, direction="high"
)


def top_events(events, n=3):
    return sorted(events, key=lambda e: e["peak"], reverse=True)[:n]


recon_elnino_top3 = top_events(recon_elnino_events, 3)
oni_elnino_top3 = top_events(oni_elnino_events, 3)

n_recon_elnino = len(recon_elnino_events)
n_recon_lanina = len(recon_lanina_events)
n_oni_elnino = len(oni_elnino_events)

recon_top3_table = [
    (e["peak_time"].year, e["peak_time"].month, round(e["peak"], 2))
    for e in recon_elnino_top3
]
oni_top3_table = [
    (e["peak_time"].year, e["peak_time"].month, round(e["peak"], 2))
    for e in oni_elnino_top3
]

# ---------------------------------------------------------------------
# 5. Compare RECON_9120 against CPC 1991-2020 Nino 3.4
# ---------------------------------------------------------------------
pair = pd.concat(
    [
        recon_9120_series.rename("RECON_9120"),
        cpc91_series.rename("CPC_9120"),
    ],
    axis=1,
).dropna()

corr = float(pair["RECON_9120"].corr(pair["CPC_9120"]))
rmse = float(
    np.sqrt(
        np.mean(
            (pair["RECON_9120"] - pair["CPC_9120"]) ** 2
        )
    )
)

diff_recon_cpc = pair["RECON_9120"] - pair["CPC_9120"]
max_dev_idx = diff_recon_cpc.abs().idxmax()
max_dev_signed = float(diff_recon_cpc.loc[max_dev_idx])
max_abs_dev = float(abs(max_dev_signed))
max_dev_year = int(max_dev_idx.year)
max_dev_month = int(max_dev_idx.month)

# ---------------------------------------------------------------------
# 6. Baseline offsets, 1950-01 to 2020-12
# ---------------------------------------------------------------------
offset_recon_all = recon_8110_series - recon_9120_series
offset_pub_all = cpc81_series - cpc91_series

offset_recon_months = offset_recon_all.loc[
    "1950-01-01":"2020-12-01"
].dropna()
offset_pub_months = offset_pub_all.loc[
    "1950-01-01":"2020-12-01"
].dropna()

offset_recon_mean = float(offset_recon_months.mean())
offset_pub_mean = float(offset_pub_months.mean())
offset_difference = float(offset_recon_mean - offset_pub_mean)

off_pair = pd.concat(
    [
        offset_recon_months.rename("recon_offset"),
        offset_pub_months.rename("published_offset"),
    ],
    axis=1,
).dropna()

offset_monthly_diff = (
    off_pair["recon_offset"] - off_pair["published_offset"]
)
offset_monthly_diff_mean = float(offset_monthly_diff.mean())
offset_monthly_diff_std = float(offset_monthly_diff.std())
offset_monthly_diff_maxabs = float(offset_monthly_diff.abs().max())
offset_pair_corr = float(
    off_pair["recon_offset"].corr(off_pair["published_offset"])
)

recon_offset_by_month = offset_recon_months.groupby(
    offset_recon_months.index.month
).mean()
pub_offset_by_month = offset_pub_months.groupby(
    offset_pub_months.index.month
).mean()

offset_recon_3dp = round(offset_recon_mean, 3)
offset_pub_3dp = round(offset_pub_mean, 3)
offset_rounded_difference = round(
    offset_recon_3dp - offset_pub_3dp, 3
)

count_explanation = (
    "RECON_9120 uses a fixed 1991-2020 monthly climatology and is smoothed here "
    "with a centred 3-month running mean before the 5-season persistence rule is "
    "applied. CPC ONI is already the official 3-month running mean and uses CPC's "
    "rolling 30-year base periods, updated every five years, rather than a fixed "
    "1991-2020 baseline. Those baseline and smoothing-history differences change "
    "how many seasons cross +/-0.5 C, so the event counts and peak values need not "
    "match even though both represent Nino 3.4."
)

offset_note = (
    "The reconstructed and published baseline offsets may differ by a "
    "near-constant amount per calendar month. That is an offset between the "
    "1981-2010 and 1991-2020 baselines, not an error in either series."
)

# Report key results (values are stored unrounded above; printing rounds only).
print("RECON_9120 vs CPC 1991-2020 (1950-01 to 2023-12)")
print("  Pearson r =", f"{corr:.6f}")
print("  RMSE      =", f"{rmse:.4f}", "degC")
print(
    "  max |dev| =",
    f"{max_abs_dev:.4f}",
    "degC at",
    f"{max_dev_year}-{max_dev_month:02d}",
    f"(signed {max_dev_signed:.4f} degC)",
)

print("\nRECON_9120 El Nino events:", n_recon_elnino)
print("  top peaks (year, month, degC):", recon_top3_table)
print("RECON_9120 La Nina events:", n_recon_lanina)

print("\nCPC ONI El Nino events (ONI thresholded directly):", n_oni_elnino)
print("  top peaks (year, month, degC):", oni_top3_table)

print("\nBaseline offsets over 1950-01 to 2020-12")
print("  reconstructed RECON_8110 - RECON_9120 =", f"{offset_recon_mean:.3f}", "degC")
print("  published CPC 1981-2010 - CPC 1991-2020 =", f"{offset_pub_mean:.3f}", "degC")
print("  difference (reconstructed - published) =", f"{offset_difference:.3f}", "degC")
print(
    "  monthly offset agreement: mean diff =",
    f"{offset_monthly_diff_mean:.3f}",
    "degC, std =",
    f"{offset_monthly_diff_std:.3f}",
    "degC, max |diff| =",
    f"{offset_monthly_diff_maxabs:.3f}",
    "degC, r =",
    f"{offset_pair_corr:.6f}",
)
print("\n", count_explanation)
print("\n", offset_note)
