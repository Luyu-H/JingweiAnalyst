import numpy as np
import pandas as pd
import xarray as xr

# -----------------------------------------------------------------------------
# Load ERSSTv5 SST and form the Nino 3.4 region mean
# -----------------------------------------------------------------------------
ersst_path = "/data1/lorraine-20419/ocean_data/ERSSTv5/sst.mnmean.nc"
cpc_root = "/data1/lorraine-20419/ocean_data/CPC"

ds = xr.open_dataset(ersst_path)
sst = ds["sst"].sel(time=slice("1950-01-01", "2023-12-31"))

# Nino 3.4: 5S-5N, 170W-120W.  ERSST lon is 0-360, so 170W=190E, 120W=240E.
lat_mask = (sst["lat"] >= -5) & (sst["lat"] <= 5)
lon360 = sst["lon"] % 360
lon_mask = (lon360 >= 190) & (lon360 <= 240)

sst_nino34 = (
    sst.where(lat_mask, drop=True)
       .where(lon_mask, drop=True)
       .mean(dim=("lat", "lon"))
       .compute()
)

def monthly_anomaly(da, clim_start, clim_end):
    clim = da.sel(time=slice(clim_start, clim_end)).groupby("time.month").mean("time")
    return da.groupby("time.month") - clim

RECON_9120 = monthly_anomaly(sst_nino34, "1991-01-01", "2020-12-31").rename("RECON_9120")
RECON_8110 = monthly_anomaly(sst_nino34, "1981-01-01", "2010-12-31").rename("RECON_8110")
RECON_9120 = RECON_9120.compute()
RECON_8110 = RECON_8110.compute()

# -----------------------------------------------------------------------------
# Read CPC monthly Nino indices
# -----------------------------------------------------------------------------
def load_cpc_monthly(path):
    raw = pd.read_csv(path, sep=r"\s+", header=None, comment="#", dtype=str)
    recs = []
    for _, r in raw.iterrows():
        if len(r) < 2:
            continue
        try:
            yr = int(float(r.iloc[0]))
            mo = int(float(r.iloc[1]))
        except Exception:
            continue
        vals = []
        for i in range(2, 10):
            if i < len(r):
                try:
                    vals.append(float(r.iloc[i]))
                except Exception:
                    vals.append(np.nan)
            else:
                vals.append(np.nan)
        recs.append([yr, mo] + vals)

    cols = [
        "YR", "MON", "NINO1+2", "NINO1+2_ANOM", "NINO3", "NINO3_ANOM",
        "NINO4", "NINO4_ANOM", "NINO3.4", "NINO3.4_ANOM"
    ]
    df = pd.DataFrame(recs, columns=cols)
    df["time"] = pd.to_datetime(dict(year=df["YR"].astype(int),
                                     month=df["MON"].astype(int),
                                     day=1))
    return df.sort_values("time").reset_index(drop=True)

cpc9120_df = load_cpc_monthly(cpc_root + "/ersst5.nino.mth.91-20.ascii")
cpc8110_df = load_cpc_monthly(cpc_root + "/ersst5.nino.mth.81-10.ascii")

cpc9120_s = cpc9120_df.set_index("time")["NINO3.4_ANOM"].sort_index()
cpc8110_s = cpc8110_df.set_index("time")["NINO3.4_ANOM"].sort_index()
cpc9120_win = cpc9120_s.loc["1950-01-01":"2023-12-31"]

# -----------------------------------------------------------------------------
# Read official CPC ONI seasonal series and map each 3-month season
# to its centre month.
# -----------------------------------------------------------------------------
season_to_month = {
    "DJF": 1, "JFM": 2, "FMA": 3, "MAM": 4, "AMJ": 5, "MJJ": 6,
    "JJA": 7, "JAS": 8, "ASO": 9, "SON": 10, "OND": 11, "NDJ": 12,
}

def load_oni(path):
    raw = pd.read_csv(path, sep=r"\s+", header=None, comment="#", dtype=str)
    recs = []
    for _, r in raw.iterrows():
        if len(r) < 4:
            continue
        seas = str(r.iloc[0]).strip().upper()
        if seas not in season_to_month:
            continue
        try:
            yr = int(float(r.iloc[1]))
            total = float(r.iloc[2])
            anom = float(r.iloc[3])
        except Exception:
            continue
        recs.append({
            "time": pd.Timestamp(year=yr, month=season_to_month[seas], day=1),
            "SEAS": seas, "YR": yr, "TOTAL": total, "ANOM": anom,
        })
    return pd.DataFrame(recs).sort_values("time").set_index("time")

oni_df = load_oni(cpc_root + "/oni.ascii.txt")
oni_anom = oni_df["ANOM"].loc["1950-01-01":"2023-12-31"]

# -----------------------------------------------------------------------------
# Monthly comparison statistics: RECON_9120 versus CPC 1991-2020 Nino3.4
# -----------------------------------------------------------------------------
recon_9120_s = RECON_9120.to_pandas()
comp = pd.concat(
    [recon_9120_s.rename("recon"), cpc9120_win.rename("cpc")],
    axis=1, join="inner"
).dropna()

diff = comp["recon"] - comp["cpc"]

corr_recon_cpc_9120 = float(np.corrcoef(comp["recon"].values, comp["cpc"].values)[0, 1])
rmse_recon_cpc_9120 = float(np.sqrt(np.mean(diff.values ** 2)))

max_abs_series = diff.abs()
max_abs_dev_recon_cpc_9120 = float(max_abs_series.max())
max_abs_dev_time = max_abs_series.idxmax()
max_abs_dev_year = int(max_abs_dev_time.year)
max_abs_dev_month = int(max_abs_dev_time.month)

corr_recon_cpc_9120_6 = round(corr_recon_cpc_9120, 6)
rmse_recon_cpc_9120_4 = round(rmse_recon_cpc_9120, 4)
max_abs_dev_recon_cpc_9120_4 = round(max_abs_dev_recon_cpc_9120, 4)

# -----------------------------------------------------------------------------
# Event detection: ONI criterion, at least 5 consecutive 3-month values
# at or above +0.5 C (or at or below -0.5 C for La Nina).
# -----------------------------------------------------------------------------
def find_runs(series, threshold, min_len=5, above=True):
    s = series.dropna().sort_index()
    if len(s) == 0:
        return []

    cond = (s >= threshold) if above else (s <= threshold)
    times = s.index
    events = []
    n = len(s)
    i = 0

    while i < n:
        if bool(cond.iloc[i]):
            j = i
            while (
                j + 1 < n
                and bool(cond.iloc[j + 1])
                and (
                    (times[j + 1].year * 12 + times[j + 1].month)
                    - (times[j].year * 12 + times[j].month)
                    == 1
                )
            ):
                j += 1

            length = j - i + 1
            if length >= min_len:
                seg = s.iloc[i:j + 1]
                if above:
                    peak_val = float(seg.max())
                    peak_time = seg.idxmax()
                else:
                    peak_val = float(seg.min())
                    peak_time = seg.idxmin()

                events.append({
                    "start": times[i],
                    "end": times[j],
                    "length": length,
                    "peak_time": peak_time,
                    "peak_value": peak_val,
                })
            i = j + 1
        else:
            i += 1

    return events

def summarize_events(events, n=3):
    top = sorted(events, key=lambda e: (-e["peak_value"], e["peak_time"]))[:n]
    return [{
        "year": int(e["peak_time"].year),
        "month": int(e["peak_time"].month),
        "peak_value": round(float(e["peak_value"]), 2),
        "start": e["start"].strftime("%Y-%m"),
        "end": e["end"].strftime("%Y-%m"),
        "length": e["length"],
    } for e in top]

recon_roll = recon_9120_s.rolling(window=3, center=True).mean()
cpc_roll = cpc9120_s.rolling(window=3, center=True).mean()

events_recon = find_runs(recon_roll, 0.5, min_len=5, above=True)
n_el_nino_recon = len(events_recon)
top3_el_nino_recon = summarize_events(events_recon, 3)

events_recon_ln = find_runs(recon_roll, -0.5, min_len=5, above=False)

events_cpc_oni = find_runs(oni_anom, 0.5, min_len=5, above=True)
n_el_nino_cpc_oni = len(events_cpc_oni)
top3_el_nino_cpc_oni = summarize_events(events_cpc_oni, 3)

if n_el_nino_recon == n_el_nino_cpc_oni:
    count_difference_explanation = (
        f"Both series contain {n_el_nino_recon} El Nino events. "
        "Small differences in event boundaries or peaks can still arise because "
        "RECON_9120 uses a fixed 1991-2020 monthly climatology and a simple area "
        "mean of ERSSTv5, while CPC ONI uses a 30-year centred baseline that is "
        "rolled forward every five years, CPC's own quality-control and missing-data "
        "treatment, and values rounded to 0.1 C."
    )
else:
    count_difference_explanation = (
        f"RECON_9120 gives {n_el_nino_recon} El Nino events, whereas the official "
        f"CPC ONI series gives {n_el_nino_cpc_oni}. The discrepancy is caused by the "
        "different climatology and processing: RECON_9120 uses a fixed 1991-2020 "
        "monthly climatology and a simple area mean of the 2-degree ERSSTv5 grid, "
        "while CPC ONI uses a 30-year centred normal rolled forward every five years, "
        "CPC's own quality-control and missing-data treatment, and rounded 0.1 C "
        "values. Near the +0.5 C threshold these differences can make a run one "
        "season shorter or longer and thus add or remove an event."
    )

# -----------------------------------------------------------------------------
# Baseline offsets: 1981-2010 anomaly minus 1991-2020 anomaly
# -----------------------------------------------------------------------------
recon_diff = (
    RECON_8110.to_pandas() - RECON_9120.to_pandas()
).loc["1950-01-01":"2020-12-31"]

offset_recon_8110_minus_9120 = float(recon_diff.mean())

cpc_offset_df = pd.concat(
    [cpc8110_s.rename("cpc8110"), cpc9120_s.rename("cpc9120")],
    axis=1, join="inner"
).loc["1950-01-01":"2020-12-31"].dropna()

offset_cpc_8110_minus_9120 = float(
    (cpc_offset_df["cpc8110"] - cpc_offset_df["cpc9120"]).mean()
)

offset_recon_8110_minus_9120_3 = round(offset_recon_8110_minus_9120, 3)
offset_cpc_8110_minus_9120_3 = round(offset_cpc_8110_minus_9120, 3)

offset_difference_cpc_minus_recon = float(
    offset_cpc_8110_minus_9120 - offset_recon_8110_minus_9120
)
offset_abs_difference = float(abs(offset_difference_cpc_minus_recon))

# Convenience printout for checking the requested numbers.
print("RECON_9120 versus CPC 1991-2020 Nino3.4:")
print(f"  Pearson r = {corr_recon_cpc_9120:.6f}")
print(f"  RMSE = {rmse_recon_cpc_9120:.4f} C")
print(f"  max |dev| = {max_abs_dev_recon_cpc_9120:.4f} C "
      f"at {max_abs_dev_year}-{max_abs_dev_month:02d}")
print(f"El Nino count RECON_9120: {n_el_nino_recon}")
print(f"  top3 = {top3_el_nino_recon}")
print(f"El Nino count CPC ONI: {n_el_nino_cpc_oni}")
print(f"  top3 = {top3_el_nino_cpc_oni}")
print(f"Baseline offset RECON (8110-9120): "
      f"{offset_recon_8110_minus_9120:.3f} C")
print(f"Baseline offset CPC   (8110-9120): "
      f"{offset_cpc_8110_minus_9120:.3f} C")
print(f"CPC minus RECON offset difference: "
      f"{offset_difference_cpc_minus_recon:.3f} C")
print(count_difference_explanation)
