import numpy as np
import pandas as pd
import xarray as xr

# ----------------------------------------------------------------------
# 1. Rebuild Niño 3.4 from the raw ERSSTv5 gridded SST field
# ----------------------------------------------------------------------
ersst_path = "/data1/lorraine-20419/ocean_data/ERSSTv5/sst.mnmean.nc"
ds = xr.open_dataset(ersst_path)
sst = ds["sst"].sel(time=slice("1950-01-01", "2023-12-31"))

# ERSSTv5 is on 0-360 longitude.  170W-120W = 190E-240E.
lat_n34 = sst["lat"].where((sst["lat"] >= -5) & (sst["lat"] <= 5), drop=True)
lon_n34 = sst["lon"].where((sst["lon"] >= 190) & (sst["lon"] <= 240), drop=True)
sst_n34 = sst.sel(lat=lat_n34, lon=lon_n34)

# Area-weighted regional mean over the Niño 3.4 box.
w = np.cos(np.deg2rad(sst_n34["lat"]))
n34_monthly = sst_n34.weighted(w).mean(dim=("lat", "lon")).load()

def monthly_anomaly(da, base_start, base_end):
    base = da.sel(time=slice(base_start, base_end))
    clim = base.groupby("time.month").mean("time")
    return da.groupby("time.month") - clim

recon_9120_da = monthly_anomaly(n34_monthly, "1991-01-01", "2020-12-31")
recon_8110_da = monthly_anomaly(n34_monthly, "1981-01-01", "2010-12-31")

def to_period_series(da, name):
    s = da.to_series()
    try:
        idx = pd.PeriodIndex(pd.to_datetime(s.index), freq="M")
    except Exception:
        idx = pd.PeriodIndex(
            year=[t.year for t in s.index],
            month=[t.month for t in s.index],
            freq="M",
        )
    s.index = idx
    s.name = name
    return s.sort_index()

RECON_9120 = to_period_series(recon_9120_da, "RECON_9120")
RECON_8110 = to_period_series(recon_8110_da, "RECON_8110")

# ----------------------------------------------------------------------
# 2. Read CPC monthly Niño indices and ONI
# ----------------------------------------------------------------------
def read_cpc_nino_monthly(path):
    rows = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if parts[0].upper() in ("YR", "YEAR"):
                continue
            if len(parts) < 10:
                continue
            try:
                vals = [float(parts[i]) for i in range(10)]
            except ValueError:
                continue
            rows.append(vals)
    cols = [
        "YR", "MON",
        "NINO1+2", "ANOM1+2",
        "NINO3", "ANOM3",
        "NINO4", "ANOM4",
        "NINO3.4", "ANOM3.4",
    ]
    df = pd.DataFrame(rows, columns=cols)
    df["time"] = pd.PeriodIndex(
        year=df["YR"].astype(int),
        month=df["MON"].astype(int),
        freq="M",
    )
    return df.set_index("time").sort_index()

cpc91 = read_cpc_nino_monthly(
    "/data1/lorraine-20419/ocean_data/CPC/ersst5.nino.mth.91-20.ascii"
)
cpc81 = read_cpc_nino_monthly(
    "/data1/lorraine-20419/ocean_data/CPC/ersst5.nino.mth.81-10.ascii"
)
cpc_9120_n34 = cpc91["ANOM3.4"].loc["1950-01":"2023-12"]
cpc_8110_n34 = cpc81["ANOM3.4"].loc["1950-01":"2020-12"]

season_center_month = {
    "DJF": 1, "JFM": 2, "FMA": 3, "MAM": 4,
    "AMJ": 5, "MJJ": 6, "JJA": 7, "JAS": 8,
    "ASO": 9, "SON": 10, "OND": 11, "NDJ": 12,
}

def read_oni(path):
    rows = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if parts[0].upper() == "SEAS":
                continue
            if len(parts) < 4:
                continue
            seas = parts[0].upper()
            if seas not in season_center_month:
                continue
            try:
                yr = int(float(parts[1]))
                total = float(parts[2])
                anom = float(parts[3])
            except ValueError:
                continue
            rows.append((seas, yr, total, anom))
    df = pd.DataFrame(rows, columns=["SEAS", "YR", "TOTAL", "ANOM"])
    df["month"] = df["SEAS"].map(season_center_month).astype(int)
    df["time"] = pd.PeriodIndex(
        year=df["YR"].astype(int),
        month=df["month"],
        freq="M",
    )
    return df.set_index("time").sort_index()

oni = read_oni("/data1/lorraine-20419/ocean_data/CPC/oni.ascii.txt")
oni_n34 = oni["ANOM"].loc["1950-01":"2023-12"]

# ----------------------------------------------------------------------
# 3. Compare reconstruction against CPC 1991-2020 monthly anomalies
# ----------------------------------------------------------------------
comp = pd.concat(
    [
        RECON_9120.rename("RECON_9120"),
        cpc_9120_n34.rename("CPC_9120"),
    ],
    axis=1,
    join="inner",
).dropna().loc["1950-01":"2023-12"]

diff_recon_cpc = comp["RECON_9120"] - comp["CPC_9120"]
corr_9120_cpc = float(comp["RECON_9120"].corr(comp["CPC_9120"]))
rmse_9120_cpc = float(np.sqrt(np.mean(diff_recon_cpc**2)))
max_abs_dev_time = diff_recon_cpc.abs().idxmax()
max_abs_dev_9120_cpc = float(diff_recon_cpc.abs().max())
max_abs_dev_signed = float(diff_recon_cpc.loc[max_abs_dev_time])
max_abs_dev_year = int(max_abs_dev_time.year)
max_abs_dev_month = int(max_abs_dev_time.month)

corr_9120_cpc_6 = round(corr_9120_cpc, 6)
rmse_9120_cpc_4 = round(rmse_9120_cpc, 4)
max_abs_dev_9120_cpc_4 = round(max_abs_dev_9120_cpc, 4)

# ----------------------------------------------------------------------
# 4. Event detection under the operational ONI threshold criterion
# ----------------------------------------------------------------------
def find_threshold_events(series, threshold, warm=True, min_len=5):
    """Return contiguous runs of length >= min_len with values >= threshold
    (warm=True) or <= threshold (warm=False).  Each run is one event; the
    event peak is the max (min for La Niña) 3-month running-mean value."""
    s = series.dropna().sort_index()
    mask = (s >= threshold) if warm else (s <= threshold)

    events = []
    run = []

    def add_run(run_periods):
        if len(run_periods) < min_len:
            return
        seg = s.loc[run_periods]
        if warm:
            peak_time = seg.idxmax()
            peak_value = float(seg.max())
        else:
            peak_time = seg.idxmin()
            peak_value = float(seg.min())
        events.append(
            {
                "start": run_periods[0],
                "end": run_periods[-1],
                "length": len(run_periods),
                "peak_time": peak_time,
                "peak_value": peak_value,
            }
        )

    for t, flag in mask.items():
        if bool(flag):
            run.append(t)
        else:
            add_run(run)
            run = []
    add_run(run)
    return events

# 3-month centred running mean of the reconstruction.
recon_9120_run3 = RECON_9120.rolling(window=3, center=True, min_periods=3).mean()

elnino_events_recon = find_threshold_events(
    recon_9120_run3, 0.5, warm=True, min_len=5
)
lanina_events_recon = find_threshold_events(
    recon_9120_run3, -0.5, warm=False, min_len=5
)
n_elnino_recon = len(elnino_events_recon)
top3_elnino_recon = sorted(
    elnino_events_recon, key=lambda e: e["peak_value"], reverse=True
)[:3]
top3_elnino_recon_table = pd.DataFrame(
    [
        {
            "year": int(e["peak_time"].year),
            "month": int(e["peak_time"].month),
            "peak_value": round(e["peak_value"], 2),
            "start": e["start"],
            "end": e["end"],
            "length": e["length"],
        }
        for e in top3_elnino_recon
    ]
)

# ONI is already a 3-month running mean.
elnino_events_oni = find_threshold_events(
    oni_n34, 0.5, warm=True, min_len=5
)
lanina_events_oni = find_threshold_events(
    oni_n34, -0.5, warm=False, min_len=5
)
n_elnino_oni = len(elnino_events_oni)
top3_elnino_oni = sorted(
    elnino_events_oni, key=lambda e: e["peak_value"], reverse=True
)[:3]
top3_elnino_oni_table = pd.DataFrame(
    [
        {
            "year": int(e["peak_time"].year),
            "month": int(e["peak_time"].month),
            "peak_value": round(e["peak_value"], 2),
            "start": e["start"],
            "end": e["end"],
            "length": e["length"],
        }
        for e in top3_elnino_oni
    ]
)

n_elnino_count_difference = n_elnino_recon - n_elnino_oni
count_explanation = (
    "RECON_9120 uses one fixed 1991-2020 monthly climatology for the whole record, "
    "while the official CPC ONI uses a 30-year centred normal that is rolled forward "
    "every five years, together with CPC's operational Nino3.4 processing.  Seasons "
    "close to the +/-0.5 degC threshold can therefore fall on different sides of the "
    "criterion, so the number of qualifying events can differ."
)

# ----------------------------------------------------------------------
# 5. Baseline offsets: 1981-2010 minus 1991-2020
# ----------------------------------------------------------------------
recon_baseline_comp = pd.concat(
    [
        RECON_8110.rename("RECON_8110"),
        RECON_9120.rename("RECON_9120"),
    ],
    axis=1,
    join="inner",
).dropna().loc["1950-01":"2020-12"]

offset_recon_8110_minus_9120 = float(
    (recon_baseline_comp["RECON_8110"]
     - recon_baseline_comp["RECON_9120"]).mean()
)
offset_recon_8110_minus_9120_3 = round(offset_recon_8110_minus_9120, 3)

cpc_baseline_comp = pd.concat(
    [
        cpc_8110_n34.rename("CPC_8110"),
        cpc_9120_n34.rename("CPC_9120"),
    ],
    axis=1,
    join="inner",
).dropna().loc["1950-01":"2020-12"]

offset_cpc_8110_minus_9120 = float(
    (cpc_baseline_comp["CPC_8110"]
     - cpc_baseline_comp["CPC_9120"]).mean()
)
offset_cpc_8110_minus_9120_3 = round(offset_cpc_8110_minus_9120, 3)

offset_agreement_abs = abs(
    offset_recon_8110_minus_9120 - offset_cpc_8110_minus_9120
)
offset_agreement_abs_3 = round(offset_agreement_abs, 3)
