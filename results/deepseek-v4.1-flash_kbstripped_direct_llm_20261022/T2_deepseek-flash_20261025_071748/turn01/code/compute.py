import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

DATA = "/data1/lorraine-20419/ocean_data"

# ==================================================================
# 1. ERSSTv5 monthly SST  ->  Nino3.4 region mean (1950-01..2023-12)
# ==================================================================
ds_sst = xr.open_dataset(f"{DATA}/ERSSTv5/sst.mnmean.nc")
sst = ds_sst["sst"]

# native grid: lon 0..358 (0-360 convention), lat descending 88..-88
# Nino3.4 = 5S-5N, 170W-120W  ==  5S-5N, 190E-240E
n34_box = sst.sel(lat=slice(5, -5), lon=slice(190, 240))
n34_da = n34_box.mean(dim=("lat", "lon")).compute()      # simple (unweighted) box mean, as CPC does

tv = n34_da["time"].values
try:
    tidx = pd.DatetimeIndex(tv)
except Exception:                                          # cftime objects
    tidx = pd.DatetimeIndex([pd.Timestamp(int(x.year), int(x.month), 1) for x in tv])

SST_N34 = pd.Series(np.asarray(n34_da.values, dtype=float), index=tidx, name="ersst5_n34")
SST_N34 = SST_N34.loc["1950-01-01":"2023-12-31"].copy()
SST_N34.index = SST_N34.index.to_period("M")
SST_N34 = SST_N34.reindex(pd.period_range("1950-01", "2023-12", freq="M"))

# --- monthly climatologies on the two baselines --------------------
def monthly_clim(s, y0, y1):
    sub = s[(s.index.year >= y0) & (s.index.year <= y1)]
    return sub.groupby(sub.index.month).mean()

clim_9120 = monthly_clim(SST_N34, 1991, 2020)
clim_8110 = monthly_clim(SST_N34, 1981, 2010)

def anomalies(s, clim):
    return pd.Series(s.values - clim.reindex(s.index.month).values, index=s.index)

RECON_9120 = anomalies(SST_N34, clim_9120).rename("RECON_9120")   # main reconstruction
RECON_8110 = anomalies(SST_N34, clim_8110).rename("RECON_8110")   # 1981-2010 control

# ==================================================================
# 2. CPC published Nino indices (1991-2020 and 1981-2010 baselines)
# ==================================================================
def read_nino_monthly(path):
    rows = []
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) < 10:
                continue
            try:
                yr, mo = int(p[0]), int(p[1])
                vals = [float(x) for x in p[2:10]]
            except ValueError:                            # header line
                continue
            rows.append([yr, mo] + vals)
    a = np.array(rows, dtype=float)
    a[a <= -90.0] = np.nan                                # -999.9 missing flags
    idx = pd.PeriodIndex(pd.to_datetime(dict(year=a[:, 0].astype(int),
                                             month=a[:, 1].astype(int), day=1)), freq="M")
    cols = ["N12", "N12_ANOM", "N3", "N3_ANOM", "N4", "N4_ANOM", "N34", "N34_ANOM"]
    return pd.DataFrame(a[:, 2:10], index=idx, columns=cols).sort_index()

CPC = read_nino_monthly(f"{DATA}/CPC/ersst5.nino.mth.91-20.ascii")
CPC_9120 = CPC["N34_ANOM"].rename("CPC_9120")
CPC81 = read_nino_monthly(f"{DATA}/CPC/ersst5.nino.mth.81-10.ascii")
CPC_8110 = CPC81["N34_ANOM"].rename("CPC_8110")

# --- CPC ONI (already a 3-month running mean, labelled by season) ---
SEASONS = ["DJF", "JFM", "FMA", "MAM", "AMJ", "MJJ",
           "JJA", "JAS", "ASO", "SON", "OND", "NDJ"]
SEAS_CENTER = {s: i + 1 for i, s in enumerate(SEASONS)}   # DJF -> Jan, JFM -> Feb, ...

rows = []
with open(f"{DATA}/CPC/oni.ascii.txt") as f:
    for line in f:
        p = line.split()
        if len(p) < 4:
            continue
        seas = p[0].upper()
        if seas not in SEAS_CENTER:
            continue
        try:
            yr, anom = int(p[1]), float(p[3])
        except ValueError:
            continue
        rows.append((yr, SEAS_CENTER[seas], anom))

o = pd.DataFrame(rows, columns=["YR", "MON", "ANOM"])
o = o[o["ANOM"] > -90]
ONI = pd.Series(o["ANOM"].values,
                index=pd.PeriodIndex(pd.to_datetime(dict(year=o.YR, month=o.MON, day=1)), freq="M"),
                name="ONI").sort_index()

# ==================================================================
# 3. RECON_9120 vs CPC 1991-2020, 1950-01..2023-12
# ==================================================================
idx = pd.period_range("1950-01", "2023-12", freq="M")
R = RECON_9120.reindex(idx)
C = CPC_9120.reindex(idx)
good = R.notna() & C.notna()
DIFF = (R - C)[good]

corr_recon_cpc = float(np.corrcoef(R[good].values, C[good].values)[0, 1])
rmse_recon_cpc = float(np.sqrt(np.mean(DIFF.values ** 2)))
absdiff = DIFF.abs()
maxdev_time = absdiff.idxmax()                 # Period (year, month) of the largest deviation
maxdev_recon_cpc = float(absdiff.max())

# ==================================================================
# 4. ONI-criterion event detection (>=5 consecutive 3-month seasons >= +0.5)
# ==================================================================
def run3(s):
    return s.rolling(3, center=True).mean()

def detect_events(series, thresh=0.5, min_len=5):
    s = series.dropna().sort_index()
    v = s.values
    n = len(v)
    out = []
    i = 0
    while i < n:
        if v[i] >= thresh:
            j = i
            while j + 1 < n and v[j + 1] >= thresh:
                j += 1
            if (j - i + 1) >= min_len:
                seg = s.iloc[i:j + 1]
                out.append({"start": s.index[i], "end": s.index[j],
                            "peak_time": seg.idxmax(), "peak": float(seg.max()),
                            "n_seasons": int(j - i + 1)})
            i = j + 1
        else:
            i += 1
    return out

RECON3 = run3(RECON_9120)                        # 3-month centred running mean of the reconstruction
EN_RECON = detect_events(RECON3, +0.5, 5)        # El Nino events
LN_RECON = detect_events(-RECON3, +0.5, 5)       # La Nina events (same rule, sign flipped)
en_count_recon = len(EN_RECON)
top3_recon = sorted(EN_RECON, key=lambda e: e["peak"], reverse=True)[:3]

# ONI on the identical window (seasons centred 1950-01..2023-12)
ONI_W = ONI[(ONI.index >= pd.Period("1950-01", "M")) & (ONI.index <= pd.Period("2023-12", "M"))]
EN_ONI = detect_events(ONI_W, +0.5, 5)
LN_ONI = detect_events(-ONI_W, +0.5, 5)
en_count_oni = len(EN_ONI)
top3_oni = sorted(EN_ONI, key=lambda e: e["peak"], reverse=True)[:3]

# ==================================================================
# 5. Baseline offsets, 1950-01..2020-12
# ==================================================================
win = pd.period_range("1950-01", "2020-12", freq="M")
offset_recon = float((RECON_8110.reindex(win) - RECON_9120.reindex(win)).mean())
offset_published = float((CPC_8110.reindex(win) - CPC_9120.reindex(win)).mean())
offset_agreement = offset_recon - offset_published

# ------------------------- reporting ------------------------------
print("RECON_9120 vs CPC 1991-2020 Nino3.4 anomaly (1950-01..2023-12, %d months)" % good.sum())
print("  Pearson r        : %.6f" % corr_recon_cpc)
print("  RMSE             : %.4f degC" % rmse_recon_cpc)
print("  max |deviation|  : %.4f degC at %04d-%02d" % (maxdev_recon_cpc, maxdev_time.year, maxdev_time.month))
print()
print("El Nino events (>=5 consecutive 3-month seasons >= +0.5 degC)")
print("  RECON_9120 : %d" % en_count_recon)
for e in top3_recon:
    print("     peak %04d-%02d  %.2f degC  (event %s..%s, %d seasons)" %
          (e["peak_time"].year, e["peak_time"].month, e["peak"], e["start"], e["end"], e["n_seasons"]))
print("  CPC ONI    : %d" % en_count_oni)
for e in top3_oni:
    print("     peak %04d-%02d  %.2f degC  (event %s..%s, %d seasons)" %
          (e["peak_time"].year, e["peak_time"].month, e["peak"], e["start"], e["end"], e["n_seasons"]))
print("  La Nina events: RECON %d, ONI %d" % (len(LN_RECON), len(LN_ONI)))
print()
print("Mean baseline offset 1981-2010 minus 1991-2020 (1950-01..2020-12)")
print("  reconstructed : %+.3f degC" % offset_recon)
print("  published CPC : %+.3f degC" % offset_published)
print("  agreement     : %+.3f degC" % offset_agreement)
