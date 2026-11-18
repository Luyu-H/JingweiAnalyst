# ============================================================
# Reconstruct the Nino 3.4 index from raw ERSSTv5 SST
# ============================================================

ersst_path = "/data1/lorraine-20419/ocean_data/ERSSTv5/sst.mnmean.nc"
cpc_dir = "/data1/lorraine-20419/ocean_data/CPC"

ds = xr.open_dataset(ersst_path, decode_times=True)
sst = ds["sst"].sel(time=slice("1950-01-01", "2023-12-31"))

# Put latitude in ascending order for ordinary slice selection.
if float(sst["lat"][0]) > float(sst["lat"][-1]):
    sst = sst.sortby("lat")

# Bring longitude to 0-360 if the file is not already there.
if float(sst["lon"].min()) < 0.0:
    sst = sst.assign_coords(lon=(sst["lon"] % 360.0)).sortby("lon")

# Nino 3.4: 5S-5N, 170W-120W = 190E-240E on a 0-360 longitude axis.
sst_box = sst.sel(lat=slice(-5.0, 5.0), lon=slice(190.0, 240.0)).load()
ds.close()


def reconstruct_nino34(sst_box, baseline_start, baseline_end):
    """Area-weighted Nino 3.4 SST anomaly from grid-cell monthly climatology."""
    base = sst_box.sel(
        time=slice(f"{baseline_start}-01-01", f"{baseline_end}-12-31")
    )
    # Calendar-month climatology at each grid cell over the requested baseline.
    clim = base.groupby("time.month").mean(dim="time", skipna=True)
    anom = sst_box.groupby("time.month") - clim

    # Cosine-latitude cell area weights, restricted to valid non-missing cells.
    w = np.cos(np.deg2rad(sst_box["lat"]))
    valid = anom.notnull()
    num = (anom * w).sum(dim=("lat", "lon"), skipna=True)
    den = (w * valid).sum(dim=("lat", "lon"))
    out = (num / den).rename("nino34_anom")

    out.attrs.update(
        {
            "description": "Nino 3.4 SST anomaly reconstructed from ERSSTv5",
            "region": "5S-5N, 170W-120W",
            "baseline": f"{baseline_start}-{baseline_end}",
            "climatology": (
                f"grid-cell calendar-month climatology over "
                f"{baseline_start}-{baseline_end}"
            ),
            "area_weighting": "cos(latitude), valid cells only",
            "units": "degC",
        }
    )
    return out


recon_9120 = reconstruct_nino34(sst_box, 1991, 2020)
recon_8110 = reconstruct_nino34(sst_box, 1981, 2010)

recon_9120_s = recon_9120.to_series().sort_index()
recon_8110_s = recon_8110.to_series().sort_index()


# ============================================================
# Read CPC monthly Nino and ONI ASCII products
# ============================================================

def read_cpc_nino(path, region="NINO3.4"):
    """Read a CPC monthly Nino index file and return the selected region anomaly."""
    header = None
    records = []

    with open(path, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()

            if header is None and any(
                p.upper() in ("YR", "YEAR") for p in parts[:3]
            ):
                header = parts
                continue

            if len(parts) < 3:
                continue

            try:
                yr = int(float(parts[0]))
                mon = int(float(parts[1]))
                vals = [float(x) for x in parts[2:]]
            except (ValueError, IndexError):
                continue

            records.append([yr, mon] + vals)

    if header is None:
        header = [
            "YR", "MON",
            "NINO1+2", "ANOM",
            "NINO3", "ANOM",
            "NINO4", "ANOM",
            "NINO3.4", "ANOM",
        ]

    header_u = [h.upper() for h in header]

    reg_idx = None
    for i, h in enumerate(header_u):
        if h == region.upper():
            reg_idx = i
            break
    if reg_idx is None:
        # Fallback for the standard CPC Nino layout.
        reg_idx = 8

    anom_idx = None
    for j in range(reg_idx + 1, len(header_u)):
        if header_u[j] == "ANOM":
            anom_idx = j
            break
    if anom_idx is None:
        anom_idx = reg_idx + 1

    arr = np.asarray(records, dtype=float)
    if arr.size == 0:
        return pd.Series(dtype=float, name=f"{region}_ANOM")

    time = pd.to_datetime(
        {
            "year": arr[:, 0].astype(int),
            "month": arr[:, 1].astype(int),
            "day": np.ones(len(arr), dtype=int),
        }
    )
    values = arr[:, anom_idx]
    s = pd.Series(values, index=pd.DatetimeIndex(time)).sort_index()
    s = s.where(s > -90.0)
    s.name = f"{region}_ANOM"
    return s


def read_oni(path):
    """Read CPC ONI seasonal ASCII and map each season to its centre month."""
    season_month = {
        "DJF": 1, "JFM": 2, "FMA": 3, "MAM": 4,
        "AMJ": 5, "MJJ": 6, "JJA": 7, "JAS": 8,
        "ASO": 9, "SON": 10, "OND": 11, "NDJ": 12,
    }
    times = []
    values = []

    with open(path, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            seas = parts[0].upper()
            if seas in ("SEAS", "SEASON"):
                continue
            if seas not in season_month:
                continue
            try:
                yr = int(float(parts[1]))
                anom = float(parts[3])
            except (ValueError, IndexError):
                continue
            times.append(
                pd.Timestamp(year=yr, month=season_month[seas], day=1)
            )
            values.append(anom)

    s = pd.Series(values, index=pd.DatetimeIndex(times)).sort_index()
    s = s.where(s > -90.0)
    s.name = "ONI_ANOM"
    return s


cpc_9120_s = read_cpc_nino(
    f"{cpc_dir}/ersst5.nino.mth.91-20.ascii", region="NINO3.4"
)
cpc_8110_s = read_cpc_nino(
    f"{cpc_dir}/ersst5.nino.mth.81-10.ascii", region="NINO3.4"
)
oni_s = read_oni(f"{cpc_dir}/oni.ascii.txt")


# ============================================================
# Compare RECON_9120 with CPC 1991-2020 Nino 3.4 anomaly
# ============================================================

comp = pd.concat(
    [
        recon_9120_s.rename("RECON_9120"),
        cpc_9120_s.rename("CPC_9120"),
    ],
    axis=1,
    join="inner",
).dropna()
comp = comp.loc["1950-01-01":"2023-12-01"]

r_vals = comp["RECON_9120"].values
c_vals = comp["CPC_9120"].values

corr = float(np.corrcoef(r_vals, c_vals)[0, 1])
rmse = float(np.sqrt(np.mean((r_vals - c_vals) ** 2)))

abs_dev = np.abs(r_vals - c_vals)
imax = int(np.argmax(abs_dev))
max_abs_dev = float(abs_dev[imax])
max_abs_dev_time = comp.index[imax]


# ============================================================
# Baseline offsets
# ============================================================

off_rec = pd.concat(
    [
        recon_8110_s.rename("RECON_8110"),
        recon_9120_s.rename("RECON_9120"),
    ],
    axis=1,
    join="inner",
).dropna()
off_rec = off_rec.loc["1950-01-01":"2020-12-01"]
recon_offset = float(
    (off_rec["RECON_8110"] - off_rec["RECON_9120"]).mean()
)

off_pub = pd.concat(
    [
        cpc_8110_s.rename("CPC_8110"),
        cpc_9120_s.rename("CPC_9120"),
    ],
    axis=1,
    join="inner",
).dropna()
off_pub = off_pub.loc["1950-01-01":"2020-12-01"]
published_offset = float(
    (off_pub["CPC_8110"] - off_pub["CPC_9120"]).mean()
)

offset_difference = published_offset - recon_offset


# ============================================================
# Event detection: contiguous runs of the 3-month running mean
# ============================================================

def detect_threshold_events(series, threshold, min_count=5, warm=True):
    """Detect contiguous threshold exceedances of at least min_count values."""
    s = series.dropna().sort_index()
    vals = s.values
    times = s.index

    if warm:
        mask = vals >= threshold
    else:
        mask = vals <= threshold

    events = []
    i = 0
    n = len(vals)

    while i < n:
        if not mask[i]:
            i += 1
            continue

        j = i
        while j + 1 < n and mask[j + 1]:
            j += 1

        length = j - i + 1
        if length >= min_count:
            seg = vals[i:j + 1]
            seg_times = times[i:j + 1]
            if warm:
                k = int(np.nanargmax(seg))
            else:
                k = int(np.nanargmin(seg))
            events.append(
                {
                    "start": times[i],
                    "end": times[j],
                    "duration": length,
                    "peak_time": seg_times[k],
                    "peak_value": float(seg[k]),
                }
            )
        i = j + 1

    return events


run3_recon = recon_9120_s.rolling(
    window=3, center=True, min_periods=3
).mean()
run3_recon = run3_recon.loc["1950-01-01":"2023-12-01"]

recon_elnino_events = detect_threshold_events(
    run3_recon, 0.5, min_count=5, warm=True
)
recon_lanina_events = detect_threshold_events(
    run3_recon, -0.5, min_count=5, warm=False
)

recon_elnino_sorted = sorted(
    recon_elnino_events, key=lambda e: e["peak_value"], reverse=True
)
recon_top3 = recon_elnino_sorted[:3]

oni_2023 = oni_s.loc["1950-01-01":"2023-12-01"]
oni_elnino_events = detect_threshold_events(
    oni_2023, 0.5, min_count=5, warm=True
)
oni_elnino_sorted = sorted(
    oni_elnino_events, key=lambda e: e["peak_value"], reverse=True
)
oni_top3 = oni_elnino_sorted[:3]

recon_cross_months = int((run3_recon >= 0.5).sum())
oni_cross_months = int((oni_2023 >= 0.5).sum())

count_note = (
    "RECON_9120 uses a fixed 1991-2020 calendar-month climatology. "
    "The official ONI uses its operational rolling 30-year base-period "
    "convention (updated every 5 years) and is distributed rounded to "
    "0.1 degC. Marginal threshold crossings therefore differ, which can "
    "change how many 5-season runs are counted."
)


# ============================================================
# Report
# ============================================================

def peak_str(ev):
    return (
        f"{ev['peak_value']:.2f} degC at "
        f"{ev['peak_time'].strftime('%Y-%m')}"
    )


print("=== RECON_9120 vs CPC 1991-2020 Nino 3.4 ===")
print(f"Pearson correlation = {corr:.6f}")
print(f"RMSE = {rmse:.4f} degC")
print(
    f"Maximum absolute deviation = {max_abs_dev:.4f} degC "
    f"at {max_abs_dev_time.strftime('%Y-%m')}"
)
print()

print("=== El Nino events, 1950-2023 ===")
print(f"RECON_9120: {len(recon_elnino_events)} events")
for ev in recon_top3:
    print(f"  {peak_str(ev)}")
print(f"CPC ONI: {len(oni_elnino_events)} events")
for ev in oni_top3:
    print(f"  {peak_str(ev)}")
print(f"Threshold-crossing months (>= +0.5 degC): "
      f"RECON_9120={recon_cross_months}, ONI={oni_cross_months}")
print()

print("=== Baseline offsets, 1950-2020 ===")
print(f"RECON_8110 minus RECON_9120 = {recon_offset:.3f} degC")
print(
    f"CPC 1981-2010 minus CPC 1991-2020 = "
    f"{published_offset:.3f} degC"
)
print(
    f"Published minus reconstructed offset = "
    f"{offset_difference:.3f} degC"
)
print()
print("Count comparison:", count_note)

result = {
    "correlation": corr,
    "rmse": rmse,
    "max_abs_dev": max_abs_dev,
    "max_abs_dev_time": max_abs_dev_time,
    "recon_elnino_count": len(recon_elnino_events),
    "recon_top3": recon_top3,
    "oni_elnino_count": len(oni_elnino_events),
    "oni_top3": oni_top3,
    "recon_offset": recon_offset,
    "published_offset": published_offset,
    "offset_difference": offset_difference,
}
