# Paths from the dataset registry
ERSST_PATH = "/data1/lorraine-20419/ocean_data/ERSSTv5/sst.mnmean.nc"
CPC_ROOT = "/data1/lorraine-20419/ocean_data/CPC"


def read_nino_monthly(path):
    """Read one CPC monthly Nino file and return the Nino3.4 anomaly series."""
    rows = []
    with open(path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 10:
                continue
            try:
                yr = int(float(parts[0]))
                mon = int(float(parts[1]))
            except (ValueError, TypeError):
                continue
            vals = []
            for p in parts[:10]:
                try:
                    vals.append(float(p))
                except ValueError:
                    vals.append(np.nan)
            rows.append(vals)

    df = pd.DataFrame(
        rows,
        columns=[
            "YR", "MON",
            "NINO1+2", "ANOM1",
            "NINO3", "ANOM2",
            "NINO4", "ANOM3",
            "NINO3.4", "ANOM4",
        ],
    )
    # The raw files use a large negative missing value; anomalies are small.
    df = df.mask(df < -90.0).dropna(subset=["YR", "MON"])
    df["time"] = pd.to_datetime(
        dict(year=df["YR"].astype(int), month=df["MON"].astype(int), day=1)
    )
    s = df.set_index("time")["ANOM4"].sort_index()
    s.name = "NINO3.4_ANOM"
    s.attrs["units"] = "degC"
    return s


def read_oni(path):
    """Read CPC seasonal ONI and map each 3-month season to its centre month."""
    season_to_month = {
        "DJF": 1, "JFM": 2, "FMA": 3, "MAM": 4, "AMJ": 5, "MJJ": 6,
        "JJA": 7, "JAS": 8, "ASO": 9, "SON": 10, "OND": 11, "NDJ": 12,
    }
    rows = []
    with open(path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            seas = parts[0].upper()
            if seas not in season_to_month:
                continue
            try:
                yr = int(float(parts[1]))
                total = float(parts[2])
                anom = float(parts[3])
            except (ValueError, TypeError):
                continue
            if yr < 0:
                continue
            if anom < -90.0:
                anom = np.nan
            rows.append((seas, yr, total, anom))

    df = pd.DataFrame(rows, columns=["SEAS", "YR", "TOTAL", "ANOM"])
    df["month"] = df["SEAS"].map(season_to_month)
    df["time"] = pd.to_datetime(
        dict(year=df["YR"].astype(int), month=df["month"].astype(int), day=1)
    )
    s = df.set_index("time")["ANOM"].sort_index()
    s.name = "ONI_ANOM"
    s.attrs["units"] = "degC"
    return s


def detect_threshold_events(series, threshold, min_duration=5, direction="above"):
    """Find contiguous threshold exceedances of at least min_duration points."""
    s = pd.Series(series).sort_index()
    if direction == "above":
        mask = s.notna() & (s >= threshold)
    elif direction == "below":
        mask = s.notna() & (s <= threshold)
    else:
        raise ValueError("direction must be 'above' or 'below'")

    m = mask.to_numpy()
    times = s.index
    n = len(s)
    events = []
    i = 0
    while i < n:
        if m[i]:
            start = i
            while i + 1 < n and m[i + 1]:
                i += 1
            end = i
            length = end - start + 1
            if length >= min_duration:
                seg = s.iloc[start:end + 1]
                if direction == "above":
                    peak_time = seg.idxmax()
                    peak_value = float(seg.max())
                else:
                    peak_time = seg.idxmin()
                    peak_value = float(seg.min())
                events.append({
                    "start": times[start],
                    "end": times[end],
                    "length": int(length),
                    "peak_time": peak_time,
                    "peak_value": peak_value,
                })
        i += 1
    return events


def events_to_mask(series, events):
    mask = pd.Series(False, index=series.index)
    for e in events:
        mask.loc[e["start"]:e["end"]] = True
    return mask


# ---------------------------------------------------------------------------
# 1. Open ERSSTv5 and build the Nino3.4 region-mean SST series
# ---------------------------------------------------------------------------
ds = xr.open_dataset(ERSST_PATH)
sst = ds["sst"].sel(time=slice("1950-01-01", "2023-12-31"))

lat = sst["lat"]
lon = sst["lon"]

# Nino3.4: 5S-5N, 170W-120W.  Native ERSST lon is 0-360, so 170W=190E,
# 120W=240E.  Use %360 in case a file variant is -180..180.
lat_mask = (lat >= -5) & (lat <= 5)
lon_360 = lon % 360
lon_mask = (lon_360 >= 190) & (lon_360 <= 240)
region_mask = lat_mask & lon_mask

sst_box = sst.where(region_mask)
weights = np.cos(np.deg2rad(lat))
valid = sst_box.notnull().astype(float)

num = (sst_box * weights).sum(dim=("lat", "lon"), skipna=True)
den = (valid * weights).sum(dim=("lat", "lon"))
nino34_sst_da = num / den
nino34_sst_da.name = "nino34_sst"
nino34_sst_da.attrs["units"] = "degC"
nino34_sst = nino34_sst_da.to_series().sort_index()

# ---------------------------------------------------------------------------
# 2. Build monthly anomalies relative to each baseline
# ---------------------------------------------------------------------------
def build_recon_anomaly(series, base_start, base_end, base_label, name):
    s = series.copy().sort_index()
    base = s.loc[base_start:base_end]
    clim = base.groupby(base.index.month).mean()
    anom_values = s.to_numpy() - clim.reindex(s.index.month).to_numpy()
    anom = pd.Series(anom_values, index=s.index, name=name)
    anom.attrs["base_period"] = base_label
    anom.attrs["climatology"] = clim.to_dict()
    anom.attrs["units"] = "degC"
    anom.attrs["long_name"] = "Nino 3.4 SST anomaly reconstructed from ERSSTv5"
    return anom, clim


recon_9120, clim_9120 = build_recon_anomaly(
    nino34_sst, "1991-01-01", "2020-12-01", "1991-2020", "RECON_9120"
)
recon_8110, clim_8110 = build_recon_anomaly(
    nino34_sst, "1981-01-01", "2010-12-01", "1981-2010", "RECON_8110"
)

# ---------------------------------------------------------------------------
# 3. Read CPC published Nino3.4 anomalies and ONI
# ---------------------------------------------------------------------------
cpc_9120 = read_nino_monthly(CPC_ROOT + "/ersst5.nino.mth.91-20.ascii")
cpc_8110 = read_nino_monthly(CPC_ROOT + "/ersst5.nino.mth.81-10.ascii")
oni = read_oni(CPC_ROOT + "/oni.ascii.txt")

cpc_9120.attrs["base_period"] = "1991-2020"
cpc_8110.attrs["base_period"] = "1981-2010"
oni.attrs["base_period"] = "published ONI rolling 30-year centred normal"

# Align reconstructed and CPC 1991-2020 series on common months, 1950-2023.
common_idx = recon_9120.index.intersection(cpc_9120.index)
common_idx = common_idx[
    (common_idx >= pd.Timestamp("1950-01-01"))
    & (common_idx <= pd.Timestamp("2023-12-01"))
]
recon_9120_common = recon_9120.reindex(common_idx)
cpc_9120_common = cpc_9120.reindex(common_idx)

valid_pair = pd.concat(
    [
        recon_9120_common.rename("recon"),
        cpc_9120_common.rename("cpc"),
    ],
    axis=1,
).dropna()

# ---------------------------------------------------------------------------
# 4. Comparison metrics: Pearson r, RMSE, largest absolute deviation
# ---------------------------------------------------------------------------
r_raw = float(valid_pair["recon"].corr(valid_pair["cpc"]))
rmse_raw = float(np.sqrt(np.mean((valid_pair["recon"] - valid_pair["cpc"]) ** 2)))

abs_diff = (valid_pair["recon"] - valid_pair["cpc"]).abs()
max_idx = abs_diff.idxmax()
max_abs_dev_raw = float(abs_diff.loc[max_idx])
max_deviation_year = int(max_idx.year)
max_deviation_month = int(max_idx.month)

correlation_r = round(r_raw, 6)
rmse_c = round(rmse_raw, 4)
max_abs_deviation_c = round(max_abs_dev_raw, 4)

diff_9120_cpc = (recon_9120_common - cpc_9120_common).dropna()
diff_9120_cpc.name = "RECON_9120_minus_CPC"

# ---------------------------------------------------------------------------
# 5. 3-month centred running mean and threshold events for RECON_9120
# ---------------------------------------------------------------------------
rm_recon = recon_9120.rolling(window=3, center=True, min_periods=3).mean()
rm_recon.name = "RECON_9120_3mo_running_mean"

elnino_events_recon = detect_threshold_events(
    rm_recon, 0.5, min_duration=5, direction="above"
)
lanina_events_recon = detect_threshold_events(
    rm_recon, -0.5, min_duration=5, direction="below"
)

n_elnino_recon = len(elnino_events_recon)
n_lanina_recon = len(lanina_events_recon)

top3_recon = sorted(
    elnino_events_recon, key=lambda e: e["peak_value"], reverse=True
)[:3]
top3_recon_report = [
    {
        "year": int(e["peak_time"].year),
        "month": int(e["peak_time"].month),
        "peak_value_c": round(float(e["peak_value"]), 2),
    }
    for e in top3_recon
]

warm_mask_recon = events_to_mask(rm_recon, elnino_events_recon)
cold_mask_recon = events_to_mask(rm_recon, lanina_events_recon)

# ---------------------------------------------------------------------------
# 6. Same event criterion applied directly to CPC ONI
# ---------------------------------------------------------------------------
elnino_events_oni = detect_threshold_events(
    oni, 0.5, min_duration=5, direction="above"
)
lanina_events_oni = detect_threshold_events(
    oni, -0.5, min_duration=5, direction="below"
)

n_elnino_oni = len(elnino_events_oni)
n_lanina_oni = len(lanina_events_oni)

top3_oni = sorted(
    elnino_events_oni, key=lambda e: e["peak_value"], reverse=True
)[:3]
top3_oni_report = [
    {
        "year": int(e["peak_time"].year),
        "month": int(e["peak_time"].month),
        "peak_value_c": round(float(e["peak_value"]), 2),
    }
    for e in top3_oni
]

oni_first_center = oni.index.min()
oni_last_center = oni.index.max()

count_explanation = (
    f"RECON_9120 gives {n_elnino_recon} El Nino events over 1950-2023 using a fixed "
    f"1991-2020 monthly climatology and a 3-month centred running mean; CPC ONI gives "
    f"{n_elnino_oni} events over its available record "
    f"({oni_first_center.strftime('%b %Y')} to {oni_last_center.strftime('%b %Y')}). "
    f"The difference is expected because (i) ONI uses its own published rolling 30-year "
    f"centred normal rather than the fixed 1991-2020 baseline of RECON_9120; "
    f"(ii) ONI is distributed at coarse precision relative to the raw ERSST-derived "
    f"reconstruction; and (iii) the centred running mean has no values at 1950-01 or "
    f"2023-12 for RECON_9120, whereas ONI starts at DJF 1950 (centred Jan 1950) and ends "
    f"at JJA 2026 (centred Jul 2026). Events straddling those record edges are therefore "
    f"truncated or absent in one series but not the other."
)

# ---------------------------------------------------------------------------
# 7. Baseline offsets over 1950-2020 (unweighted over months)
# ---------------------------------------------------------------------------
offset_period = slice("1950-01-01", "2020-12-01")

recon_offset_raw = float(
    (recon_8110 - recon_9120).loc[offset_period].dropna().mean()
)

published_pair = pd.concat(
    [
        cpc_8110.rename("cpc8110"),
        cpc_9120.rename("cpc9120"),
    ],
    axis=1,
).dropna().loc[offset_period]
published_offset_raw = float(
    (published_pair["cpc8110"] - published_pair["cpc9120"]).mean()
)

offset_difference_raw = recon_offset_raw - published_offset_raw

recon_offset_c = round(recon_offset_raw, 3)
published_offset_c = round(published_offset_raw, 3)
offset_difference_c = round(offset_difference_raw, 3)
offset_agreement_abs_c = round(abs(offset_difference_raw), 3)

# ---------------------------------------------------------------------------
# 8. Collect the requested results
# ---------------------------------------------------------------------------
result = {
    "correlation_r": correlation_r,
    "rmse_c": rmse_c,
    "max_abs_deviation_c": max_abs_deviation_c,
    "max_deviation_year": max_deviation_year,
    "max_deviation_month": max_deviation_month,
    "n_elnino_recon": n_elnino_recon,
    "n_lanina_recon": n_lanina_recon,
    "top3_recon": top3_recon_report,
    "n_elnino_oni": n_elnino_oni,
    "n_lanina_oni": n_lanina_oni,
    "top3_oni": top3_oni_report,
    "recon_offset_c": recon_offset_c,
    "published_offset_c": published_offset_c,
    "offset_difference_c": offset_difference_c,
    "offset_agreement_abs_c": offset_agreement_abs_c,
    "count_explanation": count_explanation,
}
