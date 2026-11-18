months_names = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
src_order = ["ERSSTv5", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
time_range = ["1950-01", "2014-12"]
# Slightly larger than the Nino 3.4 box; the index tool masks the exact 5S-5N, 170W-120W region itself
lat_rng = [-10, 10]
lon_rng = [-180, -110]

def to_degC(da):
    u = str(da.attrs.get("units", "")).strip().lower()
    if u in ("k", "kelvin", "degk"):
        da = transform_unit(da, from_unit="K", to_unit="degC")
    return da

raw = {}
raw["ERSSTv5"] = load_regular_data(dataset="ersst", variable="sst", time_range=time_range,
                                   lat_range=lat_rng, lon_range=lon_rng)
for m in ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]:
    raw[m] = load_regular_data(dataset="cmip6_historical", variable="tos", model=m,
                               time_range=time_range, lat_range=lat_rng, lon_range=lon_rng)

units_found = {k: str(v.attrs.get("units", "")) for k, v in raw.items()}

# Monthly Nino 3.4 anomaly, each source against its own 1950-2014 monthly climatology
nino = {}
for k in src_order:
    da = to_degC(raw[k])
    idx = compute_enso_nino_index(da, nino_region="nino3.4", anomaly_freq="month",
                                  climatology="self", base_period=time_range)
    nino[k] = idx.compute()

n_months = {k: int(nino[k].sizes["time"]) for k in src_order}
for k in src_order:
    assert n_months[k] == 780, (k, n_months[k])
    assert float(np.nanmax(np.abs(nino[k].values))) < 50, k  # guard against unmasked fill values
    assert int(np.isnan(nino[k].values).sum()) == 0, k

# 3. Overall amplitude (population std, ddof=0)
std_all = {k: float(compute_basic_stat(nino[k], op="std", dim="time").compute()) for k in src_order}

# 4. Amplitude bias and ratio
bias = {m: std_all[m] - std_all["ERSSTv5"] for m in src_order[1:]}
ratio = {m: std_all[m] / std_all["ERSSTv5"] for m in src_order[1:]}
rank_models = sorted(bias, key=lambda m: abs(bias[m]))

# 5. Seasonal cycle of amplitude: std for each calendar month (ddof=0)
monthly_std = {}
for k in src_order:
    vals = nino[k].values
    mon = nino[k]["time"].dt.month.values
    monthly_std[k] = np.array([vals[mon == mm].std(ddof=0) for mm in range(1, 13)])

seasonal_summary = {}
for k in src_order:
    s = monthly_std[k]
    ipk, imn = int(np.argmax(s)), int(np.argmin(s))
    seasonal_summary[k] = {
        "peak_month": months_names[ipk], "peak_value": float(s[ipk]),
        "min_month": months_names[imn], "min_value": float(s[imn]),
        "max_min_ratio": float(s[ipk] / s[imn]),
    }

obs_s = seasonal_summary["ERSSTv5"]
phase_locking = {}
for m in src_order[1:]:
    ss = seasonal_summary[m]
    same_peak = ss["peak_month"] == obs_s["peak_month"]
    rr = ss["max_min_ratio"] / obs_s["max_min_ratio"]
    if rr > 1.15:
        strength = "stronger phase locking than observed"
    elif rr < 0.85:
        strength = "weaker phase locking than observed"
    else:
        strength = "similar phase-locking strength to observed"
    phase_locking[m] = {"reproduces_peak_month": same_peak, "ratio_vs_obs_ratio": float(rr),
                        "assessment": strength}

# 6. Exceedance counts: months with |anomaly| > 1.0 degC out of 780
exceed = {}
for k in src_order:
    cnt = int((np.abs(nino[k].values) > 1.0).sum())
    exceed[k] = {"count": cnt, "percent": round(100.0 * cnt / 780, 1)}

# Summary table
print("Units found:", units_found)
for k in src_order:
    print(f"{k}: std={std_all[k]:.2f}, exceed={exceed[k]['count']} ({exceed[k]['percent']:.1f}%), "
          f"peak={seasonal_summary[k]['peak_month']} {seasonal_summary[k]['peak_value']:.2f}, "
          f"min={seasonal_summary[k]['min_month']} {seasonal_summary[k]['min_value']:.2f}, "
          f"ratio={seasonal_summary[k]['max_min_ratio']:.2f}")
    print("   monthly std:", np.round(monthly_std[k], 2))
for m in src_order[1:]:
    print(f"{m}: bias={bias[m]:.2f}, ratio={ratio[m]:.2f}, phase={phase_locking[m]}")
print("Rank (smallest to largest |bias|):", rank_models)
