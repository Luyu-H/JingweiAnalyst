# ---------------------------------------------------------------------------
# 1. Load the four SST fields (subset to a box containing the Nino 3.4 region)
# ---------------------------------------------------------------------------
time_range  = ["1950-01", "2014-12"]
base_period = ["1950-01", "2014-12"]

model_names = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
sources     = ["ERSSTv5"] + model_names

box_lat = [-10.0, 10.0]
box_lon = [-180.0, -110.0]      # contains 170W-120W without crossing the dateline

fields = {}
fields["ERSSTv5"] = load_regular_data(
    dataset="ersst", variable="sst", time_range=time_range,
    lat_range=box_lat, lon_range=box_lon)

for m in model_names:
    fields[m] = load_regular_data(
        dataset="cmip6_historical", variable="tos", model=m,
        time_range=time_range, lat_range=box_lat, lon_range=box_lon)

# ---------------------------------------------------------------------------
# 2. Cell areas for each source (each on its own native grid) and the total
#    valid area actually used inside 5S-5N, 170W-120W
# ---------------------------------------------------------------------------
area_fields  = {}
area_source  = {}
total_area   = {}

area_fields["ERSSTv5"] = load_areas(dataset="ersst", variable="sst")
for m in model_names:
    area_fields[m] = load_areas(dataset="cmip6_historical", model=m, variable="tos")

for src in sources:
    f = fields[src]
    lat_name = f.attrs.get("lat_coord", "lat")
    lon_name = f.attrs.get("lon_coord", "lon")

    valid = f.isel(time=0).notnull()          # cells that carry a value
    a     = area_fields[src] * valid          # xarray aligns to the valid cells

    la = valid[lat_name]
    lo = valid[lon_name]
    rmask = (la >= -5) & (la <= 5) & (lo >= -170) & (lo <= -120)

    total_area[src]  = float(a.where(rmask).sum().compute())
    area_source[src] = str(area_fields[src].attrs.get("area_source", "unknown"))

# ---------------------------------------------------------------------------
# 3. Monthly Nino 3.4 anomaly series, each against its own 1950-2014 monthly
#    climatology, area-weighted with its own native cell areas
# ---------------------------------------------------------------------------
series = {}
for src in sources:
    idx = compute_enso_nino_index(
        data=fields[src],
        nino_region="nino3.4",
        depth_level=None,
        anomaly_freq="month",
        climatology="self",
        base_period=base_period,
    )
    series[src] = idx.compute()

# ---------------------------------------------------------------------------
# 4. Statistics: all-month std, calendar-month std, peaks, exceedances
# ---------------------------------------------------------------------------
overall_std  = {}
calendar_std = {}
peak_month   = {}
peak_value   = {}
min_month    = {}
min_value    = {}
peak_ratio   = {}
exceed_count = {}
exceed_pct   = {}
n_months     = {}

for src in sources:
    s = series[src]
    v = np.asarray(s.values, dtype=float)
    n_months[src]    = int(v.size)
    overall_std[src] = float(np.std(v, ddof=0))

    mon = np.asarray(s["time"].dt.month.values, dtype=int)
    ms = [float(np.std(v[mon == m], ddof=0)) for m in range(1, 13)]
    calendar_std[src] = ms

    i_max = int(np.argmax(ms))
    i_min = int(np.argmin(ms))
    peak_month[src] = i_max + 1
    peak_value[src] = ms[i_max]
    min_month[src]  = i_min + 1
    min_value[src]  = ms[i_min]
    peak_ratio[src] = ms[i_max] / ms[i_min]

    cnt = int(np.sum(np.abs(v) > 1.0))
    exceed_count[src] = cnt
    exceed_pct[src]   = cnt / v.size * 100.0

# ---------------------------------------------------------------------------
# 5. Amplitude bias / ratio against ERSSTv5 and model ranking
# ---------------------------------------------------------------------------
obs_std = overall_std["ERSSTv5"]
amp_bias  = {m: overall_std[m] - obs_std for m in model_names}
amp_ratio = {m: overall_std[m] / obs_std for m in model_names}
ranked_models = sorted(model_names, key=lambda m: abs(amp_bias[m]))

# ---------------------------------------------------------------------------
# 6. Printed summary
# ---------------------------------------------------------------------------
print("Nino 3.4 (5S-5N, 170W-120W), 1950-01 to 2014-12, self monthly climatology")
print("-" * 78)
for src in sources:
    print(f"{src:16s} area_source={area_source[src]:12s} "
          f"valid area = {total_area[src]:.3e} m2")
print()

print(f"{'source':16s} {'N':>4s} {'std':>6s} {'peak':>14s} {'min':>14s} {'ratio':>6s} {'>1C':>6s} {'pct':>6s}")
for src in sources:
    print(f"{src:16s} {n_months[src]:4d} {overall_std[src]:6.2f} "
          f"{peak_month[src]:2d} ({peak_value[src]:.2f}) "
          f"{min_month[src]:2d} ({min_value[src]:.2f}) "
          f"{peak_ratio[src]:6.2f} {exceed_count[src]:6d} {exceed_pct[src]:6.1f}")
print()

print("Calendar-month standard deviations (Jan..Dec):")
for src in sources:
    print(f"  {src:16s} " + "  ".join(f"{x:4.2f}" for x in calendar_std[src]))
print()

print("Amplitude vs ERSSTv5:")
for m in model_names:
    print(f"  {m:16s} bias = {amp_bias[m]:+.2f} degC   ratio = {amp_ratio[m]:.2f}")
print()
print("Ranking (smallest -> largest |bias|): " + " < ".join(ranked_models))
