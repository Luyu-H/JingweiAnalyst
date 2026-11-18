"""
Category 3: reproducing published products.
"""
from __future__ import annotations

import functools
import os

import numpy as np
import pandas as pd
import xarray as xr

from src.eval.tool_verification import refs
from src.eval.tool_verification.common import case
from src.tools.basic_calculator import compute_basic_stat, detect_threshold_events, temporal_aggregate
from src.tools.data_loader import load_areas, load_climatology_data, load_regular_data
from src.tools.data_openers import get_dataset_cfg
from src.tools.scientific_calculator import compute_enso_nino_index, compute_mixed_layer_depth, compute_pdo_index

CAT = "cat3_products"
CPC_ROOT = get_dataset_cfg("cpc_indices")["root"]
END = "2025-12"          # the last months of the record are still being revised by ERSSTv5
REGIONS = {"nino1+2": ("NINO1+2", (-10, 0), (-90, -80)), "nino3": ("NINO3", (-5, 5), (-150, -90)),
           "nino3.4": ("NINO3.4", (-5, 5), (-170, -120)), "nino4": ("NINO4", (-5, 5), (160, -150))}


@functools.lru_cache(None)
def ersst(start, end):
    return load_regular_data("ersst", "sst", time_range=[start, end]).load()


def agreement(rec, label, tool, published, tol=0.01, share=0.95, tier="T2", note=""):
    """Compare two series on their common months and record the agreement statistics."""
    both = pd.concat([pd.Series(tool), pd.Series(published)], axis=1, join="inner").dropna()
    d = (both.iloc[:, 0] - both.iloc[:, 1]).values
    within = float(np.mean(np.abs(d) <= tol + 1e-9))
    rec.metric(f"{label}: share of months within {tol} degC", within, threshold=share,
               op="info" if share is None else ">=", tier=tier,
               note=f"{len(d)} months; RMSE {np.sqrt(np.mean(d ** 2)):.4f}, max |d| {np.abs(d).max():.3f} "
                    f"({both.index[int(np.argmax(np.abs(d)))].strftime('%Y-%m')}), mean d {d.mean():+.4f}"
                    + (f"; {note}" if note else ""))
    return d


def cpc_series(file, column):
    raw = refs.raw_cpc_monthly(os.path.join(CPC_ROOT, file))
    index = pd.to_datetime([f"{y}-{m:02d}-01" for (y, m) in raw])
    return pd.Series([v[column] for v in raw.values()], index=index).sort_index()


def series(da):
    return pd.Series(np.asarray(da.values, dtype="float64"), index=pd.DatetimeIndex(da["time"].values))


@case("compute_enso_nino_index", "CPC ERSSTv5 Nino indices: area-mean SST (TOTAL) and anomalies on two baselines")
def cpc_nino(rec):
    sst = ersst("1950-01", END)
    for region, (name, lat_range, lon_range) in REGIONS.items():
        box = load_regular_data("ersst", "sst", time_range=["1950-01", END], lat_range=list(lat_range),
                                lon_range=list(lon_range))
        total = series(compute_basic_stat(box, "area_weighted_mean", dim=["lat", "lon"]).load())
        published = cpc_series("ersst5.nino.mth.91-20.ascii", name)
        d_incl = agreement(rec, f"{name} TOTAL (box mean SST, edge cells included)", total, published)
        # The 2-degree centres fall exactly on the box edges. Shrinking the box by
        # 0.1 degree drops the edge rows/columns, the other reading of the convention.
        inner = load_regular_data("ersst", "sst", time_range=["1950-01", END],
                                  lat_range=[lat_range[0] + 0.1, lat_range[1] - 0.1],
                                  lon_range=[lon_range[0] + 0.1, lon_range[1] - 0.1])
        total_inner = series(compute_basic_stat(inner, "area_weighted_mean", dim=["lat", "lon"]).load())
        both = pd.concat([total_inner, published], axis=1, join="inner").dropna()
        rmse_inner = float(np.sqrt(np.mean((both.iloc[:, 0] - both.iloc[:, 1]) ** 2)))
        rec.metric(f"{name}: box-edge convention, RMSE with edges included / excluded",
                   f"{np.sqrt(np.mean(d_incl ** 2)):.4f} / {rmse_inner:.4f}", op="info", tier="T3",
                   note=f"{box.sizes['lat']}x{box.sizes['lon']} vs {inner.sizes['lat']}x{inner.sizes['lon']} cells")
        frac = series_fractional_box(sst, lat_range, lon_range)
        both = pd.concat([frac, published], axis=1, join="inner").dropna()
        d = both.iloc[:, 0] - both.iloc[:, 1]
        rec.metric(f"{name}: diagnosis - box mean weighting each cell by its fraction inside the box",
                   f"RMSE {np.sqrt(np.mean(d ** 2)):.4f}, within 0.01: {np.mean(np.abs(d) <= 0.0101):.2f}", op="info",
                   tier="T3", note="the tool gives a cell centred on the box edge full weight; CPC's values match "
                                   "partial (fractional-overlap) weighting far better")
        for file, base in (("ersst5.nino.mth.91-20.ascii", ["1991-01", "2020-12"]),
                           ("ersst5.nino.mth.81-10.ascii", ["1981-01", "2010-12"])):
            index = series(compute_enso_nino_index(sst, nino_region=region, climatology="self", base_period=base))
            agreement(rec, f"{name} ANOM, {base[0][:4]}-{base[1][:4]} base", index, cpc_series(file, f"{name}_ANOM"))


def series_fractional_box(sst, lat_range, lon_range):
    """Box-mean SST with each 2-degree cell weighted by area x the fraction of it inside the box."""
    lat, lon = sst["lat"].values, sst["lon"].values
    lo0, lo1 = lon_range
    lon_s = lon.copy()
    if lo1 < lo0:
        lon_s, lo1 = np.where(lon < 0, lon + 360, lon), lo1 + 360
    fy = np.clip(np.minimum(lat + 1, lat_range[1]) - np.maximum(lat - 1, lat_range[0]), 0, None) / 2
    fx = np.clip(np.minimum(lon_s + 1, lo1) - np.maximum(lon_s - 1, lo0), 0, None) / 2
    w = (np.sin(np.radians(np.clip(lat + 1, -90, 90))) - np.sin(np.radians(np.clip(lat - 1, -90, 90)))) * fy
    weights = w[:, None] * fx[None, :]
    values = sst.values.astype("float64")
    ok = np.isfinite(values) & (weights > 0)
    mean = np.sum(np.where(ok, values * weights, 0), axis=(1, 2)) / np.sum(np.where(ok, weights, 0), axis=(1, 2))
    return pd.Series(mean, index=pd.DatetimeIndex(sst["time"].values))


def oni_base_period(year):
    """CPC's centred 30-year base for the 5-year block holding `year`, capped at the latest complete one."""
    if year <= 1955:
        start = 1936
    else:
        block = 1956 + 5 * ((year - 1956) // 5)
        start = block - 15
    start = min(start, 1996)          # CPC added 1996-2025 at the start of 2026
    return start, start + 29


@case("temporal_aggregate", "ONI: 3-month running mean of Nino3.4, CPC rolling base and fixed 1991-2020 base")
def oni(rec):
    sst = ersst("1936-01", END)
    published_rows = refs.raw_oni(os.path.join(CPC_ROOT, "oni.ascii.txt"))
    centre = {"DJF": 1, "JFM": 2, "FMA": 3, "MAM": 4, "AMJ": 5, "MJJ": 6, "JJA": 7, "JAS": 8, "ASO": 9, "SON": 10,
              "OND": 11, "NDJ": 12}
    published = pd.Series([a for _, _, _, a in published_rows],
                          index=pd.to_datetime([f"{y}-{centre[s]:02d}-01" for s, y, _, _ in published_rows]))

    fixed = compute_enso_nino_index(sst, climatology="self", base_period=["1991-01", "2020-12"])
    oni_fixed = series(temporal_aggregate(fixed, freq="rolling", window=3, center=True))
    agreement(rec, "ONI, fixed 1991-2020 base", oni_fixed["1950":], published, tier="T3", share=None,
              note="report only: early decades sit on a cooler base")

    pieces = []
    for base in sorted({oni_base_period(y) for y in range(1950, int(END[:4]) + 1)}):
        anomaly = compute_enso_nino_index(sst, climatology="self", base_period=[f"{base[0]}-01", f"{base[1]}-12"])
        smoothed = series(temporal_aggregate(anomaly, freq="rolling", window=3, center=True))
        years = [y for y in range(1950, int(END[:4]) + 1) if oni_base_period(y) == base]
        pieces.append(smoothed[smoothed.index.year.isin(years)])
    oni_rolling = pd.concat(pieces).sort_index()
    agreement(rec, "ONI, CPC rolling 30-year bases", oni_rolling, published, tier="T3", share=None,
                  note="report only; bases 1936-65 for 1950-55, then moved every 5 years, latest 1996-2025")
    # Is the published ONI reproducible at all from CPC's own monthly ERSSTv5 Nino3.4?
    total = cpc_series("ersst5.nino.mth.91-20.ascii", "NINO3.4")
    rebuilt = []
    for base in sorted({oni_base_period(y) for y in range(1966, int(END[:4]) + 1)}):
        b = total[str(base[0]):str(base[1])]
        clim = b.groupby(b.index.month).mean()
        a = (total - clim.reindex(total.index.month).values).rolling(3, center=True).mean()
        rebuilt.append(a[[oni_base_period(y) == base and y >= 1966 for y in a.index.year]])
    both = pd.concat([pd.concat(rebuilt).sort_index(), published], axis=1, join="inner").dropna()
    dd = both.iloc[:, 0] - both.iloc[:, 1]
    rec.metric("diagnosis: published ONI rebuilt from CPC's own monthly ERSSTv5 Nino3.4 (no tool involved)",
               f"RMSE {np.sqrt(np.mean(dd ** 2)):.4f}, max {np.abs(dd).max():.3f}", op="info", tier="T3",
               note="1966-2025, documented base rule; a large residual means the ONI file rests on another SST "
                    "input or vintage, so the tool cannot be judged against it")
    both = pd.concat([oni_rolling, published], axis=1, join="inner").dropna()
    by_decade = (both.iloc[:, 0] - both.iloc[:, 1]).abs().groupby((both.index.year // 10) * 10).max()
    rec.metric("ONI rolling base: max |d| by decade", {int(k): round(float(v), 3) for k, v in by_decade.items()},
               op="info", tier="T3")


def runs(values, threshold, above, min_len=5):
    """CPC's episode rule written out directly: (start, end) index pairs of qualifying runs."""
    out, start = [], None
    for i, v in enumerate(list(values) + [np.nan]):
        hit = np.isfinite(v) and (v >= threshold if above else v <= threshold)
        if hit and start is None:
            start = i
        elif not hit and start is not None:
            if i - start >= min_len:
                out.append((start, i - 1))
            start = None
    return out


@case("detect_threshold_events", "ENSO episodes from the published ONI vs CPC's rule and well-known events")
def enso_episodes(rec):
    rows = refs.raw_oni(os.path.join(CPC_ROOT, "oni.ascii.txt"))
    centre = {"DJF": 1, "JFM": 2, "FMA": 3, "MAM": 4, "AMJ": 5, "MJJ": 6, "JJA": 7, "JAS": 8, "ASO": 9, "SON": 10,
              "OND": 11, "NDJ": 12}
    time = pd.to_datetime([f"{y}-{centre[s]:02d}-01" for s, y, _, _ in rows])
    oni = np.round(np.array([a for _, _, _, a in rows]), 1)      # CPC classifies on the one-decimal ONI
    da = xr.DataArray(oni, dims="time", coords={"time": time})
    for label, threshold, direction in (("El Nino", 0.5, "above"), ("La Nina", -0.5, "below")):
        table = detect_threshold_events(da, threshold=threshold, direction=direction, min_duration=5)
        tool = [(pd.Timestamp(s), pd.Timestamp(e)) for s, e in zip(table["start"], table["end"])]
        ref = [(time[a], time[b]) for a, b in runs(oni, threshold, direction == "above")]
        rec.metric(f"{label}: episode list identical to CPC's rule", tool == ref, op="true", tier="T0",
                   note=f"{len(tool)} episodes from the tool, {len(ref)} from the rule")
        known = {"El Nino": ["1972-10", "1982-12", "1997-12", "2015-12", "2023-12"],
                 "La Nina": ["1973-12", "1988-12", "1999-12", "2010-12", "2021-12"]}[label]
        missing = [m for m in known if not any(s <= pd.Timestamp(m) <= e for s, e in tool)]
        rec.metric(f"{label}: well-known events detected ({', '.join(known)})", not missing, op="true", tier="T2",
                   note=f"missing {missing}" if missing else "")


@case("compute_pdo_index", "NCEI PDO (ERSST v6) vs the tool on ERSSTv5")
def ncei_pdo(rec):
    sst = ersst("1900-01", END)
    index, eof, explained, regression = compute_pdo_index(sst)
    tool = series(index)
    raw = refs.raw_pdo(os.path.join(CPC_ROOT, "ersst.v6.pdo.dat"))
    published = pd.Series(list(raw.values()), index=pd.to_datetime([f"{y}-{m:02d}-01" for (y, m) in raw])).sort_index()
    both = pd.concat([tool, published], axis=1, join="inner").dropna()
    a, b = both.iloc[:, 0], both.iloc[:, 1]
    r_month = float(np.corrcoef(a, b)[0, 1])
    rec.metric("monthly correlation, 1900-2025", r_month, threshold=0.8, op=">=", tier="T3",
               note=f"{len(both)} months; explained variance {explained:.1f}%")
    annual = both.groupby(both.index.year).mean()
    rec.metric("annual-mean correlation", float(np.corrcoef(annual.iloc[:, 0], annual.iloc[:, 1])[0, 1]),
               threshold=0.85, op=">=", tier="T3")
    low = both.rolling(61, center=True).mean().dropna()
    rec.metric("5-year running-mean correlation", float(np.corrcoef(low.iloc[:, 0], low.iloc[:, 1])[0, 1]),
               op="info", tier="T3")
    rec.metric("share of months with the same sign", float(np.mean(np.sign(a) == np.sign(b))), op="info", tier="T3")
    za, zb = (a - a.mean()) / a.std(), (b - b.mean()) / b.std()
    rec.metric("RMS difference of the standardised series", float(np.sqrt(np.mean((za - zb) ** 2))), op="info",
               tier="T3")
    phases = [("1925", "1946", 1), ("1947", "1976", -1), ("1977", "1998", 1), ("1999", "2013", -1),
              ("2014", "2016", 1), ("2020", "2025", -1)]
    agree = []
    for start, end, sign in phases:
        mt, mp = float(a[start:end].mean()), float(b[start:end].mean())
        agree.append(np.sign(mt) == np.sign(mp) == sign)
        rec.metric(f"phase {start}-{end}: mean index, tool / NCEI", f"{mt:+.2f} / {mp:+.2f}", op="info", tier="T3")
    rec.metric("every multi-decadal phase has the documented sign in both", all(agree), op="true", tier="T3")


@case("compute_mixed_layer_depth", "WOA23 monthly T/S -> MLD vs the SEANOE dr003 climatology")
def seanoe_mld(rec):
    temp = load_climatology_data("woa23", "temperature", time_slice="monthly").load()
    sal = load_climatology_data("woa23", "salinity", time_slice="monthly").load()
    mld = compute_mixed_layer_depth(temp, sal, ref_depth=10.0, method="density", density_threshold=0.03)
    seanoe = load_climatology_data("seanoe_mld", "mld_dr003", version="v2023").load()
    areas = load_areas(dataset="woa23")
    tool = mld.transpose("month", "lat", "lon").values
    ref = seanoe.transpose("month", "lat", "lon").values
    ok = np.isfinite(tool) & np.isfinite(ref) & (tool > 0) & (ref > 0)
    w = np.broadcast_to(areas.values, tool.shape)[ok]
    lt, lr = np.log(tool[ok]), np.log(ref[ok])
    r = refs.weighted_pearson(lt, lr, w)
    rec.metric("area-weighted pattern correlation of log(MLD), all months", r, threshold=0.7, op=">=", tier="T3",
               note=f"{int(ok.sum())} cell-months; tool missing where WOA's 1500 m monthly column never "
                    f"reaches the threshold: {int((np.isnan(tool) & np.isfinite(ref)).sum())}")
    ratio = np.exp(np.sum(w * (lt - lr)) / np.sum(w))
    rec.metric("area-weighted geometric-mean ratio tool / SEANOE", float(ratio), op="info", tier="T3",
               note="below 1 expected: monthly-mean profiles are smoother than individual ones")
    lat = np.broadcast_to(seanoe["lat"].values[None, :, None], tool.shape)
    month = np.broadcast_to(np.arange(1, 13)[:, None, None], tool.shape)
    for label, sel in (("Northern winter (JFM, 40-70N)", (lat >= 40) & (lat <= 70) & (month <= 3)),
                       ("Southern winter (JAS, 40-70S)", (lat <= -40) & (lat >= -70) & (month >= 7) & (month <= 9)),
                       ("tropics (20S-20N, all months)", np.abs(lat) <= 20)):
        m = ok & sel
        rec.metric(f"{label}: median MLD tool / SEANOE (m)",
                   f"{np.median(tool[m]):.1f} / {np.median(ref[m]):.1f}", op="info", tier="T3")
