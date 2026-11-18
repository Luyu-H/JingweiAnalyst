"""
Harness shared by the tool verification suites.

Each suite module (cat1_analytic ... cat7_safeguards) registers its checks with
@case(tool, name). A case receives a Recorder and calls rec.metric(...) once
per number it measures. The runner executes every case in its own try block,
so one tool that raises never hides the others, and writes one JSON line per
metric to <out_dir>/<category>.jsonl.

Statuses:
  pass / fail  the metric was compared against a threshold
  info         reported without a threshold (T3: magnitude or consistency
               numbers whose acceptable range is a judgement, not a rule)
  error        the case raised before it could measure anything further

Tolerance tiers used throughout (see README.md):
  T0  implementation equivalence      float64 <= 1e-12 relative, float32 <= 1e-6
  T1  discretisation vs analytic      the method's own error bound
  T2  published product reproduction  the product's published precision
  T3  consistency / magnitude         reported, judged by a stated range
"""
from __future__ import annotations

import json
import math
import os
import sys
import time
import traceback
import warnings

import numpy as np
import xarray as xr

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

warnings.filterwarnings("ignore")

CASES = []


def case(tool: str, name: str):
    """Register a verification case for `tool`."""
    def decorator(func):
        CASES.append({"module": func.__module__, "tool": tool, "name": name, "func": func})
        return func
    return decorator


def _plain(value):
    """A JSON-safe version of a metric value."""
    if isinstance(value, xr.DataArray):
        value = value.values
    if isinstance(value, np.ndarray):
        value = value.item() if value.size == 1 else value.tolist()
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    return value


class Recorder:
    """Collects the metrics of one case."""

    def __init__(self, category, tool, case_name, sink):
        self.category, self.tool, self.case_name, self.sink = category, tool, case_name, sink
        self.records = []

    def metric(self, metric, value, *, threshold=None, op="<=", expected=None, tier="", method="",
               note=""):
        """
        Record one measured number.

        op is the comparison that must hold for a pass: '<=' and '>=' against
        `threshold`, '==' for an exact value, 'in' for a (low, high) range,
        'true' for a boolean that must be True, and 'info' to report without
        judging.
        """
        value = _plain(value)
        status = "info"
        if op != "info":
            try:
                v = float(value) if op != "true" else bool(value)
                if op == "<=":
                    ok = v <= threshold
                elif op == ">=":
                    ok = v >= threshold
                elif op == "==":
                    ok = v == threshold
                elif op == "in":
                    ok = threshold[0] <= v <= threshold[1]
                elif op == "true":
                    ok = bool(value)
                else:
                    raise ValueError(f"unknown op {op}")
                status = "pass" if ok else "fail"
            except (TypeError, ValueError):
                status = "fail"
        record = {
            "category": self.category, "tool": self.tool, "case": self.case_name, "metric": metric,
            "value": value, "expected": _plain(expected), "threshold": _plain(threshold), "op": op,
            "status": status, "tier": tier, "method": method, "note": note,
        }
        self.records.append(record)
        self.sink.write(json.dumps(record, ensure_ascii=False) + "\n")
        self.sink.flush()
        return status == "pass" or status == "info"

    def raises(self, metric, func, exc=(ValueError,), match=None, note=""):
        """
        A safeguard check: `func` must raise one of `exc`, and the message must
        contain every string in `match` (case-insensitive).
        """
        match = [match] if isinstance(match, str) else list(match or [])
        exc = exc if isinstance(exc, tuple) else (exc,)
        expected = f"raises {'/'.join(e.__name__ for e in exc)}" + (f" mentioning {match}" if match else "")
        try:
            func()
        except exc as error:
            message = str(error)
            missing = [m for m in match if m.lower() not in message.lower()]
            outcome = f"raised {type(error).__name__}" + (f", message lacks {missing}" if missing else "")
            return self.metric(metric, not missing, op="true", expected=expected, tier="safeguard",
                               note=f"{outcome}: {message[:300]}")
        except Exception as error:     # the wrong kind of failure
            return self.metric(metric, False, op="true", expected=expected, tier="safeguard",
                               note=f"raised {type(error).__name__} instead: {str(error)[:300]}")
        return self.metric(metric, False, op="true", expected=expected, tier="safeguard",
                           note="did not raise" + (f"; {note}" if note else ""))

    def error(self, exc_text):
        record = {"category": self.category, "tool": self.tool, "case": self.case_name,
                  "metric": "case raised", "value": None, "status": "error", "note": exc_text[-3000:]}
        self.records.append(record)
        self.sink.write(json.dumps(record, ensure_ascii=False) + "\n")
        self.sink.flush()


def run_cases(category, module_name, out_dir, only=None):
    """Run every case registered by `module_name`, writing <out_dir>/<category>.jsonl."""
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"{category}.jsonl")
    selected = [c for c in CASES if c["module"] == module_name
                and (only is None or any(o in f"{c['tool']}/{c['name']}" for o in only))]
    totals = {"pass": 0, "fail": 0, "info": 0, "error": 0}
    with open(path, "a") as sink:
        for index, entry in enumerate(selected, start=1):
            rec = Recorder(category, entry["tool"], entry["name"], sink)
            started = time.time()
            try:
                entry["func"](rec)
            except Exception:
                rec.error(traceback.format_exc())
            counts = {k: sum(r["status"] == k for r in rec.records) for k in totals}
            for k in totals:
                totals[k] += counts[k]
            flag = "FAIL" if counts["fail"] or counts["error"] else "ok"
            print(f"[{time.strftime('%H:%M:%S')}] {category} {index}/{len(selected)} {flag:4s} "
                  f"{entry['tool']}/{entry['name']} pass={counts['pass']} fail={counts['fail']} "
                  f"info={counts['info']} error={counts['error']} ({time.time() - started:.1f}s)",
                  flush=True)
    print(f"[{time.strftime('%H:%M:%S')}] {category} DONE {totals}", flush=True)
    return totals


# ---------------------------------------------------------------------------
# Geometry: closed-form areas on the WGS84 ellipsoid
# ---------------------------------------------------------------------------

WGS84_A = 6378137.0
WGS84_F = 1.0 / 298.257223563
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)
WGS84_E = math.sqrt(WGS84_E2)
WGS84_B = WGS84_A * (1.0 - WGS84_F)


def _q(lat_deg):
    s = np.sin(np.radians(np.asarray(lat_deg, dtype="float64")))
    return s / (1.0 - WGS84_E2 * s * s) + np.log((1.0 + WGS84_E * s) / (1.0 - WGS84_E * s)) / (2.0 * WGS84_E)


def zone_area(lat1, lat2, dlon_deg=360.0):
    """
    Area in m2 of the region between two parallels and spanning dlon degrees of
    longitude on the WGS84 ellipsoid (the quadrangle MATLAB's areaquad returns).
    From the authalic-latitude integral: A = b^2/2 * dlon * [q(lat2) - q(lat1)].
    """
    return 0.5 * WGS84_B ** 2 * np.radians(dlon_deg) * (_q(lat2) - _q(lat1))


WGS84_TOTAL_AREA = float(zone_area(-90.0, 90.0))


# ---------------------------------------------------------------------------
# Synthetic fields on the grids the tools really use
# ---------------------------------------------------------------------------

def refill(template: xr.DataArray, values, name=None, **attrs) -> xr.DataArray:
    """
    The template's coordinates, dimensions and attributes with new values.

    The weighting tools find their cell areas, volumes and layer thicknesses
    through the attrs a loader set, so a synthetic field is built on a loaded
    field rather than from scratch: the analytic answer is then tested on the
    grid, levels and provenance the tool meets in practice.
    """
    values = np.broadcast_to(np.asarray(values, dtype="float64"), template.shape).copy()
    out = xr.DataArray(values, coords=template.coords, dims=template.dims,
                       attrs=dict(template.attrs), name=name or template.name)
    out.attrs.update(attrs)
    return out


def lat_lon_like(template: xr.DataArray):
    """Latitude and longitude broadcast to the template's shape (regular or curvilinear)."""
    lat_name = template.attrs.get("lat_coord", "lat") if template.attrs.get("grid_type") == "curvilinear" else "lat"
    lon_name = template.attrs.get("lon_coord", "lon") if template.attrs.get("grid_type") == "curvilinear" else "lon"
    lat, lon, _ = xr.broadcast(template[lat_name], template[lon_name], template)
    return (lat.transpose(*template.dims).values.astype("float64"),
            lon.transpose(*template.dims).values.astype("float64"))


def coord_like(template: xr.DataArray, name: str):
    """A 1-D coordinate broadcast to the template's shape."""
    c, _ = xr.broadcast(template[name], template)
    return c.transpose(*template.dims).values.astype("float64")


def ar1(rng, n, phi, size=(), sigma=1.0):
    """AR(1) noise along the last axis, started from its stationary distribution."""
    shape = tuple(size) + (n,)
    eps = rng.standard_normal(shape) * sigma
    out = np.empty(shape)
    out[..., 0] = eps[..., 0] / math.sqrt(max(1.0 - phi ** 2, 1e-12))
    for t in range(1, n):
        out[..., t] = phi * out[..., t - 1] + eps[..., t]
    return out


def monthly_axis(start="1980-01", periods=480):
    import pandas as pd
    return pd.date_range(start, periods=periods, freq="MS")


def rel_err(a, b):
    a, b = np.asarray(a, dtype="float64"), np.asarray(b, dtype="float64")
    scale = np.maximum(np.abs(b), 1e-300)
    return float(np.nanmax(np.abs(a - b) / scale))


def max_abs(a, b):
    return float(np.nanmax(np.abs(np.asarray(a, dtype="float64") - np.asarray(b, dtype="float64"))))


def same_nan(a, b):
    return bool(np.array_equal(np.isnan(np.asarray(a, dtype="float64")),
                               np.isnan(np.asarray(b, dtype="float64"))))
