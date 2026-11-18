"""
Independent reference implementations the tools are checked against.
"""
from __future__ import annotations

import calendar as _calendar
import math
import os

import numpy as np
from scipy import stats


# ---------------------------------------------------------------------------
# Regression and trend statistics
# ---------------------------------------------------------------------------

def ols_qr(X, y):
    """Ordinary least squares with an intercept, solved by QR."""
    X = np.asarray(X, dtype="float64")
    y = np.asarray(y, dtype="float64")
    if X.ndim == 1:
        X = X[:, None]
    n, k = X.shape
    design = np.column_stack([np.ones(n), X])
    q, r = np.linalg.qr(design)
    beta = np.linalg.solve(r, q.T @ y)
    residuals = y - design @ beta
    dof = n - k - 1
    sigma2 = residuals @ residuals / dof
    r_inv = np.linalg.inv(r)
    covariance = sigma2 * (r_inv @ r_inv.T)
    se = np.sqrt(np.diag(covariance))
    t = beta / se
    p = 2 * stats.t.sf(np.abs(t), dof)
    ss_res = residuals @ residuals
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot
    f = ((ss_tot - ss_res) / k) / (ss_res / dof)
    return {"beta": beta, "se": se, "p": p, "r2": r2, "f": f, "f_p": stats.f.sf(f, k, dof),
            "residuals": residuals, "dof": dof}


def vif_from_correlation(X):
    """Variance inflation factors as diag(inverse correlation matrix)."""
    corr = np.corrcoef(np.asarray(X, dtype="float64"), rowvar=False)
    return np.diag(np.linalg.inv(np.atleast_2d(corr)))


def lag1(x):
    x = np.asarray(x, dtype="float64")
    c = x - x.mean()
    return float(np.sum(c[:-1] * c[1:]) / np.sum(c * c))


def n_eff_single(residuals):
    """AR(1) effective sample size of one series (Santer et al. 2008), r1 clipped to [0, 0.99]."""
    n = len(residuals)
    r1 = min(max(lag1(residuals), 0.0), 0.99)
    return max(n * (1 - r1) / (1 + r1), 3.0), r1


def n_eff_pair(x, y):
    """Effective sample size of a correlation (Bretherton et al. 1999), r1x*r1y clipped to [0, 0.99]."""
    n = len(x)
    rho = min(max(lag1(x) * lag1(y), 0.0), 0.99)
    return min(max(n * (1 - rho) / (1 + rho), 3.0), float(n))


def mann_kendall(y, continuity=True):
    """
    Two-sided Mann-Kendall test: (S, tau, z, p) with the tie-corrected variance.
    continuity=True is the textbook form (pymannkendall's); scipy's kendalltau
    asymptotic p-value, which the tool uses, has no continuity correction.
    """
    y = np.asarray(y, dtype="float64")
    n = y.size
    s = 0.0
    for i in range(n - 1):
        s += np.sum(np.sign(y[i + 1:] - y[i]))
    _, counts = np.unique(y, return_counts=True)
    var_s = (n * (n - 1) * (2 * n + 5) - np.sum(counts * (counts - 1) * (2 * counts + 5))) / 18.0
    correction = np.sign(s) if continuity else 0.0
    z = (s - correction) / math.sqrt(var_s) if s != 0 else 0.0
    p = 2 * stats.norm.sf(abs(z))
    tau = s / (0.5 * n * (n - 1))
    return s, tau, z, p


def theil_sen_slope(x, y):
    x = np.asarray(x, dtype="float64")
    y = np.asarray(y, dtype="float64")
    i, j = np.triu_indices(x.size, k=1)
    dx = x[j] - x[i]
    keep = dx != 0
    return float(np.median((y[j] - y[i])[keep] / dx[keep]))


def aicc(rss, n, k):
    aic = n * math.log(rss / n) + 2 * k
    return aic + 2 * k * (k + 1) / (n - k - 1)


def weighted_pearson(x, y, w):
    x, y, w = (np.asarray(a, dtype="float64").ravel() for a in (x, y, w))
    keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(w)
    x, y, w = x[keep], y[keep], w[keep]
    w = w / w.sum()
    mx, my = np.sum(w * x), np.sum(w * y)
    cov = np.sum(w * (x - mx) * (y - my))
    return cov / math.sqrt(np.sum(w * (x - mx) ** 2) * np.sum(w * (y - my) ** 2))


# ---------------------------------------------------------------------------
# EOF
# ---------------------------------------------------------------------------

def eof_svd(matrix):
    """
    Leading EOF of a (time, space) matrix by SVD, after removing each column's
    time mean. Returns (pc, eof, variance_fraction) with the PC in the same
    units as sklearn's PCA scores (U * S), so the two can be compared directly.
    """
    m = np.asarray(matrix, dtype="float64")
    m = m - m.mean(axis=0, keepdims=True)
    u, s, vt = np.linalg.svd(m, full_matrices=False)
    variance = s ** 2
    return u[:, 0] * s[0], vt[0], float(variance[0] / variance.sum())


# ---------------------------------------------------------------------------
# Calendars
# ---------------------------------------------------------------------------

def days_in_month(year, month, calendar):
    """Month length under the CF calendar rules, written out rather than borrowed from cftime."""
    calendar = calendar.lower()
    if calendar in ("noleap", "365_day"):
        return [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31][month - 1]
    if calendar in ("all_leap", "366_day"):
        return [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31][month - 1]
    if calendar == "360_day":
        return 30
    if calendar == "julian":
        leap = year % 4 == 0
        return [31, 29 if leap else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31][month - 1]
    # standard / gregorian / proleptic_gregorian (identical after 1582)
    return _calendar.monthrange(year, month)[1]


# ---------------------------------------------------------------------------
# Conservative regridding between regular grids
# ---------------------------------------------------------------------------

def _edges(centres):
    c = np.asarray(centres, dtype="float64")
    mid = (c[:-1] + c[1:]) / 2
    return np.concatenate([[c[0] - (c[1] - c[0]) / 2], mid, [c[-1] + (c[-1] - c[-2]) / 2]])


def _overlap_matrix(src_edges, tgt_edges, periodic=False):
    """Length of overlap between every target interval and every source interval."""
    w = np.zeros((tgt_edges.size - 1, src_edges.size - 1))
    shifts = (-360.0, 0.0, 360.0) if periodic else (0.0,)
    for shift in shifts:
        lo = np.maximum(tgt_edges[:-1, None], src_edges[None, :-1] + shift)
        hi = np.minimum(tgt_edges[1:, None], src_edges[None, 1:] + shift)
        w += np.clip(hi - lo, 0.0, None)
    return w


def conservative_regular(field, src_lat, src_lon, tgt_lat, tgt_lon):
    """
    First-order conservative regridding between two regular lat/lon grids on a
    sphere whose cells are bounded by parallels and meridians.
    """
    s_lat_e = np.clip(_edges(src_lat), -90, 90)
    t_lat_e = np.clip(_edges(tgt_lat), -90, 90)
    w_lat = _overlap_matrix(np.sin(np.radians(s_lat_e)), np.sin(np.radians(t_lat_e)))
    w_lon = _overlap_matrix(_edges(src_lon), _edges(tgt_lon), periodic=True)
    f = np.asarray(field, dtype="float64")
    valid = np.isfinite(f)
    f0 = np.where(valid, f, 0.0)
    num = w_lat @ f0 @ w_lon.T
    den = w_lat @ valid.astype(float) @ w_lon.T
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / den, den


def _unit_vectors(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)


def _triangle_excess(a, b, c):
    """Spherical excess of the triangle a-b-c on the unit sphere (Van Oosterom & Strackee 1983)."""
    numerator = np.abs(np.einsum("...i,...i->...", a, np.cross(b, c)))
    denominator = 1 + np.einsum("...i,...i->...", a, b) + np.einsum("...i,...i->...", b, c) \
        + np.einsum("...i,...i->...", c, a)
    return 2 * np.arctan2(numerator, denominator)


def spherical_quad_area(vertex_lat, vertex_lon):
    """
    Area on the unit sphere of quadrilaterals with great-circle edges, the cell
    geometry ESMF assumes.
    """
    v = _unit_vectors(np.asarray(vertex_lat, dtype="float64"), np.asarray(vertex_lon, dtype="float64"))
    return _triangle_excess(v[..., 0, :], v[..., 1, :], v[..., 2, :]) + \
        _triangle_excess(v[..., 0, :], v[..., 2, :], v[..., 3, :])


def regular_cell_vertices(lat_edges, lon_edges):
    """(ny, nx, 4) vertices of a regular grid's cells, anticlockwise from the south-west corner."""
    la0, la1 = np.meshgrid(lat_edges[:-1], lon_edges[:-1], indexing="ij")[0], \
        np.meshgrid(lat_edges[1:], lon_edges[:-1], indexing="ij")[0]
    lo0, lo1 = np.meshgrid(lat_edges[:-1], lon_edges[:-1], indexing="ij")[1], \
        np.meshgrid(lat_edges[:-1], lon_edges[1:], indexing="ij")[1]
    lat = np.stack([la0, la0, la1, la1], axis=-1)
    lon = np.stack([lo0, lo1, lo1, lo0], axis=-1)
    return lat, lon


def duplicate_columns(lat2d, lon2d):
    """Indices of grid columns whose every centre repeats an earlier column's (a cyclic halo)."""
    lat2d, lon2d = np.asarray(lat2d), np.asarray(lon2d)
    if lat2d.ndim == 3:            # (depth, y, x) broadcast: take one level
        lat2d, lon2d = lat2d[0], lon2d[0]
    out = []
    for j in range(lat2d.shape[1]):
        for k in range(j):
            if np.allclose(lat2d[:, j], lat2d[:, k]) and np.allclose(((lon2d[:, j] - lon2d[:, k] + 180) % 360) - 180, 0):
                out.append(j)
                break
    return out


# ---------------------------------------------------------------------------
# Raw readers (no opener in between)
# ---------------------------------------------------------------------------

def raw_nc(path, variable):
    import netCDF4
    with netCDF4.Dataset(path) as ds:
        var = ds[variable]
        var.set_auto_maskandscale(True)
        values = var[:]
        return np.ma.filled(values.astype("float64"), np.nan), dict(var.__dict__)


def raw_cpc_monthly(path):
    """The CPC monthly Nino table as {(year, month): {column: value}}, parsed by hand."""
    out = {}
    with open(path) as handle:
        header = handle.readline().split()
        names, last = [], None
        for column in header:
            if column.startswith("NINO"):
                last = column
                names.append(column)
            elif column == "ANOM":
                names.append(f"{last}_ANOM")
            else:
                names.append(column)
        for line in handle:
            parts = line.split()
            if len(parts) != len(names):
                continue
            row = dict(zip(names, parts))
            out[(int(row["YR"]), int(row["MON"]))] = {k: float(v) for k, v in row.items()
                                                      if k not in ("YR", "MON")}
    return out


def raw_oni(path):
    """The ONI table as a list of (season, year, total, anom)."""
    rows = []
    with open(path) as handle:
        handle.readline()
        for line in handle:
            parts = line.split()
            if len(parts) == 4:
                rows.append((parts[0], int(parts[1]), float(parts[2]), float(parts[3])))
    return rows


def raw_pdo(path):
    """The NCEI PDO table as {(year, month): value}, with the 99.99 sentinels dropped."""
    out = {}
    with open(path) as handle:
        handle.readline()
        handle.readline()
        for line in handle:
            parts = line.split()
            if len(parts) < 2:
                continue
            year = int(parts[0])
            for month, token in enumerate(parts[1:13], start=1):
                value = float(token)
                if value < 99:
                    out[(year, month)] = value
    return out


def raw_woa_mask(path, value_columns):
    """A WOA ASCII mask pivoted on its Latitude/Longitude columns (not filled by index)."""
    import pandas as pd
    table = pd.read_csv(path, comment="#")
    grids = {}
    for column in value_columns:
        pivot = table.pivot_table(index="Latitude", columns="Longitude", values=column, aggfunc="first")
        pivot = pivot.reindex(index=np.arange(-89.5, 90.0, 1.0), columns=np.arange(-179.5, 180.0, 1.0))
        grids[column] = pivot.to_numpy(dtype="float64")
    return grids, table


def file_exists(path):
    return os.path.exists(path)
