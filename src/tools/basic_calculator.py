# tools: compute_basic_stat, temporal_aggregate, fit_1D_trend, transform_unit, 
# fit_multivariate_relation, analyze_correlation, detect_outliers, compute_difference_metrics,
# detect_threshold_events, compute_exceedance_fraction, fit_trend_map, locate_extreme

import re
import xarray as xr
import numpy as np
import pandas as pd
from typing import Dict, Any, List, Union
from scipy import stats
from scipy.stats import f as f_dist
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from sklearn.linear_model import LinearRegression, Ridge, Lasso
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split, GroupShuffleSplit

from src.tools.registry import register_tools
from src.tools.data_loader import areas_like, layer_thickness_like, volumes_like
from src.tools.data_openers import carry_provenance, molar_basis, normalize_units


def _missing_dims_message(data, missing):
    """Why a dimension could not be found, and what to do about it."""
    message = f"Dimension(s) {missing} not found in data. Available: {list(data.dims)}."
    if {"lat", "lon"}.intersection(missing) and "source_dataset" not in data.attrs:
        message += (
            " This array carries no grid attributes, so a curvilinear grid's index dimensions cannot be "
            "addressed as 'lat'/'lon' and its cell areas cannot be located. Name the dimensions the way "
            "this array spells them, or work from a field that kept its provenance — "
            "compute_difference_metrics returns its 'difference' with the grid identity intact."
        )
    return message


def _alias_dims(data, dims):
    """
    Dimension names as the array actually spells them.

    A question is asked about latitude and longitude, but the curvilinear models
    index their grids by nlat/nlon (CESM2) or j/i (MPI-ESM1-2-LR), and their
    lat/lon are 2-D fields over those indices. Writing ['lat', 'lon'] is what
    anyone means by "over the horizontal", so it is translated here rather than
    rejected.
    """
    alias = {"lat": data.attrs.get("y_dim"), "lon": data.attrs.get("x_dim")}
    return [d if d in data.dims else (alias.get(d) or d) for d in dims]


def _horizontal_names(data):
    """Every name this array's horizontal dimensions answer to."""
    return {data.attrs.get("y_dim") or "lat", data.attrs.get("x_dim") or "lon", "lat", "lon"}


def _resolve_weighted_dims(data, dim, required, op):
    """
    Which dimensions a weighted mean should reduce over.

    Accepts the geographic names on the curvilinear models, whose grids are
    indexed by nlat/nlon or j/i, because a question is asked about latitude and
    longitude rather than about grid indices. Refuses a reduction that leaves
    out an axis the weights vary along: averaging over depth alone with area
    weights, or over longitude alone, would apply a weight that does not belong
    to the axis being collapsed.
    """
    if dim is None:
        return None
    requested = [dim] if isinstance(dim, str) else list(dim)
    alias = {"lat": data.attrs.get("y_dim"), "lon": data.attrs.get("x_dim")}
    reduce_dims = [d if d in data.dims else alias.get(d, d) for d in requested]

    missing = [d for d in reduce_dims if d not in data.dims]
    if missing:
        raise ValueError(f"Dimension(s) {missing} not found in data. Available: {list(data.dims)}.")
    if not set(required).issubset(reduce_dims):
        raise ValueError(
            f"{op} has to reduce over {list(required)}, but was asked for {reduce_dims}. "
            f"Weights along an axis that is not being collapsed do not apply; use op='mean' for those."
        )
    return reduce_dims


@register_tools(
    name="compute_basic_stat",
    description="Reduce a field along chosen dimensions with one basic statistic.\n"
        "Use a WEIGHTED mean for any physical average. A plain mean is only right when every element "
        "represents the same amount of ocean, which is rarely true: cells shrink towards the poles, and "
        "layer thicknesses vary by more than an order of magnitude with depth. The weights are taken from "
        "the grid the data came from, so the field must still be on its original grid and levels — "
        "subsetting is fine, regridding or averaging over depth beforehand is not.",
    args={
        "data": "xr.DataArray, required. Input field.",
        "op": "str, required. The statistic to apply:\n"
            "- 'mean', 'median', 'min', 'max', 'sum', 'count' (non-missing elements), 'std' (standard "
            "deviation about the mean, ddof=0).\n"
            "- 'area_weighted_mean' for a horizontal average, weighting each cell by its area.\n"
            "- 'volume_weighted_mean' whenever depth is among the dimensions averaged, weighting each cell "
            "by its volume. This is the correct operator for the mean temperature or concentration of a "
            "depth range.\n"
            "- 'argmax'/'argmin' return WHERE the extreme falls rather than what it is: the coordinate "
            "label at which it occurs, one per remaining column. They reduce a single dimension. For the "
            "one extreme of a whole field, use locate_extreme.",
        "dim": "str or list of str, optional. Dimension(s) to reduce. If None, every dimension is reduced. "
            "A weighted mean must include the dimensions its weights vary along: the horizontal pair for "
            "'area_weighted_mean', and 'depth' for 'volume_weighted_mean'. On a curvilinear grid you may "
            "write ['lat', 'lon'] and they are mapped to the model's index dimensions."
    },
    returns="xr.DataArray with the reduced dimension(s) removed. For 'argmax'/'argmin' the values are "
        "coordinate labels rather than data values. A weighted result records in attrs which weights were "
        "used."
)
def compute_basic_stat(data: xr.DataArray, op: str, dim: str | list[str] | None = None) -> xr.DataArray:
    op = op.lower()
    if op == "count":
        result = data.count(dim=dim)
    elif op == "max":
        result = data.max(dim=dim)
    elif op == "min":
        result = data.min(dim=dim)
    elif op == "mean":
        result = data.mean(dim=dim)
    elif op == "median":
        result = data.median(dim=dim)
    elif op == "sum":
        result = data.sum(dim=dim)
    elif op == "std":
        result = data.std(dim=dim)

    elif op in ("argmax", "argmin"):
        # idxmax/idxmin give the coordinate LABEL of the extreme, which is the
        # answer to "which month" or "which year"; argmax on the bare array
        # would give a position, which means nothing without the axis.
        requested = _alias_dims(data, [dim] if isinstance(dim, str) else list(dim or data.dims))
        if len(requested) != 1:
            raise ValueError(
                f"'{op}' finds where the extreme falls along ONE axis, so it needs a single dimension; it "
                f"was asked for {requested}. Reduce the other dimensions first (for instance with "
                f"op='area_weighted_mean'), then ask which step of what is left holds the extreme."
            )
        axis = requested[0]
        if axis not in data.dims:
            raise ValueError(f"Dimension '{axis}' not found in data. Available: {list(data.dims)}.")
        if axis not in data.coords:
            raise ValueError(
                f"'{axis}' has no coordinate labels, so there is nothing for '{op}' to name. Use a "
                f"dimension with coordinates, such as 'time', 'month' or 'depth'."
            )
        result = data.idxmax(dim=axis) if op == "argmax" else data.idxmin(dim=axis)
        result.attrs = dict(data.attrs)
        result.attrs["long_name"] = f"{axis} of the {'maximum' if op == 'argmax' else 'minimum'}"
        result.attrs["units"] = str(data[axis].attrs.get("units", axis))
        return result

    elif op in ("area_weighted_mean", "volume_weighted_mean"):
        # The weights come from the grid the data itself records, so ERSSTv5 is
        # weighted on its 2-degree cells and a curvilinear model with its own
        # areacello. Weighting every product with the 1-degree grid instead
        # leaves ERSSTv5 with no cell in common and silently empties the result.
        if op == "area_weighted_mean":
            full_weights = areas_like(data)
            required = tuple(d for d in full_weights.dims if d in data.dims)
            source = full_weights.attrs.get("area_source", "cell areas")
        else:
            if "depth" not in data.dims:
                raise ValueError(
                    "volume_weighted_mean needs a 'depth' dimension. For a purely horizontal average use "
                    "op='area_weighted_mean'."
                )
            full_weights = volumes_like(data)
            required = ("depth",)
            source = full_weights.attrs.get(
                "volume_source", full_weights.attrs.get("area_source", "cell volumes")
            )

        reduce_dims = _resolve_weighted_dims(data, dim, required, op)

        # Land and below-sea-floor cells carry area but no value, and the models
        # disagree on whether areacello is masked at all (CESM2's is not), so the
        # weights are masked with the field itself before they are summed.
        weights = full_weights.where(data.notnull())
        result = (data * weights).sum(dim=reduce_dims, skipna=True) / weights.sum(dim=reduce_dims)
        result.attrs = dict(data.attrs)
        result.attrs["weighting"] = f"{'area' if op == 'area_weighted_mean' else 'volume'}: {source}"

    else:
        raise ValueError(
            f"Unsupported operation '{op}'. Choose from ['count', 'max', 'min', 'mean', "
            f"'area_weighted_mean', 'volume_weighted_mean', 'median', 'sum', 'std', 'argmax', 'argmin']."
        )

    return result


# The two season conventions in use here. They are NOT interchangeable: the same
# label covers a different three months, so a composite built with one and an
# anomaly built with the other are offset by a month.
#   meteorological — DJF/MAM/JJA/SON, the climate-science default
#   woa            — calendar quarters, which is how WOA23's seasonal
#                    climatology files are defined, so compute_anomaly_vs_woa23(freq='season')
#                    necessarily uses this one
SEASON_CONVENTIONS = {
    "meteorological": {"resample": "QS-DEC",
                       "labels": {12: "DJF", 3: "MAM", 6: "JJA", 9: "SON"}},
    "woa": {"resample": "QS-JAN",
            "labels": {1: "JFM", 4: "AMJ", 7: "JAS", 10: "OND"}},
}


@register_tools(
    name="temporal_aggregate",
    description="Aggregate a time series to annual, seasonal or climatological-monthly values, or smooth "
        "it with a running window.\n"
        "Means are weighted by the number of days in each month, using the dataset's own calendar: "
        "calendar months are 28 to 31 days long, so averaging them with equal weight is not a time "
        "average. Weighting applies to op='mean' only; a median, minimum, maximum or standard deviation of "
        "unequal-length months has no weighted equivalent.\n"
        "Only a field with a real time axis can be aggregated; a climatology is already a long-term mean.",
    args={
        "data": "xr.DataArray, required. A field with a 'time' dimension.",
        "freq": "str, optional, default 'year'. What to aggregate to:\n"
            "- 'year': one value per calendar year.\n"
            "- 'season': one value per seasonal period (see season_convention).\n"
            "- 'climatology': one value per calendar month, across all years.\n"
            "- 'rolling': a running window that SMOOTHS the series instead of collapsing it — the time "
            "axis keeps every step, and each value is the statistic over the `window` steps around it.",
        "op": "str, optional, default 'mean'. 'mean' (day-weighted), 'median', 'min', 'max' or 'std'.",
        "season_convention": "str, optional, default 'meteorological'. Which three months each season "
            "covers: 'meteorological' gives DJF/MAM/JJA/SON, the usual climate definition; 'woa' gives the "
            "calendar quarters JFM/AMJ/JAS/OND used by some observational climatologies. The two are "
            "offset by a month and must not be mixed.",
        "window": "int, required for freq='rolling', ignored otherwise. Number of time steps in the "
            "window, in the units of the time axis. An odd window is centred exactly; an even one leans "
            "one step towards the past.",
        "center": "bool, optional, default True, freq='rolling' only. Place each window's value at its "
            "middle step, keeping the smoothed series in phase with the original. False places it at the "
            "last step, which is what an index computed from data already in hand does, and shifts every "
            "feature later by half a window."
    },
    returns="xr.DataArray, aggregated in time.\n"
        "- 'year': the time axis becomes an integer 'year' coordinate. attrs['incomplete_periods'] lists "
        "any year not covered by a full 12 months, which biases a trend if left in.\n"
        "- 'season': one entry per seasonal period, with 'season' and 'year' coordinates attached; a "
        "December-starting season is labelled by the year holding its January.\n"
        "- 'climatology': a 'month' coordinate (1-12) with month names attached.\n"
        "- 'rolling': the time axis is unchanged in length and labels, with the ends left missing where no "
        "complete window covers them. The running mean is unweighted, because a running window is a filter "
        "over a fixed number of steps rather than a time average."
)
def temporal_aggregate(data: xr.DataArray, freq: str = 'year', op: str = 'mean',
                       season_convention: str = 'meteorological', window: int = None,
                       center: bool = True) -> xr.DataArray:
    if 'time' not in data.dims:
        # WOA23 and SEANOE are climatologies: they have no real time axis, and
        # carry 'month' or 'season' instead. There is nothing to aggregate there.
        for climatology_dim in ('month', 'season'):
            if climatology_dim in data.dims:
                raise ValueError(
                    f"This is a climatology with a '{climatology_dim}' dimension, not a time series, so it "
                    f"cannot be aggregated in time: it is already an average over many years. Select the "
                    f"slices you want with load_climatology_data(months=...) instead, or load a dataset with "
                    f"a real time axis ('cmip6_historical', 'cmip6_ssp245', 'ersst') if you need a trend."
                )
        raise ValueError(f"Input data must have a 'time' dimension. Available: {list(data.dims)}")

    freq = freq.lower()
    op = op.lower()
    if op not in ('mean', 'median', 'min', 'max', 'std'):
        raise ValueError(f"Operation '{op}' not supported. Choose from: ['mean', 'median', 'min', 'max', 'std']")
    if season_convention not in SEASON_CONVENTIONS:
        raise ValueError(
            f"season_convention must be one of {sorted(SEASON_CONVENTIONS)}, got '{season_convention}'"
        )

    # Month length in the dataset's own calendar, so a 365-day model is weighted
    # on its own terms rather than on the Gregorian calendar's.
    days = data['time'].dt.days_in_month.astype('float64')

    def aggregate(group):
        if op != 'mean':
            return getattr(group, op)(dim='time', skipna=True)
        weights = days.sel(time=group['time'])
        # Masking the weights with the field keeps a month that is missing at one
        # cell from counting towards that cell's denominator.
        valid = weights.where(group.notnull())
        return (group * valid).sum(dim='time', skipna=True) / valid.sum(dim='time')

    if freq == 'rolling':
        if window is None:
            raise ValueError(
                "freq='rolling' needs a window: the number of time steps to average over. The ONI is a "
                "3-month running mean of monthly data, i.e. window=3."
            )
        window = int(window)
        if window < 2:
            raise ValueError(f"window must be at least 2 time steps, got {window}.")
        if window > data.sizes['time']:
            raise ValueError(
                f"A window of {window} steps does not fit in a series of {data.sizes['time']} steps."
            )
        # min_periods=window leaves the ends NaN rather than averaging a partial
        # window: a 3-month mean built from 2 months is a different statistic,
        # and CPC's own series simply starts one step in.
        rolling = data.rolling(time=window, center=center, min_periods=window)
        result = getattr(rolling, op)(skipna=True)
        result.attrs = data.attrs.copy()
        result.attrs['temporal_aggregation'] = f'rolling{window}_{op}'
        result.attrs['aggregation_operation'] = op
        result.attrs['rolling_window'] = f'{window} time steps, {"centred" if center else "trailing"}'
        result.attrs['time_weighting'] = 'unweighted (a running window is a filter over a fixed number of steps)'
        return result

    if freq == 'year':
        result = data.resample(time='1YE').map(aggregate)
        months_per_year = data['time'].dt.year.to_series().value_counts()
        incomplete = sorted(int(y) for y, n in months_per_year.items() if n != 12)
        result['time'] = result['time'].dt.year
        result = result.rename({'time': 'year'})

    elif freq == 'season':
        convention = SEASON_CONVENTIONS[season_convention]
        result = data.resample(time=convention['resample']).map(aggregate)
        labels = convention['labels']
        seasons = [labels[m] for m in result['time'].dt.month.values]
        years = result['time'].dt.year.values
        if season_convention == 'meteorological':
            # A DJF season starts in December but is named for the year holding
            # its January and February, which is the usual convention.
            years = np.where(result['time'].dt.month.values == 12, years + 1, years)
        result = result.assign_coords(season=('time', seasons), year=('time', years))
        incomplete = None

    elif freq == 'climatology':
        result = data.groupby('time.month').map(aggregate)
        month_names = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                       'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        result = result.assign_coords(month_name=('month', month_names))
        incomplete = None

    else:
        raise ValueError(
            f"Frequency '{freq}' not supported. Choose from: 'year', 'season', 'climatology', 'rolling'"
        )

    result.attrs = data.attrs.copy()
    result.attrs['temporal_aggregation'] = f'{freq}_{op}'
    result.attrs['aggregation_operation'] = op
    result.attrs['time_weighting'] = (
        'weighted by days in month' if op == 'mean' else f'unweighted ({op} has no weighted form)'
    )
    if freq == 'season':
        result.attrs['season_convention'] = (
            f"{season_convention} ({'/'.join(SEASON_CONVENTIONS[season_convention]['labels'].values())})"
        )
    if freq == 'year' and incomplete:
        result.attrs['incomplete_periods'] = (
            f"years not covered by a full 12 months: {incomplete}; their means rest on fewer months and "
            f"should usually be dropped before fitting a trend"
        )
    return result


def effective_sample_size(residuals: np.ndarray) -> tuple:
    """
    Sample size and lag-1 autocorrelation adjusted for serial correlation.

    A monthly ocean series is nowhere near independent from one step to the
    next, so the nominal sample size badly overstates the information it holds.
    Treating the residuals as AR(1) gives n_eff = n (1 - r1) / (1 + r1), the
    standard correction for trend confidence intervals in climate work
    (Santer et al. 2008, J. Geophys. Res.). With r1 = 0.9 that is a twentieth of
    the nominal count, and an uncorrected p-value is meaningless at that point.
    """
    n = residuals.size
    if n < 3:
        return float(n), 0.0
    centred = residuals - residuals.mean()
    denominator = float(np.sum(centred ** 2))
    if denominator <= 0:
        return float(n), 0.0
    r1 = float(np.sum(centred[:-1] * centred[1:]) / denominator)
    r1 = min(max(r1, 0.0), 0.99)        # negative autocorrelation is not used to inflate n
    n_eff = n * (1.0 - r1) / (1.0 + r1)
    return max(float(n_eff), 3.0), r1


def _aicc(residual_sum_squares: float, n: int, k: int) -> float:
    """
    Small-sample corrected Akaike information criterion.

    Model choice cannot be made on R2: adding a parameter can never lower it, so
    a quadratic always beats a straight line and a four-parameter sinusoid beats
    both. AICc scores fit against the number of parameters spent on it, so a
    more complex model has to earn its place.
    """
    if residual_sum_squares <= 0 or n <= k + 2:
        return np.inf
    aic = n * np.log(residual_sum_squares / n) + 2 * k
    return aic + (2 * k * (k + 1)) / (n - k - 1)


@register_tools(
    name="fit_1D_trend",
    description="Fit a trend to one series and report its magnitude with its uncertainty.\n"
        "For the usual question — is this quantity changing, and how fast — use method='linear' or "
        "method='theil_sen'. Both return the slope with a confidence interval and a p-value, and the "
        "interval is what should be quoted: a slope without one cannot be called significant.\n"
        "The significance accounts for SERIAL CORRELATION. A monthly ocean series is strongly "
        "autocorrelated, so the number of steps overstates the independent information by a large factor; "
        "the effective sample size is estimated from the lag-1 autocorrelation of the residuals and used "
        "for the interval and the p-value. The unadjusted values are returned alongside.\n"
        "On a date axis the fit runs in years, so a linear slope is a rate per YEAR whatever the sampling "
        "frequency; multiply by 10 to report it per decade. For a trend at every cell of a field, use "
        "fit_trend_map.",
    args={
        "data": "xr.DataArray, required. A strictly one-dimensional series. A 'time' coordinate is used as "
            "the x axis in years; any other numeric coordinate is used in its own units.",
        "method": "str, optional. The model to fit:\n"
            "'linear' (y = a*x + b) — the standard choice for a trend;\n"
            "'theil_sen' (median-of-slopes) — the same question, resistant to outliers and to non-normal "
            "residuals;\n"
            "'quadratic' — for an accelerating change;\n"
            "'exponential' — only attempted when every value is positive;\n"
            "'sinusoidal' — for a cyclic signal.\n"
            "If None, every applicable model is fitted and the best by AICc is returned. Note that "
            "automatic selection answers 'which curve describes this series', NOT 'what is the trend': a "
            "seasonal series will legitimately select the sinusoid. Ask for 'linear' or 'theil_sen' "
            "explicitly when the question is about a trend.",
        "alpha": "float, optional, default 0.05. Significance level, so the default gives a 95 percent "
            "confidence interval."
    },
    returns="dict:\n"
        "'trend' (for 'linear' and 'theil_sen'): {slope, slope_se, ci95_low, ci95_high, p_value, "
        "p_value_unadjusted, n, n_effective, lag1_autocorr, units, significant}. The slope is per x unit "
        "(per YEAR on a date axis). Quote it with its interval, and call it significant only if the "
        "interval excludes zero.\n"
        "'metrics': {rmse, mae, r2, adj_r2, aicc, n}. adj_r2 and aicc, unlike r2, are comparable between "
        "models with different numbers of parameters.\n"
        "'fitted_data': the fitted values on the input's coordinates.\n"
        "'equation': the fitted model written out, stating what x runs in.\n"
        "'params'/'param_std': fitted parameters and their standard errors.\n"
        "'x_unit', 'x_origin': the unit and zero point of the fitted axis.\n"
        "'method_used' and 'all_methods_aicc' when the model was selected automatically."
)
def fit_1D_trend(data: xr.DataArray, method: str = None, alpha: float = 0.05) -> Dict[str, Any]:
    if data.ndim != 1:
        raise ValueError("Input DataArray must be 1D.")

    method = method.lower() if method is not None else None

    # The slope has to carry a unit, and the array index is not one: on a
    # monthly axis a per-index slope is a per-month rate, which is what makes a
    # trend read a factor of 12 (or 120 for a decade) too small. Where the
    # dimension is a real date axis the fit is done in years since the first
    # sample, so the slope is per year whatever the sampling frequency.
    coord_name = data.dims[0]
    coord = data.coords.get(coord_name)
    if coord is not None and np.issubdtype(coord.dtype, np.datetime64):
        elapsed_days = (coord.values - coord.values[0]) / np.timedelta64(1, "D")
        x_all = elapsed_days / 365.2425          # mean Gregorian year
        x_unit = "year"
        x_origin = str(np.datetime_as_string(coord.values[0], unit="D"))
    elif coord is not None and np.issubdtype(coord.dtype, np.number):
        # An integer 'year' axis from temporal_aggregate, or depth in metres.
        x_all = coord.values.astype(float) - float(coord.values[0])
        x_unit = str(coord.attrs.get("units", coord_name))
        x_origin = str(coord.values[0])
    else:
        x_all = np.arange(len(data), dtype=float)
        x_unit = f"{coord_name} index"
        x_origin = "0"

    y_all = data.values.astype(float)

    # Fit on the finite points only: a single gap otherwise propagates through
    # every candidate model and leaves nothing to choose from.
    finite = np.isfinite(x_all) & np.isfinite(y_all)
    x, y = x_all[finite], y_all[finite]
    if x.size < 3:
        raise ValueError(
            f"Only {x.size} finite values in this series, which is too few to fit a trend. "
            f"The input has {y_all.size} points, so most of it is missing."
        )

    def linear(x, a, b):
        return a * x + b
    def quadratic(x, a, b, c):
        return a * x**2 + b * x + c
    def exponential(x, a, b):
        return a * np.exp(b * x)
    def sinusoidal(x, a, b, c, d):
        return a * np.sin(b * x + c) + d

    models = {
        "linear": (linear, [0.0, float(np.mean(y))], ["a", "b"]),
        "quadratic": (quadratic, [0.0, 0.0, float(np.mean(y))], ["a", "b", "c"]),
        # One cycle per year, expressed in whatever x runs in: 2*pi per year on a
        # date axis, 2*pi/12 on a bare monthly index.
        "sinusoidal": (sinusoidal,
                       [float(np.std(y)), 2 * np.pi if x_unit == "year" else 2 * np.pi / 12,
                        0.0, float(np.mean(y))],
                       ["a", "b", "c", "d"]),
    }
    # Fitting an exponential to a series that changes sign is meaningless, and
    # clamping the negatives to zero first would fit one curve and then score it
    # against different data.
    if np.all(y > 0):
        models["exponential"] = (exponential, [float(np.mean(y)), 0.01], ["a", "b"])

    results = {}
    for name, (func, p0, param_names) in models.items():
        try:
            params, covariance = curve_fit(func, x, y, p0=p0, maxfev=10000)
            y_fit = func(x, *params)
            residuals = y - y_fit
            rss = float(np.sum(residuals ** 2))
            k = len(params)
            with np.errstate(invalid="ignore"):
                param_std = np.sqrt(np.diag(covariance))
            results[name] = {
                "params": params, "param_names": param_names, "param_std": param_std,
                "fitted": y_fit, "residuals": residuals,
                "metrics": {
                    "rmse": float(np.sqrt(mean_squared_error(y, y_fit))),
                    "mae": float(mean_absolute_error(y, y_fit)),
                    "r2": float(r2_score(y, y_fit)),
                    "adj_r2": float(1 - (1 - r2_score(y, y_fit)) * (x.size - 1) / max(x.size - k - 1, 1)),
                    "aicc": _aicc(rss, x.size, k),
                    "n": int(x.size),
                },
            }
        except Exception:
            continue

    # Theil-Sen is not a curve_fit model; it is fitted directly when asked for,
    # and is never chosen automatically because it answers the same question as
    # 'linear' rather than describing a different shape.
    if method == "theil_sen":
        slope, intercept, lo_slope, hi_slope = stats.theilslopes(y, x, alpha=1 - alpha)
        y_fit = slope * x + intercept
        residuals = y - y_fit
        results["theil_sen"] = {
            "params": np.array([slope, intercept]), "param_names": ["a", "b"],
            "param_std": np.array([np.nan, np.nan]),
            "fitted": y_fit, "residuals": residuals,
            "theil_sen_ci": (float(lo_slope), float(hi_slope)),
            "metrics": {
                "rmse": float(np.sqrt(mean_squared_error(y, y_fit))),
                "mae": float(mean_absolute_error(y, y_fit)),
                "r2": float(r2_score(y, y_fit)),
                "adj_r2": float(1 - (1 - r2_score(y, y_fit)) * (x.size - 1) / max(x.size - 3, 1)),
                "aicc": _aicc(float(np.sum(residuals ** 2)), x.size, 2),
                "n": int(x.size),
            },
        }

    if not results:
        raise ValueError("No model could be fitted to this series; check that it is not constant.")

    all_aicc = {name: result["metrics"]["aicc"] for name, result in results.items()}
    if method is None:
        method = min(results, key=lambda m: results[m]["metrics"]["aicc"])
    if method not in results:
        raise ValueError(
            f"Method '{method}' could not be fitted. Available: {sorted(results)}. "
            f"('exponential' is only offered when every value is positive.)"
        )
    chosen = results[method]

    # Fitted values go back onto the full axis, with the gaps left as gaps.
    fitted_full = np.full(y_all.shape, np.nan)
    fitted_full[finite] = chosen["fitted"]
    fitted_data = xr.DataArray(
        fitted_full, coords=data.coords, dims=data.dims,
        attrs={**data.attrs, "fitted_by": method},
    )

    params = {name: float(value) for name, value in zip(chosen["param_names"], chosen["params"])}
    param_std = {name: float(value) for name, value in zip(chosen["param_names"], chosen["param_std"])}

    result = {
        "fitted_data": fitted_data,
        "metrics": chosen["metrics"],
        "params": params,
        "param_std": param_std,
        "x_unit": x_unit,
        "x_origin": x_origin,
        "method_used": method,
        "all_methods_aicc": all_aicc,
    }

    # A straight line is the only shape whose slope IS the trend, so that is
    # where the confidence interval and the significance test belong.
    if method in ("linear", "theil_sen"):
        slope = float(chosen["params"][0])
        residuals = chosen["residuals"]
        n = int(x.size)
        n_eff, lag1 = effective_sample_size(residuals)

        spread = float(np.sum((x - x.mean()) ** 2))
        rss = float(np.sum(residuals ** 2))
        se_nominal = np.sqrt(rss / max(n - 2, 1) / spread) if spread > 0 else np.nan
        # The same standard error with the residual degrees of freedom reduced to
        # the effective sample size, which is what serial correlation leaves.
        se_effective = np.sqrt(rss / max(n_eff - 2, 1) / spread) if spread > 0 else np.nan

        dof = max(n_eff - 2, 1)
        critical = stats.t.ppf(1 - alpha / 2, dof)
        p_adjusted = float(2 * (1 - stats.t.cdf(abs(slope / se_effective), dof))) if se_effective > 0 else np.nan
        p_nominal = float(
            2 * (1 - stats.t.cdf(abs(slope / se_nominal), max(n - 2, 1)))
        ) if se_nominal > 0 else np.nan

        if method == "theil_sen":
            ci_low, ci_high = chosen["theil_sen_ci"]
            # Kendall's tau against time is the Mann-Kendall trend test.
            tau, p_mk = stats.kendalltau(x, y)
            p_adjusted = float(p_mk)
            result["mann_kendall_tau"] = float(tau)
            # Unlike the least-squares interval above, neither the Theil-Sen
            # interval nor the Mann-Kendall p-value knows about serial
            # correlation, so both are optimistic on a monthly series. Say so
            # rather than let the two methods appear to disagree.
            result["autocorrelation_warning"] = (
                f"the Theil-Sen interval and the Mann-Kendall p-value assume independent samples, but the "
                f"lag-1 autocorrelation is {lag1:.2f}, leaving an effective sample size of {n_eff:.0f} of "
                f"{n}. Both are therefore optimistic; method='linear' returns an interval corrected for this."
            )
        else:
            ci_low = slope - critical * se_effective
            ci_high = slope + critical * se_effective

        result["trend"] = {
            "slope": slope,
            "slope_se": float(se_effective),
            "slope_se_unadjusted": float(se_nominal),
            "ci95_low": float(ci_low),
            "ci95_high": float(ci_high),
            "p_value": p_adjusted,
            "p_value_unadjusted": p_nominal,
            "n": n,
            "n_effective": float(n_eff),
            "lag1_autocorr": float(lag1),
            "units": f"{data.attrs.get('units', 'value')} per {x_unit}",
            "significant": bool(ci_low > 0 or ci_high < 0),
        }
        result["equation"] = (
            f"{method}: slope = {slope:.6g} [{ci_low:.6g}, {ci_high:.6g}] "
            f"{data.attrs.get('units', 'value')} per {x_unit} (95% CI, n_eff={n_eff:.0f} of {n}), "
            f"p = {p_adjusted:.3g}; x is {x_unit} since {x_origin}"
        )
    else:
        param_str = ", ".join(f"{value:.6g}" for value in chosen["params"])
        result["equation"] = (
            f"{method}, fit with params: [{param_str}], where x is {x_unit} since {x_origin}"
        )

    return result


def _density_in_kg_per_litre(density):
    """
    A density argument as kg/L, whichever of the three forms it arrives in.

    A field can be handed over either as the anomaly sigma that
    compute_potential_density_anomaly returns (about 25 kg/m3) or as an absolute
    density (about 1025 kg/m3). Telling them apart matters by a factor of forty,
    so it is decided from the array's own metadata first and from its magnitude
    only as a fallback, and the decision is reported rather than assumed
    silently.
    """
    if not isinstance(density, xr.DataArray):
        value = float(density)
        if not 0.9 < value < 1.2:
            raise ValueError(
                f"A constant density is read as kg/L and seawater is about 1.025; {value} is outside the "
                f"plausible range. Pass 1.025 rather than 1025, or hand over a density field."
            )
        return value, f"constant seawater density {value} kg/L ({value * 1000:g} kg/m3)"

    long_name = str(density.attrs.get("long_name", "")).lower()
    median = float(density.median(skipna=True))
    if "anomaly" in long_name or "sigma" in long_name:
        is_anomaly = True
        reason = f"its long_name, '{density.attrs.get('long_name')}'"
    elif np.isfinite(median):
        is_anomaly = median < 500.0
        reason = f"its median value of {median:.1f} kg/m3"
    else:
        raise ValueError("The density field is entirely missing; nothing can be converted with it.")

    absolute = density + 1000.0 if is_anomaly else density
    note = (
        f"density field, read as {'an anomaly (sigma, + 1000)' if is_anomaly else 'absolute density'} "
        f"from {reason}; {float(absolute.min()):.1f} to {float(absolute.max()):.1f} kg/m3 cell by cell "
        f"({density.attrs.get('reference_pressure', 'reference pressure not recorded')})"
    )
    return absolute / 1000.0, note


@register_tools(
    name="transform_unit",
    description="Convert a field from one unit to another within the same physical category.\n"
        "The conversion that matters most is between molar bases: some products report concentrations per "
        "unit MASS (umol kg-1) and others per unit VOLUME (mol m-3), and the two cannot be compared, "
        "differenced or evaluated against each other until one has been converted.\n"
        "That conversion needs seawater density. The default is a constant, which is the mean surface "
        "value and errs by about a percent in the deep ocean and in fresh polar water — enough to rival "
        "the differences such a comparison sets out to measure. Pass a density field for anything below "
        "the surface; compute_potential_density_anomaly with ref_pressure=None produces one.",
    args={
        "data": "xr.DataArray, required. The field to convert.",
        "from_unit": "str, required, case-insensitive. The unit the data is in; passing the array's own "
            "'units' attribute works, including the variant spellings different products use. Recognised "
            "categories: temperature ('degC', 'K'), concentration per mass ('umol/kg', 'mmol/kg'), "
            "concentration per volume ('mol m-3', 'mmol/m3', 'umol/L'), dimensionless ('1', 'psu', ''), "
            "length ('m', 'cm'), velocity ('m/s', 'cm/s', 'mm/s').",
        "to_unit": "str, required. The target unit, in any of the same spellings. It must belong to the "
            "same category as from_unit.",
        "density": "float or xr.DataArray, optional, default 1.025. Seawater density, used only between "
            "the mass and volume bases. A float is read as kg/L. A DataArray is aligned cell by cell, so "
            "each cell converts at its own density; both an absolute density and a density anomaly (sigma) "
            "are understood, decided from the array's metadata and magnitude and reported back in the "
            "result."
    },
    returns="xr.DataArray with converted values and its 'units' attribute updated. attrs record "
        "'original_units' and, for a mass/volume conversion, 'density_assumption' — whether a constant or "
        "a field was used, and how it was interpreted. Quote that line when converted concentrations are "
        "reported."
)
def transform_unit(data: xr.DataArray, from_unit: str, to_unit: str, density=1.025) -> xr.DataArray:
    std_from = normalize_units(from_unit)
    std_to = normalize_units(to_unit)
    density, density_note = _density_in_kg_per_litre(density)

    category_of = {
        "degC": "temperature", "K": "temperature",
        "umol kg-1": "concentration", "mmol kg-1": "concentration",
        "mol m-3": "concentration", "mmol m-3": "concentration",
        "1": "dimensionless",
        "m": "length", "cm": "length",
        "m s-1": "velocity", "cm s-1": "velocity", "mm s-1": "velocity",
    }
    for label, unit, raw in (("from_unit", std_from, from_unit), ("to_unit", std_to, to_unit)):
        if unit not in category_of:
            raise ValueError(
                f"Unsupported {label}: '{raw}'. Supported units are {sorted(category_of)} and their common "
                f"spellings (e.g. 'micromoles_per_kilogram', 'micro-mol kg-1', 'degrees_celsius', 'psu')."
            )
    if category_of[std_from] != category_of[std_to]:
        raise ValueError(
            f"'{from_unit}' is a {category_of[std_from]} unit and '{to_unit}' is a "
            f"{category_of[std_to]} unit; there is no conversion between the two."
        )

    # Factors onto each category's base unit. Concentration uses mmol m-3 as its
    # base, which is the same number as umol L-1, and the mass-basis entries are
    # the only ones that depend on density: 1 umol kg-1 x 1025 kg m-3 is
    # 1025 umol m-3, i.e. 1.025 mmol m-3.
    to_base = {
        "mmol m-3": 1.0, "mol m-3": 1000.0,
        "umol kg-1": density, "mmol kg-1": 1000.0 * density,
        "m": 1.0, "cm": 0.01,
        "m s-1": 1.0, "cm s-1": 0.01, "mm s-1": 0.001,
        "1": 1.0,
    }

    if category_of[std_from] == "temperature":
        celsius = data - 273.15 if std_from == "K" else data * 1
        converted = celsius + 273.15 if std_to == "K" else celsius
    else:
        converted = data * (to_base[std_from] / to_base[std_to])

    result = converted
    result.attrs = dict(data.attrs)
    result.attrs["units"] = std_to
    result.attrs["original_units"] = from_unit
    if molar_basis(std_from) != molar_basis(std_to):
        approximate = isinstance(density, float)
        result.attrs["density_assumption"] = density_note + (
            "; the in-situ density was not used, so this conversion is approximate"
            if approximate else "; converted at each cell's own density"
        )
    result.name = data.name
    return result


def spatial_block_labels(template: xr.DataArray, block_degrees: float = 10.0):
    """
    A block label for every cell, so that validation can hold out whole regions.

    Neighbouring ocean cells are near-duplicates of one another, so a random
    train/test split puts a cell's own neighbours in the training set and the
    validation score then measures interpolation between adjacent cells rather
    than the model's ability to generalise. Grouping cells into blocks of a few
    degrees and holding out whole blocks is the standard remedy for spatially
    autocorrelated data. Returns None when the field has no lat/lon to group by.
    """
    lat_name = template.attrs.get("lat_coord", "lat")
    lon_name = template.attrs.get("lon_coord", "lon")
    if lat_name not in template.coords or lon_name not in template.coords:
        return None
    lat, lon = xr.broadcast(template[lat_name], template[lon_name])
    lat, _ = xr.broadcast(lat, template)
    lon, _ = xr.broadcast(lon, template)
    lat_block = np.floor(lat.transpose(*template.dims).values / block_degrees)
    lon_block = np.floor(lon.transpose(*template.dims).values / block_degrees)
    return (lat_block * 10000 + lon_block).ravel()


def variance_inflation_factors(X: np.ndarray, names) -> dict:
    """
    VIF for each predictor: how much collinearity inflates its coefficient's variance.

    Ocean tracers are strongly collinear by construction — nitrate and phosphate
    track each other at close to the Redfield ratio — and collinear predictors
    give coefficients that are individually meaningless and unstable in sign
    even while the fit as a whole is excellent. Above about 5 the individual
    coefficients should not be interpreted; above 10 the collinearity is severe.
    """
    vif = {}
    for j, name in enumerate(names):
        others = np.delete(X, j, axis=1)
        if others.shape[1] == 0:
            vif[name] = 1.0
            continue
        try:
            model = LinearRegression().fit(others, X[:, j])
            r2 = r2_score(X[:, j], model.predict(others))
            vif[name] = float(1.0 / (1.0 - r2)) if r2 < 0.999999 else float("inf")
        except Exception:
            vif[name] = float("nan")
    return vif


def ols_inference(X: np.ndarray, y: np.ndarray, names) -> dict:
    """
    Full-sample ordinary-least-squares inference: per-coefficient standard
    errors, t statistics and p-values, plus the overall F test.

    These are only valid for an unpenalised least-squares fit. Ridge and lasso
    deliberately bias their coefficients, so a t or F test computed from them
    does not mean what it appears to, and none is reported for those.
    """
    n, k = X.shape
    design = np.column_stack([np.ones(n), X])
    coefficients, *_ = np.linalg.lstsq(design, y, rcond=None)
    residuals = y - design @ coefficients
    dof = n - k - 1
    if dof <= 0:
        return {}
    sigma2 = float(residuals @ residuals) / dof
    try:
        covariance = sigma2 * np.linalg.inv(design.T @ design)
    except np.linalg.LinAlgError:
        return {}
    standard_errors = np.sqrt(np.diag(covariance))
    t_values = coefficients / standard_errors
    p_values = 2 * (1 - stats.t.cdf(np.abs(t_values), dof))

    labels = ["intercept"] + list(names)
    ss_res = float(residuals @ residuals)
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    inference = {
        "coefficient_se": {l: float(v) for l, v in zip(labels, standard_errors)},
        "coefficient_p_value": {l: float(v) for l, v in zip(labels, p_values)},
    }
    if ss_tot > 0 and k >= 1:
        f_stat = ((ss_tot - ss_res) / k) / (ss_res / dof)
        inference["f_statistic"] = float(f_stat)
        inference["f_p_value"] = float(1 - f_dist.cdf(f_stat, k, dof))
    return inference


@register_tools(
    name="fit_multivariate_relation",
    description="Fit a statistical model predicting one variable from one or more others, and report the "
        "coefficients with diagnostics.\n"
        "Use it for a relationship between variables rather than a change over time: the stoichiometric "
        "slope between two nutrients, the dependence of a tracer on temperature and salinity, an empirical "
        "predictor built from several fields. With a single predictor it is an ordinary regression, "
        "returning the slope, the intercept and R2.\n"
        "Validation holds out whole spatial blocks by default. Neighbouring ocean cells are near-"
        "duplicates, so a random split scores a model against its own neighbours and reports a skill it "
        "does not have.\n"
        "Ocean tracers are often strongly collinear, which makes individual coefficients unreliable even "
        "when the fit is good; the variance inflation factors in the diagnostics say when that is "
        "happening.",
    args={
        "target": "xr.DataArray, required. The variable to be predicted.",
        "predictors": "list[xr.DataArray], required. One or more predictor fields, on the same "
            "coordinates as the target.",
        "method": "str, optional. 'linear', 'ridge', 'lasso', 'polynomial' or 'random_forest'. If None, "
            "several are fitted and the best by validation score is returned.",
        "alpha": "float, optional, default 1.0. Regularisation strength for ridge and lasso.",
        "poly_degree": "int, optional, default 2. Polynomial degree for method='polynomial'.",
        "test_size": "float, optional, default 0.2. Fraction of the data held out for validation.",
        "standardize": "bool, optional. Whether to standardise predictors. Defaults to True for ridge and "
            "lasso, whose penalties are scale-dependent.",
        "validation": "str, optional, default 'blocked'. 'blocked' holds out whole spatial regions; "
            "'random' splits cells at random, which is only appropriate for data with no spatial "
            "structure.",
        "block_degrees": "float, optional, default 10.0. Size of the held-out blocks in degrees, which "
            "should exceed the correlation length scale of the field."
    },
    returns="dict:\n"
        "'coefficients': the fitted coefficients, including the intercept, for the linear family; None for "
        "tree-based models. These are the raw coefficients and are NOT comparable between predictors.\n"
        "'metrics': {rmse, mae, r2, adj_r2, train_r2, test_r2}.\n"
        "'diagnostics': residual statistics; 'vif' per predictor with a warning when collinearity is "
        "severe; how the held-out set was chosen; and, for method='linear', coefficient standard errors, "
        "p-values and the F test.\n"
        "'feature_importance': comparable importance per predictor (standardised coefficients for the "
        "linear family, impurity importances for a forest).\n"
        "'fitted_data': predictions on the target's coordinates.\n"
        "'equation', 'method_used', and 'all_methods_scores' when the model was selected automatically."
)
def fit_multivariate_relation(target: xr.DataArray, predictors: List[xr.DataArray], method: str = None,
                              alpha: float = 1.0, poly_degree: int = 2, test_size: float = 0.2,
                              standardize: bool = None, validation: str = "blocked",
                              block_degrees: float = 10.0) -> Dict[str, Any]:
    if not predictors:
        raise ValueError("At least one predictor variable is required.")
    
    for p in predictors:
        if not p.coords.equals(target.coords):
            raise ValueError("Predictor coordinates do not match target.")


    predictor_names = []
    for i, pred in enumerate(predictors):
        if pred.name:
            predictor_names.append(str(pred.name))
        elif 'long_name' in pred.attrs:
            predictor_names.append(pred.attrs['long_name'])
        else:
            predictor_names.append(f"X{i+1}")
    
    y = target.values.flatten()
    X_list = [pred.values.flatten() for pred in predictors]
    X = np.column_stack(X_list)
    
    valid_mask = ~(np.isnan(y) | np.any(np.isnan(X), axis=1))
    X_clean = X[valid_mask]
    y_clean = y[valid_mask]
    
    n_samples = len(y_clean)
    n_features = X_clean.shape[1]
    
    if n_samples < n_features + 1:
        raise ValueError(f"Insufficient data: {n_samples} samples for {n_features} features. "
                        f"Need at least {n_features + 1} valid samples.")
    
    # Validation split. Whole spatial blocks are held out where the field has
    # lat/lon, because a random split over a gridded field validates against a
    # cell's own neighbours and reports a score that is far too optimistic.
    groups = spatial_block_labels(target, block_degrees)
    split_kind = "random"
    if groups is not None and validation == "blocked":
        groups_clean = groups[valid_mask]
        n_blocks = len(np.unique(groups_clean))
        if n_blocks >= 5:
            splitter = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=42)
            train_idx, test_idx = next(splitter.split(X_clean, y_clean, groups=groups_clean))
            X_train, X_test = X_clean[train_idx], X_clean[test_idx]
            y_train, y_test = y_clean[train_idx], y_clean[test_idx]
            split_kind = f"blocked by {block_degrees}-degree region ({n_blocks} blocks)"
        else:
            X_train, X_test, y_train, y_test = train_test_split(
                X_clean, y_clean, test_size=test_size, random_state=42)
            split_kind = f"random (only {n_blocks} spatial blocks, too few to hold any out)"
    else:
        X_train, X_test, y_train, y_test = train_test_split(
            X_clean, y_clean, test_size=test_size, random_state=42)
        if validation == "blocked":
            split_kind = "random (no lat/lon coordinates to block by)"
    
    def fit_single_method(method_name, alpha_val=alpha, standardize_val=None):
        if standardize_val is None:
            standardize_val = method_name in ["ridge", "lasso"]
        scaler = StandardScaler() if standardize_val else None
        if scaler:
            X_train_scaled = scaler.fit_transform(X_train)
            X_test_scaled = scaler.transform(X_test)
            X_clean_scaled = scaler.transform(X_clean)
        else:
            X_train_scaled = X_train
            X_test_scaled = X_test
            X_clean_scaled = X_clean
        
        if method_name == "linear":
            model = LinearRegression()
            model.fit(X_train_scaled, y_train)
            y_pred = model.predict(X_clean_scaled)
            y_train_pred = model.predict(X_train_scaled)
            y_test_pred = model.predict(X_test_scaled)
            extra_data = {"model": model, "scaler": scaler}
            
        elif method_name == "ridge":
            model = Ridge(alpha=alpha_val)
            model.fit(X_train_scaled, y_train)
            y_pred = model.predict(X_clean_scaled)
            y_train_pred = model.predict(X_train_scaled)
            y_test_pred = model.predict(X_test_scaled)
            extra_data = {"model": model, "scaler": scaler}
            
        elif method_name == "lasso":
            model = Lasso(alpha=alpha_val, max_iter=5000)
            model.fit(X_train_scaled, y_train)
            y_pred = model.predict(X_clean_scaled)
            y_train_pred = model.predict(X_train_scaled)
            y_test_pred = model.predict(X_test_scaled)
            extra_data = {"model": model, "scaler": scaler}
            
        elif method_name == "polynomial":
            poly = PolynomialFeatures(degree=poly_degree, include_bias=False)
            X_train_poly = poly.fit_transform(X_train_scaled)
            X_test_poly = poly.transform(X_test_scaled)
            X_clean_poly = poly.transform(X_clean_scaled)
            
            model = LinearRegression()
            model.fit(X_train_poly, y_train)
            
            y_pred = model.predict(X_clean_poly)
            y_train_pred = model.predict(X_train_poly)
            y_test_pred = model.predict(X_test_poly)
            
            poly_features = poly.get_feature_names_out(predictor_names)
            extra_data = {"model": model, "scaler": scaler, "poly": poly, "poly_features": poly_features}
            
        elif method_name == "random_forest":
            model = RandomForestRegressor(n_estimators=100, max_depth=10, random_state=42, n_jobs=-1)
            model.fit(X_train_scaled, y_train)
            y_pred = model.predict(X_clean_scaled)
            y_train_pred = model.predict(X_train_scaled)
            y_test_pred = model.predict(X_test_scaled)
            extra_data = {"model": model, "scaler": scaler}
            
        else:
            raise ValueError(f"Unknown method '{method_name}'. Options: linear, ridge, lasso, polynomial, random_forest")
        
        test_r2 = float(r2_score(y_test, y_test_pred))
        
        return {
            "y_pred": y_pred,
            "y_train_pred": y_train_pred,
            "y_test_pred": y_test_pred,
            "test_r2": test_r2,
            "extra_data": extra_data
        }
    
    all_methods_scores = {}
    if method is None:
        methods_to_try = ["linear", "ridge", "lasso", "polynomial", "random_forest"]
        best_method = None
        best_score = -np.inf
        best_result = None
        
        for m in methods_to_try:
            try:
                result = fit_single_method(m)
                all_methods_scores[m] = result["test_r2"]
                if result["test_r2"] > best_score:
                    best_score = result["test_r2"]
                    best_method = m
                    best_result = result
            except Exception as e:
                all_methods_scores[m] = float('nan')
                continue
        
        if best_method is None:
            raise ValueError("All methods failed. Check the input data.")
        
        method = best_method
        fit_result = best_result

    else:
        fit_result = fit_single_method(method, standardize_val=standardize)
        all_methods_scores = {method: fit_result["test_r2"]}
    
    y_pred = fit_result["y_pred"]
    y_train_pred = fit_result["y_train_pred"]
    y_test_pred = fit_result["y_test_pred"]
    extra_data = fit_result["extra_data"]
    model = extra_data["model"]
    scaler = extra_data.get("scaler")
    
    rmse = float(np.sqrt(mean_squared_error(y_clean, y_pred)))
    mae = float(mean_absolute_error(y_clean, y_pred))
    r2 = float(r2_score(y_clean, y_pred))
    train_r2 = float(r2_score(y_train, y_train_pred))
    test_r2 = float(r2_score(y_test, y_test_pred))
    
    metrics = {
        "rmse": rmse,
        "mae": mae,
        "r2": r2,
        "train_r2": train_r2,
        "test_r2": test_r2,
        "n_samples": int(n_samples),
        "n_features": int(n_features)
    }
    
    coefficients = {}
    feature_importance = {}
    
    if method in ["linear", "ridge", "lasso"]:
        coefs = model.coef_
        intercept = model.intercept_
        
        coefficients["intercept"] = float(intercept)
        for i, name in enumerate(predictor_names):
            coefficients[name] = float(coefs[i])
        
        terms = [f"{intercept:.4f}"]
        for i, name in enumerate(predictor_names):
            coef = coefs[i]
            if abs(coef) > 1e-10:  # Skip near-zero terms, especially for Lasso
                sign = "+" if coef >= 0 else ""
                terms.append(f"{sign}{coef:.4f}*{name}")
        equation = f"y = {' '.join(terms)}"
        
        # A raw coefficient is not an importance: it carries the units of its
        # predictor, so phosphate at a sixteenth of nitrate's magnitude gets a
        # coefficient sixteen times larger for the same physical effect.
        # Standardising by each predictor's spread makes them comparable.
        y_spread = float(np.std(y_clean)) or 1.0
        for i, name in enumerate(predictor_names):
            feature_importance[name] = float(abs(coefs[i]) * np.std(X_clean[:, i]) / y_spread)

    elif method == "polynomial":
        poly_features = extra_data["poly_features"]
        coefs = model.coef_
        intercept = model.intercept_
        
        coefficients["intercept"] = float(intercept)
        equation_terms = [f"{intercept:.4f}"]
        
        for i, feat_name in enumerate(poly_features):
            coef = coefs[i]
            coefficients[feat_name] = float(coef)
            if abs(coef) > 1e-10:
                sign = "+" if coef >= 0 else ""
                equation_terms.append(f"{sign}{coef:.4f}*{feat_name}")
        
        if len(equation_terms) > 10:
            equation = f"y = {' '.join(equation_terms[:10])}... (degree={poly_degree}, {len(poly_features)} terms)"
        else:
            equation = f"y = {' '.join(equation_terms)}"
        
        for name in predictor_names:
            importance = sum(abs(coefs[i]) for i, fn in enumerate(poly_features) if name in fn)
            feature_importance[name] = float(importance)
    
    elif method == "random_forest":
        importances = model.feature_importances_
        for i, name in enumerate(predictor_names):
            feature_importance[name] = float(importances[i])
        
        sorted_features = sorted(feature_importance.items(), key=lambda x: x[1], reverse=True)
        top_features = ", ".join([f"{name}({imp:.3f})" for name, imp in sorted_features[:5]])
        equation = f"Random Forest ensemble (n_estimators=100, top features: {top_features})"
    
    residuals = y_clean - y_pred
    diagnostics = {
        "residual_mean": float(np.mean(residuals)),
        "residual_std": float(np.std(residuals)),
        "residual_skewness": float(np.mean((residuals / np.std(residuals))**3) if np.std(residuals) > 0 else 0),
        "max_residual": float(np.max(np.abs(residuals))),
        "validation_split": split_kind,
        "vif": variance_inflation_factors(X_clean, predictor_names),
    }
    max_vif = max((v for v in diagnostics["vif"].values() if np.isfinite(v)), default=0.0)
    if max_vif >= 5:
        diagnostics["collinearity_warning"] = (
            f"highest VIF is {max_vif:.1f}. The predictors carry much the same information, so the "
            f"individual coefficients are unstable and should not be read as separate physical effects, "
            f"even though the fit as a whole may be good."
        )

    # Significance tests belong to an unpenalised least-squares fit only, and to
    # one fitted on the data it is being tested against, so they come from a
    # full-sample OLS rather than from the 80 percent training fit above.
    if method == "linear":
        diagnostics.update(ols_inference(X_clean, y_clean, predictor_names))
    elif method in ("ridge", "lasso"):
        diagnostics["inference_note"] = (
            f"no t or F test is reported for {method}: penalised regression biases its coefficients "
            f"deliberately, so the usual tests do not apply. Refit with method='linear' for inference."
        )
    
    y_pred_full = np.full_like(y, np.nan)
    y_pred_full[valid_mask] = y_pred
    y_pred_reshaped = y_pred_full.reshape(target.shape)

    fitted_data = xr.DataArray(
        y_pred_reshaped,
        coords=target.coords,
        dims=target.dims,
        attrs={
            **target.attrs,
            "fitted_method": method,
            "predictors": ", ".join(predictor_names),
            "r2": r2
        }
    )
    
    result = {
        "fitted_data": fitted_data,
        "equation": equation,
        "metrics": metrics,
        "coefficients": coefficients,
        "feature_importance": feature_importance,
        "diagnostics": diagnostics,
        "method_used": method
    }
    
    if all_methods_scores:
        result["all_methods_scores"] = all_methods_scores
    
    return result


def _area_weighted_correlation(var1, var2, dims, method, lag, detrend, min_periods):
    """
    Pattern correlation between two maps, weighting each cell by its area.

    An unweighted correlation over lat/lon treats every cell as one observation,
    and a 1-degree grid has as many cells in the 5 degrees around a pole as in
    5 degrees of the tropics while they cover a fortieth of the ocean. The
    weighted form is the standard pattern correlation of model evaluation.

    The p-value here is the naive one: neighbouring ocean cells are not
    independent, so the number of cells overstates the evidence just as the
    number of months does along time, and there is no single accepted correction
    for a 2-D field. Read the coefficient, not the p-value.
    """
    if method != "pearson":
        raise ValueError(
            "weighting='area' applies to the Pearson correlation only. Spearman correlates ranks, which "
            "carry no cell areas with them."
        )
    if lag != 0 or detrend:
        raise ValueError(
            "weighting='area' is for a correlation between maps, where a time lag and detrending have no "
            "meaning. Drop lag/detrend, or use weighting='none' for a correlation along time."
        )

    if not set(dims).issubset(_horizontal_names(var1)):
        raise ValueError(
            f"weighting='area' weights cells by their area, so the correlation has to be taken over the "
            f"horizontal dimensions; this one is over {dims}. Use weighting='none' for a correlation along "
            f"time or depth."
        )

    # Land, and any cell missing from either field, must not carry weight: the
    # two fields are compared over the cells they share and nothing else.
    valid = var1.notnull() & var2.notnull()
    weights = areas_like(var1).where(valid)
    total = weights.sum(dim=dims)

    def weighted_mean(field):
        return (field * weights).sum(dim=dims, skipna=True) / total

    mean1, mean2 = weighted_mean(var1), weighted_mean(var2)
    anomaly1, anomaly2 = var1 - mean1, var2 - mean2
    covariance = weighted_mean(anomaly1 * anomaly2)
    variance1, variance2 = weighted_mean(anomaly1 ** 2), weighted_mean(anomaly2 ** 2)
    correlation = covariance / np.sqrt(variance1 * variance2)

    n_samples = valid.sum(dim=dims)
    correlation = correlation.where(n_samples >= min_periods)

    def two_tailed_p(r, n):
        r = np.asarray(r, dtype="float64")
        n = np.asarray(n, dtype="float64")
        with np.errstate(divide="ignore", invalid="ignore"):
            t_stat = r * np.sqrt(np.maximum(n - 2, 1) / (1 - r ** 2))
        p = 2 * (1 - stats.t.cdf(np.abs(t_stat), np.maximum(n - 2, 1)))
        return np.where(np.isfinite(r), p, np.nan)

    pvalue = xr.apply_ufunc(two_tailed_p, correlation, n_samples, dask="parallelized",
                            output_dtypes=[float])

    def wrap(array, name, attrs):
        result = array.rename(name)
        result.attrs = attrs
        return carry_provenance(result, var1)

    area_note = weights.attrs.get("area_source", "cell areas on the fields' own grid")
    return {
        "correlation": wrap(correlation, f"{var1.name}_vs_{var2.name}_correlation", {
            "long_name": "Area-weighted pattern correlation (pearson)",
            "method": "pearson",
            "correlation_dims": list(dims),
            "weighting": f"area: {area_note}",
            "lag": 0,
            "detrended": "False",
        }),
        "pvalue": wrap(pvalue, "pvalue", {
            "long_name": "Two-tailed p-value",
            "description": "assumes independent cells; neighbouring ocean cells are not, so this is optimistic",
        }),
        "pvalue_unadjusted": wrap(pvalue, "pvalue_unadjusted", {
            "long_name": "Two-tailed p-value assuming independent samples",
            "description": "identical to 'pvalue' here: no correction exists for spatial autocorrelation",
        }),
        "n_samples": wrap(n_samples, "n_samples", {"long_name": "Number of cells common to both fields"}),
        "n_effective": wrap(n_samples.astype("float64"), "n_effective", {
            "long_name": "Effective number of independent samples",
            "description": "not estimated for a spatial correlation; equal to the cell count",
        }),
    }


@register_tools(
    name="analyze_correlation",
    description="Correlate two fields, with optional lead-lag analysis and significance testing.\n"
        "Along time the p-value is corrected for SERIAL CORRELATION: monthly ocean series are strongly "
        "autocorrelated, so an uncorrected test uses far too many degrees of freedom and declares weak "
        "correlations highly significant. The uncorrected value is returned alongside so the difference is "
        "visible.\n"
        "Set detrend=True when both series contain a long-term trend, since two independently rising "
        "series correlate through their trends alone.\n"
        "Across space, use weighting='area' for a pattern correlation between two maps: cells shrink "
        "towards the poles, and an unweighted correlation counts a polar cell as heavily as a tropical one "
        "many times its size.",
    args={
        "var1": "xr.DataArray, required. The first field.",
        "var2": "xr.DataArray, required. The second field, aligned with the first.",
        "dim": "str or list[str], optional, default 'time'. Dimension(s) to correlate over: 'time' for a "
            "temporal correlation at each cell, the horizontal pair for a pattern correlation between two "
            "maps. On a curvilinear grid you may write ['lat', 'lon'].",
        "method": "str, optional, default 'pearson'. 'pearson' for a linear relationship, 'spearman' to "
            "correlate ranks, which detects a monotonic but non-linear relationship and resists outliers.",
        "lag": "int, optional, default 0. Time lag in steps for lead-lag analysis. Positive means var2 "
            "leads var1; negative means var1 leads var2. Requires 'time' among the correlation dimensions.",
        "detrend": "bool, optional, default False. Remove a least-squares straight line from both series "
            "before correlating, so a shared trend cannot manufacture the correlation.",
        "min_periods": "int, optional, default 10. Minimum number of valid pairs required; correlations "
            "from a handful of points are unstable.",
        "weighting": "str, optional, default 'none'. 'area' weights each cell by its area, for a pattern "
            "correlation. It applies only over the horizontal dimensions, only to 'pearson', and not with "
            "a lag or with detrending; the weights come from the grid the fields are on, so put them on a "
            "common grid first."
    },
    returns="dict of xr.DataArray, each shaped by the dimensions not correlated over:\n"
        "'correlation': the coefficient, -1 to 1.\n"
        "'pvalue': two-tailed. Along time this is corrected for serial correlation and is the one to "
        "quote; across space no correction exists, so read the coefficient rather than the p-value.\n"
        "'pvalue_unadjusted': the naive value assuming every sample is independent.\n"
        "'n_samples': valid data pairs.\n"
        "'n_effective': independent samples after accounting for autocorrelation. When it is much smaller "
        "than n_samples, the correlation rests on that many independent observations."
)
def analyze_correlation(var1: xr.DataArray, var2: xr.DataArray, dim: Union[str, List[str]] = "time",
                        method: str = "pearson", lag: int = 0, detrend: bool = False,
                        min_periods: int = 10, weighting: str = "none") -> Dict[str, xr.DataArray]:
    if method not in ["pearson", "spearman"]:
        raise ValueError(f"Method must be 'pearson' or 'spearman', got '{method}'")
    if weighting not in ("none", "area"):
        raise ValueError(f"weighting must be 'none' or 'area', got '{weighting}'")
    
    var1_aligned, var2_aligned = xr.align(var1, var2, join='inner')
    
    if dim is None:
        dims_to_correlate = list(var1_aligned.dims)
    elif isinstance(dim, str):
        dims_to_correlate = [dim]
    else:
        dims_to_correlate = list(dim)
    
    # ['lat', 'lon'] on a curvilinear model means its index dimensions.
    dims_to_correlate = _alias_dims(var1_aligned, dims_to_correlate)

    for d in dims_to_correlate:
        if d not in var1_aligned.dims:
            raise ValueError(f"Dimension '{d}' not found in var1. Available: {var1_aligned.dims}")
        if d not in var2_aligned.dims:
            raise ValueError(f"Dimension '{d}' not found in var2. Available: {var2_aligned.dims}")
    
    if weighting == "area":
        return _area_weighted_correlation(var1_aligned, var2_aligned, dims_to_correlate,
                                          method, lag, detrend, min_periods)

    if lag != 0:
        if 'time' not in dims_to_correlate:
            raise ValueError("Lag can only be applied when 'time' is in correlation dimensions")

        time_len = len(var1_aligned['time'])
        if abs(lag) >= time_len:
            raise ValueError(f"Lag {lag} exceeds available time dimension length {time_len}")

        if lag > 0:
            # var2 leads var1: shift var2 backward (earlier times)
            var2_aligned = var2_aligned.isel(time=slice(None, -lag))
            var1_aligned = var1_aligned.isel(time=slice(lag, None))
        else:
            # var1 leads var2: shift var1 backward
            var1_aligned = var1_aligned.isel(time=slice(None, -abs(lag)))
            var2_aligned = var2_aligned.isel(time=slice(abs(lag), None))

        # The two windows are now paired by POSITION: by construction their time
        # labels differ by `lag` steps, and xarray's exact join inside apply_ufunc
        # below rejects arrays whose labels disagree. Put var2 on var1's time axis
        # (labels do not survive into a correlation anyway) so the pairing stands.
        stale = [c for c in var2_aligned.coords if "time" in var2_aligned[c].dims]
        var2_aligned = var2_aligned.drop_vars(stale)
        var2_aligned = var2_aligned.assign_coords(
            {c: var1_aligned[c] for c in var1_aligned.coords if "time" in var1_aligned[c].dims}
        )

    # The serial-correlation correction is only meaningful along an ordered axis
    # such as time; across a flattened spatial block the lag-1 statistic has no
    # single interpretation, so the sample size is left alone there.
    correct_serial = dims_to_correlate == ["time"]

    def correlation_func(x, y):
        """Correlate one pair of profiles, handling NaN and serial correlation."""
        x_flat = np.asarray(x, dtype="float64").ravel()
        y_flat = np.asarray(y, dtype="float64").ravel()
        valid_mask = ~(np.isnan(x_flat) | np.isnan(y_flat))
        x_valid = x_flat[valid_mask]
        y_valid = y_flat[valid_mask]

        n = len(x_valid)
        if n < max(min_periods, 3):
            return np.nan, np.nan, np.nan, n, np.nan

        if detrend:
            index = np.arange(n, dtype="float64")
            x_valid = x_valid - np.polyval(np.polyfit(index, x_valid, 1), index)
            y_valid = y_valid - np.polyval(np.polyfit(index, y_valid, 1), index)

        if np.std(x_valid) == 0 or np.std(y_valid) == 0:
            return np.nan, np.nan, np.nan, n, np.nan

        try:
            if method == "pearson":
                corr, p_naive = stats.pearsonr(x_valid, y_valid)
            else:
                corr, p_naive = stats.spearmanr(x_valid, y_valid)
        except Exception:
            return np.nan, np.nan, np.nan, n, np.nan

        if not correct_serial or abs(corr) >= 1.0:
            return float(corr), float(p_naive), float(p_naive), int(n), float(n)

        # Effective sample size from the lag-1 autocorrelation of both series.
        # Two series that each persist from month to month share far fewer
        # independent samples than they have time steps.
        def lag1(values):
            centred = values - values.mean()
            denominator = float(np.sum(centred ** 2))
            if denominator <= 0:
                return 0.0
            return float(np.sum(centred[:-1] * centred[1:]) / denominator)

        r1x, r1y = lag1(x_valid), lag1(y_valid)
        product = min(max(r1x * r1y, 0.0), 0.99)
        n_eff = min(max(n * (1 - product) / (1 + product), 3.0), float(n))

        t_stat = corr * np.sqrt(max(n_eff - 2, 1) / (1 - corr ** 2))
        p_adjusted = float(2 * (1 - stats.t.cdf(abs(t_stat), max(n_eff - 2, 1))))
        return float(corr), p_adjusted, float(p_naive), int(n), float(n_eff)

    # CMIP6 arrays arrive chunked along time (12 months per chunk), and a core
    # dimension of apply_ufunc cannot span multiple dask chunks, so the
    # correlation dimensions are gathered into one chunk first. Everything else
    # stays chunked, so a global field is still processed block by block.
    def gather_core_dims(array):
        if array.chunks is None:
            return array
        # unify_chunks first: a field's 2-D coordinates are chunked differently
        # from the field itself, and asking for chunksizes across the two raises
        # rather than answering.
        array = array.unify_chunks()
        split = [d for d in dims_to_correlate if len(array.chunksizes.get(d, ())) > 1]
        return array.chunk({d: -1 for d in split}) if split else array

    var1_aligned = gather_core_dims(var1_aligned)
    var2_aligned = gather_core_dims(var2_aligned)

    corr_values, pval_values, pval_unadjusted, n_samples, n_effective = xr.apply_ufunc(
        correlation_func,
        var1_aligned,
        var2_aligned,
        input_core_dims=[dims_to_correlate, dims_to_correlate],
        output_core_dims=[[], [], [], [], []],
        vectorize=True,
        dask='parallelized',
        output_dtypes=[float, float, float, int, float]
    )

    remaining_dims = [d for d in var1_aligned.dims if d not in dims_to_correlate]
    # Every coordinate that still fits the remaining axes, not just the ones
    # named after them: a curvilinear grid carries its latitude and longitude as
    # 2-D fields, and dropping them leaves a map that cannot be located.
    remaining_coords = {
        name: coord for name, coord in var1_aligned.coords.items()
        if set(coord.dims).issubset(remaining_dims)
    }

    def wrap(values, name, attrs):
        result = xr.DataArray(values, coords=remaining_coords, dims=remaining_dims, name=name, attrs=attrs)
        # A correlation map lives on the same grid as its inputs, so it keeps
        # their grid identity and can be area-weighted like any other field.
        return carry_provenance(result, var1_aligned)

    correlation = wrap(
        corr_values, f"{var1.name}_vs_{var2.name}_correlation",
        {
            "long_name": f"Correlation coefficient ({method})",
            "method": method,
            "correlation_dims": dims_to_correlate,
            "lag": lag,
            "detrended": str(detrend),
        },
    )
    pvalue = wrap(
        pval_values, "pvalue",
        {
            "long_name": "Two-tailed p-value",
            "description": (
                "corrected for serial correlation using the effective sample size"
                if correct_serial else
                "assumes independent samples; no serial-correlation correction applies to these dimensions"
            ),
        },
    )
    pvalue_unadjusted = wrap(
        pval_unadjusted, "pvalue_unadjusted",
        {"long_name": "Two-tailed p-value assuming independent samples",
         "description": "usually far too small for a time correlation; compare with 'pvalue'"},
    )
    n_samples_da = wrap(n_samples, "n_samples", {"long_name": "Number of valid sample pairs"})
    n_effective_da = wrap(
        n_effective, "n_effective",
        {"long_name": "Effective number of independent samples",
         "description": "from the lag-1 autocorrelation of both series"},
    )

    return {
        "correlation": correlation,
        "pvalue": pvalue,
        "pvalue_unadjusted": pvalue_unadjusted,
        "n_samples": n_samples_da,
        "n_effective": n_effective_da,
    }


@register_tools(
    name="detect_outliers",
    description="Flag values that are unusual against a distribution, for quality control or for finding "
        "extreme events.\n"
        "Run it on ANOMALIES whenever the axis being scanned is time. A raw ocean series is dominated by "
        "its seasonal cycle, so every summer maximum sits in the tail of the distribution and the "
        "'outliers' returned are simply the warm season.\n"
        "It asks whether a value is unusual, not whether a stated threshold was crossed: for a physical "
        "criterion use compute_exceedance_fraction or detect_threshold_events.",
    args={
        "data": "xr.DataArray, required. The field to scan.",
        "method": "str, optional, default 'zscore'. 'zscore' (standard deviations from the mean; assumes "
            "an approximately normal distribution), 'iqr' (interquartile range; non-parametric), 'mad' "
            "(median absolute deviation; highly robust), 'percentile' (direct percentile cutoffs).",
        "threshold": "float, optional. Strictness, interpreted per method: standard deviations for "
            "'zscore' (default 3.0), an IQR multiplier for 'iqr' (1.5), a scaled-MAD multiplier for 'mad' "
            "(3.0), and a symmetric percentile cutoff for 'percentile' (99.0, i.e. outside the 0.5th to "
            "99.5th percentiles).",
        "dim": "str or list[str] or None, optional, default 'time'. Dimension(s) along which the reference "
            "distribution is built: 'time' compares each cell against its own history, the horizontal pair "
            "compares each cell against the rest of the map, None pools everything."
    },
    returns="xr.DataArray, a boolean mask with the input's dimensions, True where a value is flagged. "
        "Cells with no data are False rather than missing."
)
def detect_outliers(data: xr.DataArray, method: str = "zscore", threshold: float = None, dim: Union[str, list, None] = "time") -> xr.DataArray:
    method = method.lower()
    valid_methods = ['zscore', 'iqr', 'mad', 'percentile']
    if method not in valid_methods:
        raise ValueError(f"Method '{method}' not supported. Choose from {valid_methods}")

    if threshold is None:
        default_thresholds = {
            'zscore': 3.0,
            'iqr': 1.5,
            'mad': 3.0,
            'percentile': 99.0
        }
        threshold = default_thresholds[method]
    
    # The parentheses matter: without them this read as
    # (method == 'percentile' and threshold <= 0) or threshold >= 100, which
    # rejected a perfectly valid z-score threshold with the percentile message.
    if method == 'percentile':
        if not 0 < threshold < 100:
            raise ValueError(
                f"For 'percentile' method, threshold must be 0 < threshold < 100, got {threshold}"
            )
    elif threshold <= 0:
        raise ValueError(f"Threshold must be positive, got {threshold}")
    
    if dim is None:
        dims_to_use = list(data.dims)
    elif isinstance(dim, str):
        dims_to_use = [dim]
    else:
        dims_to_use = list(dim)
    
    for d in dims_to_use:
        if d not in data.dims:
            raise ValueError(f"Dimension '{d}' not found in data. Available dimensions: {list(data.dims)}")
    
    if method == "zscore":
        # Standard score method: |x - μ| > threshold * σ
        mean = data.mean(dim=dims_to_use, skipna=True)
        std = data.std(dim=dims_to_use, skipna=True)
        std_safe = std.where(std > 0, np.nan)
        z_scores = xr.ufuncs.fabs((data - mean) / std_safe)
        outliers = z_scores > threshold
        
    elif method == "iqr":
        # Interquartile range method: x < Q1 - threshold*IQR or x > Q3 + threshold*IQR
        q1 = data.quantile(0.25, dim=dims_to_use, skipna=True)
        q3 = data.quantile(0.75, dim=dims_to_use, skipna=True)
        iqr = q3 - q1
        iqr_safe = iqr.where(iqr > 0, np.nan)
        lower_bound = q1 - threshold * iqr_safe
        upper_bound = q3 + threshold * iqr_safe
        outliers = (data < lower_bound) | (data > upper_bound)
        
    elif method == "mad":
        # Median Absolute Deviation: more robust to outliers than std
        median = data.median(dim=dims_to_use, skipna=True)
        mad = xr.ufuncs.fabs(data - median).median(dim=dims_to_use, skipna=True)
        # MAD to std conversion factor for normal distribution
        normalization_factor = 1.4826
        mad_safe = mad.where(mad > 0, np.nan)
        mad_scores = xr.ufuncs.fabs((data - median) / (mad_safe * normalization_factor))
        outliers = mad_scores > threshold
        
    elif method == "percentile":
        lower_percentile = (100.0 - threshold) / 2.0
        upper_percentile = 100.0 - lower_percentile
        lower_threshold = data.quantile(lower_percentile / 100.0, dim=dims_to_use, skipna=True)
        upper_threshold = data.quantile(upper_percentile / 100.0, dim=dims_to_use, skipna=True)
        outliers = (data > upper_threshold) | (data < lower_threshold)
    
    else:
        raise ValueError(f"Unknown method: {method}. Choose from 'zscore', 'iqr', 'mad', 'percentile'.")
    
    outliers = outliers.astype(bool).fillna(False)

    outliers = carry_provenance(outliers, data)
    outliers.name = f"outlier_mask_{method}"
    outliers.attrs.update({
        'detection_method': method,
        'threshold': float(threshold),
        'detection_dimensions': str(dims_to_use),
        'description': f'Outlier mask: True indicates outlier detected by {method} method',
        'input_data_name': data.name if data.name else 'unnamed',
        'input_data_dims': str(list(data.dims))
    })
    
    return outliers



def _coord_label(values):
    """
    One coordinate value, in a form worth printing.

    A time stamp becomes 'YYYY-MM' rather than the nanosecond count that
    datetime64.item() returns, which is what makes 'the largest departure
    happened in February 1998' readable; a numeric coordinate stays a number so
    a latitude can still be computed with.
    """
    array = np.asarray(values)
    if np.issubdtype(array.dtype, np.datetime64):
        return str(pd.Timestamp(array.item()))[:7]
    if array.dtype == object:               # cftime, on a model calendar
        value = array.item()
        return value.strftime("%Y-%m") if hasattr(value, "strftime") else str(value)
    if np.issubdtype(array.dtype, np.number):
        return float(array)
    return str(array.item())


@register_tools(
    name="compute_difference_metrics",
    description="Measure how far one field is from another: bias, RMSE, the largest single departure and "
        "where it happens.\n"
        "This is the arithmetic behind 'how well does this reproduce the reference', computed in one pass "
        "over the elements the two actually share, so every number describes the same sample.\n"
        "Over the horizontal the metrics are area-weighted by default, since an unweighted RMSE counts a "
        "polar cell as heavily as a tropical one many times its size. Along time they are unweighted. The "
        "weights come from the grid the fields are on, so two fields on different grids must be put on a "
        "common one first.\n"
        "The correlation is deliberately left out: use analyze_correlation, with weighting='area' for a "
        "pattern correlation.",
    args={
        "data": "xr.DataArray, required. The field being assessed. The sign of the bias follows it: "
            "bias = mean(data - reference).",
        "reference": "xr.DataArray, required. What it is measured against. The two are intersected first, "
            "so a longer series or a wider region is trimmed to the overlap rather than broadcast.",
        "dim": "str or list[str], optional. Dimension(s) to reduce. Default None, every shared dimension, "
            "giving one number per metric. Pass the horizontal pair for a per-time-step RMSE, or 'time' "
            "for a per-cell map.",
        "weighting": "str, optional, default 'auto'. 'auto' weights by cell area when the horizontal "
            "dimensions are among those reduced and leaves everything else unweighted; 'area' insists on "
            "it and fails if the grid cannot be identified; 'none' weights every element equally."
    },
    returns="dict. Scalars come back as floats, a partial reduction as xr.DataArray on the dimensions left "
        "over.\n"
        "'bias': mean(data - reference), the systematic offset.\n"
        "'mae', 'rmse': mean absolute error and root mean square error.\n"
        "'centered_rmse': RMSE with the bias removed, i.e. the disagreement in pattern alone. "
        "rmse^2 = bias^2 + centered_rmse^2, so the two together say whether a field is merely offset or "
        "actually shaped differently.\n"
        "'max_abs_error' and 'max_abs_error_at': the largest single departure and the coordinates where it "
        "occurs; the location is None when the reduction leaves dimensions behind.\n"
        "'std_data', 'std_reference', 'std_ratio': the two spreads over the same sample and their ratio — "
        "above 1 the field varies too much, below 1 too little.\n"
        "'mean_data', 'mean_reference': the weighted means over the shared sample, so the bias can be "
        "reported with the two numbers it came from. Quote these rather than means computed separately, "
        "which cover different elements.\n"
        "'n_valid', 'weighting': the sample size and what the weights were.\n"
        "'difference': the field data - reference itself, on the shared sample and carrying the grid "
        "identity, ready to map or to pass to compute_exceedance_fraction."
)
def compute_difference_metrics(data: xr.DataArray, reference: xr.DataArray,
                               dim: Union[str, List[str], None] = None,
                               weighting: str = "auto") -> Dict[str, Any]:
    if weighting not in ("auto", "area", "none"):
        raise ValueError(f"weighting must be 'auto', 'area' or 'none', got '{weighting}'")

    # Inner join: whatever the two do not share cannot be compared, and leaving
    # it in would let one field's extra years or cells into a mean of the other.
    data_aligned, reference_aligned = xr.align(data, reference, join="inner")
    if data_aligned.size == 0:
        raise ValueError(
            "The two fields have nothing in common once aligned: no shared time steps, levels or cells. "
            "Check that they are on the same grid (regrid_horizontal) and cover overlapping periods."
        )

    if dim is None:
        reduce_dims = list(data_aligned.dims)
    elif isinstance(dim, str):
        reduce_dims = [dim]
    else:
        reduce_dims = list(dim)
    # ['lat', 'lon'] on a curvilinear model means its index dimensions.
    reduce_dims = _alias_dims(data_aligned, reduce_dims)
    missing = [d for d in reduce_dims if d not in data_aligned.dims]
    if missing:
        raise ValueError(_missing_dims_message(data_aligned, missing))

    difference = data_aligned - reference_aligned
    # Binary arithmetic drops attrs, and with them the grid identity that lets
    # the areas be found and a curvilinear model's dimensions be called lat/lon.
    # The difference is handed back with them restored, so the caller never has
    # to subtract the two fields again by hand.
    difference_field = carry_provenance(difference.rename("difference"), data_aligned)
    difference_field.attrs.setdefault("units", data.attrs.get("units", ""))
    difference_field.attrs.setdefault(
        "long_name", f"{data.name or 'data'} - {reference.name or 'reference'}"
    )
    valid = difference.notnull()
    n_valid = valid.sum(dim=reduce_dims)

    wants_area = weighting == "area" or (
        weighting == "auto" and len(_horizontal_names(data_aligned).intersection(reduce_dims)) >= 2
    )
    weights, weighting_note = None, "unweighted (every element counts equally)"
    if wants_area:
        try:
            weights = areas_like(data_aligned).where(valid)
            weighting_note = f"area: {weights.attrs.get('area_source', 'cell areas on the fields own grid')}"
        except Exception as exc:
            if weighting == "area":
                raise ValueError(
                    f"Area weighting needs to know which grid these fields are on, and that could not be "
                    f"established ({exc.__class__.__name__}: {exc}). Load them with the data tools so they "
                    f"keep their attrs, or pass weighting='none'."
                ) from exc
            weighting_note = (f"unweighted: cell areas were unavailable ({exc.__class__.__name__}), so the "
                              f"horizontal metrics count every cell equally")

    if weights is None:
        def average(field):
            return field.mean(dim=reduce_dims, skipna=True)
    else:
        total = weights.sum(dim=reduce_dims)
        def average(field):
            return (field * weights).sum(dim=reduce_dims, skipna=True) / total

    bias = average(difference)
    mean_data = average(data_aligned.where(valid))
    mean_reference = average(reference_aligned.where(valid))
    mae = average(abs(difference))
    rmse = np.sqrt(average(difference ** 2))
    centered_rmse = np.sqrt(average((difference - bias) ** 2))
    std_data = np.sqrt(average((data_aligned.where(valid) - mean_data) ** 2))
    std_reference = np.sqrt(average((reference_aligned.where(valid) - mean_reference) ** 2))
    absolute = abs(difference)
    max_abs_error = absolute.max(dim=reduce_dims, skipna=True)

    # Where the worst disagreement sits is only a single place when the whole
    # field has been reduced to one number; a per-cell RMSE map has one such
    # place per cell, which is a field of its own and not this tool's business.
    location = None
    if max_abs_error.ndim == 0 and n_valid > 0:
        flattened = absolute.stack(_point=reduce_dims)
        worst = flattened.isel(_point=int(flattened.fillna(-np.inf).argmax("_point")))
        location = {d: _coord_label(worst[d].values) for d in reduce_dims if d in worst.coords}

    def scalar_or_field(value, name, long_name):
        value = value.rename(name) if isinstance(value, xr.DataArray) else value
        if isinstance(value, xr.DataArray):
            if value.ndim == 0:
                return float(value)
            value.attrs = {"long_name": long_name, "units": data.attrs.get("units", ""),
                           "weighting": weighting_note}
            return carry_provenance(value, data_aligned)
        return value

    return {
        "bias": scalar_or_field(bias, "bias", "Mean difference (data - reference)"),
        "mae": scalar_or_field(mae, "mae", "Mean absolute error"),
        "rmse": scalar_or_field(rmse, "rmse", "Root mean square error"),
        "centered_rmse": scalar_or_field(centered_rmse, "centered_rmse", "RMSE with the bias removed"),
        "max_abs_error": scalar_or_field(max_abs_error, "max_abs_error", "Largest absolute difference"),
        "max_abs_error_at": location,
        "std_data": scalar_or_field(std_data, "std_data", "Standard deviation of the data"),
        "std_reference": scalar_or_field(std_reference, "std_reference", "Standard deviation of the reference"),
        "std_ratio": scalar_or_field(std_data / std_reference, "std_ratio", "std(data) / std(reference)"),
        "mean_data": scalar_or_field(mean_data, "mean_data", "Mean of the data over the shared sample"),
        "mean_reference": scalar_or_field(mean_reference, "mean_reference",
                                          "Mean of the reference over the shared sample"),
        "n_valid": int(n_valid) if n_valid.ndim == 0 else n_valid.rename("n_valid"),
        "weighting": weighting_note,
        "difference": difference_field,
    }


@register_tools(
    name="detect_threshold_events",
    description="Find the episodes in a series that cross a threshold and stay across it, and report each "
        "one's timing, duration and peak.\n"
        "This is how events are defined operationally rather than by eye: a warm or cold phase of a "
        "climate index, a marine heatwave, a hypoxic spell, any run that has to persist to count. The "
        "persistence requirement is the point — a single step beyond the threshold is noise, and "
        "min_duration is what separates an event from a wobble.\n"
        "Missing values break a run rather than being skipped over, so a gap cannot join two separate "
        "episodes into one.",
    args={
        "data": "xr.DataArray, required. A one-dimensional series such as an index or a regional mean. "
            "Reduce a field to a series first; a threshold applied to a whole map does not define an "
            "event.",
        "threshold": "float, required. The value to be crossed, in the series' own units, with the sign "
            "that matches `direction`.",
        "direction": "str, optional, default 'above'. 'above' keeps runs at or above the threshold, "
            "'below' keeps runs at or below it.",
        "min_duration": "int, optional, default 5. How many consecutive steps the series must stay across "
            "the threshold to count as an event. 1 returns every crossing, however brief."
    },
    returns="pandas.DataFrame, one row per event in chronological order, with columns 'start', 'end', "
        "'duration' (steps), 'peak_value' (the most extreme value reached), 'peak_time' and 'mean_value'. "
        "An empty frame means no run met the criterion. To rank by intensity, sort on the absolute peak. "
        "DataFrame.attrs records the criterion used."
)
def detect_threshold_events(data: xr.DataArray, threshold: float, direction: str = "above",
                            min_duration: int = 5) -> pd.DataFrame:
    direction = str(direction).lower()
    if direction not in ("above", "below"):
        raise ValueError(f"direction must be 'above' or 'below', got '{direction}'")
    min_duration = int(min_duration)
    if min_duration < 1:
        raise ValueError(f"min_duration must be at least 1 time step, got {min_duration}")
    if data.ndim != 1:
        raise ValueError(
            f"An event is a run along one axis, so this needs a 1-D series; this array has dimensions "
            f"{list(data.dims)}. Average the spatial dimensions away first, for instance with "
            f"compute_basic_stat(op='area_weighted_mean', dim=['lat', 'lon'])."
        )

    axis = data.dims[0]
    values = np.asarray(data.values, dtype="float64")
    labels = data[axis].values if axis in data.coords else np.arange(values.size)

    across = values >= threshold if direction == "above" else values <= threshold
    # NaN compares false either way, so a gap already ends a run; this makes that
    # explicit rather than relying on it.
    across &= ~np.isnan(values)

    events = []
    start = None
    for position in range(values.size + 1):
        inside = position < values.size and across[position]
        if inside and start is None:
            start = position
        elif not inside and start is not None:
            if position - start >= min_duration:
                segment = values[start:position]
                peak = int(np.argmax(segment)) if direction == "above" else int(np.argmin(segment))
                events.append({
                    "start": labels[start],
                    "end": labels[position - 1],
                    "duration": position - start,
                    "peak_value": float(segment[peak]),
                    "peak_time": labels[start + peak],
                    "mean_value": float(np.mean(segment)),
                })
            start = None

    table = pd.DataFrame(events, columns=["start", "end", "duration", "peak_value", "peak_time",
                                          "mean_value"])
    table.attrs = {
        "series": data.name or "unnamed",
        "units": data.attrs.get("units", ""),
        "criterion": f"{direction} {threshold} for at least {min_duration} consecutive steps",
        "n_events": len(table),
    }
    return table


@register_tools(
    name="compute_exceedance_fraction",
    description="Measure what proportion of a domain lies beyond a threshold, by area, by volume or by "
        "count.\n"
        "Writing the comparison by hand gets this wrong twice: missing cells compare as False, so land "
        "lands in the denominator, and a proportion of the ocean is a proportion of its AREA, not of its "
        "cells. This excludes missing data and weights each cell by its area or volume on the grid the "
        "field is on.\n"
        "Every comparison is applied to the deviation from `reference`, which is zero unless set. Writing "
        "d = x - reference, the four are: 'greater' d > threshold, 'less' d < threshold, 'abs_greater' "
        "|d| > threshold, 'abs_less' |d| < threshold. So 'abs_less' with reference=16 and threshold=2 is "
        "the band 14 < x < 18, while 'abs_greater' with the same numbers is everything outside it. The "
        "inequalities are STRICT, and the test actually used is returned in 'criterion'.\n"
        "It answers a different question from its neighbours: detect_outliers asks whether a value is "
        "unusual against a distribution, detect_threshold_events asks whether a threshold was held long "
        "enough to count as an episode, and this asks how much of the domain is beyond it.",
    args={
        "data": "xr.DataArray, required. The field to test — a difference, an anomaly, a concentration, an "
            "index. It must still carry the attrs its loader set if the weights are to be found.",
        "threshold": "float, required. The distance from the reference that counts, in the field's own "
            "units. It must not be negative; for a one-sided test below a value use comparison='less'.",
        "comparison": "str, optional, default 'abs_greater'. Which side counts, applied to x - reference "
            "as above. 'abs_greater' for a discrepancy, where an error of either sign matters; 'abs_less' "
            "for the share that agrees to within a tolerance, or that falls inside a band; 'greater' and "
            "'less' for a one-sided physical criterion.",
        "dim": "str or list[str], optional. Dimension(s) the proportion is taken over. Default None, every "
            "dimension, giving one number. Pass the horizontal pair for an area fraction at each time "
            "step, or 'time' for the fraction of the record at each cell.",
        "weighting": "str, optional, default 'auto'. 'auto' weights by cell area when the horizontal "
            "dimensions are reduced, by cell volume when depth is among them too, and counts elements "
            "equally otherwise — which is what a fraction of time steps wants. 'area', 'volume' and 'none' "
            "force the choice.",
        "reference": "float, optional, default 0.0. The value the comparison is taken against, in the "
            "field's own units, so the threshold measures a departure from it. Use it for a band around a "
            "physical constant; leave it at zero when the field is already a difference or an anomaly."
    },
    returns="dict:\n"
        "'fraction_percent' and 'fraction': the answer, in percent and as a proportion. A float when "
        "everything is reduced, otherwise a DataArray on the dimensions left over.\n"
        "'exceeding_amount', 'total_amount': the weighted amounts behind the fraction — square metres for "
        "area weighting, cubic metres for volume, a count of elements when unweighted. Report the amount "
        "alongside the percentage when the question is about how much ocean.\n"
        "'n_exceeding', 'n_valid': element counts, with missing data already excluded.\n"
        "'mask': the test as a field of 1 where the criterion is met, 0 where it is not and NaN where "
        "there is no data — ready to stipple over a map.\n"
        "'weighting', 'criterion': what the weights were, and the exact test with its numbers filled in."
)
def compute_exceedance_fraction(data: xr.DataArray, threshold: float,
                                comparison: str = "abs_greater",
                                dim: Union[str, List[str], None] = None,
                                weighting: str = "auto",
                                reference: float = 0.0) -> Dict[str, Any]:
    # Every test is on the deviation from the reference, so a band around a
    # physical constant needs no arithmetic on the field itself — which would
    # drop the attributes the weighting depends on.
    tests = {
        "greater": (lambda deviation: deviation > threshold, "{d} > {t}", "x > {hi}"),
        "less": (lambda deviation: deviation < threshold, "{d} < {t}", "x < {hi}"),
        "abs_greater": (lambda deviation: abs(deviation) > threshold, "|{d}| > {t}",
                        "x < {lo} or x > {hi}"),
        "abs_less": (lambda deviation: abs(deviation) < threshold, "|{d}| < {t}", "{lo} < x < {hi}"),
    }
    if comparison not in tests:
        raise ValueError(f"comparison must be one of {sorted(tests)}, got '{comparison}'")
    if threshold < 0:
        raise ValueError(
            f"threshold is a distance from the reference and must not be negative, got {threshold}. For a "
            f"one-sided test below a value use comparison='less'."
        )
    if weighting not in ("auto", "area", "volume", "none"):
        raise ValueError(f"weighting must be 'auto', 'area', 'volume' or 'none', got '{weighting}'")

    if dim is None:
        reduce_dims = list(data.dims)
    elif isinstance(dim, str):
        reduce_dims = [dim]
    else:
        reduce_dims = list(dim)
    reduce_dims = _alias_dims(data, reduce_dims)
    missing = [d for d in reduce_dims if d not in data.dims]
    if missing:
        raise ValueError(_missing_dims_message(data, missing))

    valid = data.notnull()
    # The flag is NaN where the data is, so missing cells fall out of both the
    # numerator and the denominator instead of counting as 'not exceeding'.
    deviation = data - reference if reference else data
    flag = tests[comparison][0](deviation).astype("float64").where(valid)

    horizontal_reduced = len(_horizontal_names(data).intersection(reduce_dims)) >= 2
    choice = weighting
    if weighting == "auto":
        if horizontal_reduced and "depth" in reduce_dims:
            choice = "volume"
        elif horizontal_reduced:
            choice = "area"
        else:
            choice = "none"

    if choice == "none":
        weights = xr.ones_like(data.astype("float64")).where(valid)
        note, amount_units = "unweighted (every element counts equally)", "elements"
    else:
        try:
            weights = (volumes_like(data) if choice == "volume" else areas_like(data)).where(valid)
        except Exception as exc:
            if weighting != "auto":
                raise ValueError(
                    f"{choice.capitalize()} weighting needs to know which grid this field is on, and that "
                    f"could not be established ({exc.__class__.__name__}: {exc}). Load it with the data "
                    f"tools so it keeps its attrs, or pass weighting='none'."
                ) from exc
            weights = xr.ones_like(data.astype("float64")).where(valid)
            choice = "none"
            note = (f"unweighted: cell areas were unavailable ({exc.__class__.__name__}), so every element "
                    f"counts equally")
            amount_units = "elements"
        else:
            source = weights.attrs.get("volume_source" if choice == "volume" else "area_source",
                                       f"cell {choice}s on the field's own grid")
            note = f"{choice}: {source}"
            amount_units = "m3" if choice == "volume" else "m2"

    total_amount = weights.sum(dim=reduce_dims, skipna=True)
    exceeding_amount = (weights * flag).sum(dim=reduce_dims, skipna=True)
    fraction = exceeding_amount / total_amount
    n_valid = valid.sum(dim=reduce_dims)
    n_exceeding = flag.sum(dim=reduce_dims, skipna=True)

    def out(value, name, long_name, units):
        if isinstance(value, xr.DataArray) and value.ndim == 0:
            return float(value)
        value = value.rename(name)
        value.attrs = {"long_name": long_name, "units": units, "weighting": note}
        return carry_provenance(value, data)

    units = data.attrs.get("units", "")
    algebra, in_full = tests[comparison][1], tests[comparison][2]
    criterion = algebra.format(d="x" if not reference else f"x - {reference}", t=threshold)
    criterion += ", i.e. " + in_full.format(lo=reference - threshold, hi=reference + threshold)
    if units:
        criterion += f" [{units}]"
    mask = flag.rename("exceedance_mask")
    mask.attrs = {"long_name": f"1 where {criterion}, 0 where not, NaN where there is no data",
                  "criterion": criterion}

    return {
        "fraction_percent": out(fraction * 100, "fraction_percent", f"Percentage of the domain where {criterion}", "percent"),
        "fraction": out(fraction, "fraction", f"Proportion of the domain where {criterion}", "1"),
        "exceeding_amount": out(exceeding_amount, "exceeding_amount", f"Amount where {criterion}", amount_units),
        "total_amount": out(total_amount, "total_amount", "Amount tested", amount_units),
        "n_exceeding": out(n_exceeding, "n_exceeding", "Elements beyond the threshold", "count"),
        "n_valid": out(n_valid.astype("float64"), "n_valid", "Elements with data", "count"),
        "mask": carry_provenance(mask, data),
        "weighting": note,
        "criterion": criterion,
    }


@register_tools(
    name="fit_trend_map",
    description="Fit a linear trend at every cell of a field and report where it is significant.\n"
        "This is the gridded counterpart of fit_1D_trend. The model is always the straight line, because a "
        "map is read as a rate and every cell must be estimated the same way.\n"
        "The significance accounts for SERIAL CORRELATION, by the same arithmetic fit_1D_trend uses: the "
        "effective sample size comes from the lag-1 autocorrelation of each cell's own residuals. Without "
        "it a monthly field is significant almost everywhere, which says more about the persistence of the "
        "ocean than about any trend.\n"
        "On a date axis the slope is reported PER DECADE by default, while fit_1D_trend reports per YEAR; "
        "set per_decade=False when the two are quoted side by side.\n"
        "Running it on the raw field is fine over whole years, but the seasonal cycle inflates the lag-1 "
        "autocorrelation and makes the test conservative; removing the seasonal cycle first gives a "
        "sharper test and the same trend.",
    args={
        "data": "xr.DataArray, required. A field with a time axis and any number of other dimensions. "
            "Every cell is fitted independently; cells with too little data come back missing.",
        "dim": "str, optional, default 'time'. The dimension to fit along. A numeric coordinate such as "
            "'depth' also works, and the slope is then per metre.",
        "alpha": "float, optional, default 0.05. Significance level, so the default gives a 95 percent "
            "interval and marks 'significant' where that interval excludes zero.",
        "per_decade": "bool, optional, default True. Report the slope per decade rather than per year. "
            "Only applies on a date axis.",
        "min_periods": "int, optional, default 10. Cells with fewer valid steps are left missing rather "
            "than fitted."
    },
    returns="dict of xr.DataArray on the input's dimensions with the fitted one removed, all carrying the "
        "grid identity so they can be averaged or regridded directly:\n"
        "'slope': the trend, in the units named by 'units'. This is the map to draw.\n"
        "'ci95_low', 'ci95_high', 'slope_se': the confidence interval and standard error, "
        "serial-correlation corrected.\n"
        "'significant': True where the interval excludes zero — the mask to stipple with. Draw the slope "
        "everywhere and mark this, rather than blanking the insignificant cells: a coherent pattern that "
        "misses significance is still worth seeing.\n"
        "'p_value' and 'p_value_unadjusted': corrected, and assuming independence.\n"
        "'r2': variance explained at that cell. Low r2 with a significant slope is normal for a noisy "
        "field — a real drift buried in variability.\n"
        "'n', 'n_effective', 'lag1_autocorr': the sample size, how many independent steps it amounts to, "
        "and the autocorrelation behind that reduction.\n"
        "'units' (str): what the slope is in."
)
def fit_trend_map(data: xr.DataArray, dim: str = "time", alpha: float = 0.05,
                  per_decade: bool = True, min_periods: int = 10) -> Dict[str, Any]:
    if dim not in data.dims:
        raise ValueError(f"Dimension '{dim}' not found in data. Available: {list(data.dims)}.")
    if data.sizes[dim] < 3:
        raise ValueError(
            f"A trend needs at least 3 points along '{dim}'; this field has {data.sizes[dim]}."
        )
    if data.ndim == 1:
        raise ValueError(
            "This field is a single series, not a map. Use fit_1D_trend, which also chooses between "
            "straight, robust, quadratic, exponential and sinusoidal models."
        )

    # The axis has to carry a unit, or the slope is 'per time step' — a per-month
    # rate on a monthly field, which reads twelve times too small.
    coord = data.coords.get(dim)
    if coord is not None and np.issubdtype(coord.dtype, np.datetime64):
        elapsed = (coord.values - coord.values[0]) / np.timedelta64(1, "D") / 365.2425
        x_unit, x_origin = "year", str(np.datetime_as_string(coord.values[0], unit="D"))
    elif coord is not None and np.issubdtype(coord.dtype, np.number):
        elapsed = coord.values.astype(float) - float(coord.values[0])
        x_unit, x_origin = str(coord.attrs.get("units", dim)), str(coord.values[0])
    else:
        elapsed = np.arange(data.sizes[dim], dtype=float)
        x_unit, x_origin = f"{dim} index", "0"

    scale = 10.0 if (per_decade and x_unit == "year") else 1.0
    slope_unit = f"{data.attrs.get('units', 'value')} per {'decade' if scale == 10.0 else x_unit}"

    x = xr.DataArray(elapsed, dims=[dim], coords={dim: data[dim]} if coord is not None else None)

    # Ordinary least squares written as sums, so the whole map is fitted with
    # array arithmetic instead of a Python loop over a hundred thousand cells.
    # Each cell uses its own valid steps, so a gap costs that cell alone.
    valid = data.notnull()
    n = valid.sum(dim)
    x_valid = x.where(valid)
    mean_x = x_valid.mean(dim, skipna=True)
    mean_y = data.mean(dim, skipna=True)
    dx, dy = x_valid - mean_x, data - mean_y
    spread = (dx ** 2).sum(dim, skipna=True)
    slope = (dx * dy).sum(dim, skipna=True) / spread
    intercept = mean_y - slope * mean_x

    residuals = data - (slope * x_valid + intercept)
    rss = (residuals ** 2).sum(dim, skipna=True)
    sst = (dy ** 2).sum(dim, skipna=True)
    r2 = 1.0 - rss / sst.where(sst > 0)

    # Lag-1 autocorrelation of the residuals, cell by cell. Pairs that straddle a
    # gap drop out of both sums rather than being counted as zero.
    lag1 = (residuals.shift({dim: 1}) * residuals).sum(dim, skipna=True) / rss.where(rss > 0)
    lag1 = lag1.clip(0.0, 0.99)
    n_effective = (n * (1.0 - lag1) / (1.0 + lag1)).clip(3.0, None)
    n_effective = xr.where(n_effective > n, n, n_effective)

    se_effective = np.sqrt(rss / (n_effective - 2).clip(1.0, None) / spread.where(spread > 0))
    se_nominal = np.sqrt(rss / (n - 2).clip(1) / spread.where(spread > 0))

    dof = (n_effective - 2).clip(1.0, None)
    critical = xr.apply_ufunc(stats.t.ppf, 1.0 - alpha / 2.0, dof,
                              dask="parallelized", output_dtypes=[float])

    def two_tailed(statistic, degrees):
        return 2.0 * (1.0 - stats.t.cdf(np.abs(statistic), degrees))

    p_value = xr.apply_ufunc(two_tailed, slope / se_effective, dof,
                             dask="parallelized", output_dtypes=[float])
    p_unadjusted = xr.apply_ufunc(two_tailed, slope / se_nominal, (n - 2).clip(1),
                                  dask="parallelized", output_dtypes=[float])

    ci_low = (slope - critical * se_effective) * scale
    ci_high = (slope + critical * se_effective) * scale
    slope = slope * scale

    # A cell with too few points is not a weak trend, it is no trend at all.
    enough = n >= min_periods
    outputs = {
        "slope": (slope, "Linear trend", slope_unit),
        "ci95_low": (ci_low, f"Lower {100 * (1 - alpha):.0f}% confidence limit on the trend", slope_unit),
        "ci95_high": (ci_high, f"Upper {100 * (1 - alpha):.0f}% confidence limit on the trend", slope_unit),
        "slope_se": (se_effective * scale, "Standard error of the trend", slope_unit),
        "p_value": (p_value, "Two-tailed p-value, corrected for serial correlation", "1"),
        "p_value_unadjusted": (p_unadjusted, "Two-tailed p-value assuming independent samples", "1"),
        "r2": (r2, "Variance explained by the straight line", "1"),
        "n": (n.astype("float64"), f"Valid steps along {dim}", "count"),
        "n_effective": (n_effective, "Effective number of independent steps", "count"),
        "lag1_autocorr": (lag1, "Lag-1 autocorrelation of the residuals", "1"),
    }

    result = {}
    for name, (values, long_name, units) in outputs.items():
        field = values.where(enough).rename(name)
        field.attrs = {
            "long_name": long_name,
            "units": units,
            "trend_axis": f"{x_unit} since {x_origin}",
            "serial_correlation": "effective sample size from the lag-1 autocorrelation of the residuals",
        }
        result[name] = carry_provenance(field, data)

    significant = ((result["ci95_low"] > 0) | (result["ci95_high"] < 0)).where(enough)
    significant = significant.rename("significant")
    significant.attrs = {
        "long_name": f"Trend significantly different from zero at the {100 * (1 - alpha):.0f}% level",
        "criterion": "the confidence interval excludes zero, after the serial-correlation correction",
    }
    result["significant"] = carry_provenance(significant, data)
    result["units"] = slope_unit
    return result


@register_tools(
    name="locate_extreme",
    description="Find the single largest or smallest value in a field and say where it is.\n"
        "Searching across more than one axis at a time is the point: compute_basic_stat with op='max' "
        "returns the value and loses the position, and op='argmax' returns a position along one axis and "
        "so gives one answer per column. This returns one point, with every coordinate that locates it.\n"
        "On a curvilinear grid the index pair says nothing about where on Earth the extreme sits, so the "
        "latitude and longitude of that cell are read off and reported alongside.",
    args={
        "data": "xr.DataArray, required. Any field. Missing values are ignored rather than competing for "
            "the extreme.",
        "op": "str, optional, default 'max'. 'max' for the largest value, 'min' for the smallest, "
            "'abs_max' for the one furthest from zero in either direction — the one to use on a difference "
            "or an anomaly, where a large negative departure matters as much as a positive one. 'abs_max' "
            "reports the value with its original sign.",
        "dim": "str or list[str], optional. The dimensions to search over. Default None, every dimension. "
            "A subset is only accepted when it covers every dimension the field has; for a per-column "
            "answer use compute_basic_stat with op='argmax'."
    },
    returns="dict:\n"
        "'value' (float): the extreme value, with its original sign.\n"
        "'coordinates' (dict): every coordinate that locates it — latitude and longitude in degrees, depth "
        "in metres, time as 'YYYY-MM'. This is what to quote.\n"
        "'indices' (dict): the integer position along each searched dimension, for selecting the same "
        "point out of another array.\n"
        "'ties' (int): how many cells hold that same value. More than 1 means the field is flat at its "
        "extreme, and the position reported is only the first of several.\n"
        "'n_valid' (int) and 'description' (str): how many values were searched, and the answer as a "
        "sentence."
)
def locate_extreme(data: xr.DataArray, op: str = "max",
                   dim: Union[str, List[str], None] = None) -> Dict[str, Any]:
    op = str(op).lower()
    if op not in ("max", "min", "abs_max"):
        raise ValueError(f"op must be 'max', 'min' or 'abs_max', got '{op}'")

    search_dims = list(data.dims) if dim is None else ([dim] if isinstance(dim, str) else list(dim))
    search_dims = _alias_dims(data, search_dims)
    missing = [d for d in search_dims if d not in data.dims]
    if missing:
        raise ValueError(_missing_dims_message(data, missing))
    remaining = [d for d in data.dims if d not in search_dims]
    if remaining:
        raise ValueError(
            f"Searching only {search_dims} would leave {remaining} behind, so the answer would be one "
            f"extreme per {remaining} rather than a single point. Reduce those dimensions first, or use "
            f"compute_basic_stat(op='arg{'min' if op == 'min' else 'max'}') for a position along one axis "
            f"at every column."
        )

    n_valid = int(data.notnull().sum())
    if n_valid == 0:
        raise ValueError("Every value in this field is missing, so it has no extreme.")

    ordered = data.transpose(*search_dims)
    # abs_max searches the magnitude but reports the value as it stands, so a
    # record cold anomaly comes back negative rather than as its size alone.
    searched = abs(ordered) if op == "abs_max" else ordered
    values = np.asarray(searched.values, dtype="float64")

    flat_position = int(np.nanargmin(values) if op == "min" else np.nanargmax(values))
    position = np.unravel_index(flat_position, values.shape)
    indices = {d: int(p) for d, p in zip(search_dims, position)}

    point = ordered.isel(indices)
    value = float(point)
    extreme_of_searched = float(np.asarray(searched.isel(indices)))
    ties = int((values == extreme_of_searched).sum())

    # Every coordinate that survives the selection locates the point, including
    # the 2-D lat/lon a curvilinear model carries alongside its index dimensions.
    coordinates = {
        name: _coord_label(coord.values)
        for name, coord in point.coords.items() if coord.ndim == 0
    }

    # Say it in geographic terms wherever they exist: a model's index pair locates
    # nothing, and MPI-ESM1-2-LR spells its 2-D coordinates 'latitude'/'longitude'
    # rather than 'lat'/'lon', so the names come from the array's own metadata.
    geographic = (data.attrs.get("lat_coord", "lat"), data.attrs.get("lon_coord", "lon"))
    preferred = [name for name in ("time", "depth") + geographic if name in coordinates]
    spoken = preferred or list(coordinates)
    units = data.attrs.get("units", "")
    described = ", ".join(
        f"{name} {coordinates[name]:g}" if isinstance(coordinates[name], float)
        else f"{name} {coordinates[name]}" for name in spoken
    )
    label = {"max": "max", "min": "min", "abs_max": "largest departure from zero"}[op]

    return {
        "value": value,
        "coordinates": coordinates,
        "indices": indices,
        "ties": ties,
        "n_valid": n_valid,
        "description": f"{label} = {value:.6g}{' ' + units if units else ''} at {described}",
    }
