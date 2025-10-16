# tools: compute_basic_stat, fit_1D_trend

import xarray as xr
import numpy as np
import pandas as pd
from typing import Optional, Dict, Any
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

from src.tools.registry import register_tools


@register_tools(
    name="compute_basic_stat",
    description="Perform a basic statistical operation (count, max, min, mean, median, sum) " \
    "on an xarray.DataArray along specified dimension(s).",
    args={
        "data": "xr.DataArray, required. Input data array with 4 dimensions (time, depth, lat, lon).",
        "op": "str, required. The statistical operation to perform. "
            "Supported values: 'count' (number of non-NaN elements), 'max', 'min', 'mean', 'median', 'sum'.",
        "dim": "str or list of str, optional. Dimension(s) over which to perform the operation. "
            "If None, operation is applied over all dimensions."
    },
    returns="xr.DataArray. Result of the statistical operation, with reduced dimension(s)."
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
    else:
        raise ValueError(f"Unsupported operation '{op}'. "
                         f"Choose from ['count','max','min','mean','median','sum'].")

    return result


@register_tools(
    name="fit_1D_trend",
    description="Fit a trend line to a 1D xarray.DataArray using common climate trend models.",
    args={
        "data": "xr.DataArray, required. 1D data array (e.g., time series) to fit the trend.",
        "method": "str, optional. Trend fitting method. Supported values: "
        "'linear' (y = a*x + b), 'quadratic' (y = a*x^2 + b*x + c), "
        "'exponential' (y = a*exp(b*x)), 'sinusoidal' (y = a*sin(b*x + c) + d). "
        "If None, the best method is automatically chosen based on R² score."
    },
    returns="A dictionary containing: " \
    "'fitted_data' (xr.DataArray): Fitted values aligned with input." \
    "'equation' (str): String representation of the fitted function." \
    "'metrics' (dict): Difference metrics {rmse, mae, r2}." \
    "'method' (str): The fitting method used."
)
def fit_1D_trend(data: xr.DataArray, method: str = None) -> Dict[str, Any]:
    if data.ndim != 1:
        raise ValueError("Input DataArray must be 1D.")

    x = np.arange(len(data))
    y = data.values.astype(float)

    def linear(x, a, b):
        return a * x + b
    def quadratic(x, a, b, c):
        return a * x**2 + b * x + c
    def exponential(x, a, b):
        return a * np.exp(b * x)
    def sinusoidal(x, a, b, c, d):
        return a * np.sin(b * x + c) + d
    
    models = {
        "linear": (linear, [0, np.mean(y)]),
        "quadratic": (quadratic, [0, 0, np.mean(y)]),
        "exponential": (exponential, [1, 0.01]),
        "sinusoidal": (sinusoidal, [np.std(y), 2*np.pi/12, 0, np.mean(y)]),
    }

    results = {}

    for name, (func, p0) in models.items():
        try:
            params, _ = curve_fit(func, x, y, p0=p0, maxfev=10000)
            y_fit = func(x, *params)
            metrics = {
                "rmse": float(np.sqrt(mean_squared_error(y, y_fit))),
                "mae": float(mean_absolute_error(y, y_fit)),
                "r2": float(r2_score(y, y_fit)),
            }
            results[name] = {"params": params, "fitted": y_fit, "metrics": metrics}
        except Exception:
            continue
    
    if method is None:
        method = max(results.keys(), key=lambda m: results[m]["metrics"]["r2"])
    if method not in results:
        raise ValueError(f"Method '{method}' failed. Available: {list(results.keys())}")
    chosen = results[method]

    fitted_data = xr.DataArray(
        chosen["fitted"],
        coords=data.coords,
        dims=data.dims,
        attrs={**data.attrs, "fitted_by": method}
    )

    param_str = ", ".join([f"{p:.4f}" for p in chosen["params"]])
    equation = f"{method} fit with params: [{param_str}]"

    return {
        "fitted_data": fitted_data,
        "equation": equation,
        "metrics": chosen["metrics"],
    }