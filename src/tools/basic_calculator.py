# tools: compute_basic_stat, temporal_aggregate, fit_1D_trend, transform_unit, 
# fit_multivariate_relation, analyze_correlation

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
from sklearn.model_selection import train_test_split

from src.tools.registry import register_tools
from src.tools.data_loader import load_areas


@register_tools(
    name="compute_basic_stat",
    description="Perform a basic statistical operation (count, max, min, mean, area_weighted_mean, median, sum) " \
    "on an xarray.DataArray along specified dimension(s).",
    args={
        "data": "xr.DataArray, required. Input data array with 4 dimensions (time, depth, lat, lon).",
        "op": "str, required. The statistical operation to perform. "
            "Supported values: 'count' (number of non-NaN elements), 'max', 'min', 'mean', 'area_weighted_mean', 'median', 'sum'.",
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

    elif op == "area_weighted_mean":
        if 'lat' not in data.dims or 'lon' not in data.dims:
            raise ValueError("Area-weighted mean requires 'lat' and 'lon' dimensions in data.")

        lat_range = [float(data.lat.min()), float(data.lat.max())]
        lon_range = [float(data.lon.min()), float(data.lon.max())]
        areas = load_areas(lat_range, lon_range)
        areas = areas.sel(lat=data.lat, lon=data.lon, method="nearest")

        weighted_sum = (data * areas).sum(dim=dim)
        total_area = areas.where(~np.isnan(data)).sum(dim=dim)
        result = weighted_sum / total_area

    else:
        raise ValueError(f"Unsupported operation '{op}'. "
                         f"Choose from ['count','max','min','mean', 'area_weighted_mean', 'median','sum'].")

    return result


@register_tools(
    name="temporal_aggregate",
    description="Aggregate ocean data temporally at different resolutions (monthly, seasonal, annual) "
                "for trend analysis and seasonal cycle studies. Supports multiple aggregation operations.",
    args={
        "data": "xr.DataArray, required. Input data array with time dimension.",
        "freq": "str, optional. Temporal frequency: 'season' (DJF/MAM/JJA/SON), 'year' (annual), 'climatology' (multi-year monthly mean). Default: 'year'.",
        "op": "str, optional. Aggregation operation: 'mean', 'median', 'min', 'max', 'std'. Default: 'mean'."
    },
    returns="xr.DataArray. Temporally aggregated data with updated time coordinate."
)
def temporal_aggregate(data: xr.DataArray, freq: str = 'year', op: str = 'mean') -> xr.DataArray:
    if 'time' not in data.dims:
        raise ValueError("Input data must have a 'time' dimension")
    
    agg_funcs = {
        'mean': lambda x: x.mean(dim='time', skipna=True),
        'median': lambda x: x.median(dim='time', skipna=True),
        'min': lambda x: x.min(dim='time', skipna=True),
        'max': lambda x: x.max(dim='time', skipna=True),
        'std': lambda x: x.std(dim='time', skipna=True)
    }
    
    if op not in agg_funcs:
        raise ValueError(f"Operation '{op}' not supported. Choose from: {list(agg_funcs.keys())}")
    
    if freq == 'year':
        result = data.resample(time='1Y').map(agg_funcs[op])
        result['time'] = result['time'].dt.year
        result = result.rename({'time': 'year'})
        
    elif freq == 'season':
        # Seasonal aggregation (DJF, MAM, JJA, SON)
        # QS-DEC starts seasons in Dec-Jan-Feb pattern
        result = data.resample(time='QS-DEC').map(agg_funcs[op])
        
        season_map = {12: 'DJF', 3: 'MAM', 6: 'JJA', 9: 'SON'}
        seasons = [season_map[m] for m in result['time'].dt.month.values]
        years = result['time'].dt.year.values
        
        # Create multi-index or combined labels
        result = result.assign_coords(
            season=('time', seasons),
            year=('time', years)
        )
        
    elif freq == 'climatology':
        # Climatological monthly mean (average across all years for each month)
        result = data.groupby('time.month').map(agg_funcs[op])
        month_names = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                      'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        result = result.assign_coords(
            month_name=('month', month_names)
        )
        
    else:
        raise ValueError(f"Frequency '{freq}' not supported. "
                        f"Choose from: 'year', 'season', 'month', 'climatology'")
    
    # Preserve attributes
    result.attrs = data.attrs.copy()
    result.attrs['temporal_aggregation'] = f'{freq}_{op}'
    result.attrs['aggregation_operation'] = op
    
    return result


@register_tools(
    name="fit_1D_trend",
    description="Fit a trend line to a 1D xarray.DataArray using common climate trend models.",
    args={
        "data": "xr.DataArray, required. 1D data array (e.g., time series) to fit the trend.",
        "method": "str, optional. Trend fitting method. Supported values:"
            "'linear' (y = a*x + b), 'quadratic' (y = a*x^2 + b*x + c), "
            "'exponential' (y = a*exp(b*x)), 'sinusoidal' (y = a*sin(b*x + c) + d). "
            "If None, the best method is automatically chosen based on R² score."
    },
    returns="A dictionary containing: " \
    "'fitted_data' (xr.DataArray): Fitted values aligned with input." \
    "'equation' (str): String representation of the fitted function." \
    "'metrics' (dict): Difference metrics {rmse, mae, r2}."
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
        "linear, a * x + b": (linear, [0, np.mean(y)]),
        "quadratic, a * x**2 + b * x + c": (quadratic, [0, 0, np.mean(y)]),
        "exponential, a * np.exp(b * x)": (exponential, [1, 0.01]),
        "sinusoidal, a * np.sin(b * x + c) + d": (sinusoidal, [np.std(y), 2*np.pi/12, 0, np.mean(y)]),
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
    equation = f"{method}, fit with params: [{param_str}]"

    return {
        "fitted_data": fitted_data,
        "equation": equation,
        "metrics": chosen["metrics"],
    }


@register_tools(
    name="transform_unit",
    description="Transform units of oceanographic variables in an xarray.DataArray. "
                "Supports common conversions for temperature, concentration, salinity, and velocity. "
                "Uses standard seawater density (ρ = 1025 kg/m³) for mass/volume conversions.",
    args={
        "data": "xr.DataArray, required. Input data array to convert.",
        "from_unit": "str, required. Source unit (case-insensitive). "
                     "Supported: 'Kelvin', 'Celsius' (temperature); "
                     "'μmol/kg', 'μmol/L', 'mmol/m3', 'mmol/kg' (concentration); "
                     "'psu', 'PSU' (salinity); 'm/s', 'mm/s', 'cm/s' (velocity).",
        "to_unit": "str, required. Target unit (case-insensitive). "
                   "Must be compatible with from_unit category."
    },
    returns="xr.DataArray. Converted data with updated 'units' attribute and added 'original_units' metadata."
)
def transform_unit(data: xr.DataArray, from_unit: str, to_unit: str) -> xr.DataArray:
    from_unit = from_unit.strip().lower()
    to_unit = to_unit.strip().lower()

    if from_unit == to_unit:
        result = data.copy()
        result.attrs['units'] = to_unit
        return result
    
    unit_aliases = {
        'k': 'kelvin', 'kelvin': 'kelvin',
        '°c': 'celsius', 'degc': 'celsius', 'celsius': 'celsius', 'c': 'celsius',
        
        'μmol/kg': 'umol/kg', 'umol/kg': 'umol/kg', 'micromol/kg': 'umol/kg', 'μmol·kg^-1': 'umol/kg', 'umol·kg^-1': 'umol/kg', 
        
        'μmol/l': 'umol/l', 'umol/l': 'umol/l', 'micromol/l': 'umol/l', 'μmol·l^-1': 'umol/l', 'umol·l^-1': 'umol/l',
        
        'mmol/m3': 'mmol/m3', 'mmol/m^3': 'mmol/m3',
        
        'mmol/kg': 'mmol/kg',
        
        'psu': 'psu',
        
        'm/s': 'm/s', 'mm/s': 'mm/s', 'cm/s': 'cm/s',
    }

    std_from = unit_aliases.get(from_unit)
    std_to = unit_aliases.get(to_unit)
    
    if std_from is None:
        raise ValueError(f"Unsupported from_unit: '{from_unit}'")
    if std_to is None:
        raise ValueError(f"Unsupported to_unit: '{to_unit}'")

    conversion_key = (std_from, std_to)
    conversions = {
        ('kelvin', 'celsius'): lambda x: x - 273.15,
        ('celsius', 'kelvin'): lambda x: x + 273.15,
        
        # seawater density ρ ≈ 1025 kg/m³ = 1.025 kg/L
        ('umol/kg', 'umol/l'): lambda x: x * 1.025,
        ('umol/l', 'umol/kg'): lambda x: x / 1.025,
        
        ('umol/kg', 'mmol/m3'): lambda x: x * 1.025,
        ('mmol/m3', 'umol/kg'): lambda x: x / 1.025,
        
        ('umol/l', 'mmol/m3'): lambda x: x,
        ('mmol/m3', 'umol/l'): lambda x: x,
        
        ('umol/kg', 'mmol/kg'): lambda x: x / 1000,
        ('mmol/kg', 'umol/kg'): lambda x: x * 1000,
        
        ('umol/l', 'mmol/kg'): lambda x: x / 1000 / 1.025,
        ('mmol/kg', 'umol/l'): lambda x: x * 1000 * 1.025,
        
        ('mmol/m3', 'mmol/kg'): lambda x: x / 1.025,
        ('mmol/kg', 'mmol/m3'): lambda x: x * 1.025,
        
        ('m/s', 'mm/s'): lambda x: x * 1000,
        ('mm/s', 'm/s'): lambda x: x / 1000,
        ('m/s', 'cm/s'): lambda x: x * 100,
        ('cm/s', 'm/s'): lambda x: x / 100,
        ('cm/s', 'mm/s'): lambda x: x * 10,
        ('mm/s', 'cm/s'): lambda x: x / 10,
        
        ('psu', 'psu'): lambda x: x,
    }

    if conversion_key not in conversions:
        raise ValueError(f"Unsupported conversion from '{from_unit}' to '{to_unit}'.")

    convert_func = conversions[conversion_key]
    result = data.copy()
    result.values = convert_func(result.values)

    result.attrs['units'] = to_unit
    if 'original_units' not in result.attrs:
        result.attrs['original_units'] = from_unit
    
    return result


@register_tools(
    name="fit_multivariate_relation",
    description="Fit multivariate relationships between ocean variables. "
                "Supports linear, polynomial, and nonlinear models. "
                "If method is None, automatically selects the best method based on validation R².",
    args={
        "target": "xr.DataArray, required. One target variable to predict.",
        "predictors": "list[xr.DataArray], required. List of predictor variables. Must have compatible coordinates with target.",
        "method": "str, optional. Regression method. Options: "
                  "'linear', 'ridge', 'lasso', 'polynomial', 'random_forest'. "
                  "If None, automatically selects best method based on test R². Default: None.",
        "alpha": "float, optional. Regularization strength for ridge/lasso (default: 1.0).",
        "poly_degree": "int, optional. Polynomial degree for 'polynomial' method (default: 2).",
        "test_size": "float, optional. Fraction of data for validation (default: 0.2).",
        "standardize": "bool, optional. Whether to standardize predictors (default: True for ridge/lasso)."
    },
    returns="dict with keys: "
            "'fitted_data' (xr.DataArray): Predicted values with same coords as target. "
            "'equation' (str): Model equation or description. "
            "'metrics' (dict): {rmse, mae, r2, adj_r2, train_r2, test_r2}. "
            "'coefficients' (dict): Model parameters (if applicable). "
            "'feature_importance' (dict): Variable importance (for tree-based models). "
            "'diagnostics' (dict): Residual statistics and VIF (multicollinearity check)."
            "'method_used' (str): The method used (useful when method=None for auto-selection). "
            "'all_methods_scores' (dict): R² scores for all methods (when auto-selecting)."
)
def fit_multivariate_relation(target: xr.DataArray, predictors: List[xr.DataArray], method: str = None,
                              alpha: float = 1.0, poly_degree: int = 2, test_size: float = 0.2, standardize: bool = None) -> Dict[str, Any]:
    if not predictors:
        raise ValueError("At least one predictor variable is required.")

    predictor_names = []
    for i, pred in enumerate(predictors):
        if pred.name:
            predictor_names.append(str(pred.name))
        elif 'long_name' in pred.attrs:
            predictor_names.append(pred.attrs['long_name'])
        else:
            predictor_names.append(f"X{i+1}")
    
    # Align data
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
    
    # Train-test split for validation
    X_train, X_test, y_train, y_test = train_test_split(X_clean, y_clean, test_size=test_size, random_state=42)
    
    # Define function to fit a single method
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
    
    # Auto-select best method if method is None
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
    
    # Build equation string
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
        
        for i, name in enumerate(predictor_names):
            feature_importance[name] = float(abs(coefs[i]))
    
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
        
        # Truncate for readability
        if len(equation_terms) > 10:
            equation = f"y = {' '.join(equation_terms[:10])}... (degree={poly_degree}, {len(poly_features)} terms)"
        else:
            equation = f"y = {' '.join(equation_terms)}"
        
        # Aggregate importance by original features
        for name in predictor_names:
            importance = sum(abs(coefs[i]) for i, fn in enumerate(poly_features) if name in fn)
            feature_importance[name] = float(importance)
    
    elif method == "random_forest":
        importances = model.feature_importances_
        for i, name in enumerate(predictor_names):
            feature_importance[name] = float(importances[i])
        
        # Sort by importance
        sorted_features = sorted(feature_importance.items(), key=lambda x: x[1], reverse=True)
        top_features = ", ".join([f"{name}({imp:.3f})" for name, imp in sorted_features[:5]])
        equation = f"Random Forest ensemble (n_estimators=100, top features: {top_features})"
    
    residuals = y_clean - y_pred
    diagnostics = {
        "residual_mean": float(np.mean(residuals)),
        "residual_std": float(np.std(residuals)),
        "residual_skewness": float(np.mean((residuals / np.std(residuals))**3) if np.std(residuals) > 0 else 0),
        "max_residual": float(np.max(np.abs(residuals)))
    }
    
    # F-statistic for overall model significance (for linear models)
    if method in ["linear", "ridge", "lasso"]:
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y_clean - np.mean(y_clean))**2)
        if ss_tot > 0:
            f_stat = ((ss_tot - ss_res) / n_features) / (ss_res / (n_samples - n_features - 1))
            p_value = 1 - f_dist.cdf(f_stat, n_features, n_samples - n_features - 1)
            diagnostics["f_statistic"] = float(f_stat)
            diagnostics["p_value"] = float(p_value)
    
    # Reshape predictions back to original shape
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


@register_tools(
    name="analyze_correlation",
    description="Calculate correlation coefficients between two ocean variables with flexible dimension handling, lag analysis, and significance testing. ",
    args={
        "var1": "xr.DataArray, required. First variable.",
        "var2": "xr.DataArray, required. Second variable (must be broadcastable with var1).",
        "dim": "str or list[str], optional. Dimension(s) along which to compute correlation. Default: 'time'.",
        "method": "str, optional. Choose from 'pearson' and 'spearman'. Default: 'pearson'.",
        "lag": "int, optional. Time lag to apply to var2 (in units of the time dimension index). "
               "Positive lag means var2 leads var1. Example: lag=3 for 3-month lead. Default: 0.",
        "min_periods": "int, optional. Minimum number of valid (non-NaN) pairs required to compute correlation. Default: 3."
    },
    returns="dict with keys: "
            "'correlation' (xr.DataArray): Correlation coefficients in range [-1, 1]. "
            "'pvalue' (xr.DataArray): Two-tailed p-values for significance testing. "
            "'n_samples' (xr.DataArray): Number of valid sample pairs used for each correlation."
)
def analyze_correlation(var1: xr.DataArray, var2: xr.DataArray, dim: Union[str, List[str]] = "time",
                        method: str = "pearson", lag: int = 0, min_periods: int = 3) -> Dict[str, xr.DataArray]:
    if method not in ["pearson", "spearman"]:
        raise ValueError(f"Method must be 'pearson' or 'spearman', got '{method}'")
    
    var1_aligned, var2_aligned = xr.align(var1, var2, join='inner')
    
    if dim is None:
        dims_to_correlate = list(var1_aligned.dims)
    elif isinstance(dim, str):
        dims_to_correlate = [dim]
    else:
        dims_to_correlate = list(dim)
    
    for d in dims_to_correlate:
        if d not in var1_aligned.dims:
            raise ValueError(f"Dimension '{d}' not found in var1. Available: {var1_aligned.dims}")
        if d not in var2_aligned.dims:
            raise ValueError(f"Dimension '{d}' not found in var2. Available: {var2_aligned.dims}")
    
    # Apply lag if specified (only works if 'time' is in correlation dimensions)
    if lag != 0:
        if 'time' not in dims_to_correlate:
            raise ValueError("Lag can only be applied when 'time' is in correlation dimensions")
        if lag > 0:
            # var2 leads var1: shift var2 backward (earlier times)
            var2_aligned = var2_aligned.isel(time=slice(None, -lag))
            var1_aligned = var1_aligned.isel(time=slice(lag, None))
        else:
            # var1 leads var2: shift var1 backward
            var1_aligned = var1_aligned.isel(time=slice(None, lag))
            var2_aligned = var2_aligned.isel(time=slice(-lag, None))
    
    def correlation_func(x, y):
        """Calculate correlation with NaN handling"""
        x_flat = x.ravel()
        y_flat = y.ravel()
        valid_mask = ~(np.isnan(x_flat) | np.isnan(y_flat))
        x_valid = x_flat[valid_mask]
        y_valid = y_flat[valid_mask]
        
        n = len(x_valid)
        
        if n < min_periods:
            return np.nan, np.nan, n
        
        if method == "pearson":
            corr, pval = stats.pearsonr(x_valid, y_valid)
        else:  # spearman
            corr, pval = stats.spearmanr(x_valid, y_valid)
        
        return corr, pval, n
    
    corr_values, pval_values, n_samples = xr.apply_ufunc(
        correlation_func,
        var1_aligned,
        var2_aligned,
        input_core_dims=[dims_to_correlate, dims_to_correlate],
        output_core_dims=[[], [], []],
        vectorize=True,
        dask='parallelized',
        output_dtypes=[float, float, int]
    )
    
    remaining_dims = [d for d in var1_aligned.dims if d not in dims_to_correlate]
    remaining_coords = {d: var1_aligned.coords[d] for d in remaining_dims if d in var1_aligned.coords}
    
    correlation = xr.DataArray(
        corr_values,
        coords=remaining_coords,
        dims=remaining_dims,
        name=f"{var1.name}_vs_{var2.name}_correlation",
        attrs={
            "long_name": f"Correlation coefficient ({method})",
            "method": method,
            "correlation_dims": dims_to_correlate,
            "lag": lag
        }
    )
    
    pvalue = xr.DataArray(
        pval_values,
        coords=remaining_coords,
        dims=remaining_dims,
        name="pvalue",
        attrs={
            "long_name": "Two-tailed p-value",
            "description": "Probability that correlation is zero (null hypothesis)"
        }
    )
    
    n_samples_da = xr.DataArray(
        n_samples,
        coords=remaining_coords,
        dims=remaining_dims,
        name="n_samples",
        attrs={"long_name": "Number of valid sample pairs"}
    )
    
    return {
        "correlation": correlation,
        "pvalue": pvalue,
        "n_samples": n_samples_da
    }


@register_tools(
    name="detect_outliers",
    description="Detect and flag outliers in ocean data for quality control or extreme event identification. "
                "Supports multiple statistical methods suitable for oceanographic variables.",
    args={
        "data": "xr.DataArray, required. Input ocean data array with dimensions (time, depth, lat, lon) or subset.",
        "method": "str, optional. Detection method: 'zscore' (default), 'iqr', 'mad', 'percentile'. "
                  "'zscore': standard score method; 'iqr': interquartile range; 'mad': median absolute deviation; "
                  "'percentile': absolute percentile thresholds.",
        "threshold": "float, optional. Threshold value. For 'zscore'/'mad': number of std/MAD (default=3.0); "
                     "For 'iqr': IQR multiplier (default=1.5); For 'percentile': percentile value (default=99.0).",
        "dim": "str or list[str], optional. Dimension(s) along which to compute statistics. "
               "Default is 'time' for temporal anomalies. Use ['lat', 'lon'] for spatial, or None for global."
    },
    returns="xr.DataArray. Boolean array with same shape as input, True indicates outliers."
)
def detect_outliers(data: xr.DataArray, method: str = "zscore", threshold: float = None, dim: Union[str, list, None] = "time") -> xr.DataArray:
    # Set default thresholds for each method
    if threshold is None:
        default_thresholds = {
            'zscore': 3.0,
            'iqr': 1.5,
            'mad': 3.0,
            'percentile': 99.0
        }
        threshold = default_thresholds.get(method, 3.0)
    
    # Apply detection method
    if method == "zscore":
        # Standard score method: |x - μ| > threshold * σ
        mean = data.mean(dim=dim, skipna=True)
        std = data.std(dim=dim, skipna=True)
        z_scores = xr.ufuncs.fabs((data - mean) / std)
        outliers = z_scores > threshold
        
    elif method == "iqr":
        # Interquartile range method: x < Q1 - threshold*IQR or x > Q3 + threshold*IQR
        q1 = data.quantile(0.25, dim=dim, skipna=True)
        q3 = data.quantile(0.75, dim=dim, skipna=True)
        iqr = q3 - q1
        lower_bound = q1 - threshold * iqr
        upper_bound = q3 + threshold * iqr
        outliers = (data < lower_bound) | (data > upper_bound)
        
    elif method == "mad":
        # Median Absolute Deviation: more robust to outliers than std
        median = data.median(dim=dim, skipna=True)
        mad = xr.ufuncs.fabs(data - median).median(dim=dim, skipna=True)
        # MAD to std conversion factor for normal distribution
        mad_scores = xr.ufuncs.fabs((data - median) / (mad * 1.4826))
        outliers = mad_scores > threshold
        
    elif method == "percentile":
        # Values above the threshold percentile are marked as outliers
        upper_threshold = data.quantile(threshold / 100.0, dim=dim, skipna=True)
        lower_threshold = data.quantile((100 - threshold) / 100.0, dim=dim, skipna=True)
        outliers = (data > upper_threshold) | (data < lower_threshold)
    
    else:
        raise ValueError(f"Unknown method: {method}. Choose from 'zscore', 'iqr', 'mad', 'percentile'.")
    
    # Ensure output is boolean type
    outliers = outliers.astype(bool)
    
    outliers.attrs['detection_method'] = method
    outliers.attrs['threshold'] = threshold
    outliers.attrs['detection_dim'] = str(dim)
    outliers.attrs['description'] = f'Outlier mask: True indicates outlier detected by {method} method'
    
    return outliers

