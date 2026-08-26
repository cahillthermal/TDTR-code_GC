"""
Optimization and uncertainty estimation engine for TDTR model fitting.

Supports fitting arbitrary subsets of physical parameters (conductivities, heat capacities,
thicknesses, spot sizes), calculating experimental error bars / parameter uncertainties,
and performing 2D parameter correlation scans.
"""

from dataclasses import dataclass
from typing import List, Dict, Tuple, Union, Optional, Callable
import numpy as np
from scipy.optimize import minimize

from .config import Sample, ExperimentParams
from .model import tdtr_refl


@dataclass
class FitResult:
    """
    Container for parameter fitting results.

    Attributes:
        params: Dictionary mapping parameter names to optimal fitted values.
        Z_min: Sum of squared residuals at optimal fit.
        ratio_model: Fitted model ratio curve.
        Xsol: Array of optimal fitted values.
        param_names: List of fitted parameter names.
    """
    params: Dict[str, float]
    Z_min: float
    ratio_model: np.ndarray
    Xsol: np.ndarray
    param_names: List[str]


@dataclass
class UncertaintyResult:
    """
    Container for TDTR error bar / uncertainty breakdown.

    Attributes:
        param_names: List of fitted parameter names.
        Xsol: Optimal fitted values.
        C_err: Absolute errors from heat capacity uncertainties.
        k_err: Absolute errors from thermal conductivity uncertainties.
        t_err: Absolute errors from thickness uncertainties.
        r_pump_err: Absolute errors from pump radius uncertainty.
        r_probe_err: Absolute errors from probe radius uncertainty.
        phase_err: Absolute errors from phase correction uncertainty.
        total_percent_err: Total percentage uncertainty for each parameter.
        total_abs_err: Total absolute uncertainty for each parameter.
    """
    param_names: List[str]
    Xsol: np.ndarray
    C_err: np.ndarray
    k_err: np.ndarray
    t_err: np.ndarray
    r_pump_err: np.ndarray
    r_probe_err: np.ndarray
    phase_err: np.ndarray
    total_percent_err: np.ndarray
    total_abs_err: np.ndarray


def parse_param_string(param: str) -> Tuple[str, Optional[int]]:
    """
    Parses parameter string (e.g., 'k2', 'C1', 't3', 'r_pump') into parameter type and 0-based layer index.
    """
    if param in ('r_pump', 'r_probe'):
        return param, None

    param_type = param[0]
    if param_type not in ('k', 'C', 't'):
        raise ValueError(f"Unknown parameter type in string '{param}'. Must start with k, C, t, r_pump, or r_probe.")

    layer_1indexed = int(param[1:])
    return param_type, layer_1indexed - 1


def get_param_value(param: str, sample: Sample, exp_params: ExperimentParams) -> float:
    """Gets current value of parameter from sample or exp_params."""
    p_type, layer_idx = parse_param_string(param)
    if p_type == 'r_pump':
        r_arr = np.atleast_1d(exp_params.r_pump)
        return float(r_arr[0])
    elif p_type == 'r_probe':
        return float(exp_params.r_probe)
    elif p_type == 'k':
        return float(sample.layers[layer_idx].k)
    elif p_type == 'C':
        return float(sample.layers[layer_idx].C)
    elif p_type == 't':
        return float(sample.layers[layer_idx].t)
    else:
        raise ValueError(f"Unknown parameter {param}")


def apply_params_to_sample(
    X: np.ndarray,
    param_names: List[str],
    sample_base: Sample,
    exp_params_base: ExperimentParams,
) -> Tuple[Sample, ExperimentParams]:
    """
    Returns updated Sample and ExperimentParams copies with X applied to param_names.
    """
    lambda_vec = sample_base.lambda_vec.copy()
    C_vec = sample_base.C_vec.copy()
    t_vec = sample_base.t_vec.copy()
    eta_vec = sample_base.eta_vec.copy()
    r_pump = exp_params_base.r_pump
    r_probe = exp_params_base.r_probe

    for val, name in zip(X, param_names):
        p_type, layer_idx = parse_param_string(name)
        if p_type == 'k':
            lambda_vec[layer_idx] = val
        elif p_type == 'C':
            C_vec[layer_idx] = val
        elif p_type == 't':
            t_vec[layer_idx] = val
        elif p_type == 'r_pump':
            r_pump = val
        elif p_type == 'r_probe':
            r_probe = val

    # Maintain anisotropy ratio relative to original lambda_vec change
    eta_vec = sample_base.eta_vec * sample_base.lambda_vec / lambda_vec

    updated_sample = Sample.from_arrays(
        lambda_array=lambda_vec,
        C_array=C_vec,
        t_array=t_vec,
        eta_array=eta_vec,
        names=[layer.name for layer in sample_base.layers],
    )
    updated_exp_params = ExperimentParams(
        r_pump=r_pump,
        r_probe=r_probe,
        f=exp_params_base.f,
        tau_rep=exp_params_base.tau_rep,
        A_pump=exp_params_base.A_pump,
        TCR=exp_params_base.TCR,
        nnodes=exp_params_base.nnodes,
    )
    return updated_sample, updated_exp_params


from scipy.optimize import least_squares


def fit_tdtr(
    ratio_data: np.ndarray,
    tdelay_data: np.ndarray,
    sample: Sample,
    exp_params: ExperimentParams,
    fit_params: List[str],
    Xguess: Optional[List[float]] = None,
    tol: float = 1e-4,
    method: str = 'least_squares',
) -> FitResult:
    """
    Fits specified parameters of the TDTR model to experimental ratio data.

    Args:
        ratio_data: Experimental ratio curve (-Vin/Vout).
        tdelay_data: Delay time array (s).
        sample: Sample multilayer specification.
        exp_params: Experimental setup parameters.
        fit_params: List of parameter strings to fit (e.g. ['k2', 'k4']).
        Xguess: Initial parameter guesses. If None, uses current values in sample / exp_params.
        tol: Convergence tolerance for optimization.
        method: Optimization algorithm ('least_squares', 'Nelder-Mead', 'Powell', 'L-BFGS-B').

    Returns:
        FitResult object containing optimal parameter values and model ratio curve.
    """
    if Xguess is None:
        Xguess = [get_param_value(p, sample, exp_params) for p in fit_params]

    Xguess_arr = np.asarray(Xguess, dtype=float)

    if method == 'least_squares':
        def residuals(X: np.ndarray) -> np.ndarray:
            samp, exp_p = apply_params_to_sample(X, fit_params, sample, exp_params)
            _, ratio_model = tdtr_refl(
                tdelay_data,
                exp_p.TCR,
                exp_p.tau_rep,
                exp_p.f,
                samp.lambda_vec,
                samp.C_vec,
                samp.t_vec,
                samp.eta_vec,
                exp_p.r_pump,
                exp_p.r_probe,
                exp_p.A_pump,
                exp_p.nnodes,
            )
            return ratio_model - ratio_data

        bounds = (1e-6 * np.ones_like(Xguess_arr), np.inf * np.ones_like(Xguess_arr))
        res = least_squares(residuals, Xguess_arr, bounds=bounds, ftol=tol, xtol=tol)
        Xsol = res.x
        Z_min = float(np.sum(res.fun**2))
    else:
        def objective(X: np.ndarray) -> float:
            if np.any(X <= 0):
                return 1e10
            samp, exp_p = apply_params_to_sample(X, fit_params, sample, exp_params)
            _, ratio_model = tdtr_refl(
                tdelay_data,
                exp_p.TCR,
                exp_p.tau_rep,
                exp_p.f,
                samp.lambda_vec,
                samp.C_vec,
                samp.t_vec,
                samp.eta_vec,
                exp_p.r_pump,
                exp_p.r_probe,
                exp_p.A_pump,
                exp_p.nnodes,
            )
            return float(np.sum((ratio_model - ratio_data) ** 2))

        res = minimize(objective, Xguess_arr, method=method, options={'fatol': tol, 'xatol': tol})
        Xsol = res.x
        Z_min = float(res.fun)

    fitted_sample, fitted_exp_p = apply_params_to_sample(Xsol, fit_params, sample, exp_params)
    _, ratio_model_opt = tdtr_refl(
        tdelay_data,
        fitted_exp_p.TCR,
        fitted_exp_p.tau_rep,
        fitted_exp_p.f,
        fitted_sample.lambda_vec,
        fitted_sample.C_vec,
        fitted_sample.t_vec,
        fitted_sample.eta_vec,
        fitted_exp_p.r_pump,
        fitted_exp_p.r_probe,
        fitted_exp_p.A_pump,
        fitted_exp_p.nnodes,
    )

    params_dict = {p: float(val) for p, val in zip(fit_params, Xsol)}

    return FitResult(
        params=params_dict,
        Z_min=Z_min,
        ratio_model=ratio_model_opt,
        Xsol=Xsol,
        param_names=fit_params,
    )


def compute_uncertainties(
    ratio_data: np.ndarray,
    tdelay_data: np.ndarray,
    sample: Sample,
    exp_params: ExperimentParams,
    fit_params: List[str],
    Vin_data: Optional[np.ndarray] = None,
    Vout_data: Optional[np.ndarray] = None,
    C_err_perc: Optional[List[float]] = None,
    k_err_perc: Optional[List[float]] = None,
    t_err_perc: Optional[List[float]] = None,
    r_err_perc: float = 0.05,
    degphase_err: float = 0.34,  # radphase ~0.006 rad
) -> UncertaintyResult:
    """
    Calculates parameter uncertainties / error bars propagate from system parameter uncertainties.
    """
    Nlayers = sample.num_layers()
    Nfit = len(fit_params)

    if C_err_perc is None:
        C_err_perc = [0.02] * Nlayers
    if k_err_perc is None:
        k_err_perc = [0.05] * Nlayers
    if t_err_perc is None:
        t_err_perc = [0.05] * Nlayers

    # Baseline fit
    baseline_fit = fit_tdtr(ratio_data, tdelay_data, sample, exp_params, fit_params)
    Xsol = baseline_fit.Xsol

    CErr = np.zeros((Nlayers, Nfit))
    lambdaErr = np.zeros((Nlayers, Nfit))
    tErr = np.zeros((Nlayers, Nfit))
    r_probeErr = np.zeros((1, Nfit))
    r_pumpErr = np.zeros((1, Nfit))
    phaseErr = np.zeros((1, Nfit))

    # 1. Specific heat uncertainties
    for i in range(Nlayers):
        if C_err_perc[i] > 0:
            samp_pert = Sample.from_arrays(
                sample.lambda_vec,
                sample.C_vec * (1.0 + C_err_perc[i] if np.isscalar(C_err_perc[i]) else C_err_perc[i]),
                sample.t_vec,
                sample.eta_vec,
            )
            fit_p = fit_tdtr(ratio_data, tdelay_data, samp_pert, exp_params, fit_params, Xguess=Xsol)
            CErr[i, :] = np.abs(fit_p.Xsol - Xsol)

    # 2. Thermal conductivity uncertainties
    for i in range(Nlayers):
        if k_err_perc[i] > 0:
            lamb_p = sample.lambda_vec.copy()
            lamb_p[i] *= (1.0 + k_err_perc[i])
            samp_pert = Sample.from_arrays(lamb_p, sample.C_vec, sample.t_vec, sample.eta_vec)
            fit_p = fit_tdtr(ratio_data, tdelay_data, samp_pert, exp_params, fit_params, Xguess=Xsol)
            lambdaErr[i, :] = np.abs(fit_p.Xsol - Xsol)

    # 3. Layer thickness uncertainties
    for i in range(Nlayers):
        if t_err_perc[i] > 0:
            t_p = sample.t_vec.copy()
            t_p[i] *= (1.0 + t_err_perc[i])
            samp_pert = Sample.from_arrays(sample.lambda_vec, sample.C_vec, t_p, sample.eta_vec)
            fit_p = fit_tdtr(ratio_data, tdelay_data, samp_pert, exp_params, fit_params, Xguess=Xsol)
            tErr[i, :] = np.abs(fit_p.Xsol - Xsol)

    # 4. Spot size uncertainties
    if r_err_perc > 0:
        exp_probe_pert = ExperimentParams(
            r_pump=exp_params.r_pump,
            r_probe=exp_params.r_probe * (1.0 + r_err_perc),
            f=exp_params.f,
            tau_rep=exp_params.tau_rep,
            A_pump=exp_params.A_pump,
            TCR=exp_params.TCR,
            nnodes=exp_params.nnodes,
        )
        fit_p = fit_tdtr(ratio_data, tdelay_data, sample, exp_probe_pert, fit_params, Xguess=Xsol)
        r_probeErr[0, :] = np.abs(fit_p.Xsol - Xsol)

        exp_pump_pert = ExperimentParams(
            r_pump=exp_params.r_pump * (1.0 + r_err_perc),
            r_probe=exp_params.r_probe,
            f=exp_params.f,
            tau_rep=exp_params.tau_rep,
            A_pump=exp_params.A_pump,
            TCR=exp_params.TCR,
            nnodes=exp_params.nnodes,
        )
        fit_p = fit_tdtr(ratio_data, tdelay_data, sample, exp_pump_pert, fit_params, Xguess=Xsol)
        r_pumpErr[0, :] = np.abs(fit_p.Xsol - Xsol)

    # 5. Phase error
    if Vin_data is not None and Vout_data is not None and degphase_err > 0:
        radphase = np.radians(degphase_err)
        Vtemp = (Vin_data + 1j * Vout_data) * np.exp(1j * radphase)
        ratio_shifted = -np.real(Vtemp) / np.imag(Vtemp)
        fit_p = fit_tdtr(ratio_shifted, tdelay_data, sample, exp_params, fit_params, Xguess=Xsol)
        phaseErr[0, :] = np.abs(fit_p.Xsol - Xsol)

    err_components = np.vstack([CErr, lambdaErr, tErr, r_probeErr, r_pumpErr, phaseErr])
    perc_components = err_components / Xsol[None, :]

    total_percent_err = np.sqrt(np.sum(perc_components**2, axis=0))
    total_abs_err = total_percent_err * Xsol

    return UncertaintyResult(
        param_names=fit_params,
        Xsol=Xsol,
        C_err=CErr,
        k_err=lambdaErr,
        t_err=tErr,
        r_pump_err=r_pumpErr,
        r_probe_err=r_probeErr,
        phase_err=phaseErr,
        total_percent_err=total_percent_err,
        total_abs_err=total_abs_err,
    )


def parameter_correlation_scan(
    ratio_data: np.ndarray,
    tdelay_data: np.ndarray,
    sample: Sample,
    exp_params: ExperimentParams,
    fit_params: List[str],
    scan_param: str,
    refit_param: str,
    scan_fraction: float = 0.30,
    Nscan: int = 15,
) -> Dict[str, np.ndarray]:
    """
    Performs 2D parameter correlation scan (fixes scan_param across range and refits refit_param).
    """
    if len(fit_params) != 2:
        raise ValueError("parameter_correlation_scan requires exactly two fit parameters.")
    if scan_param not in fit_params or refit_param not in fit_params:
        raise ValueError("scan_param and refit_param must be in fit_params.")

    # 1. Baseline 2-parameter fit
    baseline_fit = fit_tdtr(ratio_data, tdelay_data, sample, exp_params, fit_params)
    scan_fit = baseline_fit.params[scan_param]
    refit_fit = baseline_fit.params[refit_param]

    scan_min = (1.0 - scan_fraction) * scan_fit
    scan_max = (1.0 + scan_fraction) * scan_fit
    scan_values = np.linspace(scan_min, scan_max, Nscan)

    refit_values = np.zeros(Nscan)
    Z_scan = np.zeros(Nscan)

    for i, scan_val in enumerate(scan_values):
        # Update sample/exp_params with fixed scan_val
        samp_fixed, exp_fixed = apply_params_to_sample([scan_val], [scan_param], sample, exp_params)
        fit_sub = fit_tdtr(
            ratio_data,
            tdelay_data,
            samp_fixed,
            exp_fixed,
            [refit_param],
            Xguess=[refit_fit],
        )
        refit_values[i] = fit_sub.params[refit_param]
        Z_scan[i] = fit_sub.Z_min

    return {
        "scan_param": scan_param,
        "refit_param": refit_param,
        "scan_values": scan_values,
        "refit_values": refit_values,
        "Z_scan": Z_scan,
        "scan_fit_orig": scan_fit,
        "refit_fit_orig": refit_fit,
        "Z_orig": baseline_fit.Z_min,
    }
