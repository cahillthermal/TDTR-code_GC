"""
Core thermal physics model for Time-Domain Thermoreflectance (TDTR).

Vectorized implementation of transfer matrix heat diffusion in multilayer anisotropic media,
Hankel transform integration, Fourier summation, reflectance signal calculation,
steady-state heating, and sensitivity analysis.
"""

from typing import Tuple, Dict, Union, Optional
import numpy as np
from scipy.special import roots_legendre


def lgwt(N: int, a: float, b: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Computes Legendre-Gauss quadrature nodes and weights on [a, b].

    Args:
        N: Number of quadrature nodes.
        a: Lower bound of integration.
        b: Upper bound of integration.

    Returns:
        nodes (x): Array of shape (N,) containing node locations in [a, b].
        weights (w): Array of shape (N,) containing quadrature weights.
    """
    nodes_std, weights_std = roots_legendre(N)
    # Linear map from [-1, 1] to [a, b]
    x = 0.5 * (a * (1.0 - nodes_std) + b * (1.0 + nodes_std))
    w = 0.5 * (b - a) * weights_std
    return x, w


def tdtr_temp(
    kvectin: np.ndarray,
    freq: np.ndarray,
    lambda_vec: np.ndarray,
    C_vec: np.ndarray,
    t_vec: np.ndarray,
    eta_vec: np.ndarray,
    r_pump: Union[float, np.ndarray],
    r_probe: float,
    A_pump: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Computes frequency domain thermal Green's function and Hankel integrand for periodic Gaussian beam.

    Args:
        kvectin: 1D array of wavevectors k (m^-1), shape (Nk,).
        freq: 1D array of excitation frequencies (Hz), shape (Nf,).
        lambda_vec: 1D array of layer cross-plane thermal conductivities (W/m-K), length Nlayers.
        C_vec: 1D array of layer volumetric heat capacities (J/m^3-K), length Nlayers.
        t_vec: 1D array of layer thicknesses (m), length Nlayers.
        eta_vec: 1D array of layer anisotropy ratios kx/ky, length Nlayers.
        r_pump: Pump 1/e^2 radius (m). Scalar or 1D array of length Nt.
        r_probe: Probe 1/e^2 radius (m). Scalar.
        A_pump: Pump power (W).

    Returns:
        Integrand: Array of shape (Nk, Nf) or (Nk, Nf, Nt) representing the Hankel integrand.
        G: Array of shape (Nk, Nf) representing the surface thermal response function.
    """
    kvectin = np.asarray(kvectin, dtype=float)
    freq = np.asarray(freq, dtype=float)
    lambda_vec = np.asarray(lambda_vec, dtype=float)
    C_vec = np.asarray(C_vec, dtype=float)
    t_vec = np.asarray(t_vec, dtype=float)
    eta_vec = np.asarray(eta_vec, dtype=float)

    Nk = len(kvectin)
    Nf = len(freq)
    Nlayers = len(lambda_vec)

    kvect = np.tile(kvectin[:, None], (1, Nf))  # Shape (Nk, Nf)
    kvect2 = kvect**2
    kterm2 = 4.0 * (np.pi**2) * kvect2

    alpha = lambda_vec / C_vec
    omega = 2.0 * np.pi * freq

    # Substrate layer (last layer)
    q2 = 1j * omega / alpha[-1]
    q2_mat = q2[None, :]  # Shape (1, Nf)

    un = np.sqrt(4.0 * (np.pi**2) * eta_vec[-1] * kvect2 + q2_mat)
    gamman = lambda_vec[-1] * un

    Bplus = np.zeros((Nk, Nf), dtype=complex)
    Bminus = np.ones((Nk, Nf), dtype=complex)

    # Transfer matrix recursion from bottom to top
    if Nlayers > 1:
        for n in range(Nlayers - 1, 0, -1):
            q2_n = (1j * omega / alpha[n - 1])[None, :]
            unminus = np.sqrt(eta_vec[n - 1] * kterm2 + q2_n)
            gammanminus = lambda_vec[n - 1] * unminus

            AA = gammanminus + gamman
            BB = gammanminus - gamman

            temp1 = AA * Bplus + BB * Bminus
            temp2 = BB * Bplus + AA * Bminus

            expterm = np.exp(unminus * t_vec[n - 1])

            Bplus = (0.5 / (gammanminus * expterm)) * temp1
            Bminus = (0.5 / gammanminus) * expterm * temp2

            # Numerical stability: set deep penetration to semi-infinite
            penetration_logic = (t_vec[n - 1] * np.abs(unminus)) > 100.0
            Bplus[penetration_logic] = 0.0
            Bminus[penetration_logic] = 1.0

            un = unminus
            gamman = gammanminus

    G = (Bplus + Bminus) / (Bminus - Bplus) / gamman

    # Beam profile spatial overlap Kernel
    r_pump_arr = np.atleast_1d(r_pump)
    Nt = len(r_pump_arr)

    arg1 = -0.5 * (np.pi**2) * (r_pump_arr**2 + r_probe**2)  # Shape (Nt,)

    if Nt > 1:
        expterm_3d = np.exp(kvect2[:, :, None] * arg1[None, None, :])  # (Nk, Nf, Nt)
        Kernal = 2.0 * np.pi * A_pump * expterm_3d * kvect[:, :, None]
        Integrand = G[:, :, None] * Kernal
    else:
        expterm_2d = np.exp(kvect2 * arg1[0])  # (Nk, Nf)
        Kernal = 2.0 * np.pi * A_pump * expterm_2d * kvect
        Integrand = G * Kernal

    return Integrand, G


def tdtr_refl(
    tdelay: np.ndarray,
    TCR: float,
    tau_rep: float,
    f: float,
    lambda_vec: np.ndarray,
    C_vec: np.ndarray,
    t_vec: np.ndarray,
    eta_vec: np.ndarray,
    r_pump: Union[float, np.ndarray],
    r_probe: float,
    A_pump: float,
    nnodes: int = 35,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculates the complex reflectivity signal deltaR and ratio -Re(deltaR)/Im(deltaR).

    Args:
        tdelay: 1D array of delay times (s), shape (Ntd,).
        TCR: Temperature coefficient of reflectance (1/K).
        tau_rep: Laser repetition period (s).
        f: Laser modulation frequency (Hz).
        lambda_vec: Array of layer thermal conductivities (W/m-K).
        C_vec: Array of volumetric heat capacities (J/m^3-K).
        t_vec: Array of layer thicknesses (m).
        eta_vec: Array of anisotropy ratios (kx/ky).
        r_pump: Pump spot radius (m), scalar or array of length Ntd.
        r_probe: Probe spot radius (m), scalar.
        A_pump: Pump power (W).
        nnodes: Number of Gauss-Legendre integration nodes (default 35).

    Returns:
        deltaR: 1D complex array of shape (Ntd,) representing complex reflectance fluctuation.
        ratio: 1D real array of shape (Ntd,) representing signal ratio -Re(deltaR) / Im(deltaR).
    """
    tdelay = np.asarray(tdelay, dtype=float)
    r_pump_arr = np.atleast_1d(r_pump)

    min_tdelay = np.min(np.abs(tdelay))
    fmax = 10.0 / min_tdelay
    M = int(10 * np.ceil(tau_rep / min_tdelay))
    mvect = np.arange(-M, M + 1, dtype=float)

    fudge1 = np.exp(-np.pi * (((mvect / tau_rep) + f) / fmax) ** 2)
    fudge2 = np.exp(-np.pi * (((mvect / tau_rep) - f) / fmax) ** 2)

    kmax = 1.5 / np.sqrt(r_pump_arr[-1] ** 2 + r_probe**2)
    kvect, weights = lgwt(nnodes, 0.0, kmax)

    freq1 = (mvect / tau_rep) + f
    freq2 = (mvect / tau_rep) - f

    I1, _ = tdtr_temp(
        kvect, freq1, lambda_vec, C_vec, t_vec, eta_vec, r_pump_arr, r_probe, A_pump
    )
    I2, _ = tdtr_temp(
        kvect, freq2, lambda_vec, C_vec, t_vec, eta_vec, r_pump_arr, r_probe, A_pump
    )

    # Perform Gauss-Legendre integration over k
    if I1.ndim == 3:
        # I1 shape: (Nk, Nf, Nt) -> dT1 shape: (Nf, Nt)
        dT1 = np.tensordot(weights, I1, axes=(0, 0))
        dT2 = np.tensordot(weights, I2, axes=(0, 0))
    else:
        # I1 shape: (Nk, Nf) -> dT1 shape: (Nf,)
        dT1 = np.dot(weights, I1)
        dT2 = np.dot(weights, I2)

    expterm = np.exp(1j * (2.0 * np.pi / tau_rep) * np.outer(tdelay, mvect))

    Ntd = len(tdelay)
    if len(r_pump_arr) > 1:
        # dT1 shape (Nf, Nt) -> dT1_t shape (Nt, Nf)
        dT1_t = dT1.T
        dT2_t = dT2.T
        Retemp = (
            dT1_t * fudge1[None, :] + dT2_t * fudge2[None, :]
        ) * expterm
        Imtemp = -1j * (dT1_t - dT2_t) * expterm
    else:
        # dT1 shape (Nf,)
        NNt = np.ones((Ntd, 1))
        Retemp = (
            NNt * (dT1 * fudge1 + dT2 * fudge2)[None, :]
        ) * expterm
        Imtemp = -1j * (NNt * (dT1 - dT2)[None, :]) * expterm

    Resum = np.sum(Retemp, axis=1)
    Imsum = np.sum(Imtemp, axis=1)

    deltaRm = TCR * (Resum + 1j * Imsum)
    deltaR = deltaRm * np.exp(1j * 2.0 * np.pi * f * tdelay)
    ratio = -np.real(deltaR) / np.imag(deltaR)

    return deltaR, ratio


def ss_heating(
    lambda_vec: np.ndarray,
    C_vec: np.ndarray,
    t_vec: np.ndarray,
    eta_vec: np.ndarray,
    r: float,
    absorbance: float,
    A_tot_powermeter: float,
    nnodes: int = 100,
) -> float:
    """
    Calculates steady-state temperature rise dT_SS (K) for a multilayer stack under DC illumination.

    Args:
        lambda_vec: Array of thermal conductivities (W/m-K).
        C_vec: Array of heat capacities (J/m^3-K).
        t_vec: Array of thicknesses (m).
        eta_vec: Array of anisotropy ratios.
        r: Beam 1/e^2 radius (m).
        absorbance: Fraction of optical power absorbed.
        A_tot_powermeter: Total incident optical power (W).
        nnodes: Integration nodes (default 100).

    Returns:
        dT_SS: Maximum steady state temperature rise at beam center (K).
    """
    f = 0.0
    A_abs = absorbance * A_tot_powermeter
    kmin = 1.0 / (10000.0 * r)
    kmax = 1.5 / np.sqrt(2.0 * (r**2))

    kvect, weights = lgwt(nnodes, kmin, kmax)
    Integrand, _ = tdtr_temp(
        kvect,
        np.array([f]),
        lambda_vec,
        C_vec,
        t_vec,
        eta_vec,
        r,
        r,
        A_abs,
    )
    dT_SS = float(np.real(np.sum(weights * Integrand.ravel())))
    return dT_SS


def calculate_sensitivities(
    tdelay: np.ndarray,
    TCR: float,
    tau_rep: float,
    f: float,
    lambda_vec: np.ndarray,
    C_vec: np.ndarray,
    t_vec: np.ndarray,
    eta_vec: np.ndarray,
    r_pump: Union[float, np.ndarray],
    r_probe: float,
    A_pump: float,
    nnodes: int = 35,
    sens_params: Optional[Dict[str, float]] = None,
    delta_fraction: float = 0.01,
) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """
    Calculates logarithmic sensitivities S_X = d ln(Ratio) / d ln(X) for specified parameters.

    Args:
        tdelay: Array of delay times (s).
        TCR, tau_rep, f, lambda_vec, C_vec, t_vec, eta_vec, r_pump, r_probe, A_pump, nnodes: Model params.
        sens_params: Optional dictionary mapping parameter names (e.g. 'k1', 'k2', 't1', 'r_pump')
                    to layer indices or custom targets. If None, computes sensitivities for k, C, t of all layers,
                    plus r_pump and r_probe.
        delta_fraction: Fractional perturbation size (default 0.01 -> 1%).

    Returns:
        ratio_baseline: Baseline ratio model curve.
        sensitivities: Dictionary mapping parameter keys to sensitivity arrays S_X(tdelay).
    """
    _, ratio_baseline = tdtr_refl(
        tdelay, TCR, tau_rep, f, lambda_vec, C_vec, t_vec, eta_vec, r_pump, r_probe, A_pump, nnodes
    )

    sensitivities = {}
    Nlayers = len(lambda_vec)

    # 1. Thermal conductivities k_i (or lambda_i)
    for i in range(Nlayers):
        lambda_temp = lambda_vec.copy()
        lambda_temp[i] *= (1.0 + delta_fraction)
        eta_temp = eta_vec * lambda_vec / lambda_temp
        _, ratio_pert = tdtr_refl(
            tdelay, TCR, tau_rep, f, lambda_temp, C_vec, t_vec, eta_temp, r_pump, r_probe, A_pump, nnodes
        )
        S = (np.log(ratio_pert) - np.log(ratio_baseline)) / np.log(1.0 + delta_fraction)
        sensitivities[f"k{i+1}"] = S

    # 2. Heat capacities C_i
    for i in range(Nlayers):
        C_temp = C_vec.copy()
        C_temp[i] *= (1.0 + delta_fraction)
        _, ratio_pert = tdtr_refl(
            tdelay, TCR, tau_rep, f, lambda_vec, C_temp, t_vec, eta_vec, r_pump, r_probe, A_pump, nnodes
        )
        S = (np.log(ratio_pert) - np.log(ratio_baseline)) / np.log(1.0 + delta_fraction)
        sensitivities[f"C{i+1}"] = S

    # 3. Thicknesses t_i
    for i in range(Nlayers):
        t_temp = t_vec.copy()
        t_temp[i] *= (1.0 + delta_fraction)
        _, ratio_pert = tdtr_refl(
            tdelay, TCR, tau_rep, f, lambda_vec, C_vec, t_temp, eta_vec, r_pump, r_probe, A_pump, nnodes
        )
        S = (np.log(ratio_pert) - np.log(ratio_baseline)) / np.log(1.0 + delta_fraction)
        sensitivities[f"t{i+1}"] = S

    # 4. Spot sizes
    r_pump_temp = r_pump * (1.0 + delta_fraction)
    _, ratio_pert = tdtr_refl(
        tdelay, TCR, tau_rep, f, lambda_vec, C_vec, t_vec, eta_vec, r_pump_temp, r_probe, A_pump, nnodes
    )
    sensitivities["r_pump"] = (np.log(ratio_pert) - np.log(ratio_baseline)) / np.log(1.0 + delta_fraction)

    r_probe_temp = r_probe * (1.0 + delta_fraction)
    _, ratio_pert = tdtr_refl(
        tdelay, TCR, tau_rep, f, lambda_vec, C_vec, t_vec, eta_vec, r_pump, r_probe_temp, A_pump, nnodes
    )
    sensitivities["r_probe"] = (np.log(ratio_pert) - np.log(ratio_baseline)) / np.log(1.0 + delta_fraction)

    return ratio_baseline, sensitivities
