"""
Experimental data import and signal processing routines for TDTR.

Includes phase correction (phase_jump minimization), time-zero alignment,
outlier filtering, data range extraction, multi-dataset averaging, and acoustic echo analysis.
"""

from dataclasses import dataclass
from typing import Tuple, List, Union, Optional
import numpy as np
from scipy.optimize import minimize_scalar, minimize
from scipy.interpolate import interp1d


@dataclass
class TDTRData:
    """
    Data container for raw and processed TDTR experimental signals.

    Attributes:
        time_raw: Raw delay time array (s).
        time_exp: Corrected delay time array (s, shifted so t=0 is zero delay).
        Vin_raw: Raw in-phase voltage signal.
        Vout_raw: Raw out-of-phase voltage signal.
        Vin: Phase and baseline corrected in-phase voltage.
        Vout: Phase and baseline corrected out-of-phase voltage.
        ratio: Experimental signal ratio -Vin / Vout.
        phase_angle: Applied phase correction angle (radians).
        t_zero_shift: Applied time-zero shift (s).
        acoustic_peak_ps: Detected acoustic echo peak position in ps (optional).
    """
    time_raw: np.ndarray
    time_exp: np.ndarray
    Vin_raw: np.ndarray
    Vout_raw: np.ndarray
    Vin: np.ndarray
    Vout: np.ndarray
    ratio: np.ndarray
    phase_angle: float
    t_zero_shift: float = 0.0
    acoustic_peak_ps: Optional[float] = None


def phase_jump_metric(theta: float, Vin: np.ndarray, Vout: np.ndarray, time_exp: np.ndarray) -> float:
    """
    Computes the phase jump metric across t=0 to find phase rotation angle theta.

    MATLAB implementation equivalent:
        Vin_n = Vin - Vout * theta
        Vout_n = Vout + Vin * theta
    """
    Vin_n = Vin - Vout * theta
    Vout_n = Vout + Vin * theta

    zero_indices = np.where(time_exp < 0)[0]
    if len(zero_indices) == 0:
        index = 1
    else:
        index = zero_indices[-1] + 1

    first_half = Vout_n[:index]
    second_half_end = min(2 * index, len(Vout_n))
    second_half = Vout_n[index:second_half_end]

    if len(first_half) == 0 or len(second_half) == 0:
        return 0.0

    return float(np.abs(np.mean(first_half) - np.mean(second_half)))


def phase_rotation(Vin: np.ndarray, Vout: np.ndarray, theta: float) -> Tuple[np.ndarray, np.ndarray]:
    """Applies phase rotation by angle theta (radians)."""
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    Vin_rot = Vin * cos_t - Vout * sin_t
    Vout_rot = Vout * cos_t + Vin * sin_t
    return Vin_rot, Vout_rot


def auto_correct_phase(
    time_exp: np.ndarray, Vin: np.ndarray, Vout: np.ndarray, max_iter: int = 3
) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Finds phase angle theta that minimizes phase jump across t=0, with outlier filtering.
    """
    total_theta = 0.0
    curr_Vin = Vin.copy()
    curr_Vout = Vout.copy()

    res1 = minimize_scalar(lambda th: phase_jump_metric(th, curr_Vin, curr_Vout, time_exp), bounds=(-0.5, 0.5), method='bounded')
    dth = float(res1.x)
    total_theta += dth
    curr_Vin, curr_Vout = phase_rotation(curr_Vin, curr_Vout, dth)

    # Identify zero crossing index
    abs_t = np.abs(time_exp)
    ii = int(np.argmin(abs_t))

    for _ in range(max_iter - 1):
        if ii > 4 and 2 * ii <= len(curr_Vout):
            std_val = float(np.std(curr_Vout[: max(1, ii - 4)]))
            mean_val = float(np.mean(curr_Vout[: 2 * ii]))

            outliers = np.abs(curr_Vout[: 2 * ii] - mean_val) > 2.5 * std_val
            curr_Vout[: 2 * ii][outliers] = mean_val

            # Refine phase angle after outlier cleanup
            res = minimize_scalar(lambda th: phase_jump_metric(th, curr_Vin, curr_Vout, time_exp), bounds=(-0.5, 0.5), method='bounded')
            dth = float(res.x)
            total_theta += dth
            curr_Vin, curr_Vout = phase_rotation(curr_Vin, curr_Vout, dth)

    Vin_rot, Vout_rot = phase_rotation(Vin, Vout, total_theta)
    return Vin_rot, Vout_rot, total_theta


def correct_time_zero(
    time_exp: np.ndarray, Vin: np.ndarray, baseline_window_ps: Tuple[float, float] = (-20.0, -5.0)
) -> Tuple[np.ndarray, float]:
    """
    Corrects time-zero delay by finding the half-max crossing on the rising edge of Vin.

    Returns:
        time_exp_shifted: Shifted delay time array (s).
        t_half: Time shift applied (s).
    """
    time_ps = time_exp * 1e12 if np.max(np.abs(time_exp)) < 1e-3 else time_exp

    baseline_mask = (time_ps >= baseline_window_ps[0]) & (time_ps <= baseline_window_ps[1])
    if np.any(baseline_mask):
        baseline = float(np.mean(Vin[baseline_mask]))
    else:
        baseline = 0.0

    idx_peak = int(np.argmax(Vin))
    Vin_max = float(Vin[idx_peak])
    half_val = (Vin_max + baseline) / 2.0

    rising_Vin = Vin[: max(1, idx_peak + 1)]
    rising_time = time_exp[: max(1, idx_peak + 1)]

    if len(rising_Vin) > 1:
        # Interpolate rising edge to find exact t_half
        interp_func = interp1d(rising_Vin, rising_time, kind='linear', fill_value='extrapolate')
        t_half = float(interp_func(half_val))
    else:
        t_half = float(rising_time[0])

    time_exp_shifted = time_exp - t_half
    return time_exp_shifted, t_half


def read_exp_data(
    filepath: str,
    time_correction: float = 0.0,
    time_unit: str = 'ps',
    autocorrect_phase: bool = True,
    phase_shift_deg: float = 0.0,
) -> TDTRData:
    """
    Reads experimental TDTR data file and performs lock-in amplifier phase correction
    and time-zero alignment.

    Expected file formats:
        Column 0: Dummy/index or time
        Column 1: Delay time (default ps)
        Column 2: Vin signal
        Column 3: Vout signal

    Args:
        filepath: Path to ASCII text data file.
        time_correction: Additional delay shift to add (s or ps depending on unit).
        time_unit: 'ps' or 's'. Delay time input unit. Default 'ps'.
        autocorrect_phase: Whether to perform automatic phase correction via phase jump minimization.
        phase_shift_deg: Manual phase shift offset angle (degrees) to add or apply for lock-in amplifier correction.

    Returns:
        TDTRData object containing raw and processed signals.
    """
    data = np.loadtxt(filepath)

    if data.ndim == 1 or data.shape[1] < 3:
        raise ValueError(f"Data file {filepath} must contain at least 3 columns.")

    if data.shape[1] >= 4:
        time_raw = data[:, 1]
        Vin_raw = data[:, 2]
        Vout_raw = data[:, 3]
    else:
        time_raw = data[:, 0]
        Vin_raw = data[:, 1]
        Vout_raw = data[:, 2]

    time_exp = time_raw + time_correction
    if time_unit == 'ps':
        time_exp_sec = time_exp * 1e-12
        time_raw_sec = time_raw * 1e-12
    else:
        time_exp_sec = time_exp
        time_raw_sec = time_raw

    # 1. Time-zero correction
    time_exp_shifted, t_half = correct_time_zero(time_exp_sec, Vin_raw)

    # 2. Lock-in Amplifier Phase Correction
    manual_rad = np.radians(phase_shift_deg)
    if autocorrect_phase:
        Vin_corr, Vout_corr, theta_auto = auto_correct_phase(time_exp_shifted, Vin_raw, Vout_raw)
        total_theta = theta_auto + manual_rad
        if manual_rad != 0.0:
            Vin_corr, Vout_corr = phase_rotation(Vin_raw, Vout_raw, total_theta)
    else:
        total_theta = manual_rad
        Vin_corr, Vout_corr = phase_rotation(Vin_raw, Vout_raw, total_theta)

    ratio = -Vin_corr / Vout_corr

    # 3. Acoustic Echo Peak Detection (default search window 20 - 30 ps)
    acoustic_peak_ps = None
    time_ps = time_exp_shifted * 1e12
    if np.any((time_ps >= 15.0) & (time_ps <= 35.0)):
        try:
            echo_res = find_acoustic_peaks(time_ps, Vin_corr, fit_window_ps=(15.0, 35.0))
            acoustic_peak_ps = echo_res["peak_time_ps"]
        except Exception:
            acoustic_peak_ps = None

    return TDTRData(
        time_raw=time_raw_sec,
        time_exp=time_exp_shifted,
        Vin_raw=Vin_raw,
        Vout_raw=Vout_raw,
        Vin=Vin_corr,
        Vout=Vout_corr,
        ratio=ratio,
        phase_angle=total_theta,
        t_zero_shift=t_half,
        acoustic_peak_ps=acoustic_peak_ps,
    )


def extract_interior(
    t_data: np.ndarray, f_data: np.ndarray, t_min: float, t_max: float
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extracts subset of data within range t_min <= t_data <= t_max.
    """
    mask = (t_data >= t_min) & (t_data <= t_max)
    return t_data[mask], f_data[mask]


def average_datasets(
    filepaths: List[str],
    time_ref: Optional[np.ndarray] = None,
    time_unit: str = 'ps',
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Processes and averages multiple TDTR experimental data files onto a common reference time axis.

    Returns:
        time_ref: Common delay time grid (s).
        Vin_avg: Averaged Vin signal.
        Vout_avg: Averaged Vout signal.
        ratio_avg: Signal ratio -Vin_avg / Vout_avg.
    """
    datasets = [read_exp_data(fp, time_unit=time_unit) for fp in filepaths]

    if time_ref is None:
        time_ref = datasets[0].time_exp

    Vin_matrix = np.zeros((len(time_ref), len(datasets)))
    Vout_matrix = np.zeros((len(time_ref), len(datasets)))

    for i, ds in enumerate(datasets):
        interp_vin = interp1d(ds.time_exp, ds.Vin, kind='linear', fill_value='extrapolate')
        interp_vout = interp1d(ds.time_exp, ds.Vout, kind='linear', fill_value='extrapolate')
        Vin_matrix[:, i] = interp_vin(time_ref)
        Vout_matrix[:, i] = interp_vout(time_ref)

    Vin_avg = np.mean(Vin_matrix, axis=1)
    Vout_avg = np.mean(Vout_matrix, axis=1)
    ratio_avg = -Vin_avg / Vout_avg

    return time_ref, Vin_avg, Vout_avg, ratio_avg


def find_acoustic_peaks(
    t_ps: np.ndarray,
    signal: np.ndarray,
    fit_window_ps: Tuple[float, float] = (20.0, 30.0),
    degree: int = 2,
) -> Dict[str, Union[np.ndarray, float]]:
    """
    Fits polynomial baseline over acoustic echo region and extracts baseline-subtracted residual echo.

    Args:
        t_ps: Delay time array (ps).
        signal: Vin, Vout, or Ratio signal array.
        fit_window_ps: Time range (ps) for acoustic echo fitting.
        degree: Polynomial degree (1 = linear, 2 = quadratic).

    Returns:
        Dictionary containing fitted baseline, residuals, and peak delay time.
    """
    mask = (t_ps >= fit_window_ps[0]) & (t_ps <= fit_window_ps[1])
    t_fit = t_ps[mask]
    y_fit = signal[mask]

    poly_coeffs = np.polyfit(t_fit, y_fit, deg=degree)
    baseline = np.polyval(poly_coeffs, t_ps)
    residual = signal - baseline

    residual_fit = residual[mask]
    peak_idx = int(np.argmax(np.abs(residual_fit)))
    peak_time_ps = float(t_fit[peak_idx])

    return {
        "t_fit": t_fit,
        "y_fit": y_fit,
        "baseline": baseline[mask],
        "residual": residual[mask],
        "peak_time_ps": peak_time_ps,
    }
