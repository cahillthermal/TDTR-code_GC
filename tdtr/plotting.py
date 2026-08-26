"""
Plotting and visualization routines for TDTR signals, fits, sensitivities, and correlation scans.
"""

from typing import Dict, List, Optional, Union, Tuple
import numpy as np
import matplotlib.pyplot as plt

from .config import Sample


def plot_fit_result(
    tdelay_data: np.ndarray,
    ratio_data: np.ndarray,
    ratio_model: np.ndarray,
    sample: Optional[Sample] = None,
    title: str = "TDTR Model Fit Result",
    save_path: Optional[str] = None,
    show: bool = True,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plots experimental signal ratio data and fitted thermal model curve.
    """
    t_ps = tdelay_data * 1e12 if np.max(tdelay_data) < 1e-3 else tdelay_data

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogx(t_ps, ratio_data, 'ob', label="Experiment", markersize=5, fillstyle='none')
    ax.semilogx(t_ps, ratio_model, 'r-', label="Model Fit", linewidth=2)

    ax.set_xlabel("Time Delay (ps)", fontsize=14)
    ax.set_ylabel("$-V_{in} / V_{out}$", fontsize=14)
    ax.set_title(title, fontsize=14)
    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    ax.legend(fontsize=12, loc="best")

    if sample is not None:
        lambda_str = " ".join([f"{k:.3g}" for k in sample.lambda_vec])
        C_str = " ".join([f"{c/1e6:.3g}" for c in sample.C_vec])
        t_str = " ".join([f"{t*1e9:.3g}" for t in sample.t_vec])

        annotation_text = (
            f"$\\lambda = [{lambda_str}]$ W/m-K\n"
            f"$C = [{C_str}] \\times 10^6$ J/m$^3$-K\n"
            f"$t = [{t_str}]$ nm"
        )
        ax.text(
            0.05,
            0.95,
            annotation_text,
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
    if show:
        plt.show()

    return fig, ax


def plot_sensitivities(
    tdelay: np.ndarray,
    sensitivities: Dict[str, np.ndarray],
    params_to_plot: Optional[List[str]] = None,
    title: str = "TDTR Ratio Sensitivities",
    save_path: Optional[str] = None,
    show: bool = True,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plots logarithmic sensitivities S_X = d ln(Ratio) / d ln(X) vs delay time.
    """
    t_ps = tdelay * 1e12 if np.max(tdelay) < 1e-3 else tdelay

    fig, ax = plt.subplots(figsize=(7, 5))

    if params_to_plot is None:
        params_to_plot = list(sensitivities.keys())

    markers = ['o', 's', '^', 'v', '<', '>', 'd', 'x', '+']
    for i, p in enumerate(params_to_plot):
        if p in sensitivities:
            marker = markers[i % len(markers)]
            ax.semilogx(t_ps, sensitivities[p], label=p, marker=marker, markersize=4, linewidth=1.5)

    ax.set_xlabel("Delay Time (ps)", fontsize=14)
    ax.set_ylabel("Ratio Sensitivity $S_X$", fontsize=14)
    ax.set_title(title, fontsize=14)
    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    ax.legend(fontsize=11, loc="best")

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
    if show:
        plt.show()

    return fig, ax


def plot_correlation_scan(
    scan_results: Dict[str, np.ndarray],
    save_path: Optional[str] = None,
    show: bool = True,
) -> Tuple[plt.Figure, List[plt.Axes]]:
    """
    Plots 2D parameter correlation curve (refit_param vs scan_param) and fit error Z vs scan_param.
    """
    scan_p = scan_results["scan_param"]
    refit_p = scan_results["refit_param"]
    scan_vals = scan_results["scan_values"]
    refit_vals = scan_results["refit_values"]
    Z_scan = scan_results["Z_scan"]
    scan_orig = scan_results["scan_fit_orig"]
    refit_orig = scan_results["refit_fit_orig"]
    Z_orig = scan_results["Z_orig"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Plot 1: Refit param vs Scan param
    ax1.plot(scan_vals, refit_vals, "o-b", label="Refit Scan", linewidth=1.5, markersize=6)
    ax1.plot(scan_orig, refit_orig, "sr", label="Original Fit", markersize=9)
    ax1.axhline(refit_orig, color="r", linestyle="--", alpha=0.7)
    ax1.axvline(scan_orig, color="r", linestyle="--", alpha=0.7)

    ax1.set_xlabel(f"{scan_p} (Scan Parameter)", fontsize=13)
    ax1.set_ylabel(f"{refit_p} (Refitted Parameter)", fontsize=13)
    ax1.set_title(f"Parameter Correlation: {refit_p} vs {scan_p}", fontsize=13)
    ax1.grid(True, linestyle="--", alpha=0.5)
    ax1.legend(fontsize=11)

    # Plot 2: Fit Error Z vs Scan param
    ax2.plot(scan_vals, Z_scan, "o-k", label="Fit Error Z", linewidth=1.5, markersize=6)
    ax2.plot(scan_orig, Z_orig, "sr", label="Original Fit Z", markersize=9)
    ax2.axhline(Z_orig, color="r", linestyle="--", alpha=0.7)
    ax2.axvline(scan_orig, color="r", linestyle="--", alpha=0.7)

    ax2.set_xlabel(f"{scan_p} (Scan Parameter)", fontsize=13)
    ax2.set_ylabel("Fit Residual Sum $Z$", fontsize=13)
    ax2.set_title(f"Residual Error $Z$ vs {scan_p}", fontsize=13)
    ax2.grid(True, linestyle="--", alpha=0.5)
    ax2.legend(fontsize=11)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
    if show:
        plt.show()

    return fig, [ax1, ax2]


def plot_acoustic_echo(
    t_ps: np.ndarray,
    signal: np.ndarray,
    echo_results: Dict[str, Union[np.ndarray, float]],
    save_path: Optional[str] = None,
    show: bool = True,
) -> Tuple[plt.Figure, List[plt.Axes]]:
    """
    Plots acoustic echo signal, baseline fit, and extracted residual pulse.
    """
    t_fit = echo_results["t_fit"]
    baseline = echo_results["baseline"]
    residual = echo_results["residual"]
    peak_time = echo_results["peak_time_ps"]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 7), sharex=True)

    mask = (t_ps >= t_fit[0]) & (t_ps <= t_fit[-1])

    ax1.plot(t_ps[mask], signal[mask], 'k-', label="Signal", linewidth=2)
    ax1.plot(t_fit, baseline, 'r--', label="Baseline Fit", linewidth=2)
    ax1.set_ylabel("Signal Amplitude", fontsize=13)
    ax1.set_title("Acoustic Echo Baseline Fit", fontsize=13)
    ax1.grid(True, linestyle="--", alpha=0.5)
    ax1.legend(fontsize=11)

    ax2.plot(t_fit, residual, 'b-', label="Residual Echo", linewidth=2)
    ax2.axvline(peak_time, color='r', linestyle=':', label=f"Echo Peak: {peak_time:.2f} ps")
    ax2.set_xlabel("Delay Time (ps)", fontsize=13)
    ax2.set_ylabel("Residual Echo", fontsize=13)
    ax2.set_title("Extracted Acoustic Echo Residual", fontsize=13)
    ax2.grid(True, linestyle="--", alpha=0.5)
    ax2.legend(fontsize=11)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
    if show:
        plt.show()

    return fig, [ax1, ax2]
