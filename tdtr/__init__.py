"""
TDTR: Fast, modular Python package for Time-Domain Thermoreflectance modeling, signal processing, fitting, and analysis.
"""

from .config import Layer, Sample, ExperimentParams
from .model import (
    lgwt,
    tdtr_temp,
    tdtr_refl,
    ss_heating,
    calculate_sensitivities,
)
from .data import (
    TDTRData,
    read_exp_data,
    correct_time_zero,
    auto_correct_phase,
    extract_interior,
    average_datasets,
    find_acoustic_peaks,
)
from .fitting import (
    FitResult,
    UncertaintyResult,
    fit_tdtr,
    compute_uncertainties,
    parameter_correlation_scan,
)
from .plotting import (
    plot_fit_result,
    plot_sensitivities,
    plot_correlation_scan,
    plot_acoustic_echo,
)

__all__ = [
    # Config
    "Layer",
    "Sample",
    "ExperimentParams",
    # Model
    "lgwt",
    "tdtr_temp",
    "tdtr_refl",
    "ss_heating",
    "calculate_sensitivities",
    # Data
    "TDTRData",
    "read_exp_data",
    "correct_time_zero",
    "auto_correct_phase",
    "extract_interior",
    "average_datasets",
    "find_acoustic_peaks",
    # Fitting
    "FitResult",
    "UncertaintyResult",
    "fit_tdtr",
    "compute_uncertainties",
    "parameter_correlation_scan",
    # Plotting
    "plot_fit_result",
    "plot_sensitivities",
    "plot_correlation_scan",
    "plot_acoustic_echo",
]
