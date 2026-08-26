"""
Configuration classes for TDTR physical samples and experimental setups.
"""

from dataclasses import dataclass, field
from typing import List, Union, Optional
import numpy as np


@dataclass
class Layer:
    """
    Represents a single layer in a multilayer thermal stack.

    Attributes:
        name: Name or material of the layer (e.g. 'Al', 'SiO2', 'Substrate')
        k: Thermal conductivity along cross-plane direction (W/m-K)
        C: Volumetric heat capacity (J/m^3-K)
        t: Thickness of the layer (m). For the bottom substrate, thickness is ignored in thermal response (semi-infinite).
        eta: Anisotropy ratio k_inplane / k_crossplane (kx/ky). Default is 1.0 (isotropic).
    """
    name: str
    k: float
    C: float
    t: float
    eta: float = 1.0


@dataclass
class Sample:
    """
    Represents a multilayer stack.

    Attributes:
        layers: List of Layer objects, ordered from top surface (layer 0) to substrate (last layer).
    """
    layers: List[Layer] = field(default_factory=list)

    @classmethod
    def from_arrays(
        cls,
        lambda_array: Union[List[float], np.ndarray],
        C_array: Union[List[float], np.ndarray],
        t_array: Union[List[float], np.ndarray],
        eta_array: Optional[Union[List[float], np.ndarray]] = None,
        names: Optional[List[str]] = None,
    ) -> "Sample":
        """Helper factory method to construct a Sample from parameter arrays."""
        N = len(lambda_array)
        if eta_array is None:
            eta_array = [1.0] * N
        if names is None:
            names = [f"Layer_{i+1}" for i in range(N)]

        layers = [
            Layer(
                name=names[i],
                k=float(lambda_array[i]),
                C=float(C_array[i]),
                t=float(t_array[i]),
                eta=float(eta_array[i]),
            )
            for i in range(N)
        ]
        return cls(layers=layers)

    @property
    def lambda_vec(self) -> np.ndarray:
        return np.array([layer.k for layer in self.layers], dtype=float)

    @property
    def C_vec(self) -> np.ndarray:
        return np.array([layer.C for layer in self.layers], dtype=float)

    @property
    def t_vec(self) -> np.ndarray:
        return np.array([layer.t for layer in self.layers], dtype=float)

    @property
    def eta_vec(self) -> np.ndarray:
        return np.array([layer.eta for layer in self.layers], dtype=float)

    def num_layers(self) -> int:
        return len(self.layers)


@dataclass
class ExperimentParams:
    """
    Experimental parameters for TDTR measurements.

    Attributes:
        r_pump: Pump beam 1/e^2 radius (m). Can be a scalar or 1D array matching tdelay.
        r_probe: Probe beam 1/e^2 radius (m).
        f: Laser modulation frequency (Hz).
        tau_rep: Laser repetition period (s). E.g. 1 / 74.88e6 s.
        A_pump: Pump power (W). Default 12e-3.
        TCR: Temperature coefficient of reflectance (1/K). Default 1e-4.
        nnodes: Number of Gauss-Legendre quadrature nodes for integration. Default 35.
    """
    r_pump: Union[float, np.ndarray]
    r_probe: float
    f: float
    tau_rep: float
    A_pump: float = 12e-3
    TCR: float = 1e-4
    nnodes: int = 35
