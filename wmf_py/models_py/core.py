"""Core hydrological models in Python.

This module provides dataclass containers for model parameters and state,
and a Python implementation stub for :func:`shia_v1`.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple, Any

try:  # NumPy is optional for stubs
    import numpy as np
except ModuleNotFoundError:  # pragma: no cover - executed only when numpy missing
    np = Any  # type: ignore


@dataclass
class ModelParams:
    dt: float
    h_coef: float
    h_exp: float
    v_coef: float
    v_exp: float
    max_capilar: float
    max_gravita: float
    max_aquifer: float
    drena: float
    retorno_gr: float
    storage_constant: float
    separate_fluxes: bool = False
    save_vfluxes: bool = False
    save_storage: bool = False
    verbose: bool = False


@dataclass
class ModelState:
    storage_capilar: float
    storage_gravita: float
    storage_aquifer: float


def shia_v1(
    forcings: Dict[str, np.ndarray],
    grid: Dict[str, np.ndarray],
    params: ModelParams,
    state: ModelState,
    nsteps: int,
) -> Tuple[ModelState, Dict[str, np.ndarray]]:
    """Stub for the shia_v1 model.

    Parameters
    ----------
    forcings, grid: Dict[str, np.ndarray]
        Placeholder inputs for meteorological forcings and grid data.
    params: ModelParams
        Model parameter container.
    state: ModelState
        Initial model state.
    nsteps: int
        Number of timesteps to integrate.

    Returns
    -------
    Tuple[ModelState, Dict[str, np.ndarray]]
        The updated state and an empty dictionary with diagnostic arrays.
    """
    diagnostics: Dict[str, np.ndarray] = {}
    return state, diagnostics
