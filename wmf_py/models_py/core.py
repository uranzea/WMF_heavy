"""Core hydrological models in Python.

This module provides dataclass containers for model parameters and state,
and a minimal Python implementation of :func:`shia_v1` suitable for unit
testing.  The function implements a very small water–balance model where
rainfall is added to a set of conceptual stores and any overflow is routed
directly to the outlet.  The goal is to provide a lightweight reference
implementation while the full Fortran model is ported.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Tuple

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
    """Container for the model state.

    The three storages are expressed in millimetres and are two
    dimensional arrays with shape ``(ny, nx)``.
    """

    storage_capilar: np.ndarray
    storage_gravita: np.ndarray
    storage_aquifer: np.ndarray


def shia_v1(
    forcings: Dict[str, np.ndarray],
    grid: Dict[str, np.ndarray],
    params: ModelParams,
    state: ModelState,
    nsteps: int,
) -> Tuple[ModelState, Dict[str, np.ndarray]]:
    """Very small water–balance model.

    Parameters
    ----------
    forcings : Dict[str, ndarray]
        ``rain_mm_h`` – rainfall rate in millimetres per hour with shape
        ``(nsteps, ny, nx)``.  Optionally ``et_mm_h`` or ``et_mm_d`` for
        evapotranspiration.  Time steps are converted internally using
        ``dt`` from ``params``.
    grid : Dict[str, ndarray]
        Placeholder grid information; only used for area aggregation.
    params : ModelParams
        Container with model parameters.  The maximum storage capacities
        are expressed in millimetres.
    state : ModelState
        Initial state of the three storages in millimetres.
    nsteps : int
        Number of time steps to integrate.

    Returns
    -------
    ModelState
        Updated model state after ``nsteps``.
    Dict[str, ndarray]
        Diagnostics with two keys:

        ``q_cauce_t`` – array with the same shape as ``rain_mm_h``
            representing the (yet to be implemented) discharge in each
            cell.  For now it is filled with zeros.
        ``q_outlet`` – one–dimensional array of length ``nsteps`` with the
            total surface excess (in mm) leaving the domain each step.

    Notes
    -----
    * Rainfall is converted from mm/h to mm per step using
      ``rain_step = rain_mm_h * (dt / 3600)``.
    * Evapotranspiration is handled similarly, accepting both hourly and
      daily rates.
    * Storages are clipped to ``[0, max_*]`` at every step.
    """

    rain = forcings.get("rain_mm_h")
    if rain is None:
        raise KeyError("forcings must contain 'rain_mm_h'")

    dt = params.dt
    ny, nx = rain.shape[1:]

    q_outlet = np.zeros(nsteps, dtype=float)
    q_cauce_t = np.zeros_like(rain)

    storage_cap = state.storage_capilar.copy()
    storage_grav = state.storage_gravita.copy()
    storage_aqu = state.storage_aquifer.copy()

    for t in range(nsteps):
        rain_step = rain[t] * (dt / 3600.0)

        et_step = np.zeros((ny, nx), dtype=float)
        if "et_mm_h" in forcings:
            et_step = forcings["et_mm_h"][t] * (dt / 3600.0)
        elif "et_mm_d" in forcings:
            et_step = forcings["et_mm_d"][t] * (dt / 86400.0)

        storage_cap += rain_step
        storage_cap -= et_step
        storage_cap = np.clip(storage_cap, 0.0, params.max_capilar)
        cap_excess = storage_cap - params.max_capilar
        cap_excess[cap_excess < 0] = 0.0

        storage_grav += cap_excess
        storage_grav = np.clip(storage_grav, 0.0, params.max_gravita)
        grav_excess = storage_grav - params.max_gravita
        grav_excess[grav_excess < 0] = 0.0

        storage_aqu += grav_excess
        storage_aqu = np.clip(storage_aqu, 0.0, params.max_aquifer)
        aqu_excess = storage_aqu - params.max_aquifer
        aqu_excess[aqu_excess < 0] = 0.0

        q_outlet[t] = np.sum(aqu_excess)

    new_state = ModelState(storage_cap, storage_grav, storage_aqu)
    diagnostics: Dict[str, np.ndarray] = {
        "q_cauce_t": q_cauce_t,
        "q_outlet": q_outlet,
    }
    return new_state, diagnostics
