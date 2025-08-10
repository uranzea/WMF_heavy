"""Core hydrological models in Python.

This module defines light–weight dataclasses for model parameters and
state variables together with a reference implementation of a simple
hydrological balance model (:func:`shia_v1`).  The function converts
rainfall and potential evapotranspiration into changes in three
conceptual storages (capilar, gravita and aquifer) and keeps track of
surface and baseflow discharges.  It is purposely compact and fully
vectorised to serve as an easily testable stand‑in while the legacy
Fortran code is ported.
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
    k_perc_cap_to_grav: float = 0.0
    k_perc_grav_to_aqu: float = 0.0
    k_baseflow: float = 0.0
    max_infil_mm_h: float = 1.0e6
    eta_priority: str = "capilar_first"
    separate_fluxes: bool = False
    save_vfluxes: bool = False
    save_storage: bool = False
    verbose: bool = False


@dataclass
class ModelState:
    """Container for the model state.

    The three storages are expressed in millimetres and are two
    dimensional arrays with shape ``(ny, nx)``.  Optional diagnostic
    variables keep the last computed surface and baseflow for each cell.
    """

    storage_capilar: np.ndarray
    storage_gravita: np.ndarray
    storage_aquifer: np.ndarray
    q_super_last: np.ndarray | None = None
    q_base_last: np.ndarray | None = None


def shia_v1(
    forcings: Dict[str, np.ndarray],
    grid: Dict[str, np.ndarray],
    params: ModelParams,
    state: ModelState,
    nsteps: int,
) -> Tuple[ModelState, Dict[str, np.ndarray]]:
    """Minimal hydrological balance model.

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
        Diagnostics with keys depending on selected options.  At least
        ``q_outlet`` is returned.
    """

    rain = forcings.get("rain_mm_h")
    if rain is None:
        raise KeyError("forcings must contain 'rain_mm_h'")

    dt = params.dt
    ny, nx = rain.shape[1:]

    dt_h = dt / 3600.0
    dt_d = dt / 86400.0

    q_outlet = np.zeros(nsteps, dtype=float)

    q_super_t = q_base_t = q_total_t = None
    if params.separate_fluxes:
        q_super_t = np.zeros((nsteps, ny, nx), dtype=float)
        q_base_t = np.zeros((nsteps, ny, nx), dtype=float)
        q_total_t = np.zeros((nsteps, ny, nx), dtype=float)

    storage_cap_t = storage_grav_t = storage_aqu_t = None
    if params.save_storage:
        storage_cap_t = np.zeros((nsteps, ny, nx), dtype=float)
        storage_grav_t = np.zeros((nsteps, ny, nx), dtype=float)
        storage_aqu_t = np.zeros((nsteps, ny, nx), dtype=float)

    mass_err_step = np.zeros(nsteps, dtype=float)
    mass_err_cum = np.zeros(nsteps, dtype=float)

    storage_cap = state.storage_capilar.copy()
    storage_grav = state.storage_gravita.copy()
    storage_aqu = state.storage_aquifer.copy()

    q_super_last = (
        np.zeros((ny, nx), dtype=float)
        if state.q_super_last is None
        else state.q_super_last.copy()
    )
    q_base_last = (
        np.zeros((ny, nx), dtype=float)
        if state.q_base_last is None
        else state.q_base_last.copy()
    )

    prev_storage_sum = np.sum(storage_cap + storage_grav + storage_aqu)
    cum_p = 0.0
    cum_err = 0.0

    for t in range(nsteps):
        rain_step = np.maximum(rain[t], 0.0) * dt_h

        et_step = np.zeros((ny, nx), dtype=float)
        if "et_mm_h" in forcings:
            et_step = np.maximum(forcings["et_mm_h"][t], 0.0) * dt_h
        elif "et_mm_d" in forcings:
            et_step = np.maximum(forcings["et_mm_d"][t], 0.0) * dt_d

        infil = np.minimum(rain_step, params.max_infil_mm_h * dt_h)
        exceso_super = rain_step - infil
        storage_cap += infil

        # Evapotranspiration limited by storage
        if params.eta_priority == "proportional":
            water_avail = storage_cap + storage_grav + storage_aqu
            frac_cap = np.zeros_like(water_avail)
            frac_grav = np.zeros_like(water_avail)
            frac_aqu = np.zeros_like(water_avail)
            np.divide(storage_cap, water_avail, out=frac_cap, where=water_avail > 0)
            np.divide(storage_grav, water_avail, out=frac_grav, where=water_avail > 0)
            np.divide(storage_aqu, water_avail, out=frac_aqu, where=water_avail > 0)
            ET_cap = np.minimum(et_step * frac_cap, storage_cap)
            ET_grav = np.minimum(et_step * frac_grav, storage_grav)
            ET_aqu = np.minimum(et_step * frac_aqu, storage_aqu)
        else:  # capilar_first
            ET_cap = np.minimum(et_step, storage_cap)
            rem = et_step - ET_cap
            ET_grav = np.minimum(rem, storage_grav)
            rem2 = rem - ET_grav
            ET_aqu = np.minimum(rem2, storage_aqu)

        storage_cap -= ET_cap
        storage_grav -= ET_grav
        storage_aqu -= ET_aqu
        ET_total = ET_cap + ET_grav + ET_aqu

        perc_cg = np.minimum(
            params.k_perc_cap_to_grav * storage_cap, storage_cap
        )
        storage_cap -= perc_cg
        storage_grav += perc_cg

        perc_ga = np.minimum(
            params.k_perc_grav_to_aqu * storage_grav, storage_grav
        )
        retorno = params.retorno_gr * storage_grav
        storage_grav = storage_grav - perc_ga - retorno
        storage_aqu += perc_ga

        q_base = np.minimum(params.k_baseflow * storage_aqu, storage_aqu)
        drena = params.drena * storage_aqu
        storage_aqu = storage_aqu - q_base + drena

        storage_cap = np.clip(storage_cap, 0.0, params.max_capilar)
        storage_grav = np.clip(storage_grav, 0.0, params.max_gravita)
        storage_aqu = np.clip(storage_aqu, 0.0, params.max_aquifer)

        q_super = exceso_super + retorno
        q_total = q_super + q_base

        q_super_last = q_super
        q_base_last = q_base

        if params.separate_fluxes:
            assert q_super_t is not None and q_base_t is not None and q_total_t is not None
            q_super_t[t] = q_super
            q_base_t[t] = q_base
            q_total_t[t] = q_total
        if params.save_storage:
            assert (
                storage_cap_t is not None
                and storage_grav_t is not None
                and storage_aqu_t is not None
            )
            storage_cap_t[t] = storage_cap
            storage_grav_t[t] = storage_grav
            storage_aqu_t[t] = storage_aqu

        q_outlet[t] = np.sum(q_total)

        new_storage_sum = np.sum(storage_cap + storage_grav + storage_aqu)
        deltaS = new_storage_sum - prev_storage_sum
        prev_storage_sum = new_storage_sum

        P = np.sum(rain_step)
        ETsum = np.sum(ET_total)
        Qsuper = np.sum(q_super)
        Qbase = np.sum(q_base)
        err = P - (deltaS + Qsuper + Qbase + ETsum)
        cum_p += P
        cum_err += err
        rel_step = abs(err) / P if P > 0 else 0.0
        rel_cum = abs(cum_err) / cum_p if cum_p > 0 else 0.0
        mass_err_step[t] = rel_step
        mass_err_cum[t] = rel_cum
        if rel_cum > 0.01:
            raise RuntimeError("cumulative mass balance error >1%")
        if params.verbose and ((t + 1) % max(1, nsteps // 10) == 0):
            total_q = Qsuper + Qbase
            print(
                f"step {t + 1}: P={P:.3f} ET={ETsum:.3f} Q={total_q:.3f} "
                f"dS={deltaS:.3f} err={rel_step:.3%}"
            )

    new_state = ModelState(storage_cap, storage_grav, storage_aqu, q_super_last, q_base_last)

    diagnostics: Dict[str, np.ndarray] = {"q_outlet": q_outlet}
    if params.separate_fluxes:
        assert q_super_t is not None and q_base_t is not None and q_total_t is not None
        diagnostics.update(
            {
                "q_super_t": q_super_t,
                "q_base_t": q_base_t,
                "q_total_t": q_total_t,
            }
        )
    if params.save_storage:
        assert (
            storage_cap_t is not None
            and storage_grav_t is not None
            and storage_aqu_t is not None
        )
        diagnostics.update(
            {
                "storage_capilar_t": storage_cap_t,
                "storage_gravita_t": storage_grav_t,
                "storage_aquifer_t": storage_aqu_t,
            }
        )

    diagnostics["diag"] = {
        "mass_error_step": mass_err_step,
        "mass_error_cum": mass_err_cum,
    }

    return new_state, diagnostics

