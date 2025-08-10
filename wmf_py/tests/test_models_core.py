import numpy as np

from wmf_py.models_py.core import ModelParams, ModelState, shia_v1


def test_mass_balance_simple() -> None:
    ny, nx, nt = 5, 5, 3
    dt = 3600.0
    rain = np.full((nt, ny, nx), 10.0)  # mm/h
    forc = {"rain_mm_h": rain}
    grid = {"dx": 30.0, "dy": 30.0}
    params = ModelParams(
        dt=dt,
        h_coef=1,
        h_exp=1,
        v_coef=1,
        v_exp=1,
        max_capilar=100,
        max_gravita=200,
        max_aquifer=400,
        drena=0.0,
        retorno_gr=0.0,
        storage_constant=0.0,
    )
    state = ModelState(np.zeros((ny, nx)), np.zeros((ny, nx)), np.zeros((ny, nx)))
    state2, out = shia_v1(forc, grid, params, state, nsteps=nt)

    assert "q_outlet" in out

    rain_total = np.sum(rain * (dt / 3600.0))
    delta_storage = (
        np.sum(state2.storage_capilar)
        + np.sum(state2.storage_gravita)
        + np.sum(state2.storage_aquifer)
        - (
            np.sum(state.storage_capilar)
            + np.sum(state.storage_gravita)
            + np.sum(state.storage_aquifer)
        )
    )
    mb = abs(rain_total - delta_storage - np.sum(out["q_outlet"])) / rain_total
    assert mb <= 0.005
