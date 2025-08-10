import pytest

np = pytest.importorskip("numpy")

from wmf_py.models_py.core import ModelParams, ModelState, shia_v1


def test_shia_v1_stub_runs():
    params = ModelParams(
        dt=1.0,
        h_coef=1.0,
        h_exp=1.0,
        v_coef=1.0,
        v_exp=1.0,
        max_capilar=1.0,
        max_gravita=1.0,
        max_aquifer=1.0,
        drena=0.0,
        retorno_gr=0.0,
        storage_constant=0.0,
    )
    state = ModelState(0.0, 0.0, 0.0)
    forcings = {"rain": np.zeros(1)}
    grid = {"area": np.zeros(1)}
    new_state, diag = shia_v1(forcings, grid, params, state, 1)
    assert isinstance(new_state, ModelState)
    assert isinstance(diag, dict)
