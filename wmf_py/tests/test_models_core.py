import numpy as np

from wmf_py.models_py.core import ModelParams, ModelState, shia_v1


def test_mass_balance_24_steps() -> None:
    ny, nx, nt = 5, 5, 24
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
        k_perc_cap_to_grav=0.0,
        k_perc_grav_to_aqu=0.0,
        k_baseflow=0.0,
        max_infil_mm_h=1.0e6,
        eta_priority="capilar_first",
        separate_fluxes=False,
        save_vfluxes=False,
        save_storage=False,
        verbose=False,
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
    area = grid["dx"] * grid["dy"]
    q_mm = out["q_outlet"] * dt / (area * 1e-3)
    mb = abs(rain_total - delta_storage - np.sum(q_mm)) / rain_total
    assert mb <= 0.005
    assert np.all(out["diag"]["mass_error_cum"] <= 0.005)


def test_baseflow_activation() -> None:
    dt = 3600.0
    rain = np.array([[[10.0]], [[0.0]]])  # mm/h
    forc = {"rain_mm_h": rain}
    params = ModelParams(
        dt=dt,
        h_coef=1,
        h_exp=1,
        v_coef=1,
        v_exp=1,
        max_capilar=100,
        max_gravita=100,
        max_aquifer=500,
        drena=0.0,
        retorno_gr=0.0,
        storage_constant=0.0,
        k_perc_cap_to_grav=1.0,
        k_perc_grav_to_aqu=1.0,
        k_baseflow=0.2,
        max_infil_mm_h=100.0,
        eta_priority="capilar_first",
        separate_fluxes=True,
        save_vfluxes=False,
        save_storage=False,
        verbose=False,
    )
    state = ModelState(np.zeros((1, 1)), np.zeros((1, 1)), np.zeros((1, 1)))
    state2, out = shia_v1(forc, {}, params, state, nsteps=2)

    qb = out["q_base_t"][:, 0, 0]
    assert np.isclose(qb[0], 2.0)
    assert np.isclose(qb[1], 1.6)


def test_eta_priority_capilar_first() -> None:
    dt = 3600.0
    rain = np.zeros((1, 1, 1))
    et = np.array([[[60.0]]])  # mm/h -> 60 mm/step
    forc = {"rain_mm_h": rain, "et_mm_h": et}
    params = ModelParams(
        dt=dt,
        h_coef=1,
        h_exp=1,
        v_coef=1,
        v_exp=1,
        max_capilar=200,
        max_gravita=200,
        max_aquifer=200,
        drena=0.0,
        retorno_gr=0.0,
        storage_constant=0.0,
        k_perc_cap_to_grav=0.0,
        k_perc_grav_to_aqu=0.0,
        k_baseflow=0.0,
        max_infil_mm_h=100.0,
        eta_priority="capilar_first",
        separate_fluxes=False,
        save_vfluxes=False,
        save_storage=False,
        verbose=False,
    )
    state = ModelState(np.array([[50.0]]), np.array([[20.0]]), np.array([[10.0]]))
    state2, _ = shia_v1(forc, {}, params, state, nsteps=1)
    assert np.isclose(state2.storage_capilar[0, 0], 0.0)
    assert np.isclose(state2.storage_gravita[0, 0], 10.0)
    assert np.isclose(state2.storage_aquifer[0, 0], 10.0)


def test_fluxes_and_shapes() -> None:
    dt = 3600.0
    rain = np.array([[[50.0]]])  # mm/h
    forc = {"rain_mm_h": rain}
    params = ModelParams(
        dt=dt,
        h_coef=1,
        h_exp=1,
        v_coef=1,
        v_exp=1,
        max_capilar=100,
        max_gravita=100,
        max_aquifer=100,
        drena=0.0,
        retorno_gr=0.0,
        storage_constant=0.0,
        k_perc_cap_to_grav=0.0,
        k_perc_grav_to_aqu=0.0,
        k_baseflow=0.0,
        max_infil_mm_h=10.0,
        eta_priority="capilar_first",
        separate_fluxes=True,
        save_vfluxes=False,
        save_storage=True,
        verbose=False,
    )
    state = ModelState(np.zeros((1, 1)), np.zeros((1, 1)), np.zeros((1, 1)))
    state2, out = shia_v1(forc, {}, params, state, nsteps=1)

    assert out["q_super_t"].shape == (1, 1, 1)
    assert out["q_base_t"].shape == (1, 1, 1)
    assert out["q_total_t"].shape == (1, 1, 1)
    assert out["q_outlet"].shape == (1,)
    assert out["storage_capilar_t"].shape == (1, 1, 1)
    assert out["storage_gravita_t"].shape == (1, 1, 1)
    assert out["storage_aquifer_t"].shape == (1, 1, 1)
    assert np.isclose(out["q_super_t"][0, 0, 0], 40.0)
    assert np.isclose(out["q_outlet"][0], 40.0)

    params2 = ModelParams(
        dt=dt,
        h_coef=1,
        h_exp=1,
        v_coef=1,
        v_exp=1,
        max_capilar=100,
        max_gravita=100,
        max_aquifer=100,
        drena=0.0,
        retorno_gr=0.0,
        storage_constant=0.0,
        k_perc_cap_to_grav=0.0,
        k_perc_grav_to_aqu=0.0,
        k_baseflow=0.0,
        max_infil_mm_h=10.0,
        eta_priority="capilar_first",
        separate_fluxes=False,
        save_vfluxes=False,
        save_storage=False,
        verbose=False,
    )
    _, out2 = shia_v1(forc, {}, params2, state, nsteps=1)
    assert "q_super_t" not in out2
    assert "q_base_t" not in out2
    assert "q_total_t" not in out2
    assert "storage_capilar_t" not in out2
    assert "storage_gravita_t" not in out2
    assert "storage_aquifer_t" not in out2
