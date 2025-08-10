import numpy as np
from wmf_py.cu_py.metrics import (
    hypsometric_points,
    hypsometric_curve,
    time_to_outlet,
    hand,
)

# Default parameters for time_to_outlet tests
PARAMS = {
    "n_overland": 0.05,
    "h_overland_m": 0.005,
    "n_channel": 0.035,
    "h_channel_m": 0.5,
    "S_min": 1e-4,
}


def basic_flow_grid():
    flowdir = np.array(
        [
            [2, 2, 2],
            [2, 2, 2],
            [0, 2, 4],
        ],
        dtype=int,
    )
    mask = np.ones_like(flowdir, dtype=bool)
    stream = np.zeros_like(flowdir, dtype=bool)
    outlet = (2, 1)
    return flowdir, mask, stream, outlet


def test_hypsometric_points_and_curve():
    dem = np.arange(100, dtype=float).reshape(10, 10)
    mask = np.ones_like(dem, dtype=bool)

    elev, area = hypsometric_points(dem, mask)
    assert np.all(np.diff(elev) >= 0)
    assert np.all(np.diff(area) >= 0)
    assert abs(elev[0] - dem.min()) < 1e-9
    assert abs(elev[-1] - dem.max()) < 1e-9

    curve = hypsometric_curve(dem, mask, nbins=10)
    assert len(curve["elev_bins"]) == 10
    assert len(curve["area_frac"]) == 10
    assert abs(curve["elev_bins"][0] - dem.min()) < 1e-9
    assert abs(curve["elev_bins"][-1] - dem.max()) < 1e-9
    assert np.all(np.diff(curve["elev_bins"]) >= 0)
    assert np.all(np.diff(curve["area_frac"]) >= 0)


def test_time_to_outlet_slope_effect():
    flowdir, mask, stream, outlet = basic_flow_grid()
    slope1 = np.full_like(flowdir, 0.01, dtype=float)
    slope2 = np.full_like(flowdir, 0.001, dtype=float)
    t1 = time_to_outlet(flowdir, slope1, mask, stream, 30.0, 30.0, PARAMS, outlet)
    t2 = time_to_outlet(flowdir, slope2, mask, stream, 30.0, 30.0, PARAMS, outlet)
    assert t1[0, 0] < t2[0, 0]
    assert t1[outlet] == 0.0
    assert t2[outlet] == 0.0


def test_time_to_outlet_manning_effect():
    flowdir, mask, stream, outlet = basic_flow_grid()
    slope = np.full_like(flowdir, 0.01, dtype=float)
    params_fast = PARAMS.copy()
    params_slow = PARAMS.copy()
    params_fast["n_overland"] = 0.05
    params_slow["n_overland"] = 0.1
    t_fast = time_to_outlet(flowdir, slope, mask, stream, 30.0, 30.0, params_fast, outlet)
    t_slow = time_to_outlet(flowdir, slope, mask, stream, 30.0, 30.0, params_slow, outlet)
    assert t_fast[0, 0] < t_slow[0, 0]


def test_time_to_outlet_channel_vs_overland():
    flowdir = np.array(
        [
            [0, 0, 2],
            [0, 0, 2],
            [0, 0, 6],
        ],
        dtype=int,
    )
    mask = np.ones_like(flowdir, dtype=bool)
    stream = np.zeros_like(flowdir, dtype=bool)
    stream[1, :] = True
    slope = np.full_like(flowdir, 0.01, dtype=float)
    outlet = (1, 2)
    t = time_to_outlet(flowdir, slope, mask, stream, 30.0, 30.0, PARAMS, outlet)
    assert t[1, 0] < t[0, 1]


def test_time_to_outlet_smin():
    flowdir, mask, stream, outlet = basic_flow_grid()
    slope = np.zeros_like(flowdir, dtype=float)
    t = time_to_outlet(flowdir, slope, mask, stream, 30.0, 30.0, PARAMS, outlet)
    assert np.isfinite(t[0, 0]) and t[0, 0] >= 0


def test_hand_basic():
    dem = np.array(
        [
            [110.0, 105.0, 100.0],
            [110.0, 105.0, 100.0],
            [110.0, 105.0, 100.0],
        ]
    )
    flowdir = np.array(
        [
            [0, 0, 2],
            [0, 0, 2],
            [0, 0, 2],
        ],
        dtype=int,
    )
    mask = np.ones_like(flowdir, dtype=bool)
    stream = np.zeros_like(flowdir, dtype=bool)
    stream[:, 2] = True
    h = hand(dem, flowdir, stream, mask)
    assert np.allclose(h[:, 2], 0.0, atol=1e-6)
    assert h[0, 0] > h[0, 1] > h[0, 2]


def test_hand_two_streams_flow_path():
    dem = np.array(
        [
            [5.0, 6.0, 10.0],
            [4.0, 5.0, 6.0],
            [0.0, 3.0, 5.0],
        ]
    )
    flowdir = np.array(
        [
            [2, 2, 3],
            [2, 3, 2],
            [0, 0, 0],
        ],
        dtype=int,
    )
    mask = np.ones_like(flowdir, dtype=bool)
    stream = np.zeros_like(flowdir, dtype=bool)
    stream[2, 0] = True
    stream[2, 2] = True
    h = hand(dem, flowdir, stream, mask)
    # Cell (0,2) flows to stream (2,0) even though (2,2) is closer
    assert h[0, 2] == dem[0, 2] - dem[2, 0]
    assert h[1, 2] == dem[1, 2] - dem[2, 2]
