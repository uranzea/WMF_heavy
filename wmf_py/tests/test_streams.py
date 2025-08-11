import pytest

np = pytest.importorskip("numpy")

from wmf_py.cu_py import streams


def _y_network():
    """Return a simple Y-shaped network for testing."""
    stream = np.zeros((5, 5), dtype=np.uint8)
    stream[:, 2] = 1  # main stem
    stream[2, 0] = 1
    stream[2, 1] = 1

    flowdir = np.full((5, 5), -1, dtype=int)
    for r in range(4):
        flowdir[r, 2] = 2  # south
    flowdir[4, 2] = 4  # west (outlet)
    flowdir[2, 0] = 0  # east
    flowdir[2, 1] = 0  # east

    outlet = (4, 2)
    return stream, flowdir, outlet


def test_stream_find_threshold():
    acum = np.array(
        [
            [0, 0, 0, 0, 0],
            [0, 2, 3, 4, 0],
            [0, 3, 4, 5, 0],
            [0, 4, 5, 6, 0],
            [0, 0, 0, 0, 0],
        ]
    )
    counts = []
    for thr in (2, 3, 5):
        counts.append(streams.stream_find(acum, thr).sum())
    assert counts[0] >= counts[1] >= counts[2]


def test_stream_find_to_corr_removes_island():
    stream, flowdir, outlet = _y_network()
    stream[0, 4] = 1
    stream[1, 4] = 1
    flowdir[0, 4] = 2
    flowdir[1, 4] = 2

    corrected = streams.stream_find_to_corr(stream, flowdir, outlet)
    assert corrected[0, 4] == 0 and corrected[1, 4] == 0


def test_basin_netxy_find_cut():
    stream, flowdir, outlet = _y_network()
    _, _, nodes = streams.basin_netxy_find(stream, flowdir, outlet)
    edges = streams.basin_netxy_cut(stream, nodes, flowdir)
    assert len(nodes) == 4
    assert len(edges) == 3


def test_basin_netxy_remove_tributary():
    stream, flowdir, outlet = _y_network()
    stream[2, 0] = 0
    stream[2, 1] = 0
    _, _, nodes = streams.basin_netxy_find(stream, flowdir, outlet)
    edges = streams.basin_netxy_cut(stream, nodes, flowdir)
    assert len(nodes) == 2
    assert len(edges) == 1


def test_stream_seed_from_coords_tie_break():
    acum = np.zeros((5, 5), dtype=int)
    acum[1, 2] = 10
    acum[3, 1] = 10
    xll = yll = 0.0
    dx = dy = 1.0
    # Coordinates correspond to cell (1,1)
    info = streams.stream_seed_from_coords(
        acum,
        1.5,
        1.5,
        xll,
        yll,
        dx,
        dy,
        outlet_rc=(0, 0),
        search_radius_cells=3,
    )
    assert info["seed_rc"] == (1, 2)


def test_stream_seed_from_coords_truncated_window():
    acum = np.arange(9).reshape((3, 3))
    xll = yll = 0.0
    dx = dy = 1.0
    info = streams.stream_seed_from_coords(
        acum,
        0.1,
        0.1,
        xll,
        yll,
        dx,
        dy,
        outlet_rc=(0, 0),
        search_radius_cells=4,
    )
    assert info["seed_rc"] == (2, 2)


def test_stream_seed_from_coords_mask():
    acum = np.zeros((3, 3), dtype=int)
    acum[0, 1] = 5
    acum[1, 0] = 5
    mask = np.zeros_like(acum, dtype=bool)
    mask[1, 0] = True
    xll = yll = 0.0
    dx = dy = 1.0
    info = streams.stream_seed_from_coords(
        acum,
        0.2,
        0.2,
        xll,
        yll,
        dx,
        dy,
        outlet_rc=(0, 0),
        search_radius_cells=1,
        mask_basin=mask,
    )
    assert info["seed_rc"] == (1, 0)


def test_stream_threshold_nearby_selects_highest():
    acum = np.array(
        [
            [5, 4, 3],
        ]
    )
    flowdir = np.array([[0, 0, 0]])  # east
    seed_rc = (0, 0)
    outlet_rc = (0, 2)
    thr = streams.stream_threshold_nearby(
        acum, seed_rc, flowdir, outlet_rc, [5, 4, 3]
    )
    assert thr == 3


def test_stream_threshold_nearby_min_feasible():
    acum = np.array([[5, 4, 3]])
    flowdir = np.array([[0, 0, 0]])
    seed_rc = (0, 0)
    outlet_rc = (0, 2)
    thr = streams.stream_threshold_nearby(acum, seed_rc, flowdir, outlet_rc, [6, 5])
    assert thr == 5


def test_stream_find_nearby_integration():
    acum = np.array(
        [
            [5, 4, 3],
        ]
    )
    flowdir = np.array([[0, 0, 0]])
    mask = np.ones_like(acum, dtype=bool)
    res = streams.stream_find_nearby(
        acum,
        flowdir,
        x=0.5,
        y=0.5,
        xll=0.0,
        yll=0.0,
        dx=1.0,
        dy=1.0,
        outlet_rc=(0, 2),
        mask_basin=mask,
        candidate_thresholds=[5, 4, 3],
    )
    assert res["threshold"] == 3
    assert res["stream"][0, 2]
    # ensure no islands
    corrected = streams.stream_find_to_corr(res["stream"].astype(int), flowdir, (0, 2))
    assert corrected.sum() == res["stream"].sum()

