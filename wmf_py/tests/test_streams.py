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

