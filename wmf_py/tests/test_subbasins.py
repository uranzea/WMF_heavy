import pytest

np = pytest.importorskip("numpy")

from wmf_py.cu_py import streams, metrics


def _y_network():
    stream = np.zeros((5, 5), dtype=np.uint8)
    stream[:, 2] = 1
    stream[2, 0] = 1
    stream[2, 1] = 1

    flowdir = np.full((5, 5), -1, dtype=int)
    for r in range(4):
        flowdir[r, 2] = 2
    flowdir[4, 2] = 4
    flowdir[2, 0] = 0
    flowdir[2, 1] = 0

    outlet = (4, 2)
    return stream, flowdir, outlet


def test_subbasin_find_and_cut():
    stream, flowdir, outlet = _y_network()
    _, _, nodes = streams.basin_netxy_find(stream, flowdir, outlet)
    edges = streams.basin_netxy_cut(stream, nodes, flowdir)
    submap = metrics.basin_subbasin_find(stream, nodes, edges, flowdir, outlet)
    ids = [e["id"] for e in edges]
    assert set(np.unique(submap)) >= {0, *ids}

    masks = [metrics.basin_subbasin_cut(submap, i) for i in ids]
    total = np.sum(masks, axis=0)
    assert np.all(total <= 1)
    assert total.sum() == (submap > 0).sum()


def test_subbasin_map_and_horton():
    stream, flowdir, outlet = _y_network()
    _, _, nodes = streams.basin_netxy_find(stream, flowdir, outlet)
    edges = streams.basin_netxy_cut(stream, nodes, flowdir)
    submap = metrics.basin_subbasin_find(stream, nodes, edges, flowdir, outlet)

    vals = np.ones_like(submap, dtype=float)
    areas = metrics.basin_subbasin_map2subbasin(vals, submap)
    assert sum(areas.values()) == (submap > 0).sum()

    nodemap = metrics.basin_subbasin_nod(submap, nodes, edges)
    assert set(nodemap.keys()) == {e["id"] for e in edges}

    horton = metrics.basin_subbasin_horton(stream, nodes, edges)
    edge_orders = horton["edges"].values()
    assert list(sorted(edge_orders)) == [1, 1, 2]

