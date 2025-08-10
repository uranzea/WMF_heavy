import numpy as np
import pytest
import time

from wmf_py.cu_py.basics import (
    basin_acum,
    basin_cut,
    basin_find,
    dir_reclass_opentopo,
    dir_reclass_rwatershed,
)


def test_basin_find_bounds() -> None:
    rc = basin_find(x=60, y=60, xll=0, yll=0, dx=30, dy=30, ncols=3, nrows=3)
    assert rc == (2, 2)


def test_dir_reclass_opentopo() -> None:
    src = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 1]])
    res = dir_reclass_opentopo(src)
    assert set(np.unique(res)) <= set(range(8))
    assert np.array_equal(res, dir_reclass_opentopo(res))
    with pytest.raises(ValueError):
        dir_reclass_opentopo(np.array([[9]]))


def test_dir_reclass_rwatershed_mapping() -> None:
    src = np.arange(1, 9).reshape(2, 4)
    res = dir_reclass_rwatershed(src)
    expected = np.array([0, 7, 6, 5, 4, 3, 2, 1]).reshape(2, 4)
    assert np.array_equal(res, expected)


def test_basin_acum_d8() -> None:
    flowdir = np.zeros((7, 7), dtype=int)
    flowdir[:, -1] = 2  # right column flows south
    flowdir[-1, -1] = -1  # outlet has no direction
    mask = np.ones_like(flowdir, dtype=bool)
    acc = basin_acum(flowdir, mask)
    expected = np.array(
        [
            [1, 2, 3, 4, 5, 6, 7],
            [1, 2, 3, 4, 5, 6, 14],
            [1, 2, 3, 4, 5, 6, 21],
            [1, 2, 3, 4, 5, 6, 28],
            [1, 2, 3, 4, 5, 6, 35],
            [1, 2, 3, 4, 5, 6, 42],
            [1, 2, 3, 4, 5, 6, 49],
        ],
        dtype=float,
    )
    assert np.array_equal(acc, expected)
    big_flow = np.zeros((200, 200), dtype=int)
    big_mask = np.ones_like(big_flow, dtype=bool)
    t0 = time.perf_counter()
    basin_acum(big_flow, big_mask)
    assert time.perf_counter() - t0 < 1.0


def test_basin_cut_single() -> None:
    flowdir = np.zeros((7, 7), dtype=int)
    flowdir[:, -1] = 2
    flowdir[-1, -1] = -1
    mask = np.ones_like(flowdir, dtype=bool)
    dem = np.ones_like(flowdir, dtype=float)
    res = basin_cut(dem, (6, 6), flowdir, mask)
    assert res["mask"].sum() == 49


def test_basin_cut_two_catchments() -> None:
    flowdir = np.zeros((7, 7), dtype=int)
    flowdir[:3, :] = 6  # top half drains north
    flowdir[3:, -1] = 2
    flowdir[-1, -1] = -1
    mask = np.ones_like(flowdir, dtype=bool)
    dem = np.ones_like(flowdir, dtype=float)
    res = basin_cut(dem, (6, 6), flowdir, mask)
    assert res["mask"].sum() == 28


def test_basin_cut_edge_outlet() -> None:
    flowdir = np.full((5, 5), 2, dtype=int)
    flowdir[4, :] = 0
    flowdir[4, 2:] = 4
    flowdir[4, 2] = -1
    mask = np.ones_like(flowdir, dtype=bool)
    dem = np.ones_like(flowdir, dtype=float)
    res = basin_cut(dem, (4, 2), flowdir, mask)
    assert res["mask"].sum() == 25
