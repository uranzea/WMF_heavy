import numpy as np
import time

from wmf_py.cu_py.basics import basin_find, basin_map2basin, basin_acum


def test_basin_find_bounds() -> None:
    rc = basin_find(x=60, y=60, xll=0, yll=0, dx=30, dy=30, ncols=3, nrows=3)
    assert rc == (2, 2)


def test_basin_map2basin_mask() -> None:
    v = np.arange(9).reshape(3, 3)
    m = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=bool)
    out = basin_map2basin(v, m)
    assert np.isfinite(out[m]).all()


def test_basin_acum_d8() -> None:
    flowdir = np.ones((7, 7), dtype=int)
    flowdir[:, -1] = 7  # right column flows south
    flowdir[-1, -1] = 0  # outlet
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
    # Simple performance sanity check
    big_flow = np.ones((200, 200), dtype=int)
    big_mask = np.ones_like(big_flow, dtype=bool)
    t0 = time.perf_counter()
    basin_acum(big_flow, big_mask)
    assert time.perf_counter() - t0 < 1.0
