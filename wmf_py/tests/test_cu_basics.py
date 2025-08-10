import numpy as np

from wmf_py.cu_py.basics import basin_find, basin_map2basin


def test_basin_find_bounds() -> None:
    rc = basin_find(x=60, y=60, xll=0, yll=0, dx=30, dy=30, ncols=3, nrows=3)
    assert rc == (2, 2)


def test_basin_map2basin_mask() -> None:
    v = np.arange(9).reshape(3, 3)
    m = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=bool)
    out = basin_map2basin(v, m)
    assert np.isfinite(out[m]).all()
