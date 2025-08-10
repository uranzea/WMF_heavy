import numpy as np
import pytest

from wmf_py.cu_py.basics import basin_2map, basin_map2basin
from wmf_py.utils.transform import Transform_Basin2Map


def test_idx_basin_to_map_checkerboard() -> None:
    mask = np.indices((5, 5)).sum(axis=0) % 2 == 0
    idx = np.flatnonzero(mask.ravel(order="C"))
    assert idx.size == mask.sum() == 13


def test_basin_2map_and_map2basin() -> None:
    map_shape = (50, 50)
    idx = np.array([0, 10, 101, 555, 999, 1234, 2345, 3456, 4321, 4999])
    values = np.arange(idx.size)
    nodata = -1.0
    raster = basin_2map(values, idx, map_shape, nodata)
    assert raster.shape == map_shape
    assert np.all(raster.ravel()[idx] == values)
    assert (raster != nodata).sum() == values.size
    back = basin_map2basin(raster, idx, (idx.size,))
    assert np.array_equal(back, values)


def test_transform_idempotent() -> None:
    idx = np.array([0, 2, 4, 6])
    values = np.array([[1, 2], [3, 4]])
    map_shape = (3, 3)
    nodata = -9999.0
    m = basin_2map(values, idx, map_shape, nodata)
    back = basin_map2basin(m, idx, values.shape)
    assert np.array_equal(back, values)
    # also test Transform_Basin2Map wrapper
    out = Transform_Basin2Map(values, map_shape, idx, fill=nodata)
    assert np.array_equal(out, m)
