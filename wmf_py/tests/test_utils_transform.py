import numpy as np
import pytest

from wmf_py.utils.transform import Transform_Basin2Map


def test_transform_basin2map_requires_idx() -> None:
    with pytest.raises(NotImplementedError):
        Transform_Basin2Map(np.array([1]), (1, 1), None)


def test_transform_basin2map_basic() -> None:
    values = np.array([1, 2, 3, 4])
    idx = np.array([0, 4, 6, 8])  # map to a 3x3 grid
    out = Transform_Basin2Map(values, (3, 3), idx, fill=-1)
    assert out[0, 0] == 1
    assert out[1, 1] == 2
    assert out[2, 0] == 3
    assert out[2, 2] == 4
