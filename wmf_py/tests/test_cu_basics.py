import pytest

np = pytest.importorskip("numpy")

from wmf_py.cu_py import basics


def test_dir_reclass_identity():
    arr = np.arange(4).reshape(2, 2)
    out = basics.dir_reclass_opentopo(arr)
    assert np.array_equal(out, arr)
