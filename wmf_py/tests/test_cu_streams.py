import pytest

np = pytest.importorskip("numpy")

from wmf_py.cu_py import streams


def test_stream_find_threshold():
    acum = np.array([[0, 5], [10, 1]])
    stream = streams.stream_find(acum, 3)
    assert stream.sum() == 2
