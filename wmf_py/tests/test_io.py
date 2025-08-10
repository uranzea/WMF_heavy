import numpy as np
import pandas as pd
import pytest
from pathlib import Path

from wmf_py.utils.io import (
    read_raster,
    write_geotiff,
    read_netcdf,
    write_netcdf,
    write_csv,
)

rasterio = pytest.importorskip("rasterio")
netCDF4 = pytest.importorskip("netCDF4")


def test_geotiff_roundtrip(tmp_path: Path) -> None:
    arr = np.arange(9, dtype=np.float32).reshape(3, 3)
    meta = {
        "crs": "EPSG:4326",
        "transform": rasterio.transform.from_origin(0, 3, 1, 1),
        "nodata": -9999.0,
        "dtype": "float32",
    }
    path = tmp_path / "test.tif"
    write_geotiff(path, arr, meta)
    out, out_meta = read_raster(path)
    assert np.allclose(arr, out)
    assert out_meta["crs"].to_string() == "EPSG:4326"


def test_netcdf_roundtrip(tmp_path: Path) -> None:
    arr = np.random.rand(2, 3, 4).astype("float32")
    path = tmp_path / "test.nc"
    write_netcdf(path, {"var": arr}, {"units": "m"})
    data = read_netcdf(path, ["var"])
    assert np.allclose(arr, data["var"])


def test_write_csv(tmp_path: Path) -> None:
    df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
    path = tmp_path / "out.csv"
    write_csv(path, df)
    loaded = pd.read_csv(path)
    pd.testing.assert_frame_equal(df, loaded)
