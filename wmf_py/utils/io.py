"""Utility functions for raster/NetCDF/CSV IO."""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Mapping, Tuple

import numpy as np
import pandas as pd

try:  # pragma: no cover - optional dependency
    import rasterio
except Exception:  # pragma: no cover - import-time failure
    rasterio = None  # type: ignore

try:  # pragma: no cover - optional dependency
    from netCDF4 import Dataset
except Exception:  # pragma: no cover - import-time failure
    Dataset = None  # type: ignore

Array = np.ndarray
Meta = Mapping[str, object]


class MissingDependencyError(RuntimeError):
    """Raised when an optional dependency is not available."""


# ---------------------------------------------------------------------------
# GeoTIFF helpers
# ---------------------------------------------------------------------------

def read_raster(path: str | Path) -> Tuple[Array, Dict[str, object]]:
    """Read a single-band raster returning the array and metadata.

    Parameters
    ----------
    path: str or Path
        Path to the GeoTIFF file.

    Returns
    -------
    array: np.ndarray
        Array with raster values.
    meta: dict
        Raster metadata including CRS, transform, nodata and dtype.
    """
    if rasterio is None:
        raise MissingDependencyError("rasterio is required to read rasters")

    with rasterio.open(path) as src:
        array = src.read(1, masked=False)
        meta = src.meta.copy()
    return array, meta


def write_geotiff(path: str | Path, array: Array, meta: Mapping[str, object]) -> None:
    """Write a 2D array to GeoTIFF using provided metadata.

    Handles ``np.ma.MaskedArray`` by filling with ``nodata``.
    """
    if rasterio is None:
        raise MissingDependencyError("rasterio is required to write rasters")

    data = array
    if np.ma.isMaskedArray(array):
        nodata = meta.get("nodata")
        fill_value = nodata if nodata is not None else np.nan
        data = array.filled(fill_value)

    meta_out = dict(meta)
    meta_out.update(
        {
            "driver": "GTiff",
            "height": data.shape[0],
            "width": data.shape[1],
            "count": 1,
            "dtype": str(data.dtype),
        }
    )

    with rasterio.open(path, "w", **meta_out) as dst:
        dst.write(data, 1)


# ---------------------------------------------------------------------------
# NetCDF helpers
# ---------------------------------------------------------------------------

def _require_netcdf() -> None:
    if Dataset is None:
        raise MissingDependencyError("netCDF4 is required for NetCDF operations")


def read_netcdf(paths_or_dir: str | Path | Iterable[str | Path], varnames: Iterable[str]) -> Dict[str, Array]:
    """Read variables from NetCDF files.

    Parameters
    ----------
    paths_or_dir: str, Path or iterable of paths
        Path(s) to NetCDF file(s) or directory containing them.
    varnames: iterable of str
        Variable names to load.
    """
    _require_netcdf()

    if isinstance(paths_or_dir, (str, Path)) and Path(paths_or_dir).is_dir():
        paths = [p for p in Path(paths_or_dir).glob("*.nc")]
    elif isinstance(paths_or_dir, (str, Path)):
        paths = [Path(paths_or_dir)]
    else:
        paths = [Path(p) for p in paths_or_dir]

    data: Dict[str, Array] = {}
    missing = set(varnames)

    for p in paths:
        with Dataset(p) as ds:
            for name in list(missing):
                if name in ds.variables:
                    data[name] = np.asarray(ds.variables[name][:])
                    missing.remove(name)
        if not missing:
            break

    if missing:
        raise ValueError(f"Variables not found in NetCDF files: {sorted(missing)}")

    return data


def write_netcdf(path: str | Path, data_dict: Mapping[str, Array], attrs_dict: Mapping[str, object] | None = None) -> None:
    """Write a dictionary of arrays to a NetCDF file.

    All arrays must share the same shape. Supported shapes are ``(y, x)`` or
    ``(time, y, x)``.
    """
    _require_netcdf()

    if not data_dict:
        raise ValueError("data_dict must not be empty")

    shapes = {tuple(v.shape) for v in data_dict.values()}
    if len(shapes) != 1:
        raise ValueError("All arrays in data_dict must have the same shape")
    shape = next(iter(shapes))

    if len(shape) == 2:
        dims = ("y", "x")
    elif len(shape) == 3:
        dims = ("time", "y", "x")
    else:
        raise ValueError("Arrays must be 2D or 3D")

    with Dataset(path, "w", format="NETCDF4") as ds:
        for dim, size in zip(dims, shape):
            ds.createDimension(dim, size)
        for name, array in data_dict.items():
            var = ds.createVariable(name, array.dtype, dims)
            var[:] = array
        if attrs_dict:
            for key, val in attrs_dict.items():
                ds.setncattr(key, val)


# ---------------------------------------------------------------------------
# CSV helper
# ---------------------------------------------------------------------------

def write_csv(path: str | Path, dataframe_or_records) -> None:
    """Write a CSV file from a DataFrame or an iterable of records."""
    if isinstance(dataframe_or_records, pd.DataFrame):
        dataframe_or_records.to_csv(path, index=False)
    else:
        pd.DataFrame(dataframe_or_records).to_csv(path, index=False)
