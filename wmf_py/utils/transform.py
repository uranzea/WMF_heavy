"""Utility transforms."""

from __future__ import annotations

from typing import Any

try:
    import numpy as np
except ModuleNotFoundError:  # pragma: no cover
    np = Any  # type: ignore


def Transform_Basin2Map(
    values_basin: np.ndarray,
    map_shape: tuple,
    idx_basin_to_map: np.ndarray | None,
    fill: float = np.nan,
) -> np.ndarray:
    """Map basin values to a full raster.

    Parameters
    ----------
    values_basin : ndarray
        Values defined on basin cells.
    map_shape : tuple
        Shape ``(nrows, ncols)`` of the target map.
    idx_basin_to_map : ndarray
        Linear or ``(row, col)`` indices that locate each basin cell within
        the target map.  The argument is mandatory; passing ``None`` raises
        :class:`NotImplementedError` as a reminder that real indexing is
        required for the transform to operate.
    fill : float, optional
        Value used to fill cells outside the basin (default ``np.nan``).
    """

    if idx_basin_to_map is None:
        raise NotImplementedError("idx_basin_to_map must be provided")

    out = np.full(map_shape, fill, dtype=values_basin.dtype)
    out_flat = out.ravel()
    out_flat[idx_basin_to_map] = values_basin
    return out
