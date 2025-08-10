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
    idx_basin_to_map: np.ndarray,
    fill: float = np.nan,
) -> np.ndarray:  # pragma: no cover - stub
    """Map basin values to a full raster shape."""
    out = np.full(map_shape, fill)
    out.flat[idx_basin_to_map] = values_basin
    return out
