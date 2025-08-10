"""Stream network utilities."""
from __future__ import annotations

from typing import Tuple, Any

try:
    import numpy as np
except ModuleNotFoundError:  # pragma: no cover
    np = Any  # type: ignore


def stream_find(acum: np.ndarray, threshold: int) -> np.ndarray:  # pragma: no cover - stub
    """Identify stream cells based on accumulation threshold."""
    return (acum > threshold).astype(int)


def stream_cut(stream: np.ndarray, mask_basin: np.ndarray) -> np.ndarray:  # pragma: no cover - stub
    """Mask the stream network within a basin."""
    return stream * mask_basin


def stream_find_to_corr(stream: np.ndarray, flowdir: np.ndarray) -> np.ndarray:  # pragma: no cover - stub
    """Find stream cells to correct."""
    return stream


def basin_netxy_find(flowdir: np.ndarray, stream: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:  # pragma: no cover - stub
    """Return placeholder arrays representing network coordinates."""
    idx = np.argwhere(stream > 0)
    return idx[:, 0], idx[:, 1]


def basin_netxy_cut(*args, **kwargs):  # pragma: no cover - stub
    pass


def basin_subbasin_find(*args, **kwargs):  # pragma: no cover - stub
    pass


def basin_subbasin_cut(*args, **kwargs):  # pragma: no cover - stub
    pass


def basin_subbasin_map2subbasin(*args, **kwargs):  # pragma: no cover - stub
    pass


def basin_subbasin_nod(*args, **kwargs):  # pragma: no cover - stub
    pass


def basin_subbasin_horton(*args, **kwargs):  # pragma: no cover - stub
    pass
