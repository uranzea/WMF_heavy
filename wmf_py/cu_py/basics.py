"""Basic basin processing utilities."""
from __future__ import annotations

from typing import Dict, Tuple, Any

try:  # NumPy optional for stubs
    import numpy as np
except ModuleNotFoundError:  # pragma: no cover
    np = Any  # type: ignore


def dir_reclass_opentopo(flowdir: np.ndarray) -> np.ndarray:  # pragma: no cover - stub
    """Reclassify OpenTopo flow directions."""
    return flowdir


def dir_reclass_rwatershed(flowdir: np.ndarray) -> np.ndarray:  # pragma: no cover - stub
    """Reclassify r.watershed flow directions."""
    return flowdir


def basin_acum(flowdir: np.ndarray, mask: np.ndarray) -> np.ndarray:  # pragma: no cover - stub
    """Compute upstream cell accumulation."""
    return mask


def basin_find(
    x: float,
    y: float,
    xll: float,
    yll: float,
    dx: float,
    dy: float,
    ncols: int,
    nrows: int,
) -> Tuple[int, int]:  # pragma: no cover - stub
    """Locate array indices given coordinates."""
    return 0, 0


def basin_cut(
    dem: np.ndarray,
    outlet_rc: Tuple[int, int],
    flowdir: np.ndarray,
    mask: np.ndarray,
) -> Dict[str, np.ndarray]:  # pragma: no cover - stub
    """Cut a basin mask from a DEM."""
    return {"dem": dem, "mask": mask}


def basin_2map(
    values: np.ndarray,
    nrows: int,
    ncols: int,
    nodata: float,
    xll: float,
    yll: float,
    dx: float,
    dy: float,
) -> np.ndarray:  # pragma: no cover - stub
    """Map basin values to a regular grid."""
    return values.reshape((nrows, ncols))


def basin_map2basin(values_map: np.ndarray, mask_basin: np.ndarray) -> np.ndarray:  # pragma: no cover - stub
    """Extract basin values from a map."""
    return values_map[mask_basin.astype(bool)]
