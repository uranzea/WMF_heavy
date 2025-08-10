"""Basin metrics utilities."""
from __future__ import annotations

from typing import Dict, Any

try:
    import numpy as np
except ModuleNotFoundError:  # pragma: no cover
    np = Any  # type: ignore


def basin_ppalstream_find(*args, **kwargs):  # pragma: no cover - stub
    pass


def basin_ppalstream_cut(*args, **kwargs):  # pragma: no cover - stub
    pass


def basin_ppal_hipsometric(dem_basin: np.ndarray, mask_basin: np.ndarray) -> Dict[str, np.ndarray]:  # pragma: no cover - stub
    return {"dem": dem_basin, "mask": mask_basin}


def basin_time_to_out(flowdir: np.ndarray, dx: float, dy: float, slope: np.ndarray, n_manning: np.ndarray) -> np.ndarray:  # pragma: no cover - stub
    return slope * 0 + 1.0


def geo_hand(dem: np.ndarray, stream: np.ndarray) -> np.ndarray:  # pragma: no cover - stub
    return dem - dem.min()
