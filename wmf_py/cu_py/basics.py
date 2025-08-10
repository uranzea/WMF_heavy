"""Basic basin processing utilities.

This module provides minimal implementations of a few basin-related
functions required by :mod:`wmfv2`.  The algorithms are deliberately
simple and mostly act as placeholders until the full functionality from
the original Fortran code is ported.
"""

from __future__ import annotations

from typing import Any, Dict, Tuple
from collections import deque

try:  # NumPy optional for stubs
    import numpy as np
except ModuleNotFoundError:  # pragma: no cover
    np = Any  # type: ignore

# Mapping of D8 codes to row/column offsets following a common convention
_D8 = {
    1: (0, 1),  # East
    2: (-1, 1),  # North East
    3: (-1, 0),  # North
    4: (-1, -1),  # North West
    5: (0, -1),  # West
    6: (1, -1),  # South West
    7: (1, 0),  # South
    8: (1, 1),  # South East
}


def dir_reclass_opentopo(flowdir: np.ndarray) -> np.ndarray:
    """Reclassify OpenTopo flow directions.

    Currently the function returns the input unchanged because the
    internal convention follows the D8 codes used by OpenTopo.  A full
    mapping table may be added in the future.
    """

    return flowdir.copy()


def dir_reclass_rwatershed(flowdir: np.ndarray) -> np.ndarray:
    """Reclassify ``r.watershed`` flow directions.

    The present MVP assumes the incoming array already uses the internal
    D8 numbering and simply returns a copy.  A detailed reclassification
    may be implemented later.  See :func:`dir_reclass_opentopo`.
    """

    return flowdir.copy()


def basin_acum(flowdir: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Compute upstream cell accumulation following D8 directions.

    Parameters
    ----------
    flowdir : ndarray
        Array of D8 flow directions encoded with integers 1--8.
    mask : ndarray of bool
        Basin mask where ``True`` denotes active cells.

    Returns
    -------
    ndarray
        Accumulated upstream cell counts.  Cells outside the mask return
        zero.  Cells involved in direction loops retain zero.
    """

    ny, nx = flowdir.shape
    acc = np.zeros_like(flowdir, dtype=float)
    receivers = np.full((ny, nx, 2), -1, dtype=int)
    indegree = np.zeros((ny, nx), dtype=int)

    for r in range(ny):
        for c in range(nx):
            if not mask[r, c]:
                continue
            off = _D8.get(int(flowdir[r, c]))
            if not off:
                continue
            rr, cc = r + off[0], c + off[1]
            if 0 <= rr < ny and 0 <= cc < nx and mask[rr, cc]:
                receivers[r, c] = (rr, cc)
                indegree[rr, cc] += 1

    q: deque[Tuple[int, int]] = deque()
    for r in range(ny):
        for c in range(nx):
            if mask[r, c] and indegree[r, c] == 0:
                acc[r, c] = 1.0
                q.append((r, c))

    while q:
        r, c = q.popleft()
        rr, cc = receivers[r, c]
        if rr == -1:
            continue
        acc[rr, cc] += acc[r, c]
        indegree[rr, cc] -= 1
        if indegree[rr, cc] == 0:
            if acc[rr, cc] == 0:
                acc[rr, cc] = 1.0
            q.append((rr, cc))

    return acc


def basin_find(
    x: float,
    y: float,
    xll: float,
    yll: float,
    dx: float,
    dy: float,
    ncols: int,
    nrows: int,
) -> Tuple[int, int]:
    """Locate array indices given map coordinates.

    Coordinates are given in metres and refer to the lower-left corner of
    the raster (``xll``, ``yll``).  The returned indices follow the
    ``(row, col)`` convention with the origin at the lower-left cell.  A
    :class:`ValueError` is raised if the point lies outside the raster
    bounds.
    """

    col = int((x - xll) / dx)
    row = int((y - yll) / dy)
    if not (0 <= col < ncols and 0 <= row < nrows):
        raise ValueError("point outside raster bounds")
    return row, col


def basin_cut(
    dem: np.ndarray,
    outlet_rc: Tuple[int, int],
    flowdir: np.ndarray,
    mask: np.ndarray,
) -> Dict[str, np.ndarray]:
    """Apply a basin mask to DEM and flow directions.

    This is a minimal implementation that simply clips the provided
    arrays with the boolean ``mask``.  It does **not** perform an actual
    delineation following ``flowdir``; that behaviour is left for a future
    version.
    """

    dem_cut = np.where(mask, dem, np.nan)
    flowdir_cut = np.where(mask, flowdir, 0)
    return {"dem": dem_cut, "flowdir": flowdir_cut, "mask": mask}


def basin_2map(
    values: np.ndarray,
    nrows: int,
    ncols: int,
    nodata: float,
    xll: float,
    yll: float,
    dx: float,
    dy: float,
) -> np.ndarray:
    """Map basin values to a regular grid.

    The current placeholder simply reshapes ``values`` to ``(nrows, ncols)``.
    The coordinate arguments are accepted for API compatibility and will
    be used once a proper basin↔map conversion is implemented.
    """

    return values.reshape((nrows, ncols))


def basin_map2basin(values_map: np.ndarray, mask_basin: np.ndarray) -> np.ndarray:
    """Extract basin values from a map using a boolean mask."""

    return values_map[mask_basin.astype(bool)]
