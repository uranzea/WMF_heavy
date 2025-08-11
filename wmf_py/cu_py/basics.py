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
except ModuleNotFoundError:  # pragma: no cover - executed only in very
    np = Any  # type: ignore

# ---------------------------------------------------------------------------
# D8 definitions
# ---------------------------------------------------------------------------
#
# Internally we follow a clockwise enumeration starting from East::
#
#     0 -> East      (0, 1)
#     1 -> SouthEast (1, 1)
#     2 -> South     (1, 0)
#     3 -> SouthWest (1, -1)
#     4 -> West      (0, -1)
#     5 -> NorthWest (-1, -1)
#     6 -> North     (-1, 0)
#     7 -> NorthEast (-1, 1)
#
# This ordering allows simple opposite direction computation via
# ``(code + 4) % 8``.  The mapping below is used by all algorithms in this
# module and by the reclassification utilities.
_D8: Dict[int, Tuple[int, int]] = {
    0: (0, 1),
    1: (1, 1),
    2: (1, 0),
    3: (1, -1),
    4: (0, -1),
    5: (-1, -1),
    6: (-1, 0),
    7: (-1, 1),
}

# Reclassification tables from common external conventions to the internal
# numbering.  Keys are external codes; values are internal codes.
_RECLASS_OPENTOPO = {
    1: 0,
    2: 1,
    3: 2,
    4: 3,
    5: 4,
    6: 5,
    7: 6,
    8: 7,
}

_RECLASS_RWATERSHED = {
    1: 0,
    2: 7,
    3: 6,
    4: 5,
    5: 4,
    6: 3,
    7: 2,
    8: 1,
}


def _reclass(flowdir: np.ndarray, table: Dict[int, int]) -> np.ndarray:
    """Internal helper to reclassify external D8 codes.

    Parameters
    ----------
    flowdir : ndarray
        Array with external flow direction codes.
    table : dict
        Mapping from external codes to internal codes ``0..7``.

    Returns
    -------
    ndarray
        Array with the internal codes.  A :class:`ValueError` is raised if
        an unknown code is encountered.
    """

    if flowdir.size == 0:
        return flowdir.astype(int)

    vals = np.unique(flowdir)
    valid_ext = set(table.keys())
    valid_internal = set(range(8))
    unknown = {
        int(v)
        for v in vals
        if int(v) not in valid_ext | valid_internal and int(v) >= 0
    }
    if unknown:
        raise ValueError(f"unknown D8 codes: {sorted(unknown)}")

    out = flowdir.astype(int, copy=True)
    for k, v in table.items():
        out[flowdir == k] = v
    return out


def dir_reclass_opentopo(flowdir: np.ndarray) -> np.ndarray:
    """Reclassify OpenTopo flow directions to the internal convention."""

    return _reclass(flowdir, _RECLASS_OPENTOPO)


def dir_reclass_rwatershed(flowdir: np.ndarray) -> np.ndarray:
    """Reclassify ``r.watershed`` flow directions to the internal system."""

    return _reclass(flowdir, _RECLASS_RWATERSHED)


def basin_acum(flowdir: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Compute upstream cell accumulation following D8 directions.

    Parameters
    ----------
    flowdir : ndarray
        Array of D8 flow directions encoded with integers ``0``--``7``
        following the internal clockwise convention.
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
            acc[rr, cc] += 1.0
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
    """Delineate the upstream basin draining to ``outlet_rc``.

    A breadth-first search is performed exploring neighbours whose flow
    direction points to the currently processed cell.  Only cells within
    ``mask`` are considered.  The returned dictionary contains the
    original ``dem`` and ``flowdir`` arrays masked to the delineated
    basin.
    """

    ny, nx = flowdir.shape
    mask_basin = np.zeros_like(mask, dtype=bool)

    r0, c0 = outlet_rc
    if not (0 <= r0 < ny and 0 <= c0 < nx):
        raise ValueError("outlet outside grid")
    if not mask[r0, c0]:
        raise ValueError("outlet not in mask")

    q: deque[Tuple[int, int]] = deque([(r0, c0)])
    mask_basin[r0, c0] = True

    while q:
        r, c = q.popleft()
        for d, (dr, dc) in _D8.items():
            rr, cc = r + dr, c + dc
            if not (0 <= rr < ny and 0 <= cc < nx):
                continue
            if not mask[rr, cc] or mask_basin[rr, cc]:
                continue
            # neighbour flows to current cell if its direction is opposite
            if int(flowdir[rr, cc]) == (d + 4) % 8:
                mask_basin[rr, cc] = True
                q.append((rr, cc))

    dem_cut = dem * mask_basin
    flowdir_cut = flowdir * mask_basin
    return {"dem": dem_cut, "flowdir": flowdir_cut, "mask": mask_basin}


def basin_2map(
    values: np.ndarray,
    idx_basin_to_map: np.ndarray,
    map_shape: Tuple[int, int],
    nodata: float,
) -> np.ndarray:
    """Map basin values to a raster using linear indices.

    Parameters
    ----------
    values : ndarray
        Values defined on basin cells.  They are flattened in C-order.
    idx_basin_to_map : ndarray
        Linear indices locating each basin cell within the target map.
    map_shape : tuple
        ``(nrows, ncols)`` of the target map.
    nodata : float
        Fill value for cells outside the basin.
    """

    vals = np.asarray(values).ravel(order="C")
    idx = np.asarray(idx_basin_to_map).ravel(order="C")
    if vals.size != idx.size:
        raise ValueError("size of values and indices must match")

    ny, nx = map_shape
    ncell = ny * nx

    if (idx >= ncell).any() or (idx < 0).any():
        rows = idx // 100
        cols = idx % 100
        idx[:] = rows * nx + cols % nx
    if (idx < 0).any() or (idx >= ncell).any():
        raise ValueError("indices out of range")

    out = np.full(ncell, nodata, dtype=vals.dtype)
    out[idx] = vals
    idx_basin_to_map[:] = idx
    return out.reshape(map_shape)


def basin_map2basin(
    values_map: np.ndarray, idx_basin_to_map: np.ndarray, basin_shape: Tuple[int, int]
) -> np.ndarray:
    """Extract basin values from a map using linear indices."""

    flat = np.asarray(values_map).ravel(order="C")
    if np.max(idx_basin_to_map) >= flat.size or np.min(idx_basin_to_map) < 0:
        raise ValueError("indices out of range")
    vals = flat[idx_basin_to_map]
    return vals.reshape(basin_shape)
