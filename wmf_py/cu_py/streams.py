"""Stream network utilities."""
from __future__ import annotations

from collections import deque
from typing import Any, Dict, List, Tuple

from .basics import _D8

try:  # NumPy is optional for documentation builds
    import numpy as np
except ModuleNotFoundError:  # pragma: no cover - used only when NumPy is missing
    np = Any  # type: ignore


# ---------------------------------------------------------------------------
# Basic stream operations
# ---------------------------------------------------------------------------

def stream_find(acum: np.ndarray, threshold: int) -> np.ndarray:
    """Binarise a stream network from an accumulation grid.

    Parameters
    ----------
    acum : ndarray of int
        D8 accumulation map.
    threshold : int
        Minimum accumulation value for a cell to be considered part of the
        stream network.

    Returns
    -------
    ndarray of uint8
        Boolean mask of the detected stream network.
    """

    return (acum >= threshold).astype(np.uint8)


def stream_cut(stream: np.ndarray, mask_basin: np.ndarray) -> np.ndarray:
    """Restrict a stream network to a basin mask."""

    return (stream.astype(np.uint8) & mask_basin.astype(np.uint8)).astype(np.uint8)


def stream_find_to_corr(
    stream: np.ndarray, flowdir: np.ndarray, outlet_rc: Tuple[int, int]
) -> np.ndarray:
    """Remove disconnected stream segments.

    Starting from ``outlet_rc`` the network is explored upstream following the
    inverse of the D8 directions.  Only cells connected to the outlet are
    preserved; all others are set to zero.
    """

    ny, nx = stream.shape
    r0, c0 = outlet_rc
    out = np.zeros_like(stream, dtype=np.uint8)

    if not (0 <= r0 < ny and 0 <= c0 < nx):
        raise ValueError("outlet outside grid")

    q: deque[Tuple[int, int]] = deque()
    visited = np.zeros_like(stream, dtype=bool)

    q.append((r0, c0))
    visited[r0, c0] = True

    while q:
        r, c = q.popleft()
        if stream[r, c]:
            out[r, c] = 1
        for dr, dc in _D8.values():
            rr, cc = r + dr, c + dc
            if not (0 <= rr < ny and 0 <= cc < nx):
                continue
            if visited[rr, cc] or stream[rr, cc] == 0:
                continue
            off = _D8.get(int(flowdir[rr, cc]))
            if off and rr + off[0] == r and cc + off[1] == c:
                visited[rr, cc] = True
                q.append((rr, cc))

    return out


# ---------------------------------------------------------------------------
# Stream extraction from approximate coordinates
# ---------------------------------------------------------------------------

def stream_seed_from_coords(
    acum: np.ndarray,
    x: float,
    y: float,
    xll: float,
    yll: float,
    dx: float,
    dy: float,
    outlet_rc: Tuple[int, int],
    search_radius_cells: int = 4,
    min_threshold: int = 1,
) -> Dict[str, Any]:
    """Locate a seed cell near provided map coordinates.

    The function converts the ``(x, y)`` coordinates to array indices and then
    searches within a square window of size ``(2 * R + 1)`` centred on the
    initial location.  The cell with the largest accumulation value is chosen
    as seed.  Ties are broken by selecting the cell with the smallest Manhattan
    distance to the original location.  Cells with accumulation below
    ``min_threshold`` are ignored.

    Parameters
    ----------
    acum : ndarray of int
        Accumulation grid.
    x, y : float
        Approximate coordinates of the stream in metres.
    xll, yll : float
        Lower-left corner of the grid in metres.
    dx, dy : float
        Cell resolution in ``x`` and ``y`` (metres).
    outlet_rc : tuple of int
        Unused but kept for API compatibility.
    search_radius_cells : int, optional
        Radius of the search window in cells (default ``4``).
    min_threshold : int, optional
        Minimum accumulation value to consider a cell (default ``1``).

    Returns
    -------
    dict
        ``{"seed_rc": (row, col), "acum_val": int}`` with the selected seed
        cell and its accumulation value.
    """

    from .basics import basin_find  # local import to avoid circular

    ny, nx = acum.shape
    r0, c0 = basin_find(x, y, xll, yll, dx, dy, nx, ny)

    rmin = max(r0 - search_radius_cells, 0)
    rmax = min(r0 + search_radius_cells, ny - 1)
    cmin = max(c0 - search_radius_cells, 0)
    cmax = min(c0 + search_radius_cells, nx - 1)

    best_rc = None
    best_val = min_threshold - 1
    best_dist = None

    for rr in range(rmin, rmax + 1):
        for cc in range(cmin, cmax + 1):
            val = int(acum[rr, cc])
            if val < min_threshold:
                continue
            dist = abs(rr - r0) + abs(cc - c0)
            if (
                val > best_val
                or (val == best_val and (best_dist is None or dist < best_dist))
            ):
                best_val = val
                best_dist = dist
                best_rc = (rr, cc)

    if best_rc is None:
        raise ValueError("no valid seed found near provided coordinates")

    return {"seed_rc": best_rc, "acum_val": int(best_val)}


def _connects_to_outlet(
    stream: np.ndarray,
    flowdir: np.ndarray,
    start: Tuple[int, int],
    outlet: Tuple[int, int],
) -> bool:
    """Check if ``start`` reaches ``outlet`` following ``flowdir`` within ``stream``."""

    ny, nx = stream.shape
    r, c = start
    visited = set()
    while True:
        if (r, c) == outlet:
            return True
        if (r, c) in visited:
            return False
        visited.add((r, c))
        if not stream[r, c]:
            return False
        off = _D8.get(int(flowdir[r, c]))
        if not off:
            return False
        r += off[0]
        c += off[1]
        if not (0 <= r < ny and 0 <= c < nx):
            return False


def stream_threshold_nearby(
    acum: np.ndarray,
    seed_rc: Tuple[int, int],
    flowdir: np.ndarray,
    outlet_rc: Tuple[int, int],
    candidate_thresholds: List[int],
) -> int:
    """Select the highest accumulation threshold keeping connectivity."""

    if not candidate_thresholds:
        raise ValueError("empty candidate_thresholds")

    thresholds = sorted(candidate_thresholds, reverse=True)
    feasible_min: int | None = None

    for thr in thresholds:
        stream = acum >= thr
        if not stream[seed_rc]:
            continue
        feasible_min = thr  # at least seed belongs to this threshold
        if stream[outlet_rc] and _connects_to_outlet(stream, flowdir, seed_rc, outlet_rc):
            return thr

    if feasible_min is not None:
        return feasible_min
    raise ValueError("seed cell excluded for all thresholds")


def stream_find_nearby(
    acum: np.ndarray,
    flowdir: np.ndarray,
    x: float,
    y: float,
    xll: float,
    yll: float,
    dx: float,
    dy: float,
    outlet_rc: Tuple[int, int],
    mask_basin: np.ndarray,
    search_radius_cells: int = 4,
    candidate_thresholds: List[int] | None = None,
) -> Dict[str, Any]:
    """Derive a stream network near provided coordinates."""

    if candidate_thresholds is None:
        vals = acum[mask_basin]
        if vals.size == 0:
            raise ValueError("empty basin mask")
        perc = list(range(99, 79, -1))
        candidate_thresholds = sorted(
            {int(np.percentile(vals, p)) for p in perc}, reverse=True
        )

    seed_info = stream_seed_from_coords(
        acum,
        x,
        y,
        xll,
        yll,
        dx,
        dy,
        outlet_rc,
        search_radius_cells=search_radius_cells,
    )
    seed_rc = seed_info["seed_rc"]

    thr = stream_threshold_nearby(acum, seed_rc, flowdir, outlet_rc, candidate_thresholds)
    stream = stream_find(acum, thr)
    stream = stream_cut(stream, mask_basin)
    stream = stream_find_to_corr(stream, flowdir, outlet_rc)

    return {"stream": stream.astype(bool), "threshold": thr, "seed_rc": seed_rc}


# ---------------------------------------------------------------------------
# Network nodes and segmentation
# ---------------------------------------------------------------------------

def basin_netxy_find(
    stream: np.ndarray, flowdir: np.ndarray, outlet_rc: Tuple[int, int]
) -> Tuple[np.ndarray, np.ndarray, Dict[int, Tuple[int, int]]]:
    """Detect network nodes.

    Parameters
    ----------
    stream : ndarray of bool or int
        Stream network mask.
    flowdir : ndarray of int
        D8 flow direction codes using the internal convention.
    outlet_rc : tuple of int
        Row/column indices of the basin outlet which is always included as a
        node.

    Returns
    -------
    rows, cols, nodes : tuple
        ``rows`` and ``cols`` contain the coordinates of all stream cells.
        ``nodes`` is a dictionary ``{id: (row, col)}`` with the detected
        network nodes.
    """

    ny, nx = stream.shape
    rows, cols = np.nonzero(stream)
    nodes: Dict[int, Tuple[int, int]] = {}
    node_map = np.zeros_like(stream, dtype=int)
    node_id = 1

    for r, c in zip(rows, cols):
        indeg = 0
        for dr, dc in _D8.values():
            rr, cc = r + dr, c + dc
            if 0 <= rr < ny and 0 <= cc < nx and stream[rr, cc]:
                off = _D8.get(int(flowdir[rr, cc]))
                if off and rr + off[0] == r and cc + off[1] == c:
                    indeg += 1

        off = _D8.get(int(flowdir[r, c]))
        outdeg = 0
        if off:
            rr, cc = r + off[0], c + off[1]
            if 0 <= rr < ny and 0 <= cc < nx and stream[rr, cc]:
                outdeg = 1

        deg = indeg + outdeg
        if (r, c) == outlet_rc or deg != 2:
            node_map[r, c] = node_id
            nodes[node_id] = (r, c)
            node_id += 1

    return rows, cols, nodes


def basin_netxy_cut(
    stream: np.ndarray, nodes: Dict[int, Tuple[int, int]], flowdir: np.ndarray
) -> List[Dict[str, Any]]:
    """Segment the stream network between nodes.

    Returns a list of edges where each edge is a dictionary with the keys
    ``id``, ``node_u`` (upstream node), ``node_v`` (downstream node),
    ``length`` (number of cells) and ``cells`` (list of ``(r, c)`` tuples).
    """

    ny, nx = stream.shape
    node_map = np.zeros_like(stream, dtype=int)
    for nid, (r, c) in nodes.items():
        node_map[r, c] = nid

    edges: List[Dict[str, Any]] = []
    edge_id = 1

    for nid, (r, c) in nodes.items():
        off = _D8.get(int(flowdir[r, c]))
        if not off:
            continue
        rr, cc = r + off[0], c + off[1]
        if not (0 <= rr < ny and 0 <= cc < nx) or not stream[rr, cc]:
            continue

        path = [(r, c)]
        while True:
            path.append((rr, cc))
            if node_map[rr, cc] and node_map[rr, cc] != nid:
                edges.append(
                    {
                        "id": edge_id,
                        "node_u": nid,
                        "node_v": node_map[rr, cc],
                        "length": len(path),
                        "cells": path.copy(),
                    }
                )
                edge_id += 1
                break

            off = _D8.get(int(flowdir[rr, cc]))
            if not off:
                break
            rr += off[0]
            cc += off[1]

    return edges


# Re-export sub-basin utilities from :mod:`metrics` to keep compatibility
from .metrics import (  # noqa: E402  (imported late to avoid circular refs)
    basin_subbasin_cut,
    basin_subbasin_find,
    basin_subbasin_horton,
    basin_subbasin_map2subbasin,
    basin_subbasin_nod,
)

