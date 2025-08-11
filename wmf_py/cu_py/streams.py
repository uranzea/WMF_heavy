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



def stream_seed_from_coords(*args, **kwargs):
    """Placeholder for legacy functionality."""
    raise NotImplementedError("stream_seed_from_coords is not implemented")


def stream_threshold_nearby(*args, **kwargs):
    """Placeholder for legacy functionality."""
    raise NotImplementedError("stream_threshold_nearby is not implemented")


def stream_find_nearby(*args, **kwargs):
    """Placeholder for legacy functionality."""
    raise NotImplementedError("stream_find_nearby is not implemented")


def hydro_distance_and_receiver(*args, **kwargs):
    """Placeholder for legacy functionality."""
    raise NotImplementedError("hydro_distance_and_receiver is not implemented")
