"""Basin metrics utilities."""
from __future__ import annotations

from collections import deque
from typing import Any, Dict, List, Tuple

from .basics import _D8

try:
    import numpy as np
except ModuleNotFoundError:  # pragma: no cover
    np = Any  # type: ignore


def hypsometric_points(dem: np.ndarray, mask_basin: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return sorted elevations and cumulative area fraction for a basin.

    Parameters
    ----------
    dem : ndarray of float
        Digital elevation model in metres.
    mask_basin : ndarray of bool
        Mask delineating the basin. ``True`` denotes active cells.

    Returns
    -------
    tuple of ndarray
        ``(elev_sorted, area_cum_frac)`` where ``elev_sorted`` contains the
        basin elevations in ascending order and ``area_cum_frac`` the
        corresponding cumulative area fraction in the range ``0..1``.
    """

    dem_vals = dem[mask_basin]
    dem_vals = dem_vals[~np.isnan(dem_vals)]
    if dem_vals.size == 0:
        return np.array([]), np.array([])

    elev_sorted = np.sort(dem_vals)
    n = elev_sorted.size
    area_cum_frac = np.arange(1, n + 1, dtype=float) / float(n)
    return elev_sorted, area_cum_frac


def hypsometric_curve(dem: np.ndarray, mask_basin: np.ndarray, nbins: int = 100) -> dict:
    """Compute a binned hypsometric curve for the basin.

    Parameters
    ----------
    dem : ndarray
        Digital elevation model in metres.
    mask_basin : ndarray of bool
        Basin mask where ``True`` denotes active cells.
    nbins : int, optional
        Number of bins/percentiles for the curve, by default ``100``.

    Returns
    -------
    dict
        Dictionary containing ``elev_bins`` (elevation percentiles),
        ``area_frac`` (area fraction), ``min`` and ``max`` elevations and the
        mean elevation (``mean_elev``).
    """

    dem_vals = dem[mask_basin]
    dem_vals = dem_vals[~np.isnan(dem_vals)]
    if dem_vals.size == 0:
        empty = np.array([])
        return {
            "elev_bins": empty,
            "area_frac": empty,
            "min": float("nan"),
            "max": float("nan"),
            "mean_elev": float("nan"),
        }

    elev_bins = np.percentile(dem_vals, np.linspace(0, 100, nbins))
    area_frac = np.linspace(0.0, 1.0, nbins)
    return {
        "elev_bins": elev_bins,
        "area_frac": area_frac,
        "min": float(dem_vals.min()),
        "max": float(dem_vals.max()),
        "mean_elev": float(dem_vals.mean()),
    }


def travel_velocity_manning(slope_cell: float, n: float, R: float) -> float:
    """Compute flow velocity using Manning's equation.

    Parameters
    ----------
    slope_cell : float
        Slope (``m/m``).  Values are expected to be non-negative.
    n : float
        Manning roughness coefficient.
    R : float
        Hydraulic radius (approximated by water depth) in metres.

    Returns
    -------
    float
        Velocity in ``m/s``.
    """

    return (1.0 / n) * (R ** (2.0 / 3.0)) * (slope_cell ** 0.5)


def time_to_outlet(
    flowdir: np.ndarray,
    slope: np.ndarray | None,
    mask_basin: np.ndarray,
    stream_mask: np.ndarray,
    dx: float,
    dy: float,
    params: dict,
    outlet_rc: tuple[int, int],
) -> np.ndarray:
    """Compute travel time from each cell to the basin outlet.

    The computation follows D8 directions using velocities derived from
    Manning's equation.  Times are accumulated downstream so that each cell
    reports the total travel time to the specified ``outlet_rc``.
    """

    if slope is None:
        raise ValueError("slope array is required")
    if not mask_basin[outlet_rc]:
        raise ValueError("outlet outside basin mask")

    ny, nx = flowdir.shape
    time_map = np.full((ny, nx), np.nan, dtype=float)

    donors = [[[] for _ in range(nx)] for _ in range(ny)]
    for r in range(ny):
        for c in range(nx):
            if not mask_basin[r, c]:
                continue
            off = _D8.get(int(flowdir[r, c]))
            if not off:
                continue
            rr, cc = r + off[0], c + off[1]
            if 0 <= rr < ny and 0 <= cc < nx and mask_basin[rr, cc]:
                donors[rr][cc].append((r, c))

    S_min = params.get("S_min", 1e-4)
    n_over = params.get("n_overland")
    h_over = params.get("h_overland_m")
    n_ch = params.get("n_channel")
    h_ch = params.get("h_channel_m")

    diag_len = float(np.hypot(dx, dy))

    q: deque[tuple[int, int]] = deque([outlet_rc])
    time_map[outlet_rc] = 0.0
    visited = np.zeros((ny, nx), dtype=bool)
    visited[outlet_rc] = True

    while q:
        r, c = q.popleft()
        base_time = time_map[r, c]
        for dr, dc in donors[r][c]:
            if visited[dr, dc]:
                continue
            S = max(float(slope[dr, dc]), S_min)
            if stream_mask[dr, dc]:
                v = travel_velocity_manning(S, n_ch, h_ch)
            else:
                v = travel_velocity_manning(S, n_over, h_over)
            off = _D8.get(int(flowdir[dr, dc]))
            if not off:
                continue
            if off[0] != 0 and off[1] != 0:
                L = diag_len
            elif off[0] == 0:
                L = dx
            else:
                L = dy
            t = L / v if v > 0 else np.nan
            if np.isnan(t):
                continue
            time_map[dr, dc] = base_time + t
            visited[dr, dc] = True
            q.append((dr, dc))

    time_map[~mask_basin] = np.nan
    return time_map


def hydro_distance_and_receiver(
    flowdir: np.ndarray,
    stream_mask: np.ndarray,
    mask_basin: np.ndarray,
    dx: float,
    dy: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Hydrologic distance to stream and receiving stream cell.

    Parameters
    ----------
    flowdir : ndarray of int
        D8 flow direction codes using the internal convention.
    stream_mask : ndarray of bool
        Boolean mask of stream cells.
    mask_basin : ndarray of bool
        Basin mask where ``True`` denotes active cells.
    dx, dy : float
        Cell resolution in metres.

    Returns
    -------
    tuple of ndarray
        ``(hdnd_model, a_quien)`` where ``hdnd_model`` contains the accumulated
        downstream distance to the first stream cell and ``a_quien`` the linear
        index of that receiving stream cell.  ``hdnd_model`` is ``nan`` outside
        the basin and zero on stream cells; ``a_quien`` is ``-1`` outside the
        basin and equals the cell's own index on stream cells.
    """

    ny, nx = flowdir.shape

    donors = [[[] for _ in range(nx)] for _ in range(ny)]
    for r in range(ny):
        for c in range(nx):
            if not mask_basin[r, c]:
                continue
            off = _D8.get(int(flowdir[r, c]))
            if not off:
                continue
            rr, cc = r + off[0], c + off[1]
            if 0 <= rr < ny and 0 <= cc < nx and mask_basin[rr, cc]:
                donors[rr][cc].append((r, c))

    diag_len = float(np.hypot(dx, dy))
    hdnd = np.full((ny, nx), np.nan, dtype=float)
    aquien = np.full((ny, nx), -1, dtype=int)

    q: deque[tuple[int, int]] = deque()
    visited = np.zeros((ny, nx), dtype=bool)
    for r in range(ny):
        for c in range(nx):
            if stream_mask[r, c] and mask_basin[r, c]:
                hdnd[r, c] = 0.0
                aquien[r, c] = r * nx + c
                q.append((r, c))
                visited[r, c] = True

    while q:
        r, c = q.popleft()
        for dr, dc in donors[r][c]:
            if visited[dr, dc]:
                continue
            off = _D8.get(int(flowdir[dr, dc]))
            if not off:
                continue
            if off[0] != 0 and off[1] != 0:
                L = diag_len
            elif off[0] == 0:
                L = dx
            else:
                L = dy
            hdnd[dr, dc] = hdnd[r, c] + L
            aquien[dr, dc] = aquien[r, c]
            visited[dr, dc] = True
            q.append((dr, dc))

    hdnd[~mask_basin] = np.nan
    aquien[~mask_basin] = -1
    return hdnd, aquien


def hand(
    dem: np.ndarray,
    flowdir: np.ndarray,
    stream_mask: np.ndarray,
    mask_basin: np.ndarray,
) -> np.ndarray:
    """Compute the Height Above Nearest Drainage (HAND) metric."""

    ny, nx = dem.shape
    base = np.full_like(dem, np.nan, dtype=float)
    base[stream_mask & mask_basin] = dem[stream_mask & mask_basin]

    receivers = np.full((ny, nx, 2), -1, dtype=int)
    indegree = np.zeros((ny, nx), dtype=int)
    for r in range(ny):
        for c in range(nx):
            if not mask_basin[r, c]:
                continue
            off = _D8.get(int(flowdir[r, c]))
            if not off:
                continue
            rr, cc = r + off[0], c + off[1]
            if 0 <= rr < ny and 0 <= cc < nx and mask_basin[rr, cc]:
                receivers[r, c] = (rr, cc)
                indegree[rr, cc] += 1

    order: list[tuple[int, int]] = []
    q: deque[tuple[int, int]] = deque()
    for r in range(ny):
        for c in range(nx):
            if mask_basin[r, c] and indegree[r, c] == 0:
                q.append((r, c))
    while q:
        r, c = q.popleft()
        order.append((r, c))
        rr, cc = receivers[r, c]
        if rr == -1:
            continue
        indegree[rr, cc] -= 1
        if indegree[rr, cc] == 0:
            q.append((rr, cc))

    for r, c in reversed(order):
        if not mask_basin[r, c] or stream_mask[r, c]:
            continue
        rr, cc = receivers[r, c]
        if rr != -1:
            base[r, c] = base[rr, cc]

    out = dem.astype(float) - base
    out[~mask_basin] = np.nan
    out = np.maximum(out, 0.0)
    return out


# ---------------------------------------------------------------------------
# Principal stream utilities
# ---------------------------------------------------------------------------

# ``_ppal_stream_profile`` stores the most recently computed principal stream
# profile.  This mimics the behaviour of the original Fortran routines where
# ``basin_ppalstream_find`` populated a module level array which was later
# retrieved by ``basin_ppalstream_cut``.
_ppal_stream_profile: np.ndarray | None = None


def basin_ppalstream_find(
    stream: np.ndarray,
    flowdir: np.ndarray,
    dem: np.ndarray,
    dx: float,
    dy: float,
    outlet_rc: Tuple[int, int],
) -> Tuple[int, Tuple[int, int]]:
    """Locate the principal (longest) stream in a basin.

    Parameters
    ----------
    stream : ndarray of bool or int
        Mask of the stream network.  Non-zero values denote stream cells.
    flowdir : ndarray of int
        D8 flow direction codes following the internal convention defined in
        :mod:`wmf_py.cu_py.basics`.
    dem : ndarray of float
        Elevation model used to extract the longitudinal profile.
    dx, dy : float
        Cell spacing in ``x`` and ``y`` directions respectively (metres).
    outlet_rc : tuple of int
        Row/column index of the basin outlet.

    Returns
    -------
    tuple
        ``(n_cells, start_rc)`` where ``n_cells`` is the number of cells in the
        identified principal stream and ``start_rc`` the upstream cell where the
        stream begins.

    Notes
    -----
    The function also stores internally the longitudinal profile of the
    principal stream.  It can be retrieved later using
    :func:`basin_ppalstream_cut`.
    """

    global _ppal_stream_profile

    stream_mask = stream.astype(bool)
    ny, nx = stream_mask.shape

    # ------------------------------------------------------------------
    # Build downstream pointers restricted to the stream network and
    # accumulate indegrees for topological ordering.
    receivers = np.full((ny, nx, 2), -1, dtype=int)
    indegree = np.zeros((ny, nx), dtype=int)

    for r in range(ny):
        for c in range(nx):
            if not stream_mask[r, c]:
                continue
            off = _D8.get(int(flowdir[r, c]))
            if not off:
                continue
            rr, cc = r + off[0], c + off[1]
            if 0 <= rr < ny and 0 <= cc < nx and stream_mask[rr, cc]:
                receivers[r, c] = (rr, cc)
                indegree[rr, cc] += 1

    # ------------------------------------------------------------------
    # Topological order of the stream network (upstream to downstream).
    order: List[Tuple[int, int]] = []
    q: deque[Tuple[int, int]] = deque()
    for r in range(ny):
        for c in range(nx):
            if stream_mask[r, c] and indegree[r, c] == 0:
                q.append((r, c))

    while q:
        r, c = q.popleft()
        order.append((r, c))
        rr, cc = receivers[r, c]
        if rr == -1:
            continue
        indegree[rr, cc] -= 1
        if indegree[rr, cc] == 0:
            q.append((rr, cc))

    # ------------------------------------------------------------------
    # Compute downstream distances for all stream cells using the obtained
    # order.  Distances are expressed in metres.
    dist = np.zeros((ny, nx), dtype=float)

    def _step_len(dr: int, dc: int) -> float:
        if dr == 0 or dc == 0:
            return dx if dc != 0 else dy
        return (dx**2 + dy**2) ** 0.5

    for r, c in reversed(order):
        rr, cc = receivers[r, c]
        if rr == -1:
            dist[r, c] = 0.0
            continue
        dr, dc = rr - r, cc - c
        dist[r, c] = dist[rr, cc] + _step_len(dr, dc)

    # Upstream cell corresponding to the maximum downstream distance.
    start_flat = int(np.argmax(dist))
    start_rc = divmod(start_flat, nx)

    # ------------------------------------------------------------------
    # Build longitudinal profile from ``start_rc`` to the outlet.
    prof_elev: List[float] = []
    prof_dist: List[float] = []
    prof_row: List[int] = []
    prof_col: List[int] = []

    r, c = start_rc
    cum_dist = 0.0
    while True:
        prof_elev.append(float(dem[r, c]))
        prof_dist.append(cum_dist)
        prof_row.append(r)
        prof_col.append(c)
        rr, cc = receivers[r, c]
        if rr == -1:
            break
        dr, dc = rr - r, cc - c
        cum_dist += _step_len(dr, dc)
        r, c = rr, cc

    profile = np.vstack(
        [
            np.array(prof_elev, dtype=float),
            np.array(prof_dist, dtype=float),
            np.array(prof_row, dtype=int),
            np.array(prof_col, dtype=int),
        ]
    )

    _ppal_stream_profile = profile

    return profile.shape[1], start_rc


def basin_ppalstream_cut(*args, **kwargs) -> np.ndarray:
    """Return the principal stream profile computed previously.

    Parameters are accepted for compatibility with the original Fortran
    routine but are ignored.  The function simply returns a copy of the
    profile stored by :func:`basin_ppalstream_find`.

    Returns
    -------
    ndarray
        Array with shape ``(4, N)`` containing elevation, cumulative distance
        and grid coordinates of the principal stream.
    """

    if _ppal_stream_profile is None:
        raise RuntimeError("basin_ppalstream_find must be called before cutting")
    return _ppal_stream_profile.copy()


def basin_ppal_hipsometric(
    dem_basin: np.ndarray, mask_basin: np.ndarray
) -> Dict[str, np.ndarray]:
    """Compute basic hypsometric information for a basin.

    The function wraps :func:`hypsometric_points` and :func:`hypsometric_curve`
    providing both the raw sorted elevations and a binned representation of the
    curve.

    Parameters
    ----------
    dem_basin : ndarray of float
        Elevation data of the basin.
    mask_basin : ndarray of bool
        Basin mask where ``True`` denotes active cells.

    Returns
    -------
    dict
        Dictionary with the keys ``elev_sorted`` and ``area_cum_frac`` holding
        the point-wise hypsometric data together with the output of
        :func:`hypsometric_curve` (``elev_bins``, ``area_frac``, ``min``,
        ``max`` and ``mean_elev``).
    """

    elev_sorted, area_cum_frac = hypsometric_points(dem_basin, mask_basin)
    curve = hypsometric_curve(dem_basin, mask_basin)
    out: Dict[str, np.ndarray] = {
        "elev_sorted": elev_sorted,
        "area_cum_frac": area_cum_frac,
    }
    out.update(curve)
    return out


def basin_time_to_out(
    flowdir: np.ndarray,
    dx: float,
    dy: float,
    slope: np.ndarray,
    n_manning: np.ndarray,
) -> np.ndarray:
    """Estimate travel time from each cell to the basin outlet.

    The computation follows D8 directions assuming a unit hydraulic radius in
    Manning's equation.  Cells involved in flow direction loops receive
    ``NaN`` as travel time.
    """

    ny, nx = flowdir.shape
    receivers = np.full((ny, nx, 2), -1, dtype=int)
    indegree = np.zeros((ny, nx), dtype=int)

    for r in range(ny):
        for c in range(nx):
            off = _D8.get(int(flowdir[r, c]))
            if not off:
                continue
            rr, cc = r + off[0], c + off[1]
            if 0 <= rr < ny and 0 <= cc < nx:
                receivers[r, c] = (rr, cc)
                indegree[rr, cc] += 1

    order: List[Tuple[int, int]] = []
    q: deque[Tuple[int, int]] = deque()
    for r in range(ny):
        for c in range(nx):
            if indegree[r, c] == 0:
                q.append((r, c))

    while q:
        r, c = q.popleft()
        order.append((r, c))
        rr, cc = receivers[r, c]
        if rr == -1:
            continue
        indegree[rr, cc] -= 1
        if indegree[rr, cc] == 0:
            q.append((rr, cc))

    def _cell_len(r: int, c: int, rr: int, cc: int) -> float:
        dr, dc = rr - r, cc - c
        if dr == 0 or dc == 0:
            return dx if dc != 0 else dy
        return (dx**2 + dy**2) ** 0.5

    vel = (1.0 / np.maximum(n_manning, 1e-6)) * np.sqrt(np.maximum(slope, 0.0))
    cell_time = np.full((ny, nx), np.nan, dtype=float)
    for r in range(ny):
        for c in range(nx):
            rr, cc = receivers[r, c]
            if rr == -1:
                continue
            length = _cell_len(r, c, rr, cc)
            v = vel[r, c]
            if v > 0:
                cell_time[r, c] = length / v

    time_to_out = np.full((ny, nx), np.nan, dtype=float)
    for r, c in reversed(order):
        rr, cc = receivers[r, c]
        if rr == -1:
            time_to_out[r, c] = 0.0
        elif not np.isnan(cell_time[r, c]) and not np.isnan(time_to_out[rr, cc]):
            time_to_out[r, c] = cell_time[r, c] + time_to_out[rr, cc]

    return time_to_out


def geo_hand(dem: np.ndarray, stream: np.ndarray) -> np.ndarray:
    """Compute Height Above Nearest Drainage (HAND).

    The algorithm performs a breadth-first search starting from all stream
    cells.  Each cell inherits the elevation of the closest stream cell in
    terms of grid steps (using the D8 neighbourhood).  The HAND value is the
    difference between the cell elevation and the inherited base elevation.

    Parameters
    ----------
    dem : ndarray of float
        Digital elevation model.
    stream : ndarray of bool or int
        Stream network mask where non-zero values denote stream cells.

    Returns
    -------
    ndarray of float
        Array with HAND values in metres.  Stream cells are assigned ``0``.
    """

    ny, nx = dem.shape
    stream_mask = stream.astype(bool)
    base = np.full_like(dem, np.nan, dtype=float)
    visited = np.zeros((ny, nx), dtype=bool)
    q: deque[Tuple[int, int]] = deque()

    for r in range(ny):
        for c in range(nx):
            if stream_mask[r, c]:
                base[r, c] = dem[r, c]
                visited[r, c] = True
                q.append((r, c))

    while q:
        r, c = q.popleft()
        for dr, dc in _D8.values():
            rr, cc = r + dr, c + dc
            if not (0 <= rr < ny and 0 <= cc < nx):
                continue
            if visited[rr, cc]:
                continue
            base[rr, cc] = base[r, c]
            visited[rr, cc] = True
            q.append((rr, cc))

    hand = dem.astype(float) - base
    hand[stream_mask] = 0.0
    return hand


# ---------------------------------------------------------------------------
# Sub-basin utilities
# ---------------------------------------------------------------------------

def basin_subbasin_find(
    stream: np.ndarray,
    nodes: Dict[int, Tuple[int, int]],
    edges: List[Dict[str, Any]],
    flowdir: np.ndarray,
    outlet_rc: Tuple[int, int],
) -> np.ndarray:
    """Assign sub-basin identifiers to all cells in the basin.

    The algorithm labels each cell with the identifier of the first downstream
    segment (edge) encountered along its D8 flow path.  Edge identifiers are
    taken from the ``edges`` list which is assumed to contain an ``id`` entry.
    Cells outside the basin or without a path to the outlet retain the value
    zero.
    """

    ny, nx = stream.shape
    submap = np.zeros((ny, nx), dtype=int)

    # Mark stream cells that belong to each edge (excluding the downstream
    # node which is shared with the next edge).
    for edge in edges:
        eid = int(edge.get("id"))
        for r, c in edge.get("cells", [])[:-1]:
            submap[r, c] = eid

    # Propagate the identifiers upstream for all remaining cells.
    for r in range(ny):
        for c in range(nx):
            if submap[r, c] != 0:
                continue
            path = []
            rr, cc = r, c
            eid = 0
            while True:
                path.append((rr, cc))
                if submap[rr, cc] > 0:
                    eid = submap[rr, cc]
                    break
                off = _D8.get(int(flowdir[rr, cc]))
                if not off:
                    break
                rr += off[0]
                cc += off[1]
                if not (0 <= rr < ny and 0 <= cc < nx):
                    eid = 0
                    break
            for pr, pc in path:
                submap[pr, pc] = eid

    return submap


def basin_subbasin_cut(subbasin_map: np.ndarray, sid: int) -> np.ndarray:
    """Return the mask corresponding to sub-basin ``sid``."""

    return (subbasin_map == sid).astype(np.uint8)


def basin_subbasin_map2subbasin(
    values_map: np.ndarray, subbasin_map: np.ndarray
) -> Dict[int, float]:
    """Aggregate map values per sub-basin using a simple sum."""

    out: Dict[int, float] = {}
    ids = np.unique(subbasin_map)
    for sid in ids:
        if sid == 0:
            continue
        out[int(sid)] = float(values_map[subbasin_map == sid].sum())
    return out


def basin_subbasin_nod(
    subbasin_map: np.ndarray, nodes: Dict[int, Tuple[int, int]], edges: List[Dict[str, Any]]
) -> Dict[int, int]:
    """Map sub-basin identifiers to their downstream node."""

    nodemap: Dict[int, int] = {}
    for edge in edges:
        nodemap[int(edge["id"])] = int(edge["node_v"])
    return nodemap


def basin_subbasin_horton(
    stream: np.ndarray, nodes: Dict[int, Tuple[int, int]], edges: List[Dict[str, Any]]
) -> Dict[str, Dict[int, int]]:
    """Compute Strahler (Horton) orders for the network."""

    edges_by_id = {int(e["id"]): e for e in edges}
    incoming: Dict[int, List[int]] = {nid: [] for nid in nodes}
    outgoing: Dict[int, List[int]] = {nid: [] for nid in nodes}
    for e in edges:
        eid = int(e["id"])
        u = int(e["node_u"])
        v = int(e["node_v"])
        outgoing[u].append(eid)
        incoming[v].append(eid)

    order_edge: Dict[int, int] = {}
    order_node: Dict[int, int] = {}
    upstream_orders: Dict[int, List[int]] = {nid: [] for nid in nodes}

    q: deque[int] = deque(n for n, inc in incoming.items() if len(inc) == 0)
    for n in q:
        order_node[n] = 1

    while q:
        n = q.popleft()
        for eid in outgoing.get(n, []):
            order_edge[eid] = order_node[n]
            v = edges_by_id[eid]["node_v"]
            upstream_orders[v].append(order_edge[eid])
            incoming[v].remove(eid)
            if not incoming[v]:
                vals = upstream_orders[v]
                maxo = max(vals)
                order_node[v] = maxo + 1 if vals.count(maxo) >= 2 else maxo
                q.append(v)

    return {"nodes": order_node, "edges": order_edge}

