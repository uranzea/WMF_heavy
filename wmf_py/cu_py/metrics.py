"""Basin metrics utilities."""
from __future__ import annotations

from collections import deque
from typing import Any, Dict, List, Tuple

from .basics import _D8

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

