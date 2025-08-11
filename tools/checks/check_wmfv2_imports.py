#!/usr/bin/env python3
"""Analyze wmfv2 calls to ensure they resolve to wmf_py.* modules."""
from __future__ import annotations

import ast
import json
from pathlib import Path

REQUIRED_CALLS = [
    "shia_v1",
    "dir_reclass_opentopo",
    "dir_reclass_rwatershed",
    "basin_acum",
    "basin_find",
    "basin_cut",
    "basin_2map",
    "basin_map2basin",
    "stream_find",
    "stream_cut",
    "stream_find_to_corr",
    "basin_netxy_find",
    "basin_netxy_cut",
    "basin_subbasin_find",
    "basin_subbasin_cut",
    "basin_subbasin_map2subbasin",
    "basin_subbasin_nod",
    "basin_subbasin_horton",
    "hypsometric_curve",
    "hypsometric_points",
    "time_to_outlet",
    "hand",
    "stream_seed_from_coords",
    "stream_threshold_nearby",
    "stream_find_nearby",
    "hydro_distance_and_receiver",
]

ROOT = Path(__file__).resolve().parents[2]
WMFV2_FILE = ROOT / "wmf_heavy" / "wmfv2.py"
REPORT_FILE = ROOT / "reports" / "wmfv2_calls_report.json"


def build_import_map(tree: ast.AST) -> dict[str, str]:
    import_map: dict[str, str] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            module = node.module or ""
            for alias in node.names:
                name = alias.asname or alias.name
                import_map[name] = module
        elif isinstance(node, ast.Import):
            for alias in node.names:
                name = alias.asname or alias.name
                import_map[name] = alias.name
    return import_map


def collect_calls(tree: ast.AST, import_map: dict[str, str]) -> dict[str, str | None]:
    calls: dict[str, str | None] = {}

    class Visitor(ast.NodeVisitor):
        def visit_Call(self, node: ast.Call) -> None:
            func = node.func
            name = None
            module = None
            if isinstance(func, ast.Name):
                name = func.id
                module = import_map.get(name)
            elif isinstance(func, ast.Attribute) and isinstance(func.value, ast.Name):
                name = func.attr
                module = import_map.get(func.value.id)
            if name:
                calls[name] = module
            self.generic_visit(node)

    Visitor().visit(tree)
    return calls


def main() -> None:
    tree = ast.parse(WMFV2_FILE.read_text(encoding="utf-8"))
    import_map = build_import_map(tree)
    calls = collect_calls(tree, import_map)

    called_functions = {name: module for name, module in calls.items()}
    missing = [name for name in REQUIRED_CALLS if name not in called_functions]
    unmapped = [
        name
        for name, module in called_functions.items()
        if name in REQUIRED_CALLS and (module is None or not module.startswith("wmf_py"))
    ]

    REPORT_FILE.parent.mkdir(parents=True, exist_ok=True)
    REPORT_FILE.write_text(
        json.dumps(
            {
                "called_functions": called_functions,
                "missing_required_calls": missing,
                "unmapped_required_calls": unmapped,
            },
            indent=2,
            sort_keys=True,
        )
    )

    print(REPORT_FILE.read_text())


if __name__ == "__main__":
    main()
