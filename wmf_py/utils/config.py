"""Configuration helpers for the WMF pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional

try:  # PyYAML is optional; keep tests lightweight
    import yaml
except Exception:  # pragma: no cover - fallback when yaml missing
    yaml = None  # type: ignore


@dataclass
class SeedCoords:
    x: float
    y: float
    search_radius_cells: int = 4


@dataclass
class StreamsConfig:
    candidate_thresholds: Optional[List[int]] = None


@dataclass
class OutputsConfig:
    hdnd: bool = True
    aquien: bool = True


@dataclass
class Config:
    seed_coords: Optional[SeedCoords] = None
    streams: StreamsConfig = StreamsConfig()
    outputs: OutputsConfig = OutputsConfig()


def load_config(path: str) -> Config:
    """Load configuration from a YAML file."""

    if yaml is None:
        raise RuntimeError("PyYAML is required to load configuration")

    with open(path, "r", encoding="utf8") as f:
        data: Dict[str, Any] = yaml.safe_load(f) or {}

    seed = data.get("seed_coords")
    seed_cfg = None
    if seed is not None:
        seed_cfg = SeedCoords(
            x=float(seed.get("x", 0.0)),
            y=float(seed.get("y", 0.0)),
            search_radius_cells=int(seed.get("search_radius_cells", 4)),
        )

    streams = data.get("streams", {})
    streams_cfg = StreamsConfig(
        candidate_thresholds=streams.get("candidate_thresholds")
    )

    outputs = data.get("outputs", {})
    outputs_cfg = OutputsConfig(
        hdnd=bool(outputs.get("hdnd", False)),
        aquien=bool(outputs.get("aquien", False)),
    )

    return Config(seed_coords=seed_cfg, streams=streams_cfg, outputs=outputs_cfg)


__all__ = [
    "SeedCoords",
    "StreamsConfig",
    "OutputsConfig",
    "Config",
    "load_config",
]

