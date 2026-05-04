"""Phase 1.2 — YAML config loader and Pydantic validator.

The YAML config is the single source of truth for thresholds, paths, and
backend selection. See transmission_cluster_tool_plan.md §1.2.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


class ConfigError(Exception):
    """Raised when a config file is missing, malformed, or fails validation."""


class InputsConfig(BaseModel):
    snp_matrix: Path
    mlst_profiles: Path | None = None
    epi_metadata: Path

    @field_validator("snp_matrix", "mlst_profiles", "epi_metadata")
    @classmethod
    def _expand(cls, v: Path | None) -> Path | None:
        if v is None:
            return None
        return Path(v).expanduser()


class ThresholdsConfig(BaseModel):
    snp_cutoff: int = Field(20, ge=0)
    day_cutoff: int = Field(14, ge=0)
    edge_weight_alpha: float = Field(0.4, ge=0.0, le=1.0)
    edge_weight_beta: float = Field(0.4, ge=0.0, le=1.0)
    edge_weight_gamma: float = Field(0.2, ge=0.0, le=1.0)

    @model_validator(mode="after")
    def _weights_sum(self) -> "ThresholdsConfig":
        s = self.edge_weight_alpha + self.edge_weight_beta + self.edge_weight_gamma
        if not 0.99 <= s <= 1.01:
            raise ValueError(
                f"edge_weight_alpha+beta+gamma must sum to 1.0 (got {s:.3f})"
            )
        return self


class ClusteringConfig(BaseModel):
    methods: list[Literal["snp_chains", "leiden", "spectral"]] = Field(
        default_factory=lambda: ["snp_chains", "leiden", "spectral"]
    )
    bootstrap_n: int = Field(500, ge=0)
    leiden_resolution: float | Literal["auto"] = 1.0
    random_state: int = 42


class DistanceEngineConfig(BaseModel):
    backend: Literal["auto", "cpu", "gpu", "pairsnp"] = "auto"
    n_jobs: int = -1


class StreamingConfig(BaseModel):
    enabled: bool = False
    api_port: int = Field(8000, ge=1024, le=65535)


class VisualizationConfig(BaseModel):
    static: bool = True
    dashboard: bool = False
    dashboard_port: int = Field(8050, ge=1024, le=65535)


class ClusterFlowConfig(BaseModel):
    project_name: str
    output_dir: Path
    inputs: InputsConfig
    thresholds: ThresholdsConfig = Field(default_factory=ThresholdsConfig)
    clustering: ClusteringConfig = Field(default_factory=ClusteringConfig)
    distance_engine: DistanceEngineConfig = Field(default_factory=DistanceEngineConfig)
    streaming: StreamingConfig = Field(default_factory=StreamingConfig)
    visualization: VisualizationConfig = Field(default_factory=VisualizationConfig)

    model_config = ConfigDict(extra="forbid")

    @field_validator("output_dir")
    @classmethod
    def _expand(cls, v: Path) -> Path:
        return Path(v).expanduser()

    @classmethod
    def from_yaml(cls, path: str | Path) -> "ClusterFlowConfig":
        p = Path(path).expanduser()
        if not p.exists():
            raise ConfigError(f"config file not found: {p}")
        try:
            with p.open() as fh:
                raw = yaml.safe_load(fh)
        except yaml.YAMLError as e:
            raise ConfigError(f"YAML parse error in {p}: {e}") from e
        if raw is None:
            raise ConfigError(f"config file is empty: {p}")
        cfg = cls(**raw)
        cfg.validate_paths()
        return cfg

    def validate_paths(self) -> None:
        """Confirm input files exist (called automatically by from_yaml)."""
        if not self.inputs.snp_matrix.exists():
            raise ConfigError(f"snp_matrix not found: {self.inputs.snp_matrix}")
        if not self.inputs.epi_metadata.exists():
            raise ConfigError(f"epi_metadata not found: {self.inputs.epi_metadata}")
        if (
            self.inputs.mlst_profiles is not None
            and not self.inputs.mlst_profiles.exists()
        ):
            raise ConfigError(f"mlst_profiles not found: {self.inputs.mlst_profiles}")

    def to_yaml(self, path: str | Path) -> None:
        Path(path).write_text(yaml.safe_dump(self.model_dump(mode="json"), sort_keys=False))


def diff_configs(a: ClusterFlowConfig, b: ClusterFlowConfig) -> dict[str, tuple]:
    """Return a flat dict of (key → (value_a, value_b)) for keys that differ."""
    flat_a = _flatten(a.model_dump(mode="json"))
    flat_b = _flatten(b.model_dump(mode="json"))
    keys = set(flat_a) | set(flat_b)
    return {
        k: (flat_a.get(k), flat_b.get(k))
        for k in sorted(keys)
        if flat_a.get(k) != flat_b.get(k)
    }


def _flatten(d: dict, prefix: str = "") -> dict:
    out: dict = {}
    for k, v in d.items():
        key = f"{prefix}.{k}" if prefix else k
        if isinstance(v, dict):
            out.update(_flatten(v, key))
        else:
            out[key] = v
    return out
