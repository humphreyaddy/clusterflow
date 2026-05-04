"""Phase 1.3 — shared Pydantic data models.

These are the canonical types passed between every stage of the pipeline.
See transmission_cluster_tool_plan.md §1.3.
"""

from __future__ import annotations

from datetime import date
from typing import Any

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, field_validator


class Isolate(BaseModel):
    """A single sequenced isolate with its epi metadata + (optional) MLST."""

    isolate_id: str
    collection_date: date
    facility: str
    ward: str
    sequence_type: str | None = None
    mlst_alleles: dict[str, int] | None = None

    model_config = ConfigDict(frozen=True)


class SNPMatrix(BaseModel):
    """Pairwise SNP distance matrix.

    `distances` is a square symmetric numpy array of shape (n, n) where
    `distances[i, j]` is the SNP count between `isolate_ids[i]` and
    `isolate_ids[j]`. Diagonal is zero.
    """

    isolate_ids: list[str]
    distances: np.ndarray

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @field_validator("distances")
    @classmethod
    def _square_symmetric(cls, v: np.ndarray) -> np.ndarray:
        if v.ndim != 2 or v.shape[0] != v.shape[1]:
            raise ValueError(f"distance matrix must be square, got shape {v.shape}")
        if not np.allclose(v, v.T, atol=1e-6):
            raise ValueError("distance matrix must be symmetric")
        if np.any(v < 0):
            raise ValueError("distance matrix must be non-negative")
        if not np.allclose(np.diag(v), 0, atol=1e-6):
            raise ValueError("distance matrix diagonal must be zero")
        return v

    @property
    def n(self) -> int:
        return len(self.isolate_ids)

    def index(self, isolate_id: str) -> int:
        return self.isolate_ids.index(isolate_id)

    def distance(self, a: str, b: str) -> float:
        return float(self.distances[self.index(a), self.index(b)])


class WeightedEdge(BaseModel):
    source: str
    target: str
    snp_distance: int
    day_delta: int
    mlst_mismatch: bool
    composite_weight: float


class ClusterAssignment(BaseModel):
    """Output of one clustering algorithm (or the consensus)."""

    method: str
    assignments: dict[str, int]  # isolate_id → cluster_id
    n_clusters: int
    stability_scores: dict[str, float] | None = None
    ambiguous: dict[str, bool] | None = None
    agreement_score: dict[str, float] | None = None

    @field_validator("n_clusters")
    @classmethod
    def _check_count(cls, v: int, info: Any) -> int:
        assignments: dict[str, int] | None = info.data.get("assignments")
        if assignments is not None:
            unique = len(set(assignments.values()))
            if unique != v:
                raise ValueError(
                    f"n_clusters={v} disagrees with unique assignments={unique}"
                )
        return v


class TransmissionCluster(BaseModel):
    cluster_id: int
    isolate_ids: list[str]
    sequence_types: list[str]
    wards: list[str]
    date_range: tuple[date, date]
    index_case_candidate: str | None = None
    consensus_method: str = "consensus"
    confidence: float = 0.0


class CentralityScores(BaseModel):
    """Per-isolate centrality metrics from Phase 5.2."""

    isolate_id: str
    cluster_id: int
    betweenness: float
    in_degree: int
    out_degree: int
    closeness: float
    transmission_risk_score: float


class BootstrapResult(BaseModel):
    """Per-isolate bootstrap stability — Phase 5.3."""

    isolate_id: str
    cluster_id: int
    stability_score: float
    classification: str  # "stable" | "uncertain" | "unstable"


class PipelineResult(BaseModel):
    """Top-level pipeline output."""

    project_name: str
    n_isolates: int
    isolates: dict[str, Isolate]
    snp_matrix: SNPMatrix
    cluster_assignments: dict[str, ClusterAssignment]  # method → assignment
    consensus: ClusterAssignment
    transmission_clusters: list[TransmissionCluster] = Field(default_factory=list)
    centrality: list[CentralityScores] = Field(default_factory=list)
    bootstrap: list[BootstrapResult] = Field(default_factory=list)
    metadata: dict[str, Any] = Field(default_factory=dict)

    model_config = ConfigDict(arbitrary_types_allowed=True)
