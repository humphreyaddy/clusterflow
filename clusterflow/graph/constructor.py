"""Phase 3.1 — weighted graph constructor.

Builds an undirected ``igraph.Graph`` whose vertices carry every ``Isolate``
field and whose edges carry SNP distance, day delta, MLST mismatch flag, and
the composite weight ``alpha*snp + beta*time + gamma*mlst``.
"""

from __future__ import annotations

from datetime import date
from typing import TYPE_CHECKING

import igraph as ig
import numpy as np

from clusterflow.config import ThresholdsConfig
from clusterflow.models import Isolate, SNPMatrix

if TYPE_CHECKING:  # avoid pydantic import in hot path
    pass


def _composite_weight(
    snp: float, days: int, mlst_mismatch: bool, t: ThresholdsConfig
) -> float:
    snp_term = snp / max(t.snp_cutoff, 1)
    time_term = days / max(t.day_cutoff, 1)
    mlst_term = 1.0 if mlst_mismatch else 0.0
    return (
        t.edge_weight_alpha * snp_term
        + t.edge_weight_beta * time_term
        + t.edge_weight_gamma * mlst_term
    )


class GraphConstructor:
    """Build the weighted transmission graph from SNP + epi metadata."""

    def __init__(self, thresholds: ThresholdsConfig) -> None:
        self.t = thresholds

    def build(
        self, snp: SNPMatrix, isolates: dict[str, Isolate]
    ) -> ig.Graph:
        ids = snp.isolate_ids
        n = len(ids)
        g = ig.Graph(n=n, directed=False)
        g.vs["name"] = ids
        g.vs["isolate_id"] = ids
        g.vs["collection_date"] = [
            isolates[i].collection_date.isoformat() for i in ids
        ]
        g.vs["facility"] = [isolates[i].facility for i in ids]
        g.vs["ward"] = [isolates[i].ward for i in ids]
        g.vs["sequence_type"] = [
            isolates[i].sequence_type or "" for i in ids
        ]

        edges: list[tuple[int, int]] = []
        attrs: dict[str, list] = {
            "snp_distance": [],
            "day_delta": [],
            "mlst_mismatch": [],
            "composite_weight": [],
        }
        D = snp.distances
        dates = [isolates[i].collection_date for i in ids]
        sts = [isolates[i].sequence_type for i in ids]

        for i in range(n):
            for j in range(i + 1, n):
                d_snp = float(D[i, j])
                if d_snp > self.t.snp_cutoff:
                    continue
                d_days = abs((dates[i] - dates[j]).days)
                if d_days > self.t.day_cutoff:
                    continue
                mismatch = self._mlst_mismatch(sts[i], sts[j])
                w = _composite_weight(d_snp, d_days, mismatch, self.t)
                edges.append((i, j))
                attrs["snp_distance"].append(int(round(d_snp)))
                attrs["day_delta"].append(int(d_days))
                attrs["mlst_mismatch"].append(bool(mismatch))
                attrs["composite_weight"].append(float(w))

        g.add_edges(edges)
        for k, vs in attrs.items():
            g.es[k] = vs
        return g

    @staticmethod
    def _mlst_mismatch(a: str | None, b: str | None) -> bool:
        # Both unknown → conservatively treat as match (no penalty)
        if not a or not b:
            return False
        return a != b

    def add_isolate(
        self,
        g: ig.Graph,
        new_isolate: Isolate,
        snp_distances_to_existing: dict[str, float],
    ) -> list[dict]:
        """Phase 3.3 — incremental add. Returns list of new-edge dicts."""
        existing_ids: list[str] = list(g.vs["isolate_id"])
        if new_isolate.isolate_id in existing_ids:
            raise ValueError(f"isolate {new_isolate.isolate_id} already in graph")
        new_idx = g.vcount()
        g.add_vertex(
            name=new_isolate.isolate_id,
            isolate_id=new_isolate.isolate_id,
            collection_date=new_isolate.collection_date.isoformat(),
            facility=new_isolate.facility,
            ward=new_isolate.ward,
            sequence_type=new_isolate.sequence_type or "",
        )

        new_edges: list[tuple[int, int]] = []
        new_attrs: dict[str, list] = {
            "snp_distance": [],
            "day_delta": [],
            "mlst_mismatch": [],
            "composite_weight": [],
        }
        report: list[dict] = []
        for j, other_id in enumerate(existing_ids):
            d_snp = snp_distances_to_existing.get(other_id)
            if d_snp is None or d_snp > self.t.snp_cutoff:
                continue
            other_date = date.fromisoformat(g.vs[j]["collection_date"])
            d_days = abs((new_isolate.collection_date - other_date).days)
            if d_days > self.t.day_cutoff:
                continue
            other_st = g.vs[j]["sequence_type"] or None
            mismatch = self._mlst_mismatch(new_isolate.sequence_type, other_st)
            w = _composite_weight(d_snp, d_days, mismatch, self.t)
            new_edges.append((new_idx, j))
            new_attrs["snp_distance"].append(int(round(d_snp)))
            new_attrs["day_delta"].append(int(d_days))
            new_attrs["mlst_mismatch"].append(bool(mismatch))
            new_attrs["composite_weight"].append(float(w))
            report.append(
                {
                    "source": new_isolate.isolate_id,
                    "target": other_id,
                    "snp_distance": int(round(d_snp)),
                    "day_delta": int(d_days),
                    "mlst_mismatch": bool(mismatch),
                    "composite_weight": float(w),
                }
            )

        if new_edges:
            start = g.ecount()
            g.add_edges(new_edges)
            for k, vs in new_attrs.items():
                # Pad existing edges to maintain attribute length.
                cur = g.es[k] if k in g.es.attributes() else [None] * start
                g.es[k] = list(cur)[:start] + vs
        return report


def graph_summary(g: ig.Graph) -> dict:
    n_v = g.vcount()
    n_e = g.ecount()
    if n_v < 2:
        density = 0.0
    else:
        density = n_e / (n_v * (n_v - 1) / 2)
    components = g.connected_components()
    weights = g.es["composite_weight"] if n_e else []
    snp = g.es["snp_distance"] if n_e else []
    return {
        "n_vertices": n_v,
        "n_edges": n_e,
        "density": float(density),
        "n_components": len(components),
        "largest_component_size": max(map(len, components)) if components else 0,
        "mean_degree": float(np.mean(g.degree())) if n_v else 0.0,
        "weight_min": float(min(weights)) if weights else 0.0,
        "weight_max": float(max(weights)) if weights else 0.0,
        "weight_mean": float(np.mean(weights)) if weights else 0.0,
        "snp_min": int(min(snp)) if snp else 0,
        "snp_max": int(max(snp)) if snp else 0,
    }
