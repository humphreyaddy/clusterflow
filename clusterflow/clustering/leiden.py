"""Phase 4.2 — Leiden community detection."""

from __future__ import annotations

import logging

import igraph as ig
import leidenalg
import numpy as np
from joblib import Parallel, delayed

from clusterflow.models import ClusterAssignment

log = logging.getLogger(__name__)


def _similarity_weights(g: ig.Graph) -> list[float]:
    """Convert composite-distance edge weights → similarity (lower distance →
    larger similarity). Leiden modularity expects similarity, not distance."""
    if g.ecount() == 0:
        return []
    w = np.array(g.es["composite_weight"], dtype=float)
    # Map distance ∈ [0, ~3] to similarity ∈ (0, 1]
    return list(np.exp(-w))


def _run_one(g: ig.Graph, resolution: float, sim: list[float], seed: int) -> tuple[float, list[int]]:
    partition = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        weights=sim if sim else None,
        resolution_parameter=resolution,
        n_iterations=10,
        seed=seed,
    )
    return float(partition.modularity), list(partition.membership)


def leiden_clusters(
    g: ig.Graph,
    resolution: float | str = 1.0,
    random_state: int = 42,
    n_jobs: int = -1,
) -> ClusterAssignment:
    """Run Leiden, optionally sweeping resolution to maximise modularity Q."""
    ids = list(g.vs["isolate_id"])
    sim = _similarity_weights(g)

    if g.ecount() == 0:
        # Each node is its own cluster.
        assignments = {iso: i for i, iso in enumerate(ids)}
        return ClusterAssignment(
            method="leiden", assignments=assignments, n_clusters=len(ids)
        )

    if resolution == "auto":
        resolutions = np.arange(0.5, 2.01, 0.1).round(2).tolist()
        results = Parallel(n_jobs=n_jobs, backend="loky")(
            delayed(_run_one)(g, r, sim, random_state) for r in resolutions
        )
        modularities = [r[0] for r in results]
        memberships = [r[1] for r in results]
        best = int(np.argmax(modularities))
        log.info(
            "leiden resolution sweep: best=%.2f modularity=%.4f",
            resolutions[best],
            modularities[best],
        )
        membership = memberships[best]
    else:
        _, membership = _run_one(g, float(resolution), sim, random_state)

    assignments = {ids[i]: int(c) for i, c in enumerate(membership)}
    return ClusterAssignment(
        method="leiden",
        assignments=assignments,
        n_clusters=len(set(membership)),
    )
