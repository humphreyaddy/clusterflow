"""Phase 4.4 — parallel wrapper around the three clustering algorithms."""

from __future__ import annotations

import logging

import igraph as ig
from joblib import Parallel, delayed

from clusterflow.clustering.consensus import consensus_assignment
from clusterflow.clustering.leiden import leiden_clusters
from clusterflow.clustering.snp_chains import snp_chain_clusters
from clusterflow.clustering.spectral import spectral_clusters
from clusterflow.config import ClusteringConfig, ThresholdsConfig
from clusterflow.models import ClusterAssignment

log = logging.getLogger(__name__)


def _dispatch(
    method: str,
    g: ig.Graph,
    thresholds: ThresholdsConfig,
    cfg: ClusteringConfig,
) -> ClusterAssignment:
    if method == "snp_chains":
        return snp_chain_clusters(g.copy(), thresholds.snp_cutoff, thresholds.day_cutoff)
    if method == "leiden":
        return leiden_clusters(
            g.copy(),
            resolution=cfg.leiden_resolution,
            random_state=cfg.random_state,
            n_jobs=1,  # avoid nested parallelism
        )
    if method == "spectral":
        return spectral_clusters(g.copy(), k="auto", random_state=cfg.random_state)
    raise ValueError(f"unknown clustering method: {method}")


def run_all(
    g: ig.Graph,
    thresholds: ThresholdsConfig,
    cfg: ClusteringConfig,
) -> dict[str, ClusterAssignment]:
    """Run every method in parallel + compute consensus.

    Returns a dict: ``method → ClusterAssignment`` plus ``"consensus"``.
    """
    methods = list(cfg.methods)
    log.info("running clustering methods: %s", methods)
    results = Parallel(n_jobs=min(len(methods), 3), backend="threading")(
        delayed(_dispatch)(m, g, thresholds, cfg) for m in methods
    )
    out: dict[str, ClusterAssignment] = {a.method: a for a in results}
    out["consensus"] = consensus_assignment(list(results))
    log.info(
        "clustering complete: %s",
        {k: v.n_clusters for k, v in out.items()},
    )
    return out
