"""Phase 7.1 — incremental graph + reclustering for streaming mode."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import igraph as ig

from clusterflow.clustering.leiden import leiden_clusters
from clusterflow.config import ClusterFlowConfig
from clusterflow.graph.constructor import GraphConstructor
from clusterflow.models import ClusterAssignment, Isolate

log = logging.getLogger(__name__)


@dataclass
class IncrementalState:
    config: ClusterFlowConfig
    graph: ig.Graph
    isolates: dict[str, Isolate]
    consensus: ClusterAssignment
    alert_threshold: float = 0.7
    new_edges: list[dict] = field(default_factory=list)


@dataclass
class IncrementalUpdate:
    cluster_assigned: int
    new_cluster_formed: bool
    transmission_alert: bool
    centrality_score: float
    new_edges: list[dict]


def add_isolate(
    state: IncrementalState,
    new_isolate: Isolate,
    snp_distances_to_existing: dict[str, float],
) -> IncrementalUpdate:
    constructor = GraphConstructor(state.config.thresholds)
    new_edges = constructor.add_isolate(
        state.graph, new_isolate, snp_distances_to_existing
    )
    state.new_edges.extend(new_edges)
    state.isolates[new_isolate.isolate_id] = new_isolate

    # Re-cluster (Leiden only — fast)
    new_assign = leiden_clusters(
        state.graph,
        resolution=state.config.clustering.leiden_resolution,
        random_state=state.config.clustering.random_state,
        n_jobs=1,
    )
    cluster_id = new_assign.assignments[new_isolate.isolate_id]
    pre_cluster_ids = set(state.consensus.assignments.values())
    new_cluster_formed = cluster_id not in pre_cluster_ids

    # Quick centrality on the affected component only
    component_ids = state.graph.subcomponent(state.graph.vcount() - 1)
    sub = state.graph.subgraph(component_ids)
    if sub.ecount():
        weights = [max(w, 1e-6) for w in sub.es["composite_weight"]]
        between = sub.betweenness(weights=weights)
        idx = list(sub.vs["isolate_id"]).index(new_isolate.isolate_id)
        max_b = max(between) or 1.0
        risk_score = between[idx] / max_b
    else:
        risk_score = 0.0

    state.consensus = new_assign
    alert = risk_score >= state.alert_threshold
    if alert:
        log.warning(
            "TRANSMISSION_ALERT  isolate=%s  cluster=%s  risk=%.2f",
            new_isolate.isolate_id,
            cluster_id,
            risk_score,
        )
    return IncrementalUpdate(
        cluster_assigned=int(cluster_id),
        new_cluster_formed=new_cluster_formed,
        transmission_alert=alert,
        centrality_score=float(risk_score),
        new_edges=new_edges,
    )
