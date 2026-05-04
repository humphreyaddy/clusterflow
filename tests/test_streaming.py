"""Phase 7 integration test — replay isolates one-by-one."""

from __future__ import annotations

from datetime import timedelta

from clusterflow.clustering.consensus import consensus_assignment
from clusterflow.clustering.leiden import leiden_clusters
from clusterflow.graph.constructor import GraphConstructor
from clusterflow.streaming.incremental import IncrementalState, add_isolate


def test_incremental_replay_matches_n_isolates(kp_config, kp_ingestion):
    # Seed graph with the first 25 isolates; stream the remaining 3.
    iso_ids = list(kp_ingestion.isolates)
    seed_ids = iso_ids[:25]
    stream_ids = iso_ids[25:]

    seed_isolates = {i: kp_ingestion.isolates[i] for i in seed_ids}
    snp_seed = kp_ingestion.snp.distances
    snp_idx = {i: kp_ingestion.snp.isolate_ids.index(i) for i in iso_ids}

    # Build seed-graph subgraph manually
    from clusterflow.models import SNPMatrix
    import numpy as np

    sub_idx = [kp_ingestion.snp.isolate_ids.index(i) for i in seed_ids]
    seed_snp = SNPMatrix(
        isolate_ids=seed_ids,
        distances=snp_seed[np.ix_(sub_idx, sub_idx)],
    )
    g = GraphConstructor(kp_config.thresholds).build(seed_snp, seed_isolates)
    leiden = leiden_clusters(g, resolution=1.0, random_state=42, n_jobs=1)
    state = IncrementalState(
        config=kp_config,
        graph=g,
        isolates=dict(seed_isolates),
        consensus=consensus_assignment([leiden]),
    )

    for new_id in stream_ids:
        new_iso = kp_ingestion.isolates[new_id]
        # Only feed distances under the threshold to avoid huge dicts; the
        # constructor will re-filter by threshold anyway.
        distances = {
            other: float(snp_seed[snp_idx[new_id], snp_idx[other]])
            for other in state.graph.vs["isolate_id"]
        }
        update = add_isolate(state, new_iso, distances)
        assert update.cluster_assigned is not None

    assert state.graph.vcount() == 28
