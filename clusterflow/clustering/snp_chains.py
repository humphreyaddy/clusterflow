"""Phase 4.1 — SNP chain clustering (SNP2Cluster-compatible)."""

from __future__ import annotations

from datetime import date

import igraph as ig

from clusterflow.models import ClusterAssignment


def snp_chain_clusters(
    g: ig.Graph,
    snp_cutoff: int,
    day_cutoff: int,
) -> ClusterAssignment:
    """Re-implementation of SNP2Cluster's chain logic.

    Walk isolates in collection-date order. Assign each isolate to an existing
    cluster if it is within ``snp_cutoff`` SNPs of any current member AND the
    date gap to that member is within ``day_cutoff`` days; otherwise start a
    new cluster. First qualifying assignment wins (stable over date order).
    """
    n = g.vcount()
    ids: list[str] = list(g.vs["isolate_id"])
    dates: list[date] = [date.fromisoformat(d) for d in g.vs["collection_date"]]

    # Pre-build SNP-distance lookup from edge attributes
    snp_lookup: dict[tuple[int, int], int] = {}
    for e in g.es:
        s, t = e.tuple
        d = int(e["snp_distance"])
        snp_lookup[(s, t)] = d
        snp_lookup[(t, s)] = d

    order = sorted(range(n), key=lambda i: dates[i])
    cluster_of: dict[int, int] = {}
    clusters: list[list[int]] = []

    for i in order:
        assigned = False
        for cid, members in enumerate(clusters):
            for m in members:
                snp_d = snp_lookup.get((i, m))
                if snp_d is None or snp_d > snp_cutoff:
                    continue
                day_d = abs((dates[i] - dates[m]).days)
                if day_d > day_cutoff:
                    continue
                clusters[cid].append(i)
                cluster_of[i] = cid
                assigned = True
                break
            if assigned:
                break
        if not assigned:
            cluster_of[i] = len(clusters)
            clusters.append([i])

    assignments = {ids[i]: cluster_of[i] for i in range(n)}
    return ClusterAssignment(
        method="snp_chains",
        assignments=assignments,
        n_clusters=len(clusters),
    )
