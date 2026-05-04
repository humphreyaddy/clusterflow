"""End-to-end batch pipeline orchestration.

Wires Phase 1 (ingestion) → 3 (graph) → 4 (clustering) → 5 (analysis) →
6.1 (figures) into a single :func:`run_pipeline` call.
"""

from __future__ import annotations

import json
import logging
from datetime import date
from pathlib import Path

from clusterflow.analysis import (
    bootstrap_stability,
    build_temporal_dag,
    compute_centrality,
    index_case_candidates,
    save_bootstrap,
    save_centrality,
    save_dag,
)
from clusterflow.clustering import run_all
from clusterflow.config import ClusterFlowConfig
from clusterflow.graph import GraphConstructor, save_graph
from clusterflow.ingestion import IngestionPipeline
from clusterflow.models import (
    PipelineResult,
    TransmissionCluster,
)

log = logging.getLogger(__name__)


def _build_transmission_clusters(result: PipelineResult) -> list[TransmissionCluster]:
    out: list[TransmissionCluster] = []
    consensus = result.consensus
    by_cluster: dict[int, list[str]] = {}
    for iso, cid in consensus.assignments.items():
        by_cluster.setdefault(cid, []).append(iso)

    cent_lookup = {c.isolate_id: c for c in result.centrality}
    boot_lookup = {b.isolate_id: b for b in result.bootstrap}
    for cid, iso_ids in sorted(by_cluster.items()):
        sts = sorted(
            {result.isolates[i].sequence_type for i in iso_ids if result.isolates[i].sequence_type}
        )
        wards = sorted({result.isolates[i].ward for i in iso_ids})
        dates: list[date] = sorted(result.isolates[i].collection_date for i in iso_ids)
        risks = [
            (i, cent_lookup[i].transmission_risk_score)
            for i in iso_ids
            if i in cent_lookup
        ]
        idx_case = max(risks, key=lambda x: x[1])[0] if risks else None
        confidence = (
            sum(boot_lookup[i].stability_score for i in iso_ids if i in boot_lookup)
            / max(len(iso_ids), 1)
        )
        out.append(
            TransmissionCluster(
                cluster_id=int(cid),
                isolate_ids=sorted(iso_ids),
                sequence_types=sts,
                wards=wards,
                date_range=(dates[0], dates[-1]),
                index_case_candidate=idx_case,
                consensus_method=consensus.method,
                confidence=float(confidence),
            )
        )
    return out


def run_pipeline(config: ClusterFlowConfig, *, run_viz: bool = True) -> PipelineResult:
    """Execute the full ingestion → analysis pipeline. Writes all artefacts."""
    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "graph").mkdir(parents=True, exist_ok=True)
    (out_dir / "clusters").mkdir(parents=True, exist_ok=True)
    (out_dir / "analysis").mkdir(parents=True, exist_ok=True)

    log.info("[1/6] ingestion")
    ingest = IngestionPipeline().load(config)

    log.info("[2/6] graph construction")
    g = GraphConstructor(config.thresholds).build(ingest.snp, ingest.isolates)
    save_graph(g, out_dir / "graph")

    log.info("[3/6] clustering (parallel)")
    cluster_results = run_all(g, config.thresholds, config.clustering)
    consensus = cluster_results["consensus"]
    _save_cluster_assignments(cluster_results, out_dir / "clusters")

    log.info("[4/6] temporal DAG")
    dag = build_temporal_dag(g)
    save_dag(dag, out_dir / "analysis")

    log.info("[5/6] centrality")
    centrality = compute_centrality(dag, consensus)
    candidates = index_case_candidates(
        centrality,
        {iso: ingest.isolates[iso].collection_date for iso in ingest.isolates},
    )
    save_centrality(centrality, candidates, out_dir / "analysis")

    log.info("[6/6] bootstrap stability")
    boot = bootstrap_stability(ingest.snp, ingest.isolates, consensus, config)
    save_bootstrap(boot, out_dir / "analysis")

    result = PipelineResult(
        project_name=config.project_name,
        n_isolates=len(ingest.isolates),
        isolates=ingest.isolates,
        snp_matrix=ingest.snp,
        cluster_assignments={k: v for k, v in cluster_results.items() if k != "consensus"},
        consensus=consensus,
        centrality=centrality,
        bootstrap=boot,
    )
    result.transmission_clusters = _build_transmission_clusters(result)
    _save_pipeline_summary(result, out_dir)

    if run_viz and config.visualization.static:
        from clusterflow.viz import generate_all_figures

        log.info("rendering static figures")
        generate_all_figures(result, g, dag, out_dir / "figures")

    log.info("pipeline complete: %d clusters detected", consensus.n_clusters)
    return result


def _save_cluster_assignments(
    results: dict, output_dir: Path
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    import pandas as pd

    rows: list[dict] = []
    methods = list(results.keys())
    isolate_ids: set[str] = set()
    for a in results.values():
        isolate_ids.update(a.assignments)
    for iso in sorted(isolate_ids):
        row: dict = {"isolate_id": iso}
        for m in methods:
            row[m] = results[m].assignments.get(iso)
        consensus_obj = results.get("consensus")
        if consensus_obj is not None:
            row["agreement_score"] = (
                consensus_obj.agreement_score or {}
            ).get(iso)
            row["ambiguous"] = (consensus_obj.ambiguous or {}).get(iso)
        rows.append(row)
    pd.DataFrame(rows).to_csv(output_dir / "cluster_assignments.csv", index=False)


def _save_pipeline_summary(result: PipelineResult, output_dir: Path) -> None:
    summary = {
        "project_name": result.project_name,
        "n_isolates": result.n_isolates,
        "n_clusters_consensus": result.consensus.n_clusters,
        "n_clusters_per_method": {
            m: a.n_clusters for m, a in result.cluster_assignments.items()
        },
        "transmission_clusters": [
            {
                "cluster_id": c.cluster_id,
                "n_isolates": len(c.isolate_ids),
                "sequence_types": c.sequence_types,
                "wards": c.wards,
                "date_range": [c.date_range[0].isoformat(), c.date_range[1].isoformat()],
                "index_case_candidate": c.index_case_candidate,
                "confidence": c.confidence,
            }
            for c in result.transmission_clusters
        ],
    }
    (output_dir / "pipeline_summary.json").write_text(
        json.dumps(summary, indent=2)
    )
