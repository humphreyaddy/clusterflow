"""Tests for clusterflow.analysis."""

from __future__ import annotations

from clusterflow.analysis import (
    bootstrap_stability,
    build_temporal_dag,
    cluster_stability_summary,
    compute_centrality,
    index_case_candidates,
)
from clusterflow.clustering import run_all
from clusterflow.graph import GraphConstructor


def _setup(kp_config, kp_ingestion):
    g = GraphConstructor(kp_config.thresholds).build(
        kp_ingestion.snp, kp_ingestion.isolates
    )
    results = run_all(g, kp_config.thresholds, kp_config.clustering)
    return g, results


def test_dag_is_acyclic(kp_config, kp_ingestion):
    g, _ = _setup(kp_config, kp_ingestion)
    dag = build_temporal_dag(g, uncertainty_days=2)
    # Even allowing bidirectional within ±2 days, our cycle-breaker should leave a DAG
    assert dag.is_dag()


def test_centrality_index_case(kp_config, kp_ingestion):
    g, results = _setup(kp_config, kp_ingestion)
    dag = build_temporal_dag(g)
    cent = compute_centrality(dag, results["consensus"])
    assert len(cent) == 28
    iso_dates = {iso: kp_ingestion.isolates[iso].collection_date for iso in kp_ingestion.isolates}
    cands = index_case_candidates(cent, iso_dates)
    assert len(cands) == results["consensus"].n_clusters


def test_bootstrap_stability(kp_config, kp_ingestion):
    g, results = _setup(kp_config, kp_ingestion)
    boot = bootstrap_stability(
        kp_ingestion.snp, kp_ingestion.isolates, results["consensus"], kp_config
    )
    assert len(boot) == 28
    summary = cluster_stability_summary(boot)
    assert "confidence_grade" in summary.columns


def test_bootstrap_reproducible(kp_config, kp_ingestion):
    g, results = _setup(kp_config, kp_ingestion)
    b1 = bootstrap_stability(
        kp_ingestion.snp, kp_ingestion.isolates, results["consensus"], kp_config
    )
    b2 = bootstrap_stability(
        kp_ingestion.snp, kp_ingestion.isolates, results["consensus"], kp_config
    )
    s1 = {b.isolate_id: b.stability_score for b in b1}
    s2 = {b.isolate_id: b.stability_score for b in b2}
    assert s1 == s2
