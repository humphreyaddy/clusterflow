"""Tests for clusterflow.clustering."""

from __future__ import annotations

from sklearn.metrics import adjusted_rand_score

from clusterflow.clustering import (
    consensus_assignment,
    leiden_clusters,
    run_all,
    snp_chain_clusters,
    spectral_clusters,
)
from clusterflow.graph import GraphConstructor


def _build(kp_config, kp_ingestion):
    return GraphConstructor(kp_config.thresholds).build(
        kp_ingestion.snp, kp_ingestion.isolates
    )


def test_snp_chain_clusters(kp_config, kp_ingestion, kp_truth):
    g = _build(kp_config, kp_ingestion)
    a = snp_chain_clusters(g, kp_config.thresholds.snp_cutoff, kp_config.thresholds.day_cutoff)
    assert a.method == "snp_chains"
    assert 3 <= a.n_clusters <= 8
    iso = list(a.assignments)
    pred = [a.assignments[i] for i in iso]
    truth = [kp_truth[i] for i in iso]
    ari = adjusted_rand_score(truth, pred)
    assert ari > 0.5, f"snp_chain ARI too low: {ari:.2f}"


def test_leiden_clusters(kp_config, kp_ingestion, kp_truth):
    g = _build(kp_config, kp_ingestion)
    a = leiden_clusters(g, resolution=1.0, random_state=42, n_jobs=1)
    assert a.method == "leiden"
    assert a.n_clusters >= 2
    iso = list(a.assignments)
    pred = [a.assignments[i] for i in iso]
    truth = [kp_truth[i] for i in iso]
    assert adjusted_rand_score(truth, pred) > 0.5


def test_spectral_clusters(kp_config, kp_ingestion, kp_truth):
    g = _build(kp_config, kp_ingestion)
    a = spectral_clusters(g, k="auto", random_state=42)
    assert a.method == "spectral"
    assert a.n_clusters >= 2


def test_run_all_consensus(kp_config, kp_ingestion, kp_truth):
    g = _build(kp_config, kp_ingestion)
    results = run_all(g, kp_config.thresholds, kp_config.clustering)
    assert set(results) >= {"snp_chains", "leiden", "spectral", "consensus"}
    cons = results["consensus"]
    assert cons.ambiguous is not None
    iso = list(cons.assignments)
    pred = [cons.assignments[i] for i in iso]
    truth = [kp_truth[i] for i in iso]
    assert adjusted_rand_score(truth, pred) > 0.5


def test_consensus_alignment_simple():
    a = type("C", (), {})()
    from clusterflow.models import ClusterAssignment

    a1 = ClusterAssignment(method="m1", assignments={"x": 0, "y": 0, "z": 1}, n_clusters=2)
    # Same partition but with shuffled labels
    a2 = ClusterAssignment(method="m2", assignments={"x": 5, "y": 5, "z": 9}, n_clusters=2)
    cons = consensus_assignment([a1, a2])
    # After alignment, consensus should match a1 perfectly
    assert cons.assignments == a1.assignments
    assert all(s == 1.0 for s in cons.agreement_score.values())


def test_reproducible(kp_config, kp_ingestion):
    g = _build(kp_config, kp_ingestion)
    r1 = leiden_clusters(g, resolution=1.0, random_state=42, n_jobs=1)
    r2 = leiden_clusters(g, resolution=1.0, random_state=42, n_jobs=1)
    assert r1.assignments == r2.assignments
