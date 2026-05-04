"""End-to-end pipeline test on the K. pneumoniae fixtures."""

from __future__ import annotations

import time
from pathlib import Path

from sklearn.metrics import adjusted_rand_score

from clusterflow.pipeline import run_pipeline


def test_full_pipeline_runs_under_60s(kp_config, kp_truth):
    t0 = time.perf_counter()
    result = run_pipeline(kp_config)
    elapsed = time.perf_counter() - t0

    assert result.n_isolates == 28
    assert result.consensus.n_clusters >= 2

    iso = list(result.consensus.assignments)
    pred = [result.consensus.assignments[i] for i in iso]
    truth = [kp_truth[i] for i in iso]
    ari = adjusted_rand_score(truth, pred)
    assert ari > 0.7, f"consensus ARI vs truth too low: {ari:.2f}"

    out = Path(kp_config.output_dir)
    assert (out / "graph" / "transmission_graph.graphml").exists()
    assert (out / "graph" / "graph_summary.json").exists()
    assert (out / "clusters" / "cluster_assignments.csv").exists()
    assert (out / "analysis" / "transmission_dag.graphml").exists()
    assert (out / "analysis" / "centrality_scores.csv").exists()
    assert (out / "analysis" / "bootstrap_stability.csv").exists()
    assert (out / "figures" / "snp_heatmap.png").exists()
    assert (out / "figures" / "minimum_spanning_tree.png").exists()
    assert (out / "figures" / "epi_timeline_scatter.png").exists()
    assert (out / "figures" / "cluster_comparison_grid.png").exists()
    assert (out / "figures" / "bootstrap_stability_terrain.png").exists()
    assert (out / "pipeline_summary.json").exists()

    # Performance gate from the plan: <60s on a 4-core laptop
    assert elapsed < 60, f"pipeline too slow: {elapsed:.1f}s"


def test_pipeline_reproducible(kp_config):
    a = run_pipeline(kp_config, run_viz=False)
    b = run_pipeline(kp_config, run_viz=False)
    assert a.consensus.assignments == b.consensus.assignments
