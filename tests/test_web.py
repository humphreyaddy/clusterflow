"""Tests for clusterflow.web (HTML report)."""

from __future__ import annotations

from pathlib import Path

from clusterflow.pipeline import run_pipeline
from clusterflow.web.report import render_report


def test_pipeline_emits_report(kp_config, tmp_path):
    cfg = kp_config.model_copy(update={"output_dir": tmp_path})
    run_pipeline(cfg)
    report = tmp_path / "report.html"
    assert report.exists(), "pipeline did not emit report.html"
    body = report.read_text()
    assert "<title>ClusterFlow report" in body
    assert cfg.project_name in body
    # Each figure should be embedded as a base64 PNG
    for stem in [
        "snp_heatmap",
        "minimum_spanning_tree",
        "epi_timeline_scatter",
        "cluster_comparison_grid",
        "bootstrap_stability_terrain",
    ]:
        assert f'id="fig-{stem}"' in body, f"figure {stem} missing from report"
        assert "data:image/png;base64," in body
    # Cluster table + KPIs should render
    assert "Clusters" in body
    assert "kpi" in body


def test_render_report_idempotent(kp_config, tmp_path):
    cfg = kp_config.model_copy(update={"output_dir": tmp_path})
    result = run_pipeline(cfg, run_viz=True)
    a = render_report(result, tmp_path)
    b = render_report(result, tmp_path)
    assert a == b == tmp_path / "report.html"
    # Same content on both renders
    assert a.read_text() == b.read_text()


def test_report_safe_with_no_bootstrap(kp_config, tmp_path):
    cfg = kp_config.model_copy(
        update={
            "output_dir": tmp_path,
            "clustering": kp_config.clustering.model_copy(update={"bootstrap_n": 0}),
        }
    )
    run_pipeline(cfg)
    body = (tmp_path / "report.html").read_text()
    assert "Bootstrap cluster stability" in body
