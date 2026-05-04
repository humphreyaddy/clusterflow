"""Tests for clusterflow.graph."""

from __future__ import annotations

from datetime import date

import pytest

from clusterflow.config import ThresholdsConfig
from clusterflow.graph import GraphConstructor, graph_summary, save_graph, load_graph


def test_graph_built_from_kp(kp_config, kp_ingestion):
    g = GraphConstructor(kp_config.thresholds).build(kp_ingestion.snp, kp_ingestion.isolates)
    assert g.vcount() == 28
    assert g.ecount() > 10  # within-cluster edges should exist
    for attr in ("snp_distance", "day_delta", "mlst_mismatch", "composite_weight"):
        assert attr in g.es.attributes()


def test_graph_round_trip(kp_config, kp_ingestion, tmp_path):
    g = GraphConstructor(kp_config.thresholds).build(kp_ingestion.snp, kp_ingestion.isolates)
    save_graph(g, tmp_path)
    g2 = load_graph(tmp_path / "transmission_graph.pkl")
    assert g.vcount() == g2.vcount()
    assert g.ecount() == g2.ecount()


def test_summary_stats(kp_config, kp_ingestion):
    g = GraphConstructor(kp_config.thresholds).build(kp_ingestion.snp, kp_ingestion.isolates)
    s = graph_summary(g)
    assert s["n_vertices"] == 28
    assert s["weight_max"] >= 0
