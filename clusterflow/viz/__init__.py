"""Phase 6 — static figures + interactive Dash dashboard."""

from __future__ import annotations

import logging
from pathlib import Path

import igraph as ig

from clusterflow.viz.comparison import plot_cluster_comparison
from clusterflow.viz.heatmap import plot_snp_heatmap
from clusterflow.viz.mst import plot_minimum_spanning_tree
from clusterflow.viz.terrain import plot_bootstrap_terrain
from clusterflow.viz.timeline import plot_epi_timeline

log = logging.getLogger(__name__)


def generate_all_figures(
    result, graph: ig.Graph, dag: ig.Graph, output_dir: str | Path
) -> dict[str, str]:
    """Render every static figure for a completed PipelineResult.

    ``graph`` is the undirected weighted graph; ``dag`` is the directed
    temporal DAG built in Phase 5.1.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    figs: dict[str, str] = {}
    log.info("rendering snp heatmap")
    figs["snp_heatmap"] = str(plot_snp_heatmap(result, out))
    log.info("rendering minimum spanning tree")
    figs["mst"] = str(plot_minimum_spanning_tree(result, graph, out))
    log.info("rendering epi timeline")
    figs["timeline"] = str(plot_epi_timeline(result, dag, out))
    log.info("rendering cluster comparison")
    figs["comparison"] = str(plot_cluster_comparison(result, graph, out))
    log.info("rendering bootstrap terrain")
    figs["terrain"] = str(plot_bootstrap_terrain(result, out))
    return figs


__all__ = [
    "generate_all_figures",
    "plot_cluster_comparison",
    "plot_bootstrap_terrain",
    "plot_epi_timeline",
    "plot_minimum_spanning_tree",
    "plot_snp_heatmap",
]
