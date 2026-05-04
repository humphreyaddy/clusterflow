"""Phase 6.1 Figure 2 — minimum spanning tree."""

from __future__ import annotations

from pathlib import Path

import igraph as ig
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from clusterflow.models import PipelineResult


_WARD_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]


def plot_minimum_spanning_tree(
    result: PipelineResult, g: ig.Graph, output_dir: str | Path
) -> Path:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    if g.ecount() == 0:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, "No edges to draw", ha="center", va="center")
        ax.axis("off")
        png = out / "minimum_spanning_tree.png"
        fig.savefig(png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        return png

    mst = g.spanning_tree(weights="composite_weight")
    layout = mst.layout("kk")
    coords = np.asarray(layout.coords)

    consensus = result.consensus
    cluster_palette = sns.color_palette("tab20", n_colors=max(consensus.n_clusters, 1))
    risk = {c.isolate_id: c.transmission_risk_score for c in result.centrality}
    idx_cases = {c.index_case_candidate for c in result.transmission_clusters}
    idx_cases.discard(None)

    wards = sorted({result.isolates[i].ward for i in g.vs["isolate_id"]})
    ward_marker = {w: _WARD_MARKERS[i % len(_WARD_MARKERS)] for i, w in enumerate(wards)}

    fig, ax = plt.subplots(figsize=(11, 9))

    # Edges
    for e in mst.es:
        i, j = e.tuple
        snp = e["snp_distance"]
        lw = max(0.3, 2.5 / max(1, snp))
        ax.plot(
            [coords[i, 0], coords[j, 0]],
            [coords[i, 1], coords[j, 1]],
            color="#888",
            lw=lw,
            alpha=0.7,
            zorder=1,
        )

    # Group nodes by ward to do one scatter per ward (matplotlib-safe)
    for ward in wards:
        idx = [v.index for v in mst.vs if result.isolates[v["isolate_id"]].ward == ward]
        if not idx:
            continue
        cids = [consensus.assignments[mst.vs[i]["isolate_id"]] for i in idx]
        risks = [risk.get(mst.vs[i]["isolate_id"], 0.0) for i in idx]
        ax.scatter(
            coords[idx, 0],
            coords[idx, 1],
            color=[cluster_palette[c % len(cluster_palette)] for c in cids],
            marker=ward_marker[ward],
            s=[60 + 240 * r for r in risks],
            edgecolors="black",
            linewidths=0.5,
            zorder=2,
            label=ward,
        )

    # Index-case stars (one combined scatter)
    star_idx = [v.index for v in mst.vs if v["isolate_id"] in idx_cases]
    if star_idx:
        ax.scatter(
            coords[star_idx, 0],
            coords[star_idx, 1],
            marker="*",
            s=420,
            color="red",
            edgecolors="black",
            linewidths=0.6,
            zorder=4,
            label="index case",
        )

    cluster_handles = [
        plt.Line2D(
            [0], [0], marker="o", linestyle="",
            color=cluster_palette[c % len(cluster_palette)],
            label=f"cluster {c}",
        )
        for c in sorted(set(consensus.assignments.values()))
    ]
    ward_handles = [
        plt.Line2D(
            [0], [0], marker=ward_marker[w], linestyle="", color="grey", label=w
        )
        for w in wards
    ]
    star_handle = plt.Line2D(
        [0], [0], marker="*", linestyle="", color="red",
        label="index case", markersize=12,
    )
    ax.legend(
        handles=cluster_handles + ward_handles + [star_handle],
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=8,
        ncol=1,
    )
    ax.set_title(f"{result.project_name} — Minimum Spanning Tree", fontsize=13)
    ax.set_axis_off()

    png = out / "minimum_spanning_tree.png"
    svg = out / "minimum_spanning_tree.svg"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    return png
