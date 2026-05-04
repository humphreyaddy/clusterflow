"""Phase 6.1 Figure 4 — cluster comparison grid (one panel per method)."""

from __future__ import annotations

from pathlib import Path

import igraph as ig
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from clusterflow.models import PipelineResult


def plot_cluster_comparison(
    result: PipelineResult, g: ig.Graph, output_dir: str | Path
) -> Path:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    methods = list(result.cluster_assignments.keys()) + ["consensus"]
    layout = g.layout("fr") if g.ecount() else None
    coords = (
        np.asarray(layout.coords)
        if layout is not None
        else np.random.default_rng(0).normal(size=(g.vcount(), 2))
    )

    n = len(methods)
    cols = 2
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(11, 5 * rows))
    axes = np.array(axes).reshape(-1)

    for ax, method in zip(axes, methods):
        if method == "consensus":
            assignment = result.consensus
        else:
            assignment = result.cluster_assignments[method]
        ambiguous = (
            result.consensus.ambiguous if method == "consensus" else None
        )
        n_clusters = assignment.n_clusters
        palette = sns.color_palette("tab20", n_colors=max(n_clusters, 1))

        if g.ecount():
            for e in g.es:
                i, j = e.tuple
                ax.plot(
                    [coords[i, 0], coords[j, 0]],
                    [coords[i, 1], coords[j, 1]],
                    color="lightgrey",
                    lw=0.5,
                    alpha=0.5,
                    zorder=1,
                )

        for v in g.vs:
            iso = v["isolate_id"]
            cid = assignment.assignments.get(iso, -1)
            color = palette[cid % len(palette)] if cid >= 0 else "lightgrey"
            edge = "black"
            lw = 0.4
            if ambiguous and ambiguous.get(iso, False):
                edge = "red"
                lw = 1.5
            ax.scatter(
                coords[v.index, 0],
                coords[v.index, 1],
                color=color,
                edgecolors=edge,
                linewidths=lw,
                s=55,
                zorder=2,
            )

        ax.set_title(f"{method}  (n={n_clusters})", fontsize=11)
        ax.set_axis_off()

    for ax in axes[len(methods):]:
        ax.set_axis_off()

    fig.suptitle(
        f"{result.project_name} — cluster method comparison", fontsize=13
    )
    fig.tight_layout()

    png = out / "cluster_comparison_grid.png"
    svg = out / "cluster_comparison_grid.svg"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    return png
