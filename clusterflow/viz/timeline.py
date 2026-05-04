"""Phase 6.1 Figure 3 — epi-genomic timeline scatter."""

from __future__ import annotations

from pathlib import Path

import igraph as ig
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from clusterflow.models import PipelineResult


def plot_epi_timeline(
    result: PipelineResult, dag: ig.Graph, output_dir: str | Path
) -> Path:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    consensus = result.consensus
    isolates = result.isolates

    by_cluster: dict[int, list[str]] = {}
    for iso, cid in consensus.assignments.items():
        by_cluster.setdefault(cid, []).append(iso)

    cluster_order = sorted(
        by_cluster,
        key=lambda c: min(isolates[i].collection_date for i in by_cluster[c]),
    )
    palette = sns.color_palette("tab20", n_colors=len(cluster_order))

    fig, ax = plt.subplots(figsize=(11, max(4, len(cluster_order) * 0.5)))
    cid_to_y = {c: i for i, c in enumerate(cluster_order)}
    idx_cases = {c.index_case_candidate for c in result.transmission_clusters}

    # Cluster shading bands
    for cid in cluster_order:
        members = by_cluster[cid]
        d_min = min(isolates[i].collection_date for i in members)
        d_max = max(isolates[i].collection_date for i in members)
        ax.axhspan(
            cid_to_y[cid] - 0.4, cid_to_y[cid] + 0.4,
            xmin=0, xmax=1,
            color=palette[cid_to_y[cid]],
            alpha=0.10,
        )
        ax.plot(
            [d_min, d_max],
            [cid_to_y[cid], cid_to_y[cid]],
            color=palette[cid_to_y[cid]],
            alpha=0.6,
            lw=2,
            zorder=1,
        )

    # DAG arcs (downsample if many)
    dag_edges = [(e.source, e.target) for e in dag.es]
    if len(dag_edges) <= 200:
        for s, t in dag_edges:
            iso_s = dag.vs[s]["isolate_id"]
            iso_t = dag.vs[t]["isolate_id"]
            cid_s = consensus.assignments.get(iso_s)
            cid_t = consensus.assignments.get(iso_t)
            if cid_s is None or cid_t is None or cid_s != cid_t:
                continue
            ax.plot(
                [isolates[iso_s].collection_date, isolates[iso_t].collection_date],
                [cid_to_y[cid_s], cid_to_y[cid_t]],
                color="grey",
                alpha=0.18,
                lw=0.6,
                zorder=0,
            )

    for iso, cid in consensus.assignments.items():
        ax.scatter(
            isolates[iso].collection_date,
            cid_to_y[cid],
            color=palette[cid_to_y[cid]],
            edgecolors="black" if iso in idx_cases else "none",
            linewidths=1.2 if iso in idx_cases else 0,
            s=80 if iso in idx_cases else 40,
            zorder=3,
        )

    ax.set_yticks(list(cid_to_y.values()))
    ax.set_yticklabels([f"cluster {c}" for c in cluster_order])
    ax.set_xlabel("collection date")
    ax.set_title(f"{result.project_name} — epi-genomic timeline")
    fig.autofmt_xdate()
    fig.tight_layout()

    png = out / "epi_timeline_scatter.png"
    svg = out / "epi_timeline_scatter.svg"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    return png
