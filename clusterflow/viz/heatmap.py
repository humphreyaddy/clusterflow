"""Phase 6.1 Figure 1 — SNP distance heatmap with cluster annotations."""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from clusterflow.models import PipelineResult

log = logging.getLogger(__name__)


def _palette(n: int) -> list:
    pal = sns.color_palette("tab20", n_colors=max(n, 1))
    return [tuple(c) for c in pal]


def plot_snp_heatmap(result: PipelineResult, output_dir: str | Path) -> Path:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    snp = result.snp_matrix
    consensus = result.consensus
    isolates = result.isolates

    # Order isolates by cluster, then by collection date
    order = sorted(
        snp.isolate_ids,
        key=lambda i: (
            consensus.assignments.get(i, -1),
            isolates[i].collection_date,
        ),
    )
    idx = [snp.isolate_ids.index(i) for i in order]
    D = snp.distances[np.ix_(idx, idx)]
    df = pd.DataFrame(D, index=order, columns=order)

    cluster_ids = [consensus.assignments[i] for i in order]
    wards = [isolates[i].ward for i in order]
    sts = [isolates[i].sequence_type or "?" for i in order]
    cluster_colors = dict(zip(sorted(set(cluster_ids)), _palette(len(set(cluster_ids)))))
    ward_colors = dict(zip(sorted(set(wards)), _palette(len(set(wards)))))
    st_colors = dict(zip(sorted(set(sts)), _palette(len(set(sts)))))
    row_colors = pd.DataFrame(
        {
            "cluster": [cluster_colors[c] for c in cluster_ids],
            "ward": [ward_colors[w] for w in wards],
            "ST": [st_colors[s] for s in sts],
        },
        index=order,
    )

    g = sns.clustermap(
        df,
        row_cluster=False,
        col_cluster=False,
        cmap="viridis_r",
        row_colors=row_colors,
        col_colors=row_colors,
        figsize=(10, 9),
        cbar_kws={"label": "SNP distance"},
        xticklabels=False,
        yticklabels=True,
    )
    g.ax_heatmap.set_title(
        f"{result.project_name} — pairwise SNP distances",
        fontsize=12,
        pad=14,
    )
    g.ax_heatmap.set_yticklabels(
        g.ax_heatmap.get_yticklabels(), fontsize=6, rotation=0
    )

    png_path = out / "snp_heatmap.png"
    svg_path = out / "snp_heatmap.svg"
    g.savefig(png_path, dpi=300, bbox_inches="tight")
    g.savefig(svg_path, bbox_inches="tight")
    plt.close(g.fig)
    return png_path
