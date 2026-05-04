"""Phase 6.1 Figure 5 — bootstrap stability terrain (MDS + KDE)."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.manifold import MDS

from clusterflow.models import PipelineResult


def plot_bootstrap_terrain(
    result: PipelineResult, output_dir: str | Path
) -> Path:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    snp = result.snp_matrix
    if snp.n < 3:
        # Trivial fallback
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.text(0.5, 0.5, "Too few isolates for terrain", ha="center", va="center")
        ax.axis("off")
        png = out / "bootstrap_stability_terrain.png"
        fig.savefig(png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        return png

    try:
        mds = MDS(
            n_components=2,
            dissimilarity="precomputed",
            random_state=42,
            n_init=4,
            normalized_stress="auto",
        )
    except TypeError:
        mds = MDS(
            n_components=2,
            dissimilarity="precomputed",
            random_state=42,
            n_init=4,
        )
    coords = mds.fit_transform(snp.distances)

    boot = {b.isolate_id: b.stability_score for b in result.bootstrap}
    consensus = result.consensus
    palette = sns.color_palette("tab20", n_colors=max(consensus.n_clusters, 1))

    fig, ax = plt.subplots(figsize=(10, 8))

    # KDE-style scatter weighted by stability
    for iso_idx, iso in enumerate(snp.isolate_ids):
        x, y = coords[iso_idx]
        cid = consensus.assignments.get(iso, -1)
        score = boot.get(iso, 1.0)
        ax.scatter(
            x,
            y,
            color=palette[cid % len(palette)] if cid >= 0 else "grey",
            edgecolors="black" if score < 0.7 else "none",
            linewidths=1.5 if score < 0.7 else 0,
            s=80 + 200 * score,
            alpha=0.7 + 0.3 * score,
            zorder=2,
        )

    sns.kdeplot(
        x=coords[:, 0],
        y=coords[:, 1],
        fill=True,
        cmap="YlOrRd",
        alpha=0.4,
        thresh=0.05,
        levels=8,
        ax=ax,
        zorder=0,
    )

    ax.set_title(
        f"{result.project_name} — bootstrap stability terrain (MDS + KDE)",
        fontsize=12,
    )
    ax.set_xlabel("MDS axis 1")
    ax.set_ylabel("MDS axis 2")
    fig.tight_layout()

    png = out / "bootstrap_stability_terrain.png"
    svg = out / "bootstrap_stability_terrain.svg"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    return png
