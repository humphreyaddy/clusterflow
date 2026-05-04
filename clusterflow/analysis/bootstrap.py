"""Phase 5.3 — bootstrap confidence intervals via SNP-distance perturbation."""

from __future__ import annotations

import logging
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from clusterflow.clustering.consensus import _align_labels
from clusterflow.clustering.leiden import leiden_clusters
from clusterflow.config import ClusterFlowConfig
from clusterflow.graph.constructor import GraphConstructor
from clusterflow.models import (
    BootstrapResult,
    ClusterAssignment,
    Isolate,
    SNPMatrix,
)

log = logging.getLogger(__name__)


def _classify(score: float) -> str:
    if score >= 0.90:
        return "stable"
    if score >= 0.70:
        return "uncertain"
    return "unstable"


def _one_replicate(
    snp: SNPMatrix,
    isolates: dict[str, Isolate],
    cfg: ClusterFlowConfig,
    seed: int,
    perturb_max: int = 2,
) -> dict[str, int]:
    rng = np.random.default_rng(seed)
    n = snp.n
    noise = rng.uniform(-perturb_max, perturb_max, size=(n, n))
    noise = (noise + noise.T) / 2.0
    np.fill_diagonal(noise, 0.0)
    perturbed = np.clip(snp.distances + noise, 0, None)
    np.fill_diagonal(perturbed, 0.0)
    perturbed = (perturbed + perturbed.T) / 2.0

    perturbed_snp = SNPMatrix(isolate_ids=snp.isolate_ids, distances=perturbed)
    g = GraphConstructor(cfg.thresholds).build(perturbed_snp, isolates)
    a = leiden_clusters(
        g,
        resolution=cfg.clustering.leiden_resolution,
        random_state=seed,
        n_jobs=1,
    )
    return a.assignments


def bootstrap_stability(
    snp: SNPMatrix,
    isolates: dict[str, Isolate],
    consensus: ClusterAssignment,
    cfg: ClusterFlowConfig,
) -> list[BootstrapResult]:
    n_boot = cfg.clustering.bootstrap_n
    if n_boot <= 0:
        return [
            BootstrapResult(
                isolate_id=iso,
                cluster_id=int(c),
                stability_score=1.0,
                classification="stable",
            )
            for iso, c in consensus.assignments.items()
        ]

    rng = np.random.default_rng(cfg.clustering.random_state)
    seeds = rng.integers(0, 2**31 - 1, size=n_boot).tolist()

    log.info("running %d bootstrap replicates...", n_boot)
    replicates = Parallel(n_jobs=cfg.distance_engine.n_jobs, backend="loky")(
        delayed(_one_replicate)(snp, isolates, cfg, int(s)) for s in seeds
    )

    # Align each replicate to the consensus labels via Hungarian
    aligned_assignments: list[dict[str, int]] = []
    for rep in replicates:
        mapping = _align_labels(consensus.assignments, rep)
        aligned_assignments.append({iso: mapping[c] for iso, c in rep.items()})

    out: list[BootstrapResult] = []
    for iso, gold in consensus.assignments.items():
        votes = [a.get(iso) for a in aligned_assignments if iso in a]
        if not votes:
            score = 0.0
        else:
            score = sum(1 for v in votes if v == gold) / len(votes)
        out.append(
            BootstrapResult(
                isolate_id=iso,
                cluster_id=int(gold),
                stability_score=float(score),
                classification=_classify(score),
            )
        )
    return out


def cluster_stability_summary(
    bootstrap: list[BootstrapResult],
) -> pd.DataFrame:
    if not bootstrap:
        return pd.DataFrame()
    df = pd.DataFrame([b.model_dump() for b in bootstrap])
    rows = []
    for cid, sub in df.groupby("cluster_id"):
        mean_score = float(sub["stability_score"].mean())
        if mean_score >= 0.90:
            grade = "A"
        elif mean_score >= 0.75:
            grade = "B"
        else:
            grade = "C"
        rows.append(
            {
                "cluster_id": int(cid),
                "n_isolates": len(sub),
                "mean_stability": mean_score,
                "n_stable": int((sub["classification"] == "stable").sum()),
                "n_uncertain": int((sub["classification"] == "uncertain").sum()),
                "n_unstable": int((sub["classification"] == "unstable").sum()),
                "confidence_grade": grade,
            }
        )
    return pd.DataFrame(rows).sort_values("cluster_id").reset_index(drop=True)


def save_bootstrap(
    bootstrap: list[BootstrapResult], output_dir: str | Path
) -> dict[str, str]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    a = out / "bootstrap_stability.csv"
    b = out / "cluster_stability_summary.csv"
    pd.DataFrame([r.model_dump() for r in bootstrap]).to_csv(a, index=False)
    cluster_stability_summary(bootstrap).to_csv(b, index=False)
    return {"per_isolate": str(a), "per_cluster": str(b)}
