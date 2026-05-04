"""Phase 5.2 — centrality scoring + index case identification."""

from __future__ import annotations

import logging
from datetime import date
from pathlib import Path

import igraph as ig
import numpy as np
import pandas as pd

from clusterflow.models import CentralityScores, ClusterAssignment

log = logging.getLogger(__name__)


def _normalise(arr: np.ndarray) -> np.ndarray:
    arr = np.asarray(arr, dtype=float)
    lo, hi = arr.min(), arr.max()
    if hi - lo < 1e-12:
        return np.zeros_like(arr)
    return (arr - lo) / (hi - lo)


def compute_centrality(
    dag: ig.Graph,
    consensus: ClusterAssignment,
) -> list[CentralityScores]:
    n = dag.vcount()
    if n == 0:
        return []
    ids: list[str] = list(dag.vs["isolate_id"])
    if dag.ecount() > 0:
        weights = list(dag.es["composite_weight"])
        weights = [max(w, 1e-6) for w in weights]
        between = np.asarray(
            dag.betweenness(directed=True, weights=weights), dtype=float
        )
        close = np.asarray(
            dag.closeness(mode="all", weights=weights, normalized=True),
            dtype=float,
        )
    else:
        between = np.zeros(n)
        close = np.zeros(n)
    indeg = np.asarray(dag.indegree(), dtype=int)
    outdeg = np.asarray(dag.outdegree(), dtype=int)

    # Normalise within each cluster, then compute composite risk score
    cluster_ids = np.asarray(
        [consensus.assignments.get(i, -1) for i in ids], dtype=int
    )
    risk = np.zeros(n, dtype=float)
    between_n = np.zeros(n, dtype=float)
    out_n = np.zeros(n, dtype=float)
    close_n = np.zeros(n, dtype=float)
    for cid in np.unique(cluster_ids):
        mask = cluster_ids == cid
        if mask.sum() == 0:
            continue
        between_n[mask] = _normalise(between[mask])
        out_n[mask] = _normalise(outdeg[mask].astype(float))
        close_n[mask] = _normalise(close[mask])
        risk[mask] = (
            0.5 * between_n[mask] + 0.3 * out_n[mask] + 0.2 * close_n[mask]
        )

    return [
        CentralityScores(
            isolate_id=ids[i],
            cluster_id=int(cluster_ids[i]),
            betweenness=float(between[i]),
            in_degree=int(indeg[i]),
            out_degree=int(outdeg[i]),
            closeness=float(close[i]) if not np.isnan(close[i]) else 0.0,
            transmission_risk_score=float(risk[i]),
        )
        for i in range(n)
    ]


def index_case_candidates(
    centrality: list[CentralityScores],
    isolates_dates: dict[str, date],
) -> pd.DataFrame:
    """One row per cluster with the highest-risk isolate; flag temporal mismatches."""
    if not centrality:
        return pd.DataFrame()
    df = pd.DataFrame([c.model_dump() for c in centrality])
    rows = []
    for cid, sub in df.groupby("cluster_id"):
        top = sub.sort_values("transmission_risk_score", ascending=False).iloc[0]
        cluster_isos = sub["isolate_id"].tolist()
        cluster_dates = {iso: isolates_dates[iso] for iso in cluster_isos}
        earliest = min(cluster_dates.values())
        rank = sorted(cluster_dates, key=lambda x: cluster_dates[x]).index(
            top["isolate_id"]
        )
        temporal_mismatch = rank > max(0, len(cluster_isos) // 4)
        if temporal_mismatch:
            log.warning(
                "cluster %s: index case candidate %s ranked %d/%d by date "
                "(earliest=%s) — review needed",
                cid,
                top["isolate_id"],
                rank + 1,
                len(cluster_isos),
                earliest,
            )
        rows.append(
            {
                "cluster_id": int(cid),
                "index_case_candidate": top["isolate_id"],
                "transmission_risk_score": float(top["transmission_risk_score"]),
                "betweenness": float(top["betweenness"]),
                "out_degree": int(top["out_degree"]),
                "closeness": float(top["closeness"]),
                "earliest_date": earliest.isoformat(),
                "candidate_date_rank": rank + 1,
                "n_isolates_in_cluster": len(cluster_isos),
                "temporal_mismatch": bool(temporal_mismatch),
            }
        )
    return pd.DataFrame(rows).sort_values("cluster_id").reset_index(drop=True)


def save_centrality(
    centrality: list[CentralityScores],
    candidates: pd.DataFrame,
    output_dir: str | Path,
) -> dict[str, str]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    cent_csv = out / "centrality_scores.csv"
    cand_csv = out / "index_case_candidates.csv"
    pd.DataFrame([c.model_dump() for c in centrality]).to_csv(cent_csv, index=False)
    candidates.to_csv(cand_csv, index=False)
    return {"centrality": str(cent_csv), "index_cases": str(cand_csv)}
