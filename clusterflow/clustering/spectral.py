"""Phase 4.3 — spectral clustering with eigengap-driven k selection."""

from __future__ import annotations

import logging

import igraph as ig
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

from clusterflow.models import ClusterAssignment

log = logging.getLogger(__name__)


def _affinity_matrix(g: ig.Graph) -> np.ndarray:
    n = g.vcount()
    A = np.zeros((n, n), dtype=float)
    if g.ecount() == 0:
        return A
    weights = np.exp(-np.asarray(g.es["composite_weight"], dtype=float))
    for e, w in zip(g.es, weights):
        i, j = e.tuple
        A[i, j] = w
        A[j, i] = w
    return A


def _eigengap_k(eigvals: np.ndarray, max_k: int = 12) -> int:
    """Return the k that maximises consecutive-eigenvalue gap (skip k=1)."""
    if len(eigvals) < 3:
        return max(1, len(eigvals))
    gaps = np.diff(eigvals[: max_k + 1])
    # gaps[0] = eig1-eig0 → corresponds to k=1 vs k=2; we want k>=2
    if len(gaps) <= 1:
        return 2
    return int(np.argmax(gaps[1:]) + 2)


def spectral_clusters(
    g: ig.Graph,
    k: int | str = "auto",
    random_state: int = 42,
) -> ClusterAssignment:
    """Eigengap → k → KMeans on Laplacian eigenspace, with silhouette refinement."""
    ids = list(g.vs["isolate_id"])
    n = g.vcount()
    if n == 0:
        return ClusterAssignment(method="spectral", assignments={}, n_clusters=0)
    if g.ecount() == 0:
        return ClusterAssignment(
            method="spectral",
            assignments={iso: i for i, iso in enumerate(ids)},
            n_clusters=n,
        )

    A = _affinity_matrix(g)
    d = A.sum(axis=1)
    d_safe = np.where(d > 0, d, 1.0)
    D_inv_sqrt = 1.0 / np.sqrt(d_safe)
    # Symmetric normalised Laplacian: L_sym = I - D^{-1/2} A D^{-1/2}
    L_sym = np.eye(n) - (D_inv_sqrt[:, None] * A * D_inv_sqrt[None, :])
    # Force symmetry against round-off, then eigendecompose
    L_sym = (L_sym + L_sym.T) / 2.0
    eigvals, eigvecs = np.linalg.eigh(L_sym)

    if k == "auto":
        k_try = _eigengap_k(eigvals, max_k=min(12, n - 1))
    else:
        k_try = int(k)
    k_try = max(2, min(k_try, n - 1))

    candidates = sorted({max(2, k_try - 1), k_try, min(n - 1, k_try + 1)})

    best_k = k_try
    best_score = -np.inf
    best_labels: np.ndarray | None = None
    for kc in candidates:
        emb = eigvecs[:, :kc]
        # Row-normalise (Ng-Jordan-Weiss)
        norms = np.linalg.norm(emb, axis=1, keepdims=True)
        norms = np.where(norms > 0, norms, 1.0)
        emb = emb / norms
        labels = KMeans(
            n_clusters=kc,
            random_state=random_state,
            n_init=10,
        ).fit_predict(emb)
        if len(set(labels)) < 2:
            continue
        try:
            score = float(silhouette_score(emb, labels))
        except Exception:
            score = -np.inf
        if score > best_score:
            best_score = score
            best_k = kc
            best_labels = labels

    if best_labels is None:
        best_labels = np.zeros(n, dtype=int)
        best_k = 1
    log.info("spectral: k=%d silhouette=%.3f", best_k, best_score)

    assignments = {ids[i]: int(c) for i, c in enumerate(best_labels)}
    return ClusterAssignment(
        method="spectral",
        assignments=assignments,
        n_clusters=len(set(best_labels)),
    )
