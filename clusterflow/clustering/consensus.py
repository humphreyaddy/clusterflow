"""Phase 4.5 — consensus clustering via Hungarian label alignment + majority vote."""

from __future__ import annotations

import logging
from collections import Counter

import numpy as np
from scipy.optimize import linear_sum_assignment

from clusterflow.models import ClusterAssignment

log = logging.getLogger(__name__)


def _align_labels(
    reference: dict[str, int], other: dict[str, int]
) -> dict[int, int]:
    """Solve the Hungarian label-alignment problem.

    Returns a dict mapping ``other`` cluster IDs → ``reference`` cluster IDs
    that maximises overlap. Unmapped ``other`` IDs (i.e. clusters in ``other``
    that have no good reference match) are passed through with offset to
    avoid collisions.
    """
    ref_ids = sorted(set(reference.values()))
    other_ids = sorted(set(other.values()))
    if not ref_ids or not other_ids:
        return {oid: oid for oid in other_ids}

    n_rows = len(other_ids)
    n_cols = len(ref_ids)
    overlap = np.zeros((n_rows, n_cols), dtype=int)
    ref_idx = {c: i for i, c in enumerate(ref_ids)}
    oth_idx = {c: i for i, c in enumerate(other_ids)}

    common = set(reference) & set(other)
    for iso in common:
        overlap[oth_idx[other[iso]], ref_idx[reference[iso]]] += 1

    # Pad to square if needed for linear_sum_assignment
    size = max(n_rows, n_cols)
    cost = np.zeros((size, size), dtype=int)
    cost[:n_rows, :n_cols] = -overlap  # maximise overlap → minimise -overlap
    row_ind, col_ind = linear_sum_assignment(cost)

    mapping: dict[int, int] = {}
    next_free = max(ref_ids) + 1
    for r, c in zip(row_ind, col_ind):
        if r >= n_rows:
            continue
        if c < n_cols and overlap[r, c] > 0:
            mapping[other_ids[r]] = ref_ids[c]
        else:
            mapping[other_ids[r]] = next_free
            next_free += 1
    # Any other_id we missed (shouldn't happen, but be safe)
    for oid in other_ids:
        if oid not in mapping:
            mapping[oid] = next_free
            next_free += 1
    return mapping


def consensus_assignment(
    assignments: list[ClusterAssignment],
) -> ClusterAssignment:
    """Hungarian-align all methods to the first, then majority vote per isolate."""
    if not assignments:
        raise ValueError("need at least one ClusterAssignment for consensus")
    if len(assignments) == 1:
        a = assignments[0]
        return ClusterAssignment(
            method="consensus",
            assignments=dict(a.assignments),
            n_clusters=a.n_clusters,
            ambiguous={iso: False for iso in a.assignments},
            agreement_score={iso: 1.0 for iso in a.assignments},
        )

    reference = assignments[0].assignments
    aligned = [reference]
    for a in assignments[1:]:
        mapping = _align_labels(reference, a.assignments)
        aligned.append({iso: mapping[c] for iso, c in a.assignments.items()})

    isolates = sorted({iso for a in aligned for iso in a})
    consensus: dict[str, int] = {}
    ambiguous: dict[str, bool] = {}
    agreement: dict[str, float] = {}
    n_methods = len(aligned)

    for iso in isolates:
        votes = [a[iso] for a in aligned if iso in a]
        if not votes:
            continue
        counts = Counter(votes)
        winner, top = counts.most_common(1)[0]
        score = top / n_methods
        consensus[iso] = winner
        agreement[iso] = score
        ambiguous[iso] = score < (2 / 3)

    n_full = sum(1 for s in agreement.values() if s == 1.0)
    n_amb = sum(1 for v in ambiguous.values() if v)
    log.info(
        "consensus: %d isolates fully agreed, %d ambiguous (<2/3 agreement)",
        n_full,
        n_amb,
    )

    return ClusterAssignment(
        method="consensus",
        assignments=consensus,
        n_clusters=len(set(consensus.values())),
        ambiguous=ambiguous,
        agreement_score=agreement,
    )
