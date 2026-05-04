"""Phase 2.2 — CPU backend (joblib-parallel Hamming distance)."""

from __future__ import annotations

import logging
import os
from pathlib import Path

import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm

from clusterflow.distance.backend import DistanceBackend
from clusterflow.models import SNPMatrix

log = logging.getLogger(__name__)


def _read_fasta(path: Path) -> tuple[str, str]:
    """Return (id, sequence). Single-record FASTAs only."""
    seq_chunks: list[str] = []
    iso_id: str | None = None
    with path.open() as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if iso_id is None:
                    iso_id = line[1:].split()[0]
                else:
                    break  # ignore additional records, take first
            else:
                seq_chunks.append(line.upper())
    if iso_id is None:
        iso_id = path.stem
    return iso_id, "".join(seq_chunks)


def _encode(seq: str) -> np.ndarray:
    return np.frombuffer(seq.encode("ascii"), dtype=np.uint8)


def _ambiguous_mask(arr: np.ndarray) -> np.ndarray:
    # Mask N and gap characters
    return (arr != ord("N")) & (arr != ord("-"))


def _row_distances(i: int, M: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Distances from row i to all rows >= i."""
    n = M.shape[0]
    out = np.zeros(n, dtype=np.int32)
    if i + 1 == n:
        return out
    diff = (M[i + 1 :] != M[i]) & valid[i + 1 :] & valid[i]
    out[i + 1 :] = diff.sum(axis=1)
    return out


class CPUBackend(DistanceBackend):
    name = "cpu"

    def compute(self, fasta_paths: list[Path], n_jobs: int = -1) -> SNPMatrix:
        if not fasta_paths:
            raise ValueError("no FASTA files supplied")
        records = [_read_fasta(p) for p in fasta_paths]
        ids = [r[0] for r in records]
        seqs = [r[1] for r in records]

        lens = {len(s) for s in seqs}
        if len(lens) != 1:
            raise ValueError(
                f"sequences differ in length: {sorted(lens)} — align first"
            )
        if lens.pop() == 0:
            raise ValueError("all sequences are empty")

        M = np.stack([_encode(s) for s in seqs])
        valid = _ambiguous_mask(M)

        n = M.shape[0]
        if n_jobs == -1:
            n_jobs = os.cpu_count() or 1

        rows = Parallel(n_jobs=n_jobs, backend="loky")(
            delayed(_row_distances)(i, M, valid)
            for i in tqdm(range(n - 1), desc="SNP distances", disable=n < 100)
        )

        D = np.zeros((n, n), dtype=float)
        for i, row in enumerate(rows):
            D[i, i + 1 :] = row[i + 1 :]
        D = D + D.T  # symmetrise
        return SNPMatrix(isolate_ids=ids, distances=D)
