"""Phase 2.4 — GPU backend (cupy).

Optional. Only loadable if ``cupy`` is installed and a CUDA device exists.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from clusterflow.distance.backend import DistanceBackend
from clusterflow.distance.cpu import _ambiguous_mask, _encode, _read_fasta
from clusterflow.models import SNPMatrix

log = logging.getLogger(__name__)


class GPUBackend(DistanceBackend):
    name = "gpu"

    def __init__(self) -> None:
        try:
            import cupy  # noqa: F401
        except ImportError as e:
            raise RuntimeError(
                "GPU backend requested but cupy is not installed"
            ) from e

    def compute(self, fasta_paths: list[Path], n_jobs: int = -1) -> SNPMatrix:
        import cupy as cp

        records = [_read_fasta(p) for p in fasta_paths]
        ids = [r[0] for r in records]
        seqs = [r[1] for r in records]
        lens = {len(s) for s in seqs}
        if len(lens) != 1:
            raise ValueError(f"sequences differ in length: {sorted(lens)}")

        M = cp.asarray(np.stack([_encode(s) for s in seqs]))
        valid = cp.asarray(_ambiguous_mask(np.stack([_encode(s) for s in seqs])))
        diff = (M[:, None, :] != M[None, :, :]) & valid[:, None, :] & valid[None, :, :]
        D = diff.sum(axis=2).astype(cp.float64)
        D_cpu = cp.asnumpy(D)
        D_cpu = (D_cpu + D_cpu.T) / 2.0
        np.fill_diagonal(D_cpu, 0)
        return SNPMatrix(isolate_ids=ids, distances=D_cpu)
