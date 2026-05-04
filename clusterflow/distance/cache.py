"""Phase 2.5 — SHA-256 keyed npz cache for SNP matrices."""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path

import numpy as np

from clusterflow.models import SNPMatrix

log = logging.getLogger(__name__)


def _hash_inputs(paths: list[Path]) -> str:
    h = hashlib.sha256()
    for p in sorted(paths, key=str):
        h.update(str(p.resolve()).encode())
        h.update(b"\0")
        with p.open("rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 16), b""):
                h.update(chunk)
    return h.hexdigest()[:16]


def cache_path(output_dir: Path, key: str) -> Path:
    return output_dir / "cache" / f"snp_matrix_{key}.npz"


def load_cached(output_dir: Path, fasta_paths: list[Path]) -> SNPMatrix | None:
    key = _hash_inputs(fasta_paths)
    p = cache_path(output_dir, key)
    if not p.exists():
        return None
    log.info("cache hit: %s", p)
    data = np.load(p, allow_pickle=False)
    ids = [str(x) for x in data["ids"]]
    D = data["distances"].astype(float)
    return SNPMatrix(isolate_ids=ids, distances=D)


def save_cached(output_dir: Path, fasta_paths: list[Path], snp: SNPMatrix) -> Path:
    key = _hash_inputs(fasta_paths)
    p = cache_path(output_dir, key)
    p.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        p,
        ids=np.asarray(snp.isolate_ids, dtype="<U64"),
        distances=snp.distances,
    )
    log.info("cached SNP matrix: %s", p)
    return p
