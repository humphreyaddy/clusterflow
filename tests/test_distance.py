"""Tests for clusterflow.distance (CPU backend + cache)."""

from __future__ import annotations

import numpy as np

from clusterflow.distance import CPUBackend, load_cached, save_cached, select_backend
from clusterflow.distance.backend import DistanceBackend
from clusterflow.models import SNPMatrix


def _write_fasta(tmp_path, name: str, seq: str):
    p = tmp_path / f"{name}.fasta"
    p.write_text(f">{name}\n{seq}\n")
    return p


def test_cpu_backend_correct_distances(tmp_path):
    paths = [
        _write_fasta(tmp_path, "A", "ACGTACGTACGT"),
        _write_fasta(tmp_path, "B", "ACGTAAGTACGT"),  # 1 SNP from A
        _write_fasta(tmp_path, "C", "TCGTACGTACGT"),  # 1 SNP from A
    ]
    snp = CPUBackend().compute(paths, n_jobs=1)
    assert snp.isolate_ids == ["A", "B", "C"]
    assert snp.distance("A", "B") == 1
    assert snp.distance("A", "C") == 1
    assert snp.distance("A", "A") == 0


def test_cpu_backend_handles_N(tmp_path):
    a = _write_fasta(tmp_path, "A", "ACGTACGT")
    b = _write_fasta(tmp_path, "B", "ACNTACGT")  # N at position 3 — masked
    snp = CPUBackend().compute([a, b], n_jobs=1)
    assert snp.distance("A", "B") == 0


def test_select_backend_default():
    backend = select_backend("auto")
    assert isinstance(backend, DistanceBackend)


def test_cache_round_trip(tmp_path):
    p = _write_fasta(tmp_path, "X", "AAAA")
    paths = [p]
    snp = SNPMatrix(isolate_ids=["X"], distances=np.zeros((1, 1)))
    save_cached(tmp_path, paths, snp)
    out = load_cached(tmp_path, paths)
    assert out is not None
    assert out.isolate_ids == ["X"]
