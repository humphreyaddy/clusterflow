"""Tests for clusterflow.ingestion."""

from __future__ import annotations

import pytest

from clusterflow.ingestion import IngestionError, IngestionPipeline, SNPMatrixReader


def test_loads_kp_fixtures(kp_config):
    result = IngestionPipeline().load(kp_config)
    assert len(result.isolates) == 28
    assert result.snp.n == 28
    # MLST coverage
    assert all(iso.sequence_type for iso in result.isolates.values())


def test_snp_matrix_reader_rejects_asymmetric(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("\t".join(["", "a", "b"]) + "\n")
    p.open("a").write("a\t0\t3\n")
    p.open("a").write("b\t5\t0\n")
    p.write_text("\t".join(["", "a", "b"]) + "\na\t0\t3\nb\t5\t0\n")
    with pytest.raises(IngestionError):
        SNPMatrixReader().read(p)


def test_snp_matrix_reader_rejects_missing(tmp_path):
    with pytest.raises(IngestionError):
        SNPMatrixReader().read(tmp_path / "nope.tsv")


def test_snp_matrix_reader_rejects_nonzero_diagonal(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("\t".join(["", "a", "b"]) + "\na\t1\t3\nb\t3\t0\n")
    with pytest.raises(IngestionError):
        SNPMatrixReader().read(p)
