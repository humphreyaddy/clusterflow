"""Tests for clusterflow.models."""

from __future__ import annotations

from datetime import date

import numpy as np
import pytest
from pydantic import ValidationError

from clusterflow.models import ClusterAssignment, Isolate, SNPMatrix


def test_isolate_immutable():
    iso = Isolate(
        isolate_id="A",
        collection_date=date(2024, 1, 1),
        facility="H",
        ward="W",
    )
    with pytest.raises(ValidationError):
        iso.isolate_id = "B"  # type: ignore[misc]


def test_snp_matrix_validates_symmetry():
    bad = np.array([[0.0, 1.0], [2.0, 0.0]])
    with pytest.raises(ValidationError):
        SNPMatrix(isolate_ids=["a", "b"], distances=bad)


def test_snp_matrix_validates_diagonal():
    bad = np.array([[1.0, 0.0], [0.0, 1.0]])
    with pytest.raises(ValidationError):
        SNPMatrix(isolate_ids=["a", "b"], distances=bad)


def test_snp_matrix_negative():
    bad = np.array([[0.0, -1.0], [-1.0, 0.0]])
    with pytest.raises(ValidationError):
        SNPMatrix(isolate_ids=["a", "b"], distances=bad)


def test_snp_matrix_distance_lookup():
    arr = np.array([[0.0, 5.0], [5.0, 0.0]])
    m = SNPMatrix(isolate_ids=["a", "b"], distances=arr)
    assert m.distance("a", "b") == 5.0
    assert m.n == 2


def test_cluster_assignment_count_consistency():
    with pytest.raises(ValidationError):
        ClusterAssignment(
            method="x",
            assignments={"a": 0, "b": 1},
            n_clusters=5,
        )
