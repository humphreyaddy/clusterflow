"""Tests for clusterflow.config."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from clusterflow.config import ClusterFlowConfig, ConfigError, diff_configs


def _write(tmp: Path, data: dict) -> Path:
    p = tmp / "config.yaml"
    p.write_text(yaml.safe_dump(data))
    return p


def test_loads_valid_config(fixtures_dir):
    cfg = ClusterFlowConfig.from_yaml(fixtures_dir / "kp_config.yaml")
    assert cfg.project_name == "kp_test_28"
    assert cfg.thresholds.snp_cutoff == 20


def test_missing_input_raises(tmp_path):
    p = _write(
        tmp_path,
        {
            "project_name": "x",
            "output_dir": str(tmp_path),
            "inputs": {
                "snp_matrix": str(tmp_path / "missing.tsv"),
                "epi_metadata": str(tmp_path / "missing.csv"),
            },
        },
    )
    with pytest.raises(ConfigError):
        ClusterFlowConfig.from_yaml(p)


def test_weights_must_sum_to_one(tmp_path, fixtures_dir):
    cfg = yaml.safe_load((fixtures_dir / "kp_config.yaml").read_text())
    cfg["thresholds"]["edge_weight_alpha"] = 0.9  # sum != 1
    p = tmp_path / "bad.yaml"
    p.write_text(yaml.safe_dump(cfg))
    with pytest.raises(Exception):
        ClusterFlowConfig.from_yaml(p)


def test_diff_configs(fixtures_dir, tmp_path):
    a = ClusterFlowConfig.from_yaml(fixtures_dir / "kp_config.yaml")
    b = a.model_copy(deep=True)
    b = b.model_copy(
        update={"thresholds": b.thresholds.model_copy(update={"snp_cutoff": 5})}
    )
    diff = diff_configs(a, b)
    assert "thresholds.snp_cutoff" in diff
    assert diff["thresholds.snp_cutoff"] == (20, 5)
