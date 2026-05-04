"""Shared pytest fixtures."""

from __future__ import annotations

from pathlib import Path

import pytest

from clusterflow.config import ClusterFlowConfig
from clusterflow.ingestion import IngestionPipeline


FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture(scope="session")
def fixtures_dir() -> Path:
    return FIXTURES


@pytest.fixture(scope="session")
def kp_config(tmp_path_factory: pytest.TempPathFactory) -> ClusterFlowConfig:
    """Load the bundled KP test config and redirect output to tmp dir."""
    cfg = ClusterFlowConfig.from_yaml(FIXTURES / "kp_config.yaml")
    out = tmp_path_factory.mktemp("kp_results")
    return cfg.model_copy(update={"output_dir": out})


@pytest.fixture(scope="session")
def kp_ingestion(kp_config):
    return IngestionPipeline().load(kp_config)


@pytest.fixture(scope="session")
def kp_truth(fixtures_dir) -> dict[str, int]:
    import pandas as pd

    df = pd.read_csv(fixtures_dir / "kp_truth.csv")
    return dict(zip(df["isolate_id"], df["true_cluster"]))
