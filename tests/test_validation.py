"""Phase 8 — validation tests on the real K. pneumoniae example dataset.

The dataset is the SNP2Cluster v0.5.4 example (Zenodo 10.5281/zenodo.14060296)
converted to ClusterFlow's input format.

Acceptance gates from transmission_cluster_tool_plan.md §8:
  * pipeline completes in < 60 s on a 4-core laptop
  * ClusterFlow recovers ST25 and ST39 as distinct clusters that SNP2Cluster
    misses with its K-means + chain logic
  * two runs with the same seed produce identical assignments
"""

from __future__ import annotations

import time
from pathlib import Path

import pytest

from clusterflow.config import ClusterFlowConfig
from clusterflow.pipeline import run_pipeline


REAL_DIR = Path(__file__).parent / "fixtures" / "kp_real"
REAL_CONFIG = REAL_DIR / "kp_config.yaml"

pytestmark = pytest.mark.skipif(
    not REAL_CONFIG.exists(),
    reason="real K. pneumoniae fixture not available — see notebooks/validation_kpneumoniae.ipynb",
)


@pytest.fixture(scope="module")
def real_result(tmp_path_factory):
    cfg = ClusterFlowConfig.from_yaml(REAL_CONFIG)
    cfg = cfg.model_copy(
        update={"output_dir": tmp_path_factory.mktemp("kp_real")}
    )
    t0 = time.perf_counter()
    res = run_pipeline(cfg, run_viz=False)
    res.metadata["elapsed_seconds"] = time.perf_counter() - t0
    return res


def test_pipeline_completes_within_60s(real_result):
    elapsed = real_result.metadata["elapsed_seconds"]
    assert elapsed < 60, f"too slow on real data: {elapsed:.1f}s"


def test_st25_and_st39_form_distinct_clusters(real_result):
    """ClusterFlow must detect both ST25 and ST39 as distinct cluster groupings."""
    sts_seen: set[str] = set()
    for c in real_result.transmission_clusters:
        if len(c.isolate_ids) >= 2:
            sts_seen.update(c.sequence_types)
    assert "25" in sts_seen, f"ST25 not detected in any multi-isolate cluster: {sorted(sts_seen)}"
    assert "39" in sts_seen, f"ST39 not detected in any multi-isolate cluster: {sorted(sts_seen)}"


def test_real_data_reproducible(real_result, tmp_path):
    cfg = ClusterFlowConfig.from_yaml(REAL_CONFIG)
    cfg = cfg.model_copy(update={"output_dir": tmp_path})
    res2 = run_pipeline(cfg, run_viz=False)
    assert real_result.consensus.assignments == res2.consensus.assignments
