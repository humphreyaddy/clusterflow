"""Phase 8.3 — performance benchmark across dataset sizes.

Generates synthetic outbreaks at increasing scale, runs the full pipeline,
and produces a runtime + memory report.
"""

from __future__ import annotations

import logging
import time
import tracemalloc
from datetime import date
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from clusterflow.clustering import run_all
from clusterflow.config import (
    ClusteringConfig,
    ClusterFlowConfig,
    DistanceEngineConfig,
    InputsConfig,
    StreamingConfig,
    ThresholdsConfig,
    VisualizationConfig,
)
from clusterflow.graph import GraphConstructor
from clusterflow.ingestion.epi import _parse_date  # noqa: F401 (used implicitly)
from clusterflow.models import Isolate, SNPMatrix
from clusterflow.testing.simulate import simulate_outbreak

log = logging.getLogger(__name__)


def _config_for_benchmark(out: Path) -> ClusterFlowConfig:
    """Minimal in-memory ClusterFlowConfig used for benchmark runs only."""
    return ClusterFlowConfig(
        project_name="benchmark",
        output_dir=out,
        inputs=InputsConfig(
            snp_matrix=out / "_dummy.tsv",
            mlst_profiles=None,
            epi_metadata=out / "_dummy.csv",
        ),
        thresholds=ThresholdsConfig(snp_cutoff=20, day_cutoff=30),
        clustering=ClusteringConfig(bootstrap_n=0, leiden_resolution=1.0),
        distance_engine=DistanceEngineConfig(backend="cpu", n_jobs=2),
        streaming=StreamingConfig(),
        visualization=VisualizationConfig(static=False, dashboard=False),
    )


def _bench_one(n: int, seed: int = 42) -> dict:
    sim = simulate_outbreak(
        n_isolates=n,
        n_clusters=max(2, min(10, n // 25)),
        random_state=seed,
    )
    snp_arr = sim["snp_matrix"].to_numpy(dtype=float)
    snp = SNPMatrix(
        isolate_ids=list(sim["snp_matrix"].index),
        distances=snp_arr,
    )
    isolates = {
        row["isolate_id"]: Isolate(
            isolate_id=row["isolate_id"],
            collection_date=date.fromisoformat(row["collection_date"]),
            facility=row["facility"],
            ward=row["ward"],
        )
        for _, row in sim["epi"].iterrows()
    }

    tracemalloc.start()
    timings: dict[str, float] = {}

    t0 = time.perf_counter()
    cfg = _config_for_benchmark(Path("/tmp/_bench"))
    g = GraphConstructor(cfg.thresholds).build(snp, isolates)
    timings["graph"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    run_all(g, cfg.thresholds, cfg.clustering)
    timings["clustering"] = time.perf_counter() - t0

    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return {
        "n_isolates": n,
        "n_edges": g.ecount(),
        "graph_seconds": timings["graph"],
        "clustering_seconds": timings["clustering"],
        "total_seconds": sum(timings.values()),
        "peak_memory_mb": peak / 1e6,
    }


def run_benchmark(
    sizes: list[int] | None = None,
    output_dir: str | Path = "results/benchmark",
) -> pd.DataFrame:
    sizes = sizes or [50, 200, 500, 1000]
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for n in sizes:
        log.info("benchmark n=%d", n)
        rows.append(_bench_one(n))
    df = pd.DataFrame(rows)
    df.to_csv(out / "performance_benchmark.csv", index=False)
    _plot_benchmark(df, out / "performance_plot.png")
    return df


def _plot_benchmark(df: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].plot(df["n_isolates"], df["graph_seconds"], "o-", label="graph build")
    axes[0].plot(df["n_isolates"], df["clustering_seconds"], "s-", label="clustering")
    axes[0].plot(df["n_isolates"], df["total_seconds"], "^-", label="total", lw=2)
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("n isolates")
    axes[0].set_ylabel("seconds")
    axes[0].set_title("runtime vs dataset size")
    axes[0].legend()
    axes[0].grid(True, which="both", alpha=0.3)

    axes[1].plot(df["n_isolates"], df["peak_memory_mb"], "o-", color="darkorange")
    axes[1].set_xscale("log")
    axes[1].set_xlabel("n isolates")
    axes[1].set_ylabel("peak memory (MB)")
    axes[1].set_title("peak memory vs dataset size")
    axes[1].grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
