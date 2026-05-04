"""Synthetic K. pneumoniae outbreak simulator.

Used for test fixtures and Phase 8 performance benchmarks. Produces a
realistic SNP distance matrix + epi metadata + MLST file with a tunable
number of clusters and isolates.
"""

from __future__ import annotations

from datetime import date, timedelta
from pathlib import Path

import numpy as np
import pandas as pd


WARDS = ["NICU_A", "NICU_B", "MAT", "PED", "ICU"]
FACILITIES = ["Hospital_1", "Hospital_2"]


def simulate_outbreak(
    n_isolates: int = 28,
    n_clusters: int = 4,
    snp_within_cluster: int = 8,
    snp_between_cluster: int = 80,
    base_date: date = date(2024, 1, 15),
    span_days: int = 90,
    random_state: int = 42,
) -> dict:
    """Return a dict with keys ``snp_matrix``, ``epi``, ``mlst``.

    ``snp_matrix`` is an n×n numpy array. ``epi`` and ``mlst`` are pandas
    DataFrames matching ClusterFlow's expected input format.
    """
    rng = np.random.default_rng(random_state)

    # Assign isolates to clusters (roughly balanced)
    cluster_assignments = np.repeat(
        np.arange(n_clusters), n_isolates // n_clusters
    ).tolist()
    while len(cluster_assignments) < n_isolates:
        cluster_assignments.append(rng.integers(0, n_clusters))
    rng.shuffle(cluster_assignments)
    cluster_assignments = np.asarray(cluster_assignments[:n_isolates])

    # Distance matrix: small within-cluster, large between
    D = np.zeros((n_isolates, n_isolates), dtype=float)
    for i in range(n_isolates):
        for j in range(i + 1, n_isolates):
            if cluster_assignments[i] == cluster_assignments[j]:
                d = max(0, rng.integers(0, snp_within_cluster + 1))
            else:
                d = rng.integers(snp_between_cluster - 20, snp_between_cluster + 80)
            D[i, j] = D[j, i] = d

    # Cluster-level collection date windows so within-cluster days are tight
    cluster_starts = {
        cid: rng.integers(0, max(span_days - 30, 1)) for cid in range(n_clusters)
    }

    isolate_ids: list[str] = []
    epi_rows: list[dict] = []
    mlst_rows: list[dict] = []
    cluster_to_st = {cid: f"{37 + cid}" for cid in range(n_clusters)}
    cluster_to_ward = {
        cid: WARDS[cid % len(WARDS)] for cid in range(n_clusters)
    }

    for i in range(n_isolates):
        iso = f"KP{2024}_{i + 1:03d}"
        isolate_ids.append(iso)
        cid = int(cluster_assignments[i])
        d = base_date + timedelta(
            days=int(cluster_starts[cid]) + int(rng.integers(0, 21))
        )
        # ~10% noise: occasional ward/STs swap to test robustness
        ward = (
            cluster_to_ward[cid]
            if rng.random() > 0.1
            else WARDS[rng.integers(0, len(WARDS))]
        )
        facility = FACILITIES[rng.integers(0, len(FACILITIES))]
        epi_rows.append(
            {
                "isolate_id": iso,
                "collection_date": d.isoformat(),
                "facility": facility,
                "ward": ward,
            }
        )
        st = (
            cluster_to_st[cid]
            if rng.random() > 0.05
            else str(rng.integers(40, 60))
        )
        mlst_rows.append({"isolate_id": iso, "ST": f"ST{st}"})

    df_snp = pd.DataFrame(D, index=isolate_ids, columns=isolate_ids).astype(int)
    df_epi = pd.DataFrame(epi_rows)
    df_mlst = pd.DataFrame(mlst_rows)
    return {
        "snp_matrix": df_snp,
        "epi": df_epi,
        "mlst": df_mlst,
        "true_clusters": dict(zip(isolate_ids, cluster_assignments.tolist())),
    }


def write_fixtures(
    output_dir: str | Path,
    n_isolates: int = 28,
    n_clusters: int = 4,
    random_state: int = 42,
) -> dict[str, Path]:
    """Write synthetic dataset to ``output_dir`` in ClusterFlow's expected format."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    sim = simulate_outbreak(
        n_isolates=n_isolates,
        n_clusters=n_clusters,
        random_state=random_state,
    )
    snp_path = out / "kp_snp_matrix.tsv"
    epi_path = out / "kp_epi.csv"
    mlst_path = out / "kp_mlst.csv"
    truth_path = out / "kp_truth.csv"

    sim["snp_matrix"].to_csv(snp_path, sep="\t")
    sim["epi"].to_csv(epi_path, index=False)
    sim["mlst"].to_csv(mlst_path, index=False)
    pd.DataFrame(
        [{"isolate_id": k, "true_cluster": v} for k, v in sim["true_clusters"].items()]
    ).to_csv(truth_path, index=False)

    return {
        "snp": snp_path,
        "epi": epi_path,
        "mlst": mlst_path,
        "truth": truth_path,
    }
