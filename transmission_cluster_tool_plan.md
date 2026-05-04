# Transmission Cluster Detection Tool — Phased Project Plan

> **Codename:** `ClusterFlow`
> **Stack:** Python 3.11+ (primary) · R bridge (optional, for SNP2Cluster compatibility)
> **Validation dataset:** *Klebsiella pneumoniae* neonatal unit outbreak (Kwenda et al. 2024)
> **Reference baseline:** SNP2Cluster v0.5.4 (https://github.com/stanikae/SNP2Cluster/tree/v0.5.4)

---

## Project Overview

ClusterFlow is a modular, scalable transmission cluster detection pipeline that fuses whole-genome sequencing (WGS) data with epidemiological metadata using graph-theoretic algorithms and parallel computing. It supersedes the K-means + SNP-chain approach of SNP2Cluster by introducing Leiden community detection, spectral clustering, temporal directed acyclic graphs (DAGs), centrality-based index case identification, and a real-time streaming mode with a live surveillance dashboard.

### High-Level Architecture

```
Inputs (SNP matrix · MLST · Epi metadata)
        ↓
[Phase 2] Parallel Distance Engine  (GPU/CPU pairwise SNP)
        ↓
[Phase 3] Graph Constructor         (weighted temporal graph)
        ↓ ↓ ↓ (parallel)
[Phase 4] Cluster Detection         (SNP chains · Leiden · Spectral + Consensus)
        ↓ ↓ ↓ (parallel)
[Phase 5] Analysis Layer            (Temporal DAG · Centrality · Bootstrap CIs)
        ↓
[Phase 6] Visualization Suite       (Heatmap · MST · Transmission Tree · Dashboard)
        ↓
[Phase 7] Streaming Mode            (Incremental graph · FastAPI endpoint)
        ↓
[Phase 8] Validation                (K. pneumoniae benchmark vs SNP2Cluster)
        ↓
[Phase 9] Packaging                 (CLI · Docker · Docs · PyPI)
```

---

## Phase 1 — Foundation

**Goal:** Establish the repository structure, configuration system, shared data models, and ingestion layer. Everything downstream depends on this being solid.

### 1.1 Repository Structure

**Input:** None (greenfield).

**Steps:**
1. Initialise Git repository with the following layout:
```
clusterflow/
├── clusterflow/
│   ├── __init__.py
│   ├── config.py          # YAML config loader + validator
│   ├── models.py          # Pydantic data models
│   ├── ingestion/
│   │   ├── __init__.py
│   │   ├── snp_matrix.py  # SNP distance matrix reader
│   │   ├── mlst.py        # MLST profile reader
│   │   └── epi.py         # Epidemiological metadata reader
│   ├── distance/          # Phase 2
│   ├── graph/             # Phase 3
│   ├── clustering/        # Phase 4
│   ├── analysis/          # Phase 5
│   ├── viz/               # Phase 6
│   ├── streaming/         # Phase 7
│   └── cli.py             # Phase 9
├── tests/
│   ├── fixtures/          # Synthetic + real validation datasets
│   └── test_*.py
├── configs/
│   └── example.yaml
├── notebooks/
│   └── validation_kpneumoniae.ipynb
├── docs/
├── Dockerfile
├── pyproject.toml
└── README.md
```

2. Set up `pyproject.toml` with dependency groups: `core`, `gpu` (optional), `dashboard`, `dev`.
3. Configure `pre-commit` hooks: `ruff` (linting), `black` (formatting), `mypy` (type checking).

**Output:**
- Initialised Git repository with the directory structure above.
- `pyproject.toml` with pinned core dependencies (see Section 1.4).
- `Makefile` with targets: `install`, `test`, `lint`, `docs`, `docker-build`.

---

### 1.2 Configuration System

**Input:** User-provided YAML config file (path passed via CLI or API).

**Steps:**
1. Define `ClusterFlowConfig` as a Pydantic `BaseSettings` model.
2. Support the following top-level keys (designed to be backward-compatible with SNP2Cluster's config format):

```yaml
# clusterflow example config
project_name: "KP_neonatal_2024"
output_dir: "./results"

inputs:
  snp_matrix: "data/snp_distances.tsv"   # TSV/CSV pairwise distance matrix
  mlst_profiles: "data/mlst.csv"          # CSV with columns: isolate_id, ST, allele_1..N
  epi_metadata: "data/epi.csv"            # CSV with columns: isolate_id, collection_date, facility, ward

thresholds:
  snp_cutoff: 20                          # Max SNP distance for an edge (default: 20)
  day_cutoff: 14                          # Max collection date gap in days (default: 14)
  edge_weight_alpha: 0.4                  # Weight for SNP distance term
  edge_weight_beta: 0.4                   # Weight for temporal term
  edge_weight_gamma: 0.2                  # Weight for MLST mismatch term

clustering:
  methods: ["snp_chains", "leiden", "spectral"]  # Run all three by default
  bootstrap_n: 500                        # Bootstrap replicates for CI
  leiden_resolution: 1.0                  # Leiden resolution parameter

distance_engine:
  backend: "auto"     # "auto" | "cpu" | "gpu" | "pairsnp"
  n_jobs: -1          # -1 = use all available cores

streaming:
  enabled: false
  api_port: 8000

visualization:
  static: true        # Generate heatmap, MST, scatter plots
  dashboard: false    # Launch live Dash dashboard
  dashboard_port: 8050
```

3. Implement `ConfigValidator` that checks file paths exist, thresholds are in valid ranges, and column names match expected schemas.
4. Provide a `config diff` CLI command to compare two configs and highlight changed parameters.

**Output:**
- `clusterflow/config.py` — `ClusterFlowConfig` Pydantic model with field validators.
- `configs/example.yaml` — fully annotated example config.
- `configs/snp2cluster_compat.yaml` — minimal config that mirrors SNP2Cluster's defaults exactly.
- Unit tests: `tests/test_config.py` — valid config loads, missing files raise `ConfigError`, invalid thresholds raise `ValidationError`.

---

### 1.3 Shared Data Models

**Input:** Parsed config.

**Steps:**
1. Define the following Pydantic models in `clusterflow/models.py`:

```python
class Isolate:
    isolate_id: str
    collection_date: date
    facility: str
    ward: str
    sequence_type: Optional[str]          # from MLST
    mlst_alleles: Optional[Dict[str, int]]

class SNPMatrix:
    isolate_ids: List[str]
    distances: np.ndarray                 # shape (n, n), symmetric, diagonal = 0

class WeightedEdge:
    source: str
    target: str
    snp_distance: int
    day_delta: int
    mlst_mismatch: bool
    composite_weight: float               # computed from alpha/beta/gamma

class ClusterAssignment:
    method: str                           # "snp_chains" | "leiden" | "spectral"
    assignments: Dict[str, int]           # isolate_id → cluster_id
    n_clusters: int
    stability_scores: Optional[Dict[str, float]]  # from bootstrap

class TransmissionCluster:
    cluster_id: int
    isolate_ids: List[str]
    sequence_types: List[str]
    wards: List[str]
    date_range: Tuple[date, date]
    index_case_candidate: Optional[str]   # highest betweenness isolate
    consensus_method: str                 # which method(s) agree
    confidence: float                     # bootstrap stability score
```

2. Define `PipelineResult` as the top-level output object carrying all cluster assignments, analysis results, and metadata.

**Output:**
- `clusterflow/models.py` — all Pydantic models with field descriptions.
- Unit tests: `tests/test_models.py` — model instantiation, serialisation to JSON, deserialisation.

---

### 1.4 Ingestion Layer

**Input:**
- SNP distance matrix: TSV or CSV file, symmetric, with isolate IDs as row/column headers. First column = isolate IDs.
- MLST profiles: CSV with mandatory columns `isolate_id`, `ST` and optional allele columns.
- Epi metadata: CSV with mandatory columns `isolate_id`, `collection_date` (ISO 8601), `facility`, `ward`.

**Steps:**
1. Implement `SNPMatrixReader`:
   - Accept TSV, CSV, or `snp-dists`-format output.
   - Validate: matrix is square, symmetric (within tolerance 1e-6), diagonal is zero, no negative values.
   - Raise `IngestionError` with row/column position on first violation.
   - Return `SNPMatrix` model.

2. Implement `MLSTReader`:
   - Accept CSV or TSV.
   - Normalise ST field: strip "ST" prefix if present, handle novel allele combinations ("novel", "–").
   - Return `Dict[str, str]` (isolate_id → ST string).

3. Implement `EpiMetadataReader`:
   - Accept CSV or TSV.
   - Parse `collection_date` with multiple format support (ISO 8601, DD/MM/YYYY, DD-MM-YYYY).
   - Validate: no future dates, date range is plausible (warn if span > 5 years).
   - Return `Dict[str, Isolate]`.

4. Implement `IngestionPipeline.load(config)`:
   - Calls all three readers.
   - Cross-validates: all isolate IDs in SNP matrix must appear in epi metadata; warn if MLST coverage is partial.
   - Logs summary: n isolates, date range, wards, STs detected, % missing data per field.

**Output:**
- `clusterflow/ingestion/snp_matrix.py`, `mlst.py`, `epi.py`
- `clusterflow/ingestion/__init__.py` exposing `IngestionPipeline`
- Unit tests: `tests/test_ingestion.py` — valid files load, malformed files raise descriptive errors, cross-validation catches ID mismatches.
- Test fixtures: `tests/fixtures/kp_snp_matrix.tsv`, `kp_mlst.csv`, `kp_epi.csv` (synthetic K. pneumoniae data matching the Kwenda et al. outbreak profile).

---

### Phase 1 Deliverables Summary

| Deliverable | Path | Description |
|---|---|---|
| Repo scaffold | `/` | Full directory structure, `pyproject.toml`, `Makefile` |
| Config system | `clusterflow/config.py` | Pydantic config model + YAML loader |
| Example configs | `configs/` | Annotated YAML + SNP2Cluster-compat YAML |
| Data models | `clusterflow/models.py` | All shared Pydantic models |
| Ingestion layer | `clusterflow/ingestion/` | Three readers + `IngestionPipeline` |
| Test fixtures | `tests/fixtures/` | Synthetic K. pneumoniae dataset (28 isolates) |
| Test suite (Phase 1) | `tests/test_config.py`, `test_models.py`, `test_ingestion.py` | ≥90% coverage of Phase 1 code |

---

## Phase 2 — Parallel Distance Engine

**Goal:** Compute pairwise SNP distances efficiently across all dataset scales. At small scale (≤500 isolates), use Python multiprocessing. At medium scale (≤5,000), use chunked parallel workers. At large scale (5,000+), delegate to `pairsnp` (C++ binary) or `cuGraph` (GPU). The engine should be backend-agnostic behind a single function signature.

### 2.1 Backend Interface

**Input:** List of FASTA file paths (raw sequences) OR pre-computed `SNPMatrix` (skip computation).

**Steps:**
1. Define abstract base class `DistanceBackend` with a single method:
```python
def compute(self, fasta_paths: List[Path]) -> SNPMatrix: ...
```
2. Implement `AutoBackend` that selects the best available backend at runtime:
   - If `cuGraph`/RAPIDS is importable and a CUDA device is available → `GPUBackend`
   - Else if `pairsnp` binary is on PATH → `PairSNPBackend`
   - Else → `CPUBackend`
3. Log which backend was selected and estimated runtime.

**Output:** `clusterflow/distance/backend.py` — abstract interface + `AutoBackend`.

---

### 2.2 CPU Backend

**Input:** List of FASTA file paths OR in-memory sequence strings.

**Steps:**
1. Load sequences into memory as numpy byte arrays (ASCII encoding of ACGT).
2. Partition the upper triangle of the n×n distance matrix into `n_jobs` chunks using `numpy` index arithmetic.
3. Dispatch chunks to `joblib.Parallel(backend="loky")` workers.
4. Each worker computes Hamming distances for its chunk, optionally masking ambiguous bases (N, -, gaps) based on a config flag.
5. Collect results into the full symmetric matrix.
6. Emit a progress bar via `tqdm` (suppressible via `--quiet` flag).

**Input:** FASTA sequences (list of `(id, sequence)` tuples).
**Output:** `SNPMatrix` with distances in SNP count (integer).

---

### 2.3 PairSNP Backend

**Input:** List of FASTA file paths, path to `pairsnp` binary (auto-detected from PATH).

**Steps:**
1. Verify `pairsnp` version ≥ 0.3.0 via `subprocess.run(["pairsnp", "--version"])`.
2. Build command: `pairsnp -s {snp_threshold * 2} -t {n_jobs} {fasta_path}`.
3. Stream stdout (TSV format) directly into the `SNPMatrix` builder to avoid holding the full matrix in memory for large datasets.
4. Handle non-zero exit codes with informative errors.

**Output:** `SNPMatrix`.

---

### 2.4 GPU Backend (optional)

**Input:** FASTA sequences as numpy byte arrays.

**Steps:**
1. Guard import with `try/except ImportError` — GPU backend is optional; missing RAPIDS should not break CPU path.
2. Transfer sequence matrix to GPU memory via `cupy.array()`.
3. Compute pairwise Hamming distances using `cupy` broadcasting: `(A[:, None, :] != A[None, :, :]).sum(axis=2)`.
4. Apply ambiguous base mask on GPU.
5. Transfer result back to CPU as numpy array.
6. Wrap in `SNPMatrix`.

**Output:** `SNPMatrix`.

---

### 2.5 Matrix Cache

**Steps:**
1. After computation, serialise `SNPMatrix` to a compressed `.npz` file in the output directory.
2. On subsequent runs with the same input files (detected by SHA-256 hash of all input FASTAs), load the cached matrix and skip recomputation. Log cache hit.
3. Provide `--no-cache` flag to force recomputation.

**Output:** `{output_dir}/cache/snp_matrix_{hash}.npz`

---

### Phase 2 Deliverables Summary

| Deliverable | Path | Description |
|---|---|---|
| Backend interface | `clusterflow/distance/backend.py` | Abstract class + `AutoBackend` |
| CPU backend | `clusterflow/distance/cpu.py` | `joblib`-parallel Hamming distance |
| PairSNP backend | `clusterflow/distance/pairsnp.py` | Subprocess wrapper for `pairsnp` |
| GPU backend | `clusterflow/distance/gpu.py` | Optional `cupy`-based backend |
| Matrix cache | `clusterflow/distance/cache.py` | SHA-256-keyed `.npz` cache |
| Test suite (Phase 2) | `tests/test_distance.py` | CPU/PairSNP backends; result parity check between backends; cache hit/miss |

---

## Phase 3 — Graph Constructor

**Goal:** Build a weighted, optionally directed `igraph.Graph` from the `SNPMatrix`, `Isolate` metadata, and MLST profiles. The graph is the core data structure for all downstream analyses.

### 3.1 Edge Construction

**Input:** `SNPMatrix`, `Dict[str, Isolate]`, config thresholds (`snp_cutoff`, `day_cutoff`, `alpha`, `beta`, `gamma`).

**Steps:**
1. Iterate over all pairs `(i, j)` where `i < j` (upper triangle).
2. For each pair, compute the composite edge weight:
```python
snp_term   = snp_distance[i,j] / snp_cutoff       # normalised, 0..1+
time_term  = abs((date_i - date_j).days) / day_cutoff
mlst_term  = 0 if ST_i == ST_j else 1

weight = alpha * snp_term + beta * time_term + gamma * mlst_term
```
3. Include edge if `snp_distance[i,j] <= snp_cutoff` AND `time_term <= 1.0`.
4. Store on each edge: `snp_distance`, `day_delta`, `mlst_mismatch` (bool), `composite_weight`.
5. Store on each vertex: all `Isolate` fields as vertex attributes.

**Output:** `igraph.Graph` (undirected at this stage) with vertex and edge attributes.

---

### 3.2 Graph Serialisation

**Steps:**
1. Serialise graph to `GraphML` format (lossless, human-readable, importable into Gephi/Cytoscape).
2. Also serialise to `pickle` for fast re-loading within the pipeline.
3. Write a summary JSON: `{n_vertices, n_edges, density, n_components, largest_component_size, mean_degree, edge_weight_stats}`.

**Output:**
- `{output_dir}/graph/transmission_graph.graphml`
- `{output_dir}/graph/transmission_graph.pkl`
- `{output_dir}/graph/graph_summary.json`

---

### 3.3 Incremental Graph Updates (stub for Phase 7)

**Steps:**
1. Implement `GraphConstructor.add_isolate(isolate, snp_distances_to_existing)` method.
2. Only recompute edges for the new isolate; existing edges are unchanged.
3. Return the updated graph and a list of newly added edges.

**Output:** Updated `igraph.Graph` in-place; list of `WeightedEdge` objects for new connections.

---

### Phase 3 Deliverables Summary

| Deliverable | Path | Description |
|---|---|---|
| Graph constructor | `clusterflow/graph/constructor.py` | Edge construction + composite weight computation |
| Graph serialisation | `clusterflow/graph/io.py` | GraphML + pickle serialisation + summary JSON |
| Incremental stub | `clusterflow/graph/incremental.py` | `add_isolate` method (used in Phase 7) |
| Test suite (Phase 3) | `tests/test_graph.py` | Graph built from fixtures; edge count/weight validation; serialisation round-trip |

---

## Phase 4 — Cluster Detection

**Goal:** Run three clustering algorithms in parallel on the same graph, then compute a consensus assignment. Produces four cluster assignment objects: one per algorithm plus the consensus.

### 4.1 Algorithm 1 — SNP Chain Clusters

**Input:** `igraph.Graph`, `snp_cutoff` threshold.

**Description:** A faithful re-implementation of SNP2Cluster's chain logic. Provides backward compatibility and serves as the baseline for benchmarking.

**Steps:**
1. Sort isolates by collection date.
2. Apply connected components on a subgraph where all edges have `snp_distance <= snp_cutoff`.
3. Within each component, apply SNP chain logic: walk isolates in date order; assign to current cluster if `snp_distance` to any member ≤ cutoff AND date gap ≤ `day_cutoff`; otherwise start a new cluster.
4. Each isolate is assigned to exactly one cluster (first qualifying assignment wins).

**Input:** `igraph.Graph` (with `snp_distance` edge attribute), sorted by `collection_date`.
**Output:** `ClusterAssignment(method="snp_chains", assignments={...}, n_clusters=N)`.

---

### 4.2 Algorithm 2 — Leiden Community Detection

**Input:** `igraph.Graph` with `composite_weight` edge attribute.

**Steps:**
1. Use `leidenalg` library with `ModularityVertexPartition`.
2. Set edge weights to `composite_weight` (lower weight = closer = stronger community pull; invert if needed: `1 / composite_weight`).
3. Run optimisation with `n_iterations=10` (configurable).
4. Run `leiden_resolution` sweep from 0.5 to 2.0 in 0.1 steps using `joblib.Parallel`; select resolution that maximises modularity Q.
5. Return partition with optimal Q.

**Input:** `igraph.Graph`, `leiden_resolution` (float or `"auto"`).
**Output:** `ClusterAssignment(method="leiden", assignments={...}, n_clusters=N)`.

---

### 4.3 Algorithm 3 — Spectral Clustering

**Input:** `igraph.Graph` with `composite_weight` edge attribute.

**Steps:**
1. Build the graph Laplacian matrix L from the adjacency matrix weighted by `composite_weight`.
2. Compute the `k` smallest eigenvectors of L, where `k` is determined by the eigengap heuristic (largest gap in sorted eigenvalues).
3. Project isolates into the `k`-dimensional eigenspace.
4. Run `sklearn.cluster.KMeans(n_clusters=k)` in this eigenspace.
5. Use `sklearn.metrics.silhouette_score` to validate `k`; try `k-1` and `k+1` and select best silhouette score.

**Input:** `igraph.Graph`, optional `k` override (default: `"auto"`).
**Output:** `ClusterAssignment(method="spectral", assignments={...}, n_clusters=N)`.

---

### 4.4 Parallel Execution

**Steps:**
1. Wrap all three algorithm calls in `joblib.Parallel(n_jobs=3, backend="threading")`.
2. Each algorithm receives its own copy of the graph (thread-safe `igraph.Graph.copy()`).
3. Collect results into `List[ClusterAssignment]`.

---

### 4.5 Consensus Assignment

**Input:** `List[ClusterAssignment]` (one per method).

**Steps:**
1. For each isolate, record the cluster assignment from each method.
2. Map cluster IDs to a canonical form using the Hungarian algorithm (solve the label-alignment problem across methods).
3. For each isolate, compute agreement score: fraction of methods that agree on the same (aligned) cluster.
4. Assign the majority-vote cluster as the consensus assignment.
5. Flag isolates with agreement < 2/3 as `ambiguous=True` in the output.
6. Log consensus statistics: N isolates fully agreed, N ambiguous, per-cluster stability.

**Output:** `ClusterAssignment(method="consensus", assignments={...}, n_clusters=N)` with an additional `ambiguous` flag and `agreement_score` per isolate.

---

### Phase 4 Deliverables Summary

| Deliverable | Path | Description |
|---|---|---|
| SNP chain algorithm | `clusterflow/clustering/snp_chains.py` | SNP2Cluster-compatible chain logic |
| Leiden algorithm | `clusterflow/clustering/leiden.py` | `leidenalg`-based community detection |
| Spectral algorithm | `clusterflow/clustering/spectral.py` | Eigengap-based spectral clustering |
| Parallel runner | `clusterflow/clustering/runner.py` | `joblib.Parallel` wrapper for all three |
| Consensus engine | `clusterflow/clustering/consensus.py` | Hungarian alignment + majority vote |
| Test suite (Phase 4) | `tests/test_clustering.py` | All three algorithms on K. pneumoniae fixtures; expected cluster count = 4–6; consensus agrees with SNP2Cluster on ≥80% of isolates |

---

## Phase 5 — Analysis Layer

**Goal:** Produce three analytical outputs that go beyond what SNP2Cluster can generate: (1) a directed transmission tree, (2) centrality-based index case identification, and (3) bootstrap confidence intervals for cluster assignments.

### 5.1 Temporal DAG Construction

**Input:** `igraph.Graph` (undirected), `Dict[str, Isolate]` (for collection dates), `day_cutoff`, temporal uncertainty window (default: 2 days — accounts for imprecise collection dates).

**Steps:**
1. For each edge `(A, B)` in the undirected graph:
   - If `date_A < date_B - uncertainty_window` → directed edge A → B (A precedes B).
   - If `date_B < date_A - uncertainty_window` → directed edge B → A.
   - If dates are within `±uncertainty_window` → bidirectional edge (temporal order ambiguous); retain both directions.
2. Remove any directed cycles (should be rare; log if found — indicates data quality issue).
3. Compute topological ordering of the DAG.
4. Identify source nodes (in-degree = 0 within each cluster): these are index case candidates.

**Output:**
- `igraph.Graph` (directed) — temporal transmission DAG.
- `{output_dir}/analysis/transmission_dag.graphml`
- `{output_dir}/analysis/dag_summary.json` with `{n_source_nodes, n_sink_nodes, max_chain_length, clusters_with_ambiguous_order}`

---

### 5.2 Centrality Scoring

**Input:** Temporal DAG from Step 5.1, consensus `ClusterAssignment`.

**Steps:**
1. Compute per-isolate metrics using `igraph`:
   - `betweenness_centrality(directed=True, weights="composite_weight")` — isolates on many shortest transmission paths.
   - `in_degree()` — number of probable upstream sources.
   - `out_degree()` — number of probable downstream transmissions.
   - `closeness_centrality()` — speed of potential spread from this node.
2. Normalise all scores to [0, 1] range within each cluster.
3. Compute composite `transmission_risk_score = 0.5*betweenness + 0.3*out_degree_norm + 0.2*closeness`.
4. For each cluster, rank isolates by `transmission_risk_score` and annotate the top-ranked isolate as `index_case_candidate`.
5. Cross-check: if `index_case_candidate` is not among the earliest-collected isolates in that cluster (by date rank), emit a `WARNING: temporal mismatch — early date and high centrality disagree` log entry for human review.

**Output:**
- `{output_dir}/analysis/centrality_scores.csv` — one row per isolate with all centrality metrics.
- `{output_dir}/analysis/index_case_candidates.csv` — one row per cluster with top candidate and score.

---

### 5.3 Bootstrap Confidence Intervals

**Input:** `SNPMatrix`, config `bootstrap_n` (default: 500), `ClusterFlowConfig`.

**Steps:**
1. For each bootstrap replicate `b` in `range(bootstrap_n)`:
   a. Sample isolate pairs with replacement from the upper triangle of the SNP matrix (preserving symmetry).
   b. Re-build the graph with perturbed edge weights (SNP distances ±uniform[0, 2] SNPs to model sequencing uncertainty).
   c. Re-run Leiden clustering (fastest of the three methods) on the perturbed graph.
   d. Store cluster assignment for each isolate.
2. Parallelise across `n_jobs` cores using `joblib.Parallel(backend="loky")`.
3. For each isolate, compute `stability_score` = fraction of bootstrap replicates in which it was assigned to the same (label-aligned) cluster as in the consensus run.
4. Classify isolates: `stable` (score ≥ 0.90), `uncertain` (0.70–0.90), `unstable` (< 0.70).
5. For each cluster, compute cluster-level stability as the mean stability score of its members.

**Output:**
- `{output_dir}/analysis/bootstrap_stability.csv` — per-isolate stability scores and classification.
- `{output_dir}/analysis/cluster_stability_summary.csv` — per-cluster mean stability and confidence grade (A/B/C).

---

### Phase 5 Deliverables Summary

| Deliverable | Path | Description |
|---|---|---|
| Temporal DAG | `clusterflow/analysis/dag.py` | Directed graph construction + cycle detection |
| Centrality scoring | `clusterflow/analysis/centrality.py` | All centrality metrics + index case ranking |
| Bootstrap CI | `clusterflow/analysis/bootstrap.py` | Parallel bootstrap + stability classification |
| Test suite (Phase 5) | `tests/test_analysis.py` | DAG has no cycles; K. pneumoniae known index cases recovered; bootstrap stability reproducible across seeds |

---

## Phase 6 — Visualization Suite

**Goal:** Generate all outputs in two tiers: (A) publication-ready static figures matching and extending SNP2Cluster's output, and (B) an interactive live dashboard for operational use.

### 6.1 Static Outputs (Tier A)

All static figures are generated using `matplotlib`/`seaborn` and saved as both PNG (300 DPI) and SVG (for publication editing).

#### Figure 1 — SNP Distance Heatmap

**Input:** `SNPMatrix`, consensus `ClusterAssignment`.

**Steps:**
1. Reorder isolates by cluster, then by collection date within each cluster.
2. Plot heatmap using `seaborn.clustermap` with the reordered distance matrix.
3. Add colour-coded row/column annotations for: cluster ID, ward, sequence type.
4. Draw cluster boundary boxes.
5. Add colour bar for SNP distance scale.

**Output:** `{output_dir}/figures/snp_heatmap.png`, `.svg`

---

#### Figure 2 — Minimum Spanning Tree

**Input:** `igraph.Graph`, consensus `ClusterAssignment`, centrality scores.

**Steps:**
1. Compute MST using `igraph.Graph.spanning_tree(weights="composite_weight")`.
2. Layout using Kamada-Kawai algorithm.
3. Node colour = cluster. Node size = `transmission_risk_score`. Node shape = ward.
4. Edge thickness = 1 / snp_distance (closer isolates → thicker edges).
5. Annotate index case candidates with a star marker.
6. Add legend for cluster colours, ward shapes, and centrality scale.

**Output:** `{output_dir}/figures/minimum_spanning_tree.png`, `.svg`

---

#### Figure 3 — Epi-Genomic Timeline Scatter

**Input:** Consensus `ClusterAssignment`, `Dict[str, Isolate]`.

**Steps:**
1. X-axis: collection date. Y-axis: cluster ID (categorical, sorted by first detection date).
2. Each point = one isolate, coloured by cluster, shaped by ward.
3. Draw horizontal shaded band per cluster spanning its full date range.
4. Overlay temporal DAG edges as grey arcs above the scatter.
5. Annotate index case candidates.

**Output:** `{output_dir}/figures/epi_timeline_scatter.png`, `.svg`

---

#### Figure 4 — Cluster Comparison Grid (new — no SNP2Cluster equivalent)

**Input:** All three `ClusterAssignment` objects + consensus.

**Steps:**
1. 2×2 grid of subplots: SNP chains, Leiden, Spectral, Consensus.
2. Each subplot is a force-directed network (same layout coordinates), nodes coloured by that method's cluster assignment.
3. Isolates with `ambiguous=True` in consensus are highlighted with a dashed border.
4. Title each subplot with method name + number of clusters detected.

**Output:** `{output_dir}/figures/cluster_comparison_grid.png`, `.svg`

---

#### Figure 5 — Bootstrap Stability Terrain (new)

**Input:** `SNPMatrix`, bootstrap stability scores, consensus `ClusterAssignment`.

**Steps:**
1. Project isolates into 2D using MDS (`sklearn.manifold.MDS(metric=True, dissimilarity="precomputed")`).
2. Compute a Gaussian KDE surface over the 2D space, weighted by stability scores.
3. Render as a filled contour plot (terrain-style): warm colours = high stability, cool = low.
4. Overlay isolate positions as scatter points coloured by cluster.
5. Draw cluster boundary contours at stability threshold 0.90.

**Output:** `{output_dir}/figures/bootstrap_stability_terrain.png`, `.svg`

---

### 6.2 Interactive Live Dashboard (Tier B)

**Input:** Full `PipelineResult` object, all graph objects, centrality scores.

**Framework:** `Dash` (Plotly) with `dash-cytoscape` for network rendering.

**Steps:**
1. App layout with four panels:
   - **Panel 1 (top-left):** Cytoscape network graph. Nodes: coloured by cluster, sized by centrality. Clicking a node highlights all its edges and loads details in Panel 4.
   - **Panel 2 (top-right):** Timeline scatter (same as Figure 3 but interactive — clicking a point filters the network).
   - **Panel 3 (bottom-left):** Cluster summary table (cluster ID, n isolates, STs, wards, date range, stability grade, index case candidate).
   - **Panel 4 (bottom-right):** Selected isolate detail card (all metadata, centrality scores, bootstrap stability, cluster assignment agreement across methods).
2. Add controls:
   - Method toggle (show SNP chains / Leiden / Spectral / Consensus).
   - SNP threshold slider (live re-filters edges).
   - Ward filter multi-select.
   - Date range selector.
3. In streaming mode (Phase 7), Panel 1 updates live as new isolates arrive.
4. Export button: download current view as PDF report.

**Output:**
- `clusterflow/viz/dashboard.py` — Dash app.
- Launched via `clusterflow serve --port 8050`.

---

### Phase 6 Deliverables Summary

| Deliverable | Path | Description |
|---|---|---|
| Heatmap | `clusterflow/viz/heatmap.py` | `seaborn.clustermap` with cluster annotations |
| MST figure | `clusterflow/viz/mst.py` | igraph MST + matplotlib rendering |
| Timeline scatter | `clusterflow/viz/timeline.py` | Epi-genomic scatter with DAG arcs |
| Cluster comparison grid | `clusterflow/viz/comparison.py` | 2×2 method comparison plot |
| Bootstrap terrain | `clusterflow/viz/terrain.py` | MDS + KDE stability terrain |
| Dashboard | `clusterflow/viz/dashboard.py` | Full Dash application |
| Static output runner | `clusterflow/viz/__init__.py` | `generate_all_figures(result, output_dir)` |
| Test suite (Phase 6) | `tests/test_viz.py` | All figures render without error on K. pneumoniae fixtures; dashboard app initialises |

---

## Phase 7 — Streaming / Real-Time Mode

**Goal:** Enable ClusterFlow to accept new isolates as they arrive during an active outbreak and update the transmission graph and cluster assignments incrementally — without reprocessing the full dataset.

### 7.1 Incremental Pipeline

**Input:** Existing `PipelineResult` (in-memory or loaded from disk), new `Isolate` + SNP distances to all existing isolates.

**Steps:**
1. Call `GraphConstructor.add_isolate()` (stubbed in Phase 3) to insert the new node and its edges.
2. Re-run Leiden clustering on the updated graph (`alpha=0.05` — small perturbation, fast convergence).
3. Check if the new isolate joins an existing cluster or forms a new one.
4. Recompute centrality scores for the affected connected component only (not the full graph).
5. If the new isolate's `transmission_risk_score` exceeds a configurable alert threshold, emit a `TRANSMISSION_ALERT` event.
6. Update the `PipelineResult` in-place and push the update to the dashboard via a `dcc.Interval` callback.

**Output:** Updated `PipelineResult` + `TransmissionAlert` event (if triggered).

---

### 7.2 FastAPI Endpoint

**Steps:**
1. Implement a FastAPI application in `clusterflow/streaming/api.py`.
2. Endpoints:

```
POST /isolate
  Body: { isolate_id, collection_date, ward, facility, sequence_type, snp_distances: {id: dist} }
  Response: { cluster_assigned, new_cluster_formed, transmission_alert, centrality_score }

GET /status
  Response: { n_isolates, n_clusters, n_edges, last_update }

GET /result
  Response: Full PipelineResult as JSON

GET /graph
  Response: GraphML string of current transmission graph

WebSocket /stream
  Emits: real-time update events as new isolates arrive
```

3. Mount the Dash dashboard at `/dashboard` (served from the same process).
4. Implement basic API key authentication (configurable; default: disabled for research use).

**Output:** `clusterflow/streaming/api.py` — FastAPI app.

---

### Phase 7 Deliverables Summary

| Deliverable | Path | Description |
|---|---|---|
| Incremental pipeline | `clusterflow/streaming/incremental.py` | Single-isolate graph update |
| FastAPI app | `clusterflow/streaming/api.py` | REST + WebSocket API |
| Alert system | `clusterflow/streaming/alerts.py` | Threshold-based `TRANSMISSION_ALERT` |
| Integration test | `tests/test_streaming.py` | Replay K. pneumoniae dataset isolate-by-isolate; verify final clusters match batch result |

---

## Phase 8 — Validation & Benchmarking

**Goal:** Rigorously demonstrate that ClusterFlow recovers the same clusters as SNP2Cluster on the Kwenda et al. dataset, plus the two additional clusters (ST25 and ST39) that SNP2Cluster missed. Produce a reproducible Jupyter notebook.

### 8.1 Reference Dataset Preparation

**Input:** SNP2Cluster v0.5.4 validation data (K. pneumoniae neonatal unit outbreak, Kwenda et al. 2024).

**Steps:**
1. Download SNP2Cluster v0.5.4 from Zenodo (DOI: 10.5281/zenodo.14060296).
2. Extract the SNP distance matrix, MLST profiles, and epi metadata from the R config files.
3. Convert to ClusterFlow input format (TSV/CSV).
4. Document any format conversions in `docs/data_format_guide.md`.

**Output:** `tests/fixtures/kp_real/` — real K. pneumoniae dataset in ClusterFlow format.

---

### 8.2 Cluster Recovery Benchmark

**Input:** Real K. pneumoniae dataset.

**Steps:**
1. Run SNP2Cluster v0.5.4 (via `rpy2` or subprocess) on the dataset to get the baseline result.
2. Run ClusterFlow on the same dataset.
3. Compute Adjusted Rand Index (ARI) between SNP2Cluster and ClusterFlow consensus assignments.
4. Verify: **ARI ≥ 0.85** for the four main clusters.
5. Verify: ClusterFlow detects ≥ 2 additional clusters not found by SNP2Cluster.
6. Verify: ST25 and ST39 isolates are assigned to distinct clusters in ClusterFlow output.

**Output:**
- `{output_dir}/validation/cluster_comparison.csv` — per-isolate assignment from both tools.
- `{output_dir}/validation/benchmark_summary.json` — ARI, n_clusters per method, detection of ST25/ST39.

---

### 8.3 Performance Benchmark

**Input:** Synthetic datasets of increasing size: 50, 200, 500, 1000, 5000, 10000 isolates (generated by simulation from the K. pneumoniae outbreak model).

**Steps:**
1. Generate synthetic datasets using `clusterflow.testing.simulate_outbreak(n_isolates, n_clusters, snp_rate)`.
2. For each dataset size, time: ingestion, distance computation (CPU backend), graph construction, all three clustering methods, analysis layer.
3. Measure peak memory usage via `tracemalloc`.
4. Plot: runtime vs n_isolates (log-log), memory vs n_isolates.
5. Identify the dataset size at which PairSNP backend becomes necessary (estimated: ~500 isolates for interactive use).

**Output:** `{output_dir}/validation/performance_benchmark.csv`, `performance_plot.png`.

---

### 8.4 Validation Notebook

**Steps:**
1. Create `notebooks/validation_kpneumoniae.ipynb` that:
   - Loads the real K. pneumoniae dataset.
   - Runs ClusterFlow end-to-end.
   - Produces all Phase 6 static figures inline.
   - Compares results to SNP2Cluster.
   - Documents the ST25 and ST39 cluster discovery.
   - Is fully reproducible (`papermill`-compatible, all random seeds fixed).

**Output:** `notebooks/validation_kpneumoniae.ipynb` — runnable with `jupyter nbconvert --execute`.

---

### Phase 8 Deliverables Summary

| Deliverable | Path | Description |
|---|---|---|
| Real K. pneumoniae fixtures | `tests/fixtures/kp_real/` | Converted from SNP2Cluster v0.5.4 |
| Outbreak simulator | `clusterflow/testing/simulate.py` | Synthetic dataset generator |
| Benchmark suite | `tests/test_benchmark.py` | ARI ≥ 0.85; ST25/ST39 detected |
| Performance report | `{output_dir}/validation/` | Runtime + memory by dataset size |
| Validation notebook | `notebooks/validation_kpneumoniae.ipynb` | Fully reproducible end-to-end demo |

---

## Phase 9 — Packaging & Distribution

**Goal:** Make ClusterFlow installable, documented, and deployable with a single command.

### 9.1 Command-Line Interface

**Framework:** `Typer` (type-annotated CLI, auto-generates `--help` docs).

**Commands:**

```bash
# Run full batch pipeline
clusterflow run --config config.yaml --output ./results

# Run in streaming mode (starts API + dashboard)
clusterflow serve --config config.yaml --port 8000

# Validate input files only (dry run)
clusterflow validate --config config.yaml

# Compare two result directories
clusterflow compare ./results_v1 ./results_v2

# Generate a config template
clusterflow init --output config.yaml

# Convert SNP2Cluster config to ClusterFlow format
clusterflow convert-config --snp2cluster-config r_config.R --output clusterflow_config.yaml
```

**Output:** `clusterflow/cli.py` — full Typer CLI.

---

### 9.2 Docker Image

**Steps:**
1. `Dockerfile` with multi-stage build:
   - Stage 1 (`builder`): Install all Python dependencies, compile any C extensions.
   - Stage 2 (`runtime`): Minimal image, copy installed packages, add non-root user.
2. Include `pairsnp` binary in the runtime image (compiled from source).
3. Expose ports 8000 (API) and 8050 (dashboard).
4. `docker-compose.yaml` with optional Redis service for streaming mode job queuing.

**Output:** `Dockerfile`, `docker-compose.yaml`.

---

### 9.3 Documentation

**Framework:** `MkDocs` with `Material` theme.

**Pages:**
- `index.md` — What is ClusterFlow, comparison with SNP2Cluster.
- `installation.md` — pip, conda, Docker.
- `quickstart.md` — 5-minute tutorial with the K. pneumoniae fixture dataset.
- `config_reference.md` — all config keys with descriptions and defaults.
- `data_formats.md` — input format specs, conversion guides.
- `algorithms.md` — technical descriptions of all three clustering methods + consensus.
- `api_reference.md` — auto-generated from docstrings via `mkdocstrings`.
- `streaming.md` — real-time mode setup guide.
- `validation.md` — benchmark results.
- `snp2cluster_migration.md` — migration guide from SNP2Cluster.

**Output:** `docs/` directory, `mkdocs.yaml`. Deployed via GitHub Actions to GitHub Pages.

---

### 9.4 PyPI Release

**Steps:**
1. Set up `pyproject.toml` with `[project]` metadata, version `0.1.0`.
2. Configure `hatch` for versioning (semver, auto-bumped by CI).
3. GitHub Actions workflow: on tag push `v*.*.*`, run tests → build wheel → publish to PyPI.
4. Publish to TestPyPI first; verify install; then publish to PyPI.

**Output:** `pip install clusterflow` installs a working tool.

---

### Phase 9 Deliverables Summary

| Deliverable | Path | Description |
|---|---|---|
| CLI | `clusterflow/cli.py` | Full Typer command-line interface |
| Dockerfile | `Dockerfile` | Multi-stage production image |
| Docker Compose | `docker-compose.yaml` | Full stack with optional Redis |
| Documentation | `docs/` + `mkdocs.yaml` | MkDocs site, all pages |
| CI/CD | `.github/workflows/` | Test, lint, docs-deploy, PyPI-release pipelines |
| SNP2Cluster converter | `clusterflow/cli.py → convert-config` | R config → YAML converter |

---

## Dependency Manifest

```toml
[project]
name = "clusterflow"
requires-python = ">=3.11"

[project.dependencies]
# Core
pydantic = ">=2.0"
pyyaml = ">=6.0"
numpy = ">=1.26"
pandas = ">=2.0"
scipy = ">=1.11"

# Graph
python-igraph = ">=0.11"
leidenalg = ">=0.10"

# Clustering
scikit-learn = ">=1.3"

# Parallel
joblib = ">=1.3"
tqdm = ">=4.66"

# Visualisation (static)
matplotlib = ">=3.8"
seaborn = ">=0.13"

# Visualisation (interactive)
plotly = ">=5.18"
dash = ">=2.14"
dash-cytoscape = ">=0.3"

# CLI
typer = ">=0.9"
rich = ">=13.0"            # pretty terminal output

# API
fastapi = ">=0.104"
uvicorn = ">=0.24"
websockets = ">=12.0"

[project.optional-dependencies]
gpu = [
  "cupy-cuda12x>=13.0",    # or cupy-cuda11x depending on CUDA version
]
r-bridge = [
  "rpy2>=3.5",
]
dev = [
  "pytest>=7.4",
  "pytest-cov>=4.1",
  "pytest-asyncio>=0.21",
  "ruff>=0.1",
  "black>=23.0",
  "mypy>=1.6",
  "pre-commit>=3.5",
  "papermill>=2.5",
  "mkdocs>=1.5",
  "mkdocs-material>=9.4",
  "mkdocstrings[python]>=0.24",
]
```

---

## Testing Standards

All phases must meet these standards before the next phase begins:

- **Unit test coverage:** ≥ 85% line coverage for all non-visualisation code.
- **Integration tests:** End-to-end pipeline run on K. pneumoniae fixture dataset must complete without error.
- **Type coverage:** All public functions and classes must have type annotations; `mypy --strict` must pass.
- **Performance gate:** Full batch pipeline on the 28-isolate K. pneumoniae dataset must complete in < 60 seconds on a standard 4-core laptop.
- **Reproducibility gate:** Two runs with the same config and fixed random seeds must produce bit-identical results.

---

## Suggested Build Order for a Coding Agent

If implementing phases sequentially, the recommended order is:

1. **Phase 1.1 → 1.4** — Foundation. No prior phases required.
2. **Phase 3** — Graph Constructor. Requires Phase 1 models + ingestion. Can use a pre-computed test matrix to avoid needing Phase 2.
3. **Phase 4.1** — SNP Chain algorithm only. Validates graph construction against SNP2Cluster baseline.
4. **Phase 2** — Distance Engine. Can now be validated against Phase 4.1 output.
5. **Phase 4.2, 4.3, 4.4** — Leiden, Spectral, Consensus.
6. **Phase 5** — Analysis Layer. All clustering must be complete.
7. **Phase 6.1** — Static visualisations. Produces validation figures.
8. **Phase 8** — Validation. Run before dashboard/streaming to confirm correctness.
9. **Phase 6.2** — Dashboard. Requires stable analysis outputs.
10. **Phase 7** — Streaming. Requires stable graph + clustering code.
11. **Phase 9** — Packaging. Final step.

---

## Key Design Decisions for the Coding Agent

1. **SNP2Cluster compatibility is non-negotiable for Phase 4.1.** The SNP chain algorithm must produce identical results to SNP2Cluster v0.5.4 on the same input. This is the trust anchor for the whole project.

2. **Graph is always `igraph`, never `networkx`.** `python-igraph` is 10–100× faster for community detection and centrality at scale. `networkx` is only used if a dependency forces it, and even then only in conversion utilities.

3. **All randomness must be seedable.** Leiden clustering, spectral K-means, bootstrap resampling, and synthetic data generation must all accept an integer `random_state` parameter. Default seed = `42`. This is enforced by the reproducibility gate in the test standards above.

4. **Fail loudly at ingestion, gracefully downstream.** Missing input files, malformed matrices, and ID mismatches are hard errors with informative messages. Partial MLST coverage, ambiguous cluster assignments, and bootstrap instability are warnings with actionable suggestions — not errors.

5. **The YAML config is the single source of truth.** No hardcoded thresholds anywhere in the codebase. All defaults live in the Pydantic model's `Field(default=...)` declarations and are documented in `configs/example.yaml`.

---

*Document version: 1.0 — Generated from ClusterFlow design sessions. For questions, refer to the project conversation history.*
