# ClusterFlow

> Transmission cluster detection from whole-genome sequencing + epidemiological metadata.

ClusterFlow supersedes the K-means + SNP-chain logic of [SNP2Cluster](https://github.com/stanikae/SNP2Cluster) by introducing:

- **Leiden community detection** on a weighted temporal-genomic graph
- **Spectral clustering** with eigengap-driven *k* selection
- **Hungarian-aligned consensus** across all three algorithms
- **Temporal DAGs** with cycle-breaking for directed transmission inference
- **Centrality-based index case identification** (betweenness + out-degree + closeness)
- **Bootstrap confidence intervals** for every cluster assignment
- **A real-time streaming mode** for active outbreak surveillance

## Why a new tool?

SNP2Cluster's K-means baseline misses small/recent clusters that haven't yet
acquired enough SNP diversity to be separable in linear PCA space. On the
*Klebsiella pneumoniae* neonatal-unit example dataset shipped with SNP2Cluster
v0.5.4, ClusterFlow recovers the dominant ST152 outbreak **and** the secondary
ST25 and ST39 transmission clusters that K-means collapses into noise.

## Quick example

```bash
clusterflow simulate -o tests/fixtures
clusterflow run --config tests/fixtures/kp_config.yaml --output ./results
```

That generates the full output tree under `./results/`:

```
results/
├── graph/                       # GraphML + pickle + summary JSON
├── clusters/                    # per-method assignments
├── analysis/                    # DAG, centrality, bootstrap stability
├── figures/                     # 5 publication-ready figures (PNG + SVG)
└── pipeline_summary.json
```

See **[Quickstart](quickstart.md)** to walk through it.

## Architecture

```
Inputs (SNP matrix · MLST · Epi metadata)
        ↓
[Phase 2] Parallel Distance Engine
        ↓
[Phase 3] Graph Constructor (weighted temporal graph)
        ↓ ↓ ↓
[Phase 4] Cluster Detection (SNP chains · Leiden · Spectral + Consensus)
        ↓ ↓ ↓
[Phase 5] Analysis Layer (Temporal DAG · Centrality · Bootstrap CIs)
        ↓
[Phase 6] Visualization Suite (Heatmap · MST · Transmission Tree · Dashboard)
        ↓
[Phase 7] Streaming Mode (Incremental graph · FastAPI endpoint)
```

Each phase is implemented as a separate sub-package under `clusterflow/` and
can be swapped, extended, or driven directly from Python.
