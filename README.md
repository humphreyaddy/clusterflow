# ClusterFlow

Transmission cluster detection pipeline that fuses whole-genome sequencing (WGS) data
with epidemiological metadata using graph-theoretic algorithms and parallel computing.

Supersedes the K-means + SNP-chain approach of [SNP2Cluster](https://github.com/stanikae/SNP2Cluster)
by introducing Leiden community detection, spectral clustering, temporal DAGs,
centrality-based index case identification, and a real-time streaming mode with a
live surveillance dashboard.

**Validation dataset:** *Klebsiella pneumoniae* neonatal unit outbreak (Kwenda et al. 2024).

## Status

Phase 1 scaffolding only. See [transmission_cluster_tool_plan.md](transmission_cluster_tool_plan.md)
for the full phased implementation plan.

## Layout

```
clusterflow/        # Python package (config, models, ingestion, distance, graph,
                    # clustering, analysis, viz, streaming, cli)
tests/              # Pytest suite + fixtures
configs/            # Example YAML configs
notebooks/          # Validation + demo notebooks
docs/               # MkDocs site
```

## Install (development)

```bash
make install
make test
```
