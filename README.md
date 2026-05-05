# ClusterFlow

Transmission cluster detection pipeline that fuses whole-genome sequencing (WGS) data
with epidemiological metadata using graph-theoretic algorithms and parallel computing.

We are introducing Leiden community detection, spectral clustering, temporal DAGs,
centrality-based index case identification, and a real-time streaming mode with a
live surveillance dashboard.

**Validation dataset:** *Klebsiella pneumoniae* neonatal unit outbreak (Kwenda et al. 2024).

## Status

All nine phases implemented and validated end-to-end on the SNP2Cluster v0.5.4
example *K. pneumoniae* dataset (ST25 and ST39 detected as distinct clusters).
See [transmission_cluster_tool_plan.md](transmission_cluster_tool_plan.md) for
the full design plan and [docs/](docs/) for usage.

## Three ways to run

```bash
# 1. CLI
clusterflow run --config config.yaml --output ./results

# 2. Browser (Streamlit)
pip install "clusterflow[web]"
clusterflow web

# 3. Live surveillance API + dashboard
pip install "clusterflow[streaming,dashboard]"
clusterflow serve --config config.yaml
```

Every batch run emits a self-contained `report.html` you can email or attach
to a pull request.

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
