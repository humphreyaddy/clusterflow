# Installation

ClusterFlow targets Python 3.11+ and works on Linux, macOS (Apple Silicon), and Windows.

## pip

```bash
pip install clusterflow
```

Optional extras:

```bash
pip install "clusterflow[dashboard]"    # dash + dash-cytoscape + plotly
pip install "clusterflow[streaming]"    # fastapi + uvicorn + websockets
pip install "clusterflow[gpu]"          # cupy CUDA backend (Linux/Windows only)
pip install "clusterflow[r-bridge]"     # rpy2 bridge for SNP2Cluster comparison
pip install "clusterflow[dev]"          # pytest, ruff, black, mypy, mkdocs
```

## From source

```bash
git clone https://github.com/your-org/clusterflow
cd clusterflow
pip install -e ".[dev,dashboard,streaming]"
```

## Conda / Mamba

```bash
mamba env create -f environment.yml
mamba activate clusterflow
pip install -e .
```

## Docker

```bash
docker build -t clusterflow:latest .
docker run --rm -v $PWD/data:/data clusterflow:latest \
    run --config /data/config.yaml --output /data/results
```

For the full streaming + dashboard stack:

```bash
docker compose up -d
# API:       http://localhost:8000
# Dashboard: http://localhost:8050
```

## Verifying the install

```bash
clusterflow --help
clusterflow simulate -o /tmp/kp_demo
clusterflow run --config /tmp/kp_demo/kp_config.yaml --output /tmp/kp_results
```

If the `run` step finishes with a printed `pipeline summary` table, you're ready.
