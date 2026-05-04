# Streaming / real-time mode

ClusterFlow can accept new isolates as they arrive during an active outbreak and update the transmission graph incrementally — without reprocessing the full dataset.

## Setup

```bash
pip install "clusterflow[streaming,dashboard]"
clusterflow serve --config config.yaml --port 8000
```

The `serve` command:

1. Loads the existing dataset from `config.inputs.*` and builds the initial graph + Leiden consensus.
2. Mounts the FastAPI app at `http://0.0.0.0:8000`.
3. (Optionally) mounts the Dash dashboard at port 8050.

## Endpoints

### `POST /isolate`

Add a single new isolate to the graph and re-cluster the affected component.

**Body**

```json
{
  "isolate_id": "EXM2999",
  "collection_date": "2026-01-15",
  "facility": "Hospital_1",
  "ward": "neonatal",
  "sequence_type": "152",
  "snp_distances": {
    "EXM218": 4,
    "EXM222": 7,
    "EXM254": 3
  }
}
```

You only need to supply distances to existing isolates — the API computes the new edges itself.

**Response**

```json
{
  "cluster_assigned": 2,
  "new_cluster_formed": false,
  "transmission_alert": false,
  "centrality_score": 0.21,
  "new_edges": [
    {"source": "EXM2999", "target": "EXM218", "snp_distance": 4, ...}
  ]
}
```

If the new isolate's `transmission_risk_score` ≥ `alert_threshold` (default 0.7), `transmission_alert` is `true` and a `TRANSMISSION_ALERT` is logged.

### `GET /status`

```json
{ "n_isolates": 37, "n_clusters": 13, "n_edges": 121, "last_update": "2026-01-15T10:42:08Z" }
```

### `GET /graph`

Full vertex + edge dump as JSON (every attribute preserved).

### `GET /result`

Current consensus clustering as JSON — useful for downstream tools.

## Dashboard

Open <http://localhost:8050> after starting the server. The four-panel layout:

| Panel | Content |
|---|---|
| Top-left | Cytoscape network. Click a node to load its detail. |
| Top-right | Plotly timeline scatter (isolate × date, coloured by cluster). |
| Bottom-left | Cluster summary table. |
| Bottom-right | Selected isolate detail card. |

When new isolates POST in, the network updates live.

## Programmatic incremental updates

If you'd rather drive the graph from a Python process (e.g. a Snakemake/nextflow pipeline) without the HTTP layer:

```python
from clusterflow.config import ClusterFlowConfig
from clusterflow.streaming import IncrementalState, add_isolate
from clusterflow.graph import GraphConstructor
from clusterflow.ingestion import IngestionPipeline
from clusterflow.clustering.leiden import leiden_clusters
from clusterflow.clustering.consensus import consensus_assignment

cfg = ClusterFlowConfig.from_yaml("config.yaml")
ingest = IngestionPipeline().load(cfg)
g = GraphConstructor(cfg.thresholds).build(ingest.snp, ingest.isolates)
state = IncrementalState(
    config=cfg, graph=g,
    isolates=dict(ingest.isolates),
    consensus=consensus_assignment([leiden_clusters(g, n_jobs=1)]),
)

# Loop: receive new isolates from your queue / message bus
for new_iso, snp_distances in ...:
    update = add_isolate(state, new_iso, snp_distances)
    if update.transmission_alert:
        notify_oncall(new_iso.isolate_id, update.cluster_assigned)
```

## Authentication

The default build has no auth — fine for local development and air-gapped surveillance networks. For internet-facing deployments, put a reverse proxy (nginx, Caddy, Cloudflare Access) in front, or wire FastAPI's [security dependencies](https://fastapi.tiangolo.com/tutorial/security/) into `clusterflow.streaming.api.create_app`.
