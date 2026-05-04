# Quickstart

A 5-minute tour using the bundled synthetic *K. pneumoniae* outbreak fixture (28 isolates, 4 simulated clusters).

## 1. Generate or supply input data

```bash
clusterflow simulate -o tests/fixtures
```

This writes three CSV/TSV files plus a ground-truth label file:

```
tests/fixtures/
├── kp_snp_matrix.tsv     # 28×28 pairwise SNP distances
├── kp_epi.csv            # isolate_id, collection_date, facility, ward
├── kp_mlst.csv           # isolate_id, ST
└── kp_truth.csv          # isolate_id, true_cluster   (for validation)
```

If you already have your own data, see **[Data formats](data_formats.md)**.

## 2. Write a config

```yaml title="config.yaml"
project_name: "kp_neonatal_2024"
output_dir: "./results"

inputs:
  snp_matrix: "tests/fixtures/kp_snp_matrix.tsv"
  mlst_profiles: "tests/fixtures/kp_mlst.csv"
  epi_metadata: "tests/fixtures/kp_epi.csv"

thresholds:
  snp_cutoff: 20
  day_cutoff: 30
  edge_weight_alpha: 0.4
  edge_weight_beta: 0.4
  edge_weight_gamma: 0.2

clustering:
  methods: ["snp_chains", "leiden", "spectral"]
  bootstrap_n: 100
  leiden_resolution: 1.0
  random_state: 42
```

## 3. Run

```bash
clusterflow run --config config.yaml --output ./results
```

Expected timing on a 4-core laptop: well under 10 s for the 28-isolate fixture, ~30 s for 500 isolates.

## 4. Inspect outputs

```
results/
├── graph/
│   ├── transmission_graph.graphml   # open in Gephi / Cytoscape
│   ├── transmission_graph.pkl       # fast Python re-load
│   └── graph_summary.json
├── clusters/
│   └── cluster_assignments.csv      # one row per isolate, columns per method + consensus
├── analysis/
│   ├── transmission_dag.graphml     # directed temporal DAG
│   ├── dag_summary.json
│   ├── centrality_scores.csv
│   ├── index_case_candidates.csv
│   ├── bootstrap_stability.csv
│   └── cluster_stability_summary.csv
├── figures/
│   ├── snp_heatmap.{png,svg}
│   ├── minimum_spanning_tree.{png,svg}
│   ├── epi_timeline_scatter.{png,svg}
│   ├── cluster_comparison_grid.{png,svg}
│   └── bootstrap_stability_terrain.{png,svg}
└── pipeline_summary.json
```

## 5. (Optional) Live dashboard

```bash
pip install "clusterflow[dashboard,streaming]"
clusterflow serve --config config.yaml --port 8000
```

Then open <http://localhost:8050>. POST new isolates to `http://localhost:8000/isolate` to add them to the graph in real time.

## What's next?

- **[Configuration reference](config_reference.md)** — every config key explained
- **[Algorithms](algorithms.md)** — how the three clustering methods + consensus work
- **[Validation results](validation.md)** — how ClusterFlow performs on the SNP2Cluster benchmark
