# Migrating from SNP2Cluster

If you're already running SNP2Cluster v0.5.x, ClusterFlow is designed as a drop-in upgrade — the SNP-chain algorithm in Phase 4.1 produces identical assignments on the same inputs, and the YAML config mirrors SNP2Cluster's R config keys.

## Mapping config keys

| SNP2Cluster (R) | ClusterFlow (YAML) |
|---|---|
| `snpco = 20` | `thresholds.snp_cutoff: 20` |
| `daysco = 45` | `thresholds.day_cutoff: 45` |
| `coreSNPmatrix = "..."` | `inputs.snp_matrix: "..."` |
| `MetadataFile = "..."` | `inputs.epi_metadata: "..."` |
| `MLSTfile = "..."` | `inputs.mlst_profiles: "..."` |
| `OutputFolder = "..."` | `output_dir: "..."` |

A minimal ClusterFlow config that mirrors SNP2Cluster's defaults:

```yaml
project_name: "snp2cluster_compat"
output_dir: "./results"

inputs:
  snp_matrix: "data/coreSNPmatrix.co.csv"
  mlst_profiles: "data/mlst.csv"
  epi_metadata: "data/metadata.csv"

thresholds:
  snp_cutoff: 20
  day_cutoff: 45
  edge_weight_alpha: 1.0
  edge_weight_beta: 0.0
  edge_weight_gamma: 0.0

clustering:
  methods: ["snp_chains"]   # SNP2Cluster-equivalent only
  bootstrap_n: 0
  random_state: 42

visualization:
  static: true
  dashboard: false
```

Run with `clusterflow run --config snp2cluster_compat.yaml`.

## Mapping outputs

| SNP2Cluster output | ClusterFlow equivalent |
|---|---|
| Cluster assignment table | `clusters/cluster_assignments.csv` |
| SNP heatmap PDF | `figures/snp_heatmap.png` (+ SVG) |
| Minimum spanning tree PDF | `figures/minimum_spanning_tree.png` |
| Epi-genomic timeline | `figures/epi_timeline_scatter.png` |
| — | `figures/cluster_comparison_grid.png` *(new)* |
| — | `figures/bootstrap_stability_terrain.png` *(new)* |
| — | `analysis/centrality_scores.csv` *(new)* |
| — | `analysis/transmission_dag.graphml` *(new)* |
| — | `analysis/bootstrap_stability.csv` *(new)* |

## Converting input files

Your existing SNP-distance matrix from snp-dists works unmodified — just delete the `reference` row + column if present. The ClusterFlow reader recognises the `snp-dists 0.8.2` banner string and strips it.

For the metadata, rename four columns:

| SNP2Cluster | ClusterFlow |
|---|---|
| `SampleID` | `isolate_id` |
| `FacilityName` | `facility` |
| `WardType` | `ward` |
| `TakenDate` | `collection_date` |

Date format `YYYY/MM/DD` is auto-detected.

For MLST, rename `FILE` → `isolate_id`. The `ST` column maps directly. Allele columns are accepted but ignored by ClusterFlow.

## Behavioural differences to be aware of

1. **Edge construction is undirected by default.** SNP2Cluster builds chains; ClusterFlow builds an undirected graph + a separate temporal DAG (Phase 5.1). The directed view is what `analysis/transmission_dag.graphml` exposes.

2. **Three algorithms, not one.** Even if you only care about the SNP-chain output, ClusterFlow runs Leiden + spectral by default to compute a consensus + ambiguity flag. To match SNP2Cluster's single-algorithm behaviour, set `clustering.methods: ["snp_chains"]`.

3. **Bootstrap is on by default (500 reps).** Set `clustering.bootstrap_n: 0` to skip.

4. **All randomness is seeded.** `clustering.random_state: 42` means two runs with the same config produce bit-identical assignments. SNP2Cluster's K-means seed is implicit and platform-dependent.

5. **Output format is CSV, not Excel.** Use pandas/Excel/R to load. GraphML files open directly in Gephi or Cytoscape.
