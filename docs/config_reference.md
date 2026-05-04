# Configuration reference

Every configurable knob lives in a single YAML file validated by a Pydantic model. There are no hardcoded thresholds anywhere in the codebase — defaults live in `clusterflow/config.py` and are documented here.

## Top level

| Key | Type | Default | Description |
|---|---|---|---|
| `project_name` | str | *(required)* | Used in figure titles and the summary JSON. |
| `output_dir` | path | *(required)* | Where every artefact is written. Created if missing. |
| `inputs` | section | — | See below. |
| `thresholds` | section | defaults | See below. |
| `clustering` | section | defaults | See below. |
| `distance_engine` | section | defaults | See below. |
| `streaming` | section | defaults | See below. |
| `visualization` | section | defaults | See below. |

## `inputs`

| Key | Type | Required | Description |
|---|---|---|---|
| `snp_matrix` | path | yes | TSV/CSV pairwise SNP matrix (square, symmetric, zero diagonal). |
| `mlst_profiles` | path | no | CSV with `isolate_id`, `ST`. Optional; matches by `isolate_id`. |
| `epi_metadata` | path | yes | CSV with `isolate_id`, `collection_date`, `facility`, `ward`. |

## `thresholds`

| Key | Type | Default | Description |
|---|---|---|---|
| `snp_cutoff` | int ≥ 0 | 20 | Max SNP distance for an edge. |
| `day_cutoff` | int ≥ 0 | 14 | Max collection-date gap (days) for an edge. |
| `edge_weight_alpha` | 0–1 | 0.4 | Weight on normalised SNP distance term. |
| `edge_weight_beta` | 0–1 | 0.4 | Weight on normalised temporal term. |
| `edge_weight_gamma` | 0–1 | 0.2 | Weight on MLST mismatch (0/1) term. |

The three weights must sum to 1.0; the validator rejects configs that don't.

The composite edge weight is

```
weight = α·(SNP / snp_cutoff) + β·(|Δdays| / day_cutoff) + γ·𝟙[ST_i ≠ ST_j]
```

## `clustering`

| Key | Type | Default | Description |
|---|---|---|---|
| `methods` | list | all 3 | Subset of `["snp_chains", "leiden", "spectral"]`. |
| `bootstrap_n` | int ≥ 0 | 500 | Replicates for the stability CI. Set to 0 to skip. |
| `leiden_resolution` | float \| `"auto"` | 1.0 | Resolution parameter. `"auto"` sweeps 0.5–2.0. |
| `random_state` | int | 42 | Seed for every stochastic step. |

## `distance_engine`

| Key | Type | Default | Description |
|---|---|---|---|
| `backend` | enum | `auto` | `auto` \| `cpu` \| `gpu` \| `pairsnp`. |
| `n_jobs` | int | -1 | -1 → use all available cores. |

## `streaming`

| Key | Type | Default | Description |
|---|---|---|---|
| `enabled` | bool | false | Whether `clusterflow serve` should run by default. |
| `api_port` | int | 8000 | FastAPI port. |

## `visualization`

| Key | Type | Default | Description |
|---|---|---|---|
| `static` | bool | true | Generate the 5 static figures. |
| `dashboard` | bool | false | Launch the Dash dashboard. |
| `dashboard_port` | int | 8050 | Dash port. |

## Validating a config

```bash
clusterflow validate --config config.yaml
```

Prints `[ok]  isolates=N  snp_matrix=N×N` if everything cross-validates, or a descriptive error otherwise.

## Diffing configs

```bash
clusterflow compare config_v1.yaml config_v2.yaml
```

Renders a Rich table of every key whose value differs.
