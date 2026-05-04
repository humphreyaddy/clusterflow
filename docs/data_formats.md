# Data formats

ClusterFlow's three input files all use plain CSV/TSV. Files are auto-detected by the delimiter sniffer, so `.tsv`, `.csv`, and pipe-separated all work.

## SNP distance matrix

A square symmetric matrix with isolate IDs as both row labels and column headers. The first column must contain isolate IDs; the first cell may be empty or contain a banner string (`snp-dists 0.8.2` is recognised and stripped).

```tsv
		EXM113	EXM218	EXM222
EXM113	0	1	2
EXM218	1	0	3
EXM222	2	3	0
```

**Validation rules** (errors include the offending row/column):

- Square (n × n)
- Symmetric within tolerance 1e-6
- Diagonal exactly zero
- All values non-negative
- No `NaN` or non-numeric cells

The reader symmetrises the matrix to remove float jitter before returning.

## Epidemiological metadata

CSV/TSV with the four required columns (case-insensitive headers):

```csv
isolate_id,collection_date,facility,ward
EXM113,2019-10-17,Hospital_1,neonatal
EXM218,2019-12-17,Hospital_1,neonatal
```

Date formats supported: `YYYY-MM-DD`, `DD/MM/YYYY`, `DD-MM-YYYY`, `YYYY/MM/DD` plus the generic pandas parser as a fallback. Future dates raise an error; spans > 5 years emit a warning.

## MLST profiles (optional)

CSV/TSV with `isolate_id` and `ST` columns. Allele columns are accepted but ignored. The `ST` value is normalised: `"ST405"` → `"405"`, `"-"`/blank/`"novel"` → `"novel"`.

```csv
isolate_id,ST,gapA,infB,mdh,pgi,phoE,rpoB,tonB
EXM113,17,2,1,1,1,4,4,4
EXM218,152,5,1,1,1,32,1,2
```

Partial MLST coverage emits an info log but is not an error.

## Cross-validation

`IngestionPipeline.load(config)` runs all three readers and then:

- Errors out if any isolate ID in the SNP matrix is absent from the epi metadata.
- Drops (and warns) any epi rows that don't appear in the SNP matrix.
- Logs MLST coverage as a percentage.

This is the same logic invoked by `clusterflow validate`.

## Converting from SNP2Cluster

The SNP2Cluster v0.5.4 example data ships in three files inside the Zenodo archive (DOI 10.5281/zenodo.14060296):

| SNP2Cluster file | Convert to | Note |
|---|---|---|
| `coreSNPmatrix.co.csv` | `snp_matrix.csv` | Drop the `reference` row + column. The first cell will say `snp-dists 0.8.2`; replace with empty. |
| `example_metadata.csv` | `epi.csv` | Rename: `SampleID → isolate_id`, `FacilityName → facility`, `WardType → ward`, `TakenDate → collection_date`. Reformat date as `YYYY-MM-DD`. |
| `05.mlst.xlsx` | `mlst.csv` | Rename `FILE → isolate_id`. Keep `ST` column. |

A worked example is at `tests/fixtures/kp_real/` and the conversion is reproduced in `notebooks/validation_kpneumoniae.ipynb`.
