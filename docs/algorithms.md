# Algorithms

ClusterFlow runs three clustering algorithms in parallel on the same weighted graph and combines them into a single consensus assignment via Hungarian label alignment + majority vote.

## The transmission graph

Vertices are isolates; their attributes carry every field of the `Isolate` model (collection date, facility, ward, ST). An undirected edge `(i, j)` is added iff:

```
SNP_distance(i, j) ≤ snp_cutoff   AND   |Δdays(i, j)| ≤ day_cutoff
```

Each edge carries the four attributes used downstream:

```
snp_distance     # integer SNP count
day_delta        # |date_i − date_j|
mlst_mismatch    # bool: 1 iff ST_i ≠ ST_j
composite_weight # α·(SNP/snp_cutoff) + β·(days/day_cutoff) + γ·mismatch
```

Lower composite weight = closer (more similar). Some algorithms operate on a similarity transform `exp(−composite_weight)`.

## Algorithm 1 — SNP chains

Faithful re-implementation of SNP2Cluster's chain logic. Used as the trust-anchor / SNP2Cluster-compatibility baseline.

1. Sort isolates by collection date.
2. Walk in date order. For each isolate, scan existing clusters; if any current member is within `snp_cutoff` SNPs **and** within `day_cutoff` days, join that cluster. First qualifying assignment wins.
3. Else start a new cluster.

## Algorithm 2 — Leiden community detection

Uses the [`leidenalg`](https://github.com/vtraag/leidenalg) implementation with `RBConfigurationVertexPartition` and similarity weights.

- Edge weights are `exp(−composite_weight)` so closer isolates pull harder.
- Default `leiden_resolution = 1.0`. Setting it to `"auto"` sweeps the resolution from 0.5 to 2.0 in steps of 0.1 in parallel and picks the resolution that maximises modularity Q.
- 10 optimisation iterations per run.

## Algorithm 3 — Spectral clustering with eigengap

1. Build the symmetric normalised graph Laplacian `L_sym = I − D^{−1/2} A D^{−1/2}` from the affinity matrix `A_ij = exp(−composite_weight_ij)`.
2. Eigendecompose `L_sym`. Pick *k* by the eigengap heuristic — the index where consecutive eigenvalues jump fastest.
3. Project isolates into the k smallest eigenvectors, row-normalise (Ng-Jordan-Weiss).
4. Run `KMeans(n_clusters=k, random_state=42, n_init=10)` in this eigenspace.
5. Silhouette-refine: try `k−1`, `k`, `k+1` and pick the *k* with the highest silhouette score.

## Consensus

Cluster IDs are arbitrary across methods (method A's "cluster 3" need not match method B's "cluster 3"). The consensus aligns labels with the [Hungarian algorithm](https://en.wikipedia.org/wiki/Hungarian_algorithm):

1. Compute the overlap matrix between two assignments' cluster IDs.
2. Solve the linear assignment that maximises overlap.
3. Apply the resulting permutation to canonicalise labels.
4. Per isolate, take the majority vote across all aligned methods.
5. `agreement_score(iso) = (votes_for_winner) / n_methods`.
6. Mark an isolate `ambiguous=True` if its agreement score is `< 2/3`.

## Temporal DAG (Phase 5.1)

For each undirected edge `(A, B)`:

- If `date_A < date_B − uncertainty_days` → directed edge A → B.
- If `date_B < date_A − uncertainty_days` → directed edge B → A.
- Otherwise (within ±uncertainty) keep both directions.

Bidirectional edges can introduce cycles. A simple DFS finds any cycle; we delete the edge with the highest `composite_weight` (= lowest similarity) in the cycle, then re-check. This converges in O(cycles × |E|).

## Centrality and index case ranking (Phase 5.2)

Computed on the directed DAG, weighted by `composite_weight`:

- **Betweenness** — fraction of shortest paths through this node.
- **In-degree / out-degree** — incoming/outgoing transmission events.
- **Closeness** — reciprocal of mean distance to every reachable isolate.

These are normalised within each cluster, then combined:

```
transmission_risk_score = 0.5·betweenness + 0.3·out_degree_norm + 0.2·closeness
```

For each cluster, the highest-scoring isolate is the **index case candidate**. If the candidate's collection-date rank within its cluster is in the latter 75 %, a `WARNING: temporal mismatch` is emitted because high centrality + late collection is suspicious.

## Bootstrap stability (Phase 5.3)

Cluster assignments are not point estimates — they should come with confidence. We perturb the SNP matrix and re-cluster.

1. For `bootstrap_n` replicates (default 500):
   - Add uniform noise `±2 SNPs` to every entry of the upper triangle, mirror to maintain symmetry, zero the diagonal.
   - Re-build the graph with the perturbed matrix.
   - Re-run Leiden (fastest of the three).
   - Hungarian-align labels back to the consensus.
2. For each isolate, `stability_score = fraction of replicates in which the isolate landed in the same (aligned) cluster as in the consensus`.
3. Classify:
   - `stable` (≥ 0.90)
   - `uncertain` (0.70–0.90)
   - `unstable` (< 0.70)
4. Cluster confidence grade: A (mean ≥ 0.90), B (≥ 0.75), C (< 0.75).

Replicates run in parallel via `joblib.Parallel(backend="loky")`. Two runs with the same seed produce bit-identical scores.
