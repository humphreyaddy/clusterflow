# API reference

Auto-generated from docstrings via [`mkdocstrings`](https://mkdocstrings.github.io/).

## Top-level pipeline

::: clusterflow.pipeline.run_pipeline

## Data models

::: clusterflow.models

## Configuration

::: clusterflow.config

## Ingestion

::: clusterflow.ingestion.IngestionPipeline
::: clusterflow.ingestion.SNPMatrixReader
::: clusterflow.ingestion.MLSTReader
::: clusterflow.ingestion.EpiMetadataReader

## Distance engine

::: clusterflow.distance.select_backend
::: clusterflow.distance.CPUBackend

## Graph

::: clusterflow.graph.GraphConstructor
::: clusterflow.graph.graph_summary
::: clusterflow.graph.save_graph
::: clusterflow.graph.load_graph

## Clustering

::: clusterflow.clustering.snp_chain_clusters
::: clusterflow.clustering.leiden_clusters
::: clusterflow.clustering.spectral_clusters
::: clusterflow.clustering.run_all
::: clusterflow.clustering.consensus_assignment

## Analysis

::: clusterflow.analysis.build_temporal_dag
::: clusterflow.analysis.compute_centrality
::: clusterflow.analysis.index_case_candidates
::: clusterflow.analysis.bootstrap_stability

## Visualisation

::: clusterflow.viz.generate_all_figures

## Streaming

::: clusterflow.streaming.add_isolate
::: clusterflow.streaming.IncrementalState
::: clusterflow.streaming.IncrementalUpdate

## Testing helpers

::: clusterflow.testing.simulate_outbreak
::: clusterflow.testing.run_benchmark
