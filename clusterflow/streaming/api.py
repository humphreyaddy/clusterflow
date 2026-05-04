"""Phase 7.2 — FastAPI endpoint for streaming isolates.

This module raises :class:`ImportError` at module load if FastAPI is missing,
so the CLI can detect the absence cleanly.
"""

from __future__ import annotations

import logging
from datetime import date, datetime

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

from clusterflow.config import ClusterFlowConfig
from clusterflow.graph import GraphConstructor
from clusterflow.ingestion import IngestionPipeline
from clusterflow.models import Isolate
from clusterflow.streaming.incremental import (
    IncrementalState,
    IncrementalUpdate,
    add_isolate,
)
from clusterflow.clustering.consensus import consensus_assignment
from clusterflow.clustering.leiden import leiden_clusters

log = logging.getLogger(__name__)


class IsolatePayload(BaseModel):
    isolate_id: str
    collection_date: date
    facility: str
    ward: str
    sequence_type: str | None = None
    snp_distances: dict[str, float]


class StatusResponse(BaseModel):
    n_isolates: int
    n_clusters: int
    n_edges: int
    last_update: datetime | None


def create_app(config: ClusterFlowConfig) -> FastAPI:
    """Build the FastAPI app, pre-loading the existing dataset from config."""
    app = FastAPI(title="ClusterFlow", version="0.1.0")

    ingest = IngestionPipeline().load(config)
    g = GraphConstructor(config.thresholds).build(ingest.snp, ingest.isolates)
    initial = leiden_clusters(
        g,
        resolution=config.clustering.leiden_resolution,
        random_state=config.clustering.random_state,
        n_jobs=1,
    )
    state = IncrementalState(
        config=config,
        graph=g,
        isolates=dict(ingest.isolates),
        consensus=consensus_assignment([initial]),
    )
    last_update: dict[str, datetime | None] = {"ts": None}

    @app.post("/isolate")
    def post_isolate(payload: IsolatePayload) -> IncrementalUpdate:
        if payload.isolate_id in state.isolates:
            raise HTTPException(409, f"isolate {payload.isolate_id} already present")
        iso = Isolate(
            isolate_id=payload.isolate_id,
            collection_date=payload.collection_date,
            facility=payload.facility,
            ward=payload.ward,
            sequence_type=payload.sequence_type,
        )
        update = add_isolate(state, iso, payload.snp_distances)
        last_update["ts"] = datetime.utcnow()
        return update

    @app.get("/status", response_model=StatusResponse)
    def status() -> StatusResponse:
        return StatusResponse(
            n_isolates=state.graph.vcount(),
            n_clusters=state.consensus.n_clusters,
            n_edges=state.graph.ecount(),
            last_update=last_update["ts"],
        )

    @app.get("/graph")
    def graph_dump() -> dict:
        return {
            "n_vertices": state.graph.vcount(),
            "n_edges": state.graph.ecount(),
            "vertices": [
                {k: v[k] for k in state.graph.vs.attributes()}
                for v in state.graph.vs
            ],
            "edges": [
                {
                    "source": e.source,
                    "target": e.target,
                    **{k: e[k] for k in state.graph.es.attributes()},
                }
                for e in state.graph.es
            ],
        }

    @app.get("/result")
    def result() -> dict:
        return {
            "consensus": state.consensus.model_dump(),
            "n_isolates": len(state.isolates),
            "n_edges": state.graph.ecount(),
        }

    return app
