"""Phase 3.2 — graph serialisation (GraphML + pickle + summary JSON)."""

from __future__ import annotations

import json
import pickle
from pathlib import Path

import igraph as ig

from clusterflow.graph.constructor import graph_summary


def save_graph(g: ig.Graph, output_dir: str | Path, prefix: str = "transmission_graph") -> dict:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    graphml_path = out / f"{prefix}.graphml"
    pickle_path = out / f"{prefix}.pkl"
    summary_path = out / "graph_summary.json"

    g.write_graphml(str(graphml_path))
    with pickle_path.open("wb") as fh:
        pickle.dump(g, fh)
    summary = graph_summary(g)
    summary_path.write_text(json.dumps(summary, indent=2))
    return {
        "graphml": str(graphml_path),
        "pickle": str(pickle_path),
        "summary": str(summary_path),
        **summary,
    }


def load_graph(path: str | Path) -> ig.Graph:
    p = Path(path)
    if p.suffix == ".pkl":
        with p.open("rb") as fh:
            return pickle.load(fh)
    return ig.Graph.Read_GraphML(str(p))
