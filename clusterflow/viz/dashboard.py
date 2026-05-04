"""Phase 6.2 — interactive Dash dashboard.

Importing this module raises ImportError if Dash + dash-cytoscape are not
installed. Use ``pip install dash dash-cytoscape plotly`` to enable.
"""

from __future__ import annotations

from typing import Any

import dash
import dash_cytoscape as cyto
import plotly.express as px
from dash import Input, Output, dcc, html

from clusterflow.models import PipelineResult


def build_dashboard(result: PipelineResult) -> dash.Dash:
    """Return a configured Dash app rendering the four-panel layout."""
    app = dash.Dash(__name__)
    consensus = result.consensus

    elements: list[dict[str, Any]] = []
    cluster_palette = px.colors.qualitative.Alphabet
    for iso, cid in consensus.assignments.items():
        elements.append(
            {
                "data": {
                    "id": iso,
                    "label": iso,
                    "cluster": str(cid),
                    "ward": result.isolates[iso].ward,
                },
                "classes": f"cluster_{cid % len(cluster_palette)}",
            }
        )

    timeline_df = [
        {
            "isolate_id": iso,
            "cluster": str(cid),
            "date": result.isolates[iso].collection_date,
            "ward": result.isolates[iso].ward,
        }
        for iso, cid in consensus.assignments.items()
    ]
    fig = px.scatter(
        timeline_df,
        x="date",
        y="cluster",
        color="cluster",
        hover_data=["isolate_id", "ward"],
        title="epi-genomic timeline",
    )

    app.layout = html.Div(
        [
            html.H1(f"ClusterFlow — {result.project_name}"),
            html.Div(
                [
                    cyto.Cytoscape(
                        id="network",
                        elements=elements,
                        style={"width": "100%", "height": "60vh"},
                        layout={"name": "cose"},
                    )
                ],
                style={"width": "55%", "display": "inline-block"},
            ),
            html.Div(
                [
                    dcc.Graph(figure=fig, id="timeline"),
                    html.Div(id="cluster-table"),
                    html.Div(id="isolate-detail"),
                ],
                style={"width": "44%", "display": "inline-block", "verticalAlign": "top"},
            ),
        ]
    )

    @app.callback(Output("isolate-detail", "children"), Input("network", "tapNodeData"))
    def show_detail(data):  # type: ignore[no-redef]
        if not data:
            return "Click a node to see its details."
        iso = data["id"]
        info = result.isolates.get(iso)
        if not info:
            return f"Unknown isolate {iso}"
        return html.Pre(info.model_dump_json(indent=2))

    return app
