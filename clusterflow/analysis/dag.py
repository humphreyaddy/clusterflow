"""Phase 5.1 — temporal DAG construction.

For each undirected edge ``(A, B)``:
  * ``date_A < date_B - uncertainty`` → directed A → B
  * ``date_B < date_A - uncertainty`` → directed B → A
  * within ``±uncertainty`` of each other → bidirectional (both directions retained)
Cycles are detected and broken on the lowest-weight edge.
"""

from __future__ import annotations

import json
import logging
from datetime import date
from pathlib import Path

import igraph as ig

log = logging.getLogger(__name__)


def build_temporal_dag(
    g: ig.Graph,
    uncertainty_days: int = 2,
) -> ig.Graph:
    """Return a directed igraph reflecting probable transmission direction."""
    n = g.vcount()
    dag = ig.Graph(n=n, directed=True)
    for attr in g.vs.attributes():
        dag.vs[attr] = list(g.vs[attr])

    dates: list[date] = [date.fromisoformat(d) for d in g.vs["collection_date"]]

    edges: list[tuple[int, int]] = []
    edge_attrs: dict[str, list] = {k: [] for k in g.es.attributes()}
    for e in g.es:
        i, j = e.tuple
        delta = (dates[j] - dates[i]).days
        directions: list[tuple[int, int]] = []
        if delta > uncertainty_days:
            directions = [(i, j)]
        elif -delta > uncertainty_days:
            directions = [(j, i)]
        else:
            directions = [(i, j), (j, i)]
        for s, t in directions:
            edges.append((s, t))
            for k in edge_attrs:
                edge_attrs[k].append(e[k])

    dag.add_edges(edges)
    for k, vs in edge_attrs.items():
        dag.es[k] = vs

    _break_cycles(dag)
    return dag


def _break_cycles(dag: ig.Graph) -> None:
    """Remove the lowest-similarity edge in any cycle until DAG is acyclic.

    Cycles can arise from bidirectional ambiguous edges. We delete edges with
    the highest ``composite_weight`` (= lowest similarity) first.
    """
    if dag.is_dag():
        return
    iterations = 0
    while not dag.is_dag() and dag.ecount() > 0 and iterations < 10000:
        try:
            cycle = next(_find_cycle(dag))
        except StopIteration:
            break
        # Pick the edge with the highest composite_weight in this cycle
        eids = []
        for u, v in zip(cycle, cycle[1:] + [cycle[0]]):
            try:
                eids.append(dag.get_eid(u, v))
            except Exception:
                continue
        if not eids:
            break
        weights = [dag.es[e]["composite_weight"] for e in eids]
        worst = eids[int(max(range(len(weights)), key=lambda i: weights[i]))]
        dag.delete_edges([worst])
        iterations += 1
    if not dag.is_dag():
        log.warning("could not break all cycles in DAG (data quality issue)")


def _find_cycle(dag: ig.Graph):
    """Yield one simple cycle (as a list of vertex ids) if any exists."""
    color = ["white"] * dag.vcount()
    parent: list[int] = [-1] * dag.vcount()

    def dfs(u: int):
        stack = [(u, iter(dag.successors(u)))]
        color[u] = "gray"
        while stack:
            node, it = stack[-1]
            try:
                nxt = next(it)
            except StopIteration:
                color[node] = "black"
                stack.pop()
                continue
            if color[nxt] == "gray":
                cycle = [nxt, node]
                cur = node
                while parent[cur] != -1 and parent[cur] != nxt:
                    cur = parent[cur]
                    cycle.append(cur)
                cycle.reverse()
                yield cycle
                return
            if color[nxt] == "white":
                parent[nxt] = node
                color[nxt] = "gray"
                stack.append((nxt, iter(dag.successors(nxt))))

    for v in range(dag.vcount()):
        if color[v] == "white":
            yield from dfs(v)


def save_dag(dag: ig.Graph, output_dir: str | Path) -> dict:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    graphml = out / "transmission_dag.graphml"
    dag.write_graphml(str(graphml))

    # Source/sink/chain stats
    n_source = sum(1 for v in dag.vs if dag.indegree(v) == 0)
    n_sink = sum(1 for v in dag.vs if dag.outdegree(v) == 0)
    components = dag.connected_components(mode="weak")
    longest_chain = 0
    if dag.is_dag() and dag.vcount():
        try:
            longest_chain = max(_longest_path_length(dag))
        except Exception:
            longest_chain = 0
    summary = {
        "n_vertices": dag.vcount(),
        "n_edges": dag.ecount(),
        "n_source_nodes": n_source,
        "n_sink_nodes": n_sink,
        "n_weakly_connected_components": len(components),
        "max_chain_length": longest_chain,
        "is_dag": dag.is_dag(),
    }
    (out / "dag_summary.json").write_text(json.dumps(summary, indent=2))
    return summary


def _longest_path_length(dag: ig.Graph):
    """DP over topological order; yields longest path length per vertex."""
    order = dag.topological_sorting()
    longest = [0] * dag.vcount()
    for u in order:
        for v in dag.successors(u):
            longest[v] = max(longest[v], longest[u] + 1)
    return longest
