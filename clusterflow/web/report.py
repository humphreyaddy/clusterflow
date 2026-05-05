"""Self-contained HTML report generator.

Produces a single ``report.html`` per pipeline run, with all figures embedded
as base64 PNGs so the file is shareable without any external assets.
"""

from __future__ import annotations

import base64
import logging
from datetime import datetime
from html import escape
from pathlib import Path
from typing import Iterable

import pandas as pd

from clusterflow.models import PipelineResult

log = logging.getLogger(__name__)


_FIGURES = [
    ("snp_heatmap", "Pairwise SNP distance heatmap",
     "Yellow blocks on the diagonal mark each cluster. Within-cluster distances are tight (≤ snp_cutoff); between-cluster distances are large."),
    ("minimum_spanning_tree", "Minimum spanning tree",
     "Each connected component is a cluster. Node colour = cluster, marker = ward, size scales with transmission risk score. Red stars mark candidate index cases."),
    ("epi_timeline_scatter", "Epi-genomic timeline",
     "Y-axis = cluster ID, X-axis = collection date. Black-edged circles are index case candidates. Cluster bands shade their full date range."),
    ("cluster_comparison_grid", "Method comparison",
     "Same network laid out four ways, coloured by each method's cluster assignment. Compares SNP-chains vs. Leiden vs. Spectral against the consensus."),
    ("bootstrap_stability_terrain", "Bootstrap stability terrain",
     "MDS projection of isolates over a KDE density of bootstrap stability. Warm peaks = many stable isolates clustered together."),
]


def _encode_png(path: Path) -> str | None:
    if not path.exists():
        return None
    return base64.b64encode(path.read_bytes()).decode("ascii")


def _df_to_html_table(df: pd.DataFrame, max_rows: int = 200) -> str:
    if df is None or df.empty:
        return "<p><em>(no data)</em></p>"
    truncated = len(df) > max_rows
    df = df.head(max_rows)
    rows = []
    rows.append("<thead><tr>" + "".join(f"<th>{escape(str(c))}</th>" for c in df.columns) + "</tr></thead>")
    body: list[str] = []
    for _, row in df.iterrows():
        cells = []
        for v in row:
            if isinstance(v, float):
                cells.append(f"<td>{v:.3g}</td>")
            elif pd.isna(v):
                cells.append("<td><em>—</em></td>")
            else:
                cells.append(f"<td>{escape(str(v))}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    rows.append("<tbody>" + "".join(body) + "</tbody>")
    table = '<table class="data">' + "".join(rows) + "</table>"
    if truncated:
        table += f'<p class="note">Showing first {max_rows} of {len(df)} rows.</p>'
    return table


def _cluster_summary_df(result: PipelineResult) -> pd.DataFrame:
    rows = []
    for c in result.transmission_clusters:
        rows.append(
            {
                "cluster_id": c.cluster_id,
                "n_isolates": len(c.isolate_ids),
                "sequence_types": ", ".join(c.sequence_types) or "—",
                "wards": ", ".join(c.wards),
                "date_start": c.date_range[0].isoformat(),
                "date_end": c.date_range[1].isoformat(),
                "index_case": c.index_case_candidate or "—",
                "confidence": round(c.confidence, 3),
            }
        )
    return pd.DataFrame(rows)


def _figures_html(figures_dir: Path) -> str:
    blocks: list[str] = []
    for stem, title, caption in _FIGURES:
        b64 = _encode_png(figures_dir / f"{stem}.png")
        if b64 is None:
            continue
        blocks.append(
            f"""
            <figure id="fig-{stem}">
              <h3>{escape(title)}</h3>
              <img src="data:image/png;base64,{b64}" alt="{escape(title)}" />
              <figcaption>{escape(caption)}</figcaption>
            </figure>
            """
        )
    return "\n".join(blocks)


def _per_method_counts_html(result: PipelineResult) -> str:
    rows = []
    rows.append("<thead><tr><th>method</th><th>n clusters</th></tr></thead>")
    body = []
    for m, a in result.cluster_assignments.items():
        body.append(f"<tr><td>{escape(m)}</td><td>{a.n_clusters}</td></tr>")
    body.append(f"<tr><td><strong>consensus</strong></td><td><strong>{result.consensus.n_clusters}</strong></td></tr>")
    rows.append("<tbody>" + "".join(body) + "</tbody>")
    return '<table class="data">' + "".join(rows) + "</table>"


def _agreement_summary(result: PipelineResult) -> str:
    if not result.consensus.agreement_score:
        return "<p><em>(no agreement data)</em></p>"
    s = pd.Series(result.consensus.agreement_score)
    n_full = int((s == 1.0).sum())
    n_amb = sum(1 for v in (result.consensus.ambiguous or {}).values() if v)
    return (
        f"<p><strong>{n_full}</strong> of {len(s)} isolates fully agreed across all methods. "
        f"<strong>{n_amb}</strong> isolates flagged ambiguous (&lt; 2/3 agreement).</p>"
    )


def _kpi_cards(result: PipelineResult) -> str:
    n_clusters = result.consensus.n_clusters
    n_isolates = result.n_isolates
    multi_isolate = sum(1 for c in result.transmission_clusters if len(c.isolate_ids) > 1)
    boot_stable = sum(1 for b in result.bootstrap if b.classification == "stable")
    return f"""
    <div class="kpi-grid">
      <div class="kpi"><div class="kpi-num">{n_isolates}</div><div class="kpi-label">isolates</div></div>
      <div class="kpi"><div class="kpi-num">{n_clusters}</div><div class="kpi-label">clusters (consensus)</div></div>
      <div class="kpi"><div class="kpi-num">{multi_isolate}</div><div class="kpi-label">multi-isolate clusters</div></div>
      <div class="kpi"><div class="kpi-num">{boot_stable}/{len(result.bootstrap) or 0}</div><div class="kpi-label">bootstrap-stable isolates</div></div>
    </div>
    """


_CSS = """
:root {
  --fg: #1a1a1a;
  --fg-muted: #555;
  --bg: #fdfdfd;
  --accent: #0d7377;
  --accent-soft: #e6f4f5;
  --border: #e1e1e1;
  --warn: #c45a00;
}
* { box-sizing: border-box; }
body {
  font: 15px/1.55 -apple-system, BlinkMacSystemFont, "Segoe UI", "Helvetica Neue", Arial, sans-serif;
  margin: 0;
  background: var(--bg);
  color: var(--fg);
}
header {
  background: linear-gradient(135deg, var(--accent) 0%, #14a098 100%);
  color: white;
  padding: 2rem 2.5rem;
}
header h1 { margin: 0; font-size: 1.8rem; font-weight: 600; }
header .meta { opacity: 0.85; font-size: 0.9rem; margin-top: 0.25rem; }
main { max-width: 1200px; margin: 0 auto; padding: 2rem 2.5rem; }
h2 {
  font-size: 1.3rem;
  margin: 2.5rem 0 1rem;
  padding-bottom: 0.4rem;
  border-bottom: 1px solid var(--border);
  color: var(--accent);
}
h3 { margin: 1.5rem 0 0.5rem; font-size: 1.05rem; }
p { color: var(--fg-muted); }
.kpi-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
  gap: 1rem;
  margin: 1.5rem 0;
}
.kpi {
  background: var(--accent-soft);
  border-left: 4px solid var(--accent);
  padding: 1rem 1.25rem;
  border-radius: 4px;
}
.kpi-num { font-size: 2rem; font-weight: 700; color: var(--accent); }
.kpi-label { font-size: 0.85rem; color: var(--fg-muted); text-transform: uppercase; letter-spacing: 0.04em; }
table.data {
  border-collapse: collapse;
  width: 100%;
  margin: 0.5rem 0 1rem;
  font-size: 0.9rem;
}
table.data th, table.data td {
  text-align: left;
  padding: 0.45rem 0.7rem;
  border-bottom: 1px solid var(--border);
}
table.data th {
  background: #f5f5f5;
  font-weight: 600;
  color: var(--fg);
}
table.data tbody tr:hover { background: #fafafa; }
.note { color: var(--fg-muted); font-size: 0.85rem; }
figure {
  margin: 1.5rem 0;
  padding: 1rem;
  background: white;
  border: 1px solid var(--border);
  border-radius: 6px;
}
figure img { max-width: 100%; height: auto; display: block; margin: 0 auto; }
figcaption { color: var(--fg-muted); font-size: 0.9rem; margin-top: 0.5rem; text-align: center; }
.toc {
  background: #f9f9f9;
  border: 1px solid var(--border);
  padding: 1rem 1.25rem;
  border-radius: 4px;
  margin: 1.5rem 0;
}
.toc ul { margin: 0.4rem 0 0; padding-left: 1.2rem; }
.toc a { color: var(--accent); text-decoration: none; }
.toc a:hover { text-decoration: underline; }
footer {
  border-top: 1px solid var(--border);
  margin-top: 4rem;
  padding: 1.5rem 2.5rem;
  text-align: center;
  color: var(--fg-muted);
  font-size: 0.85rem;
}
@media print {
  header { background: var(--accent) !important; -webkit-print-color-adjust: exact; }
  figure { page-break-inside: avoid; }
}
"""


def render_report(
    result: PipelineResult,
    output_dir: Path,
    figures_dir: Path | None = None,
    analysis_dir: Path | None = None,
) -> Path:
    """Render a self-contained ``report.html`` for one pipeline result.

    Parameters
    ----------
    result : PipelineResult
        Output of :func:`clusterflow.pipeline.run_pipeline`.
    output_dir : Path
        Directory where ``report.html`` is written.
    figures_dir : Path, optional
        Directory containing the five PNG figures. Defaults to
        ``output_dir / "figures"``.
    analysis_dir : Path, optional
        Directory containing centrality/bootstrap CSVs. Defaults to
        ``output_dir / "analysis"``.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    figs = figures_dir or (out / "figures")
    analysis = analysis_dir or (out / "analysis")

    cluster_df = _cluster_summary_df(result)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    centrality_html = "<p><em>(no centrality data)</em></p>"
    cent_csv = analysis / "centrality_scores.csv"
    if cent_csv.exists():
        cent_df = pd.read_csv(cent_csv).sort_values(
            "transmission_risk_score", ascending=False
        )
        centrality_html = _df_to_html_table(cent_df, max_rows=20)

    bootstrap_html = "<p><em>(no bootstrap data)</em></p>"
    boot_csv = analysis / "cluster_stability_summary.csv"
    if boot_csv.exists():
        boot_df = pd.read_csv(boot_csv)
        bootstrap_html = _df_to_html_table(boot_df, max_rows=200)

    figures_html = _figures_html(figs)
    cluster_table_html = _df_to_html_table(cluster_df)
    methods_html = _per_method_counts_html(result)
    agreement_html = _agreement_summary(result)
    kpi_html = _kpi_cards(result)

    body = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>ClusterFlow report — {escape(result.project_name)}</title>
  <style>{_CSS}</style>
</head>
<body>
<header>
  <h1>{escape(result.project_name)}</h1>
  <div class="meta">ClusterFlow report · generated {timestamp}</div>
</header>
<main>
  <section id="overview">
    <h2>Overview</h2>
    {kpi_html}
    {agreement_html}
  </section>

  <nav class="toc">
    <strong>Contents</strong>
    <ul>
      <li><a href="#clusters">Clusters</a></li>
      <li><a href="#methods">Per-method counts</a></li>
      <li><a href="#figures">Figures</a></li>
      <li><a href="#centrality">Top isolates by transmission risk</a></li>
      <li><a href="#bootstrap">Bootstrap cluster stability</a></li>
    </ul>
  </nav>

  <section id="clusters">
    <h2>Clusters</h2>
    {cluster_table_html}
  </section>

  <section id="methods">
    <h2>Per-method cluster counts</h2>
    {methods_html}
  </section>

  <section id="figures">
    <h2>Figures</h2>
    {figures_html}
  </section>

  <section id="centrality">
    <h2>Top isolates by transmission risk score</h2>
    <p>One row per isolate, ranked by composite risk score (50 % betweenness · 30 % out-degree · 20 % closeness, normalised within cluster).</p>
    {centrality_html}
  </section>

  <section id="bootstrap">
    <h2>Bootstrap cluster stability</h2>
    <p>Mean stability across {len(result.bootstrap)} bootstrap replicates; grade A ≥ 0.90, B ≥ 0.75, C otherwise.</p>
    {bootstrap_html}
  </section>
</main>
<footer>
  Generated by <a href="https://github.com/humphreyaddy/clusterflow">ClusterFlow</a>.
</footer>
</body>
</html>
"""
    path = out / "report.html"
    path.write_text(body, encoding="utf-8")
    log.info("HTML report written: %s", path)
    return path


def write_pipeline_report(result: PipelineResult, output_dir: Path) -> Path | None:
    """Convenience wrapper called from the pipeline orchestrator."""
    try:
        return render_report(result, Path(output_dir))
    except Exception:
        log.exception("failed to render HTML report")
        return None
