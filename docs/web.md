# Web UI + HTML reports

ClusterFlow ships two browser-friendly entry points so you don't have to leave the terminal-aversion zone:

1. **Streamlit app** — drag-and-drop your data, tweak thresholds with sliders, see the figures live.
2. **Self-contained HTML report** — `report.html` is emitted alongside the figures on every `clusterflow run`. Embed-by-base64 means it's a single file you can email, attach to a PR, or upload to a wiki.

## Streamlit app

### Install

```bash
pip install "clusterflow[web]"
```

### Launch

```bash
clusterflow web
# → Local URL:  http://localhost:8501
```

Or directly:

```bash
streamlit run -m clusterflow.web.streamlit_app
```

### What the UI does

- **Sidebar — data source.** Either upload your own SNP matrix + epi metadata + (optional) MLST, or pick one of the bundled samples (synthetic 28-isolate fixture or the real *K. pneumoniae* example from SNP2Cluster v0.5.4).
- **Sidebar — thresholds.** SNP cutoff, day cutoff, Leiden resolution, bootstrap replicate count.
- **Main area.** Six tabs after a run completes:
  - **Clusters** — one row per detected cluster (size, STs, wards, date range, index case, confidence)
  - **Figures** — all five static figures inline
  - **Methods** — per-method cluster counts and a histogram of cross-method agreement
  - **Centrality** — every isolate ranked by transmission risk score
  - **Bootstrap** — cluster stability summary + bar chart
  - **Downloads** — one-click download of `report.html` and a zip of the full output directory

Everything runs server-side in the same Python process — no browser-side computation. With the synthetic fixture, results appear in ~3 s; with the real 36-isolate fixture, ~3 s; for 1 k isolates expect ~30 s.

## HTML report

### Automatic emission

Every `clusterflow run --config ...` writes a `report.html` to the output directory (assuming `visualization.static: true` in the config — the default). It includes:

- KPI cards: # isolates, # clusters, # multi-isolate clusters, # bootstrap-stable
- Cluster table
- Per-method counts + cross-method agreement summary
- All five figures (PNG, base64-embedded)
- Top 20 isolates by transmission risk score
- Cluster stability summary

The file is fully self-contained (no external CSS, JS, or images), so you can:

- Email it as an attachment.
- Drop it into a Slack / Teams channel.
- Upload to a lab wiki without breaking links.
- Print it to PDF — there's a `@media print` style rule that keeps figures from breaking across pages.

### Re-rendering after the fact

Already ran the pipeline but want to regenerate the report (e.g., after a CSS tweak)?

```bash
clusterflow report ./results
```

It reads the artefacts in `./results/` (cluster_assignments.csv, centrality_scores.csv, etc.) and re-emits `report.html` without re-running clustering.

### Programmatic access

```python
from clusterflow.web import render_report
from clusterflow.pipeline import run_pipeline
from clusterflow.config import ClusterFlowConfig

cfg = ClusterFlowConfig.from_yaml("config.yaml")
result = run_pipeline(cfg)
render_report(result, cfg.output_dir)
```

`render_report` accepts a `PipelineResult` plus the output directory and returns the report path.

## Picking the right entry point

| You want to… | Use |
|---|---|
| Let a student / collaborator run their own analysis without a CLI | Streamlit app |
| Share results with a remote collaborator | HTML report |
| Drop a result into a manuscript / lab notebook | HTML report → print to PDF |
| Live outbreak surveillance with isolate streaming | [`clusterflow serve`](streaming.md) (Dash + FastAPI) |
| Batch process many samples in CI | Plain `clusterflow run` |
