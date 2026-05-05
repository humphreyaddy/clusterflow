"""Phase 9.1 — Typer command-line interface."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from clusterflow import __version__
from clusterflow.config import ClusterFlowConfig, ConfigError, diff_configs

app = typer.Typer(
    name="clusterflow",
    help="Transmission cluster detection from WGS + epi metadata.",
    no_args_is_help=True,
    add_completion=False,
)

console = Console()


def _setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )


@app.callback()
def _callback(
    version: bool = typer.Option(
        False, "--version", "-V", help="Print version and exit."
    ),
) -> None:
    if version:
        console.print(f"clusterflow {__version__}")
        raise typer.Exit()


@app.command()
def run(
    config: Path = typer.Option(..., "--config", "-c", help="YAML config path"),
    output: Path | None = typer.Option(None, "--output", "-o", help="Override output_dir"),
    no_viz: bool = typer.Option(False, "--no-viz", help="Skip static figures"),
    verbose: bool = typer.Option(False, "--verbose", "-v"),
) -> None:
    """Run the full batch pipeline."""
    _setup_logging(verbose)
    try:
        cfg = ClusterFlowConfig.from_yaml(config)
    except ConfigError as e:
        console.print(f"[bold red]config error:[/] {e}")
        raise typer.Exit(code=2) from e
    if output is not None:
        cfg = cfg.model_copy(update={"output_dir": Path(output).expanduser()})

    from clusterflow.pipeline import run_pipeline

    result = run_pipeline(cfg, run_viz=not no_viz)
    _print_summary(result)


@app.command()
def validate(
    config: Path = typer.Option(..., "--config", "-c", help="YAML config path"),
) -> None:
    """Dry-run: validate config and ingest inputs without clustering."""
    _setup_logging(False)
    try:
        cfg = ClusterFlowConfig.from_yaml(config)
    except ConfigError as e:
        console.print(f"[bold red]config error:[/] {e}")
        raise typer.Exit(code=2) from e
    from clusterflow.ingestion import IngestionPipeline

    ingest = IngestionPipeline().load(cfg)
    console.print(
        f"[green]ok[/]  isolates={len(ingest.isolates)}  "
        f"snp_matrix={ingest.snp.n}×{ingest.snp.n}"
    )


@app.command(name="init")
def init_config(
    output: Path = typer.Option(
        Path("config.yaml"), "--output", "-o", help="Where to write the template"
    ),
) -> None:
    """Write an annotated config template."""
    template = Path(__file__).resolve().parents[1] / "configs" / "example.yaml"
    if not template.exists():
        console.print(f"[red]template not found:[/] {template}")
        raise typer.Exit(code=1)
    Path(output).write_text(template.read_text())
    console.print(f"[green]wrote[/] {output}")


@app.command()
def compare(
    a: Path = typer.Argument(..., help="Config A"),
    b: Path = typer.Argument(..., help="Config B"),
) -> None:
    """Diff two YAML configs."""
    cfg_a = ClusterFlowConfig.from_yaml(a)
    cfg_b = ClusterFlowConfig.from_yaml(b)
    diff = diff_configs(cfg_a, cfg_b)
    if not diff:
        console.print("[green]configs are identical[/]")
        return
    table = Table(title="config diff")
    table.add_column("key")
    table.add_column(str(a))
    table.add_column(str(b))
    for k, (va, vb) in diff.items():
        table.add_row(k, str(va), str(vb))
    console.print(table)


@app.command()
def simulate(
    output: Path = typer.Option(..., "--output", "-o", help="Where to write fixtures"),
    n_isolates: int = typer.Option(28),
    n_clusters: int = typer.Option(4),
    seed: int = typer.Option(42),
) -> None:
    """Generate a synthetic K. pneumoniae outbreak dataset."""
    from clusterflow.testing import write_fixtures

    paths = write_fixtures(output, n_isolates=n_isolates, n_clusters=n_clusters, random_state=seed)
    for k, v in paths.items():
        console.print(f"[green]wrote[/] {k:>8}  {v}")


@app.command()
def benchmark(
    output: Path = typer.Option(Path("results/benchmark"), "--output", "-o"),
    sizes: str = typer.Option(
        "50,200,500,1000",
        help="Comma-separated dataset sizes to benchmark",
    ),
) -> None:
    """Run Phase 8.3 performance benchmark across dataset sizes."""
    _setup_logging(False)
    from clusterflow.testing import run_benchmark

    size_list = [int(s.strip()) for s in sizes.split(",") if s.strip()]
    df = run_benchmark(sizes=size_list, output_dir=output)
    table = Table(title="benchmark results")
    for col in df.columns:
        table.add_column(col)
    for _, row in df.iterrows():
        table.add_row(*[f"{v:.3f}" if isinstance(v, float) else str(v) for v in row])
    console.print(table)
    console.print(f"\n[green]wrote[/] {output}/performance_benchmark.csv")


@app.command()
def web(
    port: int = typer.Option(8501, help="Streamlit port"),
    host: str = typer.Option("localhost"),
    headless: bool = typer.Option(False, help="Run without auto-opening a browser"),
) -> None:
    """Launch the interactive Streamlit web UI."""
    _setup_logging(False)
    try:
        import streamlit  # noqa: F401
    except ImportError as e:
        console.print(
            "[red]streamlit not installed:[/] pip install 'clusterflow[web]'"
        )
        raise typer.Exit(code=3) from e

    import subprocess
    import sys

    app_path = Path(__file__).parent / "web" / "streamlit_app.py"
    cmd = [
        sys.executable,
        "-m",
        "streamlit",
        "run",
        str(app_path),
        "--server.port",
        str(port),
        "--server.address",
        host,
    ]
    if headless:
        cmd += ["--server.headless", "true"]
    console.print(f"[green]launching web UI at[/] http://{host}:{port}")
    subprocess.run(cmd, check=False)


@app.command()
def report(
    output: Path = typer.Argument(..., help="Existing output_dir from `clusterflow run`"),
) -> None:
    """Re-render report.html from a previous run's artefacts."""
    _setup_logging(False)
    import json
    from datetime import date

    from clusterflow.models import (
        BootstrapResult,
        CentralityScores,
        ClusterAssignment,
        Isolate,
        PipelineResult,
        SNPMatrix,
        TransmissionCluster,
    )
    from clusterflow.web import render_report

    # Reconstruct a minimal PipelineResult from saved artefacts
    summary_path = Path(output) / "pipeline_summary.json"
    if not summary_path.exists():
        console.print(f"[red]not a pipeline output dir:[/] {summary_path} missing")
        raise typer.Exit(2)
    summary = json.loads(summary_path.read_text())

    # Cluster assignments per method
    import pandas as pd

    assigns_csv = Path(output) / "clusters" / "cluster_assignments.csv"
    df = pd.read_csv(assigns_csv) if assigns_csv.exists() else pd.DataFrame()
    method_assigns: dict[str, ClusterAssignment] = {}
    consensus = ClusterAssignment(method="consensus", assignments={}, n_clusters=0)
    if not df.empty:
        for col in df.columns:
            if col in {"isolate_id", "agreement_score", "ambiguous"}:
                continue
            assignments = dict(zip(df["isolate_id"], df[col].astype(int)))
            ca = ClusterAssignment(
                method=col, assignments=assignments,
                n_clusters=len(set(assignments.values())),
            )
            if col == "consensus":
                if "agreement_score" in df.columns:
                    ca = ca.model_copy(
                        update={
                            "agreement_score": dict(zip(df["isolate_id"], df["agreement_score"])),
                            "ambiguous": dict(zip(df["isolate_id"], df.get("ambiguous", pd.Series([False] * len(df))).astype(bool))),
                        }
                    )
                consensus = ca
            else:
                method_assigns[col] = ca

    # Reconstruct centrality + bootstrap from CSVs
    centrality: list[CentralityScores] = []
    cent_csv = Path(output) / "analysis" / "centrality_scores.csv"
    if cent_csv.exists():
        for _, r in pd.read_csv(cent_csv).iterrows():
            centrality.append(CentralityScores(**r.to_dict()))

    bootstrap: list[BootstrapResult] = []
    boot_csv = Path(output) / "analysis" / "bootstrap_stability.csv"
    if boot_csv.exists():
        for _, r in pd.read_csv(boot_csv).iterrows():
            bootstrap.append(BootstrapResult(**r.to_dict()))

    # Cluster summaries from JSON
    transmission_clusters: list[TransmissionCluster] = []
    for c in summary.get("transmission_clusters", []):
        transmission_clusters.append(
            TransmissionCluster(
                cluster_id=c["cluster_id"],
                isolate_ids=[],  # not persisted in summary; OK for the report
                sequence_types=c.get("sequence_types", []),
                wards=c.get("wards", []),
                date_range=(date.fromisoformat(c["date_range"][0]), date.fromisoformat(c["date_range"][1])),
                index_case_candidate=c.get("index_case_candidate"),
                confidence=c.get("confidence", 0.0),
            )
        )

    # Minimal SNPMatrix + isolates so PipelineResult validates (only used in advanced sections)
    import numpy as np

    iso_ids = list(consensus.assignments) or list(df.get("isolate_id", []))
    n = max(len(iso_ids), 1)
    snp = SNPMatrix(isolate_ids=iso_ids or ["x"], distances=np.zeros((max(len(iso_ids), 1), max(len(iso_ids), 1))))
    isolates = {
        i: Isolate(isolate_id=i, collection_date=date(2000, 1, 1), facility="?", ward="?")
        for i in iso_ids
    }

    result = PipelineResult(
        project_name=summary.get("project_name", "report"),
        n_isolates=summary.get("n_isolates", n),
        isolates=isolates,
        snp_matrix=snp,
        cluster_assignments=method_assigns,
        consensus=consensus,
        transmission_clusters=transmission_clusters,
        centrality=centrality,
        bootstrap=bootstrap,
    )
    path = render_report(result, Path(output))
    console.print(f"[green]wrote[/] {path}")


@app.command()
def serve(
    config: Path = typer.Option(..., "--config", "-c"),
    port: int = typer.Option(8000, help="API port"),
    host: str = typer.Option("0.0.0.0"),
) -> None:
    """Start the streaming/dashboard server (Phase 7)."""
    _setup_logging(False)
    try:
        from clusterflow.streaming.api import create_app
    except ImportError as e:
        console.print(
            "[red]streaming requires fastapi+uvicorn:[/] pip install fastapi uvicorn"
        )
        raise typer.Exit(code=3) from e
    cfg = ClusterFlowConfig.from_yaml(config)
    app_obj = create_app(cfg)

    try:
        import uvicorn
    except ImportError:
        console.print("[red]uvicorn not installed[/]")
        raise typer.Exit(code=3)
    uvicorn.run(app_obj, host=host, port=port)


def _print_summary(result) -> None:
    table = Table(title=f"{result.project_name} — pipeline summary")
    table.add_column("cluster", justify="right")
    table.add_column("n", justify="right")
    table.add_column("STs")
    table.add_column("wards")
    table.add_column("date range")
    table.add_column("index case")
    table.add_column("conf")
    for c in result.transmission_clusters:
        table.add_row(
            str(c.cluster_id),
            str(len(c.isolate_ids)),
            ",".join(c.sequence_types[:3]) or "-",
            ",".join(c.wards[:3]),
            f"{c.date_range[0]} → {c.date_range[1]}",
            c.index_case_candidate or "-",
            f"{c.confidence:.2f}",
        )
    console.print(table)


def main() -> None:
    """Console entry point."""
    sys.exit(app())


if __name__ == "__main__":
    main()
