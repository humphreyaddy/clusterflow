"""Streamlit web UI for ClusterFlow.

Launch with::

    streamlit run -m clusterflow.web.streamlit_app

or via the CLI helper::

    clusterflow web

The user uploads a SNP matrix, epi metadata, and (optional) MLST profiles,
tweaks thresholds, hits "Run". The pipeline runs server-side; results are
shown in tabs (summary, figures, tables, HTML report download).
"""

from __future__ import annotations

import io
import logging
import shutil
import tempfile
import time
import zipfile
from pathlib import Path
from typing import Any

import pandas as pd
import streamlit as st

from clusterflow.config import (
    ClusteringConfig,
    ClusterFlowConfig,
    DistanceEngineConfig,
    InputsConfig,
    StreamingConfig,
    ThresholdsConfig,
    VisualizationConfig,
)
from clusterflow.pipeline import run_pipeline

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

REPO_ROOT = Path(__file__).resolve().parents[2]
SAMPLE_FIXTURE = REPO_ROOT / "tests" / "fixtures"
SAMPLE_REAL = SAMPLE_FIXTURE / "kp_real"


# ---------- helpers ----------


def _save_upload(upload, target: Path) -> Path:
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(upload.getvalue())
    return target


def _load_sample(name: str) -> tuple[Path, Path, Path | None]:
    """Return (snp, epi, mlst|None) paths for a bundled sample."""
    if name == "synthetic 28-isolate":
        return (
            SAMPLE_FIXTURE / "kp_snp_matrix.tsv",
            SAMPLE_FIXTURE / "kp_epi.csv",
            SAMPLE_FIXTURE / "kp_mlst.csv",
        )
    if name == "real K. pneumoniae (SNP2Cluster example, 36 isolates)":
        return (
            SAMPLE_REAL / "snp_matrix.csv",
            SAMPLE_REAL / "epi.csv",
            SAMPLE_REAL / "mlst.csv",
        )
    raise ValueError(name)


def _zip_results(out_dir: Path) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in out_dir.rglob("*"):
            if p.is_file():
                zf.write(p, p.relative_to(out_dir))
    return buf.getvalue()


def _run(
    snp_path: Path,
    epi_path: Path,
    mlst_path: Path | None,
    out_dir: Path,
    snp_cutoff: int,
    day_cutoff: int,
    bootstrap_n: int,
    leiden_resolution: float,
    methods: list[str],
    project_name: str,
):
    cfg = ClusterFlowConfig(
        project_name=project_name,
        output_dir=out_dir,
        inputs=InputsConfig(
            snp_matrix=snp_path,
            epi_metadata=epi_path,
            mlst_profiles=mlst_path,
        ),
        thresholds=ThresholdsConfig(
            snp_cutoff=snp_cutoff,
            day_cutoff=day_cutoff,
            edge_weight_alpha=0.4,
            edge_weight_beta=0.4,
            edge_weight_gamma=0.2,
        ),
        clustering=ClusteringConfig(
            methods=list(methods),
            bootstrap_n=bootstrap_n,
            leiden_resolution=leiden_resolution,
            random_state=42,
        ),
        distance_engine=DistanceEngineConfig(backend="cpu", n_jobs=2),
        streaming=StreamingConfig(),
        visualization=VisualizationConfig(static=True, dashboard=False),
    )
    cfg.validate_paths()
    return run_pipeline(cfg)


# ---------- layout ----------


def main() -> None:
    st.set_page_config(
        page_title="ClusterFlow",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded",
    )

    st.title("ClusterFlow")
    st.caption(
        "Transmission cluster detection from WGS + epi metadata. "
        "Upload your data on the left, tweak thresholds, hit run."
    )

    # ---- sidebar: inputs ----
    with st.sidebar:
        st.header("1. Data source")
        sources = ["upload my own"]
        if SAMPLE_FIXTURE.exists():
            sources.append("synthetic 28-isolate")
        if SAMPLE_REAL.exists():
            sources.append("real K. pneumoniae (SNP2Cluster example, 36 isolates)")
        source = st.radio("Data source", sources, label_visibility="collapsed")

        snp_path: Path | None = None
        epi_path: Path | None = None
        mlst_path: Path | None = None
        upload_dir: Path | None = None

        if source == "upload my own":
            snp_up = st.file_uploader(
                "SNP distance matrix (CSV/TSV)",
                type=["csv", "tsv", "txt"],
                help="Square symmetric matrix; first column = isolate IDs; column headers must match.",
            )
            epi_up = st.file_uploader(
                "Epi metadata (CSV/TSV)",
                type=["csv", "tsv", "txt"],
                help="Required columns: isolate_id, collection_date, facility, ward.",
            )
            mlst_up = st.file_uploader(
                "MLST profiles (CSV/TSV) — optional",
                type=["csv", "tsv", "txt"],
                help="Required columns: isolate_id, ST.",
            )
            if snp_up and epi_up:
                upload_dir = Path(tempfile.mkdtemp(prefix="clusterflow_upload_"))
                snp_path = _save_upload(snp_up, upload_dir / snp_up.name)
                epi_path = _save_upload(epi_up, upload_dir / epi_up.name)
                if mlst_up is not None:
                    mlst_path = _save_upload(mlst_up, upload_dir / mlst_up.name)
        else:
            snp_path, epi_path, mlst_path = _load_sample(source)
            with st.expander("preview sample inputs"):
                st.write(f"**SNP matrix:** `{snp_path.name}`")
                st.write(f"**Epi metadata:** `{epi_path.name}`")
                if mlst_path:
                    st.write(f"**MLST:** `{mlst_path.name}`")

        st.divider()
        st.header("2. Thresholds")
        snp_cutoff = st.slider("SNP cutoff", 1, 100, 20, help="Max SNP distance for an edge.")
        day_cutoff = st.slider("Day cutoff", 1, 365, 45, help="Max collection-date gap.")

        st.header("3. Clustering")
        methods = st.multiselect(
            "Methods",
            options=["snp_chains", "leiden", "spectral"],
            default=["snp_chains", "leiden", "spectral"],
        )
        leiden_res = st.slider("Leiden resolution", 0.1, 3.0, 1.0, step=0.1)
        bootstrap_n = st.slider("Bootstrap replicates", 0, 500, 100, step=10)

        st.header("4. Project")
        project_name = st.text_input("Project name", value="clusterflow_run")

        run_clicked = st.button("Run pipeline", type="primary", use_container_width=True)

    # ---- main panel ----
    if not run_clicked:
        st.info(
            "Pick a data source, tweak the thresholds on the left, and click "
            "**Run pipeline**. Results will appear here."
        )
        with st.expander("What does ClusterFlow do?"):
            st.markdown(
                """
                Each run:

                1. Builds a **weighted transmission graph** from your SNP matrix + epi metadata + MLST.
                2. Runs **three clustering algorithms in parallel** — SNP chains (SNP2Cluster-compatible),
                   Leiden community detection, and spectral clustering with eigengap k-selection.
                3. Computes a **Hungarian-aligned consensus** across the three.
                4. Builds a **temporal DAG** with cycle-breaking and identifies **index-case candidates**
                   via centrality scoring.
                5. Runs a **bootstrap stability analysis** to give every cluster a confidence grade A/B/C.
                6. Produces five **publication-ready figures** plus a **self-contained HTML report**.
                """
            )
        return

    if not snp_path or not epi_path:
        st.error("Need at least an SNP matrix and an epi-metadata file.")
        return
    if not methods:
        st.error("Select at least one clustering method.")
        return

    out_dir = Path(tempfile.mkdtemp(prefix="clusterflow_run_"))
    progress = st.progress(0.0, text="initialising")

    t0 = time.perf_counter()
    try:
        progress.progress(0.1, text="ingesting inputs")
        result = _run(
            snp_path=snp_path,
            epi_path=epi_path,
            mlst_path=mlst_path,
            out_dir=out_dir,
            snp_cutoff=snp_cutoff,
            day_cutoff=day_cutoff,
            bootstrap_n=bootstrap_n,
            leiden_resolution=leiden_res,
            methods=methods,
            project_name=project_name,
        )
        progress.progress(1.0, text="done")
    except Exception as e:  # noqa: BLE001
        progress.empty()
        st.error(f"Pipeline failed: {e}")
        st.exception(e)
        return
    elapsed = time.perf_counter() - t0

    st.success(
        f"Pipeline finished in {elapsed:.1f} s — "
        f"{result.consensus.n_clusters} clusters across {result.n_isolates} isolates."
    )

    # KPI cards
    cols = st.columns(4)
    cols[0].metric("isolates", result.n_isolates)
    cols[1].metric("clusters (consensus)", result.consensus.n_clusters)
    multi = sum(1 for c in result.transmission_clusters if len(c.isolate_ids) > 1)
    cols[2].metric("multi-isolate clusters", multi)
    stable = sum(1 for b in result.bootstrap if b.classification == "stable")
    cols[3].metric(
        "bootstrap-stable",
        f"{stable}/{len(result.bootstrap) or 0}",
    )

    tabs = st.tabs(
        ["Clusters", "Figures", "Methods", "Centrality", "Bootstrap", "Downloads"]
    )

    # ---- Clusters tab ----
    with tabs[0]:
        cluster_rows = [
            {
                "cluster_id": c.cluster_id,
                "n_isolates": len(c.isolate_ids),
                "sequence_types": ", ".join(c.sequence_types) or "—",
                "wards": ", ".join(c.wards),
                "date_start": c.date_range[0],
                "date_end": c.date_range[1],
                "index_case": c.index_case_candidate,
                "confidence": c.confidence,
            }
            for c in result.transmission_clusters
        ]
        st.dataframe(pd.DataFrame(cluster_rows), use_container_width=True, hide_index=True)

    # ---- Figures tab ----
    with tabs[1]:
        figs_dir = out_dir / "figures"
        for stem, title in [
            ("snp_heatmap", "Pairwise SNP distance heatmap"),
            ("minimum_spanning_tree", "Minimum spanning tree"),
            ("epi_timeline_scatter", "Epi-genomic timeline"),
            ("cluster_comparison_grid", "Method comparison grid"),
            ("bootstrap_stability_terrain", "Bootstrap stability terrain"),
        ]:
            png = figs_dir / f"{stem}.png"
            if png.exists():
                st.subheader(title)
                st.image(str(png), use_container_width=True)

    # ---- Methods tab ----
    with tabs[2]:
        m_rows = [
            {"method": m, "n_clusters": a.n_clusters}
            for m, a in result.cluster_assignments.items()
        ]
        m_rows.append({"method": "consensus", "n_clusters": result.consensus.n_clusters})
        st.dataframe(pd.DataFrame(m_rows), use_container_width=True, hide_index=True)
        if result.consensus.agreement_score:
            agreement = pd.Series(result.consensus.agreement_score, name="agreement_score")
            st.write(
                f"**{int((agreement == 1.0).sum())}** of {len(agreement)} "
                "isolates fully agreed across all methods."
            )
            st.bar_chart(agreement.value_counts().sort_index())

    # ---- Centrality tab ----
    with tabs[3]:
        cent_csv = out_dir / "analysis" / "centrality_scores.csv"
        if cent_csv.exists():
            df = pd.read_csv(cent_csv).sort_values(
                "transmission_risk_score", ascending=False
            )
            st.dataframe(df, use_container_width=True, hide_index=True)
        else:
            st.info("No centrality scores produced.")

    # ---- Bootstrap tab ----
    with tabs[4]:
        boot_csv = out_dir / "analysis" / "cluster_stability_summary.csv"
        if boot_csv.exists():
            df = pd.read_csv(boot_csv)
            st.dataframe(df, use_container_width=True, hide_index=True)
            st.bar_chart(df.set_index("cluster_id")["mean_stability"])
        else:
            st.info("No bootstrap data — set replicates > 0 in the sidebar.")

    # ---- Downloads tab ----
    with tabs[5]:
        st.subheader("Downloads")
        report_path = out_dir / "report.html"
        if report_path.exists():
            st.download_button(
                "📄 HTML report (self-contained)",
                report_path.read_bytes(),
                file_name=f"{project_name}_report.html",
                mime="text/html",
                use_container_width=True,
            )
        st.download_button(
            "🗜️  All results (zip)",
            _zip_results(out_dir),
            file_name=f"{project_name}_results.zip",
            mime="application/zip",
            use_container_width=True,
        )
        st.caption(f"Run output directory: `{out_dir}`")


if __name__ == "__main__":
    main()
