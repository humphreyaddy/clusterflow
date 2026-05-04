"""Phase 2.3 — pairsnp subprocess backend.

pairsnp expects a single multi-FASTA file and emits a sparse TSV of pairwise
distances. We write a temp combined FASTA, run pairsnp, then materialise the
result into a dense :class:`SNPMatrix`.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from clusterflow.distance.backend import DistanceBackend
from clusterflow.distance.cpu import _read_fasta
from clusterflow.models import SNPMatrix

log = logging.getLogger(__name__)


class PairSNPBackend(DistanceBackend):
    name = "pairsnp"

    def __init__(self, binary: str = "pairsnp") -> None:
        self.binary = shutil.which(binary) or binary

    def compute(self, fasta_paths: list[Path], n_jobs: int = -1) -> SNPMatrix:
        if not shutil.which(self.binary):
            raise RuntimeError(f"pairsnp binary not on PATH: {self.binary}")

        records = [_read_fasta(p) for p in fasta_paths]
        ids = [r[0] for r in records]

        with tempfile.TemporaryDirectory() as tmp:
            combined = Path(tmp) / "combined.fasta"
            with combined.open("w") as fh:
                for iso_id, seq in records:
                    fh.write(f">{iso_id}\n{seq}\n")
            cmd = [self.binary, "-c", str(combined)]
            log.info("running: %s", " ".join(cmd))
            res = subprocess.run(cmd, capture_output=True, text=True, check=False)
            if res.returncode != 0:
                raise RuntimeError(
                    f"pairsnp failed (exit {res.returncode}): {res.stderr.strip()}"
                )
            tsv = res.stdout
        df = pd.read_csv(
            pd.io.common.StringIO(tsv),
            sep="\t",
            index_col=0,
        )
        ids_in = [str(x) for x in df.index]
        if ids_in != list(df.columns):
            raise RuntimeError("pairsnp output rows/cols mismatched")
        D = df.to_numpy(dtype=float)
        D = (D + D.T) / 2.0
        return SNPMatrix(isolate_ids=ids_in, distances=D)
