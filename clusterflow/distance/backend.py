"""Phase 2.1 — abstract backend interface + auto-selection."""

from __future__ import annotations

import importlib
import logging
import shutil
from abc import ABC, abstractmethod
from pathlib import Path

from clusterflow.models import SNPMatrix

log = logging.getLogger(__name__)


class DistanceBackend(ABC):
    """All backends produce a :class:`SNPMatrix` from FASTA inputs."""

    name: str = "abstract"

    @abstractmethod
    def compute(
        self, fasta_paths: list[Path], n_jobs: int = -1
    ) -> SNPMatrix: ...


def _has_module(name: str) -> bool:
    try:
        importlib.import_module(name)
        return True
    except Exception:
        return False


def _has_cuda() -> bool:
    if not _has_module("cupy"):
        return False
    try:
        cupy = importlib.import_module("cupy")
        return cupy.cuda.is_available()
    except Exception:
        return False


def select_backend(preferred: str = "auto") -> DistanceBackend:
    """Return the backend instance matching the user's preference."""
    from clusterflow.distance.cpu import CPUBackend

    pref = preferred.lower()

    if pref == "cpu":
        return CPUBackend()
    if pref == "gpu":
        from clusterflow.distance.gpu import GPUBackend

        return GPUBackend()
    if pref == "pairsnp":
        from clusterflow.distance.pairsnp import PairSNPBackend

        return PairSNPBackend()

    # auto
    if _has_cuda():
        from clusterflow.distance.gpu import GPUBackend

        log.info("auto-selected GPU backend (cupy + CUDA available)")
        return GPUBackend()
    if shutil.which("pairsnp"):
        from clusterflow.distance.pairsnp import PairSNPBackend

        log.info("auto-selected pairsnp backend")
        return PairSNPBackend()
    log.info("auto-selected CPU backend (joblib parallel Hamming)")
    return CPUBackend()
