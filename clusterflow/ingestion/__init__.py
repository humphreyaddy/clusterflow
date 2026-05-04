"""Phase 1.4 — ingestion layer (SNP matrix, MLST, epi metadata)."""

from __future__ import annotations

import logging
from dataclasses import dataclass

from clusterflow.config import ClusterFlowConfig
from clusterflow.ingestion.epi import EpiMetadataReader, merge_mlst
from clusterflow.ingestion.mlst import MLSTReader
from clusterflow.ingestion.snp_matrix import IngestionError, SNPMatrixReader
from clusterflow.models import Isolate, SNPMatrix

log = logging.getLogger(__name__)


@dataclass
class IngestionResult:
    snp: SNPMatrix
    isolates: dict[str, Isolate]


class IngestionPipeline:
    """Loads SNP matrix + MLST + epi metadata, cross-validates IDs."""

    def load(self, config: ClusterFlowConfig) -> IngestionResult:
        snp = SNPMatrixReader().read(config.inputs.snp_matrix)
        epi = EpiMetadataReader().read(config.inputs.epi_metadata)

        if config.inputs.mlst_profiles is not None:
            mlst = MLSTReader().read(config.inputs.mlst_profiles)
            epi = merge_mlst(epi, mlst)
            covered = sum(
                1 for iso in epi.values() if iso.sequence_type is not None
            )
            log.info(
                "MLST coverage: %d / %d isolates (%.1f%%)",
                covered,
                len(epi),
                100 * covered / max(len(epi), 1),
            )

        # Cross-validate: every SNP-matrix ID must be in epi
        missing_in_epi = [i for i in snp.isolate_ids if i not in epi]
        if missing_in_epi:
            raise IngestionError(
                f"{len(missing_in_epi)} isolate(s) in SNP matrix missing from epi "
                f"metadata; first few: {missing_in_epi[:5]}"
            )
        # Drop epi rows that are not in the matrix (warn)
        extra = [i for i in epi if i not in snp.isolate_ids]
        if extra:
            log.warning(
                "%d epi rows not in SNP matrix — dropped: first few %s",
                len(extra),
                extra[:5],
            )
            epi = {k: v for k, v in epi.items() if k in snp.isolate_ids}

        log.info(
            "ingested %d isolates; %d wards; %d STs",
            len(epi),
            len({iso.ward for iso in epi.values()}),
            len({iso.sequence_type for iso in epi.values() if iso.sequence_type}),
        )
        return IngestionResult(snp=snp, isolates=epi)


__all__ = [
    "EpiMetadataReader",
    "IngestionError",
    "IngestionPipeline",
    "IngestionResult",
    "MLSTReader",
    "SNPMatrixReader",
]
