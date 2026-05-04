"""SNP distance matrix reader (TSV / CSV / snp-dists output).

snp-dists output uses tab separation with the first column header empty.
Generic TSV/CSV is auto-detected by sniffing the delimiter from the file.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pandas as pd

from clusterflow.models import SNPMatrix


class IngestionError(Exception):
    """Raised when an input file cannot be parsed or fails validation."""


def _sniff_delim(path: Path) -> str:
    with path.open() as fh:
        sample = fh.read(8192)
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t,;|").delimiter
    except csv.Error:
        return "\t" if "\t" in sample else ","


class SNPMatrixReader:
    """Read a square symmetric pairwise SNP distance matrix.

    Accepts TSV, CSV, or `snp-dists` output. The first column must contain
    isolate IDs; column headers must contain matching isolate IDs in the same
    order. Raises :class:`IngestionError` with a descriptive message on the
    first violation.
    """

    def read(self, path: str | Path) -> SNPMatrix:
        p = Path(path)
        if not p.exists():
            raise IngestionError(f"SNP matrix file not found: {p}")
        delim = _sniff_delim(p)
        try:
            df = pd.read_csv(p, sep=delim, index_col=0)
        except Exception as e:
            raise IngestionError(f"could not parse {p} (delim={delim!r}): {e}") from e

        if df.shape[0] == 0 or df.shape[1] == 0:
            raise IngestionError(f"SNP matrix is empty: {p}")
        if df.shape[0] != df.shape[1]:
            raise IngestionError(
                f"SNP matrix is not square ({df.shape[0]}x{df.shape[1]}): {p}"
            )

        row_ids = [str(x) for x in df.index]
        col_ids = [str(c) for c in df.columns]
        if row_ids != col_ids:
            mismatches = [
                (i, r, c) for i, (r, c) in enumerate(zip(row_ids, col_ids)) if r != c
            ]
            raise IngestionError(
                f"row/column header mismatch in {p}: first mismatch at index "
                f"{mismatches[0][0]}: row={mismatches[0][1]!r} col={mismatches[0][2]!r}"
            )

        try:
            arr = df.to_numpy(dtype=float)
        except Exception as e:
            raise IngestionError(f"non-numeric values in {p}: {e}") from e

        if np.any(np.isnan(arr)):
            i, j = np.argwhere(np.isnan(arr))[0]
            raise IngestionError(f"NaN at position ({row_ids[i]}, {col_ids[j]}) in {p}")
        if np.any(arr < 0):
            i, j = np.argwhere(arr < 0)[0]
            raise IngestionError(
                f"negative distance at ({row_ids[i]}, {col_ids[j]})={arr[i, j]} in {p}"
            )
        if not np.allclose(np.diag(arr), 0, atol=1e-6):
            i = int(np.argmax(np.abs(np.diag(arr))))
            raise IngestionError(
                f"non-zero diagonal at {row_ids[i]}: {arr[i, i]} (expected 0) in {p}"
            )
        if not np.allclose(arr, arr.T, atol=1e-6):
            d = np.abs(arr - arr.T)
            i, j = np.unravel_index(np.argmax(d), arr.shape)
            raise IngestionError(
                f"asymmetric distance at ({row_ids[i]}, {col_ids[j]}): "
                f"{arr[i, j]} vs {arr[j, i]} (delta={d[i, j]:.6f}) in {p}"
            )

        # Symmetrise to remove float jitter, round to int (SNP counts).
        arr = (arr + arr.T) / 2.0
        return SNPMatrix(isolate_ids=row_ids, distances=arr)
