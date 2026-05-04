"""MLST profile reader."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from clusterflow.ingestion.snp_matrix import IngestionError, _sniff_delim


class MLSTReader:
    """Read MLST profiles from a CSV/TSV.

    Required columns: ``isolate_id``, ``ST``. Allele columns are accepted but
    optional. The ST field is normalised: a leading ``"ST"`` prefix is
    stripped, and ``"novel"`` / ``"-"`` / ``"–"`` / blanks are mapped to
    ``"novel"``.
    """

    def read(self, path: str | Path) -> dict[str, str]:
        p = Path(path)
        if not p.exists():
            raise IngestionError(f"MLST file not found: {p}")
        delim = _sniff_delim(p)
        try:
            df = pd.read_csv(p, sep=delim, dtype=str, keep_default_na=False)
        except Exception as e:
            raise IngestionError(f"could not parse MLST {p}: {e}") from e

        cols = {c.lower(): c for c in df.columns}
        if "isolate_id" not in cols or "st" not in cols:
            raise IngestionError(
                f"MLST {p} must contain 'isolate_id' and 'ST' columns; got {list(df.columns)}"
            )
        id_col = cols["isolate_id"]
        st_col = cols["st"]

        out: dict[str, str] = {}
        for i, row in df.iterrows():
            iso = str(row[id_col]).strip()
            raw = str(row[st_col]).strip()
            out[iso] = self._normalise_st(raw)
        return out

    @staticmethod
    def _normalise_st(value: str) -> str:
        if value == "" or value.lower() in {"novel", "-", "–", "na", "nan"}:
            return "novel"
        if value.upper().startswith("ST"):
            return value[2:].strip() or "novel"
        return value
