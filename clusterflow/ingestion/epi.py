"""Epidemiological metadata reader."""

from __future__ import annotations

import logging
from datetime import date, datetime
from pathlib import Path

import pandas as pd

from clusterflow.ingestion.snp_matrix import IngestionError, _sniff_delim
from clusterflow.models import Isolate

log = logging.getLogger(__name__)


_DATE_FORMATS = ("%Y-%m-%d", "%d/%m/%Y", "%d-%m-%Y", "%Y/%m/%d")


def _parse_date(s: str) -> date:
    s = s.strip()
    for fmt in _DATE_FORMATS:
        try:
            return datetime.strptime(s, fmt).date()
        except ValueError:
            continue
    # Last resort: pandas parser
    try:
        return pd.to_datetime(s).date()
    except Exception as e:
        raise IngestionError(f"could not parse date: {s!r}") from e


class EpiMetadataReader:
    """Read epi metadata: isolate_id, collection_date, facility, ward."""

    REQUIRED = ("isolate_id", "collection_date", "facility", "ward")

    def read(self, path: str | Path) -> dict[str, Isolate]:
        p = Path(path)
        if not p.exists():
            raise IngestionError(f"epi metadata not found: {p}")
        delim = _sniff_delim(p)
        try:
            df = pd.read_csv(p, sep=delim, dtype=str, keep_default_na=False)
        except Exception as e:
            raise IngestionError(f"could not parse epi {p}: {e}") from e

        cols = {c.lower(): c for c in df.columns}
        missing = [c for c in self.REQUIRED if c not in cols]
        if missing:
            raise IngestionError(
                f"epi {p} missing required columns: {missing}; got {list(df.columns)}"
            )

        out: dict[str, Isolate] = {}
        today = date.today()
        min_date: date | None = None
        max_date: date | None = None
        for _, row in df.iterrows():
            iso = str(row[cols["isolate_id"]]).strip()
            d = _parse_date(str(row[cols["collection_date"]]))
            if d > today:
                raise IngestionError(f"future collection date for {iso}: {d}")
            min_date = d if min_date is None or d < min_date else min_date
            max_date = d if max_date is None or d > max_date else max_date
            out[iso] = Isolate(
                isolate_id=iso,
                collection_date=d,
                facility=str(row[cols["facility"]]).strip() or "unknown",
                ward=str(row[cols["ward"]]).strip() or "unknown",
            )

        if min_date is not None and max_date is not None:
            span = (max_date - min_date).days
            if span > 365 * 5:
                log.warning(
                    "epi date span %d days (>5 years) — please verify", span
                )
        return out


def merge_mlst(isolates: dict[str, Isolate], mlst: dict[str, str]) -> dict[str, Isolate]:
    """Return a new isolate dict with `sequence_type` populated from MLST."""
    out: dict[str, Isolate] = {}
    for iso_id, iso in isolates.items():
        st = mlst.get(iso_id)
        if st is None:
            out[iso_id] = iso
        else:
            out[iso_id] = iso.model_copy(update={"sequence_type": st})
    return out
