"""Phase 7 — threshold-based TRANSMISSION_ALERT events."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime


@dataclass
class TransmissionAlert:
    isolate_id: str
    cluster_id: int
    risk_score: float
    timestamp: datetime
    reason: str
