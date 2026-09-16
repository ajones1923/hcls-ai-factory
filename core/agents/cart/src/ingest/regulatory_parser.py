"""Regulatory milestone data ingest pipeline for CAR-T Intelligence Agent.

Author: Adam Jones
Date: February 2026
"""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from loguru import logger

from src.vector_collections import CARTCollectionManager
from src.models import RegulatoryRecord, RegulatoryEvent

from .base import BaseIngestPipeline



# The curated corpus records regulatory events as DISPLAY LABELS ("Breakthrough Therapy
# Designation", "Type II Variation — 2L DLBCL"), while RegulatoryEvent holds slugs. Passing the
# label straight to the enum raised ValueError for all 40 records, which the parser caught and
# logged -- so the regulatory corpus silently never loaded. Map label -> enum, keeping the
# narrower FDA/EMA distinctions, and fall back to OTHER rather than dropping a real record.
_EVENT_ALIASES = {
    "breakthrough therapy designation": RegulatoryEvent.BREAKTHROUGH_THERAPY,
    "rmat designation": RegulatoryEvent.RMAT,
    "prime designation": RegulatoryEvent.PRIME,
    "conditional approval": RegulatoryEvent.CONDITIONAL_APPROVAL,
    "manufacturing & marketing approval": RegulatoryEvent.FULL_APPROVAL,
    "new drug submission approval": RegulatoryEvent.FULL_APPROVAL,
    "odac advisory committee meeting": RegulatoryEvent.ADVISORY_COMMITTEE,
    "final guidance document": RegulatoryEvent.GUIDANCE,
    "era harmonization framework": RegulatoryEvent.GUIDANCE,
    "rwe pilot program": RegulatoryEvent.PILOT_PROGRAM,
    "rare pediatric disease prv award": RegulatoryEvent.PRV_AWARD,
    "who eml committee review": RegulatoryEvent.WHO_EML,
}


def _normalise_event(raw: str) -> "RegulatoryEvent":
    """Resolve a display label or slug to a RegulatoryEvent."""
    txt = (raw or "").strip()
    try:
        return RegulatoryEvent(txt)                      # already a slug
    except ValueError:
        pass
    key = txt.lower()
    if key in _EVENT_ALIASES:
        return _EVENT_ALIASES[key]
    # Labels that carry a scope suffix, e.g. "Supplemental BLA Approval — 2L+ Multiple Myeloma"
    head = key.split("—")[0].split(" - ")[0].strip()
    if head in _EVENT_ALIASES:
        return _EVENT_ALIASES[head]
    if head.startswith("supplemental bla"):
        return RegulatoryEvent.SUPPLEMENTAL_BLA
    if head.startswith("type i") or head.startswith("type ii"):
        return RegulatoryEvent.VARIATION
    if "boxed warning" in head:
        return RegulatoryEvent.BOXED_WARNING
    if "label update" in head:
        return RegulatoryEvent.LABEL_UPDATE
    return RegulatoryEvent.OTHER


class RegulatoryIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for CAR-T FDA regulatory milestone data."""

    COLLECTION_NAME = "cart_regulatory"

    def __init__(self, collection_manager: CARTCollectionManager, embedder: Any,
                 data_dir: Optional[Path] = None):
        super().__init__(collection_manager, embedder)
        self.data_dir = data_dir or Path(__file__).resolve().parents[2] / "data"

    def fetch(self, data_file: Optional[str] = None) -> List[Dict[str, Any]]:
        file_path = Path(data_file) if data_file else self.data_dir / "reference" / "regulatory_seed_data.json"
        with open(file_path, "r") as f:
            records = json.load(f)
        logger.info(f"Loaded {len(records)} regulatory records from {file_path}")
        return records

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[RegulatoryRecord]:
        records = []
        for data in raw_data:
            try:
                if "regulatory_event" in data and isinstance(data["regulatory_event"], str):
                    data["regulatory_event"] = _normalise_event(data["regulatory_event"])
                records.append(RegulatoryRecord(**data))
            except Exception as e:
                logger.warning(f"Failed to parse regulatory record {data.get('id', '?')}: {e}")
        return records

    def run(self, collection_name: Optional[str] = None, batch_size: int = 32,
            **fetch_kwargs) -> int:
        target = collection_name or self.COLLECTION_NAME
        raw = self.fetch(**fetch_kwargs)
        records = self.parse(raw)
        return self.embed_and_store(records, target, batch_size)
