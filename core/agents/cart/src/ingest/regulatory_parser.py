"""Regulatory milestone data ingest pipeline for CAR-T Intelligence Agent.

Author: Adam Jones
Date: February 2026
"""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from loguru import logger

from src.collections import CARTCollectionManager
from src.models import RegulatoryRecord, RegulatoryEvent

from .base import BaseIngestPipeline


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
                    data["regulatory_event"] = RegulatoryEvent(data["regulatory_event"])
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
