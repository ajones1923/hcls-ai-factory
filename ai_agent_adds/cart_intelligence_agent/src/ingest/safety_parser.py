"""Safety / pharmacovigilance data ingest pipeline for CAR-T Intelligence Agent.

Author: Adam Jones
Date: February 2026
"""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from loguru import logger

from src.collections import CARTCollectionManager
from src.models import SafetyRecord, SafetyEventType

from .base import BaseIngestPipeline


class SafetyIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for CAR-T pharmacovigilance / post-market safety data."""

    COLLECTION_NAME = "cart_safety"

    def __init__(self, collection_manager: CARTCollectionManager, embedder: Any,
                 data_dir: Optional[Path] = None):
        super().__init__(collection_manager, embedder)
        self.data_dir = data_dir or Path(__file__).resolve().parents[2] / "data"

    def fetch(self, data_file: Optional[str] = None) -> List[Dict[str, Any]]:
        file_path = Path(data_file) if data_file else self.data_dir / "reference" / "safety_seed_data.json"
        with open(file_path, "r") as f:
            records = json.load(f)
        logger.info(f"Loaded {len(records)} safety records from {file_path}")
        return records

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[SafetyRecord]:
        records = []
        for data in raw_data:
            try:
                if "event_type" in data and isinstance(data["event_type"], str):
                    data["event_type"] = SafetyEventType(data["event_type"])
                if "year" in data:
                    data["year"] = int(data["year"])
                records.append(SafetyRecord(**data))
            except Exception as e:
                logger.warning(f"Failed to parse safety record {data.get('id', '?')}: {e}")
        return records

    def run(self, collection_name: Optional[str] = None, batch_size: int = 32,
            **fetch_kwargs) -> int:
        target = collection_name or self.COLLECTION_NAME
        raw = self.fetch(**fetch_kwargs)
        records = self.parse(raw)
        return self.embed_and_store(records, target, batch_size)
