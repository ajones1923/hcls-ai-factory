"""Real-world evidence data ingest pipeline for CAR-T Intelligence Agent.

Author: Adam Jones
Date: February 2026
"""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from loguru import logger

from src.collections import CARTCollectionManager
from src.models import RealWorldRecord, RWEStudyType

from .base import BaseIngestPipeline


class RealWorldIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for CAR-T real-world evidence and outcomes data."""

    COLLECTION_NAME = "cart_realworld"

    def __init__(self, collection_manager: CARTCollectionManager, embedder: Any,
                 data_dir: Optional[Path] = None):
        super().__init__(collection_manager, embedder)
        self.data_dir = data_dir or Path(__file__).resolve().parents[2] / "data"

    def fetch(self, data_file: Optional[str] = None) -> List[Dict[str, Any]]:
        file_path = Path(data_file) if data_file else self.data_dir / "reference" / "realworld_seed_data.json"
        with open(file_path, "r") as f:
            records = json.load(f)
        logger.info(f"Loaded {len(records)} real-world evidence records from {file_path}")
        return records

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[RealWorldRecord]:
        records = []
        for data in raw_data:
            try:
                if "study_type" in data and isinstance(data["study_type"], str):
                    data["study_type"] = RWEStudyType(data["study_type"])
                if "population_size" in data:
                    data["population_size"] = int(data["population_size"])
                if "median_followup_months" in data:
                    data["median_followup_months"] = float(data["median_followup_months"])
                records.append(RealWorldRecord(**data))
            except Exception as e:
                logger.warning(f"Failed to parse RWE record {data.get('id', '?')}: {e}")
        return records

    def run(self, collection_name: Optional[str] = None, batch_size: int = 32,
            **fetch_kwargs) -> int:
        target = collection_name or self.COLLECTION_NAME
        raw = self.fetch(**fetch_kwargs)
        records = self.parse(raw)
        return self.embed_and_store(records, target, batch_size)
