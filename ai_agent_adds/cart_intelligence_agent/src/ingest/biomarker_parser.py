"""Biomarker data ingest pipeline for CAR-T Intelligence Agent.

Author: Adam Jones
Date: February 2026
"""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from loguru import logger

from src.collections import CARTCollectionManager
from src.models import BiomarkerRecord, BiomarkerType, EvidenceLevel

from .base import BaseIngestPipeline


class BiomarkerIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for CAR-T predictive and pharmacodynamic biomarker data."""

    COLLECTION_NAME = "cart_biomarkers"

    def __init__(self, collection_manager: CARTCollectionManager, embedder: Any,
                 data_dir: Optional[Path] = None):
        super().__init__(collection_manager, embedder)
        self.data_dir = data_dir or Path(__file__).resolve().parents[2] / "data"

    def fetch(self, data_file: Optional[str] = None) -> List[Dict[str, Any]]:
        file_path = Path(data_file) if data_file else self.data_dir / "reference" / "biomarker_seed_data.json"
        with open(file_path, "r") as f:
            records = json.load(f)
        logger.info(f"Loaded {len(records)} biomarker records from {file_path}")
        return records

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[BiomarkerRecord]:
        records = []
        for data in raw_data:
            try:
                if "biomarker_type" in data and isinstance(data["biomarker_type"], str):
                    data["biomarker_type"] = BiomarkerType(data["biomarker_type"])
                if "evidence_level" in data and isinstance(data["evidence_level"], str):
                    data["evidence_level"] = EvidenceLevel(data["evidence_level"])
                records.append(BiomarkerRecord(**data))
            except Exception as e:
                logger.warning(f"Failed to parse biomarker record {data.get('id', '?')}: {e}")
        return records

    def run(self, collection_name: Optional[str] = None, batch_size: int = 32,
            **fetch_kwargs) -> int:
        target = collection_name or self.COLLECTION_NAME
        raw = self.fetch(**fetch_kwargs)
        records = self.parse(raw)
        return self.embed_and_store(records, target, batch_size)
