"""Molecular / structural data ingest pipeline for CAR-T Intelligence Agent.

Author: Adam Jones
Date: February 2026
"""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from loguru import logger

from src.collections import CARTCollectionManager
from src.models import SequenceRecord

from .base import BaseIngestPipeline


class SequenceIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for CAR-T molecular and structural data (scFv, binding affinity)."""

    COLLECTION_NAME = "cart_sequences"

    def __init__(self, collection_manager: CARTCollectionManager, embedder: Any,
                 data_dir: Optional[Path] = None):
        super().__init__(collection_manager, embedder)
        self.data_dir = data_dir or Path(__file__).resolve().parents[2] / "data"

    def fetch(self, data_file: Optional[str] = None) -> List[Dict[str, Any]]:
        file_path = Path(data_file) if data_file else self.data_dir / "reference" / "sequence_seed_data.json"
        with open(file_path, "r") as f:
            records = json.load(f)
        logger.info(f"Loaded {len(records)} sequence records from {file_path}")
        return records

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[SequenceRecord]:
        records = []
        for data in raw_data:
            try:
                records.append(SequenceRecord(**data))
            except Exception as e:
                logger.warning(f"Failed to parse sequence record {data.get('id', '?')}: {e}")
        return records

    def run(self, collection_name: Optional[str] = None, batch_size: int = 32,
            **fetch_kwargs) -> int:
        target = collection_name or self.COLLECTION_NAME
        raw = self.fetch(**fetch_kwargs)
        records = self.parse(raw)
        return self.embed_and_store(records, target, batch_size)
