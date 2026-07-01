"""Manufacturing data ingest pipeline for CAR-T Intelligence Agent.

Parses manufacturing / CMC process records from CSV/JSON files into
ManufacturingRecord models and stores embeddings in the cart_manufacturing
Milvus collection.

Supports process steps: transduction, expansion, harvest, formulation,
release testing, and cryopreservation.

Author: Adam Jones
Date: February 2026
"""

import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from loguru import logger

from src.collections import CARTCollectionManager
from src.models import ManufacturingRecord, ProcessStep

from .base import BaseIngestPipeline


class ManufacturingIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for CAR-T manufacturing / CMC data.

    Reads manufacturing records from CSV or JSON files, converts them into
    ManufacturingRecord models, and stores embeddings in the cart_manufacturing
    Milvus collection.

    Expected CSV/JSON fields:
        id, text_summary, process_step, vector_type, parameter,
        parameter_value, target_spec, met_spec, batch_id, notes

    Usage:
        pipeline = ManufacturingIngestPipeline(collection_manager, embedder)
        count = pipeline.run(data_file="/path/to/manufacturing_data.json")
    """

    COLLECTION_NAME = "cart_manufacturing"

    def __init__(
        self,
        collection_manager: CARTCollectionManager,
        embedder: Any,
        data_dir: Optional[Path] = None,
    ):
        super().__init__(collection_manager, embedder)
        self.data_dir = data_dir or Path(__file__).resolve().parents[2] / "data"

    def fetch(
        self,
        data_file: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        """Fetch manufacturing data from a CSV or JSON file."""
        if data_file:
            file_path = Path(data_file)
        else:
            csv_path = self.data_dir / "manufacturing_data.csv"
            json_path = self.data_dir / "manufacturing_data.json"
            if csv_path.exists():
                file_path = csv_path
            elif json_path.exists():
                file_path = json_path
            else:
                raise FileNotFoundError(
                    f"No manufacturing data file found in {self.data_dir}. "
                    "Provide a data_file path or place manufacturing_data.csv/.json in the data directory."
                )

        suffix = file_path.suffix.lower()
        if suffix == ".csv":
            with open(file_path, "r", newline="") as f:
                reader = csv.DictReader(f)
                records = list(reader)
        elif suffix == ".json":
            with open(file_path, "r") as f:
                records = json.load(f)
        else:
            raise ValueError(f"Unsupported file format: {suffix}. Use .csv or .json.")

        logger.info(f"Loaded {len(records)} manufacturing records from {file_path}")
        return records

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[ManufacturingRecord]:
        """Parse manufacturing data dicts into ManufacturingRecord models."""
        records = []
        for data in raw_data:
            try:
                if "process_step" in data and isinstance(data["process_step"], str):
                    data["process_step"] = ProcessStep(data["process_step"])
                record = ManufacturingRecord(**data)
                records.append(record)
            except Exception as e:
                logger.warning(f"Failed to parse manufacturing record {data.get('id', '?')}: {e}")
                continue
        return records

    def run(
        self,
        collection_name: Optional[str] = None,
        batch_size: int = 32,
        **fetch_kwargs,
    ) -> int:
        target = collection_name or self.COLLECTION_NAME
        raw = self.fetch(**fetch_kwargs)
        records = self.parse(raw)
        return self.embed_and_store(records, target, batch_size)
