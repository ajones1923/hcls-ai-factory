"""Assay data ingest pipeline for CAR-T Intelligence Agent.

Parses in vitro and in vivo assay results from CSV/JSON files into
AssayResult models and stores embeddings in the cart_assays Milvus
collection.

Supports assay types: cytotoxicity, cytokine, flow cytometry,
proliferation, in vivo, persistence, and exhaustion.

Author: Adam Jones
Date: February 2026
"""

import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from loguru import logger

from src.collections import CARTCollectionManager
from src.models import AssayResult, AssayType

from .base import BaseIngestPipeline


class AssayIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for CAR-T assay result data.

    Reads assay records from CSV or JSON files, converts them into
    AssayResult models, and stores embeddings in the cart_assays
    Milvus collection.

    Expected CSV/JSON fields:
        id, text_summary, assay_type, construct_id, target_antigen,
        cell_line, effector_ratio, key_metric, metric_value, outcome, notes

    Usage:
        pipeline = AssayIngestPipeline(collection_manager, embedder)
        count = pipeline.run(data_file="/path/to/assay_results.csv")
    """

    COLLECTION_NAME = "cart_assays"

    def __init__(
        self,
        collection_manager: CARTCollectionManager,
        embedder: Any,
        data_dir: Optional[Path] = None,
    ):
        """Initialize the assay ingest pipeline.

        Args:
            collection_manager: CARTCollectionManager for Milvus operations.
            embedder: Embedding model with encode() method.
            data_dir: Directory containing assay data files.
                Defaults to the project data/ directory.
        """
        super().__init__(collection_manager, embedder)
        self.data_dir = data_dir or Path(__file__).resolve().parents[2] / "data"

    def fetch(
        self,
        data_file: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        """Fetch assay data from a CSV or JSON file.

        Supports two file formats:
          - CSV: comma-separated with header row matching AssayResult fields
          - JSON: list of objects with keys matching AssayResult fields

        Args:
            data_file: Path to a CSV or JSON file containing assay records.
                If None, looks for 'assay_results.csv' or 'assay_results.json'
                in self.data_dir.

        Returns:
            List of assay record dicts.

        """
        if data_file:
            file_path = Path(data_file)
        else:
            # Look for default files in self.data_dir
            csv_path = self.data_dir / "assay_results.csv"
            json_path = self.data_dir / "assay_results.json"
            if csv_path.exists():
                file_path = csv_path
            elif json_path.exists():
                file_path = json_path
            else:
                raise FileNotFoundError(
                    f"No assay data file found in {self.data_dir}. "
                    "Provide a data_file path or place assay_results.csv/.json in the data directory."
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

        logger.info(f"Loaded {len(records)} assay records from {file_path}")
        return records

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[AssayResult]:
        """Parse assay data dicts into AssayResult models.

        Maps assay_type strings to the AssayType enum and validates
        numeric metric values.

        Args:
            raw_data: List of dicts from fetch(), each containing fields
                matching the AssayResult model.

        Returns:
            List of validated AssayResult model instances.

        """
        records = []
        for data in raw_data:
            try:
                # Map assay_type string to AssayType enum if needed
                if "assay_type" in data and isinstance(data["assay_type"], str):
                    data["assay_type"] = AssayType(data["assay_type"])

                # Convert metric_value to float
                if "metric_value" in data:
                    data["metric_value"] = float(data["metric_value"])

                record = AssayResult(**data)
                records.append(record)
            except Exception as e:
                logger.warning(f"Failed to parse assay record {data.get('id', '?')}: {e}")
                continue
        return records

    def run(
        self,
        collection_name: Optional[str] = None,
        batch_size: int = 32,
        **fetch_kwargs,
    ) -> int:
        """Execute the full assay data ingest pipeline.

        Args:
            collection_name: Target collection (defaults to 'cart_assays').
            batch_size: Batch size for embedding and insertion.
            **fetch_kwargs: Passed to fetch() (data_file).

        Returns:
            Total number of records ingested.

        """
        target = collection_name or self.COLLECTION_NAME
        raw = self.fetch(**fetch_kwargs)
        records = self.parse(raw)
        return self.embed_and_store(records, target, batch_size)
