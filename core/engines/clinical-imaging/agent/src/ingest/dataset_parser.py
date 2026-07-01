"""Public imaging dataset seed data ingest pipeline.

Loads public imaging dataset metadata from a JSON seed file,
validates each record against the DatasetRecord Pydantic model,
and stores embeddings in the imaging_datasets Milvus collection.
"""

import json
from pathlib import Path
from typing import Any, List

from loguru import logger

from src.models import DatasetRecord

from .base import BaseIngestPipeline


class DatasetIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for public imaging dataset metadata.

    Reads dataset definitions from a reference JSON file, converts
    them into DatasetRecord models, and stores embeddings in the
    imaging_datasets collection.

    Usage:
        pipeline = DatasetIngestPipeline(collection_manager, embedder)
        count = pipeline.run()
    """

    COLLECTION_NAME = "imaging_datasets"

    def fetch(self, **kwargs) -> List[Any]:
        """Load dataset seed data from JSON file."""
        data_path = (
            Path(__file__).resolve().parent.parent.parent
            / "data"
            / "reference"
            / "dataset_seed_data.json"
        )
        logger.info(f"Loading dataset seed data from {data_path}")
        with open(data_path) as f:
            return json.load(f)

    def parse(self, raw_data: List[Any]) -> List[DatasetRecord]:
        """Parse dataset data dicts into DatasetRecord models."""
        records = []
        for item in raw_data:
            try:
                records.append(DatasetRecord(**item))
            except Exception as e:
                logger.warning(f"Skipping invalid dataset record: {e}")
        return records
