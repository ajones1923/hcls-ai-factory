"""Anatomy seed data ingest pipeline for Imaging Intelligence Agent.

Loads anatomical structure reference data from a JSON seed file, validates
each record against the AnatomyRecord Pydantic model, and stores
embeddings in the imaging_anatomy Milvus collection.

Author: Adam Jones
Date: February 2026
"""

import json
from pathlib import Path
from typing import Any, List

from loguru import logger

from src.models import AnatomyRecord

from .base import BaseIngestPipeline


class AnatomyIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for anatomy seed data.

    Reads anatomical structure definitions from a reference JSON file,
    converts them into AnatomyRecord models, and stores embeddings in
    the imaging_anatomy collection.

    Usage:
        pipeline = AnatomyIngestPipeline(collection_manager, embedder)
        count = pipeline.run()
    """

    COLLECTION_NAME = "imaging_anatomy"

    def fetch(self, **kwargs) -> List[Any]:
        """Load anatomy seed data from JSON file.

        Args:
            **kwargs: Unused; reserved for interface compatibility.

        Returns:
            List of anatomy data dicts from the seed file.
        """
        data_path = (
            Path(__file__).resolve().parent.parent.parent
            / "data"
            / "reference"
            / "anatomy_seed_data.json"
        )
        logger.info(f"Loading anatomy seed data from {data_path}")
        with open(data_path) as f:
            return json.load(f)

    def parse(self, raw_data: List[Any]) -> List[AnatomyRecord]:
        """Parse anatomy data dicts into AnatomyRecord models.

        Args:
            raw_data: List of dicts from fetch(), each containing fields
                matching the AnatomyRecord model.

        Returns:
            List of validated AnatomyRecord model instances.
        """
        records = []
        for item in raw_data:
            try:
                records.append(AnatomyRecord(**item))
            except Exception as e:
                logger.warning(f"Skipping invalid anatomy record: {e}")
        return records
