"""Guideline seed data ingest pipeline for Imaging Intelligence Agent.

Loads clinical practice guideline reference data from a JSON seed file,
validates each record against the GuidelineRecord Pydantic model, and
stores embeddings in the imaging_guidelines Milvus collection.

Author: Adam Jones
Date: February 2026
"""

import json
from pathlib import Path
from typing import Any, List

from loguru import logger

from src.models import GuidelineRecord

from .base import BaseIngestPipeline


class GuidelineIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for guideline seed data.

    Reads clinical practice guideline definitions from a reference JSON
    file, converts them into GuidelineRecord models, and stores embeddings
    in the imaging_guidelines collection.

    Usage:
        pipeline = GuidelineIngestPipeline(collection_manager, embedder)
        count = pipeline.run()
    """

    COLLECTION_NAME = "imaging_guidelines"

    def fetch(self, **kwargs) -> List[Any]:
        """Load guideline seed data from JSON file.

        Args:
            **kwargs: Unused; reserved for interface compatibility.

        Returns:
            List of guideline data dicts from the seed file.
        """
        data_path = (
            Path(__file__).resolve().parent.parent.parent
            / "data"
            / "reference"
            / "guideline_seed_data.json"
        )
        logger.info(f"Loading guideline seed data from {data_path}")
        with open(data_path) as f:
            return json.load(f)

    def parse(self, raw_data: List[Any]) -> List[GuidelineRecord]:
        """Parse guideline data dicts into GuidelineRecord models.

        Args:
            raw_data: List of dicts from fetch(), each containing fields
                matching the GuidelineRecord model.

        Returns:
            List of validated GuidelineRecord model instances.
        """
        records = []
        for item in raw_data:
            try:
                records.append(GuidelineRecord(**item))
            except Exception as e:
                logger.warning(f"Skipping invalid guideline record: {e}")
        return records
