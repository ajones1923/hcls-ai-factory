"""Imaging finding seed data ingest pipeline for Imaging Intelligence Agent.

Loads imaging finding reference data from a JSON seed file, validates
each record against the ImagingFinding Pydantic model, and stores
embeddings in the imaging_findings Milvus collection.

Author: Adam Jones
Date: February 2026
"""

import json
from pathlib import Path
from typing import Any, List

from loguru import logger

from src.models import ImagingFinding

from .base import BaseIngestPipeline


class FindingIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for imaging finding seed data.

    Reads finding definitions from a reference JSON file, converts
    them into ImagingFinding models, and stores embeddings in the
    imaging_findings collection.

    Usage:
        pipeline = FindingIngestPipeline(collection_manager, embedder)
        count = pipeline.run()
    """

    COLLECTION_NAME = "imaging_findings"

    def fetch(self, **kwargs) -> List[Any]:
        """Load imaging finding seed data from JSON file.

        Args:
            **kwargs: Unused; reserved for interface compatibility.

        Returns:
            List of finding data dicts from the seed file.
        """
        data_path = (
            Path(__file__).resolve().parent.parent.parent
            / "data"
            / "reference"
            / "finding_seed_data.json"
        )
        logger.info(f"Loading finding seed data from {data_path}")
        with open(data_path) as f:
            return json.load(f)

    def parse(self, raw_data: List[Any]) -> List[ImagingFinding]:
        """Parse finding data dicts into ImagingFinding models.

        Args:
            raw_data: List of dicts from fetch(), each containing fields
                matching the ImagingFinding model.

        Returns:
            List of validated ImagingFinding model instances.
        """
        records = []
        for item in raw_data:
            try:
                records.append(ImagingFinding(**item))
            except Exception as e:
                logger.warning(f"Skipping invalid finding record: {e}")
        return records
