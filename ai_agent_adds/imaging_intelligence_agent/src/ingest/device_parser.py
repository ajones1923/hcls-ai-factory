"""Imaging device seed data ingest pipeline for Imaging Intelligence Agent.

Loads FDA-cleared AI/ML medical device reference data from a JSON seed
file, validates each record against the ImagingDevice Pydantic model,
and stores embeddings in the imaging_devices Milvus collection.

Author: Adam Jones
Date: February 2026
"""

import json
from pathlib import Path
from typing import Any, List

from loguru import logger

from src.models import ImagingDevice

from .base import BaseIngestPipeline


class DeviceIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for imaging device seed data.

    Reads FDA-cleared AI/ML medical device definitions from a reference
    JSON file, converts them into ImagingDevice models, and stores
    embeddings in the imaging_devices collection.

    Usage:
        pipeline = DeviceIngestPipeline(collection_manager, embedder)
        count = pipeline.run()
    """

    COLLECTION_NAME = "imaging_devices"

    def fetch(self, **kwargs) -> List[Any]:
        """Load imaging device seed data from JSON file.

        Args:
            **kwargs: Unused; reserved for interface compatibility.

        Returns:
            List of device data dicts from the seed file.
        """
        data_path = (
            Path(__file__).resolve().parent.parent.parent
            / "data"
            / "reference"
            / "device_seed_data.json"
        )
        logger.info(f"Loading device seed data from {data_path}")
        with open(data_path) as f:
            return json.load(f)

    def parse(self, raw_data: List[Any]) -> List[ImagingDevice]:
        """Parse device data dicts into ImagingDevice models.

        Args:
            raw_data: List of dicts from fetch(), each containing fields
                matching the ImagingDevice model.

        Returns:
            List of validated ImagingDevice model instances.
        """
        records = []
        for item in raw_data:
            try:
                records.append(ImagingDevice(**item))
            except Exception as e:
                logger.warning(f"Skipping invalid device record: {e}")
        return records
