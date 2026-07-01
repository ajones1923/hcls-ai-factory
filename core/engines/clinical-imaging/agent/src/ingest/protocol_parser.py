"""Imaging protocol seed data ingest pipeline for Imaging Intelligence Agent.

Loads imaging protocol reference data from a JSON seed file, validates
each record against the ImagingProtocol Pydantic model, and stores
embeddings in the imaging_protocols Milvus collection.

Author: Adam Jones
Date: February 2026
"""

import json
from pathlib import Path
from typing import Any, List

from loguru import logger

from src.models import ImagingProtocol

from .base import BaseIngestPipeline


class ProtocolIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for imaging protocol seed data.

    Reads protocol definitions from a reference JSON file, converts
    them into ImagingProtocol models, and stores embeddings in the
    imaging_protocols collection.

    Usage:
        pipeline = ProtocolIngestPipeline(collection_manager, embedder)
        count = pipeline.run()
    """

    COLLECTION_NAME = "imaging_protocols"

    def fetch(self, **kwargs) -> List[Any]:
        """Load imaging protocol seed data from JSON file.

        Args:
            **kwargs: Unused; reserved for interface compatibility.

        Returns:
            List of protocol data dicts from the seed file.
        """
        data_path = (
            Path(__file__).resolve().parent.parent.parent
            / "data"
            / "reference"
            / "protocol_seed_data.json"
        )
        logger.info(f"Loading protocol seed data from {data_path}")
        with open(data_path) as f:
            return json.load(f)

    def parse(self, raw_data: List[Any]) -> List[ImagingProtocol]:
        """Parse protocol data dicts into ImagingProtocol models.

        Args:
            raw_data: List of dicts from fetch(), each containing fields
                matching the ImagingProtocol model.

        Returns:
            List of validated ImagingProtocol model instances.
        """
        records = []
        for item in raw_data:
            try:
                records.append(ImagingProtocol(**item))
            except Exception as e:
                logger.warning(f"Skipping invalid protocol record: {e}")
        return records
