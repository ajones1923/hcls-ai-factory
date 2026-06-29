"""Benchmark seed data ingest pipeline for Imaging Intelligence Agent.

Loads model benchmark / performance reference data from a JSON seed file,
validates each record against the BenchmarkRecord Pydantic model, and
stores embeddings in the imaging_benchmarks Milvus collection.

Author: Adam Jones
Date: February 2026
"""

import json
from pathlib import Path
from typing import Any, List

from loguru import logger

from src.models import BenchmarkRecord

from .base import BaseIngestPipeline


class BenchmarkIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for benchmark seed data.

    Reads model benchmark definitions from a reference JSON file,
    converts them into BenchmarkRecord models, and stores embeddings
    in the imaging_benchmarks collection.

    Usage:
        pipeline = BenchmarkIngestPipeline(collection_manager, embedder)
        count = pipeline.run()
    """

    COLLECTION_NAME = "imaging_benchmarks"

    def fetch(self, **kwargs) -> List[Any]:
        """Load benchmark seed data from JSON file.

        Args:
            **kwargs: Unused; reserved for interface compatibility.

        Returns:
            List of benchmark data dicts from the seed file.
        """
        data_path = (
            Path(__file__).resolve().parent.parent.parent
            / "data"
            / "reference"
            / "benchmark_seed_data.json"
        )
        logger.info(f"Loading benchmark seed data from {data_path}")
        with open(data_path) as f:
            return json.load(f)

    def parse(self, raw_data: List[Any]) -> List[BenchmarkRecord]:
        """Parse benchmark data dicts into BenchmarkRecord models.

        Args:
            raw_data: List of dicts from fetch(), each containing fields
                matching the BenchmarkRecord model.

        Returns:
            List of validated BenchmarkRecord model instances.
        """
        records = []
        for item in raw_data:
            try:
                records.append(BenchmarkRecord(**item))
            except Exception as e:
                logger.warning(f"Skipping invalid benchmark record: {e}")
        return records
