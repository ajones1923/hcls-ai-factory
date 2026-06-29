"""Report template seed data ingest pipeline for Imaging Intelligence Agent.

Loads structured radiology report template reference data from a JSON seed
file, validates each record against the ReportTemplate Pydantic model,
and stores embeddings in the imaging_report_templates Milvus collection.

Author: Adam Jones
Date: February 2026
"""

import json
from pathlib import Path
from typing import Any, List

from loguru import logger

from src.models import ReportTemplate

from .base import BaseIngestPipeline


class ReportTemplateIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for report template seed data.

    Reads structured radiology report template definitions from a reference
    JSON file, converts them into ReportTemplate models, and stores
    embeddings in the imaging_report_templates collection.

    Usage:
        pipeline = ReportTemplateIngestPipeline(collection_manager, embedder)
        count = pipeline.run()
    """

    COLLECTION_NAME = "imaging_report_templates"

    def fetch(self, **kwargs) -> List[Any]:
        """Load report template seed data from JSON file.

        Args:
            **kwargs: Unused; reserved for interface compatibility.

        Returns:
            List of report template data dicts from the seed file.
        """
        data_path = (
            Path(__file__).resolve().parent.parent.parent
            / "data"
            / "reference"
            / "report_template_seed_data.json"
        )
        logger.info(f"Loading report template seed data from {data_path}")
        with open(data_path) as f:
            return json.load(f)

    def parse(self, raw_data: List[Any]) -> List[ReportTemplate]:
        """Parse report template data dicts into ReportTemplate models.

        Args:
            raw_data: List of dicts from fetch(), each containing fields
                matching the ReportTemplate model.

        Returns:
            List of validated ReportTemplate model instances.
        """
        records = []
        for item in raw_data:
            try:
                records.append(ReportTemplate(**item))
            except Exception as e:
                logger.warning(f"Skipping invalid report template record: {e}")
        return records
