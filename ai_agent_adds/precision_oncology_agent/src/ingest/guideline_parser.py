"""
NCCN / ESMO / ASCO guideline ingest pipeline.

Loads curated clinical practice guideline seed data from a local JSON
file and normalises it into the ``onco_guidelines`` collection schema
for the Oncology Intelligence Agent's RAG knowledge base.

Note: Full guideline content is typically copyrighted. This pipeline
ingests curated summaries and references that have been manually
prepared as seed data.

Author: Adam Jones
Date: February 2026
"""

import json
import logging
import os
from typing import Any, Dict, List, Optional

from src.ingest.base import BaseIngestPipeline

logger = logging.getLogger(__name__)

# Default path to the seed data file (relative to project root)
DEFAULT_SEED_PATH = os.path.join(
    os.path.dirname(__file__), "..", "..", "data", "reference", "guideline_seed_data.json"
)


class GuidelineIngestPipeline(BaseIngestPipeline):
    """
    Ingest pipeline for clinical practice guideline seed data.

    Populates the ``onco_guidelines`` Milvus collection with curated
    summaries of NCCN, ESMO, and ASCO treatment guidelines for major
    cancer types and biomarker-driven therapy recommendations.

    Parameters
    ----------
    collection_manager : CollectionManager
        Milvus collection manager instance.
    embedder : object
        Embedding model with ``encode`` method.
    seed_path : str, optional
        Path to the guideline seed data JSON file.
    """

    def __init__(
        self,
        collection_manager: Any,
        embedder: Any,
        seed_path: Optional[str] = None,
    ) -> None:
        super().__init__(
            collection_manager=collection_manager,
            embedder=embedder,
            collection_name="onco_guidelines",
        )
        self.seed_path = seed_path or os.path.normpath(DEFAULT_SEED_PATH)

    # ------------------------------------------------------------------
    # Fetch
    # ------------------------------------------------------------------

    def fetch(
        self,
        query: Optional[str] = None,
        max_results: Optional[int] = None,
    ) -> List[Dict]:
        """
        Load guideline records from the curated seed data JSON file.

        Parameters
        ----------
        query : str, optional
            Not used (all seed data is loaded).
        max_results : int, optional
            Limit on records returned (default: all).

        Returns
        -------
        list of dict
            Raw guideline records from the seed file.
        """
        logger.info("Loading guideline seed data from %s", self.seed_path)

        if not os.path.exists(self.seed_path):
            logger.warning(
                "Guideline seed data file not found: %s — "
                "run seed data generation first.",
                self.seed_path,
            )
            return []

        try:
            with open(self.seed_path, "r", encoding="utf-8") as f:
                data = json.load(f)
        except (json.JSONDecodeError, OSError) as exc:
            logger.error("Failed to load guideline seed data: %s", exc)
            return []

        records = data if isinstance(data, list) else data.get("guidelines", [])
        if max_results:
            records = records[:max_results]

        logger.info("Loaded %d guideline records", len(records))
        return records

    # ------------------------------------------------------------------
    # Parse
    # ------------------------------------------------------------------

    def parse(self, raw_data: List[Dict]) -> List[Dict]:
        """
        Parse guideline seed records into ``onco_guidelines`` schema.

        Schema fields:
            - id: Unique guideline identifier
            - title: Guideline title
            - text: Full guideline summary text (embedding source)
            - source: Issuing body (NCCN, ESMO, ASCO)
            - cancer_type: Cancer type covered
            - version: Guideline version / year
            - category: Category (e.g. "treatment", "screening", "biomarker")
            - recommendation_level: Strength of recommendation
            - genes: Relevant genes
            - drugs: Relevant drugs
            - source_type: Always "guideline"

        Parameters
        ----------
        raw_data : list of dict
            Raw guideline records from ``fetch``.

        Returns
        -------
        list of dict
            Normalised records ready for embedding and insertion.
        """
        records: List[Dict] = []

        for item in raw_data:
            guideline_id = item.get("id", "")
            title = item.get("title", "")
            summary = item.get("summary", item.get("text", item.get("description", "")))

            if not summary and not title:
                continue

            text = f"{title}. {summary}" if summary else title

            records.append({
                "id": str(guideline_id) if guideline_id else f"guide_{len(records)}",
                "title": title,
                "text": text,
                "source": item.get("source", item.get("issuing_body", "")),
                "cancer_type": item.get("cancer_type", ""),
                "version": item.get("version", item.get("year", "")),
                "category": item.get("category", ""),
                "recommendation_level": item.get("recommendation_level", item.get("level", "")),
                "genes": item.get("genes", ""),
                "drugs": item.get("drugs", ""),
                "source_type": "guideline",
            })

        logger.info("Parsed %d guideline records", len(records))
        return records
