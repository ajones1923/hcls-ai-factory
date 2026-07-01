"""
Resistance mechanism ingest pipeline.

Loads curated drug resistance mechanism seed data from a local JSON file
and normalises it into the ``onco_resistance`` collection schema for the
Oncology Intelligence Agent's RAG knowledge base.

Author: Adam Jones
Date: February 2026
"""

import json
import logging
import os
from typing import Any, Dict, List, Optional

from src.ingest.base import BaseIngestPipeline

logger = logging.getLogger(__name__)

DEFAULT_SEED_PATH = os.path.join(
    os.path.dirname(__file__), "..", "..", "data", "reference", "resistance_seed_data.json"
)


class ResistanceIngestPipeline(BaseIngestPipeline):
    """
    Ingest pipeline for drug resistance mechanism seed data.

    Populates the ``onco_resistance`` Milvus collection with curated
    descriptions of known resistance mechanisms, including the drug(s)
    affected, genomic alterations driving resistance, and potential
    strategies to overcome resistance.

    Parameters
    ----------
    collection_manager : CollectionManager
        Milvus collection manager instance.
    embedder : object
        Embedding model with ``encode`` method.
    seed_path : str, optional
        Path to the resistance seed data JSON file.
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
            collection_name="onco_resistance",
        )
        self.seed_path = seed_path or os.path.normpath(DEFAULT_SEED_PATH)

    def fetch(
        self,
        query: Optional[str] = None,
        max_results: Optional[int] = None,
    ) -> List[Dict]:
        """Load resistance mechanism records from curated seed data JSON."""
        logger.info("Loading resistance seed data from %s", self.seed_path)

        if not os.path.exists(self.seed_path):
            logger.warning("Resistance seed data not found: %s", self.seed_path)
            return []

        try:
            with open(self.seed_path, "r", encoding="utf-8") as f:
                data = json.load(f)
        except (json.JSONDecodeError, OSError) as exc:
            logger.error("Failed to load resistance seed data: %s", exc)
            return []

        records = data if isinstance(data, list) else data.get("resistance_mechanisms", [])
        if max_results:
            records = records[:max_results]

        logger.info("Loaded %d resistance records", len(records))
        return records

    def parse(self, raw_data: List[Dict]) -> List[Dict]:
        """
        Parse resistance seed records into ``onco_resistance`` schema.

        Schema fields:
            - id: Resistance mechanism identifier
            - name: Mechanism name
            - text: Full description (embedding source)
            - drug: Drug(s) affected by this resistance mechanism
            - gene: Gene(s) involved in resistance
            - mutation: Specific mutation(s) driving resistance
            - cancer_type: Cancer type context
            - strategy: Strategies to overcome resistance
            - source_type: Always "resistance"
        """
        records: List[Dict] = []

        for item in raw_data:
            res_id = item.get("id", "")
            name = item.get("name", item.get("mechanism", ""))
            description = item.get("description", item.get("summary", item.get("text", "")))
            drug = item.get("drug", item.get("drugs", ""))

            if not description and not name:
                continue

            text = f"Resistance mechanism: {name}. Drug: {drug}. {description}" if description else name

            records.append({
                "id": str(res_id) if res_id else f"resistance_{len(records)}",
                "name": name,
                "text": text,
                "drug": drug if isinstance(drug, str) else ", ".join(drug) if drug else "",
                "gene": item.get("gene", item.get("genes", "")),
                "mutation": item.get("mutation", item.get("mutations", "")),
                "cancer_type": item.get("cancer_type", ""),
                "strategy": item.get("strategy", item.get("overcome", "")),
                "source_type": "resistance",
            })

        logger.info("Parsed %d resistance records", len(records))
        return records
