"""
Outcome capture ingest pipeline.

Loads curated treatment outcome seed data from a local JSON file and
normalises it into the ``onco_outcomes`` collection schema for the
Oncology Intelligence Agent's RAG knowledge base.

Outcome data enables the agent to learn from historical treatment
responses, progression patterns, and survival metrics to improve
therapy ranking and open-question generation.

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
    os.path.dirname(__file__), "..", "..", "data", "reference", "outcome_seed_data.json"
)


class OutcomeIngestPipeline(BaseIngestPipeline):
    """
    Ingest pipeline for treatment outcome seed data.

    Populates the ``onco_outcomes`` Milvus collection with curated
    outcome records capturing response rates, progression-free survival,
    overall survival, and adverse event profiles for specific
    biomarker-therapy combinations.

    Parameters
    ----------
    collection_manager : CollectionManager
        Milvus collection manager instance.
    embedder : object
        Embedding model with ``encode`` method.
    seed_path : str, optional
        Path to the outcome seed data JSON file.
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
            collection_name="onco_outcomes",
        )
        self.seed_path = seed_path or os.path.normpath(DEFAULT_SEED_PATH)

    def fetch(
        self,
        query: Optional[str] = None,
        max_results: Optional[int] = None,
    ) -> List[Dict]:
        """Load outcome records from curated seed data JSON."""
        logger.info("Loading outcome seed data from %s", self.seed_path)

        if not os.path.exists(self.seed_path):
            logger.warning("Outcome seed data not found: %s", self.seed_path)
            return []

        try:
            with open(self.seed_path, "r", encoding="utf-8") as f:
                data = json.load(f)
        except (json.JSONDecodeError, OSError) as exc:
            logger.error("Failed to load outcome seed data: %s", exc)
            return []

        records = data if isinstance(data, list) else data.get("outcomes", [])
        if max_results:
            records = records[:max_results]

        logger.info("Loaded %d outcome records", len(records))
        return records

    def parse(self, raw_data: List[Dict]) -> List[Dict]:
        """
        Parse outcome seed records into ``onco_outcomes`` schema.

        Schema fields:
            - id: Outcome record identifier
            - text: Full outcome summary (embedding source)
            - cancer_type: Cancer type
            - gene: Biomarker gene
            - variant: Specific variant (if applicable)
            - drug: Therapy / drug name
            - response_rate: Objective response rate (e.g. "42%")
            - pfs_months: Median progression-free survival (months)
            - os_months: Median overall survival (months)
            - line_of_therapy: Line of therapy (1L, 2L+, etc.)
            - trial_id: Associated clinical trial (NCT ID)
            - source: Data source reference
            - source_type: Always "outcome"
        """
        records: List[Dict] = []

        for item in raw_data:
            outcome_id = item.get("id", "")
            cancer_type = item.get("cancer_type", "")
            gene = item.get("gene", item.get("biomarker", ""))
            variant = item.get("variant", "")
            drug = item.get("drug", item.get("therapy", ""))
            description = item.get("description", item.get("summary", item.get("text", "")))

            # Build embedding text
            parts = []
            if cancer_type:
                parts.append(f"Cancer: {cancer_type}.")
            if gene:
                parts.append(f"Biomarker: {gene} {variant}." if variant else f"Biomarker: {gene}.")
            if drug:
                parts.append(f"Therapy: {drug}.")

            response_rate = item.get("response_rate", item.get("orr", ""))
            pfs = item.get("pfs_months", item.get("pfs", ""))
            os_val = item.get("os_months", item.get("os", ""))

            if response_rate:
                parts.append(f"Response rate: {response_rate}.")
            if pfs:
                parts.append(f"Median PFS: {pfs} months.")
            if os_val:
                parts.append(f"Median OS: {os_val} months.")
            if description:
                parts.append(description)

            text = " ".join(parts) if parts else str(item)

            records.append({
                "id": str(outcome_id) if outcome_id else f"outcome_{len(records)}",
                "text": text,
                "cancer_type": cancer_type,
                "gene": gene,
                "variant": variant,
                "drug": drug,
                "response_rate": str(response_rate),
                "pfs_months": str(pfs),
                "os_months": str(os_val),
                "line_of_therapy": item.get("line_of_therapy", item.get("line", "")),
                "trial_id": item.get("trial_id", item.get("nct_id", "")),
                "source": item.get("source", item.get("reference", "")),
                "source_type": "outcome",
            })

        logger.info("Parsed %d outcome records", len(records))
        return records
