"""
Signaling pathway data ingest pipeline.

Loads curated cancer signaling pathway seed data from a local JSON file
and normalises it into the ``onco_pathways`` collection schema for the
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
    os.path.dirname(__file__), "..", "..", "data", "reference", "pathway_seed_data.json"
)


class PathwayIngestPipeline(BaseIngestPipeline):
    """
    Ingest pipeline for cancer signaling pathway seed data.

    Populates the ``onco_pathways`` Milvus collection with curated
    pathway descriptions, gene membership, druggable nodes, and
    cross-talk relationships.

    Parameters
    ----------
    collection_manager : CollectionManager
        Milvus collection manager instance.
    embedder : object
        Embedding model with ``encode`` method.
    seed_path : str, optional
        Path to the pathway seed data JSON file.
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
            collection_name="onco_pathways",
        )
        self.seed_path = seed_path or os.path.normpath(DEFAULT_SEED_PATH)

    def fetch(
        self,
        query: Optional[str] = None,
        max_results: Optional[int] = None,
    ) -> List[Dict]:
        """Load pathway records from curated seed data JSON."""
        logger.info("Loading pathway seed data from %s", self.seed_path)

        if not os.path.exists(self.seed_path):
            logger.warning("Pathway seed data not found: %s", self.seed_path)
            return []

        try:
            with open(self.seed_path, "r", encoding="utf-8") as f:
                data = json.load(f)
        except (json.JSONDecodeError, OSError) as exc:
            logger.error("Failed to load pathway seed data: %s", exc)
            return []

        records = data if isinstance(data, list) else data.get("pathways", [])
        if max_results:
            records = records[:max_results]

        logger.info("Loaded %d pathway records", len(records))
        return records

    def parse(self, raw_data: List[Dict]) -> List[Dict]:
        """
        Parse pathway seed records into ``onco_pathways`` schema.

        Schema fields:
            - id: Pathway identifier
            - name: Pathway name (e.g. "RAS-MAPK", "PI3K-AKT-mTOR")
            - text: Full pathway description (embedding source)
            - genes: Member genes (comma-separated)
            - druggable_nodes: Targetable nodes in the pathway
            - cancer_types: Cancer types where pathway is frequently altered
            - cross_talk: Related / interacting pathways
            - source_type: Always "pathway"
        """
        records: List[Dict] = []

        for item in raw_data:
            pathway_id = item.get("id", "")
            name = item.get("name", item.get("pathway", ""))
            description = item.get("description", item.get("summary", item.get("text", "")))

            if not description and not name:
                continue

            text = f"{name} signaling pathway. {description}" if description else name

            records.append({
                "id": str(pathway_id) if pathway_id else f"pathway_{len(records)}",
                "name": name,
                "text": text,
                "genes": item.get("genes", ""),
                "druggable_nodes": item.get("druggable_nodes", item.get("targets", "")),
                "cancer_types": item.get("cancer_types", ""),
                "cross_talk": item.get("cross_talk", item.get("interactions", "")),
                "source_type": "pathway",
            })

        logger.info("Parsed %d pathway records", len(records))
        return records
