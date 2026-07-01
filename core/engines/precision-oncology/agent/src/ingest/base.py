"""
Abstract base class for Oncology Intelligence Agent ingest pipelines.

Provides the standard fetch -> parse -> embed_and_store orchestration
pattern used by all data-source adapters (CIViC, PubMed, ClinicalTrials.gov,
curated seed files, etc.).

Follows the same structural pattern as the CAR-T Agent ingest pipelines
for consistency across the HCLS AI Factory project.

Author: Adam Jones
Date: February 2026
"""

import logging
from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


class BaseIngestPipeline(ABC):
    """
    Abstract base class for oncology data ingest pipelines.

    Subclasses implement ``fetch`` and ``parse`` to handle source-specific
    retrieval and normalisation.  The base class provides ``embed_and_store``
    for vectorising text fields and inserting records into Milvus, as well
    as the ``run`` orchestrator that ties the three steps together.

    Parameters
    ----------
    collection_manager : CollectionManager
        Milvus collection manager instance (wraps pymilvus operations).
    embedder : object
        Embedding model with an ``encode(texts: List[str]) -> List[List[float]]``
        method (e.g. SentenceTransformer or a thin wrapper around one).
    collection_name : str
        Name of the target Milvus collection.
    batch_size : int, optional
        Number of records to embed and insert per batch (default 50).
    """

    def __init__(
        self,
        collection_manager: Any,
        embedder: Any,
        collection_name: str,
        batch_size: int = 50,
    ) -> None:
        self.collection_manager = collection_manager
        self.embedder = embedder
        self.collection_name = collection_name
        self.batch_size = batch_size

    # ------------------------------------------------------------------
    # Public orchestrator
    # ------------------------------------------------------------------

    def run(
        self,
        query: Optional[str] = None,
        max_results: Optional[int] = None,
    ) -> int:
        """
        Orchestrate the full ingest pipeline: fetch -> parse -> embed_and_store.

        Parameters
        ----------
        query : str, optional
            Search query forwarded to ``fetch``.
        max_results : int, optional
            Maximum number of raw records to retrieve.

        Returns
        -------
        int
            Total number of records inserted into the collection.
        """
        logger.info(
            "[%s] Starting ingest pipeline (collection=%s, query=%r, max_results=%s)",
            self.__class__.__name__,
            self.collection_name,
            query,
            max_results,
        )

        # Step 1 — Fetch raw data from the source
        kwargs: Dict[str, Any] = {}
        if query is not None:
            kwargs["query"] = query
        if max_results is not None:
            kwargs["max_results"] = max_results

        raw_data = self.fetch(**kwargs)
        logger.info(
            "[%s] Fetched %d raw records",
            self.__class__.__name__,
            len(raw_data),
        )

        if not raw_data:
            logger.warning("[%s] No data fetched — nothing to ingest.", self.__class__.__name__)
            return 0

        # Step 2 — Parse into collection-ready records
        records = self.parse(raw_data)
        logger.info(
            "[%s] Parsed %d records from %d raw entries",
            self.__class__.__name__,
            len(records),
            len(raw_data),
        )

        if not records:
            logger.warning("[%s] No records after parsing — nothing to store.", self.__class__.__name__)
            return 0

        # Step 3 — Embed and store
        total_inserted = self.embed_and_store(records)
        logger.info(
            "[%s] Ingest complete — %d records inserted into '%s'",
            self.__class__.__name__,
            total_inserted,
            self.collection_name,
        )
        return total_inserted

    # ------------------------------------------------------------------
    # Abstract hooks for subclasses
    # ------------------------------------------------------------------

    @abstractmethod
    def fetch(self, query: Optional[str] = None, max_results: Optional[int] = None) -> List[Dict]:
        """
        Fetch raw data from the upstream source.

        Parameters
        ----------
        query : str, optional
            Search / filter query.
        max_results : int, optional
            Cap on number of records returned.

        Returns
        -------
        list of dict
            Raw records as returned by the source API or file.
        """
        ...

    @abstractmethod
    def parse(self, raw_data: List[Dict]) -> List[Dict]:
        """
        Parse raw data into collection-ready records.

        Each returned dict must contain at minimum a ``text`` or
        ``text_chunk`` field that will be embedded, plus any metadata
        fields required by the target collection schema.

        Parameters
        ----------
        raw_data : list of dict
            Raw records from ``fetch``.

        Returns
        -------
        list of dict
            Normalised records ready for embedding and insertion.
        """
        ...

    # ------------------------------------------------------------------
    # Embedding + storage
    # ------------------------------------------------------------------

    def embed_and_store(self, records: List[Dict]) -> int:
        """
        Embed text fields and insert records into the Milvus collection
        in batches.

        The method looks for a ``text`` or ``text_chunk`` or
        ``text_summary`` field in each record to use as the embedding
        source.  The resulting vector is stored under the ``embedding``
        key.

        Parameters
        ----------
        records : list of dict
            Parsed records with text fields to embed.

        Returns
        -------
        int
            Number of records successfully inserted.
        """
        total_inserted = 0

        for batch_start in range(0, len(records), self.batch_size):
            batch = records[batch_start: batch_start + self.batch_size]

            # Resolve text field
            texts: List[str] = []
            for rec in batch:
                text = (
                    rec.get("text")
                    or rec.get("text_chunk")
                    or rec.get("text_summary")
                    or rec.get("summary")
                    or ""
                )
                texts.append(str(text))

            # Embed
            try:
                embeddings = self.embedder.encode(texts)
            except Exception:
                logger.exception(
                    "[%s] Embedding failed for batch starting at index %d",
                    self.__class__.__name__,
                    batch_start,
                )
                continue

            # Attach embeddings to records
            for rec, emb in zip(batch, embeddings):
                rec["embedding"] = emb

            # Insert into Milvus
            try:
                self.collection_manager.insert(
                    collection_name=self.collection_name,
                    records=batch,
                )
                total_inserted += len(batch)
                logger.debug(
                    "[%s] Inserted batch of %d (total so far: %d)",
                    self.__class__.__name__,
                    len(batch),
                    total_inserted,
                )
            except Exception:
                logger.exception(
                    "[%s] Insert failed for batch starting at index %d",
                    self.__class__.__name__,
                    batch_start,
                )

        return total_inserted
