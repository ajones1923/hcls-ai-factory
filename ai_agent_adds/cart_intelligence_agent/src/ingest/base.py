"""Base class for all CAR-T data ingest pipelines.

Each data source (PubMed, ClinicalTrials.gov, construct references, assay
data, manufacturing records) has its own ingest pipeline that inherits from
BaseIngestPipeline.  The base class enforces a standard fetch -> parse ->
embed_and_store workflow.

Author: Adam Jones
Date: February 2026
"""

from abc import ABC, abstractmethod
from enum import Enum
from typing import Any, List, Optional

from loguru import logger
from pydantic import BaseModel

from src.collections import CARTCollectionManager


class BaseIngestPipeline(ABC):
    """Abstract base class for CAR-T data ingest pipelines.

    Subclasses must implement:
      - fetch(**kwargs)   — retrieve raw data from the source
      - parse(raw_data)   — convert raw data into validated Pydantic models

    The base class provides:
      - embed_and_store() — embed text and insert into Milvus
      - run()             — orchestrate the full fetch -> parse -> store pipeline

    Usage:
        class MyPipeline(BaseIngestPipeline):
            def fetch(self, **kwargs): ...
            def parse(self, raw_data): ...

        pipeline = MyPipeline(collection_manager, embedder)
        pipeline.run(max_results=100)
    """

    def __init__(
        self,
        collection_manager: CARTCollectionManager,
        embedder: Any,
    ):
        """Initialize the ingest pipeline.

        Args:
            collection_manager: CARTCollectionManager instance for Milvus
                operations (insert, search, etc.).
            embedder: Embedding model or client that provides an `encode()`
                method returning List[List[float]].  Expected to be a
                SentenceTransformer or compatible wrapper using
                BGE-small-en-v1.5 (384-dim).
        """
        self.collection_manager = collection_manager
        self.embedder = embedder

    @abstractmethod
    def fetch(self, **kwargs) -> Any:
        """Fetch raw data from the upstream source.

        This method handles API calls, file reads, or any other I/O
        needed to obtain the raw data for ingestion.

        Args:
            **kwargs: Source-specific parameters (e.g. max_results, query).

        Returns:
            Raw data in the source's native format (XML, JSON, CSV rows, etc.).
        """
        ...

    @abstractmethod
    def parse(self, raw_data: Any) -> List[BaseModel]:
        """Parse raw data into validated Pydantic model instances.

        Args:
            raw_data: Output from fetch() in the source's native format.

        Returns:
            List of Pydantic model instances (CARTLiterature, ClinicalTrial,
            CARConstruct, AssayResult, or ManufacturingRecord).
        """
        ...

    def embed_and_store(
        self,
        records: List[BaseModel],
        collection_name: str,
        batch_size: int = 32,
    ) -> int:
        """Embed record text and insert into the target Milvus collection.

        Calls each record's `to_embedding_text()` method to produce the
        string that gets embedded, then inserts records in batches.

        Args:
            records: List of Pydantic model instances.  Each must have a
                `to_embedding_text() -> str` method.
            collection_name: Target Milvus collection name.
            batch_size: Number of records to embed and insert at a time.

        Returns:
            Total number of records inserted.

        Raises:
            RuntimeError: If embedding or insertion fails.
        """
        total_inserted = 0

        for i in range(0, len(records), batch_size):
            batch = records[i : i + batch_size]

            try:
                # Build embedding texts from each record
                texts = [record.to_embedding_text() for record in batch]

                # Encode texts into embedding vectors (384-dim each)
                embeddings = self.embedder.encode(texts)

                # Build dicts for insertion with embedding field added
                batch_records = []
                for record, embedding in zip(batch, embeddings):
                    record_dict = record.model_dump()
                    record_dict["embedding"] = embedding

                    # Convert any Enum values to their string .value
                    # and truncate strings to safe UTF-8 byte lengths
                    for key, value in record_dict.items():
                        if isinstance(value, Enum):
                            record_dict[key] = value.value
                        elif isinstance(value, str):
                            # Safety: truncate to Milvus VARCHAR byte limit
                            encoded = value.encode("utf-8")
                            if len(encoded) > 2990 and key in ("text_chunk", "text_summary"):
                                record_dict[key] = encoded[:2990].decode("utf-8", errors="ignore")
                            elif len(encoded) > 490 and key in ("title", "name", "known_toxicities"):
                                record_dict[key] = encoded[:490].decode("utf-8", errors="ignore")

                    batch_records.append(record_dict)

                self.collection_manager.insert_batch(collection_name, batch_records)
                total_inserted += len(batch_records)

            except Exception as exc:
                logger.error(
                    f"Failed batch {i // batch_size + 1} "
                    f"({i}-{i + len(batch)}) into '{collection_name}': {exc}"
                )
                continue

            logger.info(
                f"Inserted batch {i // batch_size + 1} "
                f"({total_inserted}/{len(records)} records) "
                f"into '{collection_name}'"
            )

        return total_inserted

    def run(
        self,
        collection_name: Optional[str] = None,
        batch_size: int = 32,
        **fetch_kwargs,
    ) -> int:
        """Execute the full ingest pipeline: fetch -> parse -> embed_and_store.

        Args:
            collection_name: Target Milvus collection.  Subclasses typically
                set a default via class attribute.
            batch_size: Batch size for embedding and insertion.
            **fetch_kwargs: Passed through to self.fetch().

        Returns:
            Total number of records ingested.

        """
        raw_data = self.fetch(**fetch_kwargs)
        records = self.parse(raw_data)
        logger.info(f"Parsed {len(records)} records for ingestion")
        total = self.embed_and_store(records, collection_name, batch_size)
        return total
