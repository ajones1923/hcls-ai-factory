"""Base ingest pipeline for Imaging Intelligence Agent.

Abstract base class that provides the fetch->parse->embed_and_store pattern.
All ingest pipelines subclass this and implement fetch() and parse().

Follows the same pattern as:
  - cart_intelligence_agent/src/ingest/base.py (CAR-T Agent)

Author: Adam Jones
Date: February 2026
"""

from abc import ABC, abstractmethod
from typing import Any, List, Optional

from loguru import logger
from pydantic import BaseModel


class BaseIngestPipeline(ABC):
    """Abstract ingest pipeline: fetch -> parse -> embed_and_store.

    Subclasses must implement:
      - fetch(**kwargs)   -- retrieve raw data from the source
      - parse(raw_data)   -- convert raw data into validated Pydantic models

    The base class provides:
      - embed_and_store() -- embed text and insert into Milvus
      - run()             -- orchestrate the full fetch -> parse -> store pipeline

    Usage:
        class MyPipeline(BaseIngestPipeline):
            COLLECTION_NAME = "imaging_xxx"
            def fetch(self, **kwargs): ...
            def parse(self, raw_data): ...

        pipeline = MyPipeline(collection_manager, embedder)
        pipeline.run(max_results=100)
    """

    COLLECTION_NAME: str = ""

    def __init__(self, collection_manager, embedder):
        """Initialize the ingest pipeline.

        Args:
            collection_manager: ImagingCollectionManager instance for Milvus
                operations (insert, search, etc.).
            embedder: Embedding model or client that provides an encode()
                method returning List[List[float]].  Expected to be a
                SentenceTransformer or compatible wrapper using
                BGE-small-en-v1.5 (384-dim).
        """
        self.collection_manager = collection_manager
        self.embedder = embedder

    @abstractmethod
    def fetch(self, **kwargs) -> Any:
        """Retrieve raw data from external source.

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
        """Parse raw data into Pydantic model instances.

        Args:
            raw_data: Output from fetch() in the source's native format.

        Returns:
            List of validated Pydantic model instances, each with a
            to_embedding_text() method.
        """
        ...

    def embed_and_store(
        self,
        records: List[BaseModel],
        collection_name: Optional[str] = None,
        batch_size: int = 32,
    ) -> int:
        """Embed record text and insert into Milvus.

        Calls each record's to_embedding_text() method to produce the
        string that gets embedded, then inserts records in batches.

        Args:
            records: List of Pydantic model instances.  Each must have a
                to_embedding_text() -> str method.
            collection_name: Target Milvus collection name.  Falls back
                to self.COLLECTION_NAME if not provided.
            batch_size: Number of records to embed and insert at a time.

        Returns:
            Total number of records inserted.
        """
        coll_name = collection_name or self.COLLECTION_NAME
        if not records:
            logger.warning(f"No records to store in {coll_name}")
            return 0

        total = 0
        for i in range(0, len(records), batch_size):
            batch = records[i : i + batch_size]

            try:
                texts = [r.to_embedding_text() for r in batch]
                embeddings = self.embedder.encode(texts, normalize_embeddings=True).tolist()

                insert_data = []
                for rec, emb in zip(batch, embeddings):
                    row = rec.model_dump()
                    # Convert enum values to strings for Milvus VARCHAR fields
                    for k, v in row.items():
                        if hasattr(v, "value"):
                            row[k] = v.value
                    row["embedding"] = emb
                    insert_data.append(row)

                count = self.collection_manager.insert_batch(coll_name, insert_data, batch_size=batch_size)
                total += count
                logger.info(f"Inserted {count} records into {coll_name} (batch {i // batch_size + 1})")

            except Exception as exc:
                logger.error(
                    f"Failed batch {i // batch_size + 1} "
                    f"({i}-{i + len(batch)}) into '{coll_name}': {exc}"
                )
                continue

        logger.info(f"Total inserted into {coll_name}: {total}")
        return total

    def run(self, collection_name: Optional[str] = None, batch_size: int = 32, **fetch_kwargs) -> int:
        """Orchestrate fetch -> parse -> embed_and_store.

        Args:
            collection_name: Target Milvus collection.  Subclasses typically
                set a default via COLLECTION_NAME class attribute.
            batch_size: Batch size for embedding and insertion.
            **fetch_kwargs: Passed through to self.fetch().

        Returns:
            Total number of records ingested.
        """
        logger.info(f"Starting ingest pipeline for {self.COLLECTION_NAME}")
        raw = self.fetch(**fetch_kwargs)
        records = self.parse(raw)
        logger.info(f"Parsed {len(records)} records")
        return self.embed_and_store(records, collection_name, batch_size)
