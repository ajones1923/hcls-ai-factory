"""Base ingest parser for the Clinical Trial Intelligence Agent.

Defines the abstract BaseIngestParser interface and the IngestRecord dataclass
that all concrete parsers (ClinicalTrials.gov, PubMed, Regulatory) implement.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


# ===================================================================
# INGEST RECORD DATACLASS
# ===================================================================


@dataclass
class IngestRecord:
    """A single record ready for embedding and insertion into Milvus.

    Attributes:
        text: The text content to be embedded.
        metadata: Key-value metadata (source, date, identifiers, etc.).
        collection_name: Target Milvus collection name.
        record_id: Optional unique identifier for deduplication.
        source: Data source identifier (e.g., "clinicaltrials", "pubmed").
    """

    text: str
    metadata: Dict[str, Any] = field(default_factory=dict)
    collection_name: str = ""
    record_id: Optional[str] = None
    source: str = ""

    def __post_init__(self) -> None:
        """Validate the record after creation."""
        if not self.text or not self.text.strip():
            raise ValueError("IngestRecord text must not be empty.")

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to a flat dictionary for storage or debugging."""
        return {
            "text": self.text,
            "metadata": self.metadata,
            "collection_name": self.collection_name,
            "record_id": self.record_id,
            "source": self.source,
        }


# ===================================================================
# INGEST STATS DATACLASS
# ===================================================================


@dataclass
class IngestStats:
    """Statistics from a single ingest run.

    Attributes:
        source: Data source identifier.
        total_fetched: Number of raw records fetched from the source.
        total_parsed: Number of records successfully parsed.
        total_validated: Number of records that passed validation.
        total_errors: Number of records that failed validation or parsing.
        duration_seconds: Total wall-clock time for the ingest run.
        error_details: List of error messages for failed records.
    """

    source: str = ""
    total_fetched: int = 0
    total_parsed: int = 0
    total_validated: int = 0
    total_errors: int = 0
    duration_seconds: float = 0.0
    error_details: List[str] = field(default_factory=list)


# ===================================================================
# BASE INGEST PARSER ABC
# ===================================================================


class BaseIngestParser(ABC):
    """Abstract base class for all ingest parsers.

    Concrete parsers must implement:
      - ``fetch()``          -- retrieve raw data from an external source
      - ``parse()``          -- transform raw data into IngestRecord objects
      - ``validate_record()`` -- validate a single IngestRecord

    The ``run()`` method orchestrates the full pipeline:
    fetch -> parse -> validate -> return validated records.

    Usage::

        class MyParser(BaseIngestParser):
            def fetch(self, **kwargs):
                return [{"id": "1", "text": "example"}]

            def parse(self, raw_data):
                return [IngestRecord(text=r["text"], source="my_source") for r in raw_data]

            def validate_record(self, record):
                return len(record.text) > 10

        parser = MyParser(source_name="my_source")
        records, stats = parser.run(max_results=100)
    """

    def __init__(
        self,
        source_name: str = "unknown",
        collection_manager: Any = None,
        embedder: Any = None,
    ) -> None:
        """Initialize the base parser.

        Args:
            source_name: Human-readable name of this data source.
            collection_manager: Optional Milvus collection manager for insertion.
            embedder: Optional embedding model with an ``encode()`` method.
        """
        self.source_name = source_name
        self.collection_manager = collection_manager
        self.embedder = embedder
        self.logger = logging.getLogger(f"{__name__}.{source_name}")

    @abstractmethod
    def fetch(self, **kwargs: Any) -> List[Dict[str, Any]]:
        """Fetch raw data from the external source.

        Args:
            **kwargs: Source-specific parameters (query, max_results, etc.).

        Returns:
            List of raw data dictionaries.
        """
        ...

    @abstractmethod
    def parse(self, raw_data: List[Dict[str, Any]]) -> List[IngestRecord]:
        """Parse raw data into IngestRecord objects.

        Args:
            raw_data: List of raw dictionaries from ``fetch()``.

        Returns:
            List of IngestRecord objects.
        """
        ...

    @abstractmethod
    def validate_record(self, record: IngestRecord) -> bool:
        """Validate a single IngestRecord.

        Args:
            record: The record to validate.

        Returns:
            True if the record is valid and should be included.
        """
        ...

    def run(self, **kwargs: Any) -> tuple:
        """Execute the full ingest pipeline: fetch -> parse -> validate.

        Args:
            **kwargs: Passed through to ``fetch()``.

        Returns:
            Tuple of (validated_records, IngestStats).
        """
        stats = IngestStats(source=self.source_name)
        start = time.time()

        try:
            # Fetch
            self.logger.info("Fetching data from %s ...", self.source_name)
            raw_data = self.fetch(**kwargs)
            stats.total_fetched = len(raw_data)
            self.logger.info("Fetched %d raw records from %s", len(raw_data), self.source_name)

            # Parse
            records = self.parse(raw_data)
            stats.total_parsed = len(records)
            self.logger.info("Parsed %d records from %s", len(records), self.source_name)

            # Validate
            validated: List[IngestRecord] = []
            for record in records:
                try:
                    if self.validate_record(record):
                        validated.append(record)
                    else:
                        stats.total_errors += 1
                        stats.error_details.append(
                            f"Validation failed for record: {record.record_id or 'unknown'}"
                        )
                except Exception as exc:
                    stats.total_errors += 1
                    stats.error_details.append(
                        f"Validation error for {record.record_id or 'unknown'}: {exc}"
                    )

            stats.total_validated = len(validated)
            self.logger.info(
                "Validated %d / %d records from %s (%d errors)",
                len(validated), len(records), self.source_name, stats.total_errors,
            )

        except Exception as exc:
            self.logger.error("Ingest pipeline failed for %s: %s", self.source_name, exc)
            stats.total_errors += 1
            stats.error_details.append(f"Pipeline error: {exc}")
            validated = []

        stats.duration_seconds = time.time() - start
        return validated, stats
