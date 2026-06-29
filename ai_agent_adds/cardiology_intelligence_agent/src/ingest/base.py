"""Base class for all Cardiology Intelligence Agent ingest parsers.

Provides the abstract interface that every domain-specific ingest parser
must implement (fetch -> parse), plus shared helpers for record validation,
batch processing, and the standard run() orchestration method.

Each parser targets one of the 12 Milvus collections defined in
``src.collections`` and produces ``IngestRecord`` dataclass instances that
downstream embedding and insertion stages consume.

Author: Adam Jones
Date: March 2026
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, List, Optional

from loguru import logger


# ═══════════════════════════════════════════════════════════════════════
# INGEST RECORD DATACLASS
# ═══════════════════════════════════════════════════════════════════════


@dataclass
class IngestRecord:
    """A single parsed record ready for embedding and Milvus insertion.

    Attributes:
        text: The primary text content to embed (abstract, description,
            guideline text, protocol narrative, etc.).  Must be non-empty.
        metadata: Arbitrary key-value metadata that maps to the target
            collection's scalar fields (e.g. ``pmid``, ``modality``,
            ``device_name``).
        collection: Target Milvus collection name (e.g.
            ``cardio_literature``).
        source: Human-readable data source identifier (e.g.
            ``"PubMed"``, ``"ClinicalTrials.gov"``, ``"ACC/AHA 2022"``).
        source_id: Optional unique identifier within the source
            (e.g. PMID, NCT ID, DOI).
    """

    text: str
    metadata: Dict[str, Any]
    collection: str
    source: str
    source_id: Optional[str] = None


# ═══════════════════════════════════════════════════════════════════════
# BASE INGEST PARSER
# ═══════════════════════════════════════════════════════════════════════


class BaseIngestParser(ABC):
    """Abstract base class for all cardiology ingest parsers.

    Subclasses must implement:
      - ``fetch(**kwargs)``  -- retrieve raw data from the upstream source
      - ``parse(raw_data)``  -- convert raw data into ``IngestRecord`` list

    The base class provides:
      - ``run(**kwargs)``          -- orchestrate fetch -> parse pipeline
      - ``validate_record(rec)``   -- basic sanity check on a record
      - ``filter_valid(records)``  -- batch filter invalid records

    Usage::

        class MyParser(BaseIngestParser):
            def fetch(self, **kwargs): ...
            def parse(self, raw_data): ...

        parser = MyParser("cardio_literature")
        records = parser.run(max_results=500)
    """

    def __init__(self, collection_name: str):
        """Initialize the parser with its target collection.

        Args:
            collection_name: Milvus collection this parser writes to
                (e.g. ``cardio_literature``, ``cardio_trials``).
        """
        self.collection = collection_name
        self.logger = logger.bind(parser=self.__class__.__name__)

    # ── Abstract interface ────────────────────────────────────────────

    @abstractmethod
    def fetch(self, **kwargs) -> List[dict]:
        """Fetch raw data from the upstream source.

        Implementations handle API calls, file reads, web scraping, or
        any other I/O needed to obtain raw records.

        Args:
            **kwargs: Source-specific parameters (e.g. ``max_results``,
                ``years_back``, ``query``).

        Returns:
            List of raw dictionaries in the source's native schema.
        """
        ...

    @abstractmethod
    def parse(self, raw_data: List[dict]) -> List[IngestRecord]:
        """Parse raw source data into validated IngestRecord instances.

        Args:
            raw_data: Output from :meth:`fetch`.

        Returns:
            List of :class:`IngestRecord` instances ready for embedding.
        """
        ...

    # ── Orchestration ─────────────────────────────────────────────────

    def run(self, **kwargs) -> List[IngestRecord]:
        """Execute the full fetch -> parse pipeline.

        Convenience method that chains :meth:`fetch` and :meth:`parse`,
        filters out invalid records, and logs a summary.

        Args:
            **kwargs: Passed through to :meth:`fetch`.

        Returns:
            List of valid :class:`IngestRecord` instances.
        """
        self.logger.info(
            f"Starting ingest for collection '{self.collection}'"
        )
        raw = self.fetch(**kwargs)
        self.logger.info(
            f"Fetched {len(raw)} raw records from upstream source"
        )

        records = self.parse(raw)
        valid = self.filter_valid(records)

        self.logger.info(
            f"Parsed {len(valid)} valid records "
            f"({len(records) - len(valid)} rejected) "
            f"for '{self.collection}'"
        )
        return valid

    # ── Validation helpers ────────────────────────────────────────────

    def validate_record(self, record: IngestRecord) -> bool:
        """Check that a record meets minimum quality requirements.

        A record is valid when:
          - ``text`` is non-empty after stripping whitespace
          - ``collection`` is set and matches this parser's target

        Args:
            record: The :class:`IngestRecord` to validate.

        Returns:
            ``True`` if the record passes validation.
        """
        if not record.text or not record.text.strip():
            return False
        if not record.collection:
            return False
        return True

    def filter_valid(self, records: List[IngestRecord]) -> List[IngestRecord]:
        """Return only records that pass :meth:`validate_record`.

        Args:
            records: List of :class:`IngestRecord` to filter.

        Returns:
            Filtered list containing only valid records.
        """
        valid = [r for r in records if self.validate_record(r)]
        rejected = len(records) - len(valid)
        if rejected > 0:
            self.logger.warning(
                f"Rejected {rejected}/{len(records)} invalid records "
                f"for '{self.collection}'"
            )
        return valid

    # ── Utility helpers ───────────────────────────────────────────────

    @staticmethod
    def truncate(text: str, max_length: int) -> str:
        """Truncate text to *max_length* characters with ellipsis.

        Args:
            text: Input string.
            max_length: Maximum allowed character length.

        Returns:
            Truncated string (with trailing ``...``) if over limit,
            otherwise the original string.
        """
        if len(text) <= max_length:
            return text
        return text[: max_length - 3] + "..."

    @staticmethod
    def safe_str(value: Any, default: str = "") -> str:
        """Safely convert a value to string, returning *default* for None.

        Args:
            value: Any value to convert.
            default: Fallback when *value* is ``None``.

        Returns:
            String representation of *value*.
        """
        if value is None:
            return default
        if isinstance(value, Enum):
            return str(value.value)
        return str(value)
