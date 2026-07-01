"""Radiology report ingest pipeline for Imaging Intelligence Agent.

Parses free-text radiology reports using RadiologyReportParser,
embeds at multiple granularities (whole-report, findings section,
impression section), and stores in the imaging_reports Milvus
collection.

Follows the same BaseIngestPipeline pattern as:
  - src/ingest/literature_parser.py (PubMedImagingIngestPipeline)
  - src/ingest/finding_parser.py (FindingIngestPipeline)

Author: Adam Jones
Date: April 2026
"""

import re
from typing import Any, Dict, List, Optional

from loguru import logger

from src.models import ReportGranularity
from src.report_parser import RadiologyReportParser

from .base import BaseIngestPipeline


def _truncate_utf8(text: str, max_bytes: int) -> str:
    """Truncate a string to fit within max_bytes when UTF-8 encoded.

    Milvus VARCHAR max_length counts bytes, not characters, so we need
    byte-aware truncation for text containing multi-byte characters.
    """
    encoded = text.encode("utf-8")
    if len(encoded) <= max_bytes:
        return text
    return encoded[:max_bytes].decode("utf-8", errors="ignore")


class RadiologyReportIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for radiology reports.

    Parses free-text reports into structured data, then embeds at
    multiple granularities for optimal retrieval:
      - whole_report: full report embedding for broad queries
      - section_findings: findings section for specific finding queries
      - section_impression: impression section for clinical summary queries

    Each granularity is stored as a separate row in the imaging_reports
    Milvus collection with a unique ID:
        {patient_id}_{study_date}_{granularity}

    Usage:
        pipeline = RadiologyReportIngestPipeline(collection_manager, embedder)
        report_id = pipeline.ingest_report(
            report_text="EXAMINATION: CT Chest...",
            patient_id="PT-001",
            study_date="2025-11-15",
            modality="ct",
        )
        count = pipeline.ingest_batch(reports)
    """

    COLLECTION_NAME = "imaging_reports"

    def __init__(self, collection_manager, embedder):
        """Initialize the radiology report ingest pipeline.

        Args:
            collection_manager: ImagingCollectionManager for Milvus operations.
            embedder: Embedding model with encode() method (BGE-small-en-v1.5).
        """
        super().__init__(collection_manager, embedder)
        self.parser = RadiologyReportParser()

    # ── Public API ──────────────────────────────────────────────────

    def ingest_report(
        self,
        report_text: str,
        patient_id: str,
        study_date: str,
        modality: str,
        report_source: str = "manual",
    ) -> str:
        """Ingest a single radiology report into Milvus.

        Parses the report, generates embeddings at multiple granularities,
        and inserts all rows into the imaging_reports collection.

        Args:
            report_text: Free-text radiology report.
            patient_id: Patient identifier.
            study_date: Study date (YYYY-MM-DD).
            modality: Imaging modality (ct, mri, xray, etc.).
            report_source: Source of the report (manual, orthanc, hl7).

        Returns:
            The base report ID (patient_id + study_date).
        """
        if not report_text or not report_text.strip():
            logger.warning(f"Empty report text for patient {patient_id}")
            return ""

        # Parse the report
        parsed = self.parser.parse(report_text)
        logger.info(
            f"Parsed report for {patient_id}/{study_date}: "
            f"{len(parsed.sections)} sections, {len(parsed.entities)} entities, "
            f"{len(parsed.measurements)} measurements"
        )

        # Use detected modality/region if not explicitly provided
        detected_modality = parsed.modality or modality
        detected_region = parsed.body_region or ""

        # Build rows for each granularity
        rows = self._build_rows(
            parsed=parsed,
            patient_id=patient_id,
            study_date=study_date,
            modality=detected_modality,
            body_region=detected_region,
            report_source=report_source,
        )

        if not rows:
            logger.warning(f"No rows generated for {patient_id}/{study_date}")
            return ""

        # Embed and insert
        try:
            texts = [r.pop("_embedding_text") for r in rows]
            embeddings = self.embedder.encode(texts, normalize_embeddings=True).tolist()

            for row, emb in zip(rows, embeddings):
                row["embedding"] = emb

            count = self.collection_manager.insert_batch(
                self.COLLECTION_NAME, rows, batch_size=len(rows)
            )
            logger.info(
                f"Ingested {count} rows for report {patient_id}/{study_date} "
                f"into {self.COLLECTION_NAME}"
            )
        except Exception as exc:
            logger.error(
                f"Failed to ingest report {patient_id}/{study_date}: {exc}"
            )
            return ""

        return f"{patient_id}_{study_date}"

    def ingest_batch(self, reports: List[Dict[str, Any]]) -> int:
        """Ingest a batch of radiology reports.

        Each report dict must contain:
          - report_text: str
          - patient_id: str
          - study_date: str
          - modality: str
        Optional:
          - report_source: str (defaults to "manual")

        Args:
            reports: List of report dicts.

        Returns:
            Total number of reports successfully ingested.
        """
        if not reports:
            logger.warning("Empty report batch provided")
            return 0

        success_count = 0
        for i, report in enumerate(reports):
            try:
                report_text = report.get("report_text", "")
                patient_id = report.get("patient_id", "")
                study_date = report.get("study_date", "")
                modality = report.get("modality", "")
                report_source = report.get("report_source", "manual")

                if not report_text or not patient_id:
                    logger.warning(f"Skipping report {i}: missing text or patient_id")
                    continue

                result = self.ingest_report(
                    report_text=report_text,
                    patient_id=patient_id,
                    study_date=study_date,
                    modality=modality,
                    report_source=report_source,
                )

                if result:
                    success_count += 1

            except Exception as exc:
                logger.error(f"Failed to ingest report {i}: {exc}")
                continue

        logger.info(
            f"Batch ingest complete: {success_count}/{len(reports)} reports ingested"
        )
        return success_count

    # ── BaseIngestPipeline interface (for compatibility) ────────────

    def fetch(self, **kwargs) -> Any:
        """Fetch is not used for report ingestion (reports are pushed).

        This method exists only to satisfy the BaseIngestPipeline
        abstract interface.  Use ingest_report() or ingest_batch()
        instead.
        """
        logger.warning(
            "fetch() is not applicable for radiology report ingestion. "
            "Use ingest_report() or ingest_batch() instead."
        )
        return []

    def parse(self, raw_data: Any) -> list:
        """Parse is not used for report ingestion.

        This method exists only to satisfy the BaseIngestPipeline
        abstract interface.  Parsing is handled internally by
        ingest_report() using RadiologyReportParser.
        """
        return []

    # ── Internal Methods ────────────────────────────────────────────

    def _build_rows(
        self,
        parsed,
        patient_id: str,
        study_date: str,
        modality: str,
        body_region: str,
        report_source: str,
    ) -> List[Dict[str, Any]]:
        """Build Milvus insert rows at multiple granularities.

        Creates up to three rows per report:
          1. whole_report — embedding of the full report summary
          2. section_findings — embedding of findings section only
          3. section_impression — embedding of impression section only

        Each row includes a special _embedding_text key that is popped
        before insertion and used for embedding generation.

        Args:
            parsed: ParsedReport from RadiologyReportParser.
            patient_id: Patient identifier.
            study_date: Study date string.
            modality: Imaging modality.
            body_region: Anatomical body region.
            report_source: Source of the report.

        Returns:
            List of row dicts ready for embedding and insertion.
        """
        rows: List[Dict[str, Any]] = []

        # Common fields for all granularities
        indication = _truncate_utf8(
            parsed.sections.get("clinical_history", ""), 490
        )
        findings_text = _truncate_utf8(
            parsed.sections.get("findings", ""), 4990
        )
        impression_text = _truncate_utf8(
            parsed.sections.get("impression", ""), 1990
        )
        comparison_date = parsed.comparison_date or ""
        findings_count = len([
            e for e in parsed.entities if e.entity_type == "finding"
        ])

        base = {
            "patient_id": _truncate_utf8(patient_id, 45),
            "study_date": _truncate_utf8(study_date, 15),
            "modality": _truncate_utf8(modality, 8),
            "body_region": _truncate_utf8(body_region, 45),
            "indication": indication,
            "findings_text": findings_text,
            "impression_text": impression_text,
            "findings_count": findings_count,
            "critical_finding": parsed.critical_finding,
            "comparison_date": _truncate_utf8(comparison_date, 15),
            "report_source": _truncate_utf8(report_source, 95),
        }

        # 1. Whole report embedding
        whole_text = self.parser.generate_embedding_text(parsed)
        if whole_text.strip():
            row = dict(base)
            row["id"] = _truncate_utf8(
                f"{patient_id}_{study_date}_{ReportGranularity.WHOLE_REPORT.value}",
                195,
            )
            row["granularity"] = ReportGranularity.WHOLE_REPORT.value
            row["_embedding_text"] = whole_text
            rows.append(row)

        # 2. Findings section embedding
        if findings_text.strip():
            row = dict(base)
            row["id"] = _truncate_utf8(
                f"{patient_id}_{study_date}_{ReportGranularity.SECTION_FINDINGS.value}",
                195,
            )
            row["granularity"] = ReportGranularity.SECTION_FINDINGS.value
            row["_embedding_text"] = f"Findings: {findings_text}"
            rows.append(row)

        # 3. Impression section embedding
        if impression_text.strip():
            row = dict(base)
            row["id"] = _truncate_utf8(
                f"{patient_id}_{study_date}_{ReportGranularity.SECTION_IMPRESSION.value}",
                195,
            )
            row["granularity"] = ReportGranularity.SECTION_IMPRESSION.value
            row["_embedding_text"] = f"Impression: {impression_text}"
            rows.append(row)

        return rows
