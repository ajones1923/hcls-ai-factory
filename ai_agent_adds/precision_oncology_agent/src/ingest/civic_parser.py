"""
CIViC (Clinical Interpretation of Variants in Cancer) ingest pipeline.

Fetches variant and evidence data from the CIViC public API and normalises
it into the ``onco_variants`` collection schema for the Oncology Intelligence
Agent's RAG knowledge base.

CIViC evidence levels:
    A — Validated association (FDA / guideline)
    B — Clinical evidence (clinical trial)
    C — Case study (case reports)
    D — Preclinical (biological rationale)
    E — Inferential association

Reference: https://civicdb.org/

Author: Adam Jones
Date: February 2026
"""

import logging
import time
from typing import Any, Dict, List, Optional

import requests

from src.ingest.base import BaseIngestPipeline

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CIVIC_API_BASE = "https://civicdb.org/api"
CIVIC_VARIANTS_ENDPOINT = f"{CIVIC_API_BASE}/variants"
CIVIC_EVIDENCE_ENDPOINT = f"{CIVIC_API_BASE}/evidence_items"

# Map CIViC evidence levels to our internal EvidenceLevel enum values
CIVIC_LEVEL_MAP = {
    "A": "level_1",   # Validated / FDA-approved
    "B": "level_2",   # Clinical evidence
    "C": "level_3",   # Case study
    "D": "level_4",   # Preclinical
    "E": "level_4",   # Inferential — closest match to preclinical
}

REQUEST_TIMEOUT = 30  # seconds
PAGE_SIZE = 50        # CIViC API page size
RATE_LIMIT_DELAY = 0.25  # seconds between requests


class CIViCIngestPipeline(BaseIngestPipeline):
    """
    Ingest pipeline for CIViC variant and evidence data.

    Populates the ``onco_variants`` Milvus collection with clinically
    interpreted cancer variants, including associated evidence items,
    drug associations, and cancer-type annotations.
    """

    def __init__(self, collection_manager: Any, embedder: Any) -> None:
        super().__init__(
            collection_manager=collection_manager,
            embedder=embedder,
            collection_name="onco_variants",
        )

    # ------------------------------------------------------------------
    # Fetch
    # ------------------------------------------------------------------

    def fetch(
        self,
        query: Optional[str] = None,
        max_results: int = 5000,
    ) -> List[Dict]:
        """
        Fetch variant records (with nested evidence) from the CIViC API.

        Uses paginated GET requests to ``/api/variants`` and enriches
        each variant with its evidence items from ``/api/evidence_items``.

        Parameters
        ----------
        query : str, optional
            Not used for CIViC (fetches all variants). Reserved for
            future gene-specific filtering.
        max_results : int
            Maximum number of variant records to retrieve (default 5000).

        Returns
        -------
        list of dict
            Raw CIViC variant records with an added ``evidence_items`` key.
        """
        variants: List[Dict] = []
        page = 1
        total_fetched = 0

        logger.info("Fetching variants from CIViC API (max_results=%d)", max_results)

        while total_fetched < max_results:
            try:
                params = {"count": PAGE_SIZE, "page": page}
                if query:
                    params["entrez_symbol"] = query

                response = requests.get(
                    CIVIC_VARIANTS_ENDPOINT,
                    params=params,
                    timeout=REQUEST_TIMEOUT,
                )
                response.raise_for_status()
                data = response.json()

            except requests.RequestException as exc:
                logger.error("CIViC API request failed on page %d: %s", page, exc)
                break

            records = data.get("records", [])
            if not records:
                logger.info("No more records on page %d — stopping.", page)
                break

            variants.extend(records)
            total_fetched += len(records)

            # Check if we've hit the last page
            total_available = data.get("_meta", {}).get("total_count", max_results)
            if total_fetched >= total_available:
                break

            page += 1
            time.sleep(RATE_LIMIT_DELAY)

        logger.info("Fetched %d variant records from CIViC", len(variants))

        # Enrich each variant with evidence items
        for variant in variants[:max_results]:
            variant["evidence_items"] = self._fetch_evidence_for_variant(
                variant.get("id")
            )
            time.sleep(RATE_LIMIT_DELAY)

        return variants[:max_results]

    def _fetch_evidence_for_variant(self, variant_id: Optional[int]) -> List[Dict]:
        """Fetch evidence items for a specific CIViC variant."""
        if variant_id is None:
            return []

        try:
            response = requests.get(
                f"{CIVIC_VARIANTS_ENDPOINT}/{variant_id}/evidence_items",
                timeout=REQUEST_TIMEOUT,
            )
            response.raise_for_status()
            return response.json() if isinstance(response.json(), list) else response.json().get("records", [])
        except requests.RequestException as exc:
            logger.warning("Failed to fetch evidence for variant %d: %s", variant_id, exc)
            return []

    # ------------------------------------------------------------------
    # Parse
    # ------------------------------------------------------------------

    def parse(self, raw_data: List[Dict]) -> List[Dict]:
        """
        Parse CIViC variant records into ``onco_variants`` collection schema.

        Schema fields:
            - id: CIViC variant ID (string)
            - gene: Gene symbol (e.g. "BRAF")
            - variant_name: Variant name (e.g. "V600E")
            - variant_type: Type of variant (e.g. "missense_variant")
            - cancer_type: Associated cancer type(s)
            - evidence_level: Mapped to our EvidenceLevel enum
            - drugs: Associated drugs / therapies
            - civic_id: CIViC identifier for cross-reference
            - text_summary: Concatenated human-readable summary for embedding
            - source_type: Always "civic"

        Parameters
        ----------
        raw_data : list of dict
            Raw CIViC variant records from ``fetch``.

        Returns
        -------
        list of dict
            Normalised records ready for embedding and insertion.
        """
        records: List[Dict] = []

        for variant in raw_data:
            civic_id = variant.get("id", "")
            gene = variant.get("entrez_name", "") or variant.get("gene", {}).get("name", "")
            variant_name = variant.get("name", "")
            variant_types = variant.get("variant_types", [])
            variant_type = variant_types[0].get("display_name", "") if variant_types else ""

            evidence_items = variant.get("evidence_items", [])

            # If no evidence items, create a single record from the variant itself
            if not evidence_items:
                summary = self._build_variant_summary(gene, variant_name, variant_type)
                records.append({
                    "id": f"civic_v{civic_id}",
                    "gene": gene,
                    "variant_name": variant_name,
                    "variant_type": variant_type,
                    "cancer_type": "",
                    "evidence_level": "level_4",
                    "drugs": "",
                    "civic_id": str(civic_id),
                    "text_summary": summary,
                    "text": summary,
                    "source_type": "civic",
                })
                continue

            # Create one record per evidence item for richer embeddings
            for ev_idx, evidence in enumerate(evidence_items):
                disease = evidence.get("disease", {})
                cancer_type = disease.get("display_name", "") if disease else ""

                civic_level = evidence.get("evidence_level", "")
                evidence_level = self._map_civic_evidence_level(civic_level)

                # Extract drugs
                drugs_list = evidence.get("drugs", [])
                drugs = ", ".join(
                    d.get("name", "") for d in drugs_list if d.get("name")
                ) if drugs_list else ""

                ev_type = evidence.get("evidence_type", "")
                ev_direction = evidence.get("evidence_direction", "")
                ev_significance = evidence.get("clinical_significance", evidence.get("significance", ""))
                ev_description = evidence.get("description", "")

                summary = self._build_evidence_summary(
                    gene=gene,
                    variant_name=variant_name,
                    variant_type=variant_type,
                    cancer_type=cancer_type,
                    drugs=drugs,
                    evidence_type=ev_type,
                    evidence_direction=ev_direction,
                    clinical_significance=ev_significance,
                    description=ev_description,
                )

                records.append({
                    "id": f"civic_v{civic_id}_e{ev_idx}",
                    "gene": gene,
                    "variant_name": variant_name,
                    "variant_type": variant_type,
                    "cancer_type": cancer_type,
                    "evidence_level": evidence_level,
                    "drugs": drugs,
                    "civic_id": str(civic_id),
                    "text_summary": summary,
                    "text": summary,
                    "source_type": "civic",
                })

        logger.info("Parsed %d records from %d CIViC variants", len(records), len(raw_data))
        return records

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _map_civic_evidence_level(civic_level: str) -> str:
        """
        Map CIViC evidence level (A-E) to our internal EvidenceLevel enum.

        CIViC levels:
            A -> level_1 (Validated / FDA)
            B -> level_2 (Clinical evidence)
            C -> level_3 (Case study)
            D -> level_4 (Preclinical)
            E -> level_4 (Inferential)

        Parameters
        ----------
        civic_level : str
            Single-character CIViC evidence level.

        Returns
        -------
        str
            Mapped EvidenceLevel enum value.
        """
        return CIVIC_LEVEL_MAP.get(civic_level.upper().strip(), "level_4")

    @staticmethod
    def _build_variant_summary(gene: str, variant_name: str, variant_type: str) -> str:
        """Build a minimal text summary for a variant without evidence."""
        parts = [f"{gene} {variant_name}"]
        if variant_type:
            parts.append(f"({variant_type})")
        parts.append("— CIViC variant with no associated clinical evidence items.")
        return " ".join(parts)

    @staticmethod
    def _build_evidence_summary(
        gene: str,
        variant_name: str,
        variant_type: str,
        cancer_type: str,
        drugs: str,
        evidence_type: str,
        evidence_direction: str,
        clinical_significance: str,
        description: str,
    ) -> str:
        """Build a rich text summary for embedding from evidence fields."""
        parts = [f"{gene} {variant_name}"]
        if variant_type:
            parts.append(f"({variant_type})")
        if cancer_type:
            parts.append(f"in {cancer_type}.")
        else:
            parts.append(".")

        if evidence_type:
            parts.append(f"Evidence type: {evidence_type}.")
        if evidence_direction:
            parts.append(f"Direction: {evidence_direction}.")
        if clinical_significance:
            parts.append(f"Clinical significance: {clinical_significance}.")
        if drugs:
            parts.append(f"Associated therapies: {drugs}.")
        if description:
            parts.append(description)

        return " ".join(parts)
