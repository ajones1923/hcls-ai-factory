"""FDA pharmacogenomic biomarker label ingest pipeline for PGx Intelligence Agent.

Parses the FDA Table of Pharmacogenomic Biomarkers in Drug Labeling,
which lists drugs with PGx information in their FDA-approved labels.
Data source: https://www.fda.gov/drugs/science-and-research-drugs/
table-pharmacogenomic-biomarkers-drug-labeling

Author: Adam Jones
Date: March 2026
"""

from typing import Any, Dict, List, Optional

import requests
from loguru import logger
from pydantic import BaseModel

from src.collections import PGxCollectionManager
from src.models import FDALabel

from .base import BaseIngestPipeline, _truncate_utf8


# FDA Table of Pharmacogenomic Biomarkers endpoint
# The FDA provides a downloadable table; we use the openFDA API for drug labels
OPENFDA_DRUG_LABEL_URL = "https://api.fda.gov/drug/label.json"
FDA_PGX_TABLE_URL = "https://api.fda.gov/drug/drugsfda.json"

# Known biomarker genes from FDA PGx table for targeted queries
FDA_PGX_GENES = [
    "CYP2D6", "CYP2C19", "CYP2C9", "CYP3A5", "CYP2B6",
    "DPYD", "TPMT", "NUDT15", "UGT1A1", "G6PD",
    "HLA-B", "HLA-A", "VKORC1", "IFNL3", "SLCO1B1",
    "CYP1A2", "NAT2", "RAS", "EGFR", "BRAF",
    "ALK", "BRCA", "PD-L1", "ERBB2", "KRAS",
]

# Map FDA labeling sections to standardized names
_SECTION_MAP = {
    "indications_and_usage": "Indications and Usage",
    "dosage_and_administration": "Dosage and Administration",
    "warnings_and_precautions": "Warnings and Precautions",
    "contraindications": "Contraindications",
    "clinical_pharmacology": "Clinical Pharmacology",
    "use_in_specific_populations": "Use in Specific Populations",
    "boxed_warning": "Boxed Warning",
}

# Keywords for classifying label type
_LABEL_TYPE_KEYWORDS = {
    "testing required": ["required", "must test", "test required", "mandatory testing"],
    "actionable PGx": ["dose adjustment", "reduce dose", "alternative", "avoid"],
    "informative PGx": ["may affect", "may influence", "polymorphism", "metabolism"],
}

_REQUIREMENT_KEYWORDS = {
    "required": ["required", "must", "mandatory", "contraindicated"],
    "recommended": ["recommended", "should", "consider testing"],
    "informational": ["may", "can", "information"],
}


class FDALabelParser(BaseIngestPipeline):
    """Ingest pipeline for FDA pharmacogenomic drug labeling.

    Queries the openFDA API for drug labels containing pharmacogenomic
    biomarker information and parses them into FDALabel models for the
    pgx_fda_labels collection.

    Usage:
        parser = FDALabelParser(collection_manager, embedder)
        count = parser.run()
    """

    COLLECTION_NAME = "pgx_fda_labels"

    def __init__(
        self,
        collection_manager: PGxCollectionManager,
        embedder: Any,
        request_timeout: int = 30,
    ):
        super().__init__(collection_manager, embedder)
        self.request_timeout = request_timeout

    def fetch(self, **kwargs) -> List[Dict[str, Any]]:
        """Fetch FDA drug labels with pharmacogenomic content from openFDA.

        Queries openFDA for labels mentioning known PGx biomarker genes.

        Returns:
            List of openFDA drug label result dicts.
        """
        genes = kwargs.get("genes", FDA_PGX_GENES)
        max_per_gene = kwargs.get("max_per_gene", 50)
        all_labels: List[Dict[str, Any]] = []
        seen_ids: set = set()

        for gene in genes:
            try:
                search_query = f'"{gene}"'
                params = {
                    "search": search_query,
                    "limit": min(max_per_gene, 100),
                }
                logger.info(f"Fetching FDA labels for {gene}...")
                resp = requests.get(
                    OPENFDA_DRUG_LABEL_URL,
                    params=params,
                    timeout=self.request_timeout,
                )
                resp.raise_for_status()
                data = resp.json()

                results = data.get("results", [])
                for result in results:
                    label_id = result.get("id", "")
                    if label_id and label_id not in seen_ids:
                        result["_query_gene"] = gene
                        all_labels.append(result)
                        seen_ids.add(label_id)

                logger.info(f"Fetched {len(results)} FDA labels for {gene}")

            except requests.RequestException as e:
                logger.error(f"Failed to fetch FDA labels for {gene}: {e}")
                continue

        logger.info(f"Fetched {len(all_labels)} unique FDA labels across {len(genes)} genes")
        return all_labels

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[FDALabel]:
        """Parse openFDA label results into FDALabel models."""
        records: List[FDALabel] = []

        for label in raw_data:
            try:
                label_id = label.get("id", "")
                query_gene = label.get("_query_gene", "")

                # Extract drug name from openfda metadata
                openfda = label.get("openfda", {})
                brand_names = openfda.get("brand_name", [])
                generic_names = openfda.get("generic_name", [])
                drug_name = (
                    brand_names[0] if brand_names
                    else generic_names[0] if generic_names
                    else "Unknown"
                )

                # Search all relevant label sections for PGx content
                for section_key, section_display in _SECTION_MAP.items():
                    section_texts = label.get(section_key, [])
                    if not section_texts:
                        continue

                    section_text = " ".join(section_texts) if isinstance(section_texts, list) else str(section_texts)

                    # Only keep sections that mention the query gene
                    if query_gene.lower() not in section_text.lower():
                        continue

                    label_type = self._classify_label_type(section_text)
                    requirement_level = self._classify_requirement(section_text)

                    # Extract last updated date from effective_time
                    effective_time = label.get("effective_time", "")
                    last_updated = ""
                    if effective_time and len(effective_time) >= 8:
                        last_updated = f"{effective_time[:4]}-{effective_time[4:6]}-{effective_time[6:8]}"

                    text_chunk = (
                        f"{drug_name} FDA label {section_display} for {query_gene}. "
                        f"{section_text}"
                    ).strip()

                    rec_id = _truncate_utf8(
                        f"fda_{label_id}_{query_gene}_{section_key}"[:95], 95
                    )

                    record = FDALabel(
                        id=rec_id,
                        drug=_truncate_utf8(drug_name, 195),
                        gene=_truncate_utf8(query_gene, 48),
                        labeling_section=_truncate_utf8(section_display, 195),
                        text_chunk=_truncate_utf8(text_chunk, 2990),
                        label_type=_truncate_utf8(label_type, 95),
                        requirement_level=_truncate_utf8(requirement_level, 95),
                        last_updated=last_updated[:10] if last_updated else "",
                    )
                    records.append(record)

            except Exception as e:
                logger.warning(f"Failed to parse FDA label {label.get('id', '?')}: {e}")
                continue

        logger.info(f"Parsed {len(records)} FDA label records")
        return records

    @staticmethod
    def _classify_label_type(text: str) -> str:
        """Classify the FDA label type from section text."""
        text_lower = text.lower()
        for label_type, keywords in _LABEL_TYPE_KEYWORDS.items():
            if any(kw in text_lower for kw in keywords):
                return label_type
        return "informative PGx"

    @staticmethod
    def _classify_requirement(text: str) -> str:
        """Classify the requirement level from section text."""
        text_lower = text.lower()
        for level, keywords in _REQUIREMENT_KEYWORDS.items():
            if any(kw in text_lower for kw in keywords):
                return level
        return "informational"
