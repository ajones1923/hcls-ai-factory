"""ClinicalTrials.gov PGx trial ingest pipeline for PGx Intelligence Agent.

Fetches pharmacogenomics-related clinical trials from the ClinicalTrials.gov
v2 API (https://clinicaltrials.gov/api/v2/studies), parses trial metadata
and summaries into PGxClinicalTrial models, and stores embeddings in the
pgx_clinical_trials Milvus collection.

API docs: https://clinicaltrials.gov/data-api/api

Author: Adam Jones
Date: March 2026
"""

import re
from typing import Any, Dict, List, Optional

import requests
from loguru import logger
from pydantic import BaseModel

from src.collections import PGxCollectionManager
from src.models import PGxClinicalTrial

from .base import BaseIngestPipeline, _truncate_utf8


# ClinicalTrials.gov v2 API
CT_GOV_BASE_URL = "https://clinicaltrials.gov/api/v2/studies"

# Default PGx search query for ClinicalTrials.gov
DEFAULT_CT_QUERY = (
    "pharmacogenomics OR pharmacogenetics OR genotype-guided OR "
    "CYP2D6 OR CYP2C19 OR CYP2C9 OR DPYD OR TPMT OR "
    "precision prescribing OR pharmacogenomic testing"
)

# Core PGx gene patterns for extraction
_PGX_GENE_REGEX = re.compile(
    r"\b(?:CYP2D6|CYP2C19|CYP2C9|CYP3A5|CYP2B6|CYP3A4|CYP1A2|"
    r"DPYD|TPMT|NUDT15|UGT1A1|SLCO1B1|VKORC1|G6PD|HLA-[AB]|"
    r"NAT2|ABCG2|IFNL3|CYP4F2)\b",
    re.IGNORECASE,
)

# Drug patterns for extraction from trial descriptions
_DRUG_REGEX = re.compile(
    r"\b(?:warfarin|clopidogrel|tamoxifen|codeine|tramadol|simvastatin|"
    r"abacavir|carbamazepine|phenytoin|fluorouracil|capecitabine|"
    r"mercaptopurine|azathioprine|tacrolimus|voriconazole|irinotecan|"
    r"allopurinol|omeprazole|sertraline|escitalopram|amitriptyline|"
    r"nortriptyline|ondansetron|atomoxetine)\b",
    re.IGNORECASE,
)

# Phase normalization
_PHASE_MAP = {
    "EARLY_PHASE1": "Early Phase 1",
    "PHASE1": "Phase 1",
    "PHASE2": "Phase 2",
    "PHASE3": "Phase 3",
    "PHASE4": "Phase 4",
    "NA": "Not Applicable",
}


class ClinicalTrialsPGxParser(BaseIngestPipeline):
    """Ingest pipeline for ClinicalTrials.gov PGx studies.

    Fetches PGx-related clinical trials from the ClinicalTrials.gov v2 API,
    parses trial metadata, extracts gene and drug mentions, and stores
    results in the pgx_clinical_trials Milvus collection.

    Usage:
        parser = ClinicalTrialsPGxParser(collection_manager, embedder)
        count = parser.run(query="CYP2D6 genotype-guided", max_results=200)
    """

    COLLECTION_NAME = "pgx_clinical_trials"

    def __init__(
        self,
        collection_manager: PGxCollectionManager,
        embedder: Any,
        request_timeout: int = 30,
    ):
        super().__init__(collection_manager, embedder)
        self.request_timeout = request_timeout

    def fetch(
        self,
        query: str = DEFAULT_CT_QUERY,
        max_results: int = 500,
        **kwargs,
    ) -> List[Dict[str, Any]]:
        """Fetch PGx clinical trials from ClinicalTrials.gov v2 API.

        Paginates through results using the nextPageToken mechanism.

        Returns:
            List of study dicts from the API response.
        """
        all_studies: List[Dict[str, Any]] = []
        page_size = min(max_results, 100)
        next_page_token: Optional[str] = None

        while len(all_studies) < max_results:
            try:
                params: Dict[str, Any] = {
                    "query.term": query,
                    "pageSize": page_size,
                    "format": "json",
                    "fields": (
                        "NCTId,BriefTitle,OfficialTitle,BriefSummary,"
                        "OverallStatus,Phase,EnrollmentCount,StartDate,"
                        "Condition,InterventionName,PrimaryOutcomeMeasure,"
                        "ResultsFirstPostDate"
                    ),
                }
                if next_page_token:
                    params["pageToken"] = next_page_token

                logger.info(f"Fetching ClinicalTrials.gov page (have {len(all_studies)} studies)...")
                resp = requests.get(
                    CT_GOV_BASE_URL,
                    params=params,
                    timeout=self.request_timeout,
                )
                resp.raise_for_status()
                data = resp.json()

                studies = data.get("studies", [])
                if not studies:
                    break

                all_studies.extend(studies)
                next_page_token = data.get("nextPageToken")

                if not next_page_token:
                    break

                logger.info(f"Fetched {len(studies)} studies (total: {len(all_studies)})")

            except requests.RequestException as e:
                logger.error(f"Failed to fetch ClinicalTrials.gov studies: {e}")
                break

        all_studies = all_studies[:max_results]
        logger.info(f"Fetched {len(all_studies)} PGx clinical trials total")
        return all_studies

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[PGxClinicalTrial]:
        """Parse ClinicalTrials.gov study dicts into PGxClinicalTrial models."""
        records: List[PGxClinicalTrial] = []

        for study in raw_data:
            try:
                # v2 API nests data under protocolSection
                protocol = study.get("protocolSection", study)
                identification = protocol.get("identificationModule", {})
                status_module = protocol.get("statusModule", {})
                design_module = protocol.get("designModule", {})
                description_module = protocol.get("descriptionModule", {})
                conditions_module = protocol.get("conditionsModule", {})
                interventions_module = protocol.get("armsInterventionsModule", {})
                outcomes_module = protocol.get("outcomesModule", {})

                nct_id = identification.get("nctId", "")
                brief_title = identification.get("briefTitle", "")
                official_title = identification.get("officialTitle", brief_title)
                brief_summary = description_module.get("briefSummary", "")

                overall_status = status_module.get("overallStatus", "")

                # Phase
                phases = design_module.get("phases", [])
                phase_raw = phases[0] if phases else ""
                phase = _PHASE_MAP.get(phase_raw, phase_raw)

                # Enrollment
                enrollment_info = design_module.get("enrollmentInfo", {})
                enrollment = 0
                if isinstance(enrollment_info, dict):
                    try:
                        enrollment = int(enrollment_info.get("count", 0))
                    except (TypeError, ValueError):
                        enrollment = 0

                # Start year
                start_date_struct = status_module.get("startDateStruct", {})
                start_date_str = start_date_struct.get("date", "") if isinstance(start_date_struct, dict) else ""
                start_year = 0
                if start_date_str:
                    try:
                        start_year = int(start_date_str[:4])
                    except (ValueError, IndexError):
                        start_year = 0

                # Extract gene and drug from title + summary
                combined_text = f"{brief_title} {official_title} {brief_summary}".strip()
                gene = self._extract_gene(combined_text)
                drug = self._extract_drug(combined_text)

                # Interventions
                interventions = interventions_module.get("interventions", [])
                intervention_names = []
                for interv in interventions:
                    if isinstance(interv, dict):
                        name = interv.get("name", "")
                        if name:
                            intervention_names.append(name)

                # Primary outcomes
                primary_outcomes = outcomes_module.get("primaryOutcomes", [])
                outcome_parts = []
                for outcome in primary_outcomes:
                    if isinstance(outcome, dict):
                        measure = outcome.get("measure", "")
                        if measure:
                            outcome_parts.append(measure)
                outcome_summary = "; ".join(outcome_parts)

                # Build text summary for embedding
                text_summary = (
                    f"{brief_title}. {brief_summary}. "
                    f"Interventions: {', '.join(intervention_names)}. "
                    f"Status: {overall_status}. Phase: {phase}."
                ).strip()

                record = PGxClinicalTrial(
                    id=_truncate_utf8(f"ct_{nct_id}", 95),
                    title=_truncate_utf8(brief_title, 490),
                    text_summary=_truncate_utf8(text_summary, 2990),
                    nct_id=_truncate_utf8(nct_id, 18),
                    phase=_truncate_utf8(phase, 48),
                    status=_truncate_utf8(overall_status, 48),
                    gene=_truncate_utf8(gene, 48),
                    drug=_truncate_utf8(drug, 195),
                    enrollment=enrollment,
                    start_year=start_year,
                    outcome_summary=_truncate_utf8(outcome_summary, 1990),
                )
                records.append(record)

            except Exception as e:
                nct = study.get("protocolSection", {}).get("identificationModule", {}).get("nctId", "?")
                logger.warning(f"Failed to parse trial {nct}: {e}")
                continue

        logger.info(f"Parsed {len(records)} PGx clinical trial records")
        return records

    @staticmethod
    def _extract_gene(text: str) -> str:
        """Extract the first PGx gene mentioned in text."""
        if not text:
            return ""
        match = _PGX_GENE_REGEX.search(text)
        return match.group(0).upper() if match else ""

    @staticmethod
    def _extract_drug(text: str) -> str:
        """Extract the first PGx-relevant drug mentioned in text."""
        if not text:
            return ""
        match = _DRUG_REGEX.search(text)
        return match.group(0).lower() if match else ""
