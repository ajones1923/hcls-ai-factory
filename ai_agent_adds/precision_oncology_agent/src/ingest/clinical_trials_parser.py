"""
ClinicalTrials.gov oncology trial ingest pipeline.

Fetches precision oncology clinical trials from the ClinicalTrials.gov v2
API and normalises them into the ``onco_trials`` collection schema for the
Oncology Intelligence Agent's RAG knowledge base.

Follows the same structural pattern as the CAR-T Agent's clinical trials
ingest for consistency across the HCLS AI Factory project.

API reference: https://clinicaltrials.gov/data-api/api

Author: Adam Jones
Date: February 2026
"""

import logging
import re
import time
from typing import Any, Dict, List, Optional

import requests

from src.ingest.base import BaseIngestPipeline

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CT_API_BASE = "https://clinicaltrials.gov/api/v2"
CT_STUDIES_ENDPOINT = f"{CT_API_BASE}/studies"

REQUEST_TIMEOUT = 30
PAGE_SIZE = 100
RATE_LIMIT_DELAY = 0.2

# Biomarker / gene patterns to extract from eligibility criteria
BIOMARKER_KEYWORDS = [
    "EGFR", "BRAF", "KRAS", "NRAS", "ALK", "ROS1", "MET", "HER2",
    "ERBB2", "RET", "NTRK", "FGFR", "IDH1", "IDH2", "BRCA", "PIK3CA",
    "TP53", "PD-L1", "TMB", "MSI", "MSI-H", "dMMR", "PTEN", "CDK",
    "FLT3", "NPM1", "JAK2", "BCR-ABL", "KIT", "PDGFRA",
]

_biomarker_pattern = re.compile(
    r"\b(" + "|".join(re.escape(b) for b in BIOMARKER_KEYWORDS) + r")\b",
    re.IGNORECASE,
)


class ClinicalTrialsIngestPipeline(BaseIngestPipeline):
    """
    Ingest pipeline for ClinicalTrials.gov oncology trials.

    Populates the ``onco_trials`` Milvus collection with trial summaries,
    eligibility criteria, biomarker requirements, and status information
    for matching against patient genomic profiles.
    """

    def __init__(self, collection_manager: Any, embedder: Any) -> None:
        super().__init__(
            collection_manager=collection_manager,
            embedder=embedder,
            collection_name="onco_trials",
        )

    # ------------------------------------------------------------------
    # Fetch
    # ------------------------------------------------------------------

    def fetch(
        self,
        query: str = "precision oncology targeted therapy biomarker",
        max_results: int = 2000,
    ) -> List[Dict]:
        """
        Fetch oncology clinical trial records from ClinicalTrials.gov v2 API.

        Parameters
        ----------
        query : str
            Search query for trial matching (default covers broad
            precision oncology terms).
        max_results : int
            Maximum number of trial records to retrieve (default 2000).

        Returns
        -------
        list of dict
            Raw study records from ClinicalTrials.gov.
        """
        studies: List[Dict] = []
        next_page_token: Optional[str] = None

        logger.info(
            "Fetching clinical trials: query=%r, max_results=%d",
            query, max_results,
        )

        while len(studies) < max_results:
            try:
                params: Dict[str, Any] = {
                    "query.term": query,
                    "filter.overallStatus": "RECRUITING,ENROLLING_BY_INVITATION,ACTIVE_NOT_RECRUITING",
                    "pageSize": min(PAGE_SIZE, max_results - len(studies)),
                    "fields": (
                        "NCTId,BriefTitle,OfficialTitle,OverallStatus,Phase,"
                        "BriefSummary,Condition,InterventionName,InterventionType,"
                        "EligibilityCriteria,LeadSponsorName,EnrollmentCount,"
                        "StartDate,PrimaryCompletionDate,StudyType"
                    ),
                }
                if next_page_token:
                    params["pageToken"] = next_page_token

                response = requests.get(
                    CT_STUDIES_ENDPOINT,
                    params=params,
                    timeout=REQUEST_TIMEOUT,
                )
                response.raise_for_status()
                data = response.json()

            except requests.RequestException as exc:
                logger.error("ClinicalTrials.gov API request failed: %s", exc)
                break

            study_list = data.get("studies", [])
            if not study_list:
                logger.info("No more studies returned — stopping.")
                break

            studies.extend(study_list)

            # Pagination
            next_page_token = data.get("nextPageToken")
            if not next_page_token:
                break

            time.sleep(RATE_LIMIT_DELAY)

        logger.info("Fetched %d trial records from ClinicalTrials.gov", len(studies))
        return studies[:max_results]

    # ------------------------------------------------------------------
    # Parse
    # ------------------------------------------------------------------

    def parse(self, raw_data: List[Dict]) -> List[Dict]:
        """
        Parse ClinicalTrials.gov study records into ``onco_trials`` schema.

        Schema fields:
            - id: NCT number (e.g. "NCT04000000")
            - title: Study title
            - text_summary: Brief summary (primary embedding source)
            - text: Same as text_summary (alias for base class)
            - phase: Study phase (e.g. "Phase 2")
            - status: Overall status (e.g. "Recruiting")
            - sponsor: Lead sponsor name
            - cancer_types: Comma-separated condition list
            - biomarker_criteria: Extracted biomarker mentions from eligibility
            - enrollment: Target enrollment count
            - start_year: Year the study started
            - interventions: Intervention names
            - source_type: Always "clinicaltrials"

        Parameters
        ----------
        raw_data : list of dict
            Raw study records from ``fetch``.

        Returns
        -------
        list of dict
            Normalised records ready for embedding and insertion.
        """
        records: List[Dict] = []

        for study in raw_data:
            # ClinicalTrials.gov v2 nests data under protocolSection
            protocol = study.get("protocolSection", study)
            identification = protocol.get("identificationModule", {})
            status_module = protocol.get("statusModule", {})
            description = protocol.get("descriptionModule", {})
            design = protocol.get("designModule", {})
            eligibility = protocol.get("eligibilityModule", {})
            sponsor_module = protocol.get("sponsorCollaboratorsModule", {})
            conditions_module = protocol.get("conditionsModule", {})
            interventions_module = protocol.get("armsInterventionsModule", {})

            # NCT ID
            nct_id = identification.get("nctId", "")
            if not nct_id:
                continue

            # Title
            title = (
                identification.get("officialTitle")
                or identification.get("briefTitle", "")
            )

            # Summary
            brief_summary = description.get("briefSummary", "")
            text_summary = f"{title}. {brief_summary}" if brief_summary else title

            # Phase
            phases = design.get("phases", [])
            phase = ", ".join(phases) if phases else ""

            # Status
            status = status_module.get("overallStatus", "")

            # Sponsor
            lead_sponsor = sponsor_module.get("leadSponsor", {})
            sponsor = lead_sponsor.get("name", "")

            # Cancer types / conditions
            conditions = conditions_module.get("conditions", [])
            cancer_types = ", ".join(conditions) if conditions else ""

            # Eligibility criteria text
            eligibility_text = eligibility.get("eligibilityCriteria", "")

            # Extract biomarker criteria
            biomarker_criteria = self._extract_biomarkers(eligibility_text)

            # Enrollment
            enrollment_info = design.get("enrollmentInfo", {})
            enrollment = enrollment_info.get("count", "")

            # Start year
            start_date = status_module.get("startDateStruct", {})
            start_year = start_date.get("date", "")
            if start_year and len(str(start_year)) >= 4:
                start_year = str(start_year)[:4]

            # Interventions
            interventions_list = interventions_module.get("interventions", [])
            interventions = ", ".join(
                i.get("name", "") for i in interventions_list if i.get("name")
            ) if interventions_list else ""

            records.append({
                "id": nct_id,
                "title": title,
                "text_summary": text_summary,
                "text": text_summary,
                "phase": phase,
                "status": status,
                "sponsor": sponsor,
                "cancer_types": cancer_types,
                "biomarker_criteria": biomarker_criteria,
                "enrollment": str(enrollment),
                "start_year": str(start_year),
                "interventions": interventions,
                "source_type": "clinicaltrials",
            })

        logger.info(
            "Parsed %d trial records from %d raw studies",
            len(records), len(raw_data),
        )
        return records

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _extract_biomarkers(eligibility_text: str) -> str:
        """Extract unique biomarker/gene mentions from eligibility criteria."""
        if not eligibility_text:
            return ""
        matches = _biomarker_pattern.findall(eligibility_text)
        unique = list(dict.fromkeys(m.upper() for m in matches))
        return ", ".join(unique) if unique else ""
