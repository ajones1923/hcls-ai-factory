"""ClinicalTrials.gov ingest pipeline for CAR-T Intelligence Agent.

Fetches CAR-T clinical trial records via the ClinicalTrials.gov API v2,
parses JSON responses into ClinicalTrial models, and stores embeddings
in the cart_trials Milvus collection.

API v2 docs: https://clinicaltrials.gov/data-api/api

Author: Adam Jones
Date: February 2026
"""

import re
import time
from typing import Any, Dict, List, Optional

import requests
from loguru import logger

from src.collections import CARTCollectionManager
from src.models import (
    CARGeneration,
    ClinicalTrial,
    TrialPhase,
    TrialStatus,
)

from .base import BaseIngestPipeline


# ClinicalTrials.gov API v2 base URL
CT_GOV_BASE_URL = "https://clinicaltrials.gov/api/v2"

# Default search parameters for CAR-T trials
DEFAULT_CONDITION = "CAR-T"
DEFAULT_INTERVENTION = "chimeric antigen receptor"

# Mapping from API phase strings to TrialPhase enum
_PHASE_MAP: Dict[str, TrialPhase] = {
    "EARLY_PHASE1": TrialPhase.EARLY_1,
    "PHASE1": TrialPhase.PHASE_1,
    "PHASE1_PHASE2": TrialPhase.PHASE_1_2,
    "PHASE2": TrialPhase.PHASE_2,
    "PHASE2_PHASE3": TrialPhase.PHASE_2_3,
    "PHASE3": TrialPhase.PHASE_3,
    "PHASE4": TrialPhase.PHASE_4,
    "NA": TrialPhase.NA,
}

# Mapping from API status strings to TrialStatus enum
_STATUS_MAP: Dict[str, TrialStatus] = {
    "RECRUITING": TrialStatus.RECRUITING,
    "ACTIVE_NOT_RECRUITING": TrialStatus.ACTIVE,
    "COMPLETED": TrialStatus.COMPLETED,
    "TERMINATED": TrialStatus.TERMINATED,
    "WITHDRAWN": TrialStatus.WITHDRAWN,
    "SUSPENDED": TrialStatus.SUSPENDED,
    "NOT_YET_RECRUITING": TrialStatus.NOT_YET,
    "UNKNOWN": TrialStatus.UNKNOWN,
}


class ClinicalTrialsIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for ClinicalTrials.gov CAR-T trials.

    Fetches trial data via the ClinicalTrials.gov API v2, parses the
    JSON response into ClinicalTrial Pydantic models, and stores
    embeddings in the cart_trials Milvus collection.

    Usage:
        pipeline = ClinicalTrialsIngestPipeline(collection_manager, embedder)
        count = pipeline.run(
            condition="CAR-T",
            intervention="chimeric antigen receptor",
            max_results=1000,
        )
    """

    COLLECTION_NAME = "cart_trials"

    def __init__(
        self,
        collection_manager: CARTCollectionManager,
        embedder: Any,
        base_url: str = CT_GOV_BASE_URL,
    ):
        """Initialize the ClinicalTrials.gov ingest pipeline.

        Args:
            collection_manager: CARTCollectionManager for Milvus operations.
            embedder: Embedding model with encode() method.
            base_url: ClinicalTrials.gov API v2 base URL.
        """
        super().__init__(collection_manager, embedder)
        self.base_url = base_url

    def fetch(
        self,
        condition: str = DEFAULT_CONDITION,
        intervention: str = DEFAULT_INTERVENTION,
        max_results: int = 1000,
        page_size: int = 100,
    ) -> List[Dict[str, Any]]:
        """Fetch CAR-T clinical trials from ClinicalTrials.gov API v2.

        Uses the GET /studies endpoint with pagination to retrieve all
        matching trials.

        API endpoint: {base_url}/studies
        Query parameters:
            query.cond  — condition/disease search
            query.intr  — intervention search
            pageSize    — results per page (max 1000)
            pageToken   — pagination cursor

        Args:
            condition: Condition search term (e.g. "CAR-T").
            intervention: Intervention search term
                (e.g. "chimeric antigen receptor").
            max_results: Maximum total number of studies to retrieve.
            page_size: Number of studies per API request (max 1000).

        Returns:
            List of study JSON objects from the API response.
        """
        url = f"{self.base_url}/studies"
        all_studies: List[Dict[str, Any]] = []
        page_token: Optional[str] = None
        page_num = 0

        while len(all_studies) < max_results:
            params: Dict[str, Any] = {
                "query.cond": condition,
                "query.intr": intervention,
                "pageSize": min(page_size, max_results - len(all_studies)),
            }
            if page_token:
                params["pageToken"] = page_token

            response = requests.get(
                url,
                params=params,
                headers={"Accept": "application/json"},
                timeout=30,
            )
            response.raise_for_status()
            data = response.json()

            studies = data.get("studies", [])
            all_studies.extend(studies)
            page_num += 1
            logger.info(
                f"Fetched page {page_num}, total {len(all_studies)} studies so far"
            )

            # Check for next page
            page_token = data.get("nextPageToken")
            if not page_token:
                break

            # Rate-limit: 1 request per second
            time.sleep(1)

        # Trim to exact max_results
        return all_studies[:max_results]

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[ClinicalTrial]:
        """Parse ClinicalTrials.gov JSON studies into ClinicalTrial models.

        Extracts key fields from the API v2 JSON structure:
          - protocolSection.identificationModule.nctId
          - protocolSection.identificationModule.officialTitle
          - protocolSection.descriptionModule.briefSummary
          - protocolSection.designModule.phases
          - protocolSection.statusModule.overallStatus
          - protocolSection.sponsorCollaboratorsModule.leadSponsor.name
          - protocolSection.designModule.enrollmentInfo.count
          - protocolSection.statusModule.startDateStruct.date

        Args:
            raw_data: List of study JSON objects from the API.

        Returns:
            List of validated ClinicalTrial model instances.
        """
        # Known CAR-T target antigens for extraction
        _ANTIGEN_PATTERN = re.compile(
            r"\b(CD19|CD20|CD22|CD30|CD33|CD38|CD123|CD138|CD171|"
            r"BCMA|B-cell maturation antigen|GD2|HER2|EGFR|EGFRvIII|"
            r"GPC3|Mesothelin|PSMA|MUC1|MUC16|ROR1|Lewis[- ]?Y|"
            r"CS1|SLAMF7|GPRC5D|FLT3|CLL-1|NKG2D|IL13R[Aa]2|"
            r"Claudin18\\.2|CLDN18)\b",
            re.IGNORECASE,
        )
        # CAR generation detection
        _GENERATION_PATTERNS = {
            CARGeneration.FOURTH: re.compile(
                r"\b(4th[- ]gen|fourth[- ]gen|armored CAR|TRUCKs?|4G\b)", re.IGNORECASE
            ),
            CARGeneration.THIRD: re.compile(
                r"\b(3rd[- ]gen|third[- ]gen|3G\b)", re.IGNORECASE
            ),
            CARGeneration.UNIVERSAL: re.compile(
                r"\b(universal|allogeneic|off[- ]the[- ]shelf)\b", re.IGNORECASE
            ),
            CARGeneration.ARMORED: re.compile(
                r"\b(armored|cytokine[- ]secreting)\b", re.IGNORECASE
            ),
        }
        # Costimulatory domain detection
        _COSTIM_PATTERN = re.compile(
            r"\b(CD28|4-1BB|CD137|ICOS|OX40|CD27)\b", re.IGNORECASE
        )

        trials: List[ClinicalTrial] = []

        for study in raw_data:
            try:
                protocol = study.get("protocolSection", {})
                id_module = protocol.get("identificationModule", {})
                desc_module = protocol.get("descriptionModule", {})
                design_module = protocol.get("designModule", {})
                status_module = protocol.get("statusModule", {})
                sponsor_module = protocol.get("sponsorCollaboratorsModule", {})
                conditions_module = protocol.get("conditionsModule", {})
                arms_module = protocol.get("armsInterventionsModule", {})

                # --- Required fields ---
                nct_id = id_module.get("nctId", "")
                if not nct_id:
                    logger.warning("Skipping study with missing nctId")
                    continue

                title = (
                    id_module.get("officialTitle")
                    or id_module.get("briefTitle")
                    or "Untitled"
                )
                brief_summary = desc_module.get("briefSummary", "")
                text_summary = f"{title}. {brief_summary}".strip()
                # Truncate to model max_length
                if len(text_summary) > 2900:
                    text_summary = text_summary[:2897] + "..."

                # --- Phase & Status ---
                phases = design_module.get("phases")
                phase = self._extract_phase(phases)

                overall_status = status_module.get("overallStatus")
                status = self._extract_status(overall_status)

                # --- Sponsor ---
                lead_sponsor = sponsor_module.get("leadSponsor", {})
                sponsor_name = lead_sponsor.get("name", "")

                # --- Enrollment ---
                enrollment_info = design_module.get("enrollmentInfo", {})
                enrollment = enrollment_info.get("count", 0) or 0

                # --- Start year ---
                start_date_struct = status_module.get("startDateStruct", {})
                start_date_str = start_date_struct.get("date", "")
                start_year = 0
                if start_date_str:
                    # Format is typically "YYYY-MM-DD" or "YYYY-MM" or "YYYY"
                    year_match = re.match(r"(\d{4})", start_date_str)
                    if year_match:
                        start_year = int(year_match.group(1))

                # --- Disease/conditions ---
                conditions = conditions_module.get("conditions", [])
                disease = "; ".join(conditions[:3]) if conditions else ""
                if len(disease) > 200:
                    disease = disease[:197] + "..."

                # --- Interventions text for downstream extraction ---
                interventions = arms_module.get("interventions", [])
                intervention_text = " ".join(
                    intv.get("description", "") + " " + intv.get("name", "")
                    for intv in interventions
                )

                # Combined searchable text for antigen/generation extraction
                searchable = f"{title} {brief_summary} {intervention_text}"

                # --- Target antigen extraction ---
                antigen_match = _ANTIGEN_PATTERN.search(searchable)
                target_antigen = antigen_match.group(0) if antigen_match else ""
                # Normalize BCMA
                if target_antigen.lower() == "b-cell maturation antigen":
                    target_antigen = "BCMA"

                # --- CAR generation extraction ---
                car_generation = CARGeneration.SECOND  # default
                for gen, pattern in _GENERATION_PATTERNS.items():
                    if pattern.search(searchable):
                        car_generation = gen
                        break

                # --- Costimulatory domain ---
                costim_matches = _COSTIM_PATTERN.findall(searchable)
                # Normalize CD137 -> 4-1BB
                costim_set = set()
                for c in costim_matches:
                    if c.upper() == "CD137":
                        costim_set.add("4-1BB")
                    else:
                        costim_set.add(c)
                costimulatory = ", ".join(sorted(costim_set)) if costim_set else ""
                if len(costimulatory) > 50:
                    costimulatory = costimulatory[:47] + "..."

                trial = ClinicalTrial(
                    id=nct_id,
                    title=title[:500] if len(title) > 500 else title,
                    text_summary=text_summary,
                    phase=phase,
                    status=status,
                    sponsor=sponsor_name[:200] if len(sponsor_name) > 200 else sponsor_name,
                    target_antigen=target_antigen[:100],
                    car_generation=car_generation,
                    costimulatory=costimulatory,
                    disease=disease,
                    enrollment=enrollment,
                    start_year=start_year if 2000 <= start_year <= 2030 else 0,
                    outcome_summary="",
                )
                trials.append(trial)

            except Exception as exc:
                nct = study.get("protocolSection", {}).get(
                    "identificationModule", {}
                ).get("nctId", "unknown")
                logger.warning(f"Failed to parse study {nct}: {exc}")
                continue

        logger.info(f"Parsed {len(trials)} ClinicalTrial records from {len(raw_data)} studies")
        return trials

    @staticmethod
    def _extract_phase(phases: Optional[List[str]]) -> TrialPhase:
        """Map ClinicalTrials.gov phase list to a TrialPhase enum.

        The API returns phases as a list (e.g. ["PHASE1", "PHASE2"]).
        This method maps the combined value to the most specific TrialPhase.

        Args:
            phases: List of phase strings from the API, or None.

        Returns:
            Corresponding TrialPhase enum value.
        """
        if not phases:
            return TrialPhase.NA

        # Join multi-phase into combined key (e.g. "PHASE1_PHASE2")
        combined = "_".join(phases)
        return _PHASE_MAP.get(combined, TrialPhase.NA)

    @staticmethod
    def _extract_status(status_str: Optional[str]) -> TrialStatus:
        """Map ClinicalTrials.gov status string to a TrialStatus enum.

        Args:
            status_str: Status string from the API (e.g. "RECRUITING").

        Returns:
            Corresponding TrialStatus enum value.
        """
        if not status_str:
            return TrialStatus.UNKNOWN
        return _STATUS_MAP.get(status_str.upper(), TrialStatus.UNKNOWN)

    def run(
        self,
        collection_name: Optional[str] = None,
        batch_size: int = 32,
        **fetch_kwargs,
    ) -> int:
        """Execute the full ClinicalTrials.gov ingest pipeline.

        Args:
            collection_name: Target collection (defaults to 'cart_trials').
            batch_size: Batch size for embedding and insertion.
            **fetch_kwargs: Passed to fetch() (condition, intervention,
                max_results, page_size).

        Returns:
            Total number of records ingested.

        """
        target = collection_name or self.COLLECTION_NAME
        logger.info(f"Starting ClinicalTrials.gov ingest pipeline -> {target}")

        raw = self.fetch(**fetch_kwargs)
        logger.info(f"Fetched {len(raw)} raw studies from ClinicalTrials.gov")

        records = self.parse(raw)
        logger.info(f"Parsed {len(records)} ClinicalTrial records")

        count = self.embed_and_store(records, target, batch_size)
        logger.info(f"Ingested {count} records into {target}")
        return count
