"""ClinicalTrials.gov ingest pipeline for Imaging Intelligence Agent.

Fetches imaging AI clinical trial records via the ClinicalTrials.gov API v2,
parses JSON responses into ImagingTrial models, and stores embeddings
in the imaging_trials Milvus collection.

API v2 docs: https://clinicaltrials.gov/data-api/api

Author: Adam Jones
Date: February 2026
"""

import re
import time
from typing import Any, Dict, List, Optional

import requests
from loguru import logger

from src.models import (
    BodyRegion,
    ImagingModality,
    ImagingTrial,
    TrialPhase,
    TrialStatus,
)

from .base import BaseIngestPipeline


# ClinicalTrials.gov API v2 base URL
CT_GOV_BASE_URL = "https://clinicaltrials.gov/api/v2"

# Default search parameters for imaging AI trials
DEFAULT_CONDITION = "imaging AI"
DEFAULT_INTERVENTION = "artificial intelligence radiology"

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

# Modality keywords for classification from trial text
_MODALITY_KEYWORDS: Dict[ImagingModality, List[str]] = {
    ImagingModality.CT: ["computed tomography", "ct scan", "ct imaging", "low-dose ct", "ldct"],
    ImagingModality.MRI: ["magnetic resonance", "mri", "mr imaging", "flair", "diffusion"],
    ImagingModality.XRAY: ["x-ray", "xray", "radiograph", "plain film"],
    ImagingModality.CXR: ["chest x-ray", "chest radiograph", "cxr"],
    ImagingModality.ULTRASOUND: ["ultrasound", "ultrasonography", "sonography", "doppler"],
    ImagingModality.PET: ["positron emission", "pet scan", "fdg-pet"],
    ImagingModality.PET_CT: ["pet/ct", "pet-ct", "hybrid pet"],
    ImagingModality.MAMMOGRAPHY: ["mammography", "mammogram", "tomosynthesis", "dbt"],
}

# Body region keywords for classification
_REGION_KEYWORDS: Dict[BodyRegion, List[str]] = {
    BodyRegion.HEAD: ["head", "intracranial", "cranial"],
    BodyRegion.BRAIN: ["brain", "cerebral", "neurological", "glioma", "stroke"],
    BodyRegion.CHEST: ["chest", "thorax", "lung", "pulmonary", "pneumonia"],
    BodyRegion.CARDIAC: ["cardiac", "heart", "coronary", "myocardial"],
    BodyRegion.BREAST: ["breast", "mammary", "bi-rads"],
    BodyRegion.ABDOMEN: ["abdominal", "abdomen", "liver", "hepatic", "renal"],
    BodyRegion.PELVIS: ["pelvic", "pelvis", "prostate"],
    BodyRegion.SPINE: ["spine", "spinal", "vertebral"],
    BodyRegion.EXTREMITY: ["extremity", "musculoskeletal", "fracture", "orthopedic"],
    BodyRegion.WHOLE_BODY: ["whole body", "whole-body", "total body"],
}


class ImagingTrialsIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for ClinicalTrials.gov imaging AI trials.

    Fetches trial data via the ClinicalTrials.gov API v2, parses the
    JSON response into ImagingTrial Pydantic models, and stores
    embeddings in the imaging_trials Milvus collection.

    Usage:
        pipeline = ImagingTrialsIngestPipeline(collection_manager, embedder)
        count = pipeline.run(
            condition="imaging AI",
            intervention="artificial intelligence radiology",
            max_results=500,
        )
    """

    COLLECTION_NAME = "imaging_trials"

    def __init__(
        self,
        collection_manager,
        embedder,
        base_url: str = CT_GOV_BASE_URL,
    ):
        """Initialize the ClinicalTrials.gov imaging ingest pipeline.

        Args:
            collection_manager: ImagingCollectionManager for Milvus operations.
            embedder: Embedding model with encode() method.
            base_url: ClinicalTrials.gov API v2 base URL.
        """
        super().__init__(collection_manager, embedder)
        self.base_url = base_url

    def fetch(
        self,
        condition: str = DEFAULT_CONDITION,
        intervention: str = DEFAULT_INTERVENTION,
        max_results: int = 500,
        page_size: int = 100,
    ) -> List[Dict[str, Any]]:
        """Fetch imaging AI clinical trials from ClinicalTrials.gov API v2.

        Uses the GET /studies endpoint with pagination to retrieve all
        matching trials.

        API endpoint: {base_url}/studies
        Query parameters:
            query.cond  -- condition/disease search
            query.intr  -- intervention search
            pageSize    -- results per page (max 1000)
            pageToken   -- pagination cursor

        Args:
            condition: Condition search term (e.g. "imaging AI").
            intervention: Intervention search term
                (e.g. "artificial intelligence radiology").
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
                f"Fetched page {page_num}, total {len(all_studies)} imaging studies so far"
            )

            # Check for next page
            page_token = data.get("nextPageToken")
            if not page_token:
                break

            # Rate-limit: 1 request per second
            time.sleep(1)

        # Trim to exact max_results
        return all_studies[:max_results]

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[ImagingTrial]:
        """Parse ClinicalTrials.gov JSON studies into ImagingTrial models.

        Extracts key fields from the API v2 JSON structure and classifies
        each trial by imaging modality, body region, and AI task.

        Args:
            raw_data: List of study JSON objects from the API.

        Returns:
            List of validated ImagingTrial model instances.
        """
        # AI task detection patterns
        _task_patterns = {
            "detection": re.compile(r"\b(detection|detect|cad|computer.aided)\b", re.IGNORECASE),
            "segmentation": re.compile(r"\b(segmentation|segment|delineation|contour)\b", re.IGNORECASE),
            "classification": re.compile(r"\b(classification|classify|diagnosis|diagnostic)\b", re.IGNORECASE),
            "triage": re.compile(r"\b(triage|prioritiz|worklist|critical finding)\b", re.IGNORECASE),
            "quantification": re.compile(r"\b(quantifi|volumetric|measurement|radiomics)\b", re.IGNORECASE),
            "reconstruction": re.compile(r"\b(reconstruction|super.resolution|denoising)\b", re.IGNORECASE),
        }

        trials: List[ImagingTrial] = []

        for study in raw_data:
            try:
                protocol = study.get("protocolSection", {})
                id_module = protocol.get("identificationModule", {})
                desc_module = protocol.get("descriptionModule", {})
                design_module = protocol.get("designModule", {})
                status_module = protocol.get("statusModule", {})
                sponsor_module = protocol.get("sponsorCollaboratorsModule", {})
                conditions_module = protocol.get("conditionsModule", {})

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
                    year_match = re.match(r"(\d{4})", start_date_str)
                    if year_match:
                        start_year = int(year_match.group(1))

                # --- Disease/conditions ---
                conditions = conditions_module.get("conditions", [])
                disease = "; ".join(conditions[:3]) if conditions else ""
                if len(disease) > 200:
                    disease = disease[:197] + "..."

                # Combined searchable text for classification
                searchable = f"{title} {brief_summary}"

                # --- Modality classification ---
                modality = self._classify_modality(searchable)

                # --- Body region classification ---
                body_region = self._classify_body_region(searchable)

                # --- AI task classification ---
                ai_task = ""
                for task_name, pattern in _task_patterns.items():
                    if pattern.search(searchable):
                        ai_task = task_name
                        break

                trial = ImagingTrial(
                    id=nct_id,
                    title=title[:500] if len(title) > 500 else title,
                    text_summary=text_summary,
                    phase=phase,
                    status=status,
                    sponsor=sponsor_name[:200] if len(sponsor_name) > 200 else sponsor_name,
                    modality=modality,
                    body_region=body_region,
                    ai_task=ai_task,
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
                logger.warning(f"Failed to parse imaging trial {nct}: {exc}")
                continue

        logger.info(f"Parsed {len(trials)} ImagingTrial records from {len(raw_data)} studies")
        return trials

    @staticmethod
    def _extract_phase(phases: Optional[List[str]]) -> TrialPhase:
        """Map ClinicalTrials.gov phase list to a TrialPhase enum.

        Args:
            phases: List of phase strings from the API, or None.

        Returns:
            Corresponding TrialPhase enum value.
        """
        if not phases:
            return TrialPhase.NA
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

    @staticmethod
    def _classify_modality(text: str) -> ImagingModality:
        """Classify a trial's imaging modality by keyword matching.

        Args:
            text: Combined trial title + summary text.

        Returns:
            The ImagingModality with the most keyword matches.
        """
        if not text:
            return ImagingModality.CT

        text_lower = text.lower()
        scores: Dict[ImagingModality, int] = {}

        for modality, keywords in _MODALITY_KEYWORDS.items():
            count = sum(1 for kw in keywords if kw in text_lower)
            scores[modality] = count

        best = max(scores, key=scores.get)  # type: ignore[arg-type]
        return best if scores[best] > 0 else ImagingModality.CT

    @staticmethod
    def _classify_body_region(text: str) -> BodyRegion:
        """Classify a trial's body region by keyword matching.

        Args:
            text: Combined trial title + summary text.

        Returns:
            The BodyRegion with the most keyword matches.
        """
        if not text:
            return BodyRegion.CHEST

        text_lower = text.lower()
        scores: Dict[BodyRegion, int] = {}

        for region, keywords in _REGION_KEYWORDS.items():
            count = sum(1 for kw in keywords if kw in text_lower)
            scores[region] = count

        best = max(scores, key=scores.get)  # type: ignore[arg-type]
        return best if scores[best] > 0 else BodyRegion.CHEST

    def run(
        self,
        collection_name: Optional[str] = None,
        batch_size: int = 32,
        **fetch_kwargs,
    ) -> int:
        """Execute the full ClinicalTrials.gov imaging ingest pipeline.

        Args:
            collection_name: Target collection (defaults to 'imaging_trials').
            batch_size: Batch size for embedding and insertion.
            **fetch_kwargs: Passed to fetch() (condition, intervention,
                max_results, page_size).

        Returns:
            Total number of records ingested.
        """
        target = collection_name or self.COLLECTION_NAME
        logger.info(f"Starting ClinicalTrials.gov imaging ingest pipeline -> {target}")

        raw = self.fetch(**fetch_kwargs)
        logger.info(f"Fetched {len(raw)} raw studies from ClinicalTrials.gov")

        records = self.parse(raw)
        logger.info(f"Parsed {len(records)} ImagingTrial records")

        count = self.embed_and_store(records, target, batch_size)
        logger.info(f"Ingested {count} records into {target}")
        return count
