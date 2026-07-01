"""ClinicalTrials.gov ingest parser for the Clinical Trial Intelligence Agent.

Fetches and parses clinical trial data from the ClinicalTrials.gov API v2,
extracting structured trial information for embedding and storage in Milvus.

Rate limiting: 3 requests/second without an API key.

API Reference: https://clinicaltrials.gov/data-api/api

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
import time
from typing import Any, Dict, List, Optional

from .base import BaseIngestParser, IngestRecord

logger = logging.getLogger(__name__)

# ===================================================================
# CONSTANTS
# ===================================================================

CTGOV_API_BASE = "https://clinicaltrials.gov/api/v2/studies"
DEFAULT_MAX_RESULTS = 200
DEFAULT_PAGE_SIZE = 20
RATE_LIMIT_DELAY = 0.34  # ~3 requests/second

# Default conditions to search when no specific query is provided
DEFAULT_CONDITIONS = [
    "cancer",
    "heart failure",
    "alzheimer",
    "diabetes",
    "rare disease",
]

# ===================================================================
# LANDMARK TRIALS
# ===================================================================

LANDMARK_TRIALS: List[Dict[str, Any]] = [
    {
        "nct_id": "NCT02477436",
        "title": "KEYNOTE-024: Pembrolizumab vs Platinum-Based Chemotherapy in PD-L1+ NSCLC",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "oncology",
        "conditions": ["Non-Small Cell Lung Cancer"],
        "interventions": ["Pembrolizumab", "Platinum-Based Chemotherapy"],
        "sponsor": "Merck Sharp & Dohme LLC",
        "significance": "Established pembrolizumab as first-line monotherapy for PD-L1-high NSCLC",
    },
    {
        "nct_id": "NCT02813135",
        "title": "KEYNOTE-189: Pembrolizumab + Chemotherapy in Metastatic NSCLC",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "oncology",
        "conditions": ["Non-Small Cell Lung Cancer"],
        "interventions": ["Pembrolizumab", "Pemetrexed", "Platinum Chemotherapy"],
        "sponsor": "Merck Sharp & Dohme LLC",
        "significance": "Established chemo-immunotherapy combination as standard of care",
    },
    {
        "nct_id": "NCT01673867",
        "title": "CheckMate 067: Nivolumab +/- Ipilimumab in Advanced Melanoma",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "oncology",
        "conditions": ["Advanced Melanoma"],
        "interventions": ["Nivolumab", "Ipilimumab"],
        "sponsor": "Bristol-Myers Squibb",
        "significance": "Demonstrated durable long-term survival with combination checkpoint inhibition",
    },
    {
        "nct_id": "NCT01820754",
        "title": "PARADIGM-HF: Sacubitril/Valsartan in Heart Failure",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "cardiology",
        "conditions": ["Heart Failure with Reduced Ejection Fraction"],
        "interventions": ["Sacubitril/Valsartan", "Enalapril"],
        "sponsor": "Novartis",
        "significance": "Established ARNI as superior to ACEi in HFrEF management",
    },
    {
        "nct_id": "NCT03036124",
        "title": "DAPA-HF: Dapagliflozin in Heart Failure",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "cardiology",
        "conditions": ["Heart Failure"],
        "interventions": ["Dapagliflozin"],
        "sponsor": "AstraZeneca",
        "significance": "First SGLT2i approved for HFrEF regardless of diabetes status",
    },
    {
        "nct_id": "NCT03057977",
        "title": "EMPEROR-Reduced: Empagliflozin in Heart Failure",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "cardiology",
        "conditions": ["Heart Failure with Reduced Ejection Fraction"],
        "interventions": ["Empagliflozin"],
        "sponsor": "Boehringer Ingelheim",
        "significance": "Confirmed SGLT2i class effect in HFrEF",
    },
    {
        "nct_id": "NCT02465268",
        "title": "FOURIER: Evolocumab Cardiovascular Outcomes",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "cardiology",
        "conditions": ["Atherosclerotic Cardiovascular Disease"],
        "interventions": ["Evolocumab"],
        "sponsor": "Amgen",
        "significance": "Demonstrated PCSK9 inhibition reduces cardiovascular events",
    },
    {
        "nct_id": "NCT02412098",
        "title": "ELOQUENT-3: Elotuzumab + Pomalidomide/Dexamethasone in Myeloma",
        "phase": "Phase II",
        "status": "Completed",
        "therapeutic_area": "hematology",
        "conditions": ["Multiple Myeloma"],
        "interventions": ["Elotuzumab", "Pomalidomide", "Dexamethasone"],
        "sponsor": "Bristol-Myers Squibb",
        "significance": "Novel immunostimulatory approach in relapsed/refractory myeloma",
    },
    {
        "nct_id": "NCT02370498",
        "title": "ELIANA: Tisagenlecleucel in Pediatric ALL",
        "phase": "Phase II",
        "status": "Completed",
        "therapeutic_area": "hematology",
        "conditions": ["Acute Lymphoblastic Leukemia"],
        "interventions": ["Tisagenlecleucel (CAR-T)"],
        "sponsor": "Novartis",
        "significance": "First FDA-approved CAR-T therapy",
    },
    {
        "nct_id": "NCT04368728",
        "title": "RECOVERY: Dexamethasone in Hospitalized COVID-19",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "infectious_disease",
        "conditions": ["COVID-19"],
        "interventions": ["Dexamethasone"],
        "sponsor": "University of Oxford",
        "significance": "First treatment proven to reduce COVID-19 mortality",
    },
    {
        "nct_id": "NCT04280705",
        "title": "ACTT-1: Remdesivir for COVID-19",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "infectious_disease",
        "conditions": ["COVID-19"],
        "interventions": ["Remdesivir"],
        "sponsor": "NIAID",
        "significance": "First FDA-approved antiviral for COVID-19",
    },
    {
        "nct_id": "NCT03548935",
        "title": "VICTORIA: Vericiguat in Heart Failure",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "cardiology",
        "conditions": ["Heart Failure"],
        "interventions": ["Vericiguat"],
        "sponsor": "Merck/Bayer",
        "significance": "Novel sGC stimulator for worsening HF",
    },
    {
        "nct_id": "NCT02507959",
        "title": "MONARCH 3: Abemaciclib + AI in HR+ Breast Cancer",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "oncology",
        "conditions": ["HR-Positive Breast Cancer"],
        "interventions": ["Abemaciclib", "Aromatase Inhibitor"],
        "sponsor": "Eli Lilly",
        "significance": "CDK4/6 inhibitor as first-line therapy in HR+/HER2- advanced breast cancer",
    },
    {
        "nct_id": "NCT03521986",
        "title": "DESTINY-Breast03: Trastuzumab Deruxtecan vs T-DM1 in HER2+ Breast Cancer",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "oncology",
        "conditions": ["HER2-Positive Breast Cancer"],
        "interventions": ["Trastuzumab Deruxtecan", "T-DM1"],
        "sponsor": "Daiichi Sankyo/AstraZeneca",
        "significance": "Next-generation ADC showed superior efficacy over T-DM1",
    },
    {
        "nct_id": "NCT03299946",
        "title": "IMpower150: Atezolizumab + Bevacizumab + Chemo in NSCLC",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "oncology",
        "conditions": ["Non-Small Cell Lung Cancer"],
        "interventions": ["Atezolizumab", "Bevacizumab", "Carboplatin", "Paclitaxel"],
        "sponsor": "Hoffmann-La Roche",
        "significance": "Four-drug regimen effective across PD-L1 expression levels",
    },
    {
        "nct_id": "NCT02684006",
        "title": "SPRINT: Intensive vs Standard Blood Pressure Control",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "cardiology",
        "conditions": ["Hypertension"],
        "interventions": ["Intensive BP Control (<120 mmHg)", "Standard BP Control (<140 mmHg)"],
        "sponsor": "NHLBI",
        "significance": "Intensive BP control reduces cardiovascular events and mortality",
    },
    {
        "nct_id": "NCT02853331",
        "title": "DECLARE-TIMI 58: Dapagliflozin Cardiovascular Outcomes in T2D",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "metabolic",
        "conditions": ["Type 2 Diabetes Mellitus"],
        "interventions": ["Dapagliflozin"],
        "sponsor": "AstraZeneca",
        "significance": "SGLT2i cardiovascular and renal benefits in type 2 diabetes",
    },
    {
        "nct_id": "NCT02564263",
        "title": "MONALEESA-7: Ribociclib in Premenopausal HR+/HER2- Breast Cancer",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "oncology",
        "conditions": ["HR-Positive Breast Cancer"],
        "interventions": ["Ribociclib", "Endocrine Therapy"],
        "sponsor": "Novartis",
        "significance": "First CDK4/6i to show OS benefit in premenopausal advanced breast cancer",
    },
    {
        "nct_id": "NCT02362594",
        "title": "PACIFIC: Durvalumab after Chemoradiotherapy in Stage III NSCLC",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "oncology",
        "conditions": ["Stage III Non-Small Cell Lung Cancer"],
        "interventions": ["Durvalumab"],
        "sponsor": "AstraZeneca",
        "significance": "Established consolidation immunotherapy after chemoradiation as SOC",
    },
    {
        "nct_id": "NCT03170960",
        "title": "CREDENCE: Canagliflozin and Renal Events in Diabetes",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "metabolic",
        "conditions": ["Diabetic Kidney Disease"],
        "interventions": ["Canagliflozin"],
        "sponsor": "Janssen",
        "significance": "First trial to show SGLT2i renal protection in diabetic nephropathy",
    },
    {
        "nct_id": "NCT03723655",
        "title": "ADAURA: Osimertinib in Early-Stage EGFR+ NSCLC",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "oncology",
        "conditions": ["EGFR-Mutant Non-Small Cell Lung Cancer"],
        "interventions": ["Osimertinib"],
        "sponsor": "AstraZeneca",
        "significance": "Adjuvant EGFR TKI dramatically improved DFS in early-stage NSCLC",
    },
    {
        "nct_id": "NCT03834506",
        "title": "HIMALAYA: Tremelimumab + Durvalumab in Hepatocellular Carcinoma",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "oncology",
        "conditions": ["Hepatocellular Carcinoma"],
        "interventions": ["Tremelimumab", "Durvalumab"],
        "sponsor": "AstraZeneca",
        "significance": "STRIDE regimen established as first-line option in unresectable HCC",
    },
    {
        "nct_id": "NCT03734029",
        "title": "DESTINY-Breast04: Trastuzumab Deruxtecan in HER2-Low Metastatic Breast Cancer",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "oncology",
        "conditions": ["HER2-Low Metastatic Breast Cancer"],
        "interventions": ["Trastuzumab Deruxtecan (T-DXd)", "Physician's Choice Chemotherapy"],
        "sponsor": "Daiichi Sankyo/AstraZeneca",
        "significance": "Established T-DXd efficacy in HER2-low mBC, PFS 9.9 vs 5.1 months",
    },
    {
        "nct_id": "NCT03036488",
        "title": "KEYNOTE-522: Pembrolizumab + Chemotherapy as Neoadjuvant Therapy in TNBC",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "oncology",
        "conditions": ["Triple-Negative Breast Cancer"],
        "interventions": ["Pembrolizumab", "Neoadjuvant Chemotherapy"],
        "sponsor": "Merck Sharp & Dohme LLC",
        "significance": "Neoadjuvant pembro+chemo achieved pCR 64.8% in TNBC, establishing immunotherapy in early-stage TNBC",
    },
    {
        "nct_id": "NCT03215706",
        "title": "CheckMate-9LA: Nivolumab + Ipilimumab + Chemotherapy in First-Line NSCLC",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "oncology",
        "conditions": ["Non-Small Cell Lung Cancer"],
        "interventions": ["Nivolumab", "Ipilimumab", "Platinum-Doublet Chemotherapy"],
        "sponsor": "Bristol-Myers Squibb",
        "significance": "Dual checkpoint inhibition + limited chemo improved OS 15.8 vs 11.0 months in first-line NSCLC",
    },
    {
        "nct_id": "NCT02511106",
        "title": "ADAURA: Osimertinib as Adjuvant Therapy in EGFR-Mutant NSCLC",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "oncology",
        "conditions": ["EGFR-Mutant Non-Small Cell Lung Cancer"],
        "interventions": ["Osimertinib"],
        "sponsor": "AstraZeneca",
        "significance": "Adjuvant osimertinib achieved DFS HR 0.17 in resected EGFR+ NSCLC, unprecedented benefit",
    },
    {
        "nct_id": "NCT03574597",
        "title": "SELECT: Semaglutide 2.4 mg Cardiovascular Outcomes in Obesity",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "metabolic",
        "conditions": ["Obesity", "Cardiovascular Disease"],
        "interventions": ["Semaglutide 2.4 mg"],
        "sponsor": "Novo Nordisk",
        "significance": "Semaglutide 2.4 mg reduced MACE by 20% in overweight/obese adults without diabetes",
    },
    {
        "nct_id": "NCT04184622",
        "title": "SURMOUNT-1: Tirzepatide for Treatment of Obesity",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "metabolic",
        "conditions": ["Obesity"],
        "interventions": ["Tirzepatide"],
        "sponsor": "Eli Lilly",
        "significance": "Tirzepatide achieved 22.5% mean weight loss, highest ever in a Phase III obesity trial",
    },
    {
        "nct_id": "NCT04437511",
        "title": "TRAILBLAZER-ALZ 2: Donanemab in Early Symptomatic Alzheimer's Disease",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "neurology",
        "conditions": ["Early Alzheimer's Disease"],
        "interventions": ["Donanemab"],
        "sponsor": "Eli Lilly",
        "significance": "Donanemab slowed cognitive and functional decline by 35% in early AD",
    },
    {
        "nct_id": "NCT04994509",
        "title": "PURPOSE 1: Lenacapavir for HIV Pre-Exposure Prophylaxis",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "infectious_disease",
        "conditions": ["HIV Prevention"],
        "interventions": ["Lenacapavir (twice-yearly injection)"],
        "sponsor": "Gilead Sciences",
        "significance": "Lenacapavir PrEP achieved 100% HIV prevention efficacy with twice-yearly dosing",
    },
    {
        "nct_id": "NCT04223856",
        "title": "EV-302: Enfortumab Vedotin + Pembrolizumab in Urothelial Carcinoma",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "oncology",
        "conditions": ["Urothelial Carcinoma"],
        "interventions": ["Enfortumab Vedotin", "Pembrolizumab"],
        "sponsor": "Astellas/Seagen",
        "significance": "EV+pembro doubled PFS (12.5 vs 6.3 months) in first-line urothelial carcinoma",
    },
    {
        "nct_id": "NCT02435849",
        "title": "ELIANA: Tisagenlecleucel in Pediatric and Young Adult ALL",
        "phase": "Phase II",
        "status": "Completed",
        "therapeutic_area": "hematology",
        "conditions": ["Pediatric Acute Lymphoblastic Leukemia"],
        "interventions": ["Tisagenlecleucel (CAR-T)"],
        "sponsor": "Novartis",
        "significance": "Tisagenlecleucel achieved 82% complete remission in pediatric r/r ALL, leading to first CAR-T approval",
    },
    {
        "nct_id": "NCT03391466",
        "title": "ZUMA-7: Axicabtagene Ciloleucel vs Standard of Care in Second-Line DLBCL",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "hematology",
        "conditions": ["Diffuse Large B-Cell Lymphoma"],
        "interventions": ["Axicabtagene Ciloleucel (Axi-cel)", "Standard of Care"],
        "sponsor": "Kite/Gilead",
        "significance": "Axi-cel EFS 8.3 vs 2.0 months compared to SOC in second-line DLBCL",
    },
    {
        "nct_id": "NCT03900078",
        "title": "MAESTRO: Resmetirom in Non-Alcoholic Steatohepatitis (NASH)",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "metabolic",
        "conditions": ["Non-Alcoholic Steatohepatitis"],
        "interventions": ["Resmetirom"],
        "sponsor": "Madrigal Pharmaceuticals",
        "significance": "Resmetirom became the first FDA-approved drug specifically for NASH",
    },
    {
        "nct_id": "NCT04194450",
        "title": "PROTECT: Teplizumab for Delay of Type 1 Diabetes Onset",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "metabolic",
        "conditions": ["Type 1 Diabetes Mellitus"],
        "interventions": ["Teplizumab"],
        "sponsor": "Provention Bio/Sanofi",
        "significance": "Teplizumab delayed clinical T1DM onset by approximately 2 years in at-risk individuals",
    },
    {
        "nct_id": "NCT04349072",
        "title": "VALOR-HCM: Mavacamten in Obstructive Hypertrophic Cardiomyopathy",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "cardiology",
        "conditions": ["Obstructive Hypertrophic Cardiomyopathy"],
        "interventions": ["Mavacamten"],
        "sponsor": "Bristol-Myers Squibb",
        "significance": "Mavacamten reduced need for septal reduction therapy in symptomatic obstructive HCM",
    },
    {
        "nct_id": "NCT04788511",
        "title": "TRILUMINATE Pivotal: TriClip for Tricuspid Regurgitation",
        "phase": "Phase III",
        "status": "Completed",
        "therapeutic_area": "cardiology",
        "conditions": ["Tricuspid Regurgitation"],
        "interventions": ["TriClip Transcatheter Tricuspid Valve Repair"],
        "sponsor": "Abbott",
        "significance": "TriClip demonstrated significant quality-of-life improvement in severe tricuspid regurgitation",
    },
]


# ===================================================================
# PARSER IMPLEMENTATION
# ===================================================================


class ClinicalTrialsParser(BaseIngestParser):
    """Ingest parser for ClinicalTrials.gov API v2.

    Fetches clinical trial study records, extracts structured fields
    (trial_id, title, phase, status, conditions, interventions,
    eligibility, sites, endpoints), and produces IngestRecord objects
    for embedding and Milvus storage.

    Rate limiting is enforced at 3 requests/second without an API key.

    Usage::

        parser = ClinicalTrialsParser()
        records, stats = parser.run(
            conditions=["breast cancer", "NSCLC"],
            max_results=100,
        )
    """

    def __init__(
        self,
        collection_manager: Any = None,
        embedder: Any = None,
        api_key: Optional[str] = None,
    ) -> None:
        """Initialize the ClinicalTrials.gov parser.

        Args:
            collection_manager: Optional Milvus collection manager.
            embedder: Optional embedding model.
            api_key: Optional ClinicalTrials.gov API key for higher rate limits.
        """
        super().__init__(
            source_name="clinicaltrials",
            collection_manager=collection_manager,
            embedder=embedder,
        )
        self.api_key = api_key
        self._last_request_time: float = 0.0

    def _rate_limit(self) -> None:
        """Enforce rate limiting between API requests."""
        elapsed = time.time() - self._last_request_time
        if elapsed < RATE_LIMIT_DELAY:
            time.sleep(RATE_LIMIT_DELAY - elapsed)
        self._last_request_time = time.time()

    def fetch(self, **kwargs: Any) -> List[Dict[str, Any]]:
        """Fetch clinical trial data from ClinicalTrials.gov API v2.

        Args:
            conditions: List of condition terms to search.
            max_results: Maximum number of trials to fetch (default 200).
            statuses: List of trial statuses to filter by.
            phases: List of trial phases to filter by.

        Returns:
            List of raw study dictionaries from the API.
        """
        try:
            import requests
        except ImportError:
            self.logger.warning(
                "requests library not available -- returning landmark trials only"
            )
            return self.seed_landmark_trials()

        conditions = kwargs.get("conditions", DEFAULT_CONDITIONS)
        max_results = kwargs.get("max_results", DEFAULT_MAX_RESULTS)
        statuses = kwargs.get("statuses", [])
        phases = kwargs.get("phases", [])

        all_studies: List[Dict[str, Any]] = []

        for condition in conditions:
            if len(all_studies) >= max_results:
                break

            params: Dict[str, Any] = {
                "query.cond": condition,
                "pageSize": min(DEFAULT_PAGE_SIZE, max_results - len(all_studies)),
                "format": "json",
                "fields": (
                    "NCTId,BriefTitle,OfficialTitle,Phase,OverallStatus,"
                    "Condition,InterventionName,InterventionType,"
                    "EligibilityCriteria,MinimumAge,MaximumAge,Gender,"
                    "LocationFacility,LocationCity,LocationCountry,"
                    "PrimaryOutcomeMeasure,SecondaryOutcomeMeasure,"
                    "LeadSponsorName,EnrollmentCount,StartDate,"
                    "PrimaryCompletionDate,StudyType"
                ),
            }

            if statuses:
                params["filter.overallStatus"] = "|".join(statuses)
            if phases:
                params["filter.phase"] = "|".join(phases)

            headers: Dict[str, str] = {}
            if self.api_key:
                headers["X-API-Key"] = self.api_key

            try:
                self._rate_limit()
                response = requests.get(
                    CTGOV_API_BASE,
                    params=params,
                    headers=headers,
                    timeout=30,
                )
                response.raise_for_status()
                data = response.json()

                studies = data.get("studies", [])
                all_studies.extend(studies)
                self.logger.info(
                    "Fetched %d studies for condition '%s'",
                    len(studies), condition,
                )

            except requests.exceptions.RequestException as exc:
                self.logger.error(
                    "Failed to fetch studies for '%s': %s", condition, exc
                )
                continue

        return all_studies[:max_results]

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[IngestRecord]:
        """Parse raw ClinicalTrials.gov study data into IngestRecord objects.

        Args:
            raw_data: List of raw study dictionaries from the API.

        Returns:
            List of IngestRecord objects.
        """
        records: List[IngestRecord] = []

        for study in raw_data:
            try:
                record = self._parse_study(study)
                if record:
                    records.append(record)
            except Exception as exc:
                self.logger.warning("Failed to parse study: %s", exc)
                continue

        return records

    def _parse_study(self, study: Dict[str, Any]) -> Optional[IngestRecord]:
        """Parse a single study dictionary into an IngestRecord.

        Handles both the API v2 nested format and the flat landmark trial format.
        """
        # Handle landmark trial format (flat dict)
        if "nct_id" in study:
            return self._parse_landmark(study)

        # API v2 format: nested under protocolSection
        protocol = study.get("protocolSection", {})
        id_module = protocol.get("identificationModule", {})
        status_module = protocol.get("statusModule", {})
        design_module = protocol.get("designModule", {})
        conditions_module = protocol.get("conditionsModule", {})
        interventions_module = protocol.get("armsInterventionsModule", {})
        eligibility_module = protocol.get("eligibilityModule", {})
        contacts_module = protocol.get("contactsLocationsModule", {})
        outcomes_module = protocol.get("outcomesModule", {})
        sponsor_module = protocol.get("sponsorCollaboratorsModule", {})

        nct_id = id_module.get("nctId", "")
        if not nct_id:
            return None

        title = id_module.get("officialTitle") or id_module.get("briefTitle", "")
        phase = design_module.get("phases", ["Not Applicable"])
        status = status_module.get("overallStatus", "Unknown")
        conditions = conditions_module.get("conditions", [])

        # Extract interventions
        interventions = []
        for arm in interventions_module.get("interventions", []):
            name = arm.get("name", "")
            itype = arm.get("type", "")
            if name:
                interventions.append(f"{name} ({itype})" if itype else name)

        # Eligibility
        eligibility_text = eligibility_module.get("eligibilityCriteria", "")
        min_age = eligibility_module.get("minimumAge", "")
        max_age = eligibility_module.get("maximumAge", "")
        sex = eligibility_module.get("sex", "All")

        # Sites
        sites = []
        for location in contacts_module.get("locations", []):
            facility = location.get("facility", "")
            city = location.get("city", "")
            country = location.get("country", "")
            if facility or city:
                sites.append(f"{facility}, {city}, {country}".strip(", "))

        # Endpoints
        primary_endpoints = []
        for outcome in outcomes_module.get("primaryOutcomes", []):
            measure = outcome.get("measure", "")
            if measure:
                primary_endpoints.append(measure)

        secondary_endpoints = []
        for outcome in outcomes_module.get("secondaryOutcomes", []):
            measure = outcome.get("measure", "")
            if measure:
                secondary_endpoints.append(measure)

        # Sponsor
        lead_sponsor = ""
        sponsor_info = sponsor_module.get("leadSponsor", {})
        if sponsor_info:
            lead_sponsor = sponsor_info.get("name", "")

        # Compose text for embedding
        text_parts = [
            f"Trial: {title}",
            f"NCT ID: {nct_id}",
            f"Phase: {', '.join(phase) if isinstance(phase, list) else phase}",
            f"Status: {status}",
            f"Conditions: {', '.join(conditions)}",
            f"Interventions: {', '.join(interventions)}",
        ]
        if primary_endpoints:
            text_parts.append(f"Primary Endpoints: {', '.join(primary_endpoints)}")
        if eligibility_text:
            # Truncate long eligibility text
            elig_preview = eligibility_text[:500]
            text_parts.append(f"Eligibility: {elig_preview}")

        text = "\n".join(text_parts)

        metadata = {
            "nct_id": nct_id,
            "title": title,
            "phase": phase if isinstance(phase, list) else [phase],
            "status": status,
            "conditions": conditions,
            "interventions": interventions,
            "eligibility_text": eligibility_text[:1000] if eligibility_text else "",
            "min_age": min_age,
            "max_age": max_age,
            "sex": sex,
            "sites": sites[:20],  # Cap at 20 sites
            "primary_endpoints": primary_endpoints,
            "secondary_endpoints": secondary_endpoints,
            "sponsor": lead_sponsor,
            "source": "clinicaltrials.gov",
        }

        return IngestRecord(
            text=text,
            metadata=metadata,
            collection_name="trial_protocols",
            record_id=nct_id,
            source="clinicaltrials",
        )

    def _parse_landmark(self, trial: Dict[str, Any]) -> IngestRecord:
        """Parse a landmark trial dictionary into an IngestRecord."""
        nct_id = trial["nct_id"]
        title = trial.get("title", "")
        phase = trial.get("phase", "")
        status = trial.get("status", "")
        conditions = trial.get("conditions", [])
        interventions = trial.get("interventions", [])
        significance = trial.get("significance", "")
        sponsor = trial.get("sponsor", "")
        therapeutic_area = trial.get("therapeutic_area", "")

        text_parts = [
            f"Landmark Trial: {title}",
            f"NCT ID: {nct_id}",
            f"Phase: {phase}",
            f"Status: {status}",
            f"Conditions: {', '.join(conditions)}",
            f"Interventions: {', '.join(interventions)}",
            f"Significance: {significance}",
        ]
        text = "\n".join(text_parts)

        metadata = {
            "nct_id": nct_id,
            "title": title,
            "phase": [phase],
            "status": status,
            "conditions": conditions,
            "interventions": interventions,
            "significance": significance,
            "sponsor": sponsor,
            "therapeutic_area": therapeutic_area,
            "is_landmark": True,
            "source": "curated_landmark",
        }

        return IngestRecord(
            text=text,
            metadata=metadata,
            collection_name="trial_protocols",
            record_id=nct_id,
            source="clinicaltrials",
        )

    def validate_record(self, record: IngestRecord) -> bool:
        """Validate a clinical trial IngestRecord.

        Checks that required fields are present and text meets minimum length.

        Args:
            record: The IngestRecord to validate.

        Returns:
            True if valid.
        """
        if not record.text or len(record.text.strip()) < 20:
            return False
        if not record.record_id:
            return False
        if not record.metadata.get("title"):
            return False
        return True

    def seed_landmark_trials(self) -> List[Dict[str, Any]]:
        """Return the curated list of landmark clinical trials.

        These are pre-defined landmark trials that represent pivotal
        studies across multiple therapeutic areas. Used for initial
        knowledge base seeding.

        Returns:
            List of landmark trial dictionaries.
        """
        return LANDMARK_TRIALS.copy()
