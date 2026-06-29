"""FDA Adverse Event Reporting System (FAERS) ingest pipeline for CAR-T Intelligence Agent.

Fetches CAR-T adverse event reports via the openFDA drug/event API,
parses JSON responses into SafetyRecord models, and stores embeddings
in the cart_safety Milvus collection.

openFDA API docs: https://open.fda.gov/apis/drug/event/
Rate limits: 240 requests per minute / 120,000 per day (without API key)

Author: Adam Jones
Date: February 2026
"""

import time
from typing import Any, Dict, List, Optional

import requests
from loguru import logger

from src.collections import CARTCollectionManager
from src.models import SafetyEventType, SafetyRecord

from .base import BaseIngestPipeline


# openFDA drug adverse event endpoint
OPENFDA_EVENT_URL = "https://api.fda.gov/drug/event.json"

# Search query for FDA-approved CAR-T products
CART_PRODUCTS_SEARCH = (
    "patient.drug.openfda.brand_name:"
    '("kymriah"+"yescarta"+"tecartus"+"breyanzi"+"abecma"+"carvykti")'
)

# Brand name -> generic/product name mapping for readability
_BRAND_TO_PRODUCT: Dict[str, str] = {
    "kymriah": "Kymriah (tisagenlecleucel)",
    "yescarta": "Yescarta (axicabtagene ciloleucel)",
    "tecartus": "Tecartus (brexucabtagene autoleucel)",
    "breyanzi": "Breyanzi (lisocabtagene maraleucel)",
    "abecma": "Abecma (idecabtagene vicleucel)",
    "carvykti": "Carvykti (ciltacabtagene autoleucel)",
}

# Map openFDA reaction preferred terms to SafetyEventType categories
_REACTION_TO_EVENT_TYPE: Dict[str, SafetyEventType] = {
    "cytokine release syndrome": SafetyEventType.CRS,
    "cytokine storm": SafetyEventType.CRS,
    "immune effector cell-associated neurotoxicity syndrome": SafetyEventType.ICANS,
    "neurotoxicity": SafetyEventType.ICANS,
    "encephalopathy": SafetyEventType.ICANS,
    "confusional state": SafetyEventType.ICANS,
    "aphasia": SafetyEventType.ICANS,
    "tremor": SafetyEventType.NEUROLOGIC,
    "seizure": SafetyEventType.NEUROLOGIC,
    "cerebral oedema": SafetyEventType.NEUROLOGIC,
    "neutropenia": SafetyEventType.CYTOPENIA,
    "thrombocytopenia": SafetyEventType.CYTOPENIA,
    "anaemia": SafetyEventType.CYTOPENIA,
    "pancytopenia": SafetyEventType.CYTOPENIA,
    "febrile neutropenia": SafetyEventType.CYTOPENIA,
    "leukopenia": SafetyEventType.CYTOPENIA,
    "lymphopenia": SafetyEventType.CYTOPENIA,
    "infection": SafetyEventType.INFECTION,
    "sepsis": SafetyEventType.INFECTION,
    "pneumonia": SafetyEventType.INFECTION,
    "bacteraemia": SafetyEventType.INFECTION,
    "fungal infection": SafetyEventType.INFECTION,
    "malignant neoplasm": SafetyEventType.SECONDARY_MALIGNANCY,
    "t-cell lymphoma": SafetyEventType.SECONDARY_MALIGNANCY,
    "myelodysplastic syndrome": SafetyEventType.SECONDARY_MALIGNANCY,
    "acute myeloid leukaemia": SafetyEventType.SECONDARY_MALIGNANCY,
    "cardiac arrest": SafetyEventType.CARDIAC,
    "cardiac failure": SafetyEventType.CARDIAC,
    "tachycardia": SafetyEventType.CARDIAC,
    "atrial fibrillation": SafetyEventType.CARDIAC,
    "hypotension": SafetyEventType.CARDIAC,
    "hepatotoxicity": SafetyEventType.ORGAN_TOXICITY,
    "renal failure": SafetyEventType.ORGAN_TOXICITY,
    "hepatic failure": SafetyEventType.ORGAN_TOXICITY,
    "multi-organ failure": SafetyEventType.ORGAN_TOXICITY,
    "tumour lysis syndrome": SafetyEventType.ORGAN_TOXICITY,
}

# Maximum retry attempts for network requests
MAX_RETRIES = 3

# Delay between requests to stay within rate limits (4 req/sec = 0.25s)
REQUEST_DELAY_SEC = 0.30


class FAERSIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for FDA Adverse Event Reporting System (FAERS) CAR-T data.

    Fetches adverse event reports for FDA-approved CAR-T products via the
    openFDA API, parses them into SafetyRecord models, and stores embeddings
    in the cart_safety Milvus collection.

    This provides real-world pharmacovigilance data complementing the static
    seed safety data from clinical trial labels.

    Usage:
        pipeline = FAERSIngestPipeline(collection_manager, embedder)
        count = pipeline.run(max_results=500)
    """

    COLLECTION_NAME = "cart_safety"

    def __init__(
        self,
        collection_manager: CARTCollectionManager,
        embedder: Any,
        api_key: Optional[str] = None,
    ):
        """Initialize the FAERS ingest pipeline.

        Args:
            collection_manager: CARTCollectionManager for Milvus operations.
            embedder: Embedding model with encode() method.
            api_key: Optional openFDA API key for higher rate limits
                (1000 req/day without key, 120,000/day with key).
        """
        super().__init__(collection_manager, embedder)
        self.api_key = api_key

    def fetch(
        self,
        search: str = CART_PRODUCTS_SEARCH,
        max_results: int = 500,
        page_size: int = 100,
    ) -> List[Dict[str, Any]]:
        """Fetch adverse event reports from the openFDA drug/event API.

        Paginates through results using the skip parameter. openFDA limits
        skip + limit to 26,000 total results.

        API endpoint: https://api.fda.gov/drug/event.json
        Query parameters:
            search — openFDA search query
            limit  — results per page (max 100)
            skip   — pagination offset
            api_key — optional API key

        Args:
            search: openFDA search query string targeting CAR-T brand names.
            max_results: Maximum total number of events to retrieve.
            page_size: Number of events per API request (max 100).

        Returns:
            List of adverse event JSON objects from the API response.
        """
        all_events: List[Dict[str, Any]] = []
        skip = 0
        page_num = 0

        while len(all_events) < max_results:
            params: Dict[str, Any] = {
                "search": search,
                "limit": min(page_size, max_results - len(all_events)),
                "skip": skip,
            }
            if self.api_key:
                params["api_key"] = self.api_key

            response = self._request_with_retry(OPENFDA_EVENT_URL, params)
            if response is None:
                logger.warning(
                    f"FAERS fetch failed after {MAX_RETRIES} retries at skip={skip}; "
                    f"returning {len(all_events)} events collected so far"
                )
                break

            data = response.json()
            results = data.get("results", [])
            if not results:
                logger.info("No more FAERS results returned; pagination complete")
                break

            all_events.extend(results)
            page_num += 1
            skip += len(results)

            logger.info(
                f"FAERS: fetched page {page_num} ({len(results)} events), "
                f"total {len(all_events)} so far"
            )

            # Rate limit: stay under 4 req/sec (240/min)
            time.sleep(REQUEST_DELAY_SEC)

        trimmed = all_events[:max_results]
        logger.info(f"FAERS fetch complete: {len(trimmed)} adverse event reports")
        return trimmed

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[SafetyRecord]:
        """Parse openFDA adverse event JSON into SafetyRecord models.

        Maps the openFDA response structure to SafetyRecord fields:
            - patient.reaction[].reactionmeddrapt -> event_type classification
            - seriousnessdeath / seriousnesshospitalization -> severity_grade
            - receivedate -> year
            - patient.drug[].openfda.brand_name -> product
            - patient.reaction descriptions -> text_summary

        Args:
            raw_data: List of adverse event JSON objects from the openFDA API.

        Returns:
            List of validated SafetyRecord model instances.
        """
        records: List[SafetyRecord] = []

        for idx, event in enumerate(raw_data):
            try:
                record = self._parse_single_event(event, idx)
                if record is not None:
                    records.append(record)
            except Exception as exc:
                safety_id = event.get("safetyreportid", f"index-{idx}")
                logger.warning(f"Failed to parse FAERS event {safety_id}: {exc}")
                continue

        logger.info(
            f"Parsed {len(records)} SafetyRecord instances "
            f"from {len(raw_data)} FAERS events"
        )
        return records

    def _parse_single_event(
        self, event: Dict[str, Any], idx: int
    ) -> Optional[SafetyRecord]:
        """Parse a single openFDA adverse event into a SafetyRecord.

        Args:
            event: Single adverse event JSON object.
            idx: Index in the batch (used for fallback ID generation).

        Returns:
            SafetyRecord instance, or None if the event lacks required data.
        """
        # --- ID ---
        safety_report_id = event.get("safetyreportid", "")
        record_id = f"FAERS-{safety_report_id}" if safety_report_id else f"FAERS-{idx:06d}"

        # --- Patient and drug info ---
        patient = event.get("patient", {})
        drugs = patient.get("drug", [])
        reactions = patient.get("reaction", [])

        if not reactions:
            return None

        # Identify the CAR-T product from the drug list
        product = self._extract_cart_product(drugs)

        # --- Reactions -> event_type and text ---
        reaction_terms = [
            r.get("reactionmeddrapt", "").strip()
            for r in reactions
            if r.get("reactionmeddrapt")
        ]
        if not reaction_terms:
            return None

        event_type = self._classify_event_type(reaction_terms)
        reaction_text = "; ".join(reaction_terms[:10])  # Limit to avoid overflow

        # --- Severity grade ---
        severity_grade = self._extract_severity(event)

        # --- Date (receivedate format: YYYYMMDD) ---
        receive_date = event.get("receivedate", "")
        year = 0
        if receive_date and len(receive_date) >= 4:
            try:
                year = int(receive_date[:4])
                if year < 2010 or year > 2030:
                    year = 0
            except ValueError:
                year = 0

        # --- Build text summary ---
        text_parts = [f"FAERS adverse event report for {product}."]
        text_parts.append(f"Reported reactions: {reaction_text}.")
        text_parts.append(f"Severity: {severity_grade}.")
        if event.get("occurcountry"):
            text_parts.append(f"Country: {event['occurcountry']}.")
        text_summary = " ".join(text_parts)
        if len(text_summary) > 2900:
            text_summary = text_summary[:2897] + "..."

        # --- Onset timing ---
        onset_timing = ""
        for drug in drugs:
            interval = drug.get("drugstartdateformat")
            if interval:
                onset_timing = f"Reported drug start format: {interval}"
                break

        # --- Outcome ---
        outcome_map = {
            "1": "recovered",
            "2": "recovering",
            "3": "not recovered",
            "4": "recovered with sequelae",
            "5": "fatal",
            "6": "unknown",
        }
        reaction_outcome = ""
        if reactions:
            outcome_code = reactions[0].get("reactionoutcome", "")
            reaction_outcome = outcome_map.get(str(outcome_code), "unknown")

        return SafetyRecord(
            id=record_id[:100],
            text_summary=text_summary,
            product=product[:200],
            event_type=event_type,
            severity_grade=severity_grade[:100],
            onset_timing=onset_timing[:100],
            incidence_rate="",  # Not available per individual report
            management_protocol="",  # Not in FAERS individual reports
            outcome=reaction_outcome[:100],
            reporting_source="FAERS",
            year=year,
        )

    @staticmethod
    def _extract_cart_product(drugs: List[Dict[str, Any]]) -> str:
        """Identify the CAR-T product from the drug list in an adverse event.

        Args:
            drugs: List of drug objects from patient.drug[].

        Returns:
            Human-readable product name, or "CAR-T product (unspecified)".
        """
        for drug in drugs:
            openfda = drug.get("openfda", {})
            brand_names = openfda.get("brand_name", [])
            for name in brand_names:
                name_lower = name.lower().strip()
                if name_lower in _BRAND_TO_PRODUCT:
                    return _BRAND_TO_PRODUCT[name_lower]

            # Fallback: check the medicinalproduct field
            med_product = drug.get("medicinalproduct", "").lower().strip()
            if med_product in _BRAND_TO_PRODUCT:
                return _BRAND_TO_PRODUCT[med_product]

        return "CAR-T product (unspecified)"

    @staticmethod
    def _classify_event_type(reaction_terms: List[str]) -> SafetyEventType:
        """Classify adverse event type from MedDRA reaction terms.

        Checks reaction terms against the known mapping. Prioritizes CRS and
        ICANS as the most clinically significant CAR-T-specific events.

        Args:
            reaction_terms: List of MedDRA preferred term strings.

        Returns:
            The most relevant SafetyEventType for this event.
        """
        # Priority order: CRS > ICANS > others
        found_types: List[SafetyEventType] = []

        for term in reaction_terms:
            term_lower = term.lower().strip()
            if term_lower in _REACTION_TO_EVENT_TYPE:
                found_types.append(_REACTION_TO_EVENT_TYPE[term_lower])

        if not found_types:
            return SafetyEventType.ORGAN_TOXICITY  # Generic fallback

        # Prioritize CRS and ICANS
        for priority in [SafetyEventType.CRS, SafetyEventType.ICANS]:
            if priority in found_types:
                return priority

        return found_types[0]

    @staticmethod
    def _extract_severity(event: Dict[str, Any]) -> str:
        """Extract severity grade from openFDA seriousness fields.

        openFDA uses binary (1/absent) flags for seriousness categories:
            - seriousnessdeath
            - seriousnesshospitalization
            - seriousnesslifethreatening
            - seriousnessdisabling
            - seriousnessother

        Args:
            event: Single adverse event JSON object.

        Returns:
            Human-readable severity string (e.g., "fatal", "hospitalization").
        """
        if event.get("seriousnessdeath") == "1":
            return "fatal (Grade 5)"
        if event.get("seriousnesslifethreatening") == "1":
            return "life-threatening (Grade 4)"
        if event.get("seriousnesshospitalization") == "1":
            return "hospitalization required (Grade 3)"
        if event.get("seriousnessdisabling") == "1":
            return "disabling (Grade 3)"
        if event.get("seriousnessother") == "1":
            return "medically significant (Grade 2)"
        if event.get("serious") == "1":
            return "serious (grade unspecified)"
        return "non-serious"

    @staticmethod
    def _request_with_retry(
        url: str,
        params: Dict[str, Any],
        max_retries: int = MAX_RETRIES,
    ) -> Optional[requests.Response]:
        """Make an HTTP GET request with retry logic.

        Args:
            url: Request URL.
            params: Query parameters.
            max_retries: Maximum number of retry attempts.

        Returns:
            requests.Response on success, None on failure.
        """
        for attempt in range(1, max_retries + 1):
            try:
                response = requests.get(
                    url,
                    params=params,
                    headers={"Accept": "application/json"},
                    timeout=30,
                )
                response.raise_for_status()
                return response
            except requests.exceptions.HTTPError as exc:
                status_code = exc.response.status_code if exc.response else 0
                if status_code == 404:
                    logger.info("FAERS API returned 404 (no results for query)")
                    return None
                if status_code == 429:
                    wait = 2 ** attempt
                    logger.warning(
                        f"FAERS rate limited (429); retrying in {wait}s "
                        f"(attempt {attempt}/{max_retries})"
                    )
                    time.sleep(wait)
                    continue
                logger.error(
                    f"FAERS HTTP error {status_code} on attempt "
                    f"{attempt}/{max_retries}: {exc}"
                )
                if attempt < max_retries:
                    time.sleep(2 ** attempt)
                    continue
                return None
            except requests.exceptions.RequestException as exc:
                logger.error(
                    f"FAERS request failed on attempt "
                    f"{attempt}/{max_retries}: {exc}"
                )
                if attempt < max_retries:
                    time.sleep(2 ** attempt)
                    continue
                return None

        return None

    def run(
        self,
        collection_name: Optional[str] = None,
        batch_size: int = 32,
        max_results: int = 500,
        **fetch_kwargs,
    ) -> int:
        """Execute the full FAERS ingest pipeline: fetch -> parse -> embed_and_store.

        Args:
            collection_name: Target Milvus collection (defaults to 'cart_safety').
            batch_size: Batch size for embedding and insertion.
            max_results: Maximum number of adverse events to fetch.
            **fetch_kwargs: Additional keyword arguments passed to fetch().

        Returns:
            Total number of records ingested.
        """
        target = collection_name or self.COLLECTION_NAME
        logger.info(f"Starting FAERS ingest pipeline -> {target}")

        raw = self.fetch(max_results=max_results, **fetch_kwargs)
        logger.info(f"Fetched {len(raw)} raw adverse event reports from FAERS")

        records = self.parse(raw)
        logger.info(f"Parsed {len(records)} SafetyRecord instances")

        count = self.embed_and_store(records, target, batch_size)
        logger.info(f"FAERS ingest complete: {count} records into {target}")
        return count
