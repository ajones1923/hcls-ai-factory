"""DailyMed / FDA drug label ingest pipeline for CAR-T Intelligence Agent.

Fetches Structured Product Labeling (SPL) data for FDA-approved CAR-T products
via the DailyMed REST API, parses JSON responses into RegulatoryRecord models,
and stores embeddings in the cart_regulatory Milvus collection.

DailyMed API docs: https://dailymed.nlm.nih.gov/dailymed/app-support-web-services.cfm

Author: Adam Jones
Date: February 2026
"""

import time
from typing import Any, Dict, List, Optional

import requests
from loguru import logger

from src.collections import CARTCollectionManager
from src.models import RegulatoryEvent, RegulatoryRecord

from .base import BaseIngestPipeline


# DailyMed SPL search endpoint
DAILYMED_SPL_URL = "https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json"

# DailyMed SPL detail endpoint (for individual label content)
DAILYMED_SPL_DETAIL_URL = "https://dailymed.nlm.nih.gov/dailymed/services/v2/spls/{set_id}.json"

# FDA-approved CAR-T products to query
CART_PRODUCT_NAMES: List[str] = [
    "kymriah",
    "yescarta",
    "tecartus",
    "breyanzi",
    "abecma",
    "carvykti",
]

# Mapping from product name to structured metadata
_PRODUCT_METADATA: Dict[str, Dict[str, str]] = {
    "kymriah": {
        "generic": "tisagenlecleucel",
        "sponsor": "Novartis",
        "target": "CD19",
        "initial_indication": "r/r B-cell ALL (pediatric/young adult), r/r DLBCL",
    },
    "yescarta": {
        "generic": "axicabtagene ciloleucel",
        "sponsor": "Kite/Gilead",
        "target": "CD19",
        "initial_indication": "r/r large B-cell lymphoma",
    },
    "tecartus": {
        "generic": "brexucabtagene autoleucel",
        "sponsor": "Kite/Gilead",
        "target": "CD19",
        "initial_indication": "r/r mantle cell lymphoma",
    },
    "breyanzi": {
        "generic": "lisocabtagene maraleucel",
        "sponsor": "Bristol Myers Squibb",
        "target": "CD19",
        "initial_indication": "r/r large B-cell lymphoma",
    },
    "abecma": {
        "generic": "idecabtagene vicleucel",
        "sponsor": "Bristol Myers Squibb",
        "target": "BCMA",
        "initial_indication": "r/r multiple myeloma",
    },
    "carvykti": {
        "generic": "ciltacabtagene autoleucel",
        "sponsor": "Janssen/Legend Biotech",
        "target": "BCMA",
        "initial_indication": "r/r multiple myeloma",
    },
}

# Static fallback seed data when the DailyMed API is unavailable
_FALLBACK_SEED_DATA: List[Dict[str, Any]] = [
    {
        "product": "Kymriah (tisagenlecleucel)",
        "set_id": "kymriah-spl-fallback",
        "title": "KYMRIAH - tisagenlecleucel suspension for intravenous infusion",
        "effective_time": "2024-03",
        "text_summary": (
            "Kymriah (tisagenlecleucel) is a CD19-directed genetically modified autologous "
            "T cell immunotherapy. FDA approved August 2017 for r/r B-cell ALL in pediatric "
            "and young adult patients, and May 2018 for r/r DLBCL after two or more lines "
            "of systemic therapy. REMS required for CRS and neurological toxicity management. "
            "Boxed warning for CRS and neurological toxicities. Label includes dosing, "
            "preparation, administration, warnings, and adverse reactions data."
        ),
        "regulatory_event": "full_approval",
        "indication": "r/r B-cell ALL (pediatric/young adult), r/r DLBCL",
    },
    {
        "product": "Yescarta (axicabtagene ciloleucel)",
        "set_id": "yescarta-spl-fallback",
        "title": "YESCARTA - axicabtagene ciloleucel suspension for intravenous infusion",
        "effective_time": "2024-06",
        "text_summary": (
            "Yescarta (axicabtagene ciloleucel) is a CD19-directed genetically modified "
            "autologous T cell immunotherapy. FDA approved October 2017 for r/r large B-cell "
            "lymphoma after two or more lines of systemic therapy. Expanded April 2022 to "
            "include second-line treatment for LBCL. REMS program required. Boxed warning "
            "for CRS and neurological toxicities including cerebral edema."
        ),
        "regulatory_event": "label_update",
        "indication": "r/r large B-cell lymphoma, 2L LBCL",
    },
    {
        "product": "Tecartus (brexucabtagene autoleucel)",
        "set_id": "tecartus-spl-fallback",
        "title": "TECARTUS - brexucabtagene autoleucel suspension for intravenous infusion",
        "effective_time": "2024-04",
        "text_summary": (
            "Tecartus (brexucabtagene autoleucel) is a CD19-directed genetically modified "
            "autologous T cell immunotherapy. FDA approved July 2020 for r/r mantle cell "
            "lymphoma. October 2021 approval expanded to r/r B-cell ALL in adults. "
            "REMS required. Boxed warning for CRS and neurological toxicities."
        ),
        "regulatory_event": "label_update",
        "indication": "r/r mantle cell lymphoma, r/r B-ALL (adult)",
    },
    {
        "product": "Breyanzi (lisocabtagene maraleucel)",
        "set_id": "breyanzi-spl-fallback",
        "title": "BREYANZI - lisocabtagene maraleucel suspension for intravenous infusion",
        "effective_time": "2024-05",
        "text_summary": (
            "Breyanzi (lisocabtagene maraleucel) is a CD19-directed genetically modified "
            "autologous T cell immunotherapy with defined CD4+ and CD8+ CAR-T cell "
            "composition. FDA approved February 2021 for r/r large B-cell lymphoma after "
            "two or more lines. June 2024 expanded to second-line LBCL. REMS required. "
            "Boxed warning for CRS and neurological toxicities."
        ),
        "regulatory_event": "label_update",
        "indication": "r/r large B-cell lymphoma, 2L LBCL",
    },
    {
        "product": "Abecma (idecabtagene vicleucel)",
        "set_id": "abecma-spl-fallback",
        "title": "ABECMA - idecabtagene vicleucel suspension for intravenous infusion",
        "effective_time": "2024-04",
        "text_summary": (
            "Abecma (idecabtagene vicleucel) is a BCMA-directed genetically modified "
            "autologous T cell immunotherapy. FDA approved March 2021 for r/r multiple "
            "myeloma after four or more prior lines of therapy including an immunomodulatory "
            "agent, a proteasome inhibitor, and an anti-CD38 monoclonal antibody. "
            "REMS required. Boxed warning for CRS and neurological toxicities. "
            "Warning for secondary T-cell malignancies."
        ),
        "regulatory_event": "full_approval",
        "indication": "r/r multiple myeloma (4L+)",
    },
    {
        "product": "Carvykti (ciltacabtagene autoleucel)",
        "set_id": "carvykti-spl-fallback",
        "title": "CARVYKTI - ciltacabtagene autoleucel suspension for intravenous infusion",
        "effective_time": "2024-07",
        "text_summary": (
            "Carvykti (ciltacabtagene autoleucel) is a BCMA-directed genetically modified "
            "autologous T cell immunotherapy with two BCMA-targeting single-domain antibodies. "
            "FDA approved February 2022 for r/r multiple myeloma after four or more prior "
            "lines. April 2024 expanded to include patients after one or more prior lines. "
            "REMS required. Boxed warning for CRS, neurological toxicities, and HLH/MAS. "
            "Warning for secondary T-cell malignancies including CAR-positive cases."
        ),
        "regulatory_event": "label_update",
        "indication": "r/r multiple myeloma (1L+, 4L+)",
    },
]

# Maximum retry attempts for network requests
MAX_RETRIES = 3

# Delay between requests to be respectful to the DailyMed API
REQUEST_DELAY_SEC = 0.5


class DailyMedIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for DailyMed FDA drug label data for CAR-T products.

    Fetches Structured Product Labeling (SPL) data for each FDA-approved
    CAR-T product via the DailyMed REST API, parses the responses into
    RegulatoryRecord models, and stores embeddings in the cart_regulatory
    Milvus collection.

    Falls back to curated static seed data if the DailyMed API is unavailable.

    Usage:
        pipeline = DailyMedIngestPipeline(collection_manager, embedder)
        count = pipeline.run()
    """

    COLLECTION_NAME = "cart_regulatory"

    def __init__(
        self,
        collection_manager: CARTCollectionManager,
        embedder: Any,
    ):
        """Initialize the DailyMed ingest pipeline.

        Args:
            collection_manager: CARTCollectionManager for Milvus operations.
            embedder: Embedding model with encode() method.
        """
        super().__init__(collection_manager, embedder)

    def fetch(
        self,
        product_names: Optional[List[str]] = None,
    ) -> List[Dict[str, Any]]:
        """Fetch SPL data for CAR-T products from the DailyMed REST API.

        Queries the DailyMed SPL search endpoint for each CAR-T product by
        name. If the API is unreachable or returns no results, falls back to
        curated static seed data.

        API endpoint: https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json
        Query parameters:
            drug_name — product name to search for

        Args:
            product_names: List of CAR-T product names to query.
                Defaults to all six FDA-approved CAR-T products.

        Returns:
            List of SPL data dicts with keys: set_id, title, product,
            effective_time, plus any additional metadata.
        """
        names = product_names or CART_PRODUCT_NAMES
        all_spls: List[Dict[str, Any]] = []
        api_available = True

        for name in names:
            if not api_available:
                break

            params = {"drug_name": name}
            response = self._request_with_retry(DAILYMED_SPL_URL, params)

            if response is None:
                logger.warning(
                    f"DailyMed API unavailable for '{name}'; "
                    "will fall back to static seed data"
                )
                api_available = False
                break

            data = response.json()
            spl_data = data.get("data", [])

            if not spl_data:
                logger.info(f"No DailyMed SPL results for '{name}'")
                continue

            # Attach the queried product name for downstream parsing
            for spl in spl_data:
                spl["_queried_product"] = name

            all_spls.extend(spl_data)
            logger.info(
                f"DailyMed: fetched {len(spl_data)} SPL records for '{name}'"
            )

            # Respectful rate limiting
            time.sleep(REQUEST_DELAY_SEC)

        if not all_spls:
            logger.info(
                "No live DailyMed results obtained; using static fallback seed data"
            )
            return _FALLBACK_SEED_DATA

        logger.info(f"DailyMed fetch complete: {len(all_spls)} total SPL records")
        return all_spls

    def parse(self, raw_data: List[Dict[str, Any]]) -> List[RegulatoryRecord]:
        """Parse DailyMed SPL data into RegulatoryRecord models.

        Handles two input formats:
          1. Live API data: contains set_id, title, effective_time fields from
             the DailyMed API response
          2. Fallback seed data: pre-structured dicts matching RegulatoryRecord
             fields directly

        Args:
            raw_data: List of SPL data dicts from fetch().

        Returns:
            List of validated RegulatoryRecord model instances.
        """
        records: List[RegulatoryRecord] = []

        for idx, spl in enumerate(raw_data):
            try:
                record = self._parse_single_spl(spl, idx)
                if record is not None:
                    records.append(record)
            except Exception as exc:
                spl_id = spl.get("set_id", spl.get("setid", f"index-{idx}"))
                logger.warning(f"Failed to parse DailyMed SPL {spl_id}: {exc}")
                continue

        logger.info(
            f"Parsed {len(records)} RegulatoryRecord instances "
            f"from {len(raw_data)} SPL entries"
        )
        return records

    def _parse_single_spl(
        self, spl: Dict[str, Any], idx: int
    ) -> Optional[RegulatoryRecord]:
        """Parse a single DailyMed SPL entry into a RegulatoryRecord.

        Args:
            spl: Single SPL data dict (either live API or fallback format).
            idx: Index in the batch for fallback ID generation.

        Returns:
            RegulatoryRecord instance, or None if the entry lacks required data.
        """
        # Check if this is fallback seed data (already has text_summary)
        if "text_summary" in spl and "product" in spl:
            return self._parse_fallback_record(spl, idx)

        # --- Live DailyMed API response parsing ---
        set_id = spl.get("setid", spl.get("set_id", ""))
        if not set_id:
            return None

        record_id = f"DAILYMED-{set_id[:80]}"
        title = spl.get("title", "").strip()
        if not title:
            return None

        # Determine product name from query or title
        queried_product = spl.get("_queried_product", "").lower()
        product = self._resolve_product_name(queried_product, title)

        # Extract date from effective_time (format varies: YYYYMMDD or YYYY-MM-DD)
        effective_time = spl.get("published_date", spl.get("effective_time", ""))
        date_str = self._normalize_date(effective_time)

        # Get product metadata for enrichment
        meta = _PRODUCT_METADATA.get(queried_product, {})
        indication = meta.get("initial_indication", "")

        # Build text summary from title and available metadata
        text_parts = [f"FDA drug label (DailyMed SPL) for {product}."]
        text_parts.append(f"Title: {title}.")
        if indication:
            text_parts.append(f"Indication(s): {indication}.")
        if meta.get("target"):
            text_parts.append(f"Target antigen: {meta['target']}.")
        if meta.get("sponsor"):
            text_parts.append(f"Sponsor: {meta['sponsor']}.")

        text_summary = " ".join(text_parts)
        if len(text_summary) > 2900:
            text_summary = text_summary[:2897] + "..."

        # Determine regulatory event type based on context
        regulatory_event = RegulatoryEvent.LABEL_UPDATE

        return RegulatoryRecord(
            id=record_id[:100],
            text_summary=text_summary,
            product=product[:200],
            regulatory_event=regulatory_event,
            date=date_str[:20],
            agency="FDA",
            indication=indication[:200],
            decision="approved",
            conditions="REMS required",
            pivotal_trial="",
        )

    def _parse_fallback_record(
        self, spl: Dict[str, Any], idx: int
    ) -> Optional[RegulatoryRecord]:
        """Parse a fallback seed data dict into a RegulatoryRecord.

        Args:
            spl: Fallback seed data dict with pre-structured fields.
            idx: Index for fallback ID generation.

        Returns:
            RegulatoryRecord instance.
        """
        set_id = spl.get("set_id", f"fallback-{idx:03d}")
        record_id = f"DAILYMED-{set_id}"

        # Map regulatory_event string to enum
        event_str = spl.get("regulatory_event", "label_update")
        try:
            regulatory_event = RegulatoryEvent(event_str)
        except ValueError:
            regulatory_event = RegulatoryEvent.LABEL_UPDATE

        date_str = spl.get("effective_time", spl.get("date", ""))

        return RegulatoryRecord(
            id=record_id[:100],
            text_summary=spl.get("text_summary", "")[:3000],
            product=spl.get("product", "")[:200],
            regulatory_event=regulatory_event,
            date=date_str[:20],
            agency="FDA",
            indication=spl.get("indication", "")[:200],
            decision="approved",
            conditions="REMS required",
            pivotal_trial="",
        )

    @staticmethod
    def _resolve_product_name(queried_product: str, title: str) -> str:
        """Resolve the full product name from the query or title.

        Args:
            queried_product: The product name that was queried (lowercase).
            title: The SPL title string.

        Returns:
            Human-readable product name with generic name.
        """
        if queried_product in _PRODUCT_METADATA:
            meta = _PRODUCT_METADATA[queried_product]
            return f"{queried_product.capitalize()} ({meta['generic']})"

        # Try to extract from title
        for name, meta in _PRODUCT_METADATA.items():
            if name.lower() in title.lower():
                return f"{name.capitalize()} ({meta['generic']})"

        return title[:200] if title else "Unknown CAR-T product"

    @staticmethod
    def _normalize_date(date_str: str) -> str:
        """Normalize date string to YYYY-MM format.

        Handles formats: YYYYMMDD, YYYY-MM-DD, YYYY-MM, YYYY.

        Args:
            date_str: Raw date string from the API.

        Returns:
            Normalized date string in YYYY-MM or YYYY format.
        """
        if not date_str:
            return ""

        cleaned = date_str.strip().replace("/", "-")

        # YYYYMMDD format
        if len(cleaned) == 8 and cleaned.isdigit():
            return f"{cleaned[:4]}-{cleaned[4:6]}"

        # Already YYYY-MM-DD or YYYY-MM
        if "-" in cleaned:
            parts = cleaned.split("-")
            if len(parts) >= 2:
                return f"{parts[0]}-{parts[1]}"
            return parts[0]

        # YYYY only
        if len(cleaned) == 4 and cleaned.isdigit():
            return cleaned

        return cleaned[:20]

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
                logger.warning(
                    f"DailyMed HTTP error {status_code} on attempt "
                    f"{attempt}/{max_retries}: {exc}"
                )
                if attempt < max_retries:
                    time.sleep(2 ** attempt)
                    continue
                return None
            except requests.exceptions.RequestException as exc:
                logger.warning(
                    f"DailyMed request failed on attempt "
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
        **fetch_kwargs,
    ) -> int:
        """Execute the full DailyMed ingest pipeline: fetch -> parse -> embed_and_store.

        Args:
            collection_name: Target Milvus collection (defaults to 'cart_regulatory').
            batch_size: Batch size for embedding and insertion.
            **fetch_kwargs: Additional keyword arguments passed to fetch().

        Returns:
            Total number of records ingested.
        """
        target = collection_name or self.COLLECTION_NAME
        logger.info(f"Starting DailyMed ingest pipeline -> {target}")

        raw = self.fetch(**fetch_kwargs)
        logger.info(f"Fetched {len(raw)} raw SPL entries from DailyMed")

        records = self.parse(raw)
        logger.info(f"Parsed {len(records)} RegulatoryRecord instances")

        count = self.embed_and_store(records, target, batch_size)
        logger.info(f"DailyMed ingest complete: {count} records into {target}")
        return count
