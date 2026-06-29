"""CPIC guideline ingest pipeline for Pharmacogenomics Intelligence Agent.

Fetches clinical pharmacogenomics guidelines from the CPIC API
(https://api.cpicpgx.org/v1/), parses gene-drug pairs, phenotype
recommendations, and dosing guidelines into DrugGuideline models,
and stores embeddings in the pgx_drug_guidelines Milvus collection.

CPIC API docs: https://cpicpgx.org/api/

Author: Adam Jones
Date: March 2026
"""

from typing import Any, Dict, List

import requests
from loguru import logger

from src.collections import PGxCollectionManager
from src.models import (
    AlertLevel,
    ClinicalAction,
    CPICLevel,
    DrugGuideline,
    GuidelineBody,
)

from .base import BaseIngestPipeline, _truncate_utf8


# CPIC API endpoints
CPIC_BASE_URL = "https://api.cpicpgx.org/v1"
CPIC_PAIRS_URL = f"{CPIC_BASE_URL}/pair"
CPIC_GUIDELINES_URL = f"{CPIC_BASE_URL}/guideline"
CPIC_RECOMMENDATIONS_URL = f"{CPIC_BASE_URL}/recommendation"
CPIC_DRUG_URL = f"{CPIC_BASE_URL}/drug"

# Map CPIC level strings to our enum
_CPIC_LEVEL_MAP = {
    "A": CPICLevel.A,
    "A/B": CPICLevel.A_B,
    "B": CPICLevel.B,
    "C": CPICLevel.C,
    "D": CPICLevel.D,
}

# Keywords for classifying clinical action from recommendation text
_ACTION_KEYWORDS = {
    ClinicalAction.CONTRAINDICATED: ["contraindicated", "do not use", "must not"],
    ClinicalAction.AVOID: ["avoid", "not recommended", "consider alternative"],
    ClinicalAction.ALTERNATIVE: ["alternative", "switch to", "use alternative"],
    ClinicalAction.DOSE_ADJUST: [
        "reduce dose", "increase dose", "dose adjustment",
        "lower dose", "higher dose", "50%", "25%", "75%",
    ],
}

_ALERT_KEYWORDS = {
    AlertLevel.CRITICAL: ["contraindicated", "do not use", "avoid", "life-threatening"],
    AlertLevel.WARNING: ["caution", "reduce dose", "monitor", "increased risk"],
}


class CPICGuidelineParser(BaseIngestPipeline):
    """Ingest pipeline for CPIC pharmacogenomics guidelines.

    Fetches gene-drug pair recommendations from the CPIC REST API,
    parses phenotype-specific dosing guidance, and stores results
    in the pgx_drug_guidelines Milvus collection.

    Usage:
        parser = CPICGuidelineParser(collection_manager, embedder)
        count = parser.run()
    """

    COLLECTION_NAME = "pgx_drug_guidelines"

    def __init__(
        self,
        collection_manager: PGxCollectionManager,
        embedder: Any,
        request_timeout: int = 30,
    ):
        super().__init__(collection_manager, embedder)
        self.request_timeout = request_timeout

    def fetch(self, **kwargs) -> Dict[str, Any]:
        """Fetch CPIC guidelines, gene-drug pairs, and recommendations.

        Returns:
            Dict with keys: pairs, recommendations, guidelines.
        """
        raw: Dict[str, Any] = {"pairs": [], "recommendations": [], "guidelines": []}

        try:
            logger.info("Fetching CPIC gene-drug pairs...")
            resp = requests.get(CPIC_PAIRS_URL, timeout=self.request_timeout)
            resp.raise_for_status()
            raw["pairs"] = resp.json()
            logger.info(f"Fetched {len(raw['pairs'])} CPIC gene-drug pairs")
        except requests.RequestException as e:
            logger.error(f"Failed to fetch CPIC pairs: {e}")

        try:
            logger.info("Fetching CPIC recommendations...")
            resp = requests.get(CPIC_RECOMMENDATIONS_URL, timeout=self.request_timeout)
            resp.raise_for_status()
            raw["recommendations"] = resp.json()
            logger.info(f"Fetched {len(raw['recommendations'])} CPIC recommendations")
        except requests.RequestException as e:
            logger.error(f"Failed to fetch CPIC recommendations: {e}")

        try:
            logger.info("Fetching CPIC guideline metadata...")
            resp = requests.get(CPIC_GUIDELINES_URL, timeout=self.request_timeout)
            resp.raise_for_status()
            raw["guidelines"] = resp.json()
            logger.info(f"Fetched {len(raw['guidelines'])} CPIC guidelines")
        except requests.RequestException as e:
            logger.error(f"Failed to fetch CPIC guidelines: {e}")

        return raw

    def parse(self, raw_data: Dict[str, Any]) -> List[DrugGuideline]:
        """Parse CPIC API responses into DrugGuideline models.

        Joins recommendations with gene-drug pair metadata to produce
        phenotype-specific guideline records.
        """
        records: List[DrugGuideline] = []

        # Build guideline version lookup
        guideline_versions: Dict[str, str] = {}
        for gl in raw_data.get("guidelines", []):
            gl_id = str(gl.get("id", ""))
            version = gl.get("version", "")
            guideline_versions[gl_id] = version

        # Build pair lookup by (gene, drug) for CPIC level
        pair_lookup: Dict[str, Dict] = {}
        for pair in raw_data.get("pairs", []):
            gene = pair.get("genesymbol", "")
            drug = pair.get("drugname", "")
            if gene and drug:
                pair_lookup[f"{gene}_{drug}"] = pair

        for rec in raw_data.get("recommendations", []):
            try:
                gene = rec.get("genesymbol", "") or ""
                drug = rec.get("drugname", "") or rec.get("drug", {}).get("name", "")
                phenotype = rec.get("phenotypes", {}).get("genesymbol", "") if isinstance(rec.get("phenotypes"), dict) else ""
                if not phenotype:
                    phenotype = rec.get("lookupkey", {}).get(gene, "") if isinstance(rec.get("lookupkey"), dict) else ""

                recommendation_text = rec.get("drugrecommendation", "") or ""
                classification = rec.get("classification", "") or ""

                # Determine CPIC level from the pair table
                pair_key = f"{gene}_{drug}"
                pair_data = pair_lookup.get(pair_key, {})
                cpic_level_str = pair_data.get("cpiclevel", "A")
                cpic_level = _CPIC_LEVEL_MAP.get(cpic_level_str, CPICLevel.A)

                # Classify clinical action from recommendation text
                clinical_action = self._classify_action(recommendation_text)

                # Classify alert level
                alert_level = self._classify_alert(recommendation_text)

                # Build text chunk for embedding
                text_chunk = (
                    f"{gene} {drug} {phenotype}. "
                    f"{recommendation_text} "
                    f"Classification: {classification}"
                ).strip()

                rec_id = _truncate_utf8(
                    f"cpic_{gene}_{drug}_{phenotype}".replace(" ", "_")[:95], 95
                )

                record = DrugGuideline(
                    id=rec_id,
                    gene=_truncate_utf8(gene, 48),
                    drug=_truncate_utf8(drug, 195),
                    phenotype=_truncate_utf8(phenotype, 95),
                    guideline_body=GuidelineBody.CPIC,
                    cpic_level=cpic_level,
                    recommendation=_truncate_utf8(recommendation_text, 1990),
                    clinical_action=clinical_action,
                    alert_level=alert_level,
                    alternative_drugs=_truncate_utf8(
                        rec.get("alternatedrugavailable", ""), 490
                    ),
                    dose_adjustment=_truncate_utf8(classification, 490),
                    evidence_pmids="",
                    guideline_version=guideline_versions.get(
                        str(pair_data.get("guidelineid", "")), ""
                    ),
                    last_updated=pair_data.get("lastmodified", "")[:10] if pair_data.get("lastmodified") else "",
                    text_chunk=_truncate_utf8(text_chunk, 2990),
                )
                records.append(record)

            except Exception as e:
                logger.warning(f"Failed to parse CPIC recommendation: {e}")
                continue

        logger.info(f"Parsed {len(records)} CPIC guideline records")
        return records

    @staticmethod
    def _classify_action(text: str) -> ClinicalAction:
        """Classify clinical action from recommendation text."""
        text_lower = text.lower()
        for action, keywords in _ACTION_KEYWORDS.items():
            if any(kw in text_lower for kw in keywords):
                return action
        return ClinicalAction.STANDARD

    @staticmethod
    def _classify_alert(text: str) -> AlertLevel:
        """Classify alert level from recommendation text."""
        text_lower = text.lower()
        for level, keywords in _ALERT_KEYWORDS.items():
            if any(kw in text_lower for kw in keywords):
                return level
        return AlertLevel.INFO
