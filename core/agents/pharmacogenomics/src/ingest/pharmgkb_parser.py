"""PharmGKB clinical annotation ingest pipeline for PGx Intelligence Agent.

Parses PharmGKB clinical annotations, variant annotations, and drug labels
into DrugInteraction and ClinicalEvidence models.  PharmGKB is the
Pharmacogenomics Knowledgebase, curating gene-drug-disease relationships
with clinical annotation levels of evidence.

PharmGKB API: https://api.pharmgkb.org/v1/data/

Author: Adam Jones
Date: March 2026
"""

from typing import Any, Dict, List, Optional

import requests
from loguru import logger
from pydantic import BaseModel

from src.collections import PGxCollectionManager
from src.models import (
    ClinicalEvidence,
    DrugInteraction,
    EvidenceLevel,
    InteractionType,
)

from .base import BaseIngestPipeline, _truncate_utf8


# PharmGKB API endpoints
PHARMGKB_BASE_URL = "https://api.pharmgkb.org/v1/data"
PHARMGKB_CLINICAL_ANN_URL = f"{PHARMGKB_BASE_URL}/clinicalAnnotation"
PHARMGKB_VARIANT_ANN_URL = f"{PHARMGKB_BASE_URL}/variantAnnotation"
PHARMGKB_DRUG_LABEL_URL = f"{PHARMGKB_BASE_URL}/drugLabel"
PHARMGKB_GUIDELINE_URL = f"{PHARMGKB_BASE_URL}/guideline"

# Map PharmGKB evidence levels to our enum
_EVIDENCE_MAP = {
    "1A": EvidenceLevel.LEVEL_1A,
    "1B": EvidenceLevel.LEVEL_1B,
    "2A": EvidenceLevel.LEVEL_2A,
    "2B": EvidenceLevel.LEVEL_2B,
    "3": EvidenceLevel.LEVEL_3,
    "4": EvidenceLevel.LEVEL_4,
}

# Keywords for classifying interaction type
_INTERACTION_KEYWORDS = {
    InteractionType.TOXICITY: ["toxicity", "adverse", "side effect", "hypersensitivity", "sjs", "ten"],
    InteractionType.EFFICACY: ["efficacy", "response", "resistance", "non-response", "treatment outcome"],
    InteractionType.PD: ["pharmacodynamic", "receptor", "target", "sensitivity"],
}


class PharmGKBParser(BaseIngestPipeline):
    """Ingest pipeline for PharmGKB clinical and variant annotations.

    Fetches gene-drug annotations from PharmGKB, parses clinical
    annotations into DrugInteraction models (pgx_drug_interactions)
    and variant-level evidence into ClinicalEvidence models
    (pgx_clinical_evidence).

    Usage:
        parser = PharmGKBParser(collection_manager, embedder)
        count = parser.run()
    """

    COLLECTION_NAME = "pgx_drug_interactions"

    def __init__(
        self,
        collection_manager: PGxCollectionManager,
        embedder: Any,
        request_timeout: int = 30,
    ):
        super().__init__(collection_manager, embedder)
        self.request_timeout = request_timeout

    def fetch(self, **kwargs) -> Dict[str, Any]:
        """Fetch clinical annotations and variant annotations from PharmGKB.

        Returns:
            Dict with keys: clinical_annotations, variant_annotations.
        """
        max_results = kwargs.get("max_results", 500)
        raw: Dict[str, Any] = {"clinical_annotations": [], "variant_annotations": []}

        # Fetch clinical annotations
        try:
            logger.info("Fetching PharmGKB clinical annotations...")
            params = {"view": "max", "offset": 0, "max": max_results}
            resp = requests.get(
                PHARMGKB_CLINICAL_ANN_URL,
                params=params,
                timeout=self.request_timeout,
            )
            resp.raise_for_status()
            data = resp.json()
            raw["clinical_annotations"] = data.get("data", data) if isinstance(data, dict) else data
            logger.info(f"Fetched {len(raw['clinical_annotations'])} clinical annotations")
        except requests.RequestException as e:
            logger.error(f"Failed to fetch PharmGKB clinical annotations: {e}")

        # Fetch variant annotations
        try:
            logger.info("Fetching PharmGKB variant annotations...")
            params = {"view": "max", "offset": 0, "max": max_results}
            resp = requests.get(
                PHARMGKB_VARIANT_ANN_URL,
                params=params,
                timeout=self.request_timeout,
            )
            resp.raise_for_status()
            data = resp.json()
            raw["variant_annotations"] = data.get("data", data) if isinstance(data, dict) else data
            logger.info(f"Fetched {len(raw['variant_annotations'])} variant annotations")
        except requests.RequestException as e:
            logger.error(f"Failed to fetch PharmGKB variant annotations: {e}")

        return raw

    def parse(self, raw_data: Dict[str, Any]) -> List[BaseModel]:
        """Parse PharmGKB data into DrugInteraction and ClinicalEvidence models.

        Clinical annotations become DrugInteraction records; variant
        annotations become ClinicalEvidence records.
        """
        records: List[BaseModel] = []

        # Parse clinical annotations -> DrugInteraction
        for ann in raw_data.get("clinical_annotations", []):
            try:
                record = self._parse_clinical_annotation(ann)
                if record:
                    records.append(record)
            except Exception as e:
                logger.warning(f"Failed to parse clinical annotation: {e}")
                continue

        # Parse variant annotations -> ClinicalEvidence
        for ann in raw_data.get("variant_annotations", []):
            try:
                record = self._parse_variant_annotation(ann)
                if record:
                    records.append(record)
            except Exception as e:
                logger.warning(f"Failed to parse variant annotation: {e}")
                continue

        logger.info(f"Parsed {len(records)} total PharmGKB records")
        return records

    def _parse_clinical_annotation(self, ann: Dict[str, Any]) -> Optional[DrugInteraction]:
        """Parse a single PharmGKB clinical annotation into a DrugInteraction."""
        ann_id = str(ann.get("id", ""))
        ann.get("location", {})
        gene_name = ""
        variant_rsid = ""

        # Extract gene and variant from relatedGenes / relatedVariants
        related_genes = ann.get("relatedGenes", [])
        if related_genes and isinstance(related_genes, list):
            first_gene = related_genes[0]
            gene_name = first_gene.get("symbol", "") if isinstance(first_gene, dict) else str(first_gene)

        related_chemicals = ann.get("relatedChemicals", [])
        drug_name = ""
        if related_chemicals and isinstance(related_chemicals, list):
            first_drug = related_chemicals[0]
            drug_name = first_drug.get("name", "") if isinstance(first_drug, dict) else str(first_drug)

        if not gene_name or not drug_name:
            return None

        related_variants = ann.get("relatedVariants", [])
        if related_variants and isinstance(related_variants, list):
            first_var = related_variants[0]
            variant_rsid = first_var.get("name", "") if isinstance(first_var, dict) else str(first_var)

        # Evidence level
        level_str = ann.get("evidenceLevel", "") or ann.get("level", "")
        evidence_level = _EVIDENCE_MAP.get(level_str, EvidenceLevel.LEVEL_2A)

        # Annotation text
        ann_text = ann.get("text", "") or ann.get("summary", "") or ""
        phenotype_cat = ann.get("phenotypeCategory", "") or ""
        significance = ann.get("significance", "") or ""

        # Classify interaction type
        interaction_type = self._classify_interaction_type(ann_text, phenotype_cat)

        text_chunk = f"{gene_name} {drug_name} {variant_rsid}. {ann_text}".strip()

        return DrugInteraction(
            id=_truncate_utf8(f"pharmgkb_ca_{ann_id}", 95),
            drug=_truncate_utf8(drug_name, 195),
            gene=_truncate_utf8(gene_name, 48),
            variant_rsid=_truncate_utf8(variant_rsid, 48),
            interaction_type=interaction_type,
            effect_description=_truncate_utf8(ann_text, 1990),
            evidence_level=evidence_level,
            clinical_significance=_truncate_utf8(significance, 490),
            pharmgkb_id=_truncate_utf8(ann_id, 48),
            affected_phenotype=_truncate_utf8(phenotype_cat, 195),
            text_chunk=_truncate_utf8(text_chunk, 2990),
        )

    def _parse_variant_annotation(self, ann: Dict[str, Any]) -> Optional[ClinicalEvidence]:
        """Parse a single PharmGKB variant annotation into a ClinicalEvidence."""
        ann_id = str(ann.get("id", ""))
        variant = ann.get("variant", {}) or {}
        variant_name = variant.get("name", "") if isinstance(variant, dict) else str(variant)

        gene_name = ""
        related_genes = ann.get("relatedGenes", ann.get("genes", []))
        if related_genes and isinstance(related_genes, list):
            first_gene = related_genes[0]
            gene_name = first_gene.get("symbol", "") if isinstance(first_gene, dict) else str(first_gene)

        drug_name = ""
        related_drugs = ann.get("relatedChemicals", ann.get("drugs", []))
        if related_drugs and isinstance(related_drugs, list):
            first_drug = related_drugs[0]
            drug_name = first_drug.get("name", "") if isinstance(first_drug, dict) else str(first_drug)

        sentence = ann.get("sentence", "") or ann.get("text", "") or ""
        study_type = ann.get("studyType", "") or ""

        # PMIDs from literature references
        pmids = []
        for lit in ann.get("literatures", ann.get("literature", [])):
            if isinstance(lit, dict):
                pmid = lit.get("pmid", "") or lit.get("pubmedId", "")
                if pmid:
                    pmids.append(str(pmid))

        title = f"{gene_name} {variant_name} {drug_name} variant annotation".strip()
        text_chunk = f"{title}. {sentence}".strip()

        return ClinicalEvidence(
            id=_truncate_utf8(f"pharmgkb_va_{ann_id}", 95),
            title=_truncate_utf8(title, 490),
            text_chunk=_truncate_utf8(text_chunk, 2990),
            study_type=_truncate_utf8(study_type, 95),
            gene=_truncate_utf8(gene_name, 48),
            drug=_truncate_utf8(drug_name, 195),
            phenotype="",
            outcome_measure=_truncate_utf8(ann.get("significance", ""), 195),
            outcome_value="",
            sample_size=0,
            pmid=_truncate_utf8(pmids[0] if pmids else "", 18),
            year=0,
            source="PharmGKB",
        )

    @staticmethod
    def _classify_interaction_type(text: str, phenotype_cat: str) -> InteractionType:
        """Classify the interaction type from annotation text and category."""
        combined = f"{text} {phenotype_cat}".lower()
        for itype, keywords in _INTERACTION_KEYWORDS.items():
            if any(kw in combined for kw in keywords):
                return itype
        return InteractionType.PK

    def run(
        self,
        collection_name: Optional[str] = None,
        batch_size: int = 32,
        **fetch_kwargs,
    ) -> int:
        """Run ingestion, storing DrugInteractions and ClinicalEvidence separately."""
        target_interactions = "pgx_drug_interactions"
        target_evidence = "pgx_clinical_evidence"

        raw_data = self.fetch(**fetch_kwargs)
        all_records = self.parse(raw_data)

        interactions = [r for r in all_records if isinstance(r, DrugInteraction)]
        evidence = [r for r in all_records if isinstance(r, ClinicalEvidence)]

        total = 0
        if interactions:
            logger.info(f"Storing {len(interactions)} DrugInteraction records")
            total += self.embed_and_store(interactions, target_interactions, batch_size)
        if evidence:
            logger.info(f"Storing {len(evidence)} ClinicalEvidence records")
            total += self.embed_and_store(evidence, target_evidence, batch_size)

        logger.info(f"PharmGKB ingestion complete: {total} records stored")
        return total
