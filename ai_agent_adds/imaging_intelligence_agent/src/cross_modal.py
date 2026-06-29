"""Cross-modal trigger: Imaging -> Genomics pipeline integration.

When imaging findings meet severity thresholds (e.g., Lung-RADS 4A+),
automatically queries the genomic_evidence collection for relevant
cancer genomics context, enriching imaging findings with molecular
insights.

This bridges the Imaging Intelligence Agent with the 3.5M-vector
genomic_evidence collection populated by the rag-chat-pipeline,
enabling a single DGX Spark to connect radiology AI findings to
precision medicine genomics -- the core value proposition of the
HCLS AI Factory.

Trigger rules:
    - Lung-RADS 4A, 4B, 4X  -> lung cancer driver mutation queries
    - CXR consolidation (critical/urgent) -> infection genomics queries
    - Brain lesion high activity -> neurological genomics queries

Author: Adam Jones
Date: February 2026
"""

import time
from typing import Any, Dict, List, Optional

from loguru import logger

from src.models import CrossModalResult, FindingSeverity, WorkflowResult


# =====================================================================
# Genomic query templates per trigger type
# =====================================================================

LUNG_CANCER_QUERIES = [
    "lung cancer driver mutations EGFR ALK ROS1 KRAS",
    "non-small cell lung cancer NSCLC targeted therapy genomics",
    "lung adenocarcinoma molecular subtypes precision medicine",
]

INFECTION_GENOMICS_QUERIES = [
    "respiratory infection pathogen genomics antimicrobial resistance",
    "pneumonia host genomic susceptibility immune response",
]

NEURO_GENOMICS_QUERIES = [
    "multiple sclerosis genomic risk factors HLA-DRB1 susceptibility",
    "neuroinflammation demyelination genomic markers treatment response",
]


# =====================================================================
# CrossModalTrigger
# =====================================================================


class CrossModalTrigger:
    """Evaluates imaging workflow results and triggers cross-modal queries.

    When a workflow result exceeds a severity threshold, the trigger
    encodes genomic-domain queries via the shared BGE-small-en-v1.5
    embedder, searches the genomic_evidence Milvus collection, and
    returns a CrossModalResult with the enrichment context.

    Args:
        collection_manager: ImagingCollectionManager with access to
            the genomic_evidence collection (read-only).
        embedder: SentenceTransformer (BGE-small-en-v1.5) for encoding
            genomic query strings to 384-dim vectors.
        enabled: Master switch; when False, evaluate() always returns
            None.  Controlled by settings.CROSS_MODAL_ENABLED.
    """

    GENOMIC_COLLECTION = "genomic_evidence"
    TOP_K_PER_QUERY = 5
    SCORE_THRESHOLD = 0.40

    def __init__(self, collection_manager, embedder, enabled: bool = True):
        self.collection_manager = collection_manager
        self.embedder = embedder
        self.enabled = enabled

    # -----------------------------------------------------------------
    # Public API
    # -----------------------------------------------------------------

    def evaluate(self, workflow_result: WorkflowResult) -> Optional[CrossModalResult]:
        """Check if a workflow result should trigger cross-modal enrichment.

        Triggers:
            - Lung-RADS 4A, 4B, 4X -> query lung cancer genomics
            - CXR consolidation + critical/urgent -> query infection genomics
            - Brain lesion highly_active -> query neurological genomics

        Args:
            workflow_result: Completed WorkflowResult from any imaging
                workflow.

        Returns:
            CrossModalResult with genomic context if triggered, or None
            if the result does not meet any threshold.
        """
        if not self.enabled:
            return None

        # Check lung nodule triggers
        if workflow_result.workflow_name == "ct_chest_lung_nodule":
            return self._evaluate_lung_nodule(workflow_result)

        # Check CXR triggers
        if workflow_result.workflow_name == "cxr_rapid_findings":
            return self._evaluate_cxr(workflow_result)

        # Check brain lesion triggers
        if workflow_result.workflow_name == "mri_brain_ms_lesion":
            return self._evaluate_brain_lesion(workflow_result)

        return None

    # -----------------------------------------------------------------
    # Per-workflow evaluators
    # -----------------------------------------------------------------

    def _evaluate_lung_nodule(
        self, result: WorkflowResult
    ) -> Optional[CrossModalResult]:
        """Trigger genomic query for Lung-RADS 4A+ findings.

        Parses the classification string (e.g. "Lung-RADS 4B") and
        triggers if the category is 4A, 4B, or 4X.
        """
        classification = result.classification or ""

        # Check if Lung-RADS is 4A, 4B, or 4X
        high_risk_categories = ["4A", "4B", "4X"]
        high_risk = any(cat in classification for cat in high_risk_categories)
        if not high_risk:
            return None

        # Determine which specific category matched for the trigger reason
        matched = next(
            (cat for cat in high_risk_categories if cat in classification),
            "4A+",
        )
        trigger_reason = f"Lung-RADS {matched} — high-risk lung nodule"

        logger.info(
            f"Cross-modal trigger fired: {trigger_reason} "
            f"(classification={classification})"
        )

        return self._query_genomics(
            queries=LUNG_CANCER_QUERIES,
            trigger_reason=trigger_reason,
        )

    def _evaluate_cxr(self, result: WorkflowResult) -> Optional[CrossModalResult]:
        """Trigger genomic query for severe CXR findings.

        Only triggers for critical or urgent severity, indicating an
        acute finding like tension pneumothorax or consolidation.
        """
        if result.severity not in (FindingSeverity.CRITICAL, FindingSeverity.URGENT):
            return None

        # Check for consolidation specifically (suggests infection)
        has_consolidation = any(
            f.get("category") == "consolidation" or f.get("class_name") == "consolidation"
            for f in result.findings
        )
        if not has_consolidation:
            return None

        trigger_reason = (
            f"CXR {result.severity.value} — consolidation detected, "
            f"querying infection genomics"
        )

        logger.info(f"Cross-modal trigger fired: {trigger_reason}")

        return self._query_genomics(
            queries=INFECTION_GENOMICS_QUERIES,
            trigger_reason=trigger_reason,
        )

    def _evaluate_brain_lesion(
        self, result: WorkflowResult
    ) -> Optional[CrossModalResult]:
        """Trigger genomic query for highly active MS lesions.

        Only triggers when the classification indicates highly active
        disease (many new/enlarging lesions).
        """
        classification = result.classification or ""
        if "highly_active" not in classification:
            return None

        trigger_reason = (
            "MS highly active — querying neurological genomics for "
            "susceptibility and treatment response markers"
        )

        logger.info(f"Cross-modal trigger fired: {trigger_reason}")

        return self._query_genomics(
            queries=NEURO_GENOMICS_QUERIES,
            trigger_reason=trigger_reason,
        )

    # -----------------------------------------------------------------
    # Genomic collection search
    # -----------------------------------------------------------------

    def _query_genomics(
        self, queries: List[str], trigger_reason: str
    ) -> CrossModalResult:
        """Search the genomic_evidence collection with multiple queries.

        Encodes each query string with BGE-small-en-v1.5 and performs a
        vector similarity search against the 3.5M-vector genomic_evidence
        collection.  De-duplicates hits by ID and builds a summary.

        Args:
            queries: List of genomic-domain query strings.
            trigger_reason: Human-readable reason the trigger fired.

        Returns:
            CrossModalResult with genomic context, hit count, and summary.
        """
        start = time.time()
        all_hits: List[Dict[str, Any]] = []
        seen_ids: set = set()

        for query in queries:
            try:
                embedding = self.embedder.encode(
                    query, normalize_embeddings=True
                )
                # Convert numpy array to list if needed
                embedding_list = (
                    embedding.tolist()
                    if hasattr(embedding, "tolist")
                    else list(embedding)
                )

                hits = self.collection_manager.search(
                    collection_name=self.GENOMIC_COLLECTION,
                    query_embedding=embedding_list,
                    top_k=self.TOP_K_PER_QUERY,
                    score_threshold=self.SCORE_THRESHOLD,
                )

                for hit in hits:
                    hit_id = hit.get("id", "")
                    if hit_id not in seen_ids:
                        seen_ids.add(hit_id)
                        all_hits.append(hit)

            except Exception as e:
                logger.warning(
                    f"Cross-modal genomic query failed for '{query}': {e}"
                )

        elapsed_ms = (time.time() - start) * 1000

        # Extract text context from hits
        genomic_context = []
        for hit in all_hits:
            text = (
                hit.get("text_chunk")
                or hit.get("text")
                or hit.get("text_summary")
                or ""
            )
            if text:
                score = hit.get("score", 0.0)
                genomic_context.append(f"[score={score:.2f}] {text}")

        # Build enrichment summary
        if all_hits:
            enrichment_summary = (
                f"Cross-modal genomics enrichment triggered by {trigger_reason}. "
                f"Retrieved {len(all_hits)} unique genomic evidence records "
                f"from {len(queries)} queries in {elapsed_ms:.0f}ms. "
                f"Top score: {max(h.get('score', 0) for h in all_hits):.2f}."
            )
        else:
            enrichment_summary = (
                f"Cross-modal genomics enrichment triggered by {trigger_reason}, "
                f"but no genomic evidence met the similarity threshold "
                f"({self.SCORE_THRESHOLD})."
            )

        logger.info(
            f"Cross-modal query complete: {len(all_hits)} hits from "
            f"{len(queries)} queries in {elapsed_ms:.0f}ms"
        )

        return CrossModalResult(
            trigger_reason=trigger_reason,
            genomic_context=genomic_context,
            genomic_hit_count=len(all_hits),
            query_count=len(queries),
            enrichment_summary=enrichment_summary,
        )
