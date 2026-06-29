"""
Oncology Intelligence Agent - Cross-Modal Triggers
Connect oncology variants to imaging intelligence and drug discovery pipelines.

When actionable variants are identified, this module queries across modalities
(genomic evidence, imaging findings) to enrich the clinical context. Follows
the cross-modal pattern established by the Imaging Intelligence Agent.

Author: Adam Jones
Date: February 2026
License: Apache 2.0
"""

import hashlib
import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from src.knowledge import ACTIONABLE_TARGETS, classify_variant_actionability

logger = logging.getLogger(__name__)


@dataclass
class CrossModalResult:
    """Result of a cross-modal trigger evaluation.

    Attributes:
        trigger_reason: Why the cross-modal trigger fired.
        actionable_variants: List of actionable variant dicts that triggered it.
        genomic_context: Retrieved genomic evidence snippets.
        imaging_context: Retrieved imaging findings (if imaging collections exist).
        genomic_hit_count: Number of genomic evidence hits.
        imaging_hit_count: Number of imaging evidence hits.
        enrichment_summary: Human-readable summary combining all cross-modal data.
    """

    trigger_reason: str
    actionable_variants: List[Dict] = field(default_factory=list)
    genomic_context: List[Dict] = field(default_factory=list)
    imaging_context: List[Dict] = field(default_factory=list)
    genomic_hit_count: int = 0
    imaging_hit_count: int = 0
    enrichment_summary: str = ""


class OncoCrossModalTrigger:
    """Cross-modal trigger connecting oncology variants to imaging and drug discovery.

    When variants with evidence level A or B are found, queries:
      1. genomic_evidence collection for variant context and literature
      2. imaging_* collections for relevant imaging correlates (e.g.,
         EGFR mutation + lung nodule imaging characteristics)

    Gracefully handles missing collections (imaging agent may not be deployed).
    """

    GENOMIC_COLLECTION = "genomic_evidence"
    IMAGING_COLLECTION_PREFIX = "imaging_"

    # Similarity threshold for vector search
    DEFAULT_THRESHOLD = 0.40

    def __init__(self, collection_manager, embedder, settings: Optional[Dict] = None):
        """Initialize the cross-modal trigger.

        Args:
            collection_manager: Milvus collection manager for vector storage.
            embedder: Embedding model wrapper (e.g., BGE-small-en-v1.5).
            settings: Optional settings dict with thresholds and configuration.
        """
        self.collection_manager = collection_manager
        self.embedder = embedder
        self.settings = settings or {}
        self.threshold = self.settings.get("cross_modal_threshold", self.DEFAULT_THRESHOLD)
        self.genomic_top_k = self.settings.get("genomic_top_k", 5)
        self.imaging_top_k = self.settings.get("imaging_top_k", 5)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def evaluate(self, case_or_variants: Dict[str, Any]) -> Optional[CrossModalResult]:
        """Evaluate whether cross-modal triggers should fire for given variants.

        Checks if any variants have evidence level A or B in ACTIONABLE_TARGETS.
        If so, queries genomic and imaging collections for enrichment context.

        Args:
            case_or_variants: Dict containing either:
              - "variants": list of variant dicts with gene, variant, actionability
              - Or a CaseSnapshot-like dict with variants, cancer_type, etc.

        Returns:
            CrossModalResult if actionable triggers found, None otherwise.
        """
        variants = case_or_variants.get("variants", [])
        cancer_type = case_or_variants.get("cancer_type", "")

        # Identify high-evidence actionable variants (level A or B)
        actionable = []
        for v in variants:
            gene = v.get("gene", "").upper().strip()
            actionability = v.get("actionability", "VUS")

            # If actionability not pre-computed, compute it
            if actionability == "VUS" and gene in ACTIONABLE_TARGETS:
                variant_str = v.get("variant", "")
                actionability = self._classify_actionability(gene, variant_str)

            if actionability in ("A", "B"):
                actionable.append({
                    "gene": gene,
                    "variant": v.get("variant", ""),
                    "evidence_level": actionability,
                    "consequence": v.get("consequence", ""),
                })

        if not actionable:
            logger.debug("No actionable A/B variants found; cross-modal trigger not fired")
            return None

        logger.info(
            "Cross-modal trigger fired: %d actionable A/B variants",
            len(actionable),
        )

        # Build genomic queries
        genomic_queries = []
        imaging_queries = []
        for av in actionable:
            gene = av["gene"]
            variant = av["variant"]
            genomic_queries.append(f"{gene} {variant} targeted therapy evidence")
            genomic_queries.append(f"{gene} mutation clinical significance")

            # Build imaging queries with cancer-type context
            if cancer_type:
                imaging_queries.append(f"{gene} mutation {cancer_type} imaging findings")
            else:
                imaging_queries.append(f"{gene} mutation imaging characteristics")

        # Query genomic evidence
        genomic_hits = self._query_genomics(genomic_queries)
        logger.info("Genomic evidence: %d hits", len(genomic_hits))

        # Query imaging (graceful failure if collections don't exist)
        imaging_hits = self._query_imaging(imaging_queries)
        logger.info("Imaging evidence: %d hits", len(imaging_hits))

        # Build enrichment summary
        enrichment_summary = self._build_enrichment_summary(
            genomic_hits, imaging_hits, actionable
        )

        # Construct trigger reason
        genes = list({av["gene"] for av in actionable})
        trigger_reason = (
            f"High-evidence actionable variants detected in {', '.join(genes)}. "
            f"Cross-modal enrichment retrieved {len(genomic_hits)} genomic and "
            f"{len(imaging_hits)} imaging evidence items."
        )

        return CrossModalResult(
            trigger_reason=trigger_reason,
            actionable_variants=actionable,
            genomic_context=genomic_hits,
            imaging_context=imaging_hits,
            genomic_hit_count=len(genomic_hits),
            imaging_hit_count=len(imaging_hits),
            enrichment_summary=enrichment_summary,
        )

    # ------------------------------------------------------------------
    # Genomic evidence retrieval
    # ------------------------------------------------------------------

    def _query_genomics(self, queries: List[str]) -> List[Dict]:
        """Embed queries and search the genomic_evidence collection.

        Args:
            queries: List of natural-language query strings.

        Returns:
            Deduplicated list of genomic evidence dicts with text, source, score.
        """
        hits = []
        seen_texts = set()

        for query in queries:
            try:
                embedding = self.embedder.embed(query)
                results = self.collection_manager.search(
                    collection_name=self.GENOMIC_COLLECTION,
                    query_vector=embedding,
                    top_k=self.genomic_top_k,
                    output_fields=["text", "source", "gene", "title"],
                )
                for r in results:
                    score = r.get("score", 0.0)
                    if score < self.threshold:
                        continue
                    text = r.get("text", "")
                    # Deduplicate by text content hash
                    text_key = hashlib.md5(text.strip().lower().encode()).hexdigest()
                    if text_key in seen_texts:
                        continue
                    seen_texts.add(text_key)
                    hits.append({
                        "text": text,
                        "source": r.get("source", r.get("title", "")),
                        "gene": r.get("gene", ""),
                        "score": round(score, 4),
                        "query": query,
                    })
            except Exception as exc:
                logger.warning("Genomic query failed for '%s': %s", query[:50], exc)

        return hits

    # ------------------------------------------------------------------
    # Imaging evidence retrieval
    # ------------------------------------------------------------------

    def _query_imaging(self, queries: List[str]) -> List[Dict]:
        """Search imaging collections for relevant findings.

        Attempts to discover imaging_* collections and search them.
        Gracefully returns empty list if no imaging collections exist
        (imaging agent may not be deployed).

        Args:
            queries: List of natural-language imaging query strings.

        Returns:
            List of imaging evidence dicts, or empty list if unavailable.
        """
        hits = []
        seen_texts = set()

        # Discover imaging collections
        imaging_collections = self._discover_imaging_collections()
        if not imaging_collections:
            logger.debug("No imaging collections found; skipping imaging cross-modal query")
            return []

        for collection_name in imaging_collections:
            for query in queries:
                try:
                    embedding = self.embedder.embed(query)
                    results = self.collection_manager.search(
                        collection_name=collection_name,
                        query_vector=embedding,
                        top_k=self.imaging_top_k,
                        output_fields=["text", "source", "modality", "finding", "title"],
                    )
                    for r in results:
                        score = r.get("score", 0.0)
                        if score < self.threshold:
                            continue
                        text = r.get("text", "")
                        text_key = text[:100].strip().lower()
                        if text_key in seen_texts:
                            continue
                        seen_texts.add(text_key)
                        hits.append({
                            "text": text,
                            "source": r.get("source", r.get("title", "")),
                            "modality": r.get("modality", ""),
                            "finding": r.get("finding", ""),
                            "score": round(score, 4),
                            "collection": collection_name,
                            "query": query,
                        })
                except Exception as exc:
                    logger.debug(
                        "Imaging query failed for collection '%s': %s",
                        collection_name, exc,
                    )

        return hits

    def _discover_imaging_collections(self) -> List[str]:
        """Discover available imaging_* collections in Milvus.

        Returns:
            List of collection names starting with 'imaging_'.
        """
        try:
            all_collections = self.collection_manager.list_collections()
            imaging = [
                c for c in all_collections
                if c.startswith(self.IMAGING_COLLECTION_PREFIX)
            ]
            return imaging
        except Exception as exc:
            logger.debug("Failed to list collections for imaging discovery: %s", exc)
            return []

    # ------------------------------------------------------------------
    # Enrichment summary
    # ------------------------------------------------------------------

    def _build_enrichment_summary(
        self,
        genomic_hits: List[Dict],
        imaging_hits: List[Dict],
        variants: List[Dict],
    ) -> str:
        """Build a human-readable enrichment summary from cross-modal results.

        Args:
            genomic_hits: List of genomic evidence dicts.
            imaging_hits: List of imaging evidence dicts.
            variants: List of actionable variant dicts.

        Returns:
            Multi-line enrichment summary string.
        """
        lines = []

        # Header
        genes = list({v["gene"] for v in variants})
        lines.append(f"Cross-Modal Enrichment Summary for: {', '.join(genes)}")
        lines.append("=" * 60)

        # Variants section
        lines.append("")
        lines.append(f"Actionable Variants ({len(variants)}):")
        for v in variants:
            lines.append(
                f"  - {v['gene']} {v['variant']} "
                f"[Evidence Level {v.get('evidence_level', '?')}]"
            )

        # Genomic evidence section
        lines.append("")
        lines.append(f"Genomic Evidence ({len(genomic_hits)} hits):")
        if genomic_hits:
            for i, hit in enumerate(genomic_hits[:10], 1):
                source = hit.get("source", "Unknown")
                score = hit.get("score", 0.0)
                text_preview = hit.get("text", "")[:120].replace("\n", " ")
                lines.append(f"  {i}. [{score:.3f}] {source}")
                lines.append(f"     {text_preview}...")
        else:
            lines.append("  No genomic evidence retrieved.")

        # Imaging evidence section
        lines.append("")
        lines.append(f"Imaging Evidence ({len(imaging_hits)} hits):")
        if imaging_hits:
            for i, hit in enumerate(imaging_hits[:10], 1):
                source = hit.get("source", "Unknown")
                modality = hit.get("modality", "")
                score = hit.get("score", 0.0)
                text_preview = hit.get("text", "")[:120].replace("\n", " ")
                modality_str = f" ({modality})" if modality else ""
                lines.append(f"  {i}. [{score:.3f}] {source}{modality_str}")
                lines.append(f"     {text_preview}...")
        else:
            lines.append("  No imaging collections available or no relevant findings.")

        lines.append("")
        lines.append("=" * 60)

        return "\n".join(lines)

    # ------------------------------------------------------------------
    # Helper
    # ------------------------------------------------------------------

    def _classify_actionability(self, gene: str, variant: str) -> str:
        """Classify variant actionability using shared knowledge function.

        Args:
            gene: Gene symbol (uppercase).
            variant: Variant string.

        Returns:
            Evidence level ("A", "B", "C", "D") or "VUS".
        """
        return classify_variant_actionability(gene, variant)
