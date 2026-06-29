"""Multi-collection RAG engine for Clinical Trial Intelligence Agent.

Searches across all 14 trial-specific Milvus collections simultaneously
using parallel ThreadPoolExecutor, synthesises findings with clinical trial
knowledge augmentation, and generates grounded LLM responses with regulatory
and evidence citations.

Extends the pattern from: rag-chat-pipeline/src/rag_engine.py

Features:
- Parallel search via ThreadPoolExecutor (13 trial + 1 shared genomic collection)
- Settings-driven weights and parameters from config/settings.py
- Workflow-based dynamic weight boosting per TrialWorkflowType
- Milvus field-based filtering (phase, status, indication, sponsor)
- Citation relevance scoring (high/medium/low) with NCT/PMID link formatting
- Cross-collection entity linking for comprehensive trial queries
- Regulatory guidance retrieval with document identifiers
- Conversation memory for multi-turn trial consultations
- Patient context injection for patient-trial matching
- Confidence scoring based on evidence quality and collection diversity

Author: Adam Jones
Date: March 2026
"""

import json
import logging
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from config.settings import settings

from .agent import (
    TRIAL_SYSTEM_PROMPT,
    WORKFLOW_COLLECTION_BOOST,
    TrialWorkflowType,
    TrialResponse,
)

logger = logging.getLogger(__name__)

# =====================================================================
# CONVERSATION PERSISTENCE HELPERS
# =====================================================================

CONVERSATION_DIR = Path(__file__).parent.parent / "data" / "cache" / "conversations"
_CONVERSATION_TTL = timedelta(hours=24)


def _save_conversation(session_id: str, history: list):
    """Persist conversation to disk as JSON."""
    try:
        CONVERSATION_DIR.mkdir(parents=True, exist_ok=True)
        path = CONVERSATION_DIR / f"{session_id}.json"
        data = {
            "session_id": session_id,
            "updated": datetime.now(timezone.utc).isoformat(),
            "messages": history,
        }
        path.write_text(json.dumps(data, indent=2))
    except Exception as exc:
        logger.warning("Failed to persist conversation %s: %s", session_id, exc)


def _load_conversation(session_id: str) -> list:
    """Load conversation from disk, respecting 24-hour TTL."""
    try:
        path = CONVERSATION_DIR / f"{session_id}.json"
        if path.exists():
            data = json.loads(path.read_text())
            updated = datetime.fromisoformat(data["updated"])
            if datetime.now(timezone.utc) - updated < _CONVERSATION_TTL:
                return data.get("messages", [])
            else:
                path.unlink(missing_ok=True)  # Expired
    except Exception as exc:
        logger.warning("Failed to load conversation %s: %s", session_id, exc)
    return []


def _cleanup_expired_conversations():
    """Remove conversation files older than 24 hours."""
    try:
        if not CONVERSATION_DIR.exists():
            return
        cutoff = datetime.now(timezone.utc) - _CONVERSATION_TTL
        for path in CONVERSATION_DIR.glob("*.json"):
            try:
                data = json.loads(path.read_text())
                updated = datetime.fromisoformat(data["updated"])
                if updated < cutoff:
                    path.unlink()
            except Exception:
                pass
    except Exception as exc:
        logger.warning("Conversation cleanup error: %s", exc)


# Allowed characters for Milvus filter expressions to prevent injection
_SAFE_FILTER_RE = re.compile(r"^[A-Za-z0-9 _.\-/\*:(),]+$")


# =====================================================================
# SEARCH RESULT DATACLASS
# =====================================================================

@dataclass
class TrialSearchResult:
    """A single search result from a Milvus collection.

    Attributes:
        collection: Source collection name (e.g. 'trial_protocols').
        record_id: Milvus record primary key.
        score: Weighted relevance score (0.0 - 1.0+).
        text: Primary text content from the record.
        metadata: Full record metadata dict from Milvus.
        relevance: Citation relevance tier ('high', 'medium', 'low').
    """
    collection: str = ""
    record_id: str = ""
    score: float = 0.0
    text: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)
    relevance: str = "low"


# =====================================================================
# COLLECTION CONFIGURATION (reads weights from settings)
# =====================================================================

COLLECTION_CONFIG: Dict[str, Dict[str, Any]] = {
    "trial_protocols": {
        "weight": settings.WEIGHT_PROTOCOLS,
        "label": "Protocol",
        "text_field": "design_summary",
        "title_field": "trial_title",
        "filterable_fields": ["phase", "status", "indication", "sponsor"],
    },
    "trial_eligibility": {
        "weight": settings.WEIGHT_ELIGIBILITY,
        "label": "Eligibility",
        "text_field": "criteria_text",
        "title_field": "trial_title",
        "filterable_fields": ["criteria_type", "indication"],
    },
    "trial_endpoints": {
        "weight": settings.WEIGHT_ENDPOINTS,
        "label": "Endpoint",
        "text_field": "endpoint_description",
        "title_field": "endpoint_name",
        "filterable_fields": ["endpoint_type", "therapeutic_area"],
    },
    "trial_sites": {
        "weight": settings.WEIGHT_SITES,
        "label": "Site",
        "text_field": "site_capabilities",
        "title_field": "site_name",
        "filterable_fields": ["country", "therapeutic_area", "performance_tier"],
    },
    "trial_investigators": {
        "weight": settings.WEIGHT_INVESTIGATORS,
        "label": "Investigator",
        "text_field": "expertise_summary",
        "title_field": "investigator_name",
        "filterable_fields": ["therapeutic_area", "country"],
    },
    "trial_results": {
        "weight": settings.WEIGHT_RESULTS,
        "label": "TrialResult",
        "text_field": "results_summary",
        "title_field": "trial_name",
        "filterable_fields": ["phase", "outcome", "indication"],
    },
    "trial_regulatory": {
        "weight": settings.WEIGHT_REGULATORY,
        "label": "Regulatory",
        "text_field": "guidance_text",
        "title_field": "document_title",
        "filterable_fields": ["agency", "document_type", "therapeutic_area"],
    },
    "trial_literature": {
        "weight": settings.WEIGHT_LITERATURE,
        "label": "Literature",
        "text_field": "abstract",
        "title_field": "title",
        "filterable_fields": ["study_type", "therapeutic_area"],
    },
    "trial_biomarkers": {
        "weight": settings.WEIGHT_BIOMARKERS,
        "label": "Biomarker",
        "text_field": "biomarker_context",
        "title_field": "biomarker_name",
        "filterable_fields": ["biomarker_type", "indication"],
    },
    "trial_safety": {
        "weight": settings.WEIGHT_SAFETY,
        "label": "Safety",
        "text_field": "safety_narrative",
        "title_field": "event_term",
        "filterable_fields": ["severity", "causality", "system_organ_class"],
    },
    "trial_rwe": {
        "weight": settings.WEIGHT_RWE,
        "label": "RWE",
        "text_field": "evidence_summary",
        "title_field": "study_title",
        "filterable_fields": ["data_source", "indication"],
    },
    "trial_adaptive": {
        "weight": settings.WEIGHT_ADAPTIVE,
        "label": "Adaptive",
        "text_field": "design_description",
        "title_field": "design_type",
        "filterable_fields": ["adaptation_type", "therapeutic_area"],
    },
    "trial_guidelines": {
        "weight": settings.WEIGHT_GUIDELINES,
        "label": "Guideline",
        "text_field": "recommendation",
        "title_field": "guideline_title",
        "filterable_fields": ["agency", "document_id", "therapeutic_area"],
    },
    "genomic_evidence": {
        "weight": settings.WEIGHT_GENOMIC,
        "label": "Genomic",
        "text_field": "text_chunk",
        "title_field": "gene",
        "filterable_fields": [],
    },
}

ALL_COLLECTION_NAMES = list(COLLECTION_CONFIG.keys())


def get_all_collection_names() -> List[str]:
    """Return all collection names."""
    return list(COLLECTION_CONFIG.keys())


# =====================================================================
# TRIAL RAG ENGINE
# =====================================================================

class TrialRAGEngine:
    """Multi-collection RAG engine for clinical trial intelligence.

    Searches across all 14 trial-specific Milvus collections plus the
    shared genomic_evidence collection. Supports workflow-specific weight
    boosting, parallel search, query expansion, patient context injection,
    and multi-turn conversation memory.

    Features:
    - Parallel search via ThreadPoolExecutor (14 collections)
    - Settings-driven weights and parameters
    - Workflow-based dynamic weight boosting (13 trial workflows)
    - Milvus field-based filtering (phase, status, indication, sponsor)
    - Citation relevance scoring (high/medium/low)
    - Cross-collection entity linking
    - Regulatory guidance retrieval with document IDs
    - Conversation memory context injection
    - Patient context for personalised trial matching
    - Confidence scoring based on evidence diversity

    Usage:
        engine = TrialRAGEngine(milvus_client, embedding_model, llm_client)
        response = engine.query("What endpoints work for NSCLC Phase 3?")
        results = engine.search("adaptive design Alzheimer's")
    """

    def __init__(
        self,
        milvus_client=None,
        embedding_model=None,
        llm_client=None,
        session_id: str = "default",
    ):
        """Initialize the TrialRAGEngine.

        Args:
            milvus_client: Connected Milvus client with access to all
                trial collections. If None, search operations will
                raise RuntimeError.
            embedding_model: Embedding model (BGE-small-en-v1.5) for query
                vectorisation. If None, embedding operations will raise.
            llm_client: LLM client (Anthropic Claude) for response synthesis.
                If None, search-only mode is available.
            session_id: Conversation session identifier for persistence
                (default: "default").
        """
        self.milvus = milvus_client
        self.embedder = embedding_model
        self.llm = llm_client
        self.session_id = session_id
        self._max_conversation_context = settings.MAX_CONVERSATION_CONTEXT

        # Load persisted conversation history (falls back to empty list)
        self._conversation_history: List[Dict[str, str]] = _load_conversation(session_id)

        # Cleanup expired conversations on startup
        _cleanup_expired_conversations()

    # ══════════════════════════════════════════════════════════════════
    # PROPERTIES
    # ══════════════════════════════════════════════════════════════════

    @property
    def conversation_history(self) -> List[Dict[str, str]]:
        """Return current conversation history."""
        return self._conversation_history

    # ══════════════════════════════════════════════════════════════════
    # PUBLIC API
    # ══════════════════════════════════════════════════════════════════

    def query(
        self,
        question: str,
        workflow: Optional[TrialWorkflowType] = None,
        top_k: int = 5,
        patient_context: Optional[dict] = None,
    ) -> TrialResponse:
        """Main query method: expand -> search -> synthesise.

        Performs the full RAG pipeline: parallel multi-collection search
        with workflow-specific weighting, result reranking, LLM synthesis
        with patient context, and confidence scoring.

        Args:
            question: Natural language clinical trial question.
            workflow: Optional TrialWorkflowType to apply domain-specific
                collection weight boosting. If None, weights are auto-detected
                or base defaults are used.
            top_k: Maximum results to return per collection.
            patient_context: Optional dict with patient-specific data
                (age, sex, diagnosis, biomarkers, prior_therapies, mutations)
                for personalised trial matching.

        Returns:
            TrialResponse with synthesised answer, search results, citations,
            confidence score, and metadata.
        """
        start = time.time()

        # Step 1: Determine collections and weights
        weights = self._get_boosted_weights(workflow)
        collections = list(weights.keys())

        # Step 2: Search across collections
        results = self.search(
            question=question,
            collections=collections,
            top_k=top_k,
        )

        # Step 3: Apply workflow-specific reranking
        results = self._rerank_results(results, question)

        # Step 4: Score citations
        results = self._score_citations(results)

        # Step 5: Score confidence
        confidence = self._score_confidence(results)

        # Step 6: Synthesise LLM response (if LLM available)
        if self.llm:
            response = self._synthesize_response(
                question=question,
                results=results,
                workflow=workflow,
                patient_context=patient_context,
            )
        else:
            response = TrialResponse(
                question=question,
                answer="[LLM not configured -- search-only mode. "
                       "See results below for retrieved evidence.]",
                results=results,
                workflow=workflow,
                confidence=confidence,
            )

        # Step 7: Extract citations
        response.citations = self._extract_citations(results)
        response.confidence = confidence
        response.search_time_ms = (time.time() - start) * 1000
        response.collections_searched = len(collections)
        response.patient_context_used = patient_context is not None

        # Step 8: Update conversation history
        self.add_conversation_context("user", question)
        if response.answer:
            self.add_conversation_context("assistant", response.answer[:500])

        return response

    def search(
        self,
        question: str,
        collections: Optional[List[str]] = None,
        top_k: int = 5,
    ) -> List[TrialSearchResult]:
        """Search across multiple collections with weighted scoring.

        Embeds the query, runs parallel Milvus searches across all specified
        collections, applies collection weights, and returns merged ranked
        results.

        Args:
            question: Natural language search query.
            collections: Optional list of collection names to search.
                If None, all 14 collections are searched.
            top_k: Maximum results per collection.

        Returns:
            List of TrialSearchResult sorted by weighted score descending.
        """
        if not self.milvus:
            raise RuntimeError(
                "Milvus client not configured. Cannot perform search."
            )

        # Embed query
        query_vector = self._embed_query(question)

        # Determine collections
        if not collections:
            collections = get_all_collection_names()

        # Get weights (base defaults for search-only calls)
        weights = {
            name: COLLECTION_CONFIG.get(name, {}).get("weight", 0.05)
            for name in collections
        }

        # Parallel search with weighting
        results = self._parallel_search(query_vector, collections, weights, top_k)

        return results

    # ══════════════════════════════════════════════════════════════════
    # EMBEDDING
    # ══════════════════════════════════════════════════════════════════

    def _embed_query(self, text: str) -> List[float]:
        """Generate embedding vector for query text.

        Uses the BGE instruction prefix for optimal retrieval performance
        with BGE-small-en-v1.5.

        Args:
            text: Query text to embed.

        Returns:
            384-dimensional float vector.

        Raises:
            RuntimeError: If embedding model is not configured.
        """
        if not self.embedder:
            raise RuntimeError(
                "Embedding model not configured. Cannot generate query vector."
            )
        prefix = "Represent this sentence for searching relevant passages: "
        return self.embedder.embed_text(prefix + text)

    # ══════════════════════════════════════════════════════════════════
    # COLLECTION SEARCH
    # ══════════════════════════════════════════════════════════════════

    def _search_collection(
        self,
        collection_name: str,
        query_vector: List[float],
        top_k: int,
        filter_expr: Optional[str] = None,
    ) -> List[dict]:
        """Search a single Milvus collection.

        Performs a vector similarity search on the specified collection
        with optional scalar field filtering.

        Args:
            collection_name: Milvus collection name.
            query_vector: 384-dimensional query embedding.
            top_k: Maximum number of results.
            filter_expr: Optional Milvus boolean filter expression
                (e.g. 'phase == "Phase 3"').

        Returns:
            List of result dicts from Milvus with score and field values.
        """
        try:
            search_params = {
                "metric_type": "COSINE",
                "params": {"nprobe": 16},
            }

            # Build search kwargs
            search_kwargs = {
                "collection_name": collection_name,
                "data": [query_vector],
                "anns_field": "embedding",
                "param": search_params,
                "limit": top_k,
                "output_fields": ["*"],
            }

            if filter_expr:
                search_kwargs["filter"] = filter_expr

            results = self.milvus.search(**search_kwargs)

            # Flatten Milvus search results
            flat_results = []
            if results and len(results) > 0:
                for hit in results[0]:
                    record = {
                        "id": str(hit.id),
                        "score": float(hit.score) if hasattr(hit, "score") else 0.0,
                    }
                    # Extract entity fields
                    if hasattr(hit, "entity"):
                        entity = hit.entity
                        if hasattr(entity, "fields"):
                            for field_name, field_value in entity.fields.items():
                                if field_name != "embedding":
                                    record[field_name] = field_value
                        elif isinstance(entity, dict):
                            for k, v in entity.items():
                                if k != "embedding":
                                    record[k] = v
                    flat_results.append(record)

            return flat_results

        except Exception as exc:
            logger.warning(
                "Search failed for collection '%s': %s", collection_name, exc,
            )
            return []

    def _parallel_search(
        self,
        query_vector: List[float],
        collections: List[str],
        weights: Dict[str, float],
        top_k: int,
    ) -> List[TrialSearchResult]:
        """Search multiple collections in parallel with weighted scoring.

        Uses ThreadPoolExecutor for concurrent Milvus searches across
        all specified collections. Applies collection-specific weights
        to raw similarity scores for unified ranking.

        Args:
            query_vector: 384-dimensional query embedding.
            collections: List of collection names to search.
            weights: Dict mapping collection name to weight multiplier.
            top_k: Maximum results per collection.

        Returns:
            List of TrialSearchResult sorted by weighted score descending.
        """
        all_results: List[TrialSearchResult] = []
        max_workers = min(len(collections), 8)

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_collection = {
                executor.submit(
                    self._search_collection, coll, query_vector, top_k,
                ): coll
                for coll in collections
            }

            for future in as_completed(future_to_collection):
                coll_name = future_to_collection[future]
                try:
                    raw_results = future.result(timeout=30)
                except Exception as exc:
                    logger.warning(
                        "Parallel search failed for '%s': %s", coll_name, exc,
                    )
                    continue

                cfg = COLLECTION_CONFIG.get(coll_name, {})
                label = cfg.get("label", coll_name)
                weight = weights.get(coll_name, 0.05)
                text_field = cfg.get("text_field", "text_chunk")
                title_field = cfg.get("title_field", "")

                for record in raw_results:
                    raw_score = record.get("score", 0.0)
                    weighted_score = raw_score * (1.0 + weight)

                    # Citation relevance tier
                    if raw_score >= settings.CITATION_HIGH_THRESHOLD:
                        relevance = "high"
                    elif raw_score >= settings.CITATION_MEDIUM_THRESHOLD:
                        relevance = "medium"
                    else:
                        relevance = "low"

                    # Extract text content
                    text = record.get(text_field, "")
                    if not text and title_field:
                        text = record.get(title_field, "")
                    if not text:
                        # Fallback: try common text fields
                        for fallback in ("abstract", "content", "recommendation",
                                         "guidance_text", "results_summary",
                                         "design_summary", "criteria_text"):
                            text = record.get(fallback, "")
                            if text:
                                break

                    # Build metadata (exclude embedding vector)
                    metadata = {
                        k: v for k, v in record.items()
                        if k not in ("embedding",)
                    }
                    metadata["relevance"] = relevance
                    metadata["collection_label"] = label
                    metadata["weight_applied"] = weight

                    result = TrialSearchResult(
                        collection=coll_name,
                        record_id=str(record.get("id", "")),
                        score=weighted_score,
                        text=text,
                        metadata=metadata,
                        relevance=relevance,
                    )
                    all_results.append(result)

        # Sort by weighted score descending
        all_results.sort(key=lambda r: r.score, reverse=True)

        # Deduplicate by record_id
        seen_ids: set = set()
        unique_results: List[TrialSearchResult] = []
        for result in all_results:
            dedup_key = f"{result.collection}:{result.record_id}"
            if dedup_key not in seen_ids:
                seen_ids.add(dedup_key)
                unique_results.append(result)

        # Cap at reasonable limit
        return unique_results[:top_k * len(collections)]

    # ══════════════════════════════════════════════════════════════════
    # RERANKING
    # ══════════════════════════════════════════════════════════════════

    def _rerank_results(
        self,
        results: List[TrialSearchResult],
        query: str,
    ) -> List[TrialSearchResult]:
        """Rerank results based on relevance to original query.

        Applies heuristic boosts for:
        - Regulatory results matching query agencies (FDA, EMA)
        - Results from protocol/endpoint collections for design queries
        - Results with high citation relevance
        - NCT/PMID-bearing results (evidence-based)
        - Results matching detected biomarker or drug terms

        Args:
            results: Raw search results to rerank.
            query: Original query text for relevance matching.

        Returns:
            Reranked list of TrialSearchResult.
        """
        query_lower = query.lower()
        query_terms = set(query_lower.split())

        for result in results:
            boost = 0.0

            # Boost regulatory results when query mentions regulatory bodies
            if result.collection == "trial_regulatory":
                agency = result.metadata.get("agency", "").lower()
                if agency and agency in query_lower:
                    boost += 0.15
                # Boost ICH guidelines
                doc_type = result.metadata.get("document_type", "").lower()
                if "ich" in query_lower and "ich" in doc_type:
                    boost += 0.10

            # Boost guideline results
            if result.collection == "trial_guidelines":
                boost += 0.05

            # Boost results with NCT IDs
            nct_id = result.metadata.get("nct_id", "")
            if nct_id:
                boost += 0.05

            # Boost results with PMIDs
            pmid = result.metadata.get("pmid", "")
            if pmid:
                boost += 0.05

            # Boost results with high relevance
            if result.relevance == "high":
                boost += 0.10
            elif result.relevance == "medium":
                boost += 0.05

            # Boost protocol results for design queries
            design_terms = {"protocol", "design", "phase", "randomized",
                            "sample size", "power", "endpoint"}
            if result.collection == "trial_protocols":
                if query_terms & design_terms:
                    boost += 0.10

            # Boost endpoint results for endpoint queries
            endpoint_terms = {"endpoint", "primary", "secondary", "surrogate",
                              "outcome", "overall survival", "pfs", "orr"}
            if result.collection == "trial_endpoints":
                if query_terms & endpoint_terms:
                    boost += 0.10

            # Boost safety results for safety queries
            safety_terms = {"safety", "adverse", "sae", "dsmb", "toxicity",
                            "susar", "dlt", "crs"}
            if result.collection == "trial_safety":
                if query_terms & safety_terms:
                    boost += 0.10

            # Boost biomarker results for precision queries
            biomarker_terms = {"biomarker", "companion", "diagnostic",
                               "ctdna", "genomic", "molecular", "precision"}
            if result.collection == "trial_biomarkers":
                if query_terms & biomarker_terms:
                    boost += 0.10

            # Boost eligibility results for criteria queries
            eligibility_terms = {"eligibility", "inclusion", "exclusion",
                                 "criteria", "enroll", "screen"}
            if result.collection == "trial_eligibility":
                if query_terms & eligibility_terms:
                    boost += 0.10

            # Apply boost
            result.score += boost

        # Re-sort after boosting
        results.sort(key=lambda r: r.score, reverse=True)
        return results

    # ══════════════════════════════════════════════════════════════════
    # CITATION SCORING
    # ══════════════════════════════════════════════════════════════════

    def _score_citations(
        self,
        results: List[TrialSearchResult],
    ) -> List[TrialSearchResult]:
        """Score and label results with citation relevance tiers.

        Assigns high/medium/low relevance based on raw similarity score
        thresholds from settings.

        Args:
            results: Search results to score.

        Returns:
            Same list with updated relevance fields.
        """
        for result in results:
            raw_score = result.metadata.get("score", result.score)
            if raw_score >= settings.CITATION_HIGH_THRESHOLD:
                result.relevance = "high"
            elif raw_score >= settings.CITATION_MEDIUM_THRESHOLD:
                result.relevance = "medium"
            else:
                result.relevance = "low"
            result.metadata["relevance"] = result.relevance
        return results

    # ══════════════════════════════════════════════════════════════════
    # LLM SYNTHESIS
    # ══════════════════════════════════════════════════════════════════

    def _synthesize_response(
        self,
        question: str,
        results: List[TrialSearchResult],
        workflow: Optional[TrialWorkflowType] = None,
        patient_context: Optional[dict] = None,
    ) -> TrialResponse:
        """Use LLM to synthesise search results into a trial intelligence response.

        Builds a structured prompt with retrieved evidence, patient context,
        conversation history, and workflow-specific instructions. Generates
        a grounded answer via the configured LLM.

        Args:
            question: Original trial question.
            results: Ranked search results for context.
            workflow: Optional workflow for instruction tuning.
            patient_context: Optional patient-specific data dict.

        Returns:
            TrialResponse with synthesised answer and metadata.
        """
        context = self._build_context(results, patient_context)
        patient_section = self._format_patient_context(patient_context)
        conversation_section = self._format_conversation_history()
        workflow_section = self._format_workflow_instructions(workflow)

        prompt = (
            f"## Retrieved Evidence\n\n{context}\n\n"
            f"{patient_section}"
            f"{conversation_section}"
            f"{workflow_section}"
            f"---\n\n"
            f"## Question\n\n{question}\n\n"
            f"Please provide a comprehensive, evidence-based clinical trial "
            f"intelligence answer grounded in the retrieved evidence above. "
            f"Follow the system prompt instructions for trial citation format, "
            f"severity badges, statistical context, and structured output sections.\n\n"
            f"Cite sources using clickable markdown links where NCT IDs are available: "
            f"[NCT01234567](https://clinicaltrials.gov/study/NCT01234567). "
            f"For PubMed evidence, use [PMID:12345678](https://pubmed.ncbi.nlm.nih.gov/12345678/). "
            f"For collection-sourced evidence, use [Collection:record-id]. "
            f"Prioritise [high relevance] citations."
        )

        answer = self.llm.generate(
            prompt=prompt,
            system_prompt=TRIAL_SYSTEM_PROMPT,
            max_tokens=2048,
            temperature=0.7,
        )

        return TrialResponse(
            question=question,
            answer=answer,
            results=results,
            workflow=workflow,
        )

    def _build_context(
        self,
        results: List[TrialSearchResult],
        patient_context: Optional[dict] = None,
    ) -> str:
        """Build context string from search results for LLM prompt.

        Organises results by collection, formatting each with its
        citation reference, relevance tag, score, and text excerpt.

        Args:
            results: Ranked search results to format.
            patient_context: Optional patient context (used for additional
                context augmentation).

        Returns:
            Formatted evidence context string for the LLM prompt.
        """
        if not results:
            return "No evidence found in the knowledge base."

        # Group results by collection
        by_collection: Dict[str, List[TrialSearchResult]] = {}
        for result in results:
            label = result.metadata.get("collection_label", result.collection)
            if label not in by_collection:
                by_collection[label] = []
            by_collection[label].append(result)

        sections: List[str] = []
        for label, coll_results in by_collection.items():
            section_lines = [f"### Evidence from {label}"]
            for i, result in enumerate(coll_results[:5], 1):
                citation = self._format_citation_link(result)
                relevance_tag = (
                    f" [{result.relevance} relevance]"
                    if result.relevance else ""
                )
                text_excerpt = result.text[:500] if result.text else "(no text)"
                section_lines.append(
                    f"{i}. {citation}{relevance_tag} "
                    f"(score={result.score:.3f}) {text_excerpt}"
                )
            sections.append("\n".join(section_lines))

        return "\n\n".join(sections)

    def _format_citation_link(self, result: TrialSearchResult) -> str:
        """Format a citation with clickable URL where possible.

        Args:
            result: Search result to format citation for.

        Returns:
            Markdown-formatted citation string.
        """
        label = result.metadata.get("collection_label", result.collection)
        record_id = result.record_id

        # ClinicalTrials.gov
        nct_id = result.metadata.get("nct_id", "")
        if nct_id:
            return (
                f"[{label}:{nct_id}]"
                f"(https://clinicaltrials.gov/study/{nct_id})"
            )

        # PubMed literature
        pmid = result.metadata.get("pmid", "")
        if pmid:
            return (
                f"[{label}:PMID {pmid}]"
                f"(https://pubmed.ncbi.nlm.nih.gov/{pmid}/)"
            )

        # FDA guidance
        fda_doc = result.metadata.get("fda_document_id", "")
        if fda_doc:
            return f"[{label}:FDA {fda_doc}]"

        return f"[{label}:{record_id}]"

    def _format_patient_context(self, patient_context: Optional[dict]) -> str:
        """Format patient context for LLM prompt injection.

        Used primarily for patient-trial matching workflows.

        Args:
            patient_context: Optional patient data dict with keys like
                age, sex, diagnosis, biomarkers, prior_therapies, mutations,
                performance_status, stage, comorbidities.

        Returns:
            Formatted patient context section or empty string.
        """
        if not patient_context:
            return ""

        lines = ["### Patient Context\n"]

        field_labels = {
            "age": "Age",
            "sex": "Sex",
            "diagnosis": "Primary Diagnosis",
            "stage": "Disease Stage",
            "histology": "Histology",
            "performance_status": "ECOG Performance Status",
            "biomarkers": "Biomarker Status",
            "mutations": "Genomic Alterations",
            "prior_therapies": "Prior Lines of Therapy",
            "prior_lines": "Number of Prior Lines",
            "comorbidities": "Comorbidities",
            "medications": "Current Medications",
            "organ_function": "Organ Function",
            "labs": "Recent Labs",
            "imaging": "Imaging Findings",
            "allergies": "Allergies",
            "geographic_region": "Geographic Region",
            "language": "Preferred Language",
            "travel_distance_km": "Max Travel Distance (km)",
        }

        for key, label in field_labels.items():
            value = patient_context.get(key)
            if value is not None:
                if isinstance(value, list):
                    value = ", ".join(str(v) for v in value)
                elif isinstance(value, dict):
                    value = "; ".join(f"{k}: {v}" for k, v in value.items())
                lines.append(f"- **{label}:** {value}")

        lines.append("\n")
        return "\n".join(lines)

    def _format_conversation_history(self) -> str:
        """Format recent conversation history for multi-turn context.

        Returns:
            Formatted conversation history section or empty string.
        """
        if not self._conversation_history:
            return ""

        # Use only the most recent exchanges
        recent = self._conversation_history[-self._max_conversation_context * 2:]

        lines = ["### Conversation History\n"]
        for entry in recent:
            role = entry.get("role", "unknown").capitalize()
            content = entry.get("content", "")[:300]
            lines.append(f"**{role}:** {content}")

        lines.append("\n")
        return "\n".join(lines)

    def _format_workflow_instructions(
        self,
        workflow: Optional[TrialWorkflowType],
    ) -> str:
        """Format workflow-specific instructions for the LLM prompt.

        Args:
            workflow: Optional workflow type for tailored instructions.

        Returns:
            Workflow instruction section or empty string.
        """
        if not workflow:
            return ""

        instructions = {
            TrialWorkflowType.PROTOCOL_DESIGN: (
                "### Workflow: Protocol Design\n"
                "Focus on: study design rationale, randomization scheme, blinding strategy, "
                "sample size justification, statistical analysis plan, key design elements "
                "(primary endpoint, treatment duration, follow-up), and regulatory precedent "
                "for similar designs. Reference ICH E8(R1) and ICH E9(R1) where applicable.\n\n"
            ),
            TrialWorkflowType.PATIENT_MATCHING: (
                "### Workflow: Patient-Trial Matching\n"
                "Focus on: eligibility criteria alignment with patient profile, biomarker "
                "requirements, prior therapy restrictions, washout periods, performance status "
                "requirements, geographic accessibility, and trial availability. Highlight "
                "potential barriers and suggest alternative trials if primary match fails.\n\n"
            ),
            TrialWorkflowType.SITE_SELECTION: (
                "### Workflow: Site Selection\n"
                "Focus on: therapeutic area expertise, patient access and demographics, "
                "enrollment track record, regulatory infrastructure, investigator experience, "
                "geographic diversity requirements, decentralized trial capabilities, and "
                "site performance metrics.\n\n"
            ),
            TrialWorkflowType.ELIGIBILITY_ANALYSIS: (
                "### Workflow: Eligibility Criteria Analysis\n"
                "Focus on: criteria stringency vs enrollment feasibility, impact on "
                "generalizability, screen failure rate prediction, diversity implications, "
                "biomarker testing requirements, and comparison with similar trial criteria. "
                "Suggest optimizations to balance scientific rigor with enrollment speed.\n\n"
            ),
            TrialWorkflowType.ENDPOINT_STRATEGY: (
                "### Workflow: Endpoint Strategy\n"
                "Focus on: primary endpoint selection and regulatory acceptance, surrogate "
                "endpoint validation status, composite endpoint construction, PRO instruments, "
                "effect size estimation, multiplicity considerations, and precedent endpoints "
                "in successful programs. Reference ICH E9(R1) estimand framework.\n\n"
            ),
            TrialWorkflowType.REGULATORY_STRATEGY: (
                "### Workflow: Regulatory Strategy\n"
                "Focus on: applicable regulatory pathways (Breakthrough, Fast Track, "
                "Accelerated Approval, Priority Review, PRIME, Orphan), pre-IND/pre-BLA "
                "meeting strategy, Type A/B/C meeting requests, regulatory precedent analysis, "
                "labeling strategy, and post-marketing commitments. Cite specific FDA guidance "
                "documents and ICH guidelines.\n\n"
            ),
            TrialWorkflowType.COMPETITIVE_INTELLIGENCE: (
                "### Workflow: Competitive Intelligence\n"
                "Focus on: competitive landscape mapping, mechanism differentiation, "
                "trial design comparison, enrollment timeline analysis, interim results "
                "interpretation, regulatory pathway comparison, and strategic positioning. "
                "Provide structured comparison tables where possible.\n\n"
            ),
            TrialWorkflowType.SAFETY_MONITORING: (
                "### Workflow: Safety Monitoring\n"
                "Focus on: DSMB charter provisions, stopping rules (efficacy and futility), "
                "safety signal detection methodology, SUSAR reporting requirements, risk "
                "management plan, REMS considerations, dose modification algorithms, and "
                "safety database review timeline. Reference ICH E2A/E2B/E6(R3).\n\n"
            ),
            TrialWorkflowType.ADAPTIVE_DESIGN: (
                "### Workflow: Adaptive Design\n"
                "Focus on: adaptation type selection (sample size re-estimation, interim "
                "analysis, dose-finding, enrichment, seamless Phase 2/3), statistical "
                "methodology (Bayesian, frequentist, hybrid), type I error control, "
                "operational complexity, and regulatory acceptance. Reference ICH E20 and "
                "FDA Adaptive Design Guidance (2019).\n\n"
            ),
            TrialWorkflowType.BIOMARKER_STRATEGY: (
                "### Workflow: Biomarker Strategy\n"
                "Focus on: biomarker validation level (exploratory, probable, known valid), "
                "companion diagnostic development timeline, assay standardization, "
                "enrichment design vs all-comers, co-development with therapeutic, "
                "regulatory pathway for CDx, and real-world testing feasibility.\n\n"
            ),
            TrialWorkflowType.RWE_ANALYSIS: (
                "### Workflow: Real-World Evidence Analysis\n"
                "Focus on: RWE study design (retrospective, prospective, hybrid), data "
                "source selection (claims, EHR, registries), endpoint definitions in RWD, "
                "confounding control methodology, regulatory acceptability of RWE, and "
                "integration with clinical trial evidence. Reference FDA RWE Framework.\n\n"
            ),
            TrialWorkflowType.RECRUITMENT_OPTIMIZATION: (
                "### Workflow: Recruitment Optimization\n"
                "Focus on: enrollment rate projections, diversity and inclusion strategy, "
                "digital recruitment channels, site-level enrollment targets, screen failure "
                "mitigation, retention strategies, decentralized trial elements, and "
                "community engagement approaches.\n\n"
            ),
        }

        return instructions.get(workflow, "")

    # ══════════════════════════════════════════════════════════════════
    # CITATIONS & CONFIDENCE
    # ══════════════════════════════════════════════════════════════════

    def _extract_citations(
        self,
        results: List[TrialSearchResult],
    ) -> List[dict]:
        """Extract and format citations from search results.

        Generates a structured citation list from all results, including
        NCT links, PMID links, and regulatory document references.

        Args:
            results: Search results to extract citations from.

        Returns:
            List of citation dicts with keys: source, id, title, url,
            relevance, score.
        """
        citations: List[dict] = []
        seen: set = set()

        for result in results:
            cite = {
                "source": result.metadata.get("collection_label", result.collection),
                "id": result.record_id,
                "title": "",
                "url": "",
                "relevance": result.relevance,
                "score": round(result.score, 4),
            }

            # Extract title from metadata
            cfg = COLLECTION_CONFIG.get(result.collection, {})
            title_field = cfg.get("title_field", "")
            if title_field:
                cite["title"] = result.metadata.get(title_field, "")

            # Generate URL for known reference types
            nct_id = result.metadata.get("nct_id", "")
            if nct_id:
                cite["url"] = f"https://clinicaltrials.gov/study/{nct_id}"
                cite["id"] = nct_id

            pmid = result.metadata.get("pmid", "")
            if pmid:
                cite["url"] = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                cite["id"] = f"PMID:{pmid}"

            doi = result.metadata.get("doi", "")
            if doi and not cite["url"]:
                cite["url"] = f"https://doi.org/{doi}"

            # Deduplicate
            dedup_key = cite["id"] or f"{cite['source']}:{result.record_id}"
            if dedup_key not in seen:
                seen.add(dedup_key)
                citations.append(cite)

        return citations

    def _score_confidence(
        self,
        results: List[TrialSearchResult],
    ) -> float:
        """Score overall confidence based on result quality.

        Confidence is based on:
        - Number of high-relevance results
        - Collection diversity
        - Average similarity score
        - Presence of regulatory/guideline evidence

        Args:
            results: Search results to evaluate.

        Returns:
            Confidence score between 0.0 and 1.0.
        """
        if not results:
            return 0.0

        # Factor 1: High-relevance ratio (0-0.3)
        high_count = sum(1 for r in results if r.relevance == "high")
        relevance_score = min(high_count / max(len(results), 1), 1.0) * 0.3

        # Factor 2: Collection diversity (0-0.3)
        unique_collections = len(set(r.collection for r in results))
        len(COLLECTION_CONFIG)
        diversity_score = min(unique_collections / 4, 1.0) * 0.3

        # Factor 3: Average score of top results (0-0.2)
        top_scores = [r.score for r in results[:5]]
        avg_score = sum(top_scores) / max(len(top_scores), 1)
        quality_score = min(avg_score, 1.0) * 0.2

        # Factor 4: Regulatory/guideline evidence present (0-0.2)
        has_regulatory = any(
            r.collection in ("trial_regulatory", "trial_guidelines")
            for r in results
        )
        regulatory_score = 0.2 if has_regulatory else 0.0

        confidence = relevance_score + diversity_score + quality_score + regulatory_score
        return round(min(confidence, 1.0), 3)

    # ══════════════════════════════════════════════════════════════════
    # ENTITY & TRIAL SEARCH
    # ══════════════════════════════════════════════════════════════════

    def find_related(
        self,
        entity: str,
        entity_type: str = "condition",
        top_k: int = 5,
    ) -> List[TrialSearchResult]:
        """Find related entities across collections.

        Searches all collections for evidence related to a clinical
        entity (condition, drug, biomarker, regulatory body). Useful
        for building entity profiles and cross-referencing.

        Args:
            entity: Entity name (e.g. 'NSCLC', 'pembrolizumab', 'PD-L1').
            entity_type: Entity category for targeted search:
                'condition', 'drug', 'biomarker', 'regulatory', 'investigator'.
            top_k: Maximum results per collection.

        Returns:
            List of TrialSearchResult from all relevant collections.
        """
        type_collection_map = {
            "condition": [
                "trial_protocols", "trial_results", "trial_literature",
                "trial_eligibility", "trial_guidelines",
            ],
            "drug": [
                "trial_results", "trial_protocols", "trial_literature",
                "trial_safety", "trial_regulatory",
            ],
            "biomarker": [
                "trial_biomarkers", "genomic_evidence", "trial_eligibility",
                "trial_literature", "trial_endpoints",
            ],
            "regulatory": [
                "trial_regulatory", "trial_guidelines", "trial_protocols",
            ],
            "investigator": [
                "trial_investigators", "trial_sites", "trial_literature",
            ],
        }

        collections = type_collection_map.get(entity_type, get_all_collection_names())
        return self.search(entity, collections=collections, top_k=top_k)

    def get_trial_details(self, trial_id: str) -> Optional[dict]:
        """Retrieve details for a specific trial by NCT ID.

        Searches the trial_protocols collection for a specific trial
        using scalar filtering on the nct_id field.

        Args:
            trial_id: ClinicalTrials.gov NCT identifier (e.g. 'NCT01234567').

        Returns:
            Trial details dict or None if not found.
        """
        if not self.milvus:
            raise RuntimeError("Milvus client not configured.")

        # Sanitize input
        safe_id = trial_id.strip()
        if not _SAFE_FILTER_RE.match(safe_id):
            logger.warning("Rejected unsafe trial ID: %r", safe_id)
            return None

        try:
            # Use a generic query vector for filtered search
            query_vector = self._embed_query(f"clinical trial {trial_id}")
            filter_expr = f'nct_id == "{safe_id}"'

            raw_results = self._search_collection(
                "trial_protocols", query_vector, top_k=1,
                filter_expr=filter_expr,
            )

            if raw_results:
                return raw_results[0]
            return None

        except Exception as exc:
            logger.warning("Failed to get trial details for %s: %s", trial_id, exc)
            return None

    def search_eligibility(
        self,
        trial_id: str,
        patient_profile: dict,
    ) -> dict:
        """Assess patient eligibility for a specific trial.

        Retrieves the trial's eligibility criteria and evaluates them
        against the provided patient profile. Returns a structured
        assessment with matched/unmatched criteria.

        Args:
            trial_id: ClinicalTrials.gov NCT identifier.
            patient_profile: Patient data dict with keys like age, sex,
                diagnosis, biomarkers, prior_therapies, performance_status.

        Returns:
            Dict with keys: trial_id, eligible (bool), matched_criteria (list),
            unmatched_criteria (list), uncertain_criteria (list), assessment (str).
        """
        if not self.milvus:
            raise RuntimeError("Milvus client not configured.")

        assessment = {
            "trial_id": trial_id,
            "eligible": False,
            "matched_criteria": [],
            "unmatched_criteria": [],
            "uncertain_criteria": [],
            "assessment": "",
        }

        # Get trial eligibility criteria
        safe_id = trial_id.strip()
        if not _SAFE_FILTER_RE.match(safe_id):
            assessment["assessment"] = "Invalid trial ID format."
            return assessment

        try:
            query_vector = self._embed_query(f"eligibility criteria {trial_id}")
            filter_expr = f'nct_id == "{safe_id}"'

            criteria_results = self._search_collection(
                "trial_eligibility", query_vector, top_k=10,
                filter_expr=filter_expr,
            )

            if not criteria_results:
                assessment["assessment"] = f"No eligibility criteria found for {trial_id}."
                return assessment

            # If LLM available, use it for assessment
            if self.llm:
                criteria_text = "\n".join(
                    r.get("criteria_text", str(r)) for r in criteria_results
                )
                patient_text = self._format_patient_context(patient_profile)

                prompt = (
                    f"## Trial Eligibility Criteria ({trial_id})\n\n"
                    f"{criteria_text}\n\n"
                    f"{patient_text}\n\n"
                    f"Evaluate whether this patient meets the eligibility criteria. "
                    f"For each criterion, indicate: MET, NOT MET, or UNCERTAIN. "
                    f"Provide an overall eligibility assessment."
                )

                answer = self.llm.generate(
                    prompt=prompt,
                    system_prompt=TRIAL_SYSTEM_PROMPT,
                    max_tokens=1024,
                    temperature=0.3,
                )
                assessment["assessment"] = answer
                assessment["eligible"] = "NOT MET" not in answer.upper()[:200]
            else:
                assessment["assessment"] = (
                    "LLM not configured. Manual review of criteria required."
                )

        except Exception as exc:
            logger.warning(
                "Eligibility assessment failed for %s: %s", trial_id, exc,
            )
            assessment["assessment"] = f"Assessment error: {exc}"

        return assessment

    # ══════════════════════════════════════════════════════════════════
    # CONVERSATION MEMORY
    # ══════════════════════════════════════════════════════════════════

    def add_conversation_context(
        self,
        role: str,
        content: str,
        session_id: Optional[str] = None,
    ):
        """Add to conversation history for multi-turn context.

        Maintains a rolling window of recent conversation exchanges
        for follow-up query context injection. Persists to disk so
        history survives restarts.

        Args:
            role: Message role ('user' or 'assistant').
            content: Message content text.
            session_id: Optional override; defaults to self.session_id.
        """
        self._conversation_history.append({
            "role": role,
            "content": content,
            "timestamp": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        })

        # Trim to max context window
        max_entries = self._max_conversation_context * 2
        if len(self._conversation_history) > max_entries:
            self._conversation_history = self._conversation_history[-max_entries:]

        # Persist to disk
        _save_conversation(session_id or self.session_id, self._conversation_history)

    def clear_conversation(self, session_id: Optional[str] = None):
        """Clear conversation history.

        Resets the multi-turn context and removes the persisted file.
        Useful when starting a new consultation or switching topics.

        Args:
            session_id: Optional override; defaults to self.session_id.
        """
        self._conversation_history.clear()
        sid = session_id or self.session_id
        try:
            path = CONVERSATION_DIR / f"{sid}.json"
            if path.exists():
                path.unlink()
        except Exception as exc:
            logger.warning("Failed to remove conversation file %s: %s", sid, exc)

    # ══════════════════════════════════════════════════════════════════
    # WEIGHT COMPUTATION
    # ══════════════════════════════════════════════════════════════════

    def _get_boosted_weights(
        self,
        workflow: Optional[TrialWorkflowType] = None,
    ) -> Dict[str, float]:
        """Compute collection weights with optional workflow boosting.

        When a workflow is specified, applies boost multipliers from
        WORKFLOW_COLLECTION_BOOST on top of the base weights from
        settings. Weights are then renormalized to sum to ~1.0.

        Args:
            workflow: Optional TrialWorkflowType for boosting.

        Returns:
            Dict mapping collection name to adjusted weight.
        """
        # Start with base weights
        base_weights = {
            name: cfg.get("weight", 0.05)
            for name, cfg in COLLECTION_CONFIG.items()
        }

        if not workflow or workflow not in WORKFLOW_COLLECTION_BOOST:
            return base_weights

        # Apply boost multipliers
        boosts = WORKFLOW_COLLECTION_BOOST[workflow]
        boosted = {}
        for name, base_w in base_weights.items():
            multiplier = boosts.get(name, 1.0)
            boosted[name] = base_w * multiplier

        # Renormalize to sum to ~1.0
        total = sum(boosted.values())
        if total > 0:
            boosted = {name: w / total for name, w in boosted.items()}

        return boosted

    # ══════════════════════════════════════════════════════════════════
    # HEALTH CHECK
    # ══════════════════════════════════════════════════════════════════

    def health_check(self) -> dict:
        """Check Milvus connection and collection status.

        Verifies connectivity to the Milvus server and checks that
        all expected trial collections exist and are loaded.

        Returns:
            Dict with keys: status ('healthy'/'degraded'/'unhealthy'),
            milvus_connected (bool), collections_available (list),
            collections_missing (list), embedding_model (str),
            llm_configured (bool).
        """
        health = {
            "status": "unhealthy",
            "milvus_connected": False,
            "collections_available": [],
            "collections_missing": [],
            "embedding_model": settings.EMBEDDING_MODEL,
            "llm_configured": self.llm is not None,
            "timestamp": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        }

        if not self.milvus:
            health["error"] = "Milvus client not configured"
            return health

        try:
            available_collections = []
            expected_names = get_all_collection_names()

            for coll_name in expected_names:
                try:
                    has_collection = self.milvus.has_collection(coll_name)
                    if has_collection:
                        available_collections.append(coll_name)
                    else:
                        health["collections_missing"].append(coll_name)
                except Exception:
                    health["collections_missing"].append(coll_name)

            health["milvus_connected"] = True
            health["collections_available"] = available_collections

            total_expected = len(expected_names)
            total_available = len(available_collections)

            if total_available == total_expected:
                health["status"] = "healthy"
            elif total_available >= total_expected * 0.5:
                health["status"] = "degraded"
            else:
                health["status"] = "unhealthy"

        except Exception as exc:
            health["error"] = str(exc)
            health["status"] = "unhealthy"

        return health
