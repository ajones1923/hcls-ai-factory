"""Multi-collection RAG engine for Cardiology Intelligence Agent.

Searches across all 12 cardiology-specific Milvus collections simultaneously
using parallel ThreadPoolExecutor, synthesises findings with clinical knowledge
augmentation, and generates grounded LLM responses with guideline citations.

Extends the pattern from: rag-chat-pipeline/src/rag_engine.py

Features:
- Parallel search via ThreadPoolExecutor (12 cardiology + 1 shared genomic collection)
- Settings-driven weights and parameters from config/settings.py
- Workflow-based dynamic weight boosting per CardioWorkflowType
- Milvus field-based filtering (modality, hf_type, valve, severity, condition)
- Citation relevance scoring (high/medium/low) with PMID link formatting
- Cross-collection entity linking for comprehensive cardiovascular queries
- Guideline retrieval with Class of Recommendation / Level of Evidence
- Conversation memory for multi-turn clinical consultations
- Patient context injection for personalised recommendations
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
from typing import Any, Dict, Generator, List, Optional

from config.settings import settings

from .models import (
    CardioWorkflowType,
)
from .collections import (
    get_all_collection_names,
)
from .agent import (
    CARDIO_SYSTEM_PROMPT,
    WORKFLOW_COLLECTION_BOOST,
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
class CardioSearchResult:
    """A single search result from a Milvus collection.

    Attributes:
        collection: Source collection name (e.g. 'cardio_guidelines').
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


@dataclass
class CardioResponse:
    """Complete response from the RAG engine.

    Attributes:
        question: Original query text.
        answer: LLM-synthesised answer text.
        results: Ranked search results used for synthesis.
        workflow: Clinical workflow that was applied.
        confidence: Overall confidence score (0.0 - 1.0).
        citations: Formatted citation list.
        search_time_ms: Total search time in milliseconds.
        collections_searched: Number of collections queried.
        patient_context_used: Whether patient context was injected.
        timestamp: ISO 8601 timestamp of response generation.
    """
    question: str = ""
    answer: str = ""
    results: List[CardioSearchResult] = field(default_factory=list)
    workflow: Optional[CardioWorkflowType] = None
    confidence: float = 0.0
    citations: List[Dict[str, str]] = field(default_factory=list)
    search_time_ms: float = 0.0
    collections_searched: int = 0
    patient_context_used: bool = False
    timestamp: str = field(
        default_factory=lambda: datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    )


# =====================================================================
# COLLECTION CONFIGURATION (reads weights from settings)
# =====================================================================

COLLECTION_CONFIG: Dict[str, Dict[str, Any]] = {
    "cardio_literature": {
        "weight": settings.WEIGHT_LITERATURE,
        "label": "Literature",
        "text_field": "abstract",
        "title_field": "title",
        "filterable_fields": ["subspecialty", "study_type"],
    },
    "cardio_trials": {
        "weight": settings.WEIGHT_TRIALS,
        "label": "Trial",
        "text_field": "key_findings",
        "title_field": "trial_name",
        "filterable_fields": ["condition", "status", "phase"],
    },
    "cardio_imaging": {
        "weight": settings.WEIGHT_IMAGING,
        "label": "Imaging",
        "text_field": "clinical_significance",
        "title_field": "finding",
        "filterable_fields": ["modality"],
    },
    "cardio_electrophysiology": {
        "weight": settings.WEIGHT_ELECTROPHYSIOLOGY,
        "label": "Electrophysiology",
        "text_field": "clinical_significance",
        "title_field": "finding",
        "filterable_fields": ["category", "urgency"],
    },
    "cardio_heart_failure": {
        "weight": settings.WEIGHT_HEART_FAILURE,
        "label": "HeartFailure",
        "text_field": "content",
        "title_field": "topic",
        "filterable_fields": ["hf_type", "nyha_class", "acc_stage", "drug_class"],
    },
    "cardio_valvular": {
        "weight": settings.WEIGHT_VALVULAR,
        "label": "Valvular",
        "text_field": "indication",
        "title_field": "valve",
        "filterable_fields": ["valve", "pathology", "severity"],
    },
    "cardio_prevention": {
        "weight": settings.WEIGHT_PREVENTION,
        "label": "Prevention",
        "text_field": "content",
        "title_field": "topic",
        "filterable_fields": ["risk_factor", "evidence_class"],
    },
    "cardio_interventional": {
        "weight": settings.WEIGHT_INTERVENTIONAL,
        "label": "Interventional",
        "text_field": "outcomes",
        "title_field": "procedure_name",
        "filterable_fields": [],
    },
    "cardio_oncology": {
        "weight": settings.WEIGHT_ONCOLOGY,
        "label": "CardioOncology",
        "text_field": "management",
        "title_field": "chemotherapy_agent",
        "filterable_fields": ["risk_level", "cardiotoxicity_type"],
    },
    "cardio_devices": {
        "weight": settings.WEIGHT_DEVICES,
        "label": "Devices",
        "text_field": "clinical_application",
        "title_field": "device_name",
        "filterable_fields": ["device_type", "fda_status"],
    },
    "cardio_guidelines": {
        "weight": settings.WEIGHT_GUIDELINES,
        "label": "Guideline",
        "text_field": "recommendation",
        "title_field": "guideline_title",
        "filterable_fields": ["society", "class_of_rec", "evidence_level", "condition"],
    },
    "cardio_hemodynamics": {
        "weight": settings.WEIGHT_HEMODYNAMICS,
        "label": "Hemodynamics",
        "text_field": "clinical_significance",
        "title_field": "parameter_name",
        "filterable_fields": ["cathlab_context"],
    },
    "genomic_evidence": {
        "weight": settings.WEIGHT_GENOMIC,
        "label": "Genomic",
        "text_field": "text_chunk",
        "title_field": "gene",
        "filterable_fields": [],
    },
}


# =====================================================================
# CARDIO RAG ENGINE
# =====================================================================

class CardioRAGEngine:
    """Multi-collection RAG engine for cardiovascular clinical intelligence.

    Searches across all 12 cardiology-specific Milvus collections plus the
    shared genomic_evidence collection. Supports workflow-specific weight
    boosting, parallel search, query expansion, patient context injection,
    and multi-turn conversation memory.

    Features:
    - Parallel search via ThreadPoolExecutor (13 collections)
    - Settings-driven weights and parameters
    - Workflow-based dynamic weight boosting (8 clinical workflows)
    - Milvus field-based filtering (modality, hf_type, valve, condition, etc.)
    - Citation relevance scoring (high/medium/low)
    - Cross-collection entity linking
    - Guideline retrieval with COR/LOE
    - Conversation memory context injection
    - Patient context for personalised recommendations
    - Confidence scoring based on evidence diversity

    Usage:
        engine = CardioRAGEngine(milvus_client, embedding_model, llm_client)
        response = engine.query("What is GDMT for HFrEF?")
        results = engine.search("aortic stenosis criteria")
    """

    def __init__(
        self,
        milvus_client=None,
        embedding_model=None,
        llm_client=None,
        query_expander=None,
        session_id: str = "default",
    ):
        """Initialize the CardioRAGEngine.

        Args:
            milvus_client: Connected Milvus client with access to all
                cardiology collections. If None, search operations will
                raise RuntimeError.
            embedding_model: Embedding model (BGE-small-en-v1.5) for query
                vectorisation. If None, embedding operations will raise.
            llm_client: LLM client (Anthropic Claude) for response synthesis.
                If None, search-only mode is available.
            query_expander: Optional QueryExpander instance for synonym and
                term expansion.
            session_id: Conversation session identifier for persistence
                (default: "default").
        """
        self.milvus = milvus_client
        self.embedder = embedding_model
        self.llm = llm_client
        self.query_expander = query_expander
        self.session_id = session_id
        self._max_conversation_context = settings.MAX_CONVERSATION_CONTEXT

        # Load persisted conversation history (falls back to empty list)
        self.conversation_history: List[Dict[str, str]] = _load_conversation(session_id)

        # Cleanup expired conversations on startup
        _cleanup_expired_conversations()

    # ══════════════════════════════════════════════════════════════════
    # PUBLIC API
    # ══════════════════════════════════════════════════════════════════

    def query(
        self,
        question: str,
        workflow: Optional[CardioWorkflowType] = None,
        top_k: int = 5,
        patient_context: Optional[dict] = None,
    ) -> CardioResponse:
        """Main query method: expand -> search -> synthesise.

        Performs the full RAG pipeline: query expansion, parallel multi-collection
        search with workflow-specific weighting, result reranking, LLM synthesis
        with patient context, and confidence scoring.

        Args:
            question: Natural language clinical question.
            workflow: Optional CardioWorkflowType to apply domain-specific
                collection weight boosting. If None, weights are auto-detected
                or base defaults are used.
            top_k: Maximum results to return per collection.
            patient_context: Optional dict with patient-specific data
                (age, sex, comorbidities, medications, labs, imaging)
                for personalised recommendations.

        Returns:
            CardioResponse with synthesised answer, search results, citations,
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

        # Step 4: Score confidence
        confidence = self._score_confidence(results)

        # Step 5: Synthesise LLM response (if LLM available)
        if self.llm:
            response = self._synthesize(
                question=question,
                results=results,
                workflow=workflow,
                patient_context=patient_context,
            )
        else:
            response = CardioResponse(
                question=question,
                answer="[LLM not configured -- search-only mode. "
                       "See results below for retrieved evidence.]",
                results=results,
                workflow=workflow,
                confidence=confidence,
            )

        # Step 6: Extract citations
        response.citations = self._extract_citations(results)
        response.confidence = confidence
        response.search_time_ms = (time.time() - start) * 1000
        response.collections_searched = len(collections)
        response.patient_context_used = patient_context is not None

        # Step 7: Update conversation history
        self.add_conversation_context("user", question)
        if response.answer:
            self.add_conversation_context("assistant", response.answer[:500])

        return response

    def search(
        self,
        question: str,
        collections: Optional[List[str]] = None,
        top_k: int = 5,
    ) -> List[CardioSearchResult]:
        """Search across multiple collections with weighted scoring.

        Embeds the query, runs parallel Milvus searches across all specified
        collections, applies collection weights, and returns merged ranked
        results.

        Args:
            question: Natural language search query.
            collections: Optional list of collection names to search.
                If None, all 13 collections are searched.
            top_k: Maximum results per collection.

        Returns:
            List of CardioSearchResult sorted by weighted score descending.
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
            # Add shared genomic collection
            if "genomic_evidence" not in collections:
                collections.append("genomic_evidence")

        # Get weights (base defaults for search-only calls)
        weights = {
            name: COLLECTION_CONFIG.get(name, {}).get("weight", 0.05)
            for name in collections
        }

        # Parallel search with weighting
        results = self._parallel_search(collections, query_vector, top_k, weights)

        return results

    def query_stream(
        self,
        question: str,
        workflow: Optional[CardioWorkflowType] = None,
        top_k: int = 5,
        patient_context: Optional[dict] = None,
    ) -> Generator[Dict, None, None]:
        """Streaming RAG query -- yields evidence then token chunks.

        First yields the search results, then streams LLM tokens as they
        are generated. Useful for real-time UI rendering.

        Args:
            question: Natural language clinical question.
            workflow: Optional workflow for weight boosting.
            top_k: Maximum results per collection.
            patient_context: Optional patient-specific context dict.

        Yields:
            Dicts with 'type' key: 'evidence', 'token', or 'done'.
        """
        # Step 1: Search
        weights = self._get_boosted_weights(workflow)
        collections = list(weights.keys())
        results = self.search(question, collections, top_k)
        results = self._rerank_results(results, question)

        yield {"type": "evidence", "content": results}

        # Step 2: Build prompt and stream LLM response
        if not self.llm:
            yield {
                "type": "done",
                "content": "[LLM not configured -- search-only mode]",
            }
            return

        context = self._build_context(results)
        patient_section = self._format_patient_context(patient_context)
        conversation_section = self._format_conversation_history()

        prompt = (
            f"## Retrieved Evidence\n\n{context}\n\n"
            f"{patient_section}"
            f"{conversation_section}"
            f"---\n\n"
            f"## Question\n\n{question}\n\n"
            f"Please provide a comprehensive, evidence-based cardiovascular "
            f"medicine answer grounded in the retrieved evidence. "
            f"Cite sources using the format [Collection:record-id] or "
            f"[PMID:12345678](https://pubmed.ncbi.nlm.nih.gov/12345678/). "
            f"Follow ACC/AHA/ESC guideline citation format with Class/LOE."
        )

        full_answer = ""
        for token in self.llm.generate_stream(
            prompt=prompt,
            system_prompt=CARDIO_SYSTEM_PROMPT,
            max_tokens=2048,
            temperature=0.7,
        ):
            full_answer += token
            yield {"type": "token", "content": token}

        yield {"type": "done", "content": full_answer}

        # Update conversation history
        self.add_conversation_context("user", question)
        self.add_conversation_context("assistant", full_answer[:500])

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
                (e.g. 'modality == "echocardiography"').

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
        collections: List[str],
        query_vector: List[float],
        top_k: int,
        weights: Dict[str, float],
    ) -> List[CardioSearchResult]:
        """Search multiple collections in parallel with weighted scoring.

        Uses ThreadPoolExecutor for concurrent Milvus searches across
        all specified collections. Applies collection-specific weights
        to raw similarity scores for unified ranking.

        Args:
            collections: List of collection names to search.
            query_vector: 384-dimensional query embedding.
            top_k: Maximum results per collection.
            weights: Dict mapping collection name to weight multiplier.

        Returns:
            List of CardioSearchResult sorted by weighted score descending.
        """
        all_results: List[CardioSearchResult] = []
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
                                         "clinical_significance", "key_findings"):
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

                    result = CardioSearchResult(
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

        # Deduplicate by record_id (same record may appear via expansion)
        seen_ids: set = set()
        unique_results: List[CardioSearchResult] = []
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
        results: List[CardioSearchResult],
        query: str,
    ) -> List[CardioSearchResult]:
        """Rerank results based on relevance to original query.

        Applies heuristic boosts for:
        - Guideline results matching query conditions
        - Results from critical-flagged collections
        - Results with high citation relevance
        - PMID-bearing results (evidence-based)
        - Results matching detected biomarker or drug terms

        Args:
            results: Raw search results to rerank.
            query: Original query text for relevance matching.

        Returns:
            Reranked list of CardioSearchResult.
        """
        query_lower = query.lower()
        query_terms = set(query_lower.split())

        for result in results:
            boost = 0.0

            # Boost guideline results when query mentions conditions
            if result.collection == "cardio_guidelines":
                condition = result.metadata.get("condition", "").lower()
                if condition and condition in query_lower:
                    boost += 0.15
                # Boost Class I recommendations
                cor = result.metadata.get("class_of_rec", "")
                if cor == "I":
                    boost += 0.10
                elif cor == "IIa":
                    boost += 0.05

            # Boost literature with PMIDs
            pmid = result.metadata.get("pmid", "")
            if pmid:
                boost += 0.05

            # Boost results with high relevance
            if result.relevance == "high":
                boost += 0.10
            elif result.relevance == "medium":
                boost += 0.05

            # Boost heart failure results when GDMT terms present
            gdmt_terms = {"gdmt", "arni", "sglt2", "beta-blocker", "mra",
                          "sacubitril", "dapagliflozin", "empagliflozin",
                          "spironolactone", "eplerenone"}
            if result.collection == "cardio_heart_failure":
                if query_terms & gdmt_terms:
                    boost += 0.10

            # Boost imaging results when modality matches
            if result.collection == "cardio_imaging":
                modality = result.metadata.get("modality", "").lower()
                if modality and modality in query_lower:
                    boost += 0.10

            # Boost EP results for arrhythmia queries
            arrhythmia_terms = {"afib", "atrial fibrillation", "flutter",
                                "tachycardia", "bradycardia", "arrhythmia",
                                "ecg", "ekg", "rhythm"}
            if result.collection == "cardio_electrophysiology":
                if query_terms & arrhythmia_terms:
                    boost += 0.10

            # Apply boost
            result.score += boost

        # Re-sort after boosting
        results.sort(key=lambda r: r.score, reverse=True)
        return results

    # ══════════════════════════════════════════════════════════════════
    # LLM SYNTHESIS
    # ══════════════════════════════════════════════════════════════════

    def _synthesize(
        self,
        question: str,
        results: List[CardioSearchResult],
        workflow: Optional[CardioWorkflowType] = None,
        patient_context: Optional[dict] = None,
    ) -> CardioResponse:
        """Use LLM to synthesise search results into a clinical response.

        Builds a structured prompt with retrieved evidence, patient context,
        conversation history, and workflow-specific instructions. Generates
        a grounded clinical answer via the configured LLM.

        Args:
            question: Original clinical question.
            results: Ranked search results for context.
            workflow: Optional workflow for instruction tuning.
            patient_context: Optional patient-specific data dict.

        Returns:
            CardioResponse with synthesised answer and metadata.
        """
        context = self._build_context(results)
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
            f"Please provide a comprehensive, evidence-based cardiovascular medicine answer "
            f"grounded in the retrieved evidence above. Follow the system prompt instructions "
            f"for guideline citation format, severity badges, measurements with normal ranges, "
            f"and structured output sections.\n\n"
            f"Cite sources using clickable markdown links where PMIDs are available: "
            f"[PMID:12345678](https://pubmed.ncbi.nlm.nih.gov/12345678/). "
            f"For collection-sourced evidence, use [Collection:record-id]. "
            f"Prioritise [high relevance] citations."
        )

        answer = self.llm.generate(
            prompt=prompt,
            system_prompt=CARDIO_SYSTEM_PROMPT,
            max_tokens=2048,
            temperature=0.7,
        )

        return CardioResponse(
            question=question,
            answer=answer,
            results=results,
            workflow=workflow,
        )

    def _build_context(self, results: List[CardioSearchResult]) -> str:
        """Build context string from search results for LLM prompt.

        Organises results by collection, formatting each with its
        citation reference, relevance tag, score, and text excerpt.

        Args:
            results: Ranked search results to format.

        Returns:
            Formatted evidence context string for the LLM prompt.
        """
        if not results:
            return "No evidence found in the knowledge base."

        # Group results by collection
        by_collection: Dict[str, List[CardioSearchResult]] = {}
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

    def _format_citation_link(self, result: CardioSearchResult) -> str:
        """Format a citation with clickable URL where possible.

        Args:
            result: Search result to format citation for.

        Returns:
            Markdown-formatted citation string.
        """
        label = result.metadata.get("collection_label", result.collection)
        record_id = result.record_id

        # PubMed literature
        pmid = result.metadata.get("pmid", "")
        if pmid:
            return (
                f"[{label}:PMID {pmid}]"
                f"(https://pubmed.ncbi.nlm.nih.gov/{pmid}/)"
            )

        # Clinical trials
        nct_id = result.metadata.get("nct_id", "")
        if nct_id:
            return (
                f"[{label}:{nct_id}]"
                f"(https://clinicaltrials.gov/study/{nct_id})"
            )

        return f"[{label}:{record_id}]"

    def _format_patient_context(self, patient_context: Optional[dict]) -> str:
        """Format patient context for LLM prompt injection.

        Args:
            patient_context: Optional patient data dict with keys like
                age, sex, comorbidities, medications, labs, imaging.

        Returns:
            Formatted patient context section or empty string.
        """
        if not patient_context:
            return ""

        lines = ["### Patient Context\n"]

        field_labels = {
            "age": "Age",
            "sex": "Sex",
            "weight_kg": "Weight (kg)",
            "height_cm": "Height (cm)",
            "bmi": "BMI",
            "race": "Race/Ethnicity",
            "comorbidities": "Comorbidities",
            "medications": "Current Medications",
            "allergies": "Allergies",
            "labs": "Recent Labs",
            "vitals": "Vital Signs",
            "imaging": "Imaging Findings",
            "ecg": "ECG Findings",
            "lvef": "LVEF (%)",
            "nyha_class": "NYHA Class",
            "hf_stage": "HF Stage",
            "cha2ds2_vasc": "CHA2DS2-VASc Score",
            "has_bled": "HAS-BLED Score",
            "surgical_history": "Surgical History",
            "family_history": "Family History",
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
        if not self.conversation_history:
            return ""

        # Use only the most recent exchanges
        recent = self.conversation_history[-self._max_conversation_context * 2:]

        lines = ["### Conversation History\n"]
        for entry in recent:
            role = entry.get("role", "unknown").capitalize()
            content = entry.get("content", "")[:300]
            lines.append(f"**{role}:** {content}")

        lines.append("\n")
        return "\n".join(lines)

    def _format_workflow_instructions(
        self,
        workflow: Optional[CardioWorkflowType],
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
            CardioWorkflowType.CAD_ASSESSMENT: (
                "### Workflow: Coronary Artery Disease Assessment\n"
                "Focus on: pre-test probability, diagnostic pathway (anatomic vs functional), "
                "revascularization indications, optimal medical therapy, and secondary prevention. "
                "Reference CAD-RADS classification if CCTA results are available.\n\n"
            ),
            CardioWorkflowType.HEART_FAILURE: (
                "### Workflow: Heart Failure Evaluation\n"
                "Systematically address: HF phenotype (HFrEF/HFmrEF/HFpEF), ACC/AHA stage, "
                "NYHA class, all four GDMT pillars (ARNI, beta-blocker, MRA, SGLT2i) with "
                "titration status and target doses, device therapy eligibility (ICD, CRT), "
                "and advanced therapies if applicable.\n\n"
            ),
            CardioWorkflowType.VALVULAR_DISEASE: (
                "### Workflow: Valvular Heart Disease\n"
                "Address: severity grading with quantitative criteria, timing of intervention, "
                "surgical vs transcatheter approach, prosthetic valve selection, "
                "and long-term follow-up protocol.\n\n"
            ),
            CardioWorkflowType.ARRHYTHMIA: (
                "### Workflow: Arrhythmia / Electrophysiology\n"
                "Address: rhythm identification, hemodynamic stability assessment, "
                "rate vs rhythm control strategy, anticoagulation decision (CHA2DS2-VASc), "
                "ablation candidacy, and device therapy indications.\n\n"
            ),
            CardioWorkflowType.CARDIAC_MRI: (
                "### Workflow: Cardiac MRI Interpretation\n"
                "Address: tissue characterisation findings (T1/T2/LGE pattern), "
                "volumetric and functional parameters with normal ranges, "
                "differential diagnosis based on LGE distribution, "
                "and clinical implications.\n\n"
            ),
            CardioWorkflowType.STRESS_TEST: (
                "### Workflow: Stress Testing\n"
                "Address: test modality selection (exercise vs pharmacologic, "
                "with/without imaging), interpretation of results (Duke treadmill score, "
                "perfusion defects, wall motion abnormalities), and downstream management "
                "recommendations.\n\n"
            ),
            CardioWorkflowType.PREVENTIVE_RISK: (
                "### Workflow: Preventive Cardiology / Risk Assessment\n"
                "Address: ASCVD risk calculation, risk-enhancing factors, "
                "treatment thresholds, statin benefit groups, blood pressure targets, "
                "lifestyle modifications, and screening recommendations.\n\n"
            ),
            CardioWorkflowType.CARDIO_ONCOLOGY: (
                "### Workflow: Cardio-Oncology Consultation\n"
                "Address: baseline cardiovascular risk assessment, "
                "cancer-therapy cardiotoxicity risk stratification, "
                "monitoring protocol (biomarkers + imaging), "
                "cardioprotective strategies, and criteria for holding/modifying therapy.\n\n"
            ),
        }

        return instructions.get(workflow, "")

    # ══════════════════════════════════════════════════════════════════
    # CITATIONS & CONFIDENCE
    # ══════════════════════════════════════════════════════════════════

    def _extract_citations(
        self,
        results: List[CardioSearchResult],
    ) -> List[dict]:
        """Extract and format citations from search results.

        Generates a structured citation list from all results, including
        PMID links, trial references, and guideline citations.

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
            pmid = result.metadata.get("pmid", "")
            if pmid:
                cite["url"] = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                cite["id"] = f"PMID:{pmid}"

            nct_id = result.metadata.get("nct_id", "")
            if nct_id:
                cite["url"] = f"https://clinicaltrials.gov/study/{nct_id}"
                cite["id"] = nct_id

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
        results: List[CardioSearchResult],
    ) -> float:
        """Score overall confidence based on result quality.

        Confidence is based on:
        - Number of high-relevance results
        - Collection diversity
        - Average similarity score
        - Presence of guideline evidence

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

        # Factor 4: Guideline evidence present (0-0.2)
        has_guidelines = any(
            r.collection == "cardio_guidelines" for r in results
        )
        guideline_score = 0.2 if has_guidelines else 0.0

        confidence = relevance_score + diversity_score + quality_score + guideline_score
        return round(min(confidence, 1.0), 3)

    # ══════════════════════════════════════════════════════════════════
    # ENTITY & GUIDELINE SEARCH
    # ══════════════════════════════════════════════════════════════════

    def find_related(
        self,
        entity: str,
        entity_type: str = "condition",
    ) -> List[CardioSearchResult]:
        """Find related entities across collections.

        Searches all collections for evidence related to a clinical
        entity (condition, drug, gene, procedure). Useful for building
        entity profiles and cross-referencing.

        Args:
            entity: Entity name (e.g. 'heart failure', 'amiodarone', 'MYH7').
            entity_type: Entity category for targeted search:
                'condition', 'drug', 'gene', 'procedure', 'imaging'.

        Returns:
            List of CardioSearchResult from all relevant collections.
        """
        # Choose collections based on entity type
        type_collection_map = {
            "condition": [
                "cardio_guidelines", "cardio_literature", "cardio_trials",
                "cardio_heart_failure", "cardio_imaging",
            ],
            "drug": [
                "cardio_guidelines", "cardio_trials", "cardio_literature",
                "cardio_heart_failure", "cardio_prevention",
            ],
            "gene": [
                "genomic_evidence", "cardio_literature", "cardio_guidelines",
            ],
            "procedure": [
                "cardio_interventional", "cardio_guidelines", "cardio_trials",
            ],
            "imaging": [
                "cardio_imaging", "cardio_guidelines", "cardio_hemodynamics",
            ],
        }

        collections = type_collection_map.get(entity_type, get_all_collection_names())
        return self.search(entity, collections=collections, top_k=5)

    def get_guideline(
        self,
        condition: str,
        recommendation_class: Optional[str] = None,
    ) -> List[dict]:
        """Retrieve guideline recommendations for a condition.

        Searches the cardio_guidelines collection with optional filtering
        by class of recommendation.

        Args:
            condition: Clinical condition to search guidelines for.
            recommendation_class: Optional COR filter (e.g. 'I', 'IIa', 'IIb').

        Returns:
            List of guideline recommendation dicts with fields:
            society, title, year, recommendation, class_of_rec, evidence_level.
        """
        if not self.milvus:
            raise RuntimeError("Milvus client not configured.")

        query_vector = self._embed_query(condition)

        # Build filter expression
        filter_expr = None
        if recommendation_class:
            safe_class = recommendation_class.strip()
            if _SAFE_FILTER_RE.match(safe_class):
                filter_expr = f'class_of_rec == "{safe_class}"'
            else:
                logger.warning(
                    "Rejected unsafe filter value: %r", safe_class,
                )

        raw_results = self._search_collection(
            "cardio_guidelines", query_vector, top_k=10,
            filter_expr=filter_expr,
        )

        guidelines = []
        for record in raw_results:
            guideline = {
                "society": record.get("society", ""),
                "title": record.get("guideline_title", ""),
                "year": record.get("year", 0),
                "recommendation": record.get("recommendation", ""),
                "class_of_rec": record.get("class_of_rec", ""),
                "evidence_level": record.get("evidence_level", ""),
                "condition": record.get("condition", ""),
                "section": record.get("section", ""),
                "score": record.get("score", 0.0),
            }
            guidelines.append(guideline)

        return guidelines

    # ══════════════════════════════════════════════════════════════════
    # CONVERSATION MEMORY
    # ══════════════════════════════════════════════════════════════════

    def add_conversation_context(self, role: str, content: str, session_id: Optional[str] = None):
        """Add to conversation history for multi-turn context.

        Maintains a rolling window of recent conversation exchanges
        for follow-up query context injection.  Persists to disk so
        history survives restarts.

        Args:
            role: Message role ('user' or 'assistant').
            content: Message content text.
            session_id: Optional override; defaults to self.session_id.
        """
        self.conversation_history.append({
            "role": role,
            "content": content,
            "timestamp": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        })

        # Trim to max context window
        max_entries = self._max_conversation_context * 2
        if len(self.conversation_history) > max_entries:
            self.conversation_history = self.conversation_history[-max_entries:]

        # Persist to disk
        _save_conversation(session_id or self.session_id, self.conversation_history)

    def clear_conversation(self, session_id: Optional[str] = None):
        """Clear conversation history.

        Resets the multi-turn context and removes the persisted file.
        Useful when starting a new clinical consultation or switching
        topics.

        Args:
            session_id: Optional override; defaults to self.session_id.
        """
        self.conversation_history.clear()
        sid = session_id or self.session_id
        try:
            path = CONVERSATION_DIR / f"{sid}.json"
            if path.exists():
                path.unlink()
        except Exception as exc:
            logger.warning("Failed to remove conversation file %s: %s", sid, exc)

    # ══════════════════════════════════════════════════════════════════
    # HEALTH CHECK
    # ══════════════════════════════════════════════════════════════════

    def health_check(self) -> dict:
        """Check Milvus connection and collection status.

        Verifies connectivity to the Milvus server and checks that
        all expected cardiology collections exist and are loaded.

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
            # Check Milvus connectivity
            # Try listing collections to verify connection
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

            # Determine overall status
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

    # ══════════════════════════════════════════════════════════════════
    # WEIGHT COMPUTATION
    # ══════════════════════════════════════════════════════════════════

    def _get_boosted_weights(
        self,
        workflow: Optional[CardioWorkflowType] = None,
    ) -> Dict[str, float]:
        """Compute collection weights with optional workflow boosting.

        When a workflow is specified, applies boost multipliers from
        WORKFLOW_COLLECTION_BOOST on top of the base weights from
        settings. Weights are then renormalized to sum to ~1.0.

        Args:
            workflow: Optional CardioWorkflowType for boosting.

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


# =====================================================================
# FACTORY FUNCTION
# =====================================================================

def create_rag_engine(
    milvus_client=None,
    embedding_model=None,
    llm_client=None,
    query_expander=None,
) -> CardioRAGEngine:
    """Factory function to create a configured CardioRAGEngine.

    Convenience function that creates and returns a CardioRAGEngine
    with the provided components. If no components are provided,
    returns an engine in search-only mode (will raise on search
    without Milvus client).

    Args:
        milvus_client: Connected Milvus client. If None, attempts to
            create one from settings.
        embedding_model: Embedding model instance. If None, attempts to
            load BGE-small-en-v1.5 from settings.
        llm_client: LLM client instance. If None, LLM synthesis is
            disabled (search-only mode).
        query_expander: Optional QueryExpander instance.

    Returns:
        Configured CardioRAGEngine instance.
    """
    # Attempt auto-configuration if components not provided
    if milvus_client is None:
        try:
            from pymilvus import MilvusClient
            milvus_client = MilvusClient(
                uri=f"http://{settings.MILVUS_HOST}:{settings.MILVUS_PORT}",
            )
            logger.info(
                "Auto-connected to Milvus at %s:%d",
                settings.MILVUS_HOST, settings.MILVUS_PORT,
            )
        except Exception as exc:
            logger.warning(
                "Could not auto-connect to Milvus: %s. "
                "Search operations will fail.", exc,
            )

    if embedding_model is None:
        try:
            from sentence_transformers import SentenceTransformer

            class _EmbeddingWrapper:
                """Lightweight wrapper to provide embed_text() interface."""
                def __init__(self, model_name: str):
                    self.model = SentenceTransformer(model_name)

                def embed_text(self, text: str) -> List[float]:
                    return self.model.encode(text, normalize_embeddings=True).tolist()

            embedding_model = _EmbeddingWrapper(settings.EMBEDDING_MODEL)
            logger.info("Loaded embedding model: %s", settings.EMBEDDING_MODEL)
        except Exception as exc:
            logger.warning(
                "Could not load embedding model '%s': %s. "
                "Embedding operations will fail.",
                settings.EMBEDDING_MODEL, exc,
            )

    if llm_client is None and settings.ANTHROPIC_API_KEY:
        try:
            import anthropic

            class _LLMWrapper:
                """Lightweight wrapper to provide generate() and generate_stream() interface."""
                def __init__(self, api_key: str, model: str):
                    self.client = anthropic.Anthropic(api_key=api_key)
                    self.model = model

                def generate(
                    self,
                    prompt: str,
                    system_prompt: str = "",
                    max_tokens: int = 2048,
                    temperature: float = 0.7,
                ) -> str:
                    response = self.client.messages.create(
                        model=self.model,
                        max_tokens=max_tokens,
                        temperature=temperature,
                        system=system_prompt,
                        messages=[{"role": "user", "content": prompt}],
                    )
                    return response.content[0].text

                def generate_stream(
                    self,
                    prompt: str,
                    system_prompt: str = "",
                    max_tokens: int = 2048,
                    temperature: float = 0.7,
                ) -> Generator[str, None, None]:
                    with self.client.messages.stream(
                        model=self.model,
                        max_tokens=max_tokens,
                        temperature=temperature,
                        system=system_prompt,
                        messages=[{"role": "user", "content": prompt}],
                    ) as stream:
                        for text in stream.text_stream:
                            yield text

            llm_client = _LLMWrapper(settings.ANTHROPIC_API_KEY, settings.LLM_MODEL)
            logger.info("Configured LLM: %s", settings.LLM_MODEL)
        except Exception as exc:
            logger.warning(
                "Could not configure LLM client: %s. "
                "LLM synthesis disabled (search-only mode).", exc,
            )

    if query_expander is None:
        try:
            from .query_expansion import QueryExpander
            query_expander = QueryExpander()
            logger.info("Loaded query expander")
        except Exception as exc:
            logger.debug("Query expander not available: %s", exc)

    engine = CardioRAGEngine(
        milvus_client=milvus_client,
        embedding_model=embedding_model,
        llm_client=llm_client,
        query_expander=query_expander,
    )

    logger.info(
        "CardioRAGEngine created: milvus=%s, embedder=%s, llm=%s, expander=%s",
        milvus_client is not None,
        embedding_model is not None,
        llm_client is not None,
        query_expander is not None,
    )

    return engine


# =====================================================================
# DISPLAY HELPER
# =====================================================================

def format_search_results(results: List[CardioSearchResult]) -> str:
    """Format search results for human-readable display.

    Produces a structured text output suitable for terminal or log
    display, with collection grouping, relevance tags, and score.

    Args:
        results: List of CardioSearchResult to format.

    Returns:
        Formatted multi-line string.
    """
    if not results:
        return "No results found."

    lines: List[str] = [
        f"Found {len(results)} results across "
        f"{len(set(r.collection for r in results))} collections:\n",
    ]

    # Group by collection
    by_collection: Dict[str, List[CardioSearchResult]] = {}
    for result in results:
        label = result.metadata.get("collection_label", result.collection)
        if label not in by_collection:
            by_collection[label] = []
        by_collection[label].append(result)

    for label, coll_results in by_collection.items():
        lines.append(f"--- {label} ({len(coll_results)} results) ---")
        for i, result in enumerate(coll_results[:5], 1):
            relevance_tag = f"[{result.relevance}]" if result.relevance else ""
            text_preview = result.text[:120].replace("\n", " ") if result.text else "(no text)"
            lines.append(
                f"  {i}. {relevance_tag} (score={result.score:.3f}) "
                f"{text_preview}"
            )

            # Show PMID or NCT if available
            pmid = result.metadata.get("pmid", "")
            nct = result.metadata.get("nct_id", "")
            if pmid:
                lines.append(f"     PMID: {pmid}")
            if nct:
                lines.append(f"     NCT: {nct}")

        lines.append("")

    return "\n".join(lines)
