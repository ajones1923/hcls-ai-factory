"""Multi-collection RAG engine for Neurology Intelligence Agent.

Searches across all 14 neurology-specific Milvus collections simultaneously
using parallel ThreadPoolExecutor, synthesises findings with neurological
knowledge augmentation, and generates grounded LLM responses with clinical
evidence citations.

Extends the pattern from: rag-chat-pipeline/src/rag_engine.py

Features:
- Parallel search via ThreadPoolExecutor (13 neuro + 1 shared genomic collection)
- Settings-driven weights and parameters from config/settings.py
- Workflow-based dynamic weight boosting per NeuroWorkflowType
- Milvus field-based filtering (modality, condition, severity, region)
- Citation relevance scoring (high/medium/low) with PMID/DOI link formatting
- Cross-collection entity linking for comprehensive neurological queries
- Guideline retrieval with AAN/EAN/ILAE/MDS document identifiers
- Conversation memory for multi-turn clinical consultations
- Patient context injection for personalised clinical decision support
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
    NEURO_SYSTEM_PROMPT,
    WORKFLOW_COLLECTION_BOOST,
    NeuroWorkflowType,
    NeuroResponse,
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
class NeuroSearchResult:
    """A single search result from a Milvus collection.

    Attributes:
        collection: Source collection name (e.g. 'neuro_cerebrovascular').
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
    "neuro_literature": {
        "weight": settings.WEIGHT_LITERATURE,
        "label": "Literature",
        "text_field": "abstract",
        "title_field": "title",
        "filterable_fields": ["study_type", "condition", "journal"],
    },
    "neuro_trials": {
        "weight": settings.WEIGHT_TRIALS,
        "label": "ClinicalTrial",
        "text_field": "trial_summary",
        "title_field": "trial_title",
        "filterable_fields": ["phase", "status", "condition", "intervention"],
    },
    "neuro_imaging": {
        "weight": settings.WEIGHT_IMAGING,
        "label": "Imaging",
        "text_field": "imaging_findings",
        "title_field": "modality",
        "filterable_fields": ["modality", "condition", "brain_region", "scoring_system"],
    },
    "neuro_electrophysiology": {
        "weight": settings.WEIGHT_ELECTROPHYSIOLOGY,
        "label": "Electrophysiology",
        "text_field": "findings_summary",
        "title_field": "study_type",
        "filterable_fields": ["study_type", "pattern", "condition", "localization"],
    },
    "neuro_degenerative": {
        "weight": settings.WEIGHT_DEGENERATIVE,
        "label": "Neurodegeneration",
        "text_field": "clinical_description",
        "title_field": "condition",
        "filterable_fields": ["condition", "stage", "biomarker", "gene"],
    },
    "neuro_cerebrovascular": {
        "weight": settings.WEIGHT_CEREBROVASCULAR,
        "label": "Cerebrovascular",
        "text_field": "clinical_description",
        "title_field": "condition",
        "filterable_fields": ["stroke_type", "territory", "etiology", "severity"],
    },
    "neuro_epilepsy": {
        "weight": settings.WEIGHT_EPILEPSY,
        "label": "Epilepsy",
        "text_field": "seizure_description",
        "title_field": "epilepsy_type",
        "filterable_fields": ["classification", "etiology", "eeg_pattern", "syndrome"],
    },
    "neuro_oncology": {
        "weight": settings.WEIGHT_ONCOLOGY,
        "label": "NeuroOncology",
        "text_field": "tumor_description",
        "title_field": "tumor_type",
        "filterable_fields": ["who_grade", "molecular_markers", "location", "treatment"],
    },
    "neuro_ms": {
        "weight": settings.WEIGHT_MS,
        "label": "MS_Neuroimmunology",
        "text_field": "clinical_description",
        "title_field": "condition",
        "filterable_fields": ["ms_type", "activity_status", "dmt", "antibody_status"],
    },
    "neuro_movement": {
        "weight": settings.WEIGHT_MOVEMENT,
        "label": "MovementDisorder",
        "text_field": "clinical_description",
        "title_field": "condition",
        "filterable_fields": ["disorder_type", "phenotype", "gene", "treatment"],
    },
    "neuro_headache": {
        "weight": settings.WEIGHT_HEADACHE,
        "label": "Headache",
        "text_field": "headache_description",
        "title_field": "headache_type",
        "filterable_fields": ["ichd3_code", "frequency", "treatment_type", "severity"],
    },
    "neuro_neuromuscular": {
        "weight": settings.WEIGHT_NEUROMUSCULAR,
        "label": "Neuromuscular",
        "text_field": "clinical_description",
        "title_field": "condition",
        "filterable_fields": ["disorder_type", "antibody", "gene", "emg_pattern"],
    },
    "neuro_guidelines": {
        "weight": settings.WEIGHT_GUIDELINES,
        "label": "Guideline",
        "text_field": "recommendation",
        "title_field": "guideline_title",
        "filterable_fields": ["society", "evidence_class", "recommendation_level",
                              "condition"],
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
# NEURO RAG ENGINE
# =====================================================================

class NeuroRAGEngine:
    """Multi-collection RAG engine for neurology intelligence.

    Searches across all 14 neurology-specific Milvus collections plus the
    shared genomic_evidence collection. Supports workflow-specific weight
    boosting, parallel search, query expansion, patient context injection,
    and multi-turn conversation memory.

    Features:
    - Parallel search via ThreadPoolExecutor (14 collections)
    - Settings-driven weights and parameters
    - Workflow-based dynamic weight boosting (9 neuro workflows)
    - Milvus field-based filtering (modality, condition, severity, region)
    - Citation relevance scoring (high/medium/low)
    - Cross-collection entity linking
    - Guideline retrieval with AAN/EAN/ILAE/MDS document IDs
    - Conversation memory context injection
    - Patient context for personalised clinical decision support
    - Confidence scoring based on evidence diversity

    Usage:
        engine = NeuroRAGEngine(milvus_client, embedding_model, llm_client)
        response = engine.query("68yo acute stroke, NIHSS 14, 2 hours onset")
        results = engine.search("Dravet syndrome SCN1A fenfluramine")
    """

    def __init__(
        self,
        milvus_client=None,
        embedding_model=None,
        llm_client=None,
        session_id: str = "default",
    ):
        """Initialize the NeuroRAGEngine.

        Args:
            milvus_client: Connected Milvus client with access to all
                neuro collections. If None, search operations will
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

    # ==================================================================
    # PROPERTIES
    # ==================================================================

    @property
    def conversation_history(self) -> List[Dict[str, str]]:
        """Return current conversation history."""
        return self._conversation_history

    # ==================================================================
    # PUBLIC API
    # ==================================================================

    def query(
        self,
        question: str,
        workflow: Optional[NeuroWorkflowType] = None,
        top_k: int = 5,
        patient_context: Optional[dict] = None,
    ) -> NeuroResponse:
        """Main query method: expand -> search -> synthesise.

        Performs the full RAG pipeline: parallel multi-collection search
        with workflow-specific weighting, result reranking, LLM synthesis
        with patient context, and confidence scoring.

        Args:
            question: Natural language neurology question.
            workflow: Optional NeuroWorkflowType to apply domain-specific
                collection weight boosting. If None, weights are auto-detected
                or base defaults are used.
            top_k: Maximum results to return per collection.
            patient_context: Optional dict with patient-specific data
                (age, sex, symptoms, imaging, labs, medications, history)
                for personalised clinical decision support.

        Returns:
            NeuroResponse with synthesised answer, search results, citations,
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
            response = NeuroResponse(
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
    ) -> List[NeuroSearchResult]:
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
            List of NeuroSearchResult sorted by weighted score descending.
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

    # ==================================================================
    # EMBEDDING
    # ==================================================================

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

    # ==================================================================
    # COLLECTION SEARCH
    # ==================================================================

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
                (e.g. 'condition == "ischemic stroke"').

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
    ) -> List[NeuroSearchResult]:
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
            List of NeuroSearchResult sorted by weighted score descending.
        """
        all_results: List[NeuroSearchResult] = []
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
                                         "clinical_description", "findings_summary",
                                         "imaging_findings", "seizure_description",
                                         "tumor_description", "headache_description",
                                         "trial_summary", "text_chunk"):
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

                    result = NeuroSearchResult(
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
        unique_results: List[NeuroSearchResult] = []
        for result in all_results:
            dedup_key = f"{result.collection}:{result.record_id}"
            if dedup_key not in seen_ids:
                seen_ids.add(dedup_key)
                unique_results.append(result)

        # Cap at reasonable limit
        return unique_results[:top_k * len(collections)]

    # ==================================================================
    # RERANKING
    # ==================================================================

    def _rerank_results(
        self,
        results: List[NeuroSearchResult],
        query: str,
    ) -> List[NeuroSearchResult]:
        """Rerank results based on relevance to original query.

        Applies heuristic boosts for:
        - Guideline results matching query societies (AAN, EAN, ILAE, MDS)
        - Results from imaging/electrophysiology collections for diagnostic queries
        - Results with high citation relevance
        - PMID/DOI-bearing results (evidence-based)
        - Results matching detected biomarker or drug terms

        Args:
            results: Raw search results to rerank.
            query: Original query text for relevance matching.

        Returns:
            Reranked list of NeuroSearchResult.
        """
        query_lower = query.lower()
        query_terms = set(query_lower.split())

        for result in results:
            boost = 0.0

            # Boost guideline results when query mentions guideline bodies
            if result.collection == "neuro_guidelines":
                society = result.metadata.get("society", "").lower()
                if society and society in query_lower:
                    boost += 0.15
                # Boost AAN/EAN/ILAE/MDS guidelines
                for body in ("aan", "ean", "ilae", "mds", "ihs"):
                    if body in query_lower:
                        boost += 0.08
                        break

            # Boost guideline results generally
            if result.collection == "neuro_guidelines":
                boost += 0.05

            # Boost results with PMIDs
            pmid = result.metadata.get("pmid", "")
            if pmid:
                boost += 0.05

            # Boost results with DOIs
            doi = result.metadata.get("doi", "")
            if doi:
                boost += 0.03

            # Boost results with NCT IDs (clinical trial evidence)
            nct_id = result.metadata.get("nct_id", "")
            if nct_id:
                boost += 0.05

            # Boost results with high relevance
            if result.relevance == "high":
                boost += 0.10
            elif result.relevance == "medium":
                boost += 0.05

            # Boost imaging results for imaging-related queries
            imaging_terms = {"mri", "ct", "pet", "flair", "dwi", "swi",
                             "angiography", "dat-spect", "datscan", "imaging",
                             "scan", "aspects", "mra", "cta", "atrophy"}
            if result.collection == "neuro_imaging":
                if query_terms & imaging_terms:
                    boost += 0.10

            # Boost electrophysiology results for EEG/EMG queries
            ephys_terms = {"eeg", "emg", "ncs", "nerve conduction",
                           "spike-wave", "epileptiform", "electroencephalography",
                           "electromyography", "evoked potential"}
            if result.collection == "neuro_electrophysiology":
                if query_terms & ephys_terms:
                    boost += 0.10

            # Boost cerebrovascular results for stroke queries
            stroke_terms = {"stroke", "tpa", "thrombectomy", "hemorrhage",
                            "nihss", "ischemic", "lvo", "sah", "tia",
                            "infarction", "hematoma", "carotid"}
            if result.collection == "neuro_cerebrovascular":
                if query_terms & stroke_terms:
                    boost += 0.10

            # Boost degenerative results for dementia queries
            dementia_terms = {"dementia", "alzheimer", "cognitive", "mci",
                              "moca", "amyloid", "tau", "frontotemporal",
                              "lewy", "neurodegeneration", "memory"}
            if result.collection == "neuro_degenerative":
                if query_terms & dementia_terms:
                    boost += 0.10

            # Boost epilepsy results for seizure queries
            epilepsy_terms = {"epilepsy", "seizure", "convulsion", "status",
                              "antiseizure", "asm", "aed", "ilae",
                              "spike", "focal", "generalized", "absence"}
            if result.collection == "neuro_epilepsy":
                if query_terms & epilepsy_terms:
                    boost += 0.10

            # Boost MS results for MS/neuroimmunology queries
            ms_terms = {"multiple sclerosis", "ms", "rrms", "ppms", "spms",
                        "nmosd", "mogad", "oligoclonal", "mcdonald",
                        "dmt", "relapse", "demyelination"}
            if result.collection == "neuro_ms":
                if query_terms & ms_terms:
                    boost += 0.10

            # Boost neuromuscular results for NMJ/motor neuron queries
            nm_terms = {"myasthenia", "als", "motor neuron", "neuropathy",
                        "gbs", "cidp", "weakness", "fasciculation",
                        "sma", "muscular dystrophy", "emg"}
            if result.collection == "neuro_neuromuscular":
                if query_terms & nm_terms:
                    boost += 0.10

            # Boost genomic results for genetic queries
            genetic_terms = {"gene", "mutation", "variant", "genetic",
                             "apoe", "gba", "lrrk2", "scn1a", "sod1",
                             "c9orf72", "htt", "smn1", "smn2"}
            if result.collection == "genomic_evidence":
                if query_terms & genetic_terms:
                    boost += 0.10

            # Apply boost
            result.score += boost

        # Re-sort after boosting
        results.sort(key=lambda r: r.score, reverse=True)
        return results

    # ==================================================================
    # CITATION SCORING
    # ==================================================================

    def _score_citations(
        self,
        results: List[NeuroSearchResult],
    ) -> List[NeuroSearchResult]:
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

    # ==================================================================
    # LLM SYNTHESIS
    # ==================================================================

    def _synthesize_response(
        self,
        question: str,
        results: List[NeuroSearchResult],
        workflow: Optional[NeuroWorkflowType] = None,
        patient_context: Optional[dict] = None,
    ) -> NeuroResponse:
        """Use LLM to synthesise search results into a neurology response.

        Builds a structured prompt with retrieved evidence, patient context,
        conversation history, and workflow-specific instructions. Generates
        a grounded answer via the configured LLM.

        Args:
            question: Original neurology question.
            results: Ranked search results for context.
            workflow: Optional workflow for instruction tuning.
            patient_context: Optional patient-specific data dict.

        Returns:
            NeuroResponse with synthesised answer and metadata.
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
            f"Please provide a comprehensive, evidence-based neurology clinical "
            f"decision support answer grounded in the retrieved evidence above. "
            f"Follow the system prompt instructions for clinical citation format, "
            f"severity badges, clinical scale scoring, and structured output sections.\n\n"
            f"Cite sources using clickable markdown links where PMIDs are available: "
            f"[PMID:12345678](https://pubmed.ncbi.nlm.nih.gov/12345678/). "
            f"For clinical trials, use [NCT01234567](https://clinicaltrials.gov/study/NCT01234567). "
            f"For collection-sourced evidence, use [Collection:record-id]. "
            f"Prioritise [high relevance] citations and guideline references."
        )

        answer = self.llm.generate(
            prompt=prompt,
            system_prompt=NEURO_SYSTEM_PROMPT,
            max_tokens=2048,
            temperature=0.7,
        )

        return NeuroResponse(
            question=question,
            answer=answer,
            results=results,
            workflow=workflow,
        )

    def _build_context(
        self,
        results: List[NeuroSearchResult],
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
        by_collection: Dict[str, List[NeuroSearchResult]] = {}
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

    def _format_citation_link(self, result: NeuroSearchResult) -> str:
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

        # ClinicalTrials.gov
        nct_id = result.metadata.get("nct_id", "")
        if nct_id:
            return (
                f"[{label}:{nct_id}]"
                f"(https://clinicaltrials.gov/study/{nct_id})"
            )

        # DOI
        doi = result.metadata.get("doi", "")
        if doi:
            return f"[{label}:DOI {doi}](https://doi.org/{doi})"

        # AAN/EAN guideline reference
        guideline_id = result.metadata.get("guideline_id", "")
        if guideline_id:
            return f"[{label}:{guideline_id}]"

        return f"[{label}:{record_id}]"

    def _format_patient_context(self, patient_context: Optional[dict]) -> str:
        """Format patient context for LLM prompt injection.

        Used for patient-specific clinical decision support.

        Args:
            patient_context: Optional patient data dict with keys like
                age, sex, symptoms, onset_time, nihss, gcs, moca,
                imaging_findings, eeg_findings, csf_results, medications,
                prior_therapies, comorbidities, family_history, genomic_data.

        Returns:
            Formatted patient context section or empty string.
        """
        if not patient_context:
            return ""

        lines = ["### Patient Context\n"]

        field_labels = {
            "age": "Age",
            "sex": "Sex",
            "chief_complaint": "Chief Complaint",
            "onset_time": "Symptom Onset",
            "onset_duration": "Duration",
            "symptoms": "Presenting Symptoms",
            "neurological_exam": "Neurological Examination",
            "nihss": "NIHSS Score",
            "gcs": "GCS Score",
            "moca": "MoCA Score",
            "updrs": "UPDRS Part III Score",
            "edss": "EDSS Score",
            "alsfrs_r": "ALSFRS-R Score",
            "mrs": "Modified Rankin Scale",
            "hit6": "HIT-6 Score",
            "imaging_findings": "Imaging Findings",
            "eeg_findings": "EEG Findings",
            "emg_findings": "EMG/NCS Findings",
            "csf_results": "CSF Results",
            "biomarkers": "Biomarkers",
            "genomic_data": "Genomic Data",
            "medications": "Current Medications",
            "prior_therapies": "Prior Treatments",
            "allergies": "Allergies",
            "comorbidities": "Comorbidities",
            "family_history": "Family History",
            "surgical_history": "Surgical History",
            "functional_status": "Functional Status",
            "occupation": "Occupation",
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
        workflow: Optional[NeuroWorkflowType],
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
            NeuroWorkflowType.STROKE_ACUTE: (
                "### Workflow: Acute Stroke Evaluation\n"
                "Focus on: time from onset/LKW, tPA eligibility (< 4.5h, NIHSS, "
                "contraindications), thrombectomy eligibility (LVO, ASPECTS >= 6, "
                "favorable perfusion), blood pressure management, imaging interpretation "
                "(CT/CTA/CTP or MRI DWI/FLAIR/MRA), stroke etiology (TOAST classification), "
                "and acute complications (hemorrhagic transformation, malignant edema, "
                "seizures). Reference AHA/ASA 2019 acute stroke guidelines and "
                "ESO thrombolysis/thrombectomy guidelines.\n\n"
            ),
            NeuroWorkflowType.STROKE_PREVENTION: (
                "### Workflow: Stroke Prevention\n"
                "Focus on: stroke etiology and TOAST classification, antiplatelet vs "
                "anticoagulation selection, dual antiplatelet duration (POINT/CHANCE trials), "
                "carotid intervention criteria (NASCET/ECST), PFO closure criteria "
                "(RESPECT/CLOSE/DEFENSE PFO), intracranial stenosis management "
                "(SAMMPRIS/CASSISS), and vascular risk factor optimization (blood pressure, "
                "lipids, diabetes, lifestyle). Reference AHA/ASA secondary prevention "
                "guidelines.\n\n"
            ),
            NeuroWorkflowType.DEMENTIA_EVALUATION: (
                "### Workflow: Dementia Evaluation\n"
                "Focus on: cognitive screening (MoCA, MMSE), neuropsychological testing "
                "domains, NIA-AA ATN biomarker framework (amyloid/tau/neurodegeneration), "
                "CSF biomarkers (Abeta42, p-tau181, p-tau217, NfL, t-tau), blood-based "
                "biomarkers (p-tau217, GFAP, NfL), amyloid and tau PET interpretation, "
                "structural MRI (hippocampal atrophy, Scheltens visual rating, Fazekas), "
                "differential diagnosis (AD vs FTD vs DLB vs VaD vs NPH), reversible "
                "causes (B12, TSH, RPR, HIV), and treatment options (cholinesterase "
                "inhibitors, memantine, anti-amyloid antibodies). Reference NIA-AA 2024 "
                "clinical criteria and DLB Consortium fourth report.\n\n"
            ),
            NeuroWorkflowType.EPILEPSY_CLASSIFICATION: (
                "### Workflow: Epilepsy Classification & Management\n"
                "Focus on: ILAE 2017 seizure classification (focal vs generalized vs "
                "unknown onset, aware vs impaired awareness, motor vs non-motor), "
                "epilepsy syndrome identification, EEG interpretation (interictal "
                "epileptiform discharges, ictal patterns, background), MRI lesion "
                "correlation, first-line ASM selection by seizure type (focal: "
                "lamotrigine/levetiracetam, generalized: valproate/lamotrigine), "
                "drug-resistant epilepsy evaluation (failed 2 appropriate ASMs), "
                "presurgical workup (video-EEG, MRI 3T epilepsy protocol, PET, SEEG), "
                "and non-pharmacological options (VNS, RNS, DBS, dietary therapy). "
                "Reference ILAE 2017 classification and AAN/AES treatment guidelines.\n\n"
            ),
            NeuroWorkflowType.MS_MANAGEMENT: (
                "### Workflow: MS & Neuroimmunology Management\n"
                "Focus on: McDonald 2017 diagnostic criteria (DIS/DIT, CSF OCB "
                "substitution), MRI lesion characterization (periventricular, "
                "juxtacortical, infratentorial, spinal cord; Dawson fingers), "
                "disease activity assessment (relapses, EDSS progression, MRI "
                "activity, NfL), DMT escalation algorithm (platform therapies vs "
                "high-efficacy therapies), JCV risk stratification for natalizumab, "
                "NMOSD vs MOGAD differentiation (AQP4-IgG, MOG-IgG), acute relapse "
                "management (IV methylprednisolone, PLEX), and pregnancy planning. "
                "Reference AAN 2018 DMT guidelines and ECTRIMS/EAN 2024 treatment "
                "guidelines.\n\n"
            ),
            NeuroWorkflowType.MOVEMENT_DISORDER: (
                "### Workflow: Movement Disorder Evaluation\n"
                "Focus on: MDS clinical diagnostic criteria for PD (bradykinesia "
                "plus rigidity or tremor, supportive features, red flags, exclusion "
                "criteria), atypical parkinsonism differentiation (PSP: vertical gaze "
                "palsy, early falls; MSA: autonomic failure, cerebellar ataxia; CBD: "
                "alien limb, apraxia, cortical sensory loss), DaT-SPECT interpretation, "
                "PD treatment algorithm (levodopa, dopamine agonists, MAO-B inhibitors, "
                "COMT inhibitors, DBS), essential tremor vs PD tremor, dystonia "
                "classification and botulinum toxin injection, and genetic parkinsonism "
                "(GBA, LRRK2). Reference MDS diagnostic criteria and AAN PD treatment "
                "guidelines.\n\n"
            ),
            NeuroWorkflowType.HEADACHE_DIAGNOSIS: (
                "### Workflow: Headache Diagnosis & Management\n"
                "Focus on: ICHD-3 diagnostic criteria for migraine (with/without aura), "
                "cluster headache, and tension-type headache, headache red flags "
                "(SNNOOP10 mnemonic: thunderclap, papilledema, fever, focal deficits, "
                "age > 50, cancer history, immunosuppression, positional, progressive), "
                "acute treatment (triptans, gepants, ditans, NSAIDs), preventive "
                "therapy (CGRP monoclonal antibodies, beta-blockers, topiramate, "
                "amitriptyline), medication overuse headache management, and IIH "
                "evaluation (opening pressure, acetazolamide, visual monitoring). "
                "Reference IHS ICHD-3 classification and AAN migraine treatment "
                "guidelines.\n\n"
            ),
            NeuroWorkflowType.NEUROMUSCULAR_EVAL: (
                "### Workflow: Neuromuscular Evaluation\n"
                "Focus on: electrodiagnostic interpretation (EMG: spontaneous activity, "
                "MUP morphology, recruitment; NCS: demyelinating vs axonal patterns, "
                "conduction block, temporal dispersion), MG diagnostics (AChR/MuSK "
                "antibodies, RNS decrement, single-fiber EMG jitter), ALS diagnostic "
                "criteria (El Escorial/Awaji: UMN + LMN signs in multiple regions), "
                "GBS evaluation (areflexia, ascending weakness, FVC monitoring for "
                "respiratory failure, CSF albuminocytologic dissociation), CIDP "
                "criteria (EFNS/PNS), genetic neuromuscular disease testing (SMN1, "
                "dystrophin, SOD1, C9orf72), and treatment selection (IVIG, PLEX, "
                "FcRn inhibitors, complement inhibitors, gene therapy). Reference "
                "AAN evidence-based guidelines for each condition.\n\n"
            ),
            NeuroWorkflowType.NEURO_ONCOLOGY: (
                "### Workflow: Neuro-Oncology Evaluation\n"
                "Focus on: WHO 2021 CNS tumor classification (integrated diagnosis "
                "with molecular markers), glioma molecular stratification (IDH status, "
                "1p/19q codeletion, MGMT methylation, TERT promoter, H3K27), RANO "
                "criteria for treatment response assessment, standard of care "
                "(Stupp protocol for GBM: maximal safe resection + concurrent "
                "temozolomide/radiation + adjuvant temozolomide), TTFields, "
                "brain metastasis management (SRS vs WBRT, immunotherapy "
                "considerations, GPA scoring), and clinical trial options. "
                "Reference NCCN CNS cancer guidelines and EANO guidelines.\n\n"
            ),
        }

        return instructions.get(workflow, "")

    # ==================================================================
    # CITATIONS & CONFIDENCE
    # ==================================================================

    def _extract_citations(
        self,
        results: List[NeuroSearchResult],
    ) -> List[dict]:
        """Extract and format citations from search results.

        Generates a structured citation list from all results, including
        PMID links, NCT links, DOI links, and guideline references.

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

            # Guideline identifiers
            guideline_id = result.metadata.get("guideline_id", "")
            if guideline_id and not cite["url"]:
                cite["id"] = guideline_id

            # Deduplicate
            dedup_key = cite["id"] or f"{cite['source']}:{result.record_id}"
            if dedup_key not in seen:
                seen.add(dedup_key)
                citations.append(cite)

        return citations

    def _score_confidence(
        self,
        results: List[NeuroSearchResult],
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
        diversity_score = min(unique_collections / 4, 1.0) * 0.3

        # Factor 3: Average score of top results (0-0.2)
        top_scores = [r.score for r in results[:5]]
        avg_score = sum(top_scores) / max(len(top_scores), 1)
        quality_score = min(avg_score, 1.0) * 0.2

        # Factor 4: Guideline evidence present (0-0.2)
        has_guidelines = any(
            r.collection == "neuro_guidelines" for r in results
        )
        guideline_score = 0.2 if has_guidelines else 0.0

        confidence = relevance_score + diversity_score + quality_score + guideline_score
        return round(min(confidence, 1.0), 3)

    # ==================================================================
    # ENTITY & CONDITION SEARCH
    # ==================================================================

    def find_related(
        self,
        entity: str,
        entity_type: str = "condition",
        top_k: int = 5,
    ) -> List[NeuroSearchResult]:
        """Find related entities across collections.

        Searches all collections for evidence related to a neurological
        entity (condition, drug, biomarker, imaging modality). Useful
        for building entity profiles and cross-referencing.

        Args:
            entity: Entity name (e.g. 'Parkinson disease', 'lecanemab', 'NfL').
            entity_type: Entity category for targeted search:
                'condition', 'drug', 'biomarker', 'imaging', 'gene'.
            top_k: Maximum results per collection.

        Returns:
            List of NeuroSearchResult from all relevant collections.
        """
        type_collection_map = {
            "condition": [
                "neuro_literature", "neuro_trials", "neuro_guidelines",
                "neuro_degenerative", "neuro_cerebrovascular",
                "neuro_epilepsy", "neuro_ms", "neuro_movement",
                "neuro_headache", "neuro_neuromuscular",
            ],
            "drug": [
                "neuro_trials", "neuro_literature", "neuro_guidelines",
                "neuro_cerebrovascular", "neuro_degenerative",
                "neuro_epilepsy", "neuro_ms",
            ],
            "biomarker": [
                "neuro_degenerative", "neuro_imaging", "neuro_electrophysiology",
                "neuro_literature", "neuro_guidelines", "neuro_ms",
                "neuro_neuromuscular",
            ],
            "imaging": [
                "neuro_imaging", "neuro_literature", "neuro_guidelines",
                "neuro_cerebrovascular", "neuro_degenerative",
                "neuro_oncology",
            ],
            "gene": [
                "genomic_evidence", "neuro_degenerative", "neuro_movement",
                "neuro_neuromuscular", "neuro_epilepsy", "neuro_literature",
            ],
        }

        collections = type_collection_map.get(entity_type, get_all_collection_names())
        return self.search(entity, collections=collections, top_k=top_k)

    def get_guideline(
        self,
        condition: str,
        society: Optional[str] = None,
    ) -> List[NeuroSearchResult]:
        """Retrieve clinical guidelines for a condition.

        Searches the neuro_guidelines collection with optional society
        filtering (AAN, EAN, ILAE, MDS, IHS).

        Args:
            condition: Neurological condition to search guidelines for.
            society: Optional guideline body filter (e.g. 'AAN', 'ILAE').

        Returns:
            List of NeuroSearchResult from guidelines collection.
        """
        query = f"clinical guideline recommendation {condition}"
        if society:
            query = f"{society} {query}"

        return self.search(query, collections=["neuro_guidelines"], top_k=10)

    def search_imaging(
        self,
        modality: str,
        findings: str,
        top_k: int = 10,
    ) -> List[NeuroSearchResult]:
        """Search imaging-specific evidence.

        Targeted search combining modality and findings for neuroimaging
        interpretation support.

        Args:
            modality: Imaging modality (e.g. 'MRI brain', 'CT angiography').
            findings: Imaging findings description.
            top_k: Maximum results.

        Returns:
            List of NeuroSearchResult from imaging and related collections.
        """
        query = f"{modality} {findings}"
        collections = ["neuro_imaging", "neuro_literature", "neuro_guidelines"]
        return self.search(query, collections=collections, top_k=top_k)

    # ==================================================================
    # CONVERSATION MEMORY
    # ==================================================================

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

    # ==================================================================
    # WEIGHT COMPUTATION
    # ==================================================================

    def _get_boosted_weights(
        self,
        workflow: Optional[NeuroWorkflowType] = None,
    ) -> Dict[str, float]:
        """Compute collection weights with optional workflow boosting.

        When a workflow is specified, applies boost multipliers from
        WORKFLOW_COLLECTION_BOOST on top of the base weights from
        settings. Weights are then renormalized to sum to ~1.0.

        Args:
            workflow: Optional NeuroWorkflowType for boosting.

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

    # ==================================================================
    # HEALTH CHECK
    # ==================================================================

    def health_check(self) -> dict:
        """Check Milvus connection and collection status.

        Verifies connectivity to the Milvus server and checks that
        all expected neuro collections exist and are loaded.

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
