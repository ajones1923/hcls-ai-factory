"""Multi-collection RAG engine for Rare Disease Diagnostic Agent.

Searches across all 14 rare-disease-specific Milvus collections simultaneously
using parallel ThreadPoolExecutor, synthesises findings with clinical genetics
knowledge augmentation, and generates grounded LLM responses with OMIM, HPO,
and ACMG citations.

Extends the pattern from: rag-chat-pipeline/src/rag_engine.py

Features:
- Parallel search via ThreadPoolExecutor (13 rare disease + 1 shared genomic collection)
- Settings-driven weights and parameters from config/settings.py
- Workflow-based dynamic weight boosting per DiagnosticWorkflowType
- Milvus field-based filtering (gene, inheritance, phenotype, OMIM)
- Citation relevance scoring (high/medium/low) with OMIM/HPO link formatting
- Cross-collection entity linking for comprehensive diagnostic queries
- Phenotype-to-disease association retrieval
- Conversation memory for multi-turn diagnostic consultations
- Patient context injection for personalised diagnostic assessment
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

from .agent import (  # noqa: F401 -- re-exported for convenience
    RARE_DISEASE_SYSTEM_PROMPT,
    WORKFLOW_COLLECTION_BOOST,
    RARE_DISEASE_CONDITIONS,
    RARE_DISEASE_GENES,
    RARE_DISEASE_PHENOTYPES,
    DiagnosticWorkflowType,
    DiagnosticResult,
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
class RareDiseaseSearchResult:
    """A single search result from a Milvus collection.

    Attributes:
        collection: Source collection name (e.g. 'rd_phenotypes').
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
    "rd_phenotypes": {
        "weight": settings.WEIGHT_PHENOTYPES,
        "label": "Phenotype",
        "text_field": "phenotype_description",
        "title_field": "hpo_term",
        "filterable_fields": ["hpo_id", "category", "onset"],
    },
    "rd_diseases": {
        "weight": settings.WEIGHT_DISEASES,
        "label": "Disease",
        "text_field": "disease_description",
        "title_field": "disease_name",
        "filterable_fields": ["omim_id", "orpha_id", "inheritance", "prevalence"],
    },
    "rd_genes": {
        "weight": settings.WEIGHT_GENES,
        "label": "Gene",
        "text_field": "gene_summary",
        "title_field": "gene_symbol",
        "filterable_fields": ["gene_symbol", "inheritance", "clingen_validity"],
    },
    "rd_variants": {
        "weight": settings.WEIGHT_VARIANTS,
        "label": "Variant",
        "text_field": "variant_interpretation",
        "title_field": "variant_hgvs",
        "filterable_fields": ["gene", "acmg_class", "clinvar_significance", "review_status"],
    },
    "rd_literature": {
        "weight": settings.WEIGHT_LITERATURE,
        "label": "Literature",
        "text_field": "abstract",
        "title_field": "title",
        "filterable_fields": ["study_type", "disease_area", "year"],
    },
    "rd_trials": {
        "weight": settings.WEIGHT_TRIALS,
        "label": "ClinicalTrial",
        "text_field": "trial_summary",
        "title_field": "trial_title",
        "filterable_fields": ["phase", "status", "condition", "intervention_type"],
    },
    "rd_therapies": {
        "weight": settings.WEIGHT_THERAPIES,
        "label": "Therapy",
        "text_field": "therapy_description",
        "title_field": "therapy_name",
        "filterable_fields": ["therapy_type", "approval_status", "condition"],
    },
    "rd_case_reports": {
        "weight": settings.WEIGHT_CASE_REPORTS,
        "label": "CaseReport",
        "text_field": "case_narrative",
        "title_field": "case_title",
        "filterable_fields": ["diagnosis", "gene", "phenotypes"],
    },
    "rd_guidelines": {
        "weight": settings.WEIGHT_GUIDELINES,
        "label": "Guideline",
        "text_field": "recommendation",
        "title_field": "guideline_title",
        "filterable_fields": ["issuing_body", "disease_area", "year"],
    },
    "rd_pathways": {
        "weight": settings.WEIGHT_PATHWAYS,
        "label": "Pathway",
        "text_field": "pathway_description",
        "title_field": "pathway_name",
        "filterable_fields": ["pathway_type", "metabolites"],
    },
    "rd_registries": {
        "weight": settings.WEIGHT_REGISTRIES,
        "label": "Registry",
        "text_field": "registry_description",
        "title_field": "registry_name",
        "filterable_fields": ["disease", "country", "enrollment"],
    },
    "rd_natural_history": {
        "weight": settings.WEIGHT_NATURAL_HISTORY,
        "label": "NaturalHistory",
        "text_field": "natural_history_summary",
        "title_field": "disease_name",
        "filterable_fields": ["disease", "age_of_onset", "progression"],
    },
    "rd_newborn_screening": {
        "weight": settings.WEIGHT_NEWBORN_SCREENING,
        "label": "NewbornScreening",
        "text_field": "screening_description",
        "title_field": "condition_name",
        "filterable_fields": ["analyte", "rusp_status", "confirmatory_test"],
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
# RARE DISEASE RAG ENGINE
# =====================================================================

class RareDiseaseRAGEngine:
    """Multi-collection RAG engine for rare disease diagnostics.

    Searches across all 14 rare-disease-specific Milvus collections plus the
    shared genomic_evidence collection. Supports workflow-specific weight
    boosting, parallel search, query expansion, patient context injection,
    and multi-turn conversation memory.

    Features:
    - Parallel search via ThreadPoolExecutor (14 collections)
    - Settings-driven weights and parameters
    - Workflow-based dynamic weight boosting (13 diagnostic workflows)
    - Milvus field-based filtering (gene, inheritance, HPO, OMIM)
    - Citation relevance scoring (high/medium/low)
    - Cross-collection entity linking (gene-disease-phenotype)
    - Phenotype overlap scoring for differential diagnosis
    - Conversation memory context injection
    - Patient context for personalised diagnostic assessment
    - Confidence scoring based on evidence diversity

    Usage:
        engine = RareDiseaseRAGEngine(milvus_client, embedding_model, llm_client)
        response = engine.query("Infant with hypotonia and cardiomyopathy")
        results = engine.search("CFTR F508del variant interpretation")
    """

    def __init__(
        self,
        milvus_client=None,
        embedding_model=None,
        llm_client=None,
        session_id: str = "default",
    ):
        """Initialize the RareDiseaseRAGEngine.

        Args:
            milvus_client: Connected Milvus client with access to all
                rare disease collections. If None, search operations will
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
        workflow: Optional[DiagnosticWorkflowType] = None,
        top_k: int = 5,
        patient_context: Optional[dict] = None,
    ) -> DiagnosticResult:
        """Main query method: expand -> search -> synthesise.

        Performs the full RAG pipeline: parallel multi-collection search
        with workflow-specific weighting, result reranking, LLM synthesis
        with patient context, and confidence scoring.

        Args:
            question: Natural language rare disease diagnostic question.
            workflow: Optional DiagnosticWorkflowType to apply domain-specific
                collection weight boosting. If None, weights are auto-detected
                or base defaults are used.
            top_k: Maximum results to return per collection.
            patient_context: Optional dict with patient-specific data
                (age, sex, phenotypes, hpo_terms, variants, gene_panel_results,
                family_history, consanguinity, ethnicity, newborn_screen_results)
                for personalised diagnostic assessment.

        Returns:
            DiagnosticResult with synthesised answer, search results, citations,
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
            response = DiagnosticResult(
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
    ) -> List[RareDiseaseSearchResult]:
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
            List of RareDiseaseSearchResult sorted by weighted score descending.
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
                (e.g. 'gene == "CFTR"').

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
    ) -> List[RareDiseaseSearchResult]:
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
            List of RareDiseaseSearchResult sorted by weighted score descending.
        """
        all_results: List[RareDiseaseSearchResult] = []
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
                                         "disease_description", "gene_summary",
                                         "variant_interpretation", "case_narrative",
                                         "therapy_description", "pathway_description",
                                         "screening_description", "text_chunk"):
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

                    result = RareDiseaseSearchResult(
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
        unique_results: List[RareDiseaseSearchResult] = []
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
        results: List[RareDiseaseSearchResult],
        query: str,
    ) -> List[RareDiseaseSearchResult]:
        """Rerank results based on relevance to original query.

        Applies heuristic boosts for:
        - Disease results matching query gene symbols
        - Variant results with ACMG classification data
        - Phenotype results matching HPO terms in query
        - Results with OMIM/Orphanet identifiers (evidence-based)
        - Gene results matching detected gene symbols
        - Therapy results for actionable conditions
        - Guideline results from recognized bodies (ACMG, ACOG, NORD)

        Args:
            results: Raw search results to rerank.
            query: Original query text for relevance matching.

        Returns:
            Reranked list of RareDiseaseSearchResult.
        """
        query_lower = query.lower()
        query_upper = query.upper()
        query_terms = set(query_lower.split())

        for result in results:
            boost = 0.0

            # Boost disease results when query mentions specific diseases
            if result.collection == "rd_diseases":
                omim_id = result.metadata.get("omim_id", "")
                if omim_id:
                    boost += 0.05
                # Boost if gene symbol present in query
                gene = result.metadata.get("gene_symbol", "")
                if gene and gene in query_upper:
                    boost += 0.15

            # Boost variant results with ACMG classification
            if result.collection == "rd_variants":
                acmg = result.metadata.get("acmg_class", "")
                if acmg:
                    boost += 0.05
                # Boost pathogenic/likely pathogenic
                if acmg in ("pathogenic", "likely_pathogenic"):
                    boost += 0.10
                # Boost if variant notation present in query
                hgvs = result.metadata.get("variant_hgvs", "")
                if hgvs and hgvs.lower() in query_lower:
                    boost += 0.15

            # Boost phenotype results when HPO terms are in query
            if result.collection == "rd_phenotypes":
                hpo_id = result.metadata.get("hpo_id", "")
                if hpo_id and hpo_id in query:
                    boost += 0.15
                boost += 0.05  # General phenotype boost

            # Boost gene results when gene symbols are in query
            if result.collection == "rd_genes":
                gene_sym = result.metadata.get("gene_symbol", "")
                if gene_sym and gene_sym in query_upper:
                    boost += 0.15
                # ClinGen validity boost
                validity = result.metadata.get("clingen_validity", "")
                if validity in ("definitive", "strong"):
                    boost += 0.10

            # Boost guideline results
            if result.collection == "rd_guidelines":
                boost += 0.05
                issuer = result.metadata.get("issuing_body", "").upper()
                acmg_bodies = {"ACMG", "ACOG", "NORD", "NICE", "ESHG"}
                if any(body in issuer for body in acmg_bodies):
                    boost += 0.10

            # Boost therapy results for therapy-related queries
            therapy_terms = {"treatment", "therapy", "drug", "ert", "gene therapy",
                             "orphan drug", "medication", "enzyme replacement"}
            if result.collection == "rd_therapies":
                if query_terms & therapy_terms:
                    boost += 0.10

            # Boost trial results for trial queries
            trial_terms = {"trial", "clinical trial", "nct", "recruiting",
                           "study", "investigational"}
            if result.collection == "rd_trials":
                if query_terms & trial_terms:
                    boost += 0.10

            # Boost newborn screening results for screening queries
            nbs_terms = {"newborn", "screening", "nbs", "dried blood spot",
                         "tandem mass", "rusp"}
            if result.collection == "rd_newborn_screening":
                if query_terms & nbs_terms:
                    boost += 0.10

            # Boost metabolic pathway results for metabolic queries
            metabolic_terms = {"metabolic", "ammonia", "lactate", "acidosis",
                               "amino acid", "organic acid", "urea cycle"}
            if result.collection == "rd_pathways":
                if query_terms & metabolic_terms:
                    boost += 0.10

            # Boost case reports for phenotype-driven queries
            if result.collection == "rd_case_reports":
                if len(query_terms & {"presenting", "phenotype", "case",
                                      "patient", "infant", "child"}) > 0:
                    boost += 0.05

            # Boost results with high relevance
            if result.relevance == "high":
                boost += 0.10
            elif result.relevance == "medium":
                boost += 0.05

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
        results: List[RareDiseaseSearchResult],
    ) -> List[RareDiseaseSearchResult]:
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
        results: List[RareDiseaseSearchResult],
        workflow: Optional[DiagnosticWorkflowType] = None,
        patient_context: Optional[dict] = None,
    ) -> DiagnosticResult:
        """Use LLM to synthesise search results into a diagnostic response.

        Builds a structured prompt with retrieved evidence, patient context,
        conversation history, and workflow-specific instructions. Generates
        a grounded answer via the configured LLM.

        Args:
            question: Original diagnostic question.
            results: Ranked search results for context.
            workflow: Optional workflow for instruction tuning.
            patient_context: Optional patient-specific data dict.

        Returns:
            DiagnosticResult with synthesised answer and metadata.
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
            f"Please provide a comprehensive, evidence-based rare disease "
            f"diagnostic assessment grounded in the retrieved evidence above. "
            f"Follow the system prompt instructions for HPO term format, "
            f"OMIM/Orphanet citations, ACMG variant classification, severity "
            f"badges, and structured output sections.\n\n"
            f"Cite sources using clickable markdown links where OMIM numbers "
            f"are available: [OMIM:219700](https://omim.org/entry/219700). "
            f"For HPO terms, use [HP:0001250](https://hpo.jax.org/browse/term/HP:0001250). "
            f"For PubMed evidence, use [PMID:12345678](https://pubmed.ncbi.nlm.nih.gov/12345678/). "
            f"For collection-sourced evidence, use [Collection:record-id]. "
            f"Prioritise [high relevance] citations and pathogenic/likely pathogenic variants."
        )

        answer = self.llm.generate(
            prompt=prompt,
            system_prompt=RARE_DISEASE_SYSTEM_PROMPT,
            max_tokens=2048,
            temperature=0.7,
        )

        return DiagnosticResult(
            question=question,
            answer=answer,
            results=results,
            workflow=workflow,
        )

    def _build_context(
        self,
        results: List[RareDiseaseSearchResult],
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
        by_collection: Dict[str, List[RareDiseaseSearchResult]] = {}
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

    def _format_citation_link(self, result: RareDiseaseSearchResult) -> str:
        """Format a citation with clickable URL where possible.

        Args:
            result: Search result to format citation for.

        Returns:
            Markdown-formatted citation string.
        """
        label = result.metadata.get("collection_label", result.collection)
        record_id = result.record_id

        # OMIM references
        omim_id = result.metadata.get("omim_id", "")
        if omim_id:
            omim_num = omim_id.replace("OMIM:", "").strip()
            return (
                f"[{label}:OMIM {omim_num}]"
                f"(https://omim.org/entry/{omim_num})"
            )

        # Orphanet references
        orpha_id = result.metadata.get("orpha_id", "")
        if orpha_id:
            orpha_num = orpha_id.replace("ORPHA:", "").strip()
            return (
                f"[{label}:ORPHA {orpha_num}]"
                f"(https://www.orpha.net/consor/cgi-bin/OC_Exp.php?lng=en&Expert={orpha_num})"
            )

        # HPO references
        hpo_id = result.metadata.get("hpo_id", "")
        if hpo_id:
            return (
                f"[{label}:{hpo_id}]"
                f"(https://hpo.jax.org/browse/term/{hpo_id})"
            )

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

        # Gene symbol
        gene_symbol = result.metadata.get("gene_symbol", "")
        if gene_symbol:
            return f"[{label}:{gene_symbol}]"

        return f"[{label}:{record_id}]"

    def _format_patient_context(self, patient_context: Optional[dict]) -> str:
        """Format patient context for LLM prompt injection.

        Used primarily for phenotype-driven diagnostic workflows.

        Args:
            patient_context: Optional patient data dict with keys like
                age, sex, phenotypes, hpo_terms, variants, gene_panel_results,
                family_history, consanguinity, ethnicity, newborn_screen_results.

        Returns:
            Formatted patient context section or empty string.
        """
        if not patient_context:
            return ""

        lines = ["### Patient Context\n"]

        field_labels = {
            "age": "Age",
            "sex": "Sex",
            "gestational_age": "Gestational Age",
            "ethnicity": "Ethnicity",
            "consanguinity": "Parental Consanguinity",
            "phenotypes": "Clinical Phenotypes",
            "hpo_terms": "HPO Terms",
            "variants": "Genomic Variants",
            "gene_panel_results": "Gene Panel Results",
            "exome_results": "Exome/Genome Results",
            "family_history": "Family History",
            "pedigree_notes": "Pedigree Notes",
            "newborn_screen_results": "Newborn Screen Results",
            "metabolic_labs": "Metabolic Laboratories",
            "enzyme_assays": "Enzyme Assay Results",
            "imaging": "Imaging Findings",
            "medications": "Current Medications",
            "prior_diagnoses": "Prior Diagnoses",
            "developmental_milestones": "Developmental Milestones",
            "growth_parameters": "Growth Parameters",
            "dietary_restrictions": "Dietary Restrictions",
            "organ_involvement": "Organ System Involvement",
        }

        for key, label in field_labels.items():
            value = patient_context.get(key)
            if value is not None:
                if isinstance(value, list):
                    value = ", ".join(str(v) for v in value)
                elif isinstance(value, dict):
                    value = "; ".join(f"{k}: {v}" for k, v in value.items())
                elif isinstance(value, bool):
                    value = "Yes" if value else "No"
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
        workflow: Optional[DiagnosticWorkflowType],
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
            DiagnosticWorkflowType.PHENOTYPE_DRIVEN: (
                "### Workflow: Phenotype-Driven Diagnosis\n"
                "Focus on: HPO term analysis, phenotype-gene associations, organ system "
                "categorization, syndromic pattern recognition, dysmorphology assessment, "
                "and ordered differential diagnosis based on phenotype overlap score. "
                "Prioritize conditions with available therapies (treatability bias). "
                "Recommend confirmatory diagnostic testing.\n\n"
            ),
            DiagnosticWorkflowType.VARIANT_INTERPRETATION: (
                "### Workflow: Variant Interpretation (ACMG/AMP)\n"
                "Focus on: ACMG/AMP 2015 criteria application (PVS, PS, PM, PP, BS, BP, BA), "
                "ClinVar review status and assertion counts, population frequency databases "
                "(gnomAD, ClinGen), computational predictions (REVEL, CADD, SpliceAI), "
                "functional studies, segregation data, and clinical correlation. Provide "
                "final ACMG classification with evidence summary.\n\n"
            ),
            DiagnosticWorkflowType.DIFFERENTIAL_DIAGNOSIS: (
                "### Workflow: Differential Diagnosis\n"
                "Focus on: ranked differential list with OMIM numbers, phenotype overlap "
                "analysis using HPO semantic similarity, distinguishing clinical features, "
                "inheritance pattern matching, gene-disease validity (ClinGen), diagnostic "
                "confirmation methods, and recommended testing cascade (tier 1: targeted, "
                "tier 2: panel, tier 3: exome/genome, tier 4: functional studies).\n\n"
            ),
            DiagnosticWorkflowType.GENE_THERAPY_ELIGIBILITY: (
                "### Workflow: Gene Therapy Eligibility Assessment\n"
                "Focus on: approved gene therapies and active trials for the condition, "
                "eligibility criteria (age, genotype, disease severity, organ function, "
                "pre-existing AAV antibodies), vector type and route of administration, "
                "efficacy and safety data from clinical trials, regulatory designations "
                "(RMAT, Breakthrough, Orphan), and alternative therapeutic options.\n\n"
            ),
            DiagnosticWorkflowType.NEWBORN_SCREENING: (
                "### Workflow: Newborn Screening Triage\n"
                "Focus on: RUSP inclusion status, screening analyte and methodology, "
                "confirmatory testing protocol (timing, specimen requirements), "
                "false positive/negative rates, genotype-phenotype correlations for "
                "screening-detected patients, early intervention protocols, and "
                "long-term follow-up recommendations.\n\n"
            ),
            DiagnosticWorkflowType.METABOLIC_WORKUP: (
                "### Workflow: Metabolic Crisis Workup\n"
                "Focus on: acute stabilization priorities (ammonia scavengers, glucose, "
                "dialysis indications), biochemical pathway analysis, first-tier metabolic "
                "labs (ammonia, lactate, amino acids, acylcarnitines, urine organic acids), "
                "second-tier confirmatory testing, differential by biochemical pattern, "
                "and acute vs chronic management strategies.\n\n"
            ),
            DiagnosticWorkflowType.CARRIER_SCREENING: (
                "### Workflow: Carrier Screening Assessment\n"
                "Focus on: carrier frequency by ethnicity/population, recommended "
                "screening panels (ACOG, ACMG), detection rate and residual risk "
                "calculations, partner testing strategy, reproductive implications, "
                "and genetic counseling considerations for identified carriers.\n\n"
            ),
            DiagnosticWorkflowType.PRENATAL_DIAGNOSIS: (
                "### Workflow: Prenatal Diagnosis\n"
                "Focus on: available prenatal testing options (CVS, amniocentesis, NIPT, "
                "PGT-M), optimal timing for each modality, sensitivity and specificity, "
                "fetal imaging findings suggestive of the condition, reproductive decision "
                "support (non-directive), and postnatal management planning.\n\n"
            ),
            DiagnosticWorkflowType.NATURAL_HISTORY: (
                "### Workflow: Natural History Assessment\n"
                "Focus on: disease progression timeline, major clinical milestones, "
                "organ-specific complications and their typical age of onset, survival "
                "statistics, functional outcome measures, genotype-phenotype correlations "
                "affecting prognosis, surveillance guidelines, and quality of life data.\n\n"
            ),
            DiagnosticWorkflowType.THERAPY_SELECTION: (
                "### Workflow: Therapy Selection\n"
                "Focus on: approved therapies (ERT, SRT, ASO, gene therapy, chaperone), "
                "evidence level for each therapeutic option, treatment initiation criteria, "
                "monitoring parameters, expected treatment response, adverse effects, "
                "FDA/EMA approval status, orphan drug and breakthrough designations, "
                "and emerging therapies in clinical trials.\n\n"
            ),
            DiagnosticWorkflowType.CLINICAL_TRIAL_MATCHING: (
                "### Workflow: Clinical Trial Matching\n"
                "Focus on: active recruiting trials for the condition, eligibility criteria "
                "alignment with patient profile, therapeutic modality (gene therapy, ASO, "
                "ERT, small molecule, cell therapy), trial phase and primary endpoints, "
                "trial site locations, and expanded access / compassionate use programs.\n\n"
            ),
            DiagnosticWorkflowType.GENETIC_COUNSELING: (
                "### Workflow: Genetic Counseling\n"
                "Focus on: inheritance pattern explanation (AD, AR, XL, mito, digenic), "
                "recurrence risk calculations, variable expressivity and reduced penetrance, "
                "anticipation in repeat disorders, cascade testing recommendations, "
                "reproductive options (PGT-M, prenatal testing, donor gametes), "
                "psychosocial support resources, and patient advocacy organizations.\n\n"
            ),
        }

        return instructions.get(workflow, "")

    # ==================================================================
    # CITATIONS & CONFIDENCE
    # ==================================================================

    def _extract_citations(
        self,
        results: List[RareDiseaseSearchResult],
    ) -> List[dict]:
        """Extract and format citations from search results.

        Generates a structured citation list from all results, including
        OMIM links, HPO links, PMID links, and gene references.

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
            omim_id = result.metadata.get("omim_id", "")
            if omim_id:
                omim_num = omim_id.replace("OMIM:", "").strip()
                cite["url"] = f"https://omim.org/entry/{omim_num}"
                cite["id"] = omim_id

            orpha_id = result.metadata.get("orpha_id", "")
            if orpha_id:
                orpha_num = orpha_id.replace("ORPHA:", "").strip()
                cite["url"] = (
                    f"https://www.orpha.net/consor/cgi-bin/OC_Exp.php"
                    f"?lng=en&Expert={orpha_num}"
                )
                cite["id"] = orpha_id

            hpo_id = result.metadata.get("hpo_id", "")
            if hpo_id and not cite["url"]:
                cite["url"] = f"https://hpo.jax.org/browse/term/{hpo_id}"
                cite["id"] = hpo_id

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
        results: List[RareDiseaseSearchResult],
    ) -> float:
        """Score overall confidence based on result quality.

        Confidence is based on:
        - Number of high-relevance results
        - Collection diversity
        - Average similarity score
        - Presence of guideline/gene evidence

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

        # Factor 4: Gene/guideline evidence present (0-0.2)
        has_authoritative = any(
            r.collection in ("rd_guidelines", "rd_genes", "rd_diseases")
            for r in results
        )
        authoritative_score = 0.2 if has_authoritative else 0.0

        confidence = (
            relevance_score + diversity_score + quality_score + authoritative_score
        )
        return round(min(confidence, 1.0), 3)

    # ==================================================================
    # ENTITY & DISEASE SEARCH
    # ==================================================================

    def find_related(
        self,
        entity: str,
        entity_type: str = "disease",
        top_k: int = 5,
    ) -> List[RareDiseaseSearchResult]:
        """Find related entities across collections.

        Searches relevant collections for evidence related to a clinical
        entity (disease, gene, phenotype, therapy). Useful for building
        entity profiles and cross-referencing.

        Args:
            entity: Entity name (e.g. 'Gaucher disease', 'CFTR', 'hypotonia').
            entity_type: Entity category for targeted search:
                'disease', 'gene', 'phenotype', 'therapy', 'variant'.
            top_k: Maximum results per collection.

        Returns:
            List of RareDiseaseSearchResult from all relevant collections.
        """
        type_collection_map = {
            "disease": [
                "rd_diseases", "rd_genes", "rd_phenotypes",
                "rd_natural_history", "rd_literature",
            ],
            "gene": [
                "rd_genes", "rd_variants", "rd_diseases",
                "genomic_evidence", "rd_literature",
            ],
            "phenotype": [
                "rd_phenotypes", "rd_diseases", "rd_case_reports",
                "rd_natural_history", "rd_literature",
            ],
            "therapy": [
                "rd_therapies", "rd_trials", "rd_guidelines",
                "rd_diseases", "rd_literature",
            ],
            "variant": [
                "rd_variants", "rd_genes", "genomic_evidence",
                "rd_diseases", "rd_literature",
            ],
        }

        collections = type_collection_map.get(entity_type, get_all_collection_names())
        return self.search(entity, collections=collections, top_k=top_k)

    def get_disease_details(self, omim_id: str) -> Optional[dict]:
        """Retrieve details for a specific disease by OMIM ID.

        Searches the rd_diseases collection for a specific disease
        using scalar filtering on the omim_id field.

        Args:
            omim_id: OMIM MIM number (e.g. 'OMIM:219700' or '219700').

        Returns:
            Disease details dict or None if not found.
        """
        if not self.milvus:
            raise RuntimeError("Milvus client not configured.")

        # Sanitize input
        safe_id = omim_id.strip()
        if not _SAFE_FILTER_RE.match(safe_id):
            logger.warning("Rejected unsafe OMIM ID: %r", safe_id)
            return None

        try:
            # Use a generic query vector for filtered search
            query_vector = self._embed_query(f"rare disease OMIM {omim_id}")
            filter_expr = f'omim_id == "{safe_id}"'

            raw_results = self._search_collection(
                "rd_diseases", query_vector, top_k=1,
                filter_expr=filter_expr,
            )

            if raw_results:
                return raw_results[0]
            return None

        except Exception as exc:
            logger.warning("Failed to get disease details for %s: %s", omim_id, exc)
            return None

    def get_variant_interpretation(
        self,
        variant_hgvs: str,
        gene: Optional[str] = None,
    ) -> List[RareDiseaseSearchResult]:
        """Retrieve variant interpretation evidence.

        Searches variant and gene collections for interpretation evidence
        on a specific genomic variant, optionally filtered by gene.

        Args:
            variant_hgvs: HGVS variant notation
                (e.g. 'NM_000492.4(CFTR):c.1521_1523del').
            gene: Optional gene symbol for filtering (e.g. 'CFTR').

        Returns:
            List of RareDiseaseSearchResult with variant evidence.
        """
        collections = ["rd_variants", "rd_genes", "genomic_evidence",
                        "rd_diseases", "rd_literature"]
        query = f"variant interpretation {variant_hgvs}"
        if gene:
            query += f" gene {gene}"
        return self.search(query, collections=collections, top_k=10)

    def search_by_phenotypes(
        self,
        hpo_terms: List[str],
        top_k: int = 10,
    ) -> List[RareDiseaseSearchResult]:
        """Search for diseases matching a set of HPO phenotype terms.

        Constructs a phenotype-based query from multiple HPO terms to
        find diseases with the highest phenotype overlap.

        Args:
            hpo_terms: List of HPO term IDs or names
                (e.g. ['HP:0001252', 'HP:0001638', 'HP:0001433']).
            top_k: Maximum results per collection.

        Returns:
            List of RareDiseaseSearchResult prioritizing phenotype
            and disease collections.
        """
        query = "differential diagnosis for patient presenting with: " + ", ".join(hpo_terms)
        collections = ["rd_phenotypes", "rd_diseases", "rd_genes",
                        "rd_case_reports", "rd_natural_history"]
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
        workflow: Optional[DiagnosticWorkflowType] = None,
    ) -> Dict[str, float]:
        """Compute collection weights with optional workflow boosting.

        When a workflow is specified, applies boost multipliers from
        WORKFLOW_COLLECTION_BOOST on top of the base weights from
        settings. Weights are then renormalized to sum to ~1.0.

        Args:
            workflow: Optional DiagnosticWorkflowType for boosting.

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
        all expected rare disease collections exist and are loaded.

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
