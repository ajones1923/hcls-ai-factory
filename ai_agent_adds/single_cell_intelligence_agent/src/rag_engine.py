"""Multi-collection RAG engine for Single-Cell Intelligence Agent.

Searches across all 12 single-cell-specific Milvus collections simultaneously
using parallel ThreadPoolExecutor, synthesises findings with single-cell
knowledge augmentation, and generates grounded LLM responses with evidence
citations.

Extends the pattern from: rag-chat-pipeline/src/rag_engine.py

Features:
- Parallel search via ThreadPoolExecutor (11 SC + 1 shared genomic collection)
- Settings-driven weights and parameters from config/settings.py
- Workflow-based dynamic weight boosting per SCWorkflowType
- Milvus field-based filtering (modality, condition, cell_type, tissue)
- Citation relevance scoring (high/medium/low) with PMID/DOI link formatting
- Cross-collection entity linking for comprehensive single-cell queries
- Cell Ontology and atlas reference retrieval
- Conversation memory for multi-turn analytical consultations
- Sample context injection for personalised analysis support
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
    SC_SYSTEM_PROMPT,
    WORKFLOW_COLLECTION_BOOST,
    SCWorkflowType,
    SCResponse,
)
from .knowledge import (
    CELL_TYPE_ATLAS,
    TME_PROFILES,
    DRUG_SENSITIVITY_DATABASE,
    LIGAND_RECEPTOR_PAIRS,
    CANCER_TME_ATLAS,
    IMMUNE_SIGNATURES,
    MARKER_GENE_DATABASE,
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
class SCSearchResult:
    """A single search result from a Milvus collection.

    Attributes:
        collection: Source collection name (e.g. 'sc_cell_types').
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
    "sc_cell_types": {
        "weight": settings.WEIGHT_CELL_TYPES,
        "label": "CellType",
        "text_field": "cell_description",
        "title_field": "cell_type",
        "filterable_fields": ["cell_ontology_id", "tissue", "condition", "domain"],
    },
    "sc_markers": {
        "weight": settings.WEIGHT_MARKERS,
        "label": "Marker",
        "text_field": "marker_description",
        "title_field": "gene_symbol",
        "filterable_fields": ["cell_type", "tissue", "source_db", "specificity"],
    },
    "sc_spatial": {
        "weight": settings.WEIGHT_SPATIAL,
        "label": "Spatial",
        "text_field": "spatial_findings",
        "title_field": "tissue_region",
        "filterable_fields": ["technology", "tissue", "spatial_domain", "condition"],
    },
    "sc_tme": {
        "weight": settings.WEIGHT_TME,
        "label": "TME",
        "text_field": "tme_description",
        "title_field": "tme_subtype",
        "filterable_fields": ["classification", "cancer_type", "immune_phenotype",
                              "therapy_response"],
    },
    "sc_drug_response": {
        "weight": settings.WEIGHT_DRUG_RESPONSE,
        "label": "DrugResponse",
        "text_field": "response_description",
        "title_field": "drug_name",
        "filterable_fields": ["cancer_type", "mechanism", "cell_type", "sensitivity"],
    },
    "sc_literature": {
        "weight": settings.WEIGHT_LITERATURE,
        "label": "Literature",
        "text_field": "abstract",
        "title_field": "title",
        "filterable_fields": ["study_type", "technology", "tissue", "year"],
    },
    "sc_methods": {
        "weight": settings.WEIGHT_METHODS,
        "label": "Method",
        "text_field": "method_description",
        "title_field": "method_name",
        "filterable_fields": ["category", "scalability", "benchmark_rank"],
    },
    "sc_datasets": {
        "weight": settings.WEIGHT_DATASETS,
        "label": "Dataset",
        "text_field": "dataset_description",
        "title_field": "dataset_name",
        "filterable_fields": ["technology", "tissue", "species", "cell_count"],
    },
    "sc_trajectories": {
        "weight": settings.WEIGHT_TRAJECTORIES,
        "label": "Trajectory",
        "text_field": "trajectory_description",
        "title_field": "lineage",
        "filterable_fields": ["tissue", "method", "start_cell", "end_cell"],
    },
    "sc_pathways": {
        "weight": settings.WEIGHT_PATHWAYS,
        "label": "Pathway",
        "text_field": "pathway_description",
        "title_field": "pathway_name",
        "filterable_fields": ["cell_type_a", "cell_type_b", "ligand", "receptor"],
    },
    "sc_clinical": {
        "weight": settings.WEIGHT_CLINICAL,
        "label": "Clinical",
        "text_field": "clinical_description",
        "title_field": "clinical_context",
        "filterable_fields": ["cancer_type", "treatment", "response", "timepoint"],
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
# SINGLE-CELL RAG ENGINE
# =====================================================================

class SingleCellRAGEngine:
    """Multi-collection RAG engine for single-cell intelligence.

    Searches across all 12 single-cell-specific Milvus collections plus the
    shared genomic_evidence collection. Supports workflow-specific weight
    boosting, parallel search, query expansion, sample context injection,
    and multi-turn conversation memory.

    Features:
    - Parallel search via ThreadPoolExecutor (12 collections)
    - Settings-driven weights and parameters
    - Workflow-based dynamic weight boosting (11 SC workflows)
    - Milvus field-based filtering (cell_type, tissue, technology, condition)
    - Citation relevance scoring (high/medium/low)
    - Cross-collection entity linking
    - Cell Ontology and atlas reference retrieval
    - Conversation memory context injection
    - Sample context for personalised analysis support
    - Confidence scoring based on evidence diversity

    Usage:
        engine = SingleCellRAGEngine(milvus_client, embedding_model, llm_client)
        response = engine.query("Classify TME in NSCLC scRNA-seq dataset")
        results = engine.search("CD8 T cell exhaustion markers PD-1 TIM-3")
    """

    def __init__(
        self,
        milvus_client=None,
        embedding_model=None,
        llm_client=None,
        session_id: str = "default",
    ):
        """Initialize the SingleCellRAGEngine.

        Args:
            milvus_client: Connected Milvus client with access to all
                SC collections. If None, search operations will
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
        workflow: Optional[SCWorkflowType] = None,
        top_k: int = 5,
        patient_context: Optional[dict] = None,
    ) -> SCResponse:
        """Main query method: expand -> search -> synthesise.

        Performs the full RAG pipeline: parallel multi-collection search
        with workflow-specific weighting, result reranking, LLM synthesis
        with sample context, and confidence scoring.

        Args:
            question: Natural language single-cell analysis question.
            workflow: Optional SCWorkflowType to apply domain-specific
                collection weight boosting. If None, weights are auto-detected
                or base defaults are used.
            top_k: Maximum results to return per collection.
            patient_context: Optional dict with sample/patient-specific data
                (tissue, technology, cell_count, condition, treatment,
                timepoint, markers_of_interest, genomic_data)
                for personalised analysis support.

        Returns:
            SCResponse with synthesised answer, search results, citations,
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
            response = SCResponse(
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
    ) -> List[SCSearchResult]:
        """Search across multiple collections with weighted scoring.

        Embeds the query, runs parallel Milvus searches across all specified
        collections, applies collection weights, and returns merged ranked
        results.

        Args:
            question: Natural language search query.
            collections: Optional list of collection names to search.
                If None, all 12 collections are searched.
            top_k: Maximum results per collection.

        Returns:
            List of SCSearchResult sorted by weighted score descending.
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
                (e.g. 'cell_type == "CD8+ T cell"').

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
    ) -> List[SCSearchResult]:
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
            List of SCSearchResult sorted by weighted score descending.
        """
        all_results: List[SCSearchResult] = []
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
                        for fallback in ("abstract", "content", "cell_description",
                                         "marker_description", "spatial_findings",
                                         "tme_description", "response_description",
                                         "method_description", "dataset_description",
                                         "trajectory_description", "pathway_description",
                                         "clinical_description", "text_chunk"):
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

                    result = SCSearchResult(
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
        unique_results: List[SCSearchResult] = []
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
        results: List[SCSearchResult],
        query: str,
    ) -> List[SCSearchResult]:
        """Rerank results based on relevance to original query.

        Applies heuristic boosts for:
        - Cell type results matching query cell types
        - Results from spatial collections for spatial queries
        - Results with high citation relevance
        - PMID/DOI-bearing results (evidence-based)
        - Results matching detected biomarker or cell type terms

        Args:
            results: Raw search results to rerank.
            query: Original query text for relevance matching.

        Returns:
            Reranked list of SCSearchResult.
        """
        query_lower = query.lower()
        query_terms = set(query_lower.split())

        for result in results:
            boost = 0.0

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

            # Boost cell type results for annotation queries
            annotation_terms = {"cell type", "annotation", "annotate", "cluster",
                                "identity", "marker", "canonical", "cellmarker",
                                "panglaodb", "celltypist", "azimuth"}
            if result.collection == "sc_cell_types":
                if query_terms & annotation_terms:
                    boost += 0.10

            # Boost marker results for biomarker queries
            marker_terms = {"marker", "biomarker", "gene", "signature",
                            "expression", "differential", "deg",
                            "surface marker", "canonical marker"}
            if result.collection == "sc_markers":
                if query_terms & marker_terms:
                    boost += 0.10

            # Boost spatial results for spatial queries
            spatial_terms = {"spatial", "visium", "merfish", "slide-seq",
                             "niche", "architecture", "colocalization",
                             "neighborhood", "tissue", "region", "domain"}
            if result.collection == "sc_spatial":
                if query_terms & spatial_terms:
                    boost += 0.10

            # Boost TME results for microenvironment queries
            tme_terms = {"tme", "microenvironment", "immune", "hot", "cold",
                         "excluded", "immunosuppressive", "infiltrate",
                         "tumor", "immune desert", "inflamed"}
            if result.collection == "sc_tme":
                if query_terms & tme_terms:
                    boost += 0.10

            # Boost drug response results for therapy queries
            drug_terms = {"drug", "therapy", "response", "resistance",
                          "sensitivity", "ic50", "treatment",
                          "depmap", "chemotherapy", "immunotherapy",
                          "checkpoint", "target"}
            if result.collection == "sc_drug_response":
                if query_terms & drug_terms:
                    boost += 0.10

            # Boost trajectory results for differentiation queries
            traj_terms = {"trajectory", "pseudotime", "lineage",
                          "differentiation", "fate", "velocity",
                          "monocle", "scvelo", "paga", "cellrank",
                          "progenitor", "stem"}
            if result.collection == "sc_trajectories":
                if query_terms & traj_terms:
                    boost += 0.10

            # Boost pathway results for signaling queries
            pathway_terms = {"ligand", "receptor", "signaling", "pathway",
                             "communication", "cellchat", "liana",
                             "nichenet", "interaction", "crosstalk",
                             "paracrine", "cytokine", "chemokine"}
            if result.collection == "sc_pathways":
                if query_terms & pathway_terms:
                    boost += 0.10

            # Boost clinical results for treatment queries
            clinical_terms = {"clinical", "patient", "treatment",
                              "response", "relapse", "mrd",
                              "longitudinal", "pre-treatment",
                              "post-treatment", "cart", "car-t"}
            if result.collection == "sc_clinical":
                if query_terms & clinical_terms:
                    boost += 0.10

            # Boost method results for analytical queries
            method_terms = {"method", "algorithm", "pipeline",
                            "benchmark", "tool", "software",
                            "clustering", "integration", "batch",
                            "harmony", "scvi", "rapids", "gpu"}
            if result.collection == "sc_methods":
                if query_terms & method_terms:
                    boost += 0.10

            # Boost genomic results for genetic queries
            genetic_terms = {"gene", "mutation", "variant", "cnv",
                             "copy number", "snv", "infercnv",
                             "copybat", "clonal", "subclone",
                             "allele", "genotype"}
            if result.collection == "genomic_evidence":
                if query_terms & genetic_terms:
                    boost += 0.10

            # Boost literature results generally (evidence-based)
            if result.collection == "sc_literature":
                boost += 0.03

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
        results: List[SCSearchResult],
    ) -> List[SCSearchResult]:
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
        results: List[SCSearchResult],
        workflow: Optional[SCWorkflowType] = None,
        patient_context: Optional[dict] = None,
    ) -> SCResponse:
        """Use LLM to synthesise search results into a single-cell response.

        Builds a structured prompt with retrieved evidence, sample context,
        conversation history, and workflow-specific instructions. Generates
        a grounded answer via the configured LLM.

        Args:
            question: Original single-cell analysis question.
            results: Ranked search results for context.
            workflow: Optional workflow for instruction tuning.
            patient_context: Optional sample/patient-specific data dict.

        Returns:
            SCResponse with synthesised answer and metadata.
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
            f"Please provide a comprehensive, evidence-based single-cell genomics "
            f"analysis grounded in the retrieved evidence above. "
            f"Follow the system prompt instructions for Cell Ontology citation format, "
            f"severity badges, cell type annotation scoring, TME classification, "
            f"and structured output sections.\n\n"
            f"Cite sources using clickable markdown links where PMIDs are available: "
            f"[PMID:12345678](https://pubmed.ncbi.nlm.nih.gov/12345678/). "
            f"For clinical trials, use [NCT01234567](https://clinicaltrials.gov/study/NCT01234567). "
            f"For Cell Ontology references, use [CL:0000084](https://www.ebi.ac.uk/ols4/ontologies/cl/classes/CL_0000084). "
            f"For collection-sourced evidence, use [Collection:record-id]. "
            f"Prioritise [high relevance] citations and atlas references."
        )

        answer = self.llm.generate(
            prompt=prompt,
            system_prompt=SC_SYSTEM_PROMPT,
            max_tokens=2048,
            temperature=0.7,
        )

        return SCResponse(
            question=question,
            answer=answer,
            results=results,
            workflow=workflow,
        )

    def _build_context(
        self,
        results: List[SCSearchResult],
        patient_context: Optional[dict] = None,
    ) -> str:
        """Build context string from search results for LLM prompt.

        Organises results by collection, formatting each with its
        citation reference, relevance tag, score, and text excerpt.

        Args:
            results: Ranked search results to format.
            patient_context: Optional sample context (used for additional
                context augmentation).

        Returns:
            Formatted evidence context string for the LLM prompt.
        """
        if not results:
            return "No evidence found in the knowledge base."

        # Group results by collection
        by_collection: Dict[str, List[SCSearchResult]] = {}
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

    def _format_citation_link(self, result: SCSearchResult) -> str:
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

        # Cell Ontology ID
        cl_id = result.metadata.get("cell_ontology_id", "")
        if cl_id:
            cl_class = cl_id.replace(":", "_")
            return (
                f"[{label}:{cl_id}]"
                f"(https://www.ebi.ac.uk/ols4/ontologies/cl/classes/{cl_class})"
            )

        return f"[{label}:{record_id}]"

    def _format_patient_context(self, patient_context: Optional[dict]) -> str:
        """Format sample/patient context for LLM prompt injection.

        Used for sample-specific analysis support.

        Args:
            patient_context: Optional sample data dict with keys like
                tissue, technology, cell_count, condition, treatment,
                timepoint, markers_of_interest, genomic_data, species,
                patient_id, sample_id, batch, quality_metrics.

        Returns:
            Formatted sample context section or empty string.
        """
        if not patient_context:
            return ""

        lines = ["### Sample / Patient Context\n"]

        field_labels = {
            "tissue": "Tissue / Organ",
            "condition": "Condition / Disease",
            "technology": "Sequencing Technology",
            "cell_count": "Cell Count",
            "species": "Species",
            "patient_id": "Patient ID",
            "sample_id": "Sample ID",
            "treatment": "Treatment",
            "timepoint": "Timepoint",
            "batch": "Batch",
            "markers_of_interest": "Markers of Interest",
            "genomic_data": "Genomic Data",
            "quality_metrics": "Quality Metrics",
            "clustering_resolution": "Clustering Resolution",
            "n_clusters": "Number of Clusters",
            "integration_method": "Integration Method",
            "reference_atlas": "Reference Atlas",
            "age": "Patient Age",
            "sex": "Patient Sex",
            "prior_therapies": "Prior Treatments",
            "comorbidities": "Comorbidities",
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
        workflow: Optional[SCWorkflowType],
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
            SCWorkflowType.CELL_TYPE_ANNOTATION: (
                "### Workflow: Cell Type Annotation\n"
                "Focus on: multi-strategy annotation consensus (reference-based mapping "
                "via Azimuth/CellTypist/scArches, marker-based scoring via CellMarker/"
                "PanglaoDB, LLM-assisted annotation), Cell Ontology ID assignment, "
                "canonical marker gene validation, cluster purity assessment, novel or "
                "ambiguous cell state identification, doublet detection impact, and "
                "annotation confidence scoring. Reference Human Cell Atlas datasets "
                "and Tabula Sapiens for normal tissue baselines.\n\n"
            ),
            SCWorkflowType.TME_CLASSIFICATION: (
                "### Workflow: Tumor Microenvironment Classification\n"
                "Focus on: four-category TME classification (hot/immune-inflamed, "
                "cold/immune-desert, excluded, immunosuppressive), immune cell "
                "composition quantification (CD8+ T cells, Tregs, TAMs, MDSCs, NK cells, "
                "DCs), spatial distribution of immune infiltrates (core vs. margin vs. "
                "stroma), gene signature scoring (IFN-gamma, TGF-beta, cytolytic "
                "activity), PD-L1 expression by cell type, treatment implications "
                "(checkpoint inhibitor eligibility, combination strategies), and "
                "comparison to published TME subtypes (Bagaev et al. 2021).\n\n"
            ),
            SCWorkflowType.DRUG_RESPONSE: (
                "### Workflow: Drug Response Prediction\n"
                "Focus on: cell-type-specific drug sensitivity scores from scRNA-seq, "
                "resistance mechanism identification (target mutation, pathway bypass, "
                "lineage switch, microenvironment-mediated), DepMap genetic dependency "
                "cross-reference, pharmacogenomic signature scoring, resistant subclone "
                "detection and frequency estimation, actionable target nomination, "
                "and combination strategy rationale. Distinguish computational "
                "predictions from experimentally validated findings.\n\n"
            ),
            SCWorkflowType.SUBCLONAL_ARCHITECTURE: (
                "### Workflow: Subclonal Architecture Detection\n"
                "Focus on: CNV inference from scRNA-seq (inferCNV, CopyKAT), clonal "
                "hierarchy reconstruction, founder vs. subclonal mutation assignment, "
                "clonal frequency estimation, phylogenetic tree inference, integration "
                "with bulk WGS/WES data, therapy-relevant subclone identification, "
                "cancer stem cell population characterization, and clonal dynamics "
                "between timepoints. Reference TCGA and PCAWG for known CNV patterns.\n\n"
            ),
            SCWorkflowType.SPATIAL_ANALYSIS: (
                "### Workflow: Spatial Transcriptomics Analysis\n"
                "Focus on: technology-appropriate analysis (Visium: spot deconvolution "
                "via cell2location/Tangram/RCTD; MERFISH/Slide-seq: single-cell "
                "resolution), spatial domain identification (SpatialDE, BayesSpace, "
                "CellCharter), spatial variable gene detection (Moran's I, SpatialDE, "
                "SPARK), cell-cell proximity analysis (Squidpy), niche identification "
                "and composition, ligand-receptor interactions in spatial context "
                "(Commot, LIANA+), and tissue architecture characterization. Note "
                "technology-specific resolution limits.\n\n"
            ),
            SCWorkflowType.TRAJECTORY_INFERENCE: (
                "### Workflow: Trajectory & Lineage Inference\n"
                "Focus on: trajectory inference method selection (Monocle3, PAGA, "
                "CellRank, Slingshot), RNA velocity analysis (scVelo: stochastic, "
                "dynamical, or steady-state model), pseudotime ordering, cell fate "
                "bifurcation identification, driver gene detection along trajectories, "
                "CytoTRACE stemness scoring, lineage tracing data integration, and "
                "benchmarking against Saelens et al. 2019 framework. Specify topology "
                "assumptions (linear, bifurcating, tree, graph).\n\n"
            ),
            SCWorkflowType.LIGAND_RECEPTOR: (
                "### Workflow: Ligand-Receptor Interaction Mapping\n"
                "Focus on: cell-cell communication inference (CellChat, LIANA+, "
                "NicheNet, Commot), ligand-receptor database selection (CellChatDB, "
                "OmniPath, CellPhoneDB), signaling pathway enrichment, sender-receiver "
                "cell type identification, spatial context integration, druggable "
                "interaction prioritization, autocrine vs. paracrine distinction, and "
                "comparison across conditions (tumor vs. normal, pre vs. post treatment). "
                "Note computational assumptions and validation requirements.\n\n"
            ),
            SCWorkflowType.BIOMARKER_DISCOVERY: (
                "### Workflow: Cell-Type-Specific Biomarker Discovery\n"
                "Focus on: differential expression analysis (Wilcoxon, MAST, "
                "pseudobulk DESeq2), gene module discovery (NMF, cNMF, Hotspot), "
                "surface marker identification for therapeutic targeting, prognostic "
                "signature construction and validation, AUROC calculation for marker "
                "specificity, cross-dataset validation, and integration with bulk "
                "RNA-seq cohorts (TCGA, GEO) for survival analysis. Distinguish "
                "discovery from validation findings.\n\n"
            ),
            SCWorkflowType.CART_TARGET_VALIDATION: (
                "### Workflow: CAR-T Target Validation\n"
                "Focus on: on-tumor expression quantification (percentage of tumor cells, "
                "expression level distribution, comparison across tumor subtypes), "
                "off-tumor expression profiling (Tabula Sapiens, Human Protein Atlas, "
                "GTEx for normal tissue expression), antigen escape risk assessment "
                "(heterogeneous expression, antigen-low subclones), co-expression with "
                "other surface antigens (dual-targeting opportunities), DepMap dependency "
                "for target essentiality, CITE-seq protein-level validation, and safety "
                "profile assessment (critical organ expression). Flag any expression in "
                "heart, lung, CNS, or liver as potential on-target/off-tumor toxicity.\n\n"
            ),
            SCWorkflowType.TREATMENT_MONITORING: (
                "### Workflow: Treatment Response Monitoring\n"
                "Focus on: longitudinal cell population dynamics (composition shifts "
                "between timepoints), clonal evolution under therapy pressure, immune "
                "cell state changes (exhaustion, activation, memory formation), MRD "
                "detection sensitivity, resistance clone emergence tracking, TME "
                "remodeling assessment, treatment-specific signatures (checkpoint "
                "inhibitor: T cell reinvigoration; CAR-T: expansion/contraction; "
                "chemotherapy: sensitive clone depletion), and response prediction "
                "biomarker identification.\n\n"
            ),
        }

        return instructions.get(workflow, "")

    # ==================================================================
    # CITATIONS & CONFIDENCE
    # ==================================================================

    def _extract_citations(
        self,
        results: List[SCSearchResult],
    ) -> List[dict]:
        """Extract and format citations from search results.

        Generates a structured citation list from all results, including
        PMID links, NCT links, DOI links, and Cell Ontology references.

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

            # Cell Ontology identifiers
            cl_id = result.metadata.get("cell_ontology_id", "")
            if cl_id and not cite["url"]:
                cl_class = cl_id.replace(":", "_")
                cite["url"] = f"https://www.ebi.ac.uk/ols4/ontologies/cl/classes/{cl_class}"
                cite["id"] = cl_id

            # Deduplicate
            dedup_key = cite["id"] or f"{cite['source']}:{result.record_id}"
            if dedup_key not in seen:
                seen.add(dedup_key)
                citations.append(cite)

        return citations

    def _score_confidence(
        self,
        results: List[SCSearchResult],
    ) -> float:
        """Score overall confidence based on result quality.

        Confidence is based on:
        - Number of high-relevance results
        - Collection diversity
        - Average similarity score
        - Presence of literature evidence

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

        # Factor 4: Literature evidence present (0-0.2)
        has_literature = any(
            r.collection == "sc_literature" for r in results
        )
        literature_score = 0.2 if has_literature else 0.0

        confidence = relevance_score + diversity_score + quality_score + literature_score
        return round(min(confidence, 1.0), 3)

    # ==================================================================
    # ENTITY & CONDITION SEARCH
    # ==================================================================

    def find_related(
        self,
        entity: str,
        entity_type: str = "cell_type",
        top_k: int = 5,
    ) -> List[SCSearchResult]:
        """Find related entities across collections.

        Searches all collections for evidence related to a single-cell
        entity (cell type, condition, biomarker, method). Useful
        for building entity profiles and cross-referencing.

        Args:
            entity: Entity name (e.g. 'CD8+ T cell', 'NSCLC', 'RNA velocity').
            entity_type: Entity category for targeted search:
                'cell_type', 'condition', 'biomarker', 'method', 'gene'.
            top_k: Maximum results per collection.

        Returns:
            List of SCSearchResult from all relevant collections.
        """
        type_collection_map = {
            "cell_type": [
                "sc_cell_types", "sc_markers", "sc_tme",
                "sc_spatial", "sc_trajectories", "sc_literature",
            ],
            "condition": [
                "sc_literature", "sc_clinical", "sc_tme",
                "sc_drug_response", "sc_cell_types", "sc_datasets",
            ],
            "biomarker": [
                "sc_markers", "sc_cell_types", "sc_clinical",
                "sc_literature", "sc_drug_response", "sc_tme",
            ],
            "method": [
                "sc_methods", "sc_literature", "sc_datasets",
                "sc_trajectories", "sc_spatial",
            ],
            "gene": [
                "genomic_evidence", "sc_markers", "sc_drug_response",
                "sc_pathways", "sc_literature", "sc_cell_types",
            ],
        }

        collections = type_collection_map.get(entity_type, get_all_collection_names())
        return self.search(entity, collections=collections, top_k=top_k)

    def search_cell_type(
        self,
        cell_type: str,
        tissue: Optional[str] = None,
    ) -> List[SCSearchResult]:
        """Search for cell type annotation evidence.

        Targeted search for cell type markers, references, and ontology data.

        Args:
            cell_type: Cell type name (e.g. 'CD8+ T cell', 'macrophage').
            tissue: Optional tissue context for tissue-specific markers.

        Returns:
            List of SCSearchResult from cell type and marker collections.
        """
        query = f"cell type annotation markers {cell_type}"
        if tissue:
            query = f"{tissue} {query}"

        return self.search(
            query,
            collections=["sc_cell_types", "sc_markers", "sc_literature"],
            top_k=10,
        )

    def search_spatial(
        self,
        tissue: str,
        features: str,
        top_k: int = 10,
    ) -> List[SCSearchResult]:
        """Search spatial transcriptomics evidence.

        Targeted search combining tissue and spatial features for spatial
        analysis support.

        Args:
            tissue: Tissue type (e.g. 'breast tumor', 'lung').
            features: Spatial feature description.
            top_k: Maximum results.

        Returns:
            List of SCSearchResult from spatial and related collections.
        """
        query = f"{tissue} {features}"
        collections = ["sc_spatial", "sc_tme", "sc_literature"]
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
        Useful when starting a new analysis session or switching topics.

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
        workflow: Optional[SCWorkflowType] = None,
    ) -> Dict[str, float]:
        """Compute collection weights with optional workflow boosting.

        When a workflow is specified, applies boost multipliers from
        WORKFLOW_COLLECTION_BOOST on top of the base weights from
        settings. Weights are then renormalized to sum to ~1.0.

        Args:
            workflow: Optional SCWorkflowType for boosting.

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
    # WORKFLOW METHODS (called from api/routes/sc_clinical.py)
    # ==================================================================

    def annotate_cell_types(
        self,
        gene_expression: Optional[dict] = None,
        marker_genes: Optional[List[str]] = None,
        tissue_type: Optional[str] = None,
        strategy: str = "marker_based",
        top_n: int = 5,
    ) -> dict:
        """Annotate cell types using RAG search augmented with knowledge base.

        Args:
            gene_expression: Gene expression profile (gene -> value).
            marker_genes: Detected marker genes for annotation.
            tissue_type: Source tissue (e.g., lung, breast, PBMC).
            strategy: Annotation strategy: reference_based | marker_based | llm_based.
            top_n: Max cell types to return.

        Returns:
            Dict with keys matching AnnotateResponse schema.
        """
        # RAG search for context
        query_parts = ["cell type annotation"]
        if marker_genes:
            query_parts.append(" ".join(marker_genes[:10]))
        if tissue_type:
            query_parts.append(tissue_type)
        query = " ".join(query_parts)

        evidence = []
        try:
            results = self.search(query, top_k=5)
            evidence = [
                {"text": r.text[:300], "collection": r.collection, "score": round(r.score, 3)}
                for r in results
            ]
        except Exception as exc:
            logger.warning("RAG search for annotate_cell_types failed: %s", exc)

        # Knowledge base matching
        cell_types = []
        marker_evidence = []
        if marker_genes:
            input_markers = set(g.upper() for g in marker_genes)
            for ct_name, ct_data in CELL_TYPE_ATLAS.items():
                ct_markers = set(m.upper() for m in ct_data.get("markers", []))
                overlap = input_markers & ct_markers
                if overlap:
                    score = len(overlap) / max(len(ct_markers), 1)
                    cell_types.append({
                        "cell_type": ct_name.replace("_", " "),
                        "confidence": round(score, 2),
                        "markers": list(overlap),
                        "cell_ontology_id": ct_data.get("cell_ontology_id", ""),
                        "source": "knowledge_base",
                    })
                    marker_evidence.append({
                        "cell_type": ct_name,
                        "matched_markers": list(overlap),
                        "total_markers": len(ct_markers),
                    })
            cell_types.sort(key=lambda x: x["confidence"], reverse=True)
            cell_types = cell_types[:top_n]

        # LLM synthesis if available
        recommendations = []
        if self.llm and (marker_genes or gene_expression):
            try:
                prompt = (
                    f"Annotate cell types. Markers: {marker_genes}, "
                    f"Tissue: {tissue_type or 'unknown'}, Strategy: {strategy}.\n"
                    f"Evidence: {evidence[:3]}\n"
                    f"Provide top cell type annotations with reasoning."
                )
                answer = self.llm.generate(prompt, system_prompt=SC_SYSTEM_PROMPT, max_tokens=1024)
                recommendations.append(answer)
            except Exception:
                pass

        if not cell_types:
            cell_types = [{"cell_type": "unknown", "confidence": 0.0, "markers": marker_genes or []}]

        return {
            "cell_types": cell_types,
            "cluster_id": None,
            "tissue_type": tissue_type,
            "strategy_used": strategy,
            "num_cells": None,
            "confidence": cell_types[0]["confidence"] if cell_types else 0.0,
            "marker_evidence": marker_evidence,
            "recommendations": recommendations,
        }

    def profile_tme(
        self,
        cell_type_proportions: Optional[dict] = None,
        immune_markers: Optional[dict] = None,
        tumor_type: Optional[str] = None,
    ) -> dict:
        """Profile the tumor microenvironment using RAG and knowledge base.

        Args:
            cell_type_proportions: Cell type -> proportion mapping.
            immune_markers: Immune marker expression levels.
            tumor_type: Cancer type (e.g., NSCLC, melanoma).

        Returns:
            Dict with keys matching TMEProfileResponse schema.
        """
        query = f"tumor microenvironment classification {tumor_type or ''} immune infiltration"
        evidence = []
        try:
            results = self.search(query, top_k=5)
            evidence = [
                {"text": r.text[:300], "collection": r.collection, "score": round(r.score, 3)}
                for r in results
            ]
        except Exception as exc:
            logger.warning("RAG search for profile_tme failed: %s", exc)

        # Knowledge base: match tumor type against CANCER_TME_ATLAS
        tme_class = "unknown"
        therapy_prediction = {}
        kb_evidence = []
        if tumor_type:
            tumor_key = tumor_type.lower().replace(" ", "_")
            for atlas_key, atlas_data in CANCER_TME_ATLAS.items():
                if tumor_key in atlas_key.lower():
                    tme_class = atlas_data.get("dominant_tme", "unknown")
                    therapy_prediction = atlas_data.get("therapy_implications", {})
                    kb_evidence.append({"source": "cancer_tme_atlas", "entry": atlas_key})
                    break

        for tme_key, tme_data in TME_PROFILES.items():
            kb_evidence.append({
                "source": "tme_profiles",
                "class": tme_key,
                "description": tme_data.get("description", ""),
            })

        # Compute immune/stromal scores from proportions
        immune_score = 0.0
        stromal_score = 0.0
        if cell_type_proportions:
            immune_types = {"t_cell", "cd8_t", "cd4_t", "nk_cell", "b_cell", "macrophage",
                            "dendritic_cell", "treg", "mast_cell", "neutrophil"}
            stromal_types = {"fibroblast", "endothelial", "pericyte", "caf"}
            for ct, prop in cell_type_proportions.items():
                ct_lower = ct.lower().replace(" ", "_")
                if ct_lower in immune_types:
                    immune_score += prop
                elif ct_lower in stromal_types:
                    stromal_score += prop

        # LLM synthesis
        recommendations = []
        if self.llm:
            try:
                prompt = (
                    f"Classify the tumor microenvironment.\n"
                    f"Cell proportions: {cell_type_proportions}\n"
                    f"Tumor type: {tumor_type or 'unknown'}\n"
                    f"Evidence: {evidence[:3]}\n"
                    f"Classify as hot/cold/excluded/immunosuppressive with therapy recommendations."
                )
                answer = self.llm.generate(prompt, system_prompt=SC_SYSTEM_PROMPT, max_tokens=1024)
                recommendations.append(answer)
            except Exception:
                pass

        return {
            "tme_classification": tme_class,
            "immune_score": round(immune_score, 3),
            "stromal_score": round(stromal_score, 3),
            "cell_composition": cell_type_proportions or {},
            "immune_phenotype": tme_class,
            "checkpoint_expression": {},
            "therapy_prediction": therapy_prediction,
            "recommendations": recommendations,
            "evidence": (kb_evidence + evidence)[:10],
            "guidelines_cited": [],
        }

    def predict_drug_response(
        self,
        gene_expression: Optional[dict] = None,
        cell_type: Optional[str] = None,
        drug_name: Optional[str] = None,
        tumor_type: Optional[str] = None,
    ) -> dict:
        """Predict drug response using RAG and knowledge base.

        Args:
            gene_expression: Gene expression profile.
            cell_type: Target cell type.
            drug_name: Drug/compound to evaluate.
            tumor_type: Cancer type.

        Returns:
            Dict with keys matching DrugResponseResponse schema.
        """
        query = f"drug response {drug_name or ''} {cell_type or ''} {tumor_type or ''} sensitivity resistance"
        evidence = []
        try:
            results = self.search(query, top_k=5)
            evidence = [
                {"text": r.text[:300], "collection": r.collection, "score": round(r.score, 3)}
                for r in results
            ]
        except Exception as exc:
            logger.warning("RAG search for predict_drug_response failed: %s", exc)

        # Knowledge base matching
        predictions = []
        biomarkers = []
        resistance_mechanisms = []
        if drug_name:
            drug_key = drug_name.lower().replace(" ", "_")
            for db_key, db_data in DRUG_SENSITIVITY_DATABASE.items():
                if drug_key in db_key.lower() or db_key.lower() in drug_key:
                    predictions.append({
                        "drug": db_key,
                        "drug_class": db_data.get("drug_class", ""),
                        "targets": db_data.get("targets", []),
                        "sensitive_cell_types": db_data.get("sensitive_cell_types", []),
                        "source": "knowledge_base",
                    })
                    biomarkers.extend(db_data.get("biomarkers", []))
                    resistance_mechanisms.extend(db_data.get("resistance_mechanisms", []))

        # LLM synthesis
        recommendations = []
        if self.llm:
            try:
                prompt = (
                    f"Predict drug response.\nDrug: {drug_name or 'not specified'}\n"
                    f"Cell type: {cell_type or 'not specified'}\n"
                    f"Tumor type: {tumor_type or 'unknown'}\n"
                    f"Evidence: {evidence[:3]}\n"
                    f"Provide sensitivity prediction and resistance mechanisms."
                )
                answer = self.llm.generate(prompt, system_prompt=SC_SYSTEM_PROMPT, max_tokens=1024)
                recommendations.append(answer)
            except Exception:
                pass

        return {
            "predictions": predictions[:5],
            "sensitivity_score": 0.0,
            "resistance_mechanisms": list(set(resistance_mechanisms))[:5],
            "biomarkers": list(set(biomarkers))[:10],
            "recommendations": recommendations,
            "evidence": evidence[:5],
            "guidelines_cited": [],
        }

    def analyze_subclones(
        self,
        cnv_profile: Optional[dict] = None,
        mutation_data: Optional[dict] = None,
        target_antigen: Optional[str] = None,
    ) -> dict:
        """Analyze subclonal architecture using RAG and knowledge base.

        Args:
            cnv_profile: Copy number variation data.
            mutation_data: Somatic mutation data.
            target_antigen: CAR-T/therapy target antigen.

        Returns:
            Dict with keys matching SubclonalResponse schema.
        """
        query = f"subclonal architecture CNV inference clonal heterogeneity {target_antigen or ''}"
        evidence = []
        try:
            results = self.search(query, top_k=5)
            evidence = [
                {"text": r.text[:300], "collection": r.collection, "score": round(r.score, 3)}
                for r in results
            ]
        except Exception as exc:
            logger.warning("RAG search for analyze_subclones failed: %s", exc)

        # Estimate clone count from mutation data
        clones = []
        num_clones = 0
        if mutation_data:
            num_clones = min(len(mutation_data), 5)
            for i, (gene, info) in enumerate(list(mutation_data.items())[:5]):
                clones.append({
                    "clone_id": f"clone_{i+1}",
                    "driver_gene": gene,
                    "frequency": info if isinstance(info, (int, float)) else 0.0,
                    "source": "input_data",
                })

        # LLM synthesis
        recommendations = []
        if self.llm:
            try:
                prompt = (
                    f"Analyze subclonal architecture.\n"
                    f"CNV profile: {cnv_profile}\nMutation data: {mutation_data}\n"
                    f"Target antigen: {target_antigen or 'N/A'}\n"
                    f"Evidence: {evidence[:3]}\n"
                    f"Provide clonal hierarchy and escape risk assessment."
                )
                answer = self.llm.generate(prompt, system_prompt=SC_SYSTEM_PROMPT, max_tokens=1024)
                recommendations.append(answer)
            except Exception:
                pass

        return {
            "clones": clones,
            "num_clones": num_clones,
            "dominant_clone": clones[0]["clone_id"] if clones else None,
            "heterogeneity_index": min(num_clones / 5.0, 1.0) if num_clones else 0.0,
            "escape_risk": "high" if num_clones >= 3 else ("moderate" if num_clones >= 1 else "unknown"),
            "recommendations": recommendations,
            "evidence": evidence[:5],
        }

    def map_spatial_niches(
        self,
        spatial_coordinates: Optional[dict] = None,
        cell_types: Optional[dict] = None,
        platform: str = "visium",
    ) -> dict:
        """Map spatial niches from spatial transcriptomics data.

        Args:
            spatial_coordinates: Cell coordinate data.
            cell_types: Cell type assignments per spot/cell.
            platform: Spatial platform (visium, merfish, cosmx, etc.).

        Returns:
            Dict with keys matching SpatialNicheResponse schema.
        """
        query = f"spatial niche mapping {platform} cell neighborhoods tissue architecture"
        evidence = []
        try:
            results = self.search(query, top_k=5)
            evidence = [
                {"text": r.text[:300], "collection": r.collection, "score": round(r.score, 3)}
                for r in results
            ]
        except Exception as exc:
            logger.warning("RAG search for map_spatial_niches failed: %s", exc)

        # Knowledge base: get cell-cell interactions from LIGAND_RECEPTOR_PAIRS
        interactions = []
        if cell_types:
            input_types = set(ct.lower() for ct in (cell_types.values() if isinstance(cell_types, dict) else []))
            for lr_key, lr_data in LIGAND_RECEPTOR_PAIRS.items():
                source_types = set(ct.lower() for ct in lr_data.get("source_cell_types", []))
                target_types = set(ct.lower() for ct in lr_data.get("target_cell_types", []))
                if input_types & (source_types | target_types):
                    interactions.append({
                        "ligand": lr_data.get("ligand", lr_key),
                        "receptor": lr_data.get("receptor", ""),
                        "source": "knowledge_base",
                    })

        # LLM synthesis
        recommendations = []
        if self.llm:
            try:
                prompt = (
                    f"Map spatial niches.\nPlatform: {platform}\n"
                    f"Cell types: {cell_types}\n"
                    f"Evidence: {evidence[:3]}\n"
                    f"Identify spatial niches and cell-cell interactions."
                )
                answer = self.llm.generate(prompt, system_prompt=SC_SYSTEM_PROMPT, max_tokens=1024)
                recommendations.append(answer)
            except Exception:
                pass

        return {
            "niches": [],
            "num_niches": 0,
            "spatial_statistics": {},
            "cell_cell_interactions": interactions[:10],
            "platform": platform,
            "recommendations": recommendations,
            "evidence": evidence[:5],
        }

    def infer_trajectory(
        self,
        cell_types: Optional[List[str]] = None,
        root_cell_type: Optional[str] = None,
        method: str = "monocle3",
    ) -> dict:
        """Infer cell trajectories and pseudotime ordering.

        Args:
            cell_types: Cell types in trajectory.
            root_cell_type: Root/stem cell type.
            method: Trajectory inference method.

        Returns:
            Dict with keys matching TrajectoryResponse schema.
        """
        query = f"trajectory inference {method} pseudotime {root_cell_type or ''} differentiation lineage"
        evidence = []
        try:
            results = self.search(query, top_k=5)
            evidence = [
                {"text": r.text[:300], "collection": r.collection, "score": round(r.score, 3)}
                for r in results
            ]
        except Exception as exc:
            logger.warning("RAG search for infer_trajectory failed: %s", exc)

        # Knowledge base: provide driver genes from MARKER_GENE_DATABASE
        driver_genes = []
        if cell_types:
            for ct in cell_types:
                ct_key = ct.lower().replace(" ", "_")
                for mk_key, mk_data in MARKER_GENE_DATABASE.items():
                    if ct_key in mk_key.lower():
                        for gene in mk_data.get("markers", [])[:3]:
                            driver_genes.append({
                                "gene": gene,
                                "cell_type": mk_key,
                                "source": "knowledge_base",
                            })

        # LLM synthesis
        recommendations = []
        if self.llm:
            try:
                prompt = (
                    f"Infer cell trajectory.\nMethod: {method}\n"
                    f"Cell types: {cell_types}\nRoot: {root_cell_type or 'N/A'}\n"
                    f"Evidence: {evidence[:3]}\n"
                    f"Provide trajectory analysis with branch points and driver genes."
                )
                answer = self.llm.generate(prompt, system_prompt=SC_SYSTEM_PROMPT, max_tokens=1024)
                recommendations.append(answer)
            except Exception:
                pass

        return {
            "trajectory": {},
            "branch_points": [],
            "driver_genes": driver_genes[:10],
            "pseudotime_range": [],
            "method_used": method,
            "recommendations": recommendations,
            "evidence": evidence[:5],
        }

    def analyze_ligand_receptor(
        self,
        cell_types: Optional[List[str]] = None,
        source_cell_type: Optional[str] = None,
        target_cell_type: Optional[str] = None,
        database: str = "cellphonedb",
    ) -> dict:
        """Analyze ligand-receptor interactions between cell types.

        Args:
            cell_types: Cell types to analyze.
            source_cell_type: Ligand source cell type.
            target_cell_type: Receptor target cell type.
            database: Database: cellphonedb | nichenet | cellchat | celltalkdb.

        Returns:
            Dict with keys matching LigandReceptorResponse schema.
        """
        query = f"ligand receptor interaction {database} {source_cell_type or ''} {target_cell_type or ''} cell communication"
        evidence = []
        try:
            results = self.search(query, top_k=5)
            evidence = [
                {"text": r.text[:300], "collection": r.collection, "score": round(r.score, 3)}
                for r in results
            ]
        except Exception as exc:
            logger.warning("RAG search for analyze_ligand_receptor failed: %s", exc)

        # Knowledge base matching
        interactions = []
        query_types = set()
        if cell_types:
            query_types.update(ct.lower() for ct in cell_types)
        if source_cell_type:
            query_types.add(source_cell_type.lower())
        if target_cell_type:
            query_types.add(target_cell_type.lower())

        if query_types:
            for lr_key, lr_data in LIGAND_RECEPTOR_PAIRS.items():
                source_types = set(ct.lower() for ct in lr_data.get("source_cell_types", []))
                target_types = set(ct.lower() for ct in lr_data.get("target_cell_types", []))
                if query_types & (source_types | target_types):
                    interactions.append({
                        "ligand": lr_data.get("ligand", lr_key),
                        "receptor": lr_data.get("receptor", ""),
                        "pathway": lr_data.get("pathway", ""),
                        "source_cell_types": lr_data.get("source_cell_types", []),
                        "target_cell_types": lr_data.get("target_cell_types", []),
                        "source": "knowledge_base",
                    })

        # LLM synthesis
        recommendations = []
        if self.llm:
            try:
                prompt = (
                    f"Analyze ligand-receptor interactions.\n"
                    f"Database: {database}\nCell types: {cell_types}\n"
                    f"Source: {source_cell_type}, Target: {target_cell_type}\n"
                    f"Evidence: {evidence[:3]}\n"
                    f"Provide interaction network analysis and therapeutic targets."
                )
                answer = self.llm.generate(prompt, system_prompt=SC_SYSTEM_PROMPT, max_tokens=1024)
                recommendations.append(answer)
            except Exception:
                pass

        return {
            "interactions": interactions[:15],
            "num_significant": len(interactions),
            "pathways_enriched": [],
            "network_summary": {},
            "recommendations": recommendations,
            "evidence": evidence[:5],
        }

    def discover_biomarkers(
        self,
        differential_expression: Optional[dict] = None,
        cell_type: Optional[str] = None,
        min_log2fc: float = 1.0,
        max_pval: float = 0.05,
    ) -> dict:
        """Discover biomarker candidates using RAG and knowledge base.

        Args:
            differential_expression: DE results (gene -> stats).
            cell_type: Cell type of interest.
            min_log2fc: Minimum log2 fold change threshold.
            max_pval: Maximum adjusted p-value threshold.

        Returns:
            Dict with keys matching BiomarkerResponse schema.
        """
        query = f"biomarker discovery {cell_type or ''} differential expression marker specificity"
        evidence = []
        try:
            results = self.search(query, top_k=5)
            evidence = [
                {"text": r.text[:300], "collection": r.collection, "score": round(r.score, 3)}
                for r in results
            ]
        except Exception as exc:
            logger.warning("RAG search for discover_biomarkers failed: %s", exc)

        # Knowledge base matching
        biomarkers = []
        top_markers = []
        if cell_type:
            ct_key = cell_type.lower().replace(" ", "_")
            for mk_key, mk_data in MARKER_GENE_DATABASE.items():
                if ct_key in mk_key.lower():
                    for gene in mk_data.get("markers", []):
                        biomarkers.append({
                            "gene": gene,
                            "cell_type": mk_key,
                            "source": "marker_gene_database",
                        })
                        top_markers.append(gene)
            # Check immune signatures
            for sig_key, sig_data in IMMUNE_SIGNATURES.items():
                sig_genes = sig_data.get("genes", [])
                if any(ct_key in g.lower() for g in sig_data.get("cell_types", [])):
                    for gene in sig_genes[:3]:
                        biomarkers.append({
                            "gene": gene,
                            "signature": sig_key,
                            "source": "immune_signatures",
                        })
                        top_markers.append(gene)

        # Filter DE results if provided
        if differential_expression:
            for gene, stats in differential_expression.items():
                if isinstance(stats, dict):
                    log2fc = abs(stats.get("log2fc", 0))
                    pval = stats.get("pval_adj", stats.get("pval", 1.0))
                    if log2fc >= min_log2fc and pval <= max_pval:
                        biomarkers.append({
                            "gene": gene,
                            "log2fc": log2fc,
                            "pval_adj": pval,
                            "source": "differential_expression",
                        })
                        top_markers.append(gene)

        # LLM synthesis
        recommendations = []
        if self.llm:
            try:
                prompt = (
                    f"Discover biomarkers.\nCell type: {cell_type or 'N/A'}\n"
                    f"Thresholds: log2FC >= {min_log2fc}, p-adj <= {max_pval}\n"
                    f"Evidence: {evidence[:3]}\n"
                    f"Identify top biomarker candidates with validation suggestions."
                )
                answer = self.llm.generate(prompt, system_prompt=SC_SYSTEM_PROMPT, max_tokens=1024)
                recommendations.append(answer)
            except Exception:
                pass

        return {
            "biomarkers": biomarkers[:15],
            "num_candidates": len(biomarkers),
            "top_markers": list(dict.fromkeys(top_markers))[:10],
            "validation_suggestions": [],
            "recommendations": recommendations,
            "evidence": evidence[:5],
        }

    def validate_cart_target(
        self,
        target_gene: str = "",
        tumor_type: Optional[str] = None,
        expression_data: Optional[dict] = None,
        tme_data: Optional[dict] = None,
    ) -> dict:
        """Validate a CAR-T therapy target using RAG and knowledge base.

        Args:
            target_gene: Target antigen gene symbol.
            tumor_type: Cancer type.
            expression_data: Expression across cell types.
            tme_data: TME composition data.

        Returns:
            Dict with keys matching CARTValidateResponse schema.
        """
        query = f"CAR-T target validation {target_gene} {tumor_type or ''} on-tumor off-tumor expression"
        evidence = []
        try:
            results = self.search(query, top_k=5)
            evidence = [
                {"text": r.text[:300], "collection": r.collection, "score": round(r.score, 3)}
                for r in results
            ]
        except Exception as exc:
            logger.warning("RAG search for validate_cart_target failed: %s", exc)

        # Knowledge base: check CELL_TYPE_ATLAS for off-tumor expression
        kb_evidence = []
        off_tumor_risk = {}
        for ct_name, ct_data in CELL_TYPE_ATLAS.items():
            markers = [m.upper() for m in ct_data.get("markers", [])]
            if target_gene.upper() in markers:
                kb_evidence.append({
                    "cell_type": ct_name,
                    "expresses_target": True,
                    "source": "cell_type_atlas",
                })
                tissues = ct_data.get("tissues", [])
                if tissues:
                    off_tumor_risk[ct_name] = tissues

        # Check CANCER_TME_ATLAS for tumor context
        if tumor_type:
            tumor_key = tumor_type.lower().replace(" ", "_")
            for atlas_key, atlas_data in CANCER_TME_ATLAS.items():
                if tumor_key in atlas_key.lower():
                    kb_evidence.append({"source": "cancer_tme_atlas", "entry": atlas_key})
                    break

        # LLM synthesis
        recommendations = []
        if self.llm:
            try:
                prompt = (
                    f"Evaluate CAR-T target {target_gene} for "
                    f"{tumor_type or 'solid tumor'}.\n"
                    f"Off-tumor risk: {off_tumor_risk}\n"
                    f"Evidence: {evidence[:3]}\n"
                    f"Assess on-tumor expression, off-tumor toxicity, TME compatibility, "
                    f"and antigen escape risk."
                )
                answer = self.llm.generate(prompt, system_prompt=SC_SYSTEM_PROMPT, max_tokens=1024)
                recommendations.append(answer)
            except Exception:
                pass

        return {
            "target_gene": target_gene,
            "on_tumor_pct": 0.0,
            "off_tumor_risk": off_tumor_risk,
            "tme_compatibility": 0.0,
            "escape_risk": "unknown",
            "therapeutic_index": 0.0,
            "subclonal_heterogeneity": None,
            "recommendations": recommendations,
            "evidence": (kb_evidence + evidence)[:10],
            "guidelines_cited": [],
        }

    def monitor_treatment(
        self,
        timepoints: Optional[List[dict]] = None,
        treatment: Optional[str] = None,
        baseline_composition: Optional[dict] = None,
        current_composition: Optional[dict] = None,
    ) -> dict:
        """Monitor treatment response using longitudinal single-cell data.

        Args:
            timepoints: Longitudinal data per timepoint.
            treatment: Treatment regimen.
            baseline_composition: Pre-treatment cell composition.
            current_composition: Current cell composition.

        Returns:
            Dict with keys matching TreatmentMonitorResponse schema.
        """
        query = f"treatment monitoring {treatment or ''} longitudinal clonal dynamics resistance"
        evidence = []
        try:
            results = self.search(query, top_k=5)
            evidence = [
                {"text": r.text[:300], "collection": r.collection, "score": round(r.score, 3)}
                for r in results
            ]
        except Exception as exc:
            logger.warning("RAG search for monitor_treatment failed: %s", exc)

        # Knowledge base: resistance mechanisms from DRUG_SENSITIVITY_DATABASE
        resistance_indicators = []
        if treatment:
            treatment_lower = treatment.lower()
            for drug_key, drug_data in DRUG_SENSITIVITY_DATABASE.items():
                if drug_key.lower() in treatment_lower or treatment_lower in drug_key.lower():
                    resistance_indicators.extend([
                        {"mechanism": m, "source": "drug_sensitivity_database"}
                        for m in drug_data.get("resistance_mechanisms", [])[:3]
                    ])

        # Compute composition changes
        immune_shift = {}
        if baseline_composition and current_composition:
            for ct in baseline_composition:
                baseline_val = baseline_composition.get(ct, 0)
                current_val = current_composition.get(ct, 0)
                if baseline_val > 0:
                    immune_shift[ct] = {
                        "baseline": baseline_val,
                        "current": current_val,
                        "change_pct": round((current_val - baseline_val) / baseline_val * 100, 1),
                    }

        # Assess response
        response_assessment = "insufficient_data"
        if immune_shift:
            increases = sum(1 for v in immune_shift.values() if v["change_pct"] > 10)
            decreases = sum(1 for v in immune_shift.values() if v["change_pct"] < -10)
            if increases > decreases:
                response_assessment = "responding"
            elif decreases > increases:
                response_assessment = "progressing"
            else:
                response_assessment = "stable"

        # LLM synthesis
        recommendations = []
        if self.llm:
            try:
                prompt = (
                    f"Monitor treatment response.\nTreatment: {treatment or 'N/A'}\n"
                    f"Composition changes: {immune_shift}\n"
                    f"Evidence: {evidence[:3]}\n"
                    f"Assess response, resistance indicators, and immune shifts."
                )
                answer = self.llm.generate(prompt, system_prompt=SC_SYSTEM_PROMPT, max_tokens=1024)
                recommendations.append(answer)
            except Exception:
                pass

        return {
            "response_assessment": response_assessment,
            "composition_changes": immune_shift,
            "resistance_indicators": resistance_indicators[:5],
            "clone_dynamics": [],
            "immune_shift": immune_shift,
            "recommendations": recommendations,
            "evidence": evidence[:5],
        }
