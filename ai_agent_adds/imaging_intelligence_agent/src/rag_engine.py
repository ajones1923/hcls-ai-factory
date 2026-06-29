"""Multi-collection RAG engine for Imaging Intelligence Agent.

Searches 10 imaging collections + 1 genomic collection, augments results
with domain knowledge, and synthesizes answers via LLM.
"""

import time
from typing import Any, Dict, Generator, List, Optional, Set

from loguru import logger

from src.models import (
    AgentResponse,
    ComparativeResult,
    CrossCollectionResult,
    SearchHit,
)
from src.knowledge import (
    IMAGING_ANATOMY,
    IMAGING_MODALITIES,
    IMAGING_PATHOLOGIES,
    get_anatomy_context,
    get_modality_context,
    get_pathology_context,
    get_comparison_context,
    resolve_comparison_entity,
)
from src.query_expansion import expand_query


# Collection configuration with weights and metadata
COLLECTION_CONFIG = {
    "imaging_literature": {"weight": 0.18, "label": "Literature", "has_modality": True, "year_field": "year"},
    "imaging_trials": {"weight": 0.12, "label": "Trial", "has_modality": True, "year_field": "start_year"},
    "imaging_findings": {"weight": 0.15, "label": "Finding", "has_modality": True, "year_field": None},
    "imaging_protocols": {"weight": 0.08, "label": "Protocol", "has_modality": True, "year_field": None},
    "imaging_devices": {"weight": 0.08, "label": "Device", "has_modality": True, "year_field": None},
    "imaging_anatomy": {"weight": 0.06, "label": "Anatomy", "has_modality": False, "year_field": None},
    "imaging_benchmarks": {"weight": 0.08, "label": "Benchmark", "has_modality": True, "year_field": None},
    "imaging_guidelines": {"weight": 0.10, "label": "Guideline", "has_modality": True, "year_field": "year"},
    "imaging_report_templates": {"weight": 0.05, "label": "ReportTemplate", "has_modality": True, "year_field": None},
    "imaging_datasets": {"weight": 0.06, "label": "Dataset", "has_modality": True, "year_field": None},
    "genomic_evidence": {"weight": 0.04, "label": "Genomic", "has_modality": False, "year_field": None},
}

SYSTEM_PROMPT = """You are an expert medical imaging intelligence assistant with deep knowledge in:

1. CT Analysis — head hemorrhage triage, chest lung nodule tracking, abdominal pathology
2. MRI Interpretation — brain MS lesions, tumor characterization, cardiac imaging
3. Chest X-ray — findings classification, pneumonia, cardiac silhouette, pneumothorax
4. Imaging AI Models — MONAI, VISTA-3D, nnU-Net, SwinUNETR, DenseNet architectures
5. Clinical Guidelines — ACR Appropriateness Criteria, Lung-RADS, BI-RADS, TI-RADS, LI-RADS
6. Imaging Protocols — CT/MRI acquisition parameters, contrast agents, radiation dose
7. FDA-Cleared Devices — 510(k), De Novo AI devices for medical imaging
8. Radiology Reporting — structured reporting, RadLex terminology, DICOM SR
9. Public Datasets — RSNA, TCIA, NIH, LIDC-IDRI, BraTS, CheXpert, MIMIC-CXR
10. Quantitative Imaging — volumetrics, RECIST criteria, volume doubling time
11. NVIDIA NIMs — VISTA-3D segmentation, MAISI synthetic CT, VILA-M3 VLM, Llama3 clinical reasoning

Always cite evidence from the provided context. Use clinical terminology appropriately.
When discussing AI models, mention their architecture, training data, and validation metrics.
For clinical findings, reference relevant classification systems (Lung-RADS, BI-RADS, etc.).

IMPORTANT: This system is for research purposes only. All outputs require clinician review.
Do not provide definitive clinical diagnoses."""


class ImagingRAGEngine:
    """Multi-collection RAG engine with knowledge augmentation."""

    def __init__(self, collection_manager, embedder, llm_client, nim_service_manager=None):
        self.collection_manager = collection_manager
        self.embedder = embedder
        self.llm_client = llm_client
        self.nim_manager = nim_service_manager
        self.system_prompt = SYSTEM_PROMPT

    def _embed_query(self, text: str) -> List[float]:
        return self.embedder.encode(text, normalize_embeddings=True).tolist()

    def _get_knowledge_context(self, query: str) -> str:
        """Check query for pathology/modality/anatomy mentions and return context."""
        query_lower = query.lower()
        contexts = []

        for key in IMAGING_PATHOLOGIES:
            if key.replace("_", " ") in query_lower:
                contexts.append(get_pathology_context(key))
                break

        for key in IMAGING_MODALITIES:
            if key in query_lower:
                contexts.append(get_modality_context(key))
                break

        for key in IMAGING_ANATOMY:
            if key in query_lower:
                contexts.append(get_anatomy_context(key))
                break

        return "\n\n".join(contexts) if contexts else ""

    def _is_comparative(self, query: str) -> bool:
        comparative_keywords = [" vs ", " versus ", " compared to ", " comparison ", " better than ",
                                " differ ", " difference ", " advantages ", " disadvantages "]
        return any(kw in query.lower() for kw in comparative_keywords)

    def _build_prompt(self, question: str, evidence: CrossCollectionResult, conversation_context: str = "") -> List[Dict]:
        """Build LLM prompt with system message + evidence + question."""
        evidence_text = ""
        for hit in evidence.hits[:20]:  # Limit to top 20 hits
            config = COLLECTION_CONFIG.get(hit.collection, {})
            label = config.get("label", hit.collection)
            evidence_text += f"\n[{label}] (score: {hit.score:.3f}) {hit.text}\n"

        knowledge = evidence.knowledge_context

        user_content = ""
        if knowledge:
            user_content += f"## Domain Knowledge\n{knowledge}\n\n"
        if evidence_text:
            user_content += f"## Retrieved Evidence ({evidence.hit_count} results from {evidence.total_collections_searched} collections)\n{evidence_text}\n\n"
        if conversation_context:
            user_content += f"## Conversation Context\n{conversation_context}\n\n"
        user_content += f"## Question\n{question}"

        return [
            {"role": "system", "content": self.system_prompt},
            {"role": "user", "content": user_content},
        ]

    def retrieve(
        self,
        query: str,
        top_k_per_collection: int = 5,
        collections_filter: Optional[List[str]] = None,
        year_min: Optional[int] = None,
        year_max: Optional[int] = None,
        modality_filter: Optional[str] = None,
        body_region_filter: Optional[str] = None,
    ) -> CrossCollectionResult:
        """Retrieve evidence from multiple collections."""
        start = time.time()

        # Expand query for better recall
        expanded_terms = expand_query(query)
        search_text = query
        if expanded_terms:
            search_text = f"{query} {' '.join(list(expanded_terms)[:10])}"

        query_embedding = self._embed_query(search_text)
        knowledge_context = self._get_knowledge_context(query)

        # Determine which collections to search
        collections = collections_filter or list(COLLECTION_CONFIG.keys())

        # Build per-collection filters
        filter_exprs = {}
        # Collections that have a body_region field
        _has_body_region = {"imaging_literature", "imaging_trials", "imaging_findings",
                           "imaging_protocols", "imaging_devices", "imaging_anatomy",
                           "imaging_benchmarks", "imaging_guidelines",
                           "imaging_report_templates", "imaging_datasets"}
        for coll in collections:
            config = COLLECTION_CONFIG.get(coll, {})
            filters = []
            if modality_filter and config.get("has_modality"):
                filters.append(f'modality == "{modality_filter}"')
            if body_region_filter and coll in _has_body_region:
                filters.append(f'body_region == "{body_region_filter}"')
            if year_min and config.get("year_field"):
                filters.append(f'{config["year_field"]} >= {int(year_min)}')
            if year_max and config.get("year_field"):
                filters.append(f'{config["year_field"]} <= {int(year_max)}')
            if filters:
                filter_exprs[coll] = " and ".join(filters)

        # Search all collections
        all_hits = []
        for coll in collections:
            try:
                results = self.collection_manager.search(
                    coll, query_embedding,
                    top_k=top_k_per_collection,
                    filter_expr=filter_exprs.get(coll),
                )
                weight = COLLECTION_CONFIG.get(coll, {}).get("weight", 0.05)
                for r in results:
                    all_hits.append(SearchHit(
                        collection=coll,
                        id=r.get("id", ""),
                        score=r.get("score", 0.0) * weight,
                        text=r.get("text_chunk", r.get("text_summary", "")),
                        metadata={k: v for k, v in r.items() if k not in ("embedding", "text_chunk", "text_summary")},
                    ))
            except Exception as e:
                logger.warning(f"Search failed for {coll}: {e}")

        # Sort by weighted score
        all_hits.sort(key=lambda h: h.score, reverse=True)

        return CrossCollectionResult(
            query=query,
            hits=all_hits,
            knowledge_context=knowledge_context,
            total_collections_searched=len(collections),
            search_time_ms=(time.time() - start) * 1000,
        )

    def query(self, question: str, conversation_context: str = "", **kwargs) -> str:
        """Full RAG query: retrieve + synthesize."""
        evidence = self.retrieve(question, **kwargs)
        messages = self._build_prompt(question, evidence, conversation_context)
        return self.llm_client.generate(messages)

    def query_stream(self, question: str, conversation_context: str = "", **kwargs) -> Generator[str, None, None]:
        """Streaming RAG query."""
        evidence = self.retrieve(question, **kwargs)
        messages = self._build_prompt(question, evidence, conversation_context)
        yield from self.llm_client.generate_stream(messages)

    def find_related(self, entity: str, top_k: int = 5) -> Dict[str, List[SearchHit]]:
        """Find related evidence across collections for an entity."""
        embedding = self._embed_query(entity)
        results = {}
        for coll in COLLECTION_CONFIG:
            try:
                hits = self.collection_manager.search(coll, embedding, top_k=top_k)
                if hits:
                    results[coll] = [
                        SearchHit(collection=coll, id=h.get("id", ""), score=h.get("score", 0.0),
                                  text=h.get("text_chunk", h.get("text_summary", "")),
                                  metadata={k: v for k, v in h.items() if k not in ("embedding",)})
                        for h in hits
                    ]
            except Exception as exc:
                logger.warning(f"find_related search failed for {coll}: {exc}")
        return results

    def retrieve_comparative(self, question: str, **kwargs) -> ComparativeResult:
        """Retrieve evidence for a comparative query (e.g., CT vs MRI)."""
        start = time.time()
        # Try to extract two entities from the query
        entity_a_str, entity_b_str = self._parse_comparison_entities(question)

        evidence_a = self.retrieve(f"{question} {entity_a_str}", **kwargs)
        evidence_b = self.retrieve(f"{question} {entity_b_str}", **kwargs)

        entity_a_dict = resolve_comparison_entity(entity_a_str)
        entity_b_dict = resolve_comparison_entity(entity_b_str)
        comparison_context = ""
        if entity_a_dict and entity_b_dict:
            comparison_context = get_comparison_context(entity_a_dict, entity_b_dict)

        return ComparativeResult(
            query=question,
            entity_a=entity_a_str,
            entity_b=entity_b_str,
            evidence_a=evidence_a,
            evidence_b=evidence_b,
            comparison_context=comparison_context,
            total_search_time_ms=(time.time() - start) * 1000,
        )

    def _parse_comparison_entities(self, query: str) -> tuple:
        """Extract two entities from a comparative query."""
        for sep in [" vs ", " versus ", " compared to ", " vs. "]:
            if sep in query.lower():
                parts = query.lower().split(sep, 1)
                return parts[0].strip(), parts[1].strip().rstrip("?.")
        return query, ""
