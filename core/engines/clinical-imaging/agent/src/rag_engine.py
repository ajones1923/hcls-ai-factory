"""Multi-collection RAG engine for Imaging Intelligence Agent.

Searches 10 imaging collections + 1 genomic collection, augments results
with domain knowledge, and synthesizes answers via LLM.

Includes a full comparative analysis pipeline for head-to-head queries
(e.g., "CT vs MRI for lung nodule detection").
"""

import re
import time
from typing import Any, Dict, Generator, List, Optional, Set, Tuple

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


# ═══════════════════════════════════════════════════════════════════════════════
# Comparative analysis regex and entity map
# ═══════════════════════════════════════════════════════════════════════════════

_COMPARATIVE_RE = re.compile(
    r"\b(compare|compared to|vs\.?|versus|difference between|head.to.head|better than|advantages|disadvantages)\b",
    re.IGNORECASE,
)

# Patterns for extracting the two entities from a comparative query
_ENTITY_PATTERNS = [
    # "A vs B", "A vs. B"
    re.compile(r"(.+?)\s+vs\.?\s+(.+)", re.IGNORECASE),
    # "A versus B"
    re.compile(r"(.+?)\s+versus\s+(.+)", re.IGNORECASE),
    # "compare A and B", "compare A with B", "compare A to B"
    re.compile(r"compare\s+(.+?)\s+(?:and|with|to)\s+(.+)", re.IGNORECASE),
    # "A compared to B"
    re.compile(r"(.+?)\s+compared\s+to\s+(.+)", re.IGNORECASE),
    # "difference between A and B"
    re.compile(r"difference(?:s)?\s+between\s+(.+?)\s+and\s+(.+)", re.IGNORECASE),
    # "A head to head B" / "A head-to-head B"
    re.compile(r"(.+?)\s+head[\s-]to[\s-]head\s+(.+)", re.IGNORECASE),
    # "advantages of A over B"
    re.compile(r"advantages\s+of\s+(.+?)\s+over\s+(.+)", re.IGNORECASE),
    # "A better than B"
    re.compile(r"(.+?)\s+better\s+than\s+(.+)", re.IGNORECASE),
]

IMAGING_ENTITY_MAP: Dict[str, Dict[str, Any]] = {
    # ── Modalities ──
    "ct": {
        "type": "modality",
        "full_name": "Computed Tomography",
        "aliases": ["CT", "cat scan", "computed tomography"],
    },
    "mri": {
        "type": "modality",
        "full_name": "Magnetic Resonance Imaging",
        "aliases": ["MRI", "MR", "magnetic resonance"],
    },
    "pet": {
        "type": "modality",
        "full_name": "Positron Emission Tomography",
        "aliases": ["PET", "PET/CT", "PET-CT"],
    },
    "ultrasound": {
        "type": "modality",
        "full_name": "Ultrasound",
        "aliases": ["US", "sonography", "echo"],
    },
    "xray": {
        "type": "modality",
        "full_name": "X-ray",
        "aliases": ["radiograph", "plain film", "CXR"],
    },
    "mammography": {
        "type": "modality",
        "full_name": "Mammography",
        "aliases": ["mammo", "DBT", "tomosynthesis"],
    },
    # ── AI architectures ──
    "cnn": {
        "type": "architecture",
        "full_name": "Convolutional Neural Network",
        "aliases": ["CNN", "ConvNet"],
    },
    "transformer": {
        "type": "architecture",
        "full_name": "Vision Transformer",
        "aliases": ["ViT", "transformer", "attention"],
    },
    "unet": {
        "type": "architecture",
        "full_name": "U-Net",
        "aliases": ["U-Net", "nnU-Net", "segmentation network"],
    },
    # ── Techniques ──
    "dlir": {
        "type": "technique",
        "full_name": "Deep Learning Image Reconstruction",
        "aliases": ["DLIR", "DL reconstruction", "AI reconstruction"],
    },
    "iterative": {
        "type": "technique",
        "full_name": "Iterative Reconstruction",
        "aliases": ["IR", "MBIR", "model-based"],
    },
    # ── Tasks ──
    "detection": {
        "type": "task",
        "full_name": "AI Detection/CADe",
        "aliases": ["CADe", "detection", "screening"],
    },
    "segmentation": {
        "type": "task",
        "full_name": "AI Segmentation",
        "aliases": ["segmentation", "contouring", "delineation"],
    },
    "classification": {
        "type": "task",
        "full_name": "AI Classification/CADx",
        "aliases": ["CADx", "classification", "diagnosis"],
    },
    # ── Technologies ──
    "photon_counting": {
        "type": "technology",
        "full_name": "Photon-Counting CT",
        "aliases": ["PCCT", "photon counting", "photon-counting detector"],
    },
    "dual_energy": {
        "type": "technology",
        "full_name": "Dual-Energy CT",
        "aliases": ["DECT", "dual energy", "spectral CT"],
    },
    "conventional_ct": {
        "type": "technology",
        "full_name": "Conventional CT",
        "aliases": ["standard CT", "conventional", "energy-integrating"],
    },
}


def _resolve_imaging_entity(text: str) -> Dict[str, Any]:
    """Fuzzy-match user input against IMAGING_ENTITY_MAP.

    Tries exact key match, then alias substring matching, then falls back
    to the knowledge.py ``resolve_comparison_entity`` for pathologies,
    modalities, and anatomy lookups.

    Returns:
        Dict with at least ``key``, ``type``, ``full_name``, or empty dict
        if nothing matched.
    """
    cleaned = text.strip().lower().replace("-", "_").replace("/", "_")

    # 1. Exact key match
    if cleaned in IMAGING_ENTITY_MAP:
        entry = IMAGING_ENTITY_MAP[cleaned]
        return {"key": cleaned, **entry}

    # 2. Alias match (case-insensitive substring)
    for key, entry in IMAGING_ENTITY_MAP.items():
        for alias in entry.get("aliases", []):
            if alias.lower() == cleaned or cleaned == alias.lower().replace("-", "_"):
                return {"key": key, **entry}

    # 3. Fuzzy substring in aliases
    for key, entry in IMAGING_ENTITY_MAP.items():
        for alias in entry.get("aliases", []):
            if cleaned in alias.lower() or alias.lower() in cleaned:
                return {"key": key, **entry}

    # 4. Fall back to knowledge.py entity resolution
    kg_result = resolve_comparison_entity(text)
    if kg_result:
        return {
            "key": kg_result["canonical"],
            "type": kg_result["type"],
            "full_name": kg_result["canonical"].replace("_", " ").title(),
            "aliases": [],
        }

    # 5. Not found — return with the raw text so downstream still works
    return {
        "key": cleaned,
        "type": "unknown",
        "full_name": text.strip(),
        "aliases": [],
    }


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

COMPARATIVE_SYSTEM_PROMPT = """You are an expert medical imaging intelligence assistant performing
a structured comparative analysis. You must produce a comprehensive, evidence-based comparison
using the following format:

## Comparison: {entity_a} vs {entity_b}

### 1. Technical Specifications
| Feature | {entity_a} | {entity_b} |
|---------|-----------|-----------|
| ... | ... | ... |

### 2. Clinical Performance
Compare sensitivity, specificity, AUC, PPV, NPV where evidence is available.
| Metric | {entity_a} | {entity_b} |
|--------|-----------|-----------|
| ... | ... | ... |

### 3. Radiation Dose / Safety
Compare dose, contraindications, and safety considerations (if applicable).

### 4. AI Integration Capabilities
Compare AI model availability, FDA-cleared devices, and workflow integration.

### 5. Clinical Guidelines & Recommendations
What do ACR, RSNA, and other bodies recommend?

### 6. Cost-Effectiveness Considerations
Relative cost, availability, throughput.

### 7. Summary Recommendation
Provide a balanced summary with appropriate caveats. State when each option is preferred.

Rules:
- Use markdown tables for structured comparisons.
- Cite specific evidence from the retrieved context with [Source] tags.
- If evidence is missing for a comparison axis, say so explicitly.
- Always note that clinical decisions require physician oversight.
- Be balanced — avoid bias toward either entity unless evidence strongly favors one."""


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

    # ═══════════════════════════════════════════════════════════════════════
    # Comparative analysis helpers
    # ═══════════════════════════════════════════════════════════════════════

    def _is_comparative(self, query: str) -> bool:
        """Detect whether a query is asking for a comparison."""
        return bool(_COMPARATIVE_RE.search(query))

    def _parse_comparison_entities(self, query: str) -> Tuple[str, str]:
        """Extract two entities from a comparative query using regex patterns.

        Tries multiple patterns in order and returns the first match.
        Falls back to splitting on the first comparative keyword.
        """
        # Strip trailing punctuation for cleaner matching
        clean_query = query.strip().rstrip("?.")

        for pattern in _ENTITY_PATTERNS:
            m = pattern.search(clean_query)
            if m:
                a = m.group(1).strip().rstrip("?.,")
                b = m.group(2).strip().rstrip("?.,")
                # Remove leading filler words
                for prefix in ["compare ", "what is the ", "what are the ", "is ", "are ",
                               "how does ", "how do ", "which is ", "when to use "]:
                    if a.lower().startswith(prefix):
                        a = a[len(prefix):]
                    if b.lower().startswith(prefix):
                        b = b[len(prefix):]
                # Strip trailing context phrases from entity B
                for suffix_pat in [r"\s+for\s+.+", r"\s+in\s+.+", r"\s+during\s+.+"]:
                    b = re.sub(suffix_pat, "", b, flags=re.IGNORECASE).strip()
                return a.strip().lower(), b.strip().lower()

        # Fallback: split on the first comparative keyword found
        for sep in [" vs ", " vs. ", " versus ", " compared to "]:
            if sep in clean_query.lower():
                idx = clean_query.lower().index(sep)
                a = clean_query[:idx].strip()
                b = clean_query[idx + len(sep):].strip()
                return a.lower(), b.lower()

        return clean_query, ""

    def _find_shared_evidence(
        self,
        hits_a: List[SearchHit],
        hits_b: List[SearchHit],
    ) -> List[SearchHit]:
        """Find evidence that appears in both result sets (same record id)."""
        ids_a = {h.id for h in hits_a}
        ids_b = {h.id for h in hits_b}
        shared_ids = ids_a & ids_b
        if not shared_ids:
            return []

        # Return the higher-scored version of each shared hit
        id_to_hit: Dict[str, SearchHit] = {}
        for h in hits_a + hits_b:
            if h.id in shared_ids:
                existing = id_to_hit.get(h.id)
                if existing is None or h.score > existing.score:
                    id_to_hit[h.id] = h
        shared = list(id_to_hit.values())
        shared.sort(key=lambda h: h.score, reverse=True)
        return shared

    def _build_comparative_prompt(
        self,
        question: str,
        result: ComparativeResult,
        conversation_context: str = "",
    ) -> List[Dict[str, str]]:
        """Build a prompt instructing the LLM to produce structured comparison."""
        entity_a_label = result.entity_a_resolved.get("full_name", result.entity_a)
        entity_b_label = result.entity_b_resolved.get("full_name", result.entity_b)

        # Format the system prompt with entity names
        sys_prompt = COMPARATIVE_SYSTEM_PROMPT.replace("{entity_a}", entity_a_label).replace("{entity_b}", entity_b_label)

        # Evidence sections
        sections = []

        # Entity A evidence
        evidence_a_text = ""
        for hit in result.evidence_a.hits[:15]:
            config = COLLECTION_CONFIG.get(hit.collection, {})
            label = config.get("label", hit.collection)
            evidence_a_text += f"\n[{label}] (score: {hit.score:.3f}) {hit.text}\n"
        if evidence_a_text:
            sections.append(f"## Evidence for {entity_a_label}\n{evidence_a_text}")

        # Entity B evidence
        evidence_b_text = ""
        for hit in result.evidence_b.hits[:15]:
            config = COLLECTION_CONFIG.get(hit.collection, {})
            label = config.get("label", hit.collection)
            evidence_b_text += f"\n[{label}] (score: {hit.score:.3f}) {hit.text}\n"
        if evidence_b_text:
            sections.append(f"## Evidence for {entity_b_label}\n{evidence_b_text}")

        # Shared evidence
        if result.shared_evidence:
            shared_text = ""
            for hit in result.shared_evidence[:10]:
                config = COLLECTION_CONFIG.get(hit.collection, {})
                label = config.get("label", hit.collection)
                shared_text += f"\n[{label}] (score: {hit.score:.3f}) {hit.text}\n"
            sections.append(f"## Shared / Cross-Referenced Evidence\n{shared_text}")

        # Knowledge graph context
        if result.comparison_context:
            sections.append(f"## Domain Knowledge\n{result.comparison_context}")

        # Conversation context
        if conversation_context:
            sections.append(f"## Conversation Context\n{conversation_context}")

        user_content = "\n\n".join(sections)
        user_content += f"\n\n## Question\n{question}"

        return [
            {"role": "system", "content": sys_prompt},
            {"role": "user", "content": user_content},
        ]

    def retrieve_comparative(self, question: str, **kwargs) -> ComparativeResult:
        """Retrieve evidence for a comparative query with dual retrieval.

        Steps:
        1. Parse two entities from the query.
        2. Resolve each entity against IMAGING_ENTITY_MAP.
        3. Run dual retrieval (one search per entity).
        4. Find shared evidence (same record_id in both result sets).
        5. Build knowledge-graph comparison context.
        6. Return structured ComparativeResult.
        """
        start = time.time()

        # 1. Parse entities
        entity_a_str, entity_b_str = self._parse_comparison_entities(question)
        logger.info(f"Comparative parse: '{entity_a_str}' vs '{entity_b_str}'")

        # 2. Resolve entities
        entity_a_resolved = _resolve_imaging_entity(entity_a_str)
        entity_b_resolved = _resolve_imaging_entity(entity_b_str)
        logger.info(
            f"Entity A resolved: {entity_a_resolved.get('key', '?')} "
            f"({entity_a_resolved.get('type', 'unknown')})"
        )
        logger.info(
            f"Entity B resolved: {entity_b_resolved.get('key', '?')} "
            f"({entity_b_resolved.get('type', 'unknown')})"
        )

        # 3. Dual retrieval — search with entity-augmented queries
        search_a = f"{question} {entity_a_resolved.get('full_name', entity_a_str)}"
        search_b = f"{question} {entity_b_resolved.get('full_name', entity_b_str)}"
        evidence_a = self.retrieve(search_a, **kwargs)
        evidence_b = self.retrieve(search_b, **kwargs)

        # 4. Find shared evidence
        shared = self._find_shared_evidence(evidence_a.hits, evidence_b.hits)

        # 5. Knowledge-graph comparison context
        kg_a = resolve_comparison_entity(entity_a_str)
        kg_b = resolve_comparison_entity(entity_b_str)
        comparison_context = ""
        if kg_a and kg_b:
            comparison_context = get_comparison_context(kg_a, kg_b)

        return ComparativeResult(
            query=question,
            entity_a=entity_a_str,
            entity_b=entity_b_str,
            entity_a_resolved=entity_a_resolved,
            entity_b_resolved=entity_b_resolved,
            evidence_a=evidence_a,
            evidence_b=evidence_b,
            shared_evidence=shared,
            comparison_context=comparison_context,
            total_search_time_ms=(time.time() - start) * 1000,
        )

    def _handle_comparative(
        self,
        question: str,
        conversation_context: str = "",
        **kwargs,
    ) -> str:
        """End-to-end comparative analysis: retrieve -> prompt -> LLM."""
        comp_result = self.retrieve_comparative(question, **kwargs)
        messages = self._build_comparative_prompt(question, comp_result, conversation_context)
        synthesis = self.llm_client.generate(messages)

        # Store the synthesis on the result for callers that want both
        comp_result.comparative_synthesis = synthesis

        # Also stash the latest comparative result so the UI can access it
        self._last_comparative_result = comp_result

        return synthesis

    def _handle_comparative_stream(
        self,
        question: str,
        conversation_context: str = "",
        **kwargs,
    ) -> Generator[str, None, None]:
        """Streaming comparative analysis."""
        comp_result = self.retrieve_comparative(question, **kwargs)
        messages = self._build_comparative_prompt(question, comp_result, conversation_context)
        self._last_comparative_result = comp_result
        yield from self.llm_client.generate_stream(messages)

    # ═══════════════════════════════════════════════════════════════════════
    # Standard prompt builder
    # ═══════════════════════════════════════════════════════════════════════

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

    # ═══════════════════════════════════════════════════════════════════════
    # Core retrieval
    # ═══════════════════════════════════════════════════════════════════════

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

    # ═══════════════════════════════════════════════════════════════════════
    # Main query entry points (with comparative branch)
    # ═══════════════════════════════════════════════════════════════════════

    def query(self, question: str, conversation_context: str = "", **kwargs) -> str:
        """Full RAG query: retrieve + synthesize.

        Automatically detects comparative queries and routes them through
        the comparative analysis pipeline.
        """
        if self._is_comparative(question):
            logger.info(f"Comparative query detected: {question}")
            return self._handle_comparative(question, conversation_context, **kwargs)

        evidence = self.retrieve(question, **kwargs)
        messages = self._build_prompt(question, evidence, conversation_context)
        return self.llm_client.generate(messages)

    def query_stream(self, question: str, conversation_context: str = "", **kwargs) -> Generator[str, None, None]:
        """Streaming RAG query."""
        if self._is_comparative(question):
            logger.info(f"Comparative query (streaming) detected: {question}")
            yield from self._handle_comparative_stream(question, conversation_context, **kwargs)
            return

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

    def get_last_comparative_result(self) -> Optional[ComparativeResult]:
        """Return the most recent ComparativeResult (set by query/query_stream)."""
        return getattr(self, "_last_comparative_result", None)
