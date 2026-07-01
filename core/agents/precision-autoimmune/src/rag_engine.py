"""
Precision Autoimmune Agent — RAG Engine

Multi-collection retrieval-augmented generation for autoimmune disease analysis.
Performs weighted parallel search across 14 Milvus collections, applies
domain-specific reranking, and synthesizes responses via Claude.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import hashlib
import re
import threading
import time
from collections import deque
from dataclasses import dataclass, field
from typing import Any, Dict, Generator, List, Optional

from loguru import logger

# ── Data classes ──────────────────────────────────────────────────────────

@dataclass
class SearchHit:
    collection: str
    id: str
    score: float
    text: str
    metadata: Dict[str, Any] = field(default_factory=dict)
    relevance: str = "low"  # high / medium / low


@dataclass
class CrossCollectionResult:
    query: str
    hits: List[SearchHit] = field(default_factory=list)
    knowledge_context: str = ""
    total_collections_searched: int = 0
    search_time_ms: float = 0.0


# ── Input validation ─────────────────────────────────────────────────────

_SAFE_FILTER_RE = re.compile(r"^[A-Za-z0-9 _\-\.]+$")


def _sanitize_filter_value(value: str, max_len: int = 64) -> Optional[str]:
    """Validate and sanitize values used in Milvus filter expressions."""
    if not value or len(value) > max_len:
        return None
    if not _SAFE_FILTER_RE.match(value):
        logger.warning(f"Rejected unsafe filter value: {value!r}")
        return None
    return value


# ── System prompt ─────────────────────────────────────────────────────────

SYSTEM_PROMPT = """You are the Precision Autoimmune Intelligence Agent, an expert clinical \
decision-support system specializing in autoimmune disease analysis. You are part of the \
HCLS AI Factory platform running on NVIDIA DGX Spark.

Your expertise spans:
- **Autoantibody interpretation**: ANA patterns (homogeneous, speckled, nucleolar, centromere), \
anti-dsDNA, anti-CCP, anti-SSA/SSB, anti-Scl-70, anti-Jo-1, AChR, anti-tTG, and 14+ \
autoantibody types with sensitivity/specificity data
- **HLA-disease associations**: 50+ HLA allele-disease associations with odds ratios \
(e.g., HLA-B*27:05 → AS OR=87.4, HLA-DRB1*04:01 → RA OR=4.2, HLA-C*06:02 → Psoriasis OR=10.0)
- **Disease activity scoring**: DAS28-CRP/ESR for RA, SLEDAI-2K for SLE, CDAI for RA, \
BASDAI for AS — with remission/low/moderate/high thresholds
- **Flare prediction**: Biomarker pattern recognition for RA (CRP/ESR/IL-6/MMP-3), \
SLE (anti-dsDNA/complement/lymphocytes), IBD (calprotectin/CRP/albumin)
- **Biologic therapy selection**: TNF inhibitors, anti-CD20, IL-6R, IL-17A, BLyS, \
IL-12/23, T-cell co-stimulation, JAK inhibitors — with PGx considerations, \
contraindications, monitoring requirements
- **Diagnostic odyssey analysis**: Identifying patterns across fragmented clinical \
records spanning years of multi-specialist visits (avg 4.6 years for SLE, 5-7 for POTS)
- **ACR/EULAR classification criteria**: 2010 RA, 2019 SLE, ASAS axSpA, and more
- **Cross-disease overlap syndromes**: MCTD, rhupus, POTS/hEDS/MCAS triad, \
shared pathogenic mechanisms
- **Pharmacogenomics**: CYP2C19/CYP3A4 metabolism of JAK inhibitors, \
FCGR3A polymorphisms affecting rituximab ADCC, HLA-based drug response prediction

When answering:
1. Always cite evidence sources using these formats:
   - Autoantibodies: [AutoAb:name] (e.g., [AutoAb:anti-CCP])
   - HLA alleles: [HLA:allele] (e.g., [HLA:B*27:05])
   - Activity scores: [Activity:score_name] (e.g., [Activity:DAS28-CRP])
   - Therapies: [Therapy:drug] (e.g., [Therapy:Adalimumab])
   - Literature: [Literature:PMID](https://pubmed.ncbi.nlm.nih.gov/PMID/)
   - Clinical trials: [Trial:NCT_ID](https://clinicaltrials.gov/study/NCT_ID)
2. Distinguish confirmed diagnoses from differential considerations
3. Flag critical alerts: HIGH DISEASE ACTIVITY, IMMINENT FLARE, strong HLA associations (OR>5)
4. Consider pharmacogenomic implications for therapy recommendations
5. Note the diagnostic odyssey timeline when analyzing patient records
6. Provide actionable, evidence-based recommendations with specific next steps
7. When data is insufficient, clearly state what additional tests would be informative
8. Cross-reference across domains (e.g., HLA + autoantibodies + biomarkers) for integrated assessment

IMPORTANT: This is a clinical decision-support tool. All recommendations should be \
reviewed by a qualified healthcare provider before clinical action."""


# ── Disease area detection ────────────────────────────────────────────────

DISEASE_KEYWORDS = {
    "rheumatoid_arthritis": [
        "rheumatoid", " ra ", "anti-ccp", "das28", "joint swelling", "synovitis",
        "morning stiffness", "erosive arthritis", "citrullinated",
    ],
    "systemic_lupus": [
        "lupus", " sle ", "anti-dsdna", "sledai", "butterfly rash", "nephritis",
        "complement", "photosensitivity", "malar rash", "discoid",
    ],
    "multiple_sclerosis": [
        "multiple sclerosis", " ms ", "demyelinating", "oligoclonal",
        "relapsing remitting", "interferon beta",
    ],
    "ankylosing_spondylitis": [
        "ankylosing", "spondylitis", "basdai", "hla-b27", "sacroiliitis",
        "axial spondyloarthritis", "bamboo spine", "uveitis",
    ],
    "sjogrens": [
        "sjogren", "sjögren", "dry eyes", "dry mouth", "schirmer",
        "anti-ssa", "anti-ssb", "anti-ro", "anti-la", "sicca",
    ],
    "systemic_sclerosis": [
        "scleroderma", "systemic sclerosis", "raynaud", "anti-scl-70",
        "anti-centromere", "skin thickening", "pulmonary fibrosis", "crest",
    ],
    "inflammatory_bowel": [
        "crohn", "colitis", "ibd", "calprotectin", "ulcerative",
        "inflammatory bowel", "fistula",
    ],
    "psoriasis": [
        "psoriasis", "psoriatic", "pasi", "il-17", "plaque", "enthesitis",
    ],
    "myasthenia_gravis": [
        "myasthenia", "achr", "ptosis", "diplopia", "thymoma",
        "acetylcholine receptor", "lambert-eaton",
    ],
    "celiac": [
        "celiac", "coeliac", "anti-ttg", "villous atrophy", "gluten",
        "tissue transglutaminase", "endomysial",
    ],
    "thyroid_autoimmune": [
        "graves", "hashimoto", "thyroid", "tsi", "anti-tpo",
        "thyroid peroxidase", "thyroglobulin",
    ],
    "pots_eds_mcas": [
        "pots", "dysautonomia", "hypermobility", "eds", "mcas",
        "mast cell", "tilt table", "beighton", "tryptase",
    ],
}

# ── Comparative query detection ───────────────────────────────────────────

_COMPARATIVE_RE = re.compile(
    r"\b(compare|vs\.?|versus|difference between|head[\s-]to[\s-]head|"
    r"which is better|pros and cons)\b",
    re.IGNORECASE,
)


class AutoimmuneRAGEngine:
    """Multi-collection RAG engine for autoimmune disease analysis."""

    def __init__(
        self,
        collection_manager,
        embedder,
        llm_client,
        settings=None,
        knowledge=None,
    ):
        self.cm = collection_manager
        self.embedder = embedder
        self.llm = llm_client
        self.knowledge = knowledge

        if settings is None:
            from config.settings import settings as _settings
            settings = _settings
        self.settings = settings
        self._conversation_lock = threading.Lock()
        self._conversation_history: deque = deque(maxlen=self.settings.CONVERSATION_MEMORY_SIZE if self.settings else 3)
        self._embed_cache: Dict[str, List[float]] = {}
        self._embed_cache_max = 256

    # ── Service readiness ─────────────────────────────────────────────
    @property
    def is_ready(self) -> bool:
        """Check if all required services are available."""
        return self.embedder is not None and self.llm is not None

    @property
    def can_search(self) -> bool:
        """Check if vector search is available (embedder required)."""
        return self.embedder is not None

    # ── Embedding ─────────────────────────────────────────────────────
    def _embed(self, text: str) -> List[float]:
        if self.embedder is None:
            raise RuntimeError(
                "Embedder not initialized — cannot encode queries. "
                "Ensure sentence-transformers is installed and the model is loaded."
            )

        # Check cache
        cache_key = text[:512]
        if cache_key in self._embed_cache:
            return self._embed_cache[cache_key]

        instruction = self.settings.BGE_INSTRUCTION
        vec = self.embedder.encode(instruction + text).tolist()

        # Evict oldest if cache full
        if len(self._embed_cache) >= self._embed_cache_max:
            oldest_key = next(iter(self._embed_cache))
            del self._embed_cache[oldest_key]
        self._embed_cache[cache_key] = vec

        return vec

    # ── Disease area detection ────────────────────────────────────────
    def _detect_disease_areas(self, query: str) -> List[str]:
        query_lower = query.lower()
        detected = []
        for area, keywords in DISEASE_KEYWORDS.items():
            if any(kw in query_lower for kw in keywords):
                detected.append(area)
        return detected

    def _is_comparative(self, query: str) -> bool:
        """Detect if query is a comparison between two entities."""
        return bool(_COMPARATIVE_RE.search(query))

    # ── Knowledge augmentation ────────────────────────────────────────
    def _build_knowledge_context(self, query: str, disease_areas: List[str]) -> str:
        if not self.knowledge:
            return ""

        parts = []

        try:
            from .knowledge import (
                AUTOANTIBODY_DISEASE_MAP,
                BIOLOGIC_THERAPIES,
                FLARE_BIOMARKER_PATTERNS,
                HLA_DISEASE_ASSOCIATIONS,
            )
        except ImportError:
            return ""

        query_lower = query.lower()

        # Add relevant HLA associations
        for allele, assocs in HLA_DISEASE_ASSOCIATIONS.items():
            if allele.lower() in query_lower:
                for a in assocs:
                    parts.append(
                        f"[HLA Knowledge] {allele} → {a['disease'].replace('_', ' ').title()} "
                        f"(OR={a['odds_ratio']}, PMID:{a.get('pmid', 'N/A')}). "
                        f"{a.get('note', '')}"
                    )

        # Add relevant autoantibody mappings
        for ab, assocs in AUTOANTIBODY_DISEASE_MAP.items():
            if ab.lower() in query_lower:
                for a in assocs:
                    parts.append(
                        f"[Autoantibody Knowledge] {ab} → {a['disease'].replace('_', ' ').title()} "
                        f"(sensitivity={a['sensitivity']}, specificity={a['specificity']}). "
                        f"{a.get('note', '')}"
                    )

        # Add relevant biologic therapy context
        for therapy in BIOLOGIC_THERAPIES:
            drug_lower = therapy["drug_name"].lower()
            class_lower = therapy["drug_class"].lower()
            if drug_lower in query_lower or class_lower in query_lower:
                diseases = ", ".join(d.replace("_", " ").title() for d in therapy["indicated_diseases"])
                pgx = "; ".join(therapy.get("pgx_considerations", []))
                parts.append(
                    f"[Therapy Knowledge] {therapy['drug_name']} ({therapy['drug_class']}): "
                    f"{therapy.get('mechanism', '')}. Indicated for: {diseases}. PGx: {pgx}"
                )

        # Add flare pattern context
        for disease, pattern in FLARE_BIOMARKER_PATTERNS.items():
            if disease.replace("_", " ") in query_lower or disease in query_lower:
                markers = ", ".join(pattern["early_warning_biomarkers"])
                parts.append(
                    f"[Flare Knowledge] {disease.replace('_', ' ').title()}: "
                    f"Early warning biomarkers: {markers}. "
                    f"Protective signals: {', '.join(pattern.get('protective_signals', []))}"
                )

        return "\n".join(parts[:25])  # Cap at 25 knowledge items

    # ── Retrieval ─────────────────────────────────────────────────────
    def retrieve(
        self,
        query: str,
        top_k_per_collection: Optional[int] = None,
        collections_filter: Optional[List[str]] = None,
        year_min: Optional[int] = None,
        year_max: Optional[int] = None,
        patient_id: Optional[str] = None,
    ) -> CrossCollectionResult:
        """Retrieve evidence from multiple Milvus collections."""
        start = time.time()

        if not self.can_search:
            logger.warning("RAG retrieve called but embedder not available")
            return CrossCollectionResult(query=query, search_time_ms=0)

        if top_k_per_collection is None:
            top_k_per_collection = self.settings.TOP_K_PER_COLLECTION

        try:
            query_embedding = self._embed(query)
        except Exception as exc:
            logger.error(f"Embedding failed: {exc}")
            return CrossCollectionResult(query=query, search_time_ms=0)

        disease_areas = self._detect_disease_areas(query)

        # Build filter expressions with injection prevention
        filter_exprs = {}
        if patient_id:
            safe_pid = _sanitize_filter_value(patient_id)
            if safe_pid:
                for coll_name in (self.settings.COLL_CLINICAL_DOCUMENTS,
                                  self.settings.COLL_PATIENT_LABS,
                                  self.settings.COLL_PATIENT_TIMELINES):
                    filter_exprs[coll_name] = f'patient_id == "{safe_pid}"'

        if year_min and self.settings.COLL_LITERATURE in (collections_filter or self.settings.all_collection_names):
            yr_filter = f"year >= {int(year_min)}"
            if year_max:
                yr_filter += f" and year <= {int(year_max)}"
            filter_exprs[self.settings.COLL_LITERATURE] = yr_filter

        # Search
        try:
            raw_results = self.cm.search_all(
                query_embedding=query_embedding,
                top_k_per_collection=top_k_per_collection,
                collections=collections_filter,
                filter_exprs=filter_exprs,
                score_threshold=self.settings.SCORE_THRESHOLD,
            )
        except Exception as exc:
            logger.error(f"Milvus search_all failed: {exc}")
            raw_results = {}

        # Convert to SearchHit objects with weighted scoring
        config = self.settings.collection_config
        hits = []
        seen_ids = set()
        seen_texts = set()

        for coll_name, coll_hits in raw_results.items():
            weight = config.get(coll_name, {}).get("weight", 0.05)
            config.get(coll_name, {}).get("label", coll_name)

            for h in coll_hits:
                hit_id = h["id"]
                if hit_id in seen_ids:
                    continue
                seen_ids.add(hit_id)

                text = h.get("text_chunk", h.get("text_summary", ""))
                # Deduplicate by text content hash
                text_hash = hashlib.md5(text[:300].encode()).hexdigest()
                if text_hash in seen_texts:
                    continue
                seen_texts.add(text_hash)

                min(h["score"] * (1 + weight), 1.0)
                relevance = self._score_relevance(h["score"])

                hits.append(SearchHit(
                    collection=coll_name,
                    id=hit_id,
                    score=h["score"],
                    text=text,
                    metadata={k: v for k, v in h.items() if k not in ("id", "score", "text_chunk", "embedding")},
                    relevance=relevance,
                ))

        # Sort by raw score descending, cap at max_evidence
        hits.sort(key=lambda x: x.score, reverse=True)
        hits = hits[: self.settings.MAX_EVIDENCE_ITEMS]

        # Build knowledge context
        knowledge_ctx = self._build_knowledge_context(query, disease_areas)

        elapsed = (time.time() - start) * 1000

        logger.info(
            f"RAG retrieve: {len(hits)} hits from {len(raw_results)} collections "
            f"in {elapsed:.0f}ms (query: {query[:80]}...)"
        )

        return CrossCollectionResult(
            query=query,
            hits=hits,
            knowledge_context=knowledge_ctx,
            total_collections_searched=len(raw_results),
            search_time_ms=round(elapsed, 1),
        )

    def _score_relevance(self, score: float) -> str:
        if score >= self.settings.CITATION_HIGH:
            return "high"
        elif score >= self.settings.CITATION_MEDIUM:
            return "medium"
        return "low"

    # ── Search (evidence only) ────────────────────────────────────────
    def search(self, question: str, **kwargs) -> List[SearchHit]:
        """Return evidence hits without LLM synthesis."""
        result = self.retrieve(question, **kwargs)
        return result.hits

    # ── Full RAG query ────────────────────────────────────────────────
    def query(
        self,
        question: str,
        patient_context: Optional[str] = None,
        **kwargs,
    ) -> str:
        """Full RAG: retrieve evidence, synthesize with Claude."""
        if not self.is_ready:
            return "Error: RAG engine not fully initialized. Ensure embedder and LLM client are available."
        result = self.retrieve(question, **kwargs)
        return self._synthesize(question, result, patient_context)

    def query_stream(
        self,
        question: str,
        patient_context: Optional[str] = None,
        **kwargs,
    ) -> Generator[str, None, None]:
        """Streaming RAG: yield tokens as they arrive."""
        if not self.is_ready:
            yield "Error: RAG engine not fully initialized."
            return
        result = self.retrieve(question, **kwargs)
        yield from self._synthesize_stream(question, result, patient_context)

    # ── LLM synthesis ─────────────────────────────────────────────────
    def _build_evidence_block(self, result: CrossCollectionResult) -> str:
        parts = []
        config = self.settings.collection_config
        for i, hit in enumerate(result.hits, 1):
            label = config.get(hit.collection, {}).get("label", hit.collection)
            relevance_tag = f"[{hit.relevance} relevance]"

            # Build citation reference
            citation = f"[{label}:{hit.id}]"
            pmid = hit.metadata.get("pmid", "")
            if pmid:
                citation = f"[{label}:{hit.id}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)"
            nct_id = hit.metadata.get("nct_id", "")
            if nct_id:
                citation = f"[{label}:{nct_id}](https://clinicaltrials.gov/study/{nct_id})"

            source = f"{citation} (score={hit.score:.3f}) {relevance_tag}"
            parts.append(f"--- Evidence {i} {source} ---\n{hit.text[:1500]}")
        return "\n\n".join(parts)

    def _build_messages(
        self,
        question: str,
        result: CrossCollectionResult,
        patient_context: Optional[str] = None,
    ) -> List[Dict[str, str]]:
        evidence_block = self._build_evidence_block(result)

        user_content = f"**Question:** {question}\n\n"
        if patient_context:
            user_content += f"**Patient Context:**\n{patient_context}\n\n"
        if result.knowledge_context:
            user_content += f"**Knowledge Base Context:**\n{result.knowledge_context}\n\n"
        user_content += (
            f"**Retrieved Evidence ({len(result.hits)} items from "
            f"{result.total_collections_searched} collections, "
            f"search time: {result.search_time_ms:.0f}ms):**\n\n{evidence_block}"
        )

        # Add final instruction
        user_content += (
            "\n\n---\n"
            "Cite sources using the [Label:ID] format above. "
            "Prioritize [high relevance] evidence. "
            "Integrate across HLA, autoantibody, biomarker, and therapy domains. "
            "Flag any critical clinical alerts."
        )

        # Include conversation history (limited context)
        messages = []
        with self._conversation_lock:
            history = list(self._conversation_history)
        for entry in history:
            messages.append({"role": "user", "content": entry["question"][:200]})
            messages.append({"role": "assistant", "content": entry["answer"][:800]})
        messages.append({"role": "user", "content": user_content})

        return messages

    def _synthesize(
        self,
        question: str,
        result: CrossCollectionResult,
        patient_context: Optional[str] = None,
    ) -> str:
        if self.llm is None:
            return "Error: LLM client not available. Cannot synthesize response."

        messages = self._build_messages(question, result, patient_context)

        try:
            response = self.llm.messages.create(
                model=self.settings.LLM_MODEL,
                max_tokens=self.settings.LLM_MAX_TOKENS,
                system=SYSTEM_PROMPT,
                messages=messages,
            )
            answer = response.content[0].text

            with self._conversation_lock:
                self._conversation_history.append({"question": question, "answer": answer})
            return answer

        except Exception as exc:
            logger.error(f"LLM synthesis failed: {exc}", exc_info=True)
            # Return evidence summary as fallback
            if result.hits:
                fallback = f"*LLM synthesis unavailable ({type(exc).__name__}). Showing top evidence:*\n\n"
                for i, hit in enumerate(result.hits[:5], 1):
                    fallback += f"**{i}.** [{hit.collection}] (score={hit.score:.3f})\n{hit.text[:400]}\n\n"
                return fallback
            return f"Error generating response: {type(exc).__name__}: {exc}"

    def _synthesize_stream(
        self,
        question: str,
        result: CrossCollectionResult,
        patient_context: Optional[str] = None,
    ) -> Generator[str, None, None]:
        if self.llm is None:
            yield "Error: LLM client not available."
            return

        messages = self._build_messages(question, result, patient_context)
        full_answer = []

        try:
            with self.llm.messages.stream(
                model=self.settings.LLM_MODEL,
                max_tokens=self.settings.LLM_MAX_TOKENS,
                system=SYSTEM_PROMPT,
                messages=messages,
            ) as stream:
                for text in stream.text_stream:
                    full_answer.append(text)
                    yield text

            with self._conversation_lock:
                self._conversation_history.append({
                    "question": question,
                    "answer": "".join(full_answer),
                })
        except Exception as exc:
            logger.error(f"LLM streaming failed: {exc}", exc_info=True)
            yield f"\n\n*Streaming interrupted: {type(exc).__name__}*"

    # ── Cross-collection entity search ────────────────────────────────
    def find_related(self, entity: str, top_k: int = 5) -> Dict[str, List[SearchHit]]:
        """Find mentions of an entity across all collections."""
        result = self.retrieve(entity, top_k_per_collection=top_k)
        grouped: Dict[str, List[SearchHit]] = {}
        for hit in result.hits:
            grouped.setdefault(hit.collection, []).append(hit)
        return grouped

    def clear_history(self) -> None:
        """Reset conversation memory."""
        with self._conversation_lock:
            self._conversation_history.clear()
