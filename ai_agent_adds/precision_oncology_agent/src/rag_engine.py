"""
Multi-Collection RAG Engine for Oncology Intelligence
=======================================================

Provides the ``OncoRAGEngine`` class that orchestrates retrieval-augmented
generation across 11 oncology-specific Milvus collections, weighted scoring,
citation formatting, knowledge injection, query expansion, comparative
retrieval, and streaming LLM output.

Author: Adam Jones
Date:   February 2026
"""

from __future__ import annotations

import logging
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Dict, Generator, List, Optional, Set, Tuple

from config.settings import settings
from src.models import AgentQuery, CrossCollectionResult, SearchHit

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# System prompt
# ---------------------------------------------------------------------------
ONCO_SYSTEM_PROMPT = """\
You are an **Oncology Intelligence Agent** — an expert AI assistant
purpose-built for clinical and translational oncology decision support.

Your core competencies include:
* **Molecular profiling** — somatic/germline variant interpretation, tumor
  mutational burden, microsatellite instability, copy-number alterations,
  and gene fusions.
* **Variant interpretation** — CIViC and OncoKB evidence levels (Tier I–IV),
  AMP/ASCO/CAP classification, clinical significance assessment.
* **Therapy selection** — NCCN and ESMO guideline-concordant treatment
  recommendations, FDA-approved indications, and emerging therapies.
* **Clinical trial matching** — eligibility assessment against active
  ClinicalTrials.gov registrations, basket/umbrella trial awareness.
* **Resistance mechanisms** — on-target mutations, bypass signaling, lineage
  plasticity, and actionable resistance biomarkers.
* **Biomarker assessment** — TMB, MSI, PD-L1, HRD scoring; companion
  diagnostic requirements and assay considerations.
* **Outcomes monitoring** — RECIST response criteria, survival endpoints,
  minimal residual disease (MRD) tracking, ctDNA dynamics.
* **Cross-modal integration** — linking genomic findings to imaging,
  pathology, and drug-discovery pipelines (MolMIM / DiffDock).

Behavioral instructions:
1. **Cite evidence** — Reference source documents with clickable PubMed or
   ClinicalTrials.gov links wherever possible.
2. **Think cross-functionally** — Connect genomic variants to downstream
   therapy options, trials, and resistance patterns.
3. **Highlight resistance & contraindications** — Proactively note known
   resistance mechanisms and relevant safety concerns.
4. **Reference guidelines** — Cite NCCN, ESMO, or ASCO guidelines when
   making treatment recommendations.
5. **Acknowledge uncertainty** — Clearly state evidence gaps, limited data,
   or situations requiring multidisciplinary tumor board review.
"""

# ---------------------------------------------------------------------------
# Collection configuration
# ---------------------------------------------------------------------------
COLLECTION_CONFIG: Dict[str, Dict[str, Any]] = {
    "onco_variants": {
        "weight": settings.WEIGHT_VARIANTS,
        "label": "Variant",
        "filter_field": "gene",
        "year_field": None,
    },
    "onco_literature": {
        "weight": settings.WEIGHT_LITERATURE,
        "label": "Literature",
        "filter_field": "gene",
        "year_field": "year",
    },
    "onco_therapies": {
        "weight": settings.WEIGHT_THERAPIES,
        "label": "Therapy",
        "filter_field": None,
        "year_field": None,
    },
    "onco_guidelines": {
        "weight": settings.WEIGHT_GUIDELINES,
        "label": "Guideline",
        "filter_field": None,
        "year_field": "year",
    },
    "onco_trials": {
        "weight": settings.WEIGHT_TRIALS,
        "label": "Trial",
        "filter_field": None,
        "year_field": "start_year",
    },
    "onco_biomarkers": {
        "weight": settings.WEIGHT_BIOMARKERS,
        "label": "Biomarker",
        "filter_field": None,
        "year_field": None,
    },
    "onco_resistance": {
        "weight": settings.WEIGHT_RESISTANCE,
        "label": "Resistance",
        "filter_field": "gene",
        "year_field": None,
    },
    "onco_pathways": {
        "weight": settings.WEIGHT_PATHWAYS,
        "label": "Pathway",
        "filter_field": None,
        "year_field": None,
    },
    "onco_outcomes": {
        "weight": settings.WEIGHT_OUTCOMES,
        "label": "Outcome",
        "filter_field": None,
        "year_field": None,
    },
    "onco_cases": {
        "weight": settings.WEIGHT_CASES,
        "label": "Case",
        "filter_field": None,
        "year_field": None,
    },
    "genomic_evidence": {
        "weight": settings.WEIGHT_GENOMIC,
        "label": "Genomic",
        "filter_field": None,
        "year_field": None,
    },
}

# Maximum evidence items forwarded to the LLM prompt
_MAX_EVIDENCE = 30

# BGE-small-en-v1.5 instruction prefix for retrieval queries
_BGE_INSTRUCTION = "Represent this sentence for searching relevant passages: "

# Comparative trigger words
_COMPARATIVE_RE = re.compile(
    r"\b(compare|vs\.?|versus|difference between|head.to.head)\b",
    re.IGNORECASE,
)


class OncoRAGEngine:
    """Multi-collection RAG engine for precision oncology queries.

    Parameters
    ----------
    collection_manager
        Backend that exposes ``.search(collection, vector, top_k, filters)``
        returning a list of ``SearchHit`` objects.
    embedder
        Callable / object with ``.encode(text) -> List[float]``.
    llm_client
        LLM wrapper exposing ``.chat(messages) -> str`` and
        ``.chat_stream(messages) -> Generator[str, None, None]``.
    knowledge : optional
        Domain knowledge store providing contextual look-ups.
    query_expander : optional
        Callable ``(str) -> List[str]`` returning expansion terms.
    """

    def __init__(
        self,
        collection_manager,
        embedder,
        llm_client=None,
        knowledge=None,
        query_expander=None,
        settings=None,
    ):
        self.collection_manager = collection_manager
        self.embedder = embedder
        self.llm_client = llm_client
        self.knowledge = knowledge
        self.query_expander = query_expander
        self.settings = settings

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def retrieve(
        self,
        query: AgentQuery,
        top_k: int = 10,
        collections_filter: Optional[List[str]] = None,
        year_min: Optional[int] = None,
        year_max: Optional[int] = None,
        conversation_context: Optional[str] = None,
        **kwargs,
    ) -> CrossCollectionResult:
        """Retrieve evidence across all configured collections.

        Parameters
        ----------
        query : AgentQuery
            Structured query object (contains ``.text`` and optional filters).
        top_k : int
            Per-collection hit limit.
        collections_filter : list[str], optional
            Restrict search to these collection names.
        year_min, year_max : int, optional
            Publication / start-year range filters.
        conversation_context : str, optional
            Previous conversation turns for contextual embedding.

        Returns
        -------
        CrossCollectionResult
            Merged, ranked evidence with per-hit citation metadata.
        """
        embed_text = query.text
        if conversation_context:
            embed_text = f"{conversation_context} {query.text}"

        query_vector = self._embed_query(embed_text)

        # Primary search across collections
        all_hits = self._search_all_collections(
            query_vector=query_vector,
            query=query,
            top_k=top_k,
            collections_filter=collections_filter,
            year_min=year_min,
            year_max=year_max,
        )

        # Expanded search (if expander available)
        if self.query_expander:
            expanded_hits = self._expanded_search(
                query=query,
                top_k=max(top_k // 2, 3),
                collections_filter=collections_filter,
                year_min=year_min,
                year_max=year_max,
            )
            all_hits.extend(expanded_hits)

        ranked = self._merge_and_rank(all_hits)

        target_collections = (
            [c for c in collections_filter if c in COLLECTION_CONFIG]
            if collections_filter
            else list(COLLECTION_CONFIG.keys())
        )

        return CrossCollectionResult(
            query=query.text,
            hits=ranked,
            total_collections_searched=len(target_collections),
        )

    def cross_collection_search(self, query: AgentQuery) -> List[SearchHit]:
        """Search all collections and return ranked hits.

        This is the entry point used by the OncoIntelligenceAgent.

        Parameters
        ----------
        query : AgentQuery
            Structured query object.

        Returns
        -------
        list[SearchHit]
            Merged, ranked evidence hits.
        """
        result = self.retrieve(query)
        return result.hits

    def search(self, question: str, **kwargs) -> List[SearchHit]:
        """Evidence-only search (no LLM generation).

        Parameters
        ----------
        question : str
            Natural-language oncology question.
        **kwargs
            Forwarded to :meth:`retrieve`.

        Returns
        -------
        list[SearchHit]
            Ranked evidence hits.
        """
        agent_query = AgentQuery(question=question)
        result = self.retrieve(agent_query, **kwargs)
        return result.hits

    def synthesize(self, question: str, evidence: list, plan=None) -> "AgentResponse":  # noqa: F821
        """Synthesize an answer from pre-retrieved evidence.

        Used by the OncoIntelligenceAgent after gathering evidence.

        Parameters
        ----------
        question : str
            The original question.
        evidence : list
            Pre-retrieved evidence items (SearchHit or CrossCollectionResult).
        plan : SearchPlan, optional
            The search plan used.

        Returns
        -------
        AgentResponse
        """
        from src.models import AgentResponse, CrossCollectionResult

        # Flatten evidence if it's a list of CrossCollectionResult
        flat_hits = []
        for item in evidence:
            if hasattr(item, "hits"):
                flat_hits.extend(item.hits)
            elif hasattr(item, "score"):
                flat_hits.append(item)

        # Build prompt and get LLM answer
        prompt = self._build_prompt(question, flat_hits[:_MAX_EVIDENCE])
        messages = [
            {"role": "system", "content": ONCO_SYSTEM_PROMPT},
            {"role": "user", "content": prompt},
        ]

        answer = ""
        if self.llm_client is not None:
            try:
                answer = self.llm_client.chat(messages)
            except Exception as exc:
                logger.warning("LLM synthesis failed: %s", exc)
                answer = "Unable to synthesize answer — LLM unavailable."
        else:
            answer = "LLM client not configured. Evidence retrieved but synthesis unavailable."

        cross_result = CrossCollectionResult(
            query=question,
            hits=flat_hits,
            total_collections_searched=len(COLLECTION_CONFIG),
        )

        response = AgentResponse(
            question=question,
            answer=answer,
            evidence=cross_result,
        )
        # Attach plan if provided
        if plan is not None:
            response.plan = plan

        return response

    def query(self, question: str, **kwargs) -> str:
        """Full RAG pipeline: retrieve evidence then generate an answer.

        Parameters
        ----------
        question : str
            Natural-language oncology question.
        **kwargs
            Forwarded to :meth:`retrieve`.

        Returns
        -------
        str
            LLM-generated answer grounded in retrieved evidence.
        """
        agent_query = AgentQuery(question=question)

        if self._is_comparative(question):
            return self._handle_comparative(question, **kwargs)

        result = self.retrieve(agent_query, **kwargs)
        prompt = self._build_prompt(question, result.hits)
        messages = [
            {"role": "system", "content": ONCO_SYSTEM_PROMPT},
            {"role": "user", "content": prompt},
        ]
        return self.llm_client.chat(messages)

    def query_stream(self, question: str, **kwargs) -> Generator[str, None, None]:
        """Streaming variant of :meth:`query`.

        Yields
        ------
        str
            Incremental text chunks from the LLM.
        """
        agent_query = AgentQuery(question=question)

        if self._is_comparative(question):
            yield from self._handle_comparative_stream(question, **kwargs)
            return

        result = self.retrieve(agent_query, **kwargs)
        prompt = self._build_prompt(question, result.hits)
        messages = [
            {"role": "system", "content": ONCO_SYSTEM_PROMPT},
            {"role": "user", "content": prompt},
        ]
        yield from self.llm_client.chat_stream(messages)

    def find_related(
        self,
        entity: str,
        top_k: int = 5,
    ) -> Dict[str, List[SearchHit]]:
        """Find hits related to *entity* across every collection.

        Parameters
        ----------
        entity : str
            Gene name, drug name, pathway, or other oncology entity.
        top_k : int
            Per-collection hit limit.

        Returns
        -------
        dict[str, list[SearchHit]]
            Mapping of collection name to ranked hits.
        """
        vector = self._embed_query(entity)
        results: Dict[str, List[SearchHit]] = {}
        for coll_name in COLLECTION_CONFIG:
            try:
                hits = self.collection_manager.search(
                    collection=coll_name,
                    vector=vector,
                    top_k=top_k,
                )
                if hits:
                    results[coll_name] = hits
            except Exception:
                logger.warning("find_related failed for %s", coll_name, exc_info=True)
        return results

    # ------------------------------------------------------------------
    # Embedding
    # ------------------------------------------------------------------

    def _embed_query(self, text: str) -> List[float]:
        """Embed *text* with BGE instruction prefix for retrieval."""
        return self.embedder.encode(f"{_BGE_INSTRUCTION}{text}")

    # ------------------------------------------------------------------
    # Collection search helpers
    # ------------------------------------------------------------------

    def _search_all_collections(
        self,
        query_vector: List[float],
        query: AgentQuery,
        top_k: int,
        collections_filter: Optional[List[str]],
        year_min: Optional[int],
        year_max: Optional[int],
    ) -> List[SearchHit]:
        """Search all (or filtered) collections in parallel, apply weights."""
        target_collections = (
            [c for c in collections_filter if c in COLLECTION_CONFIG]
            if collections_filter
            else list(COLLECTION_CONFIG.keys())
        )

        all_hits: List[SearchHit] = []

        def _search_one(coll_name: str) -> List[SearchHit]:
            cfg = COLLECTION_CONFIG[coll_name]
            filters: Dict[str, Any] = {}

            # Gene filter
            if cfg["filter_field"] and hasattr(query, "gene") and query.gene:
                filters[cfg["filter_field"]] = query.gene

            # Year range filter
            year_field = cfg.get("year_field")
            if year_field:
                if year_min is not None:
                    filters[f"{year_field}__gte"] = year_min
                if year_max is not None:
                    filters[f"{year_field}__lte"] = year_max

            try:
                raw_hits = self.collection_manager.search(
                    collection=coll_name,
                    vector=query_vector,
                    top_k=top_k,
                    filters=filters if filters else None,
                )
            except Exception:
                logger.warning("Search failed for %s", coll_name, exc_info=True)
                return []

            weight = cfg["weight"]
            label = cfg["label"]

            # Convert raw dicts from collection manager to SearchHit objects
            search_hits: List[SearchHit] = []
            for raw in raw_hits:
                hit = SearchHit(
                    collection=coll_name,
                    id=str(raw.get("id", "")),
                    score=float(raw.get("_distance", 0.0)) * weight,
                    text=raw.get("text_summary", raw.get("text_chunk", "")),
                    metadata={k: v for k, v in raw.items() if k not in ("_distance", "_collection")},
                    label=label,
                )
                hit.citation = self._format_citation(coll_name, hit.record_id)
                hit.relevance = self._score_relevance(hit.score)
                search_hits.append(hit)
            return search_hits

        with ThreadPoolExecutor(max_workers=min(len(target_collections), 8)) as pool:
            futures = {
                pool.submit(_search_one, name): name
                for name in target_collections
            }
            for future in as_completed(futures):
                try:
                    all_hits.extend(future.result())
                except Exception:
                    logger.warning(
                        "Parallel search error for %s",
                        futures[future],
                        exc_info=True,
                    )

        return all_hits

    def _expanded_search(
        self,
        query: AgentQuery,
        top_k: int,
        collections_filter: Optional[List[str]],
        year_min: Optional[int],
        year_max: Optional[int],
    ) -> List[SearchHit]:
        """Run additional searches using query expansion terms."""
        expansion_terms = self.query_expander(query.text)
        if not expansion_terms:
            return []

        expanded_text = f"{query.text} {' '.join(expansion_terms)}"
        expanded_vector = self._embed_query(expanded_text)

        return self._search_all_collections(
            query_vector=expanded_vector,
            query=query,
            top_k=top_k,
            collections_filter=collections_filter,
            year_min=year_min,
            year_max=year_max,
        )

    # ------------------------------------------------------------------
    # Merge, rank, and de-duplicate
    # ------------------------------------------------------------------

    def _merge_and_rank(self, hits: List[SearchHit]) -> List[SearchHit]:
        """De-duplicate by ID, sort descending by score, cap at 30."""
        seen_ids: Set[str] = set()
        unique: List[SearchHit] = []
        for hit in hits:
            if hit.record_id not in seen_ids:
                seen_ids.add(hit.record_id)
                unique.append(hit)
        unique.sort(key=lambda h: h.score, reverse=True)
        return unique[:_MAX_EVIDENCE]

    @staticmethod
    def _score_relevance(score: float) -> str:
        """Classify a weighted score into high / medium / low relevance."""
        if score >= 0.85:
            return "high"
        if score >= 0.65:
            return "medium"
        return "low"

    # ------------------------------------------------------------------
    # Knowledge context injection
    # ------------------------------------------------------------------

    def _get_knowledge_context(self, query: AgentQuery) -> str:
        """Pull domain-knowledge snippets that match the query.

        Checks for gene mentions, therapy mentions, resistance mentions,
        pathway mentions, and biomarker mentions in the knowledge store.
        """
        if self.knowledge is None:
            return ""

        sections: List[str] = []
        query_lower = query.text.lower()

        # Gene mentions
        try:
            gene_context = self.knowledge.lookup_gene(query_lower)
            if gene_context:
                sections.append(f"[Gene Knowledge]\n{gene_context}")
        except AttributeError:
            pass  # Knowledge module doesn't implement this method
        except Exception as exc:
            logger.debug("Gene knowledge lookup failed: %s", exc)

        # Therapy mentions
        try:
            therapy_context = self.knowledge.lookup_therapy(query_lower)
            if therapy_context:
                sections.append(f"[Therapy Knowledge]\n{therapy_context}")
        except AttributeError:
            pass
        except Exception as exc:
            logger.debug("Therapy knowledge lookup failed: %s", exc)

        # Resistance mentions
        try:
            resistance_context = self.knowledge.lookup_resistance(query_lower)
            if resistance_context:
                sections.append(f"[Resistance Knowledge]\n{resistance_context}")
        except AttributeError:
            pass
        except Exception as exc:
            logger.debug("Resistance knowledge lookup failed: %s", exc)

        # Pathway mentions
        try:
            pathway_context = self.knowledge.lookup_pathway(query_lower)
            if pathway_context:
                sections.append(f"[Pathway Knowledge]\n{pathway_context}")
        except AttributeError:
            pass
        except Exception as exc:
            logger.debug("Pathway knowledge lookup failed: %s", exc)

        # Biomarker mentions
        try:
            biomarker_context = self.knowledge.lookup_biomarker(query_lower)
            if biomarker_context:
                sections.append(f"[Biomarker Knowledge]\n{biomarker_context}")
        except AttributeError:
            pass
        except Exception as exc:
            logger.debug("Biomarker knowledge lookup failed: %s", exc)

        return "\n\n".join(sections)

    # ------------------------------------------------------------------
    # Citation formatting
    # ------------------------------------------------------------------

    @staticmethod
    def _format_citation(collection: str, record_id: str) -> str:
        """Generate a clickable citation string for *record_id*.

        - PubMed IDs (PMID:*) link to pubmed.ncbi.nlm.nih.gov.
        - NCT IDs link to clinicaltrials.gov.
        - All others use a bracketed reference.
        """
        if record_id.startswith("PMID:"):
            pmid = record_id.replace("PMID:", "")
            return f"[PubMed {pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)"

        if record_id.upper().startswith("NCT"):
            nct = record_id.upper()
            return (
                f"[{nct}](https://clinicaltrials.gov/study/{nct})"
            )

        label = COLLECTION_CONFIG.get(collection, {}).get("label", collection)
        return f"[{label}: {record_id}]"

    # ------------------------------------------------------------------
    # Prompt construction
    # ------------------------------------------------------------------

    def _build_prompt(
        self,
        question: str,
        evidence: List[SearchHit],
    ) -> str:
        """Assemble the user prompt from evidence hits + knowledge context."""
        parts: List[str] = []

        # Knowledge context
        knowledge_ctx = self._get_knowledge_context(AgentQuery(question=question))
        if knowledge_ctx:
            parts.append("=== Domain Knowledge ===")
            parts.append(knowledge_ctx)
            parts.append("")

        # Retrieved evidence
        parts.append("=== Retrieved Evidence ===")
        for idx, hit in enumerate(evidence, 1):
            relevance_tag = f"[{hit.relevance}]" if hasattr(hit, "relevance") else ""
            label = getattr(hit, "label", None) or f"[{hit.collection}:{hit.id}]"
            citation = getattr(hit, "citation", None) or hit.id
            parts.append(
                f"{idx}. {label} {relevance_tag} "
                f"(score {hit.score:.3f}) — {citation}\n"
                f"   {hit.text}"
            )
        parts.append("")

        # User question
        parts.append("=== Question ===")
        parts.append(question)
        parts.append("")
        parts.append(
            "Using the evidence above, provide a thorough, well-cited answer. "
            "Include clickable reference links. If evidence is insufficient, "
            "state what is known and what remains uncertain."
        )
        return "\n".join(parts)

    # ------------------------------------------------------------------
    # Comparative retrieval
    # ------------------------------------------------------------------

    @staticmethod
    def _is_comparative(question: str) -> bool:
        """Return True when *question* contains comparative language."""
        return bool(_COMPARATIVE_RE.search(question))

    @staticmethod
    def _parse_comparison_entities(question: str) -> Tuple[str, str]:
        """Extract the two entities being compared.

        Handles patterns like:
        - ``A vs B``
        - ``A versus B``
        - ``compare A and B``
        - ``difference between A and B``
        """
        # "A vs B" / "A versus B"
        vs_match = re.search(
            r"(.+?)\s+(?:vs\.?|versus)\s+(.+?)(?:\?|$)",
            question,
            re.IGNORECASE,
        )
        if vs_match:
            return vs_match.group(1).strip(), vs_match.group(2).strip()

        # "compare A and B"
        cmp_match = re.search(
            r"compare\s+(.+?)\s+and\s+(.+?)(?:\?|$)",
            question,
            re.IGNORECASE,
        )
        if cmp_match:
            return cmp_match.group(1).strip(), cmp_match.group(2).strip()

        # "difference between A and B"
        diff_match = re.search(
            r"difference between\s+(.+?)\s+and\s+(.+?)(?:\?|$)",
            question,
            re.IGNORECASE,
        )
        if diff_match:
            return diff_match.group(1).strip(), diff_match.group(2).strip()

        return question, ""

    def retrieve_comparative(
        self,
        question: str,
        top_k: int = 10,
        collections_filter: Optional[List[str]] = None,
        year_min: Optional[int] = None,
        year_max: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Dual-entity retrieval for comparative questions.

        Returns
        -------
        dict
            Keys: ``entity_a``, ``entity_b``, ``hits_a``, ``hits_b``,
            ``shared_hits``.
        """
        entity_a, entity_b = self._parse_comparison_entities(question)

        query_a = AgentQuery(question=entity_a)
        query_b = AgentQuery(question=entity_b) if entity_b else None

        result_a = self.retrieve(
            query_a,
            top_k=top_k,
            collections_filter=collections_filter,
            year_min=year_min,
            year_max=year_max,
        )

        result_b = None
        if query_b:
            result_b = self.retrieve(
                query_b,
                top_k=top_k,
                collections_filter=collections_filter,
                year_min=year_min,
                year_max=year_max,
            )

        # Identify shared hits
        ids_a = {h.record_id for h in result_a.hits}
        ids_b = {h.record_id for h in (result_b.hits if result_b else [])}
        shared_ids = ids_a & ids_b
        shared_hits = [h for h in result_a.hits if h.record_id in shared_ids]

        return {
            "entity_a": entity_a,
            "entity_b": entity_b,
            "hits_a": result_a.hits,
            "hits_b": result_b.hits if result_b else [],
            "shared_hits": shared_hits,
        }

    def _build_comparative_prompt(
        self,
        question: str,
        comp: Dict[str, Any],
    ) -> str:
        """Build an LLM prompt with structured comparison sections."""
        parts: List[str] = []

        # Knowledge context
        knowledge_ctx = self._get_knowledge_context(AgentQuery(question=question))
        if knowledge_ctx:
            parts.append("=== Domain Knowledge ===")
            parts.append(knowledge_ctx)
            parts.append("")

        # Entity A evidence
        parts.append(f"=== Evidence for: {comp['entity_a']} ===")
        for idx, hit in enumerate(comp["hits_a"], 1):
            parts.append(
                f"{idx}. {hit.label} (score {hit.score:.3f}) — {hit.citation}\n"
                f"   {hit.text}"
            )
        parts.append("")

        # Entity B evidence
        if comp["entity_b"]:
            parts.append(f"=== Evidence for: {comp['entity_b']} ===")
            for idx, hit in enumerate(comp["hits_b"], 1):
                parts.append(
                    f"{idx}. {hit.label} (score {hit.score:.3f}) — {hit.citation}\n"
                    f"   {hit.text}"
                )
            parts.append("")

        # Shared evidence
        if comp["shared_hits"]:
            parts.append("=== Shared / Head-to-Head Evidence ===")
            for idx, hit in enumerate(comp["shared_hits"], 1):
                parts.append(
                    f"{idx}. {hit.label} (score {hit.score:.3f}) — {hit.citation}\n"
                    f"   {hit.text}"
                )
            parts.append("")

        # Comparison instructions
        parts.append("=== Question ===")
        parts.append(question)
        parts.append("")
        parts.append(
            "Provide a structured comparison addressing:\n"
            "1. Mechanism of action differences\n"
            "2. Efficacy data (ORR, PFS, OS where available)\n"
            "3. Safety / toxicity profile comparison\n"
            "4. Biomarker or patient-selection considerations\n"
            "5. Resistance mechanisms unique to each\n"
            "6. Guideline recommendations (NCCN/ESMO)\n"
            "7. Clinical trial evidence (cite specific trials)\n"
            "8. Summary recommendation with caveats\n\n"
            "Cite all evidence with clickable links. Acknowledge uncertainty "
            "where head-to-head data is lacking."
        )
        return "\n".join(parts)

    # ------------------------------------------------------------------
    # Comparative query handlers (non-streaming and streaming)
    # ------------------------------------------------------------------

    def _handle_comparative(self, question: str, **kwargs) -> str:
        """Execute a comparative RAG query and return the answer."""
        comp = self.retrieve_comparative(question, **kwargs)
        prompt = self._build_comparative_prompt(question, comp)
        messages = [
            {"role": "system", "content": ONCO_SYSTEM_PROMPT},
            {"role": "user", "content": prompt},
        ]
        return self.llm_client.chat(messages)

    def _handle_comparative_stream(
        self, question: str, **kwargs
    ) -> Generator[str, None, None]:
        """Streaming variant of :meth:`_handle_comparative`."""
        comp = self.retrieve_comparative(question, **kwargs)
        prompt = self._build_comparative_prompt(question, comp)
        messages = [
            {"role": "system", "content": ONCO_SYSTEM_PROMPT},
            {"role": "user", "content": prompt},
        ]
        yield from self.llm_client.chat_stream(messages)
