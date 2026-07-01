"""Multi-collection RAG engine for Pharmacogenomics Intelligence Agent.

Searches across all 15 collections simultaneously using parallel ThreadPoolExecutor,
synthesizes findings with pharmacogenomic knowledge augmentation (gene references,
drug guidelines, HLA screening, phenoconversion, dosing), and generates grounded
LLM responses.

Extends the pattern from: core/engines/precision-intelligence/src/rag_engine.py

Author: Adam Jones
Date: March 2026
"""

import logging
import re
import time
from typing import Dict, Generator, List, Optional

from config.settings import settings

from .models import (
    AgentQuery,
    CrossCollectionResult,
    PGxWorkflowType,
    SearchHit,
)

logger = logging.getLogger(__name__)

# Allowed characters for Milvus filter expressions to prevent injection
_SAFE_FILTER_RE = re.compile(r"^[A-Za-z0-9 _.\-/\*:]+$")

# ═══════════════════════════════════════════════════════════════════════
# SYSTEM PROMPT
# ═══════════════════════════════════════════════════════════════════════

PGX_SYSTEM_PROMPT = """You are a pharmacogenomics intelligence agent with deep expertise in:

1. **Pharmacogene Interpretation** — star allele nomenclature, haplotype resolution, defining variants, activity scores, PharmVar database curation
2. **Drug-Gene Interaction Analysis** — PharmGKB clinical annotations, variant-level drug response, PK/PD impact assessment, evidence level classification
3. **Star Allele Nomenclature** — CYP450 allele definitions (CYP2D6, CYP2C19, CYP2C9, CYP3A5), UGT, DPYD, TPMT, NUDT15 allele systems
4. **Diplotype-to-Phenotype Translation** — activity score summation, metabolizer status assignment (PM, IM, NM, RM, UM), CPIC standardized terminology
5. **CPIC/DPWG/FDA Guideline Application** — evidence-based prescribing recommendations, therapeutic alternatives, dose adjustments, alert levels
6. **HLA-Mediated Hypersensitivity Screening** — HLA-B*57:01 (abacavir), HLA-B*15:02 (carbamazepine), HLA-A*31:01, mandatory pre-treatment testing
7. **Phenoconversion Detection** — CYP inhibitor/inducer effects on genotype-predicted phenotype, clinical drug-drug-gene interactions
8. **Multi-Gene Interaction Modeling** — polygenic pharmacogenomic profiles, combined CYP2D6+CYP2C19 effects, multi-pathway metabolism
9. **Population Pharmacogenetics** — ancestry-specific allele frequencies, health equity implications, population-adjusted phenotype predictions
10. **Dosing Algorithm Application** — warfarin dosing algorithms (IWPC), tacrolimus CYP3A5-guided dosing, genotype-adjusted dose calculations
11. **Clinical Implementation** — preemptive vs reactive testing, EHR clinical decision support integration, pharmacist-led PGx programs, IGNITE network

You have access to evidence from MULTIPLE data sources spanning the entire pharmacogenomics knowledge domain,
from molecular variant definitions through clinical implementation programs.

When answering questions:
- **Cite evidence using clickable markdown links** provided in the evidence. Use the exact
  link format from the evidence, e.g. [Evidence:PMID 12345678](https://pubmed.ncbi.nlm.nih.gov/12345678/)
  or [Trial:NCT12345678](https://clinicaltrials.gov/study/NCT12345678). For GeneReference,
  DrugGuideline, DrugInteraction, HLA, Phenoconversion, DosingAlgorithm, PopulationData,
  FDALabel, DrugAlternative, PatientProfile, Implementation, and Education sources, use
  the format [Collection:record-id] (no URL needed).
- **Think cross-functionally** — connect insights across pharmacogene biology, clinical guidelines,
  population data, and implementation science
  (e.g., how allele frequencies affect guideline applicability across populations)
- **Highlight actionable findings** — distinguish required testing from informational PGx
- **Apply CPIC evidence levels** — clearly state the strength of recommendations (Level A, B, C, D)
- **Be specific** — cite gene-drug pairs (CYP2D6-codeine, CYP2C19-clopidogrel), star alleles (*1, *2, *17),
  and quantitative data (activity scores, allele frequencies) when available
- **Include regulatory context** — reference FDA labeling requirements and boxed warnings
- **Consider phenoconversion** — flag potential drug-drug-gene interactions that alter predicted phenotype
- **Acknowledge uncertainty** — distinguish validated PGx associations from emerging evidence

Your goal is to provide unified pharmacogenomic intelligence that enables precision prescribing
by integrating variant data, clinical guidelines, population genetics, and implementation science."""

# ═══════════════════════════════════════════════════════════════════════
# COLLECTION CONFIGURATION (reads weights from settings)
# ═══════════════════════════════════════════════════════════════════════

COLLECTION_CONFIG = {
    "pgx_gene_reference":       {"weight": settings.WEIGHT_GENE_REFERENCE,       "label": "GeneReference",    "has_gene_field": True,  "has_drug_field": False},
    "pgx_drug_guidelines":      {"weight": settings.WEIGHT_DRUG_GUIDELINES,      "label": "DrugGuideline",    "has_gene_field": True,  "has_drug_field": True},
    "pgx_drug_interactions":    {"weight": settings.WEIGHT_DRUG_INTERACTIONS,    "label": "DrugInteraction",  "has_gene_field": True,  "has_drug_field": True},
    "pgx_hla_hypersensitivity": {"weight": settings.WEIGHT_HLA_HYPERSENSITIVITY, "label": "HLA",             "has_gene_field": False, "has_drug_field": True},
    "pgx_phenoconversion":      {"weight": settings.WEIGHT_PHENOCONVERSION,      "label": "Phenoconversion",  "has_gene_field": False, "has_drug_field": False},
    "pgx_dosing_algorithms":    {"weight": settings.WEIGHT_DOSING_ALGORITHMS,    "label": "DosingAlgorithm",  "has_gene_field": False, "has_drug_field": True},
    "pgx_clinical_evidence":    {"weight": settings.WEIGHT_CLINICAL_EVIDENCE,    "label": "Evidence",         "has_gene_field": True,  "has_drug_field": True},
    "pgx_population_data":      {"weight": settings.WEIGHT_POPULATION_DATA,      "label": "PopulationData",   "has_gene_field": True,  "has_drug_field": False},
    "pgx_clinical_trials":      {"weight": settings.WEIGHT_CLINICAL_TRIALS,      "label": "Trial",            "has_gene_field": True,  "has_drug_field": True},
    "pgx_fda_labels":           {"weight": settings.WEIGHT_FDA_LABELS,           "label": "FDALabel",         "has_gene_field": True,  "has_drug_field": True},
    "pgx_drug_alternatives":    {"weight": settings.WEIGHT_DRUG_ALTERNATIVES,    "label": "DrugAlternative",  "has_gene_field": True,  "has_drug_field": False},
    "pgx_patient_profiles":     {"weight": settings.WEIGHT_PATIENT_PROFILES,     "label": "PatientProfile",   "has_gene_field": True,  "has_drug_field": False},
    "pgx_implementation":       {"weight": settings.WEIGHT_IMPLEMENTATION,       "label": "Implementation",   "has_gene_field": False, "has_drug_field": False},
    "pgx_education":            {"weight": settings.WEIGHT_EDUCATION,            "label": "Education",        "has_gene_field": False, "has_drug_field": False},
    "genomic_evidence":         {"weight": settings.WEIGHT_GENOMIC,              "label": "Genomic",          "has_gene_field": False, "has_drug_field": False},
}

# ═══════════════════════════════════════════════════════════════════════
# WORKFLOW → COLLECTION BOOST MAPPING
# ═══════════════════════════════════════════════════════════════════════

WORKFLOW_COLLECTION_BOOST: Dict[PGxWorkflowType, List[str]] = {
    PGxWorkflowType.GENE_QUERY: [
        "pgx_gene_reference", "pgx_drug_guidelines",
    ],
    PGxWorkflowType.DRUG_QUERY: [
        "pgx_drug_guidelines", "pgx_drug_interactions", "pgx_drug_alternatives",
    ],
    PGxWorkflowType.PROFILE_QUERY: [
        "pgx_gene_reference", "pgx_drug_guidelines", "pgx_patient_profiles",
        "pgx_drug_interactions", "pgx_hla_hypersensitivity", "pgx_phenoconversion",
        "pgx_dosing_algorithms", "pgx_clinical_evidence", "pgx_population_data",
        "pgx_clinical_trials", "pgx_fda_labels", "pgx_drug_alternatives",
        "pgx_implementation", "pgx_education", "genomic_evidence",
    ],
    PGxWorkflowType.INTERACTION_QUERY: [
        "pgx_phenoconversion", "pgx_drug_interactions",
    ],
    PGxWorkflowType.DOSING_QUERY: [
        "pgx_dosing_algorithms", "pgx_drug_guidelines",
    ],
    PGxWorkflowType.HLA_SCREEN: [
        "pgx_hla_hypersensitivity",
    ],
}

# Known pharmacogenes for filter matching
_KNOWN_GENES = {
    "CYP2D6", "CYP2C19", "CYP2C9", "CYP3A5", "CYP3A4", "CYP2B6",
    "CYP1A2", "CYP2A6", "CYP2E1", "CYP4F2",
    "DPYD", "TPMT", "NUDT15", "UGT1A1", "SLCO1B1", "VKORC1",
    "NAT2", "G6PD", "IFNL3", "CFTR", "RYR1", "CACNA1S",
    "HLA-A", "HLA-B", "HLA-C", "HLA-DRB1",
}


class PGxRAGEngine:
    """Multi-collection RAG engine for pharmacogenomics cross-functional queries.

    Searches across all PGx collections simultaneously using parallel
    ThreadPoolExecutor, merges results with knowledge graph context, and
    generates grounded LLM responses.

    Features:
    - Parallel search via ThreadPoolExecutor (all 15 collections)
    - Settings-driven weights and parameters
    - Workflow-based dynamic weight boosting (GENE_QUERY, DRUG_QUERY, etc.)
    - Gene/drug field-based filtering
    - Citation relevance scoring (high/medium/low)
    - Cross-collection entity linking
    - Comparative analysis support
    - Conversation memory context injection
    """

    def __init__(self, collection_manager, embedder, llm_client,
                 knowledge=None, query_expander=None):
        self.collections = collection_manager
        self.embedder = embedder
        self.llm = llm_client
        self.knowledge = knowledge
        self.expander = query_expander

    def _compute_boosted_weights(self, workflows: List[PGxWorkflowType]) -> Dict[str, float]:
        """Compute adjusted collection weights based on relevant PGx workflows.

        Collections associated with the given workflows receive a 1.5x weight
        boost. All weights are then renormalized to sum to 1.0.

        Args:
            workflows: List of PGxWorkflowType values identified by query classification.

        Returns:
            Dict mapping collection name to adjusted weight (summing to ~1.0).
        """
        # Start with base weights from COLLECTION_CONFIG
        weights = {
            name: cfg["weight"]
            for name, cfg in COLLECTION_CONFIG.items()
        }

        # Determine which collections to boost
        boosted_collections: set = set()
        for workflow in workflows:
            for coll in WORKFLOW_COLLECTION_BOOST.get(workflow, []):
                boosted_collections.add(coll)

        # Apply 1.5x boost
        for coll in boosted_collections:
            if coll in weights:
                weights[coll] *= 1.5

        # Renormalize to sum to 1.0
        total = sum(weights.values())
        if total > 0:
            weights = {name: w / total for name, w in weights.items()}

        return weights

    def _classify_workflow(self, query: str) -> List[PGxWorkflowType]:
        """Classify a query into one or more PGxWorkflowType categories based on keywords.

        Args:
            query: The user's natural language question.

        Returns:
            List of PGxWorkflowType values matching the query intent.
        """
        q_upper = query.upper()
        workflows = []

        # HLA screening keywords
        hla_keywords = ["HLA", "HYPERSENSITIVITY", "SJS", "TEN", "DRESS",
                         "ABACAVIR", "CARBAMAZEPINE", "ALLOPURINOL",
                         "HLA-B*57:01", "HLA-B*15:02", "HLA-A*31:01"]
        if any(kw in q_upper for kw in hla_keywords):
            workflows.append(PGxWorkflowType.HLA_SCREEN)

        # Dosing algorithm keywords
        dosing_keywords = ["DOSE", "DOSING", "ALGORITHM", "DOSE ADJUST",
                           "DOSE CALCULATION", "DOSE REDUCTION", "WARFARIN DOSE",
                           "TACROLIMUS DOSE", "IWPC"]
        if any(kw in q_upper for kw in dosing_keywords):
            workflows.append(PGxWorkflowType.DOSING_QUERY)

        # Interaction / phenoconversion keywords
        interaction_keywords = ["PHENOCONVERSION", "DRUG-DRUG-GENE", "DDI",
                                "DRUG INTERACTION", "INHIBITOR", "INDUCER",
                                "PRECIPITANT", "SUBSTRATE"]
        if any(kw in q_upper for kw in interaction_keywords):
            workflows.append(PGxWorkflowType.INTERACTION_QUERY)

        # Drug-centric keywords
        drug_keywords = ["DRUG", "MEDICATION", "PRESCRI", "ALTERNATIVE",
                         "SUBSTITUTE", "THERAPEUTIC", "CONTRAINDIC"]
        if any(kw in q_upper for kw in drug_keywords):
            workflows.append(PGxWorkflowType.DRUG_QUERY)

        # Gene-centric keywords (check for known gene symbols)
        for gene in _KNOWN_GENES:
            if gene in q_upper:
                workflows.append(PGxWorkflowType.GENE_QUERY)
                break

        # Star allele pattern (*1, *2, *17, etc.)
        if re.search(r'\*\d+', query):
            if PGxWorkflowType.GENE_QUERY not in workflows:
                workflows.append(PGxWorkflowType.GENE_QUERY)

        # Profile-level keywords
        profile_keywords = ["PATIENT", "PROFILE", "DIPLOTYPE", "PHENOTYPE",
                            "METABOLIZER", "COMPREHENSIVE", "FULL PANEL",
                            "PANEL RESULT"]
        if any(kw in q_upper for kw in profile_keywords):
            workflows.append(PGxWorkflowType.PROFILE_QUERY)

        # Default to profile query if nothing matched
        if not workflows:
            workflows.append(PGxWorkflowType.PROFILE_QUERY)

        return workflows

    def retrieve(self, query: AgentQuery,
                 top_k_per_collection: int = None,
                 collections_filter: List[str] = None,
                 conversation_context: str = None,
                 workflows: List[PGxWorkflowType] = None) -> CrossCollectionResult:
        """Retrieve evidence from collections for a query.

        Args:
            query: The agent query with question and optional filters
            top_k_per_collection: Max results per collection (default from settings)
            collections_filter: Optional list of collection names to search
            conversation_context: Optional prior conversation context for follow-ups
            workflows: Optional list of PGxWorkflowType values for dynamic weight boosting
        """
        top_k = top_k_per_collection or settings.TOP_K_PER_COLLECTION
        start = time.time()

        # Optionally prepend conversation context for follow-up queries
        search_text = query.question
        if conversation_context:
            search_text = f"{conversation_context}\n\nCurrent question: {query.question}"

        # Step 1: Embed query
        query_embedding = self._embed_query(search_text)

        # Step 2: Determine collections to search
        collections_to_search = collections_filter or list(COLLECTION_CONFIG.keys())

        # Step 3: Build per-collection filters
        filter_exprs = {}
        for coll in collections_to_search:
            parts = []
            cfg = COLLECTION_CONFIG.get(coll, {})
            if query.gene and cfg.get("has_gene_field"):
                safe_gene = query.gene.strip()
                if _SAFE_FILTER_RE.match(safe_gene):
                    parts.append(f'gene == "{safe_gene}"')
                else:
                    logger.warning("Rejected unsafe gene filter value: %r", safe_gene)
            if query.drug and cfg.get("has_drug_field"):
                safe_drug = query.drug.strip()
                if _SAFE_FILTER_RE.match(safe_drug):
                    parts.append(f'drug == "{safe_drug}"')
                else:
                    logger.warning("Rejected unsafe drug filter value: %r", safe_drug)
            if parts:
                filter_exprs[coll] = " and ".join(parts)

        # Step 3b: Auto-classify workflow if not provided
        if not workflows:
            workflows = self._classify_workflow(query.question)

        # Step 3c: Compute boosted weights
        boosted_weights = self._compute_boosted_weights(workflows)

        # Step 4: Parallel search across all collections
        all_hits = self._search_all_collections(
            query_embedding, collections_to_search, top_k, filter_exprs,
            weight_overrides=boosted_weights,
        )

        # Step 5: Query expansion (semantic search, not field-filter)
        if self.expander:
            expanded_hits = self._expanded_search(
                query.question, query_embedding, collections_to_search, top_k,
            )
            all_hits.extend(expanded_hits)

        # Step 6: Deduplicate, score citations, rank
        hits = self._merge_and_rank(all_hits)

        # Step 7: Knowledge graph augmentation
        knowledge_context = ""
        if self.knowledge:
            knowledge_context = self._get_knowledge_context(query.question)

        elapsed = (time.time() - start) * 1000

        return CrossCollectionResult(
            query=query.question,
            hits=hits,
            knowledge_context=knowledge_context,
            total_collections_searched=len(collections_to_search),
            search_time_ms=elapsed,
        )

    def query(self, question: str, **kwargs) -> str:
        """Full RAG query: retrieve evidence + generate LLM response."""
        agent_query = AgentQuery(question=question, **kwargs)
        evidence = self.retrieve(agent_query)
        prompt = self._build_prompt(agent_query.question, evidence)
        return self.llm.generate(
            prompt=prompt,
            system_prompt=PGX_SYSTEM_PROMPT,
            max_tokens=2048,
            temperature=0.7,
        )

    def query_stream(self, question: str,
                     **kwargs) -> Generator[Dict, None, None]:
        """Streaming RAG query -- yields evidence then token chunks."""
        agent_query = AgentQuery(question=question, **kwargs)
        evidence = self.retrieve(agent_query)
        yield {"type": "evidence", "content": evidence}

        prompt = self._build_prompt(agent_query.question, evidence)
        full_answer = ""
        for token in self.llm.generate_stream(
            prompt=prompt,
            system_prompt=PGX_SYSTEM_PROMPT,
            max_tokens=2048,
            temperature=0.7,
        ):
            full_answer += token
            yield {"type": "token", "content": token}
        yield {"type": "done", "content": full_answer}

    # ── Cross-Collection Entity Linking ─────────────────────────────

    def find_related(self, entity: str, top_k: int = 5) -> Dict[str, List[SearchHit]]:
        """Find all evidence related to an entity across all 15 collections.

        Enables "show me everything about CYP2D6" spanning all collections.

        Args:
            entity: Gene name, drug name, star allele, etc.
            top_k: Max results per collection

        Returns:
            Dict of collection_name -> List[SearchHit]
        """
        embedding = self._embed_query(entity)
        results = {}

        all_results = self.collections.search_all(
            embedding, top_k_per_collection=top_k,
            score_threshold=settings.SCORE_THRESHOLD,
        )
        for coll_name, hits in all_results.items():
            label = COLLECTION_CONFIG.get(coll_name, {}).get("label", coll_name)
            search_hits = []
            for r in hits:
                hit = SearchHit(
                    collection=label,
                    id=r.get("id", ""),
                    score=r.get("score", 0.0),
                    text=r.get("text_summary", r.get("text_chunk", "")),
                    metadata=r,
                )
                search_hits.append(hit)
            if search_hits:
                results[coll_name] = search_hits
        return results

    # ── Private Methods ──────────────────────────────────────────────

    def _embed_query(self, text: str):
        """Embed query text with BGE instruction prefix."""
        prefix = "Represent this sentence for searching relevant passages: "
        return self.embedder.embed_text(prefix + text)

    def _search_all_collections(
        self, query_embedding, collections: List[str],
        top_k: int, filter_exprs: Dict[str, str],
        weight_overrides: Dict[str, float] = None,
    ) -> List[SearchHit]:
        """Search all collections in parallel via ThreadPoolExecutor."""
        all_hits = []

        # Use the parallel search_all method from PGxCollectionManager
        parallel_results = self.collections.search_all(
            query_embedding,
            top_k_per_collection=top_k,
            filter_exprs=filter_exprs,
            score_threshold=settings.SCORE_THRESHOLD,
        )

        for coll_name, results in parallel_results.items():
            if coll_name not in collections:
                continue
            cfg = COLLECTION_CONFIG.get(coll_name, {})
            if weight_overrides and coll_name in weight_overrides:
                weight = weight_overrides[coll_name]
            else:
                weight = cfg.get("weight", 0.1)
            label = cfg.get("label", coll_name)

            for r in results:
                raw_score = r.get("score", 0.0)
                weighted_score = raw_score * (1 + weight)

                # Citation relevance scoring
                if raw_score >= settings.CITATION_HIGH_THRESHOLD:
                    relevance = "high"
                elif raw_score >= settings.CITATION_MEDIUM_THRESHOLD:
                    relevance = "medium"
                else:
                    relevance = "low"

                metadata = {k: v for k, v in r.items() if k not in ("embedding",)}
                metadata["relevance"] = relevance

                hit = SearchHit(
                    collection=label,
                    id=r.get("id", ""),
                    score=weighted_score,
                    text=r.get("text_summary", r.get("text_chunk", "")),
                    metadata=metadata,
                )
                all_hits.append(hit)

        return all_hits

    def _expanded_search(
        self, query: str, query_embedding,
        collections: List[str], top_k: int,
    ) -> List[SearchHit]:
        """Use query expansion for additional coverage.

        Expansion terms that are pharmacogenes use field filters.
        Non-gene terms are re-embedded for semantic search across all collections.
        """
        if not self.expander:
            return []

        from .query_expansion import expand_query
        expanded_terms = expand_query(query)

        additional_hits = []
        for term in expanded_terms[:5]:
            term_upper = term.upper().replace("-", "").replace(" ", "")

            # Check if this term is a known pharmacogene
            is_gene = any(
                term_upper == g.upper().replace("-", "").replace(" ", "")
                for g in _KNOWN_GENES
            )

            if is_gene:
                # Use as field filter on gene-capable collections
                safe_term = term.strip()
                if not _SAFE_FILTER_RE.match(safe_term):
                    continue
                for coll_name in collections:
                    if not COLLECTION_CONFIG.get(coll_name, {}).get("has_gene_field"):
                        continue
                    try:
                        filter_expr = f'gene == "{safe_term}"'
                        results = self.collections.search(
                            coll_name, query_embedding, min(3, top_k), filter_expr,
                        )
                        label = COLLECTION_CONFIG.get(coll_name, {}).get("label", coll_name)
                        for r in results:
                            additional_hits.append(SearchHit(
                                collection=label,
                                id=r.get("id", ""),
                                score=r.get("score", 0.0) * 0.8,
                                text=r.get("text_summary", r.get("text_chunk", "")),
                                metadata=r,
                            ))
                    except Exception as exc:
                        logger.warning("Expanded gene search failed for %s/%s: %s", coll_name, safe_term, exc)
            else:
                # Semantic search: re-embed the expansion term and search all collections
                try:
                    term_embedding = self._embed_query(term)
                    term_results = self.collections.search_all(
                        term_embedding, top_k_per_collection=2,
                        score_threshold=settings.SCORE_THRESHOLD,
                    )
                    for coll_name, results in term_results.items():
                        if coll_name not in collections:
                            continue
                        label = COLLECTION_CONFIG.get(coll_name, {}).get("label", coll_name)
                        for r in results:
                            additional_hits.append(SearchHit(
                                collection=label,
                                id=r.get("id", ""),
                                score=r.get("score", 0.0) * 0.7,
                                text=r.get("text_summary", r.get("text_chunk", "")),
                                metadata=r,
                            ))
                except Exception as exc:
                    logger.warning("Expanded semantic search failed for '%s': %s", term[:50], exc)

        return additional_hits

    def _merge_and_rank(self, hits: List[SearchHit]) -> List[SearchHit]:
        """Deduplicate by ID, sort by score descending, cap at 30."""
        seen = set()
        unique = []
        for hit in hits:
            if hit.id not in seen:
                seen.add(hit.id)
                unique.append(hit)
        unique.sort(key=lambda h: h.score, reverse=True)
        return unique[:30]

    def _get_knowledge_context(self, query: str) -> str:
        """Extract knowledge graph context from PGx domains."""
        if not self.knowledge:
            return ""

        from .knowledge import (
            get_gene_context,
            get_drug_context,
            get_drug_category_context,
            get_hla_context,
            get_inhibitor_context,
        )

        context_parts = []
        query_upper = query.upper()

        # Check for pharmacogene mentions
        for gene in _KNOWN_GENES:
            if gene in query_upper:
                ctx = get_gene_context(gene)
                if ctx:
                    context_parts.append(ctx)

        # Check for HLA mentions
        hla_keywords = {
            "HLA-B*57:01": "HLA-B*57:01", "ABACAVIR": "HLA-B*57:01",
            "HLA-B*15:02": "HLA-B*15:02", "CARBAMAZEPINE": "HLA-B*15:02",
            "HLA-A*31:01": "HLA-A*31:01", "OXCARBAZEPINE": "HLA-B*15:02",
            "HLA-B*58:01": "HLA-B*58:01", "ALLOPURINOL": "HLA-B*58:01",
        }
        for keyword, hla_allele in hla_keywords.items():
            if keyword in query_upper:
                ctx = get_hla_context(hla_allele)
                if ctx:
                    context_parts.append(ctx)
                break

        # Check for phenoconversion mentions
        phenoconv_keywords = {
            "PHENOCONVERSION": "CYP2D6",
            "FLUOXETINE": "CYP2D6",
            "PAROXETINE": "CYP2D6",
            "BUPROPION": "CYP2B6",
            "INHIBITOR": "CYP2D6",
            "INDUCER": "CYP3A4",
        }
        for keyword, enzyme in phenoconv_keywords.items():
            if keyword in query_upper:
                ctx = get_inhibitor_context(enzyme)
                if ctx:
                    context_parts.append(ctx)
                break

        # Check for drug-specific guideline mentions
        drug_keywords = {
            "CODEINE": "codeine", "TRAMADOL": "tramadol",
            "CLOPIDOGREL": "clopidogrel", "VORICONAZOLE": "voriconazole",
            "WARFARIN": "warfarin", "TAMOXIFEN": "tamoxifen",
            "FLUOROURACIL": "fluorouracil", "5-FU": "fluorouracil",
            "CAPECITABINE": "capecitabine", "THIOPURINE": "thiopurine",
            "AZATHIOPRINE": "azathioprine", "MERCAPTOPURINE": "mercaptopurine",
            "TACROLIMUS": "tacrolimus", "SIMVASTATIN": "simvastatin",
            "IRINOTECAN": "irinotecan", "ATOMOXETINE": "atomoxetine",
        }
        for keyword, drug in drug_keywords.items():
            if keyword in query_upper:
                ctx = get_drug_context(drug)
                if ctx:
                    context_parts.append(ctx)

        # Check for general guideline terms
        if any(term in query_upper for term in ["CPIC", "DPWG", "GUIDELINE", "RECOMMENDATION"]):
            if not any("Guideline" in p for p in context_parts):
                ctx = get_drug_category_context("CPIC")
                if ctx:
                    context_parts.append(ctx)

        return "\n\n".join(context_parts)

    @staticmethod
    def _format_citation(collection: str, record_id: str) -> str:
        """Format a citation with clickable URL where possible."""
        if collection == "Evidence" and record_id.isdigit():
            return (
                f"[Evidence:PMID {record_id}]"
                f"(https://pubmed.ncbi.nlm.nih.gov/{record_id}/)"
            )
        if collection == "Trial" and record_id.upper().startswith("NCT"):
            return (
                f"[Trial:{record_id}]"
                f"(https://clinicaltrials.gov/study/{record_id})"
            )
        return f"[{collection}:{record_id}]"

    def _build_prompt(self, question: str,
                      evidence: CrossCollectionResult) -> str:
        """Build LLM prompt with evidence, knowledge context, and relevance tags."""
        sections = []
        by_coll = evidence.hits_by_collection()

        for coll_name, hits in by_coll.items():
            section_lines = [f"### Evidence from {coll_name}"]
            for i, hit in enumerate(hits[:5], 1):
                citation = self._format_citation(hit.collection, hit.id)
                relevance = hit.metadata.get("relevance", "")
                relevance_tag = f" [{relevance} relevance]" if relevance else ""
                section_lines.append(
                    f"{i}. {citation}{relevance_tag} "
                    f"(score={hit.score:.3f}) {hit.text[:500]}"
                )
            sections.append("\n".join(section_lines))

        evidence_text = "\n\n".join(sections) if sections else "No evidence found."

        knowledge_text = ""
        if evidence.knowledge_context:
            knowledge_text = (
                f"\n\n### Knowledge Graph Context\n"
                f"{evidence.knowledge_context}"
            )

        return (
            f"## Retrieved Evidence\n\n"
            f"{evidence_text}"
            f"{knowledge_text}\n\n"
            f"---\n\n"
            f"## Question\n\n"
            f"{question}\n\n"
            f"Please provide a comprehensive answer grounded in the evidence above. "
            f"Cite sources using the clickable markdown links provided in each evidence item. "
            f"Prioritize [high relevance] citations. "
            f"Consider cross-functional insights across pharmacogene biology, clinical guidelines, "
            f"population genetics, and implementation science."
        )

    # ── Comparative Analysis Methods ────────────────────────────────

    def _is_comparative(self, question: str) -> bool:
        q_upper = question.upper()
        return ("COMPARE" in q_upper or " VS " in q_upper
                or "VERSUS" in q_upper or "COMPARING" in q_upper)

    def retrieve_comparative(self, question: str,
                             entity_a_name: str,
                             entity_b_name: str,
                             collections_filter: List[str] = None) -> Optional['CrossCollectionResult']:
        """Run comparative retrieval with two separate searches.

        Useful for comparisons like "CYP2D6 vs CYP2C19" or "codeine vs tramadol".

        Args:
            question: The comparative question text.
            entity_a_name: First entity (gene or drug name).
            entity_b_name: Second entity (gene or drug name).
            collections_filter: Optional list of collections to limit search.

        Returns:
            ComparativeResult with evidence for both entities, or None if parsing fails.
        """
        from .models import ComparativeResult

        start = time.time()

        # Determine if entities are genes or drugs for filtering
        gene_a = entity_a_name.upper() if entity_a_name.upper() in _KNOWN_GENES else None
        gene_b = entity_b_name.upper() if entity_b_name.upper() in _KNOWN_GENES else None

        query_a = AgentQuery(question=question, gene=gene_a)
        query_b = AgentQuery(question=question, gene=gene_b)

        evidence_a = self.retrieve(query_a, collections_filter=collections_filter)
        evidence_b = self.retrieve(query_b, collections_filter=collections_filter)

        comparison_context = ""
        if self.knowledge:
            from .knowledge import get_comparison_context, resolve_comparison_entity
            entity_a_resolved = resolve_comparison_entity(entity_a_name)
            entity_b_resolved = resolve_comparison_entity(entity_b_name)
            if entity_a_resolved and entity_b_resolved:
                comparison_context = get_comparison_context(entity_a_resolved, entity_b_resolved)

        elapsed = (time.time() - start) * 1000

        return ComparativeResult(
            query=question,
            entity_a=entity_a_name,
            entity_b=entity_b_name,
            evidence_a=evidence_a,
            evidence_b=evidence_b,
            comparison_context=comparison_context,
            total_search_time_ms=elapsed,
        )

    def _build_comparative_prompt(self, question: str, comp) -> str:
        def _fmt(label: str, evidence) -> str:
            sections = []
            by_coll = evidence.hits_by_collection()
            for coll_name, hits in by_coll.items():
                lines = [f"#### {coll_name}"]
                for i, hit in enumerate(hits[:4], 1):
                    citation = self._format_citation(hit.collection, hit.id)
                    lines.append(
                        f"{i}. {citation} "
                        f"(score={hit.score:.3f}) {hit.text[:400]}"
                    )
                sections.append("\n".join(lines))
            if not sections:
                return f"### Evidence for {label}\nNo evidence found."
            return f"### Evidence for {label}\n\n" + "\n\n".join(sections)

        evidence_a_text = _fmt(comp.entity_a, comp.evidence_a)
        evidence_b_text = _fmt(comp.entity_b, comp.evidence_b)

        knowledge_text = ""
        if comp.comparison_context:
            knowledge_text = (
                f"\n\n### Knowledge Graph Comparison Data\n"
                f"{comp.comparison_context}"
            )

        return (
            f"## Comparative Analysis Evidence\n\n"
            f"{evidence_a_text}\n\n"
            f"---\n\n"
            f"{evidence_b_text}"
            f"{knowledge_text}\n\n"
            f"---\n\n"
            f"## Question\n\n{question}\n\n"
            f"## Instructions\n\n"
            f"Provide a structured comparison of **{comp.entity_a}** vs "
            f"**{comp.entity_b}**. Your response MUST include:\n\n"
            f"1. A **comparison table** in markdown format with key dimensions "
            f"as rows and the two entities as columns.\n"
            f"2. **Key pharmacogenomic differences** (bulleted list).\n"
            f"3. **Clinical implications** of each (bulleted list).\n"
            f"4. A **clinical context** paragraph explaining when each is most "
            f"relevant for pharmacogenomic testing.\n\n"
            f"Cite sources using the clickable markdown links provided in "
            f"the evidence above."
        )
