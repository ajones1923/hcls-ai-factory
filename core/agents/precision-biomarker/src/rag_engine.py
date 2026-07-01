"""Multi-collection RAG engine for Biomarker Intelligence Agent.

Searches across all 14 collections simultaneously using parallel ThreadPoolExecutor,
synthesizes findings with full knowledge graph augmentation (biomarker domains, PGx,
PhenoAge, cross-modal links), and generates grounded LLM responses.

Extends the pattern from: cart_intelligence_agent/src/rag_engine.py

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
    SearchHit,
)

logger = logging.getLogger(__name__)

# Allowed characters for Milvus filter expressions to prevent injection
_SAFE_FILTER_RE = re.compile(r"^[A-Za-z0-9 _\-]+$")

# =====================================================================
# SYSTEM PROMPT
# =====================================================================

BIOMARKER_SYSTEM_PROMPT = """You are a Biomarker Intelligence Agent with deep expertise in:

1. **Biological Aging** -- PhenoAge (Levine 2018), GrimAge (Lu 2019), epigenetic clocks, and aging biomarker interpretation
2. **Pre-Symptomatic Disease Detection** -- trajectory analysis for diabetes, cardiovascular, liver, thyroid, iron, and nutritional conditions years before clinical diagnosis
3. **Pharmacogenomic Drug-Gene Interactions** -- CPIC guidelines, CYP2D6/CYP2C19/CYP2C9/VKORC1/SLCO1B1/DPYD/CYP3A5 star allele interpretation, metabolizer phenotype classification
4. **Genotype-Adjusted Reference Ranges** -- MTHFR, APOE, PNPLA3, HFE, DIO2, VDR, FADS1 impact on biomarker interpretation
5. **Nutritional Genomics** -- MTHFR/methylfolate, FADS1/omega-3 conversion, VDR/vitamin D, BCMO1/vitamin A, FUT2/B12 absorption
6. **Cardiovascular Risk Stratification** -- Lp(a), ApoB, APOE genotyping, PCSK9 variants, homocysteine/MTHFR
7. **Liver Health Assessment** -- PNPLA3 I148M, TM6SF2, FIB-4 index, NAFLD/NASH risk stratification
8. **Iron Metabolism** -- HFE C282Y/H63D genotyping, hereditary hemochromatosis, TMPRSS6 and iron-refractory anemia
9. **Ashkenazi Jewish Carrier Screening** -- 10-gene AJ panel (BRCA1/2, GBA, HEXA, FANCC, ASPA, BLM, SMPD1, IKBKAP/ELP1, MCOLN1), compound risks (GBA+APOE E4 for Parkinson's), reproductive counseling

You have access to evidence from MULTIPLE data sources spanning:
- Biomarker reference databases
- Genetic variant-biomarker interactions
- Pharmacogenomic dosing rules (CPIC)
- Disease trajectory models
- Clinical literature evidence
- Nutritional genomics guidelines
- Drug interaction databases
- Aging clock marker data
- Genotype-adjusted reference ranges
- Monitoring protocols
- Patient genomic evidence (shared VCF-derived collection)
- AJ carrier screening genetic data

When answering questions:
- **Cite evidence using the collection labels** provided in the evidence, e.g. [BiomarkerRef:marker-id],
  [ClinicalEvidence:PMID 12345678](https://pubmed.ncbi.nlm.nih.gov/12345678/),
  [PGxRule:gene-drug-id], etc.
- **Always specify units** when discussing biomarker values (mg/dL, ng/mL, %, etc.)
- **Provide genotype-specific interpretation** when patient genotype data is available
- **Highlight critical findings** (PGx drug safety, severe iron overload, accelerated aging) prominently
- **Recommend actionable interventions** grounded in evidence (CPIC level, clinical trial data)
- **Explain pre-symptomatic disease trajectories** with estimated timelines
- **Integrate across domains** -- connect biomarker, genetic, and pharmacogenomic findings
- **Acknowledge uncertainty** -- distinguish validated associations from emerging data
- **Flag cross-modal triggers** when findings warrant imaging, oncology, or genomic pipeline follow-up

Your goal is to transform raw biomarker and genetic data into actionable precision health
intelligence that enables truly personalized medicine."""

# =====================================================================
# COLLECTION CONFIGURATION (reads weights from settings)
# =====================================================================

COLLECTION_CONFIG = {
    "biomarker_reference":        {"weight": settings.WEIGHT_BIOMARKER_REF,        "label": "BiomarkerRef",       "has_disease_area": False, "year_field": None},
    "biomarker_genetic_variants": {"weight": settings.WEIGHT_GENETIC_VARIANTS,     "label": "GeneticVariant",     "has_disease_area": True,  "year_field": None},
    "biomarker_pgx_rules":        {"weight": settings.WEIGHT_PGX_RULES,            "label": "PGxRule",            "has_disease_area": False, "year_field": None},
    "biomarker_disease_trajectories": {"weight": settings.WEIGHT_DISEASE_TRAJECTORIES, "label": "DiseaseTrajectory", "has_disease_area": True,  "year_field": None},
    "biomarker_clinical_evidence": {"weight": settings.WEIGHT_CLINICAL_EVIDENCE,   "label": "ClinicalEvidence",   "has_disease_area": True,  "year_field": "year"},
    "biomarker_nutrition":        {"weight": settings.WEIGHT_NUTRITION,             "label": "Nutrition",          "has_disease_area": False, "year_field": None},
    "biomarker_drug_interactions": {"weight": settings.WEIGHT_DRUG_INTERACTIONS,   "label": "DrugInteraction",    "has_disease_area": False, "year_field": None},
    "biomarker_aging_markers":    {"weight": settings.WEIGHT_AGING_MARKERS,        "label": "AgingMarker",        "has_disease_area": False, "year_field": None},
    "biomarker_genotype_adjustments": {"weight": settings.WEIGHT_GENOTYPE_ADJUSTMENTS, "label": "GenotypeAdj",   "has_disease_area": False, "year_field": None},
    "biomarker_monitoring":       {"weight": settings.WEIGHT_MONITORING,            "label": "Monitoring",         "has_disease_area": False, "year_field": None},
    "biomarker_critical_values":     {"weight": settings.WEIGHT_CRITICAL_VALUES,     "label": "CriticalValue",    "has_disease_area": False, "year_field": None},
    "biomarker_discordance_rules":   {"weight": settings.WEIGHT_DISCORDANCE_RULES,   "label": "DiscordanceRule",  "has_disease_area": False, "year_field": None},
    "biomarker_aj_carrier_screening": {"weight": settings.WEIGHT_AJ_CARRIER_SCREENING, "label": "AJCarrierScreen", "has_disease_area": False, "year_field": None},
    "genomic_evidence":           {"weight": settings.WEIGHT_GENOMIC_EVIDENCE,      "label": "Genomic",            "has_disease_area": False, "year_field": None},
}

# Maximum merged results returned after deduplication and ranking
MAX_MERGED_RESULTS = 30

# Maximum evidence items per collection included in the LLM prompt
MAX_PROMPT_EVIDENCE = 5

# Known disease areas for filter expressions
_KNOWN_DISEASE_AREAS = {
    "diabetes", "cardiovascular", "liver", "thyroid", "iron", "nutritional",
    "metabolic", "renal", "hematologic", "oncology",
}


class BiomarkerRAGEngine:
    """Multi-collection RAG engine for precision biomarker queries.

    Searches across all biomarker collections simultaneously using parallel
    ThreadPoolExecutor, merges results with knowledge graph context, and
    generates grounded LLM responses.

    Features:
    - Parallel search via ThreadPoolExecutor (delegated to collection manager)
    - Settings-driven weights and parameters
    - Full knowledge graph augmentation (domains, PGx, PhenoAge, cross-modal)
    - Citation relevance scoring (high/medium/low)
    - Disease-area filtering on applicable collections
    - Temporal date-range filtering for clinical evidence
    - Cross-collection entity linking
    - Patient profile context injection in prompts
    - Conversation memory context injection
    """

    def __init__(self, collection_manager, embedder, llm_client,
                 knowledge=None):
        """Initialize the RAG engine.

        Args:
            collection_manager: Milvus collection manager with search_all() method.
            embedder: Embedding model with embed_text() method.
            llm_client: LLM client with generate() and generate_stream() methods.
            knowledge: Knowledge module (src.knowledge) for context augmentation.
        """
        self.collections = collection_manager
        self.embedder = embedder
        self.llm = llm_client
        self.knowledge = knowledge

    def retrieve(self, query: AgentQuery,
                 top_k_per_collection: int = None,
                 collections_filter: List[str] = None,
                 year_min: int = None,
                 year_max: int = None,
                 conversation_context: str = None) -> CrossCollectionResult:
        """Retrieve evidence from collections for a query.

        Args:
            query: The agent query with question and optional patient profile.
            top_k_per_collection: Max results per collection (default from settings).
            collections_filter: Optional list of collection names to search.
            year_min: Optional minimum year filter for clinical evidence.
            year_max: Optional maximum year filter for clinical evidence.
            conversation_context: Optional prior conversation context for follow-ups.
        """
        top_k = top_k_per_collection or settings.TOP_K_PER_COLLECTION
        start = time.time()

        # Optionally prepend conversation context for follow-up queries
        search_text = query.question
        if conversation_context:
            # Limit context to prevent token overflow (rough estimate: 4 chars ≈ 1 token)
            max_context_chars = 2000
            if len(conversation_context) > max_context_chars:
                conversation_context = conversation_context[-max_context_chars:]
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

            # Disease area filter (applied to collections that have it)
            disease_area = self._detect_disease_area(query.question)
            if disease_area and cfg.get("has_disease_area"):
                safe_area = disease_area.strip()
                if _SAFE_FILTER_RE.match(safe_area) and len(safe_area) <= 50:
                    parts.append(f'disease_area == "{safe_area}"')

            # Year filter for clinical evidence
            year_field = cfg.get("year_field")
            if year_field:
                if year_min:
                    parts.append(f'{year_field} >= {int(year_min)}')
                if year_max:
                    parts.append(f'{year_field} <= {int(year_max)}')

            if parts:
                filter_exprs[coll] = " and ".join(parts)

        # Step 4: Parallel search across all collections
        all_hits = self._search_all_collections(
            query_embedding, collections_to_search, top_k, filter_exprs,
        )

        # Step 5: Deduplicate, score citations, rank
        hits = self._merge_and_rank(all_hits)

        # Step 6: Full knowledge graph augmentation
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

    def query(self, question: str, patient_profile=None, **kwargs) -> str:
        """Full RAG query: retrieve evidence + generate LLM response.

        Args:
            question: Natural language question about biomarkers.
            patient_profile: Optional PatientProfile for context injection.
            **kwargs: Additional query parameters.

        Returns:
            LLM-generated response string.
        """
        if not self.embedder:
            raise RuntimeError("Embedding model not initialized. Cannot run query.")
        if not self.llm:
            raise RuntimeError("LLM client not initialized. Cannot generate response.")
        agent_query = AgentQuery(question=question, patient_profile=patient_profile)
        evidence = self.retrieve(agent_query, **kwargs)
        prompt = self._build_prompt(agent_query.question, evidence, patient_profile)
        return self.llm.generate(
            prompt=prompt,
            system_prompt=BIOMARKER_SYSTEM_PROMPT,
            max_tokens=2048,
            temperature=0.7,
        )

    def query_stream(self, question: str, patient_profile=None,
                     **kwargs) -> Generator[Dict, None, None]:
        """Streaming RAG query -- yields evidence then token chunks.

        Args:
            question: Natural language question about biomarkers.
            patient_profile: Optional PatientProfile for context injection.
            **kwargs: Additional query parameters.

        Yields:
            Dicts with type='evidence', type='token', or type='done'.
        """
        if not self.embedder:
            raise RuntimeError("Embedding model not initialized. Cannot run query.")
        if not self.llm:
            raise RuntimeError("LLM client not initialized. Cannot generate response.")
        agent_query = AgentQuery(question=question, patient_profile=patient_profile)
        evidence = self.retrieve(agent_query, **kwargs)
        yield {"type": "evidence", "content": evidence}

        prompt = self._build_prompt(agent_query.question, evidence, patient_profile)
        full_answer = ""
        for token in self.llm.generate_stream(
            prompt=prompt,
            system_prompt=BIOMARKER_SYSTEM_PROMPT,
            max_tokens=2048,
            temperature=0.7,
        ):
            full_answer += token
            yield {"type": "token", "content": token}
        yield {"type": "done", "content": full_answer}

    # -- Cross-Collection Entity Linking -----------------------------------

    def find_related(self, entity: str, top_k: int = 5) -> Dict[str, List[SearchHit]]:
        """Find all evidence related to an entity across all collections.

        Enables queries like "show me everything about MTHFR" or
        "find all CYP2D6 drug interactions" spanning all 11 collections.

        Args:
            entity: Biomarker name, gene symbol, drug name, etc.
            top_k: Max results per collection.

        Returns:
            Dict of collection_name -> List[SearchHit].
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

    # -- Private Methods ---------------------------------------------------

    def _embed_query(self, text: str):
        """Embed query text with BGE instruction prefix."""
        prefix = "Represent this sentence for searching relevant passages: "
        return self.embedder.embed_text(prefix + text)

    def _detect_disease_area(self, question: str) -> Optional[str]:
        """Detect disease area from the question for collection filtering.

        Args:
            question: The user's question.

        Returns:
            Disease area string or None if no specific area detected.
        """
        q_lower = question.lower()
        disease_keywords = {
            "diabetes": ["diabetes", "hba1c", "insulin", "glucose", "metabolic syndrome", "homa"],
            "cardiovascular": ["cardiovascular", "heart", "cardiac", "lipid", "cholesterol",
                               "lp(a)", "apob", "ldl", "hdl", "atherosclerosis", "cvd"],
            "liver": ["liver", "hepatic", "nafld", "nash", "masld", "alt", "ast",
                       "ggt", "fib-4", "fibrosis", "steatosis"],
            "thyroid": ["thyroid", "tsh", "t3", "t4", "hashimoto", "hypothyroid", "hyperthyroid"],
            "iron": ["iron", "ferritin", "hemochromatosis", "hfe", "transferrin", "hepcidin"],
            "nutritional": ["vitamin", "folate", "b12", "omega", "nutrient", "mthfr", "supplement"],
        }
        for area, keywords in disease_keywords.items():
            if any(kw in q_lower for kw in keywords):
                return area
        return None

    def _search_all_collections(
        self, query_embedding, collections: List[str],
        top_k: int, filter_exprs: Dict[str, str],
    ) -> List[SearchHit]:
        """Search all collections in parallel via collection manager."""
        all_hits = []

        # Use the parallel search_all method from the collection manager
        parallel_results = self.collections.search_all(
            query_embedding,
            top_k_per_collection=top_k,
            filter_exprs=filter_exprs,
            score_threshold=settings.SCORE_THRESHOLD,
        )

        for coll_name, results in parallel_results.items():
            if coll_name not in [c for c in collections]:
                continue
            cfg = COLLECTION_CONFIG.get(coll_name, {})
            weight = cfg.get("weight", 0.1)
            label = cfg.get("label", coll_name)

            for r in results:
                raw_score = r.get("score", 0.0)
                weighted_score = min(raw_score * (1 + weight), 1.0)

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

    def _merge_and_rank(self, hits: List[SearchHit]) -> List[SearchHit]:
        """Deduplicate by ID and text content, sort by score descending, cap at 30."""
        seen_ids = set()
        seen_texts = set()
        unique = []
        for hit in hits:
            text_key = hit.text[:200].strip().lower()
            if hit.id not in seen_ids and text_key not in seen_texts:
                seen_ids.add(hit.id)
                seen_texts.add(text_key)
                unique.append(hit)
        unique.sort(key=lambda h: h.score, reverse=True)
        return unique[:MAX_MERGED_RESULTS]

    def _get_knowledge_context(self, query: str) -> str:
        """Extract knowledge graph context from ALL domains.

        Scans the query for disease domain mentions, PGx gene mentions,
        biomarker names, and PhenoAge/aging references to build context.
        """
        if not self.knowledge:
            return ""

        from .knowledge import (
            PGX_KNOWLEDGE,
            get_biomarker_context,
            get_domain_context,
            get_pgx_context,
        )

        context_parts = []
        q_lower = query.lower()
        q_upper = query.upper()

        # Check for disease domain mentions
        domain_keywords = {
            "diabetes": ["diabetes", "hba1c", "glucose", "insulin", "homa", "metabolic"],
            "cardiovascular": ["cardiovascular", "cardiac", "heart", "lipid", "cholesterol",
                               "lp(a)", "apob", "atherosclerosis"],
            "liver": ["liver", "hepatic", "nafld", "nash", "masld", "fibrosis", "fib-4"],
            "thyroid": ["thyroid", "tsh", "hashimoto", "hypothyroid"],
            "iron": ["iron", "ferritin", "hemochromatosis", "hfe", "transferrin"],
            "nutritional": ["vitamin", "folate", "b12", "omega-3", "nutrient", "supplement"],
        }
        for domain, keywords in domain_keywords.items():
            if any(kw in q_lower for kw in keywords):
                ctx = get_domain_context(domain)
                if ctx:
                    context_parts.append(ctx)
                break  # One domain context is usually sufficient

        # Check for PGx gene mentions
        for gene in PGX_KNOWLEDGE.keys():
            if gene in q_upper or gene.replace("CYP", "").lower() in q_lower:
                ctx = get_pgx_context(gene)
                if ctx:
                    context_parts.append(ctx)

        # Check for specific PGx-related terms
        pgx_keywords = ["PHARMACOGENOMIC", "PGX", "METABOLIZER", "CPIC",
                         "STAR ALLELE", "DRUG INTERACTION"]
        if any(kw in q_upper for kw in pgx_keywords):
            if not any("Pharmacogene" in p for p in context_parts):
                ctx = get_pgx_context("CYP2D6")  # Most common; provide as reference
                if ctx:
                    context_parts.append(ctx)

        # Check for biomarker-specific mentions
        biomarker_keywords = {
            "ALBUMIN": "albumin", "CREATININE": "creatinine", "CRP": "hs_CRP",
            "RDW": "rdw", "MCV": "mcv", "WBC": "wbc", "GLUCOSE": "glucose",
            "HBA1C": "HbA1c", "FERRITIN": "ferritin", "TSH": "TSH",
            "HOMOCYSTEINE": "homocysteine", "LP(A)": "Lp(a)",
            "APOB": "ApoB", "VITAMIN D": "vitamin_D_25OH",
        }
        for keyword, biomarker in biomarker_keywords.items():
            if keyword in q_upper:
                ctx = get_biomarker_context(biomarker)
                if ctx:
                    context_parts.append(ctx)
                break  # One biomarker context is enough

        # Check for aging / PhenoAge mentions
        aging_keywords = ["PHENOAGE", "GRIMAGE", "BIOLOGICAL AGE", "AGING",
                          "EPIGENETIC CLOCK", "AGE ACCELERATION"]
        if any(kw in q_upper for kw in aging_keywords):
            from .knowledge import PHENOAGE_KNOWLEDGE
            interp = PHENOAGE_KNOWLEDGE.get("clinical_interpretation", {})
            aging_lines = [
                "PhenoAge Clock Context:",
                PHENOAGE_KNOWLEDGE.get("description", ""),
                f"Reference: {PHENOAGE_KNOWLEDGE.get('reference', '')}",
            ]
            if interp.get("actionable_drivers"):
                aging_lines.append(f"Actionable Drivers: {interp['actionable_drivers']}")
            context_parts.append("\n".join(aging_lines))

        return "\n\n".join(context_parts)

    @staticmethod
    def _format_citation(collection: str, record_id: str) -> str:
        """Format a citation with clickable URL where possible."""
        if collection == "ClinicalEvidence" and record_id.isdigit():
            url = f"https://pubmed.ncbi.nlm.nih.gov/{record_id}/"
            return f"[ClinicalEvidence:PMID {record_id}]({url})"
        return f"[{collection}:{record_id}]"

    def _build_prompt(self, question: str,
                      evidence: CrossCollectionResult,
                      patient_profile=None) -> str:
        """Build LLM prompt with evidence, knowledge context, patient profile, and relevance tags.

        Args:
            question: The user's question.
            evidence: CrossCollectionResult with retrieved hits.
            patient_profile: Optional PatientProfile for context injection.
        """
        sections = []
        by_coll = evidence.hits_by_collection()

        for coll_name, hits in by_coll.items():
            section_lines = [f"### Evidence from {coll_name}"]
            for i, hit in enumerate(hits[:MAX_PROMPT_EVIDENCE], 1):
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

        patient_text = ""
        if patient_profile:
            patient_lines = [
                "\n\n### Patient Profile Context",
                f"Age: {patient_profile.age}, Sex: {patient_profile.sex}",
            ]
            if patient_profile.biomarkers:
                marker_strs = [f"{k}: {v}" for k, v in list(patient_profile.biomarkers.items())[:15]]
                patient_lines.append(f"Biomarkers: {', '.join(marker_strs)}")
            if patient_profile.genotypes:
                geno_strs = [f"{k}: {v}" for k, v in list(patient_profile.genotypes.items())[:10]]
                patient_lines.append(f"Genotypes: {', '.join(geno_strs)}")
            if patient_profile.star_alleles:
                star_strs = [f"{k}: {v}" for k, v in patient_profile.star_alleles.items()]
                patient_lines.append(f"Star Alleles: {', '.join(star_strs)}")
            patient_text = "\n".join(patient_lines)

        return (
            f"## Retrieved Evidence\n\n"
            f"{evidence_text}"
            f"{knowledge_text}"
            f"{patient_text}\n\n"
            f"---\n\n"
            f"## Question\n\n"
            f"{question}\n\n"
            f"Please provide a comprehensive answer grounded in the evidence above. "
            f"Cite sources using the labels provided in each evidence item. "
            f"Prioritize [high relevance] citations. "
            f"When patient profile data is available, provide genotype-specific interpretation. "
            f"Integrate insights across biomarker domains, genetics, and pharmacogenomics."
        )
