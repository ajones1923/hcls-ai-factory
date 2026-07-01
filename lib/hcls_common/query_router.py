"""
Cross-Collection Query Router for the HCLS AI Factory.

Routes user queries across the genomic_evidence collection (RAG pipeline) and
the 11 CAR-T Intelligence Agent collections.  Implements Section 10.3 Gap 5:
cross-collection query routing with intent classification, parallel execution,
weighted result merging, and cross-collection joins.

Pipeline context:
    Patient DNA -> VCF (Parabricks) -> Milvus RAG (genomic_evidence) -> Drug Targets
    CAR-T Agent maintains 11 specialized collections for cell-therapy intelligence.

Collections routed:
    genomic_evidence    RAG pipeline -- variants, genes, clinical significance
    cart_literature      CAR-T research papers and reviews
    cart_trials          Clinical trial protocols and results
    cart_constructs      CAR construct designs (scFv, hinge, costim, etc.)
    cart_assays          In-vitro / in-vivo efficacy assay data
    cart_manufacturing   GMP manufacturing protocols and QC
    cart_safety          Adverse events, CRS, ICANS, neurotoxicity
    cart_biomarkers      Predictive / pharmacodynamic biomarkers
    cart_regulatory      FDA / EMA regulatory filings and guidance
    cart_sequences       Antibody and scFv sequences (heavy/light chain)
    cart_realworld       Real-world evidence, registry data, post-market

Hardware target: NVIDIA DGX Spark (GB10 GPU, 128 GB unified memory, $4,699).
"""

from __future__ import annotations

import os
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from enum import Enum, unique
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from loguru import logger

try:
    from prometheus_client import Counter, Histogram

    QUERIES_ROUTED = Counter(
        "hcls_queries_routed_total",
        "Total queries routed across collections",
        ["intent"],
    )
    ROUTE_LATENCY = Histogram(
        "hcls_query_route_seconds",
        "End-to-end query routing latency",
        buckets=[0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0],
    )
    COLLECTION_SEARCH_LATENCY = Histogram(
        "hcls_collection_search_seconds",
        "Per-collection search latency",
        ["collection"],
        buckets=[0.05, 0.1, 0.25, 0.5, 1.0, 2.5],
    )
except ImportError:
    QUERIES_ROUTED = None
    ROUTE_LATENCY = None
    COLLECTION_SEARCH_LATENCY = None


# ---------------------------------------------------------------------------
# Intent taxonomy
# ---------------------------------------------------------------------------

@unique
class QueryIntent(str, Enum):
    """Semantic intent categories for cross-collection routing."""

    GENOMIC_VARIANT = "genomic_variant"
    DRUG_TARGET = "drug_target"
    CART_THERAPY = "cart_therapy"
    CLINICAL_TRIAL = "clinical_trial"
    SAFETY_TOXICITY = "safety_toxicity"
    MANUFACTURING = "manufacturing"
    BIOMARKER = "biomarker"
    REGULATORY = "regulatory"
    DRUG_DESIGN = "drug_design"
    GENERAL = "general"

    @property
    def display_name(self) -> str:
        return self.value.replace("_", " ").title()


# ---------------------------------------------------------------------------
# Keyword maps -- used by the fast keyword-based classifier
# ---------------------------------------------------------------------------

INTENT_KEYWORDS: dict[QueryIntent, list[str]] = {
    QueryIntent.GENOMIC_VARIANT: [
        "variant", "mutation", "snp", "snv", "indel", "vcf", "pathogenic",
        "likely pathogenic", "benign", "clinvar", "alphamissense", "missense",
        "nonsense", "frameshift", "splice", "synonymous", "heterozygous",
        "homozygous", "genotype", "allele", "hgvs", "rsid", "rs", "dbsnp",
        "chromosome", "chr", "exome", "genome", "germline", "somatic",
        "loss of function", "gain of function", "lof", "gof",
        "copy number", "cnv", "structural variant", "sv",
        "variant of uncertain significance", "vus",
    ],
    QueryIntent.DRUG_TARGET: [
        "drug target", "druggable", "therapeutic target", "protein target",
        "target identification", "target validation", "target hypothesis",
        "binding pocket", "active site", "allosteric", "pocket",
        "pdb", "structure", "crystal structure", "cryo-em",
        "protein function", "pathway", "mechanism of action",
        "kinase", "receptor", "enzyme", "transporter", "ion channel",
        "protease", "phosphatase", "ligase", "atpase",
        "vcp", "p97", "brca1", "brca2", "tp53", "egfr", "kras",
        "braf", "alk", "her2", "erbb2", "pik3ca",
    ],
    QueryIntent.CART_THERAPY: [
        "car-t", "cart", "car t", "chimeric antigen receptor",
        "cell therapy", "cell-therapy", "adoptive cell",
        "cd19", "cd22", "bcma", "cd20", "cd30", "cd33", "cd38",
        "cd7", "gpc3", "her2", "mesothelin", "egfrviii",
        "scfv", "single-chain", "costimulatory", "cd28", "4-1bb",
        "cd3 zeta", "tcr", "t cell receptor",
        "kymriah", "yescarta", "tecartus", "breyanzi", "abecma",
        "carvykti", "tisagenlecleucel", "axicabtagene",
        "lentivirus", "retrovirus", "transduction",
        "expansion", "t cell expansion", "lymphodepletion",
        "construct", "car construct", "vector", "payload",
    ],
    QueryIntent.CLINICAL_TRIAL: [
        "clinical trial", "trial", "phase i", "phase ii", "phase iii",
        "phase 1", "phase 2", "phase 3", "phase iv", "phase 4",
        "nct", "clinicaltrials.gov", "enrollment", "eligibility",
        "inclusion criteria", "exclusion criteria",
        "primary endpoint", "secondary endpoint", "overall survival",
        "progression-free survival", "pfs", "os", "orr",
        "complete response", "partial response", "cr", "pr",
        "dose escalation", "dose limiting toxicity", "dlt",
        "investigator", "sponsor", "multi-center", "randomized",
        "double-blind", "placebo-controlled", "open-label",
    ],
    QueryIntent.SAFETY_TOXICITY: [
        "safety", "toxicity", "adverse event", "adverse effect",
        "side effect", "cytokine release syndrome", "crs",
        "neurotoxicity", "icans", "immune effector cell-associated",
        "grade 3", "grade 4", "grade 5", "serious adverse",
        "dose limiting", "dlt", "hepatotoxicity", "cardiotoxicity",
        "nephrotoxicity", "myelosuppression", "neutropenia",
        "thrombocytopenia", "anemia", "infection",
        "gvhd", "graft versus host", "off-target",
        "on-target off-tumor", "tumor lysis syndrome",
        "b-cell aplasia", "hypogammaglobulinemia",
        "tocilizumab", "dexamethasone", "siltuximab",
        "management", "monitoring", "risk mitigation",
    ],
    QueryIntent.MANUFACTURING: [
        "manufacturing", "gmp", "good manufacturing practice",
        "cleanroom", "grade a", "grade b", "grade c",
        "leukapheresis", "apheresis", "cell processing",
        "transduction efficiency", "vector copy number", "vcn",
        "cell viability", "cell count", "potency assay",
        "release criteria", "release testing", "qc",
        "quality control", "quality assurance", "coa",
        "certificate of analysis", "batch record",
        "cryopreservation", "thaw", "formulation",
        "vein-to-vein time", "turnaround time",
        "decentralized manufacturing", "point of care",
        "scale-up", "scale-out", "bioreactor",
        "closed system", "automated manufacturing",
    ],
    QueryIntent.BIOMARKER: [
        "biomarker", "predictive biomarker", "prognostic biomarker",
        "pharmacodynamic", "pd biomarker", "companion diagnostic",
        "response biomarker", "resistance biomarker",
        "ctdna", "circulating tumor dna", "liquid biopsy",
        "mrd", "minimal residual disease",
        "cd4/cd8", "t cell phenotype", "t cell subset",
        "cytokine profile", "il-6", "il-2", "ifn-gamma",
        "tnf-alpha", "ferritin", "crp", "c-reactive protein",
        "tumor mutational burden", "tmb", "pd-l1",
        "microsatellite instability", "msi", "mmr",
        "car-t persistence", "b-cell recovery",
        "peak expansion", "auc", "cmax",
    ],
    QueryIntent.REGULATORY: [
        "regulatory", "fda", "ema", "pmda", "nmpa",
        "bla", "biologics license", "ind", "investigational new drug",
        "nda", "approval", "accelerated approval",
        "breakthrough therapy", "fast track", "priority review",
        "rems", "risk evaluation", "black box warning",
        "labeling", "prescribing information", "uspi",
        "advisory committee", "adcom", "cber", "cder",
        "post-marketing", "pharmacovigilance", "rmp",
        "conditional marketing authorization", "cma",
        "orphan drug", "rare pediatric disease",
        "regenerative medicine advanced therapy", "rmat",
    ],
    QueryIntent.DRUG_DESIGN: [
        "molecule", "compound", "drug design", "drug candidate",
        "small molecule", "smiles", "inchi", "molecular weight",
        "docking", "binding affinity", "binding energy",
        "diffdock", "molmim", "bionemo",
        "cb-5083", "nms-873", "dbeq",
        "lead optimization", "lead compound", "hit-to-lead",
        "sar", "structure-activity", "structure activity",
        "admet", "adme", "pharmacokinetics", "pk",
        "lipinski", "rule of five", "drug-likeness", "qed",
        "synthetic accessibility", "sa score", "retrosynthesis",
        "conformer", "3d structure", "pharmacophore",
        "ic50", "ki", "kd", "ec50", "potency",
        "selectivity", "off-target", "toxicophore",
        "rdkit", "scaffold", "fragment", "linker",
        "molecule generation", "de novo design", "generative",
    ],
}

# ---------------------------------------------------------------------------
# Collection configuration
# ---------------------------------------------------------------------------

# All collections managed by the platform
ALL_COLLECTIONS: list[str] = [
    "genomic_evidence",
    "cart_literature",
    "cart_trials",
    "cart_constructs",
    "cart_assays",
    "cart_manufacturing",
    "cart_safety",
    "cart_biomarkers",
    "cart_regulatory",
    "cart_sequences",
    "cart_realworld",
]

# Default weights for scoring results from each collection.
# Higher weight = results from this collection rank higher in merged output.
DEFAULT_COLLECTION_WEIGHTS: dict[str, float] = {
    "genomic_evidence": 1.0,
    "cart_literature": 0.85,
    "cart_trials": 0.90,
    "cart_constructs": 0.80,
    "cart_assays": 0.75,
    "cart_manufacturing": 0.70,
    "cart_safety": 0.95,
    "cart_biomarkers": 0.85,
    "cart_regulatory": 0.80,
    "cart_sequences": 0.70,
    "cart_realworld": 0.75,
}

# Maps each intent to the collections it should query (primary + secondary).
INTENT_COLLECTION_MAP: dict[QueryIntent, dict[str, list[str]]] = {
    QueryIntent.GENOMIC_VARIANT: {
        "primary": ["genomic_evidence"],
        "secondary": ["cart_biomarkers", "cart_literature"],
    },
    QueryIntent.DRUG_TARGET: {
        "primary": ["genomic_evidence", "cart_literature"],
        "secondary": ["cart_constructs", "cart_trials", "cart_assays"],
    },
    QueryIntent.CART_THERAPY: {
        "primary": [
            "cart_literature", "cart_constructs", "cart_trials",
        ],
        "secondary": [
            "cart_assays", "cart_manufacturing", "cart_safety",
            "cart_biomarkers", "genomic_evidence",
        ],
    },
    QueryIntent.CLINICAL_TRIAL: {
        "primary": ["cart_trials"],
        "secondary": [
            "cart_safety", "cart_biomarkers", "cart_literature",
            "cart_realworld",
        ],
    },
    QueryIntent.SAFETY_TOXICITY: {
        "primary": ["cart_safety"],
        "secondary": [
            "cart_trials", "cart_biomarkers", "cart_realworld",
            "cart_literature",
        ],
    },
    QueryIntent.MANUFACTURING: {
        "primary": ["cart_manufacturing"],
        "secondary": [
            "cart_constructs", "cart_regulatory", "cart_literature",
        ],
    },
    QueryIntent.BIOMARKER: {
        "primary": ["cart_biomarkers", "genomic_evidence"],
        "secondary": [
            "cart_trials", "cart_safety", "cart_literature",
        ],
    },
    QueryIntent.REGULATORY: {
        "primary": ["cart_regulatory"],
        "secondary": [
            "cart_trials", "cart_safety", "cart_literature",
        ],
    },
    QueryIntent.DRUG_DESIGN: {
        "primary": ["genomic_evidence"],
        "secondary": [
            "cart_literature", "cart_constructs", "cart_trials",
        ],
    },
    QueryIntent.GENERAL: {
        "primary": [
            "genomic_evidence", "cart_literature",
        ],
        "secondary": [
            "cart_trials", "cart_safety", "cart_biomarkers",
        ],
    },
}

# Cross-collection join keys.  When results from both sides of a join are
# present, the router can merge records on the shared key for richer context.
CROSS_COLLECTION_JOINS: list[dict[str, Any]] = [
    {
        "name": "genomic_variant_to_trial",
        "left_collection": "genomic_evidence",
        "right_collection": "cart_trials",
        "left_key": "gene",
        "right_key": "target_antigen",
        "description": (
            "Join patient variants with CAR-T trials targeting the same gene / antigen."
        ),
    },
    {
        "name": "genomic_variant_to_biomarker",
        "left_collection": "genomic_evidence",
        "right_collection": "cart_biomarkers",
        "left_key": "gene",
        "right_key": "gene",
        "description": (
            "Join patient variants with known CAR-T biomarkers for the same gene."
        ),
    },
    {
        "name": "trial_to_safety",
        "left_collection": "cart_trials",
        "right_collection": "cart_safety",
        "left_key": "nct_id",
        "right_key": "nct_id",
        "description": (
            "Join trial protocols with safety / adverse event data."
        ),
    },
    {
        "name": "construct_to_assay",
        "left_collection": "cart_constructs",
        "right_collection": "cart_assays",
        "left_key": "construct_id",
        "right_key": "construct_id",
        "description": (
            "Join CAR construct designs with their efficacy assay results."
        ),
    },
    {
        "name": "construct_to_sequence",
        "left_collection": "cart_constructs",
        "right_collection": "cart_sequences",
        "left_key": "construct_id",
        "right_key": "construct_id",
        "description": (
            "Join CAR construct designs with heavy/light chain sequences."
        ),
    },
]


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class IntentClassification:
    """Result of query intent classification."""

    primary_intent: QueryIntent
    confidence: float  # 0.0 -- 1.0
    secondary_intents: list[QueryIntent] = field(default_factory=list)
    matched_keywords: list[str] = field(default_factory=list)
    method: str = "keyword"  # "keyword" or "llm"
    raw_query: str = ""


@dataclass
class CollectionSearchSpec:
    """Specification for searching a single collection."""

    collection_name: str
    tier: str  # "primary" or "secondary"
    weight: float = 1.0
    top_k: int = 10
    filter_expr: str | None = None
    output_fields: list[str] | None = None


@dataclass
class QueryPlan:
    """
    Execution plan mapping an intent classification to collection searches.

    The router builds a QueryPlan, then hands it to ``execute_plan()`` which
    runs the searches in parallel via ThreadPoolExecutor.
    """

    query: str
    intent: IntentClassification
    searches: list[CollectionSearchSpec] = field(default_factory=list)
    join_specs: list[dict[str, Any]] = field(default_factory=list)
    created_at: float = field(default_factory=time.time)
    max_workers: int = 4
    timeout_seconds: float = 30.0

    @property
    def collection_names(self) -> list[str]:
        return [s.collection_name for s in self.searches]

    @property
    def primary_collections(self) -> list[str]:
        return [s.collection_name for s in self.searches if s.tier == "primary"]

    @property
    def secondary_collections(self) -> list[str]:
        return [s.collection_name for s in self.searches if s.tier == "secondary"]

    def to_dict(self) -> dict[str, Any]:
        return {
            "query": self.query,
            "intent": self.intent.primary_intent.value,
            "confidence": self.intent.confidence,
            "collections": self.collection_names,
            "join_specs": [j["name"] for j in self.join_specs],
            "created_at": self.created_at,
        }


@dataclass
class MergedResult:
    """A single result after cross-collection merge and scoring."""

    id: str
    collection: str
    score: float  # weighted similarity score
    raw_score: float  # original similarity score before weighting
    data: dict[str, Any] = field(default_factory=dict)
    joined_data: dict[str, Any] | None = None
    tier: str = "primary"


# ---------------------------------------------------------------------------
# Intent classifier
# ---------------------------------------------------------------------------

class IntentClassifier:
    """
    Two-stage intent classifier: fast keyword matching with optional LLM fallback.

    Stage 1 (keyword): Counts keyword hits per intent category.  If the top
    category has >= ``min_keyword_confidence`` score ratio over the runner-up,
    return immediately.

    Stage 2 (LLM): If keyword matching is ambiguous (ratio < threshold), send
    the query to Claude for explicit intent classification.
    """

    def __init__(
        self,
        llm_client: Any | None = None,
        min_keyword_confidence: float = 0.35,
        llm_fallback: bool = True,
    ):
        self._llm = llm_client
        self._min_keyword_confidence = min_keyword_confidence
        self._llm_fallback = llm_fallback and (llm_client is not None)

        # Pre-compile keyword patterns for faster matching
        self._compiled_patterns: dict[QueryIntent, list[re.Pattern]] = {}
        for intent, keywords in INTENT_KEYWORDS.items():
            self._compiled_patterns[intent] = [
                re.compile(r"\b" + re.escape(kw) + r"\b", re.IGNORECASE)
                for kw in keywords
            ]

    # ---- public API --------------------------------------------------------

    def classify(self, query: str) -> IntentClassification:
        """
        Classify the intent of a user query.

        Returns an IntentClassification with primary intent, confidence,
        secondary intents, and matched keywords.
        """
        if not query or not query.strip():
            return IntentClassification(
                primary_intent=QueryIntent.GENERAL,
                confidence=0.0,
                raw_query=query or "",
            )

        # Stage 1: keyword matching
        kw_result = self._classify_keywords(query)

        if kw_result.confidence >= self._min_keyword_confidence:
            return kw_result

        # Stage 2: LLM fallback for ambiguous queries
        if self._llm_fallback:
            try:
                return self._classify_llm(query, kw_result)
            except Exception as exc:
                logger.warning(f"LLM intent classification failed: {exc}")
                # Fall through to keyword result

        return kw_result

    # ---- keyword classifier ------------------------------------------------

    def _classify_keywords(self, query: str) -> IntentClassification:
        """Score each intent by keyword hit count."""
        scores: dict[QueryIntent, float] = {}
        matches: dict[QueryIntent, list[str]] = {}

        query_lower = query.lower()

        for intent, patterns in self._compiled_patterns.items():
            intent_matches: list[str] = []
            for pattern in patterns:
                if pattern.search(query_lower):
                    intent_matches.append(pattern.pattern)
            if intent_matches:
                scores[intent] = len(intent_matches)
                matches[intent] = intent_matches

        if not scores:
            return IntentClassification(
                primary_intent=QueryIntent.GENERAL,
                confidence=0.1,
                raw_query=query,
                method="keyword",
            )

        # Normalize scores to 0..1 range
        total_hits = sum(scores.values())
        normalized: dict[QueryIntent, float] = {
            intent: count / total_hits for intent, count in scores.items()
        }

        # Sort by score descending
        ranked = sorted(normalized.items(), key=lambda kv: kv[1], reverse=True)
        primary_intent = ranked[0][0]
        primary_score = ranked[0][1]

        secondary_intents = [intent for intent, _ in ranked[1:4]]

        all_matched: list[str] = []
        for kws in matches.values():
            all_matched.extend(kws)

        return IntentClassification(
            primary_intent=primary_intent,
            confidence=primary_score,
            secondary_intents=secondary_intents,
            matched_keywords=list(set(all_matched)),
            method="keyword",
            raw_query=query,
        )

    # ---- LLM classifier ----------------------------------------------------

    _LLM_INTENT_PROMPT = """You are an intent classifier for a precision medicine genomics platform.

The platform has these query intent categories:
- genomic_variant: Questions about DNA variants, mutations, SNPs, clinical significance
- drug_target: Questions about protein targets, druggability, binding sites
- cart_therapy: Questions about CAR-T cell therapy, constructs, adoptive cell therapy
- clinical_trial: Questions about clinical trials, phases, endpoints, enrollment
- safety_toxicity: Questions about adverse events, CRS, neurotoxicity, side effects
- manufacturing: Questions about GMP manufacturing, cell processing, quality control
- biomarker: Questions about biomarkers, MRD, ctDNA, response prediction
- regulatory: Questions about FDA/EMA approvals, BLA, IND, regulatory guidance
- drug_design: Questions about small molecule design, docking, SMILES, lead optimization
- general: Broad or unclear queries

Classify the following query into EXACTLY ONE primary intent.
Respond with ONLY the intent name, nothing else.

Query: {query}"""

    def _classify_llm(
        self,
        query: str,
        keyword_fallback: IntentClassification,
    ) -> IntentClassification:
        """Use Claude (or another LLM) to classify ambiguous queries."""
        prompt = self._LLM_INTENT_PROMPT.format(query=query)

        response_text: str = self._llm.generate(
            prompt=prompt,
            system_prompt="You are a precise intent classifier. Respond with only the intent category name.",
            max_tokens=32,
            temperature=0.0,
        )

        # Parse the LLM response
        response_clean = response_text.strip().lower().replace("-", "_").replace(" ", "_")

        # Map to QueryIntent
        intent_map = {intent.value: intent for intent in QueryIntent}
        if response_clean in intent_map:
            primary = intent_map[response_clean]
            return IntentClassification(
                primary_intent=primary,
                confidence=0.80,
                secondary_intents=keyword_fallback.secondary_intents,
                matched_keywords=keyword_fallback.matched_keywords,
                method="llm",
                raw_query=query,
            )

        # If LLM returned something unexpected, fall back to keyword result
        logger.warning(
            f"LLM returned unrecognized intent '{response_clean}', "
            f"falling back to keyword: {keyword_fallback.primary_intent.value}"
        )
        return keyword_fallback


# ---------------------------------------------------------------------------
# Query router
# ---------------------------------------------------------------------------

class QueryRouter:
    """
    Cross-collection query router for the HCLS AI Factory.

    Orchestrates query intent classification, plan building, parallel
    execution across Milvus collections, result merging with weighted
    scoring, and cross-collection joins.

    Usage::

        from hcls_common.query_router import QueryRouter

        router = QueryRouter(milvus_client=milvus, embedder=embedder)
        results = router.route_query("What pathogenic VCP variants exist?")
        for r in results:
            print(r.collection, r.score, r.data)

    Parameters
    ----------
    milvus_client
        A Milvus client instance that exposes ``search()`` and ``query()``
        methods.  Compatible with both the per-pipeline MilvusClient and
        the unified ``UnifiedMilvusClient`` from hcls_common.
    embedder
        An embedding provider that exposes ``embed_query(text) -> ndarray``.
    llm_client : optional
        LLM client for ambiguous intent classification (Claude fallback).
    collection_weights : optional
        Override default weights per collection.
    max_workers : int
        Maximum parallel threads for cross-collection search.
    top_k : int
        Default number of results per collection search.
    score_threshold : float
        Minimum similarity score to include a result.
    """

    def __init__(
        self,
        milvus_client: Any | None = None,
        embedder: Any | None = None,
        llm_client: Any | None = None,
        collection_weights: dict[str, float] | None = None,
        max_workers: int = 4,
        top_k: int = 10,
        score_threshold: float = 0.40,
    ):
        self._milvus = milvus_client
        self._embedder = embedder
        self._classifier = IntentClassifier(
            llm_client=llm_client,
            llm_fallback=(llm_client is not None),
        )
        self._weights = collection_weights or dict(DEFAULT_COLLECTION_WEIGHTS)
        self._max_workers = max_workers
        self._default_top_k = top_k
        self._score_threshold = score_threshold

        # Cache for query embeddings within a session
        self._embedding_cache: dict[str, Any] = {}

        logger.info(
            f"QueryRouter initialized: {len(self._weights)} collections, "
            f"max_workers={max_workers}, top_k={top_k}"
        )

    # ---- public API --------------------------------------------------------

    def classify_intent(self, query: str) -> IntentClassification:
        """Classify the intent of a user query."""
        return self._classifier.classify(query)

    def build_plan(
        self,
        query: str,
        intent: IntentClassification | None = None,
        include_secondary: bool = True,
        top_k: int | None = None,
        filter_expr: str | None = None,
    ) -> QueryPlan:
        """
        Build a QueryPlan from a query string and optional pre-classified intent.

        Parameters
        ----------
        query : str
            The user's natural language query.
        intent : IntentClassification, optional
            Pre-classified intent.  If None, ``classify_intent`` is called.
        include_secondary : bool
            Whether to include secondary collections in the plan.
        top_k : int, optional
            Override per-collection result limit.
        filter_expr : str, optional
            Milvus filter expression applied to all collections.
        """
        if intent is None:
            intent = self.classify_intent(query)

        k = top_k or self._default_top_k
        mapping = INTENT_COLLECTION_MAP.get(
            intent.primary_intent, INTENT_COLLECTION_MAP[QueryIntent.GENERAL]
        )

        searches: list[CollectionSearchSpec] = []

        # Primary collections
        for coll_name in mapping["primary"]:
            searches.append(
                CollectionSearchSpec(
                    collection_name=coll_name,
                    tier="primary",
                    weight=self._weights.get(coll_name, 1.0),
                    top_k=k,
                    filter_expr=filter_expr,
                )
            )

        # Secondary collections (reduced top_k)
        if include_secondary:
            secondary_k = max(k // 2, 3)
            for coll_name in mapping.get("secondary", []):
                searches.append(
                    CollectionSearchSpec(
                        collection_name=coll_name,
                        tier="secondary",
                        weight=self._weights.get(coll_name, 0.7) * 0.8,
                        top_k=secondary_k,
                        filter_expr=filter_expr,
                    )
                )

        # Identify applicable cross-collection joins
        plan_collections = {s.collection_name for s in searches}
        applicable_joins = [
            j for j in CROSS_COLLECTION_JOINS
            if j["left_collection"] in plan_collections
            and j["right_collection"] in plan_collections
        ]

        plan = QueryPlan(
            query=query,
            intent=intent,
            searches=searches,
            join_specs=applicable_joins,
            max_workers=self._max_workers,
        )

        logger.info(
            f"Built QueryPlan: intent={intent.primary_intent.value}, "
            f"collections={plan.collection_names}, "
            f"joins={[j['name'] for j in applicable_joins]}"
        )
        return plan

    def route_query(
        self,
        query: str,
        top_k: int | None = None,
        include_secondary: bool = True,
        filter_expr: str | None = None,
    ) -> list[MergedResult]:
        """
        End-to-end query routing: classify -> plan -> execute -> merge.

        This is the primary entry point for cross-collection queries.

        Parameters
        ----------
        query : str
            Natural language query.
        top_k : int, optional
            Results per collection.
        include_secondary : bool
            Include secondary collections.
        filter_expr : str, optional
            Milvus filter expression.

        Returns
        -------
        list[MergedResult]
            Deduplicated, weighted, ranked results across all collections.
        """
        t0 = time.time()

        plan = self.build_plan(
            query,
            include_secondary=include_secondary,
            top_k=top_k,
            filter_expr=filter_expr,
        )

        results = self.execute_plan(plan)

        elapsed = time.time() - t0

        if QUERIES_ROUTED is not None:
            QUERIES_ROUTED.labels(intent=plan.intent.primary_intent.value).inc()
        if ROUTE_LATENCY is not None:
            ROUTE_LATENCY.observe(elapsed)

        logger.info(
            f"route_query completed: {len(results)} results in {elapsed:.3f}s "
            f"(intent={plan.intent.primary_intent.value})"
        )
        return results

    def execute_plan(self, plan: QueryPlan) -> list[MergedResult]:
        """
        Execute a QueryPlan: run all collection searches in parallel,
        apply cross-collection joins, merge and deduplicate results.
        """
        if self._milvus is None or self._embedder is None:
            logger.warning(
                "QueryRouter.execute_plan called without milvus_client or embedder. "
                "Returning empty results."
            )
            return []

        # Embed the query (cached)
        query_embedding = self._get_embedding(plan.query)

        # Parallel search across all collections
        all_raw_results: dict[str, list[dict[str, Any]]] = {}

        with ThreadPoolExecutor(max_workers=plan.max_workers) as executor:
            future_to_spec: dict[Any, CollectionSearchSpec] = {}

            for spec in plan.searches:
                future = executor.submit(
                    self._search_collection,
                    spec=spec,
                    query_embedding=query_embedding,
                )
                future_to_spec[future] = spec

            for future in as_completed(future_to_spec, timeout=plan.timeout_seconds):
                spec = future_to_spec[future]
                try:
                    results = future.result()
                    all_raw_results[spec.collection_name] = results
                    logger.debug(
                        f"Collection '{spec.collection_name}': {len(results)} results"
                    )
                except Exception as exc:
                    logger.error(
                        f"Search failed for collection '{spec.collection_name}': {exc}"
                    )
                    all_raw_results[spec.collection_name] = []

        # Apply cross-collection joins
        joined_data = self._apply_joins(all_raw_results, plan.join_specs)

        # Merge and deduplicate
        merged = self._merge_results(all_raw_results, plan.searches, joined_data)

        return merged

    # ---- private helpers ---------------------------------------------------

    def _get_embedding(self, text: str) -> Any:
        """Get or compute embedding for a text string (with caching)."""
        if text in self._embedding_cache:
            return self._embedding_cache[text]

        embedding = self._embedder.embed_query(text)
        self._embedding_cache[text] = embedding

        # Keep cache bounded
        if len(self._embedding_cache) > 100:
            oldest_key = next(iter(self._embedding_cache))
            del self._embedding_cache[oldest_key]

        return embedding

    def _search_collection(
        self,
        spec: CollectionSearchSpec,
        query_embedding: Any,
    ) -> list[dict[str, Any]]:
        """
        Search a single Milvus collection.

        Adapts to whatever Milvus client interface is available:
        - UnifiedMilvusClient.search(collection_name, ...)
        - Per-pipeline MilvusClient.search(...)
        """
        t0 = time.time()

        try:
            # Try unified client interface first (collection_name as param)
            if hasattr(self._milvus, "search_collection"):
                results = self._milvus.search_collection(
                    collection_name=spec.collection_name,
                    query_embedding=query_embedding,
                    top_k=spec.top_k,
                    score_threshold=self._score_threshold,
                    filter_expr=spec.filter_expr,
                    output_fields=spec.output_fields,
                )
            elif hasattr(self._milvus, "search"):
                # Per-pipeline client -- may need collection switching
                results = self._milvus.search(
                    query_embedding=query_embedding,
                    top_k=spec.top_k,
                    score_threshold=self._score_threshold,
                    filter_expr=spec.filter_expr,
                )
            else:
                logger.warning(
                    f"Milvus client has no search interface for '{spec.collection_name}'"
                )
                results = []

            # Tag each result with its source collection
            for r in results:
                r["_collection"] = spec.collection_name
                r["_tier"] = spec.tier
                r["_weight"] = spec.weight

        except Exception as exc:
            logger.error(f"Error searching '{spec.collection_name}': {exc}")
            results = []

        elapsed = time.time() - t0
        if COLLECTION_SEARCH_LATENCY is not None:
            COLLECTION_SEARCH_LATENCY.labels(collection=spec.collection_name).observe(elapsed)

        return results

    def _apply_joins(
        self,
        raw_results: dict[str, list[dict[str, Any]]],
        join_specs: list[dict[str, Any]],
    ) -> dict[str, dict[str, Any]]:
        """
        Apply cross-collection joins.

        For each join spec, index results from the right collection by the
        join key, then attach matching records to left-side results.

        Returns a mapping of ``result_id -> joined_data`` for use during merge.
        """
        joined: dict[str, dict[str, Any]] = {}

        for join_spec in join_specs:
            left_coll = join_spec["left_collection"]
            right_coll = join_spec["right_collection"]
            left_key = join_spec["left_key"]
            right_key = join_spec["right_key"]
            join_name = join_spec["name"]

            left_results = raw_results.get(left_coll, [])
            right_results = raw_results.get(right_coll, [])

            if not left_results or not right_results:
                continue

            # Build right-side index
            right_index: dict[str, list[dict[str, Any]]] = {}
            for r in right_results:
                key_val = r.get(right_key, "")
                if key_val:
                    key_normalized = str(key_val).upper().strip()
                    if key_normalized not in right_index:
                        right_index[key_normalized] = []
                    right_index[key_normalized].append(r)

            # Match left-side results
            join_count = 0
            for lr in left_results:
                left_val = lr.get(left_key, "")
                if not left_val:
                    continue

                left_normalized = str(left_val).upper().strip()
                if left_normalized in right_index:
                    result_id = lr.get("id", "")
                    if result_id:
                        joined[result_id] = {
                            "join_name": join_name,
                            "joined_from": right_coll,
                            "join_key": left_normalized,
                            "joined_records": right_index[left_normalized],
                        }
                        join_count += 1

            if join_count > 0:
                logger.info(
                    f"Cross-collection join '{join_name}': "
                    f"{join_count} records joined ({left_coll} <-> {right_coll})"
                )

        return joined

    def _merge_results(
        self,
        raw_results: dict[str, list[dict[str, Any]]],
        search_specs: list[CollectionSearchSpec],
        joined_data: dict[str, dict[str, Any]],
    ) -> list[MergedResult]:
        """
        Merge results from all collections: deduplicate by ID, apply
        collection weights, attach join data, sort by weighted score.
        """
        # Build weight lookup from specs
        weight_by_collection: dict[str, float] = {}
        tier_by_collection: dict[str, str] = {}
        for spec in search_specs:
            weight_by_collection[spec.collection_name] = spec.weight
            tier_by_collection[spec.collection_name] = spec.tier

        seen_ids: set[str] = set()
        merged: list[MergedResult] = []

        for collection_name, results in raw_results.items():
            weight = weight_by_collection.get(collection_name, 0.7)
            tier = tier_by_collection.get(collection_name, "secondary")

            for r in results:
                result_id = r.get("id", "")

                # Deduplicate: if the same variant/record appears in multiple
                # collections, keep the one with the higher weighted score.
                dedup_key = f"{collection_name}:{result_id}" if result_id else ""

                if dedup_key and dedup_key in seen_ids:
                    continue
                if dedup_key:
                    seen_ids.add(dedup_key)

                raw_score = r.get("score", 0.0)
                weighted_score = raw_score * weight

                # Clean internal routing metadata from the user-facing data
                data = {
                    k: v for k, v in r.items()
                    if not k.startswith("_")
                }

                mr = MergedResult(
                    id=result_id or f"anon-{len(merged)}",
                    collection=collection_name,
                    score=weighted_score,
                    raw_score=raw_score,
                    data=data,
                    joined_data=joined_data.get(result_id),
                    tier=tier,
                )
                merged.append(mr)

        # Sort by weighted score descending
        merged.sort(key=lambda m: m.score, reverse=True)

        logger.debug(f"Merged {len(merged)} results from {len(raw_results)} collections")
        return merged

    # ---- convenience methods -----------------------------------------------

    def get_available_collections(self) -> list[str]:
        """Return the list of all configured collection names."""
        return list(ALL_COLLECTIONS)

    def get_collection_weights(self) -> dict[str, float]:
        """Return current collection weight configuration."""
        return dict(self._weights)

    def set_collection_weight(self, collection_name: str, weight: float) -> None:
        """Update the weight for a specific collection."""
        if weight < 0.0 or weight > 2.0:
            raise ValueError(f"Weight must be between 0.0 and 2.0, got {weight}")
        self._weights[collection_name] = weight
        logger.info(f"Updated weight for '{collection_name}': {weight}")

    def explain_plan(self, query: str) -> dict[str, Any]:
        """
        Return a human-readable explanation of how a query would be routed.

        Useful for debugging and transparency.
        """
        intent = self.classify_intent(query)
        plan = self.build_plan(query, intent=intent)

        return {
            "query": query,
            "intent": {
                "primary": intent.primary_intent.value,
                "confidence": round(intent.confidence, 3),
                "secondary": [i.value for i in intent.secondary_intents],
                "matched_keywords": intent.matched_keywords,
                "method": intent.method,
            },
            "plan": {
                "primary_collections": plan.primary_collections,
                "secondary_collections": plan.secondary_collections,
                "joins": [
                    {
                        "name": j["name"],
                        "left": j["left_collection"],
                        "right": j["right_collection"],
                        "on": f"{j['left_key']} = {j['right_key']}",
                    }
                    for j in plan.join_specs
                ],
                "total_searches": len(plan.searches),
            },
        }

    def clear_embedding_cache(self) -> None:
        """Clear the in-memory embedding cache."""
        self._embedding_cache.clear()
        logger.debug("Embedding cache cleared")


# ---------------------------------------------------------------------------
# Module-level convenience functions
# ---------------------------------------------------------------------------

def classify_intent(query: str, llm_client: Any = None) -> IntentClassification:
    """
    Standalone intent classification without a full QueryRouter instance.

    Useful for quick keyword-based classification in notebooks or scripts.
    """
    classifier = IntentClassifier(llm_client=llm_client)
    return classifier.classify(query)


def get_collections_for_intent(intent: QueryIntent) -> dict[str, list[str]]:
    """
    Return the primary and secondary collections for a given intent.
    """
    return INTENT_COLLECTION_MAP.get(intent, INTENT_COLLECTION_MAP[QueryIntent.GENERAL])


def create_query_router(
    milvus_client: Any | None = None,
    embedder: Any | None = None,
    llm_client: Any | None = None,
    **kwargs: Any,
) -> QueryRouter:
    """
    Factory function to create a configured QueryRouter.

    Parameters are forwarded to the QueryRouter constructor.
    """
    return QueryRouter(
        milvus_client=milvus_client,
        embedder=embedder,
        llm_client=llm_client,
        **kwargs,
    )
