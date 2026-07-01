"""
Meta-agent orchestration for the HCLS AI Factory precision medicine platform.

Uses Claude tool-use to coordinate across specialist sub-agents:
  - Genomic variant analysis (Milvus genomic_evidence collection)
  - CAR-T cell therapy evaluation (11 CAR-T Milvus collections)
  - Small molecule drug discovery (MolMIM + DiffDock pipeline results)

The meta-agent receives clinical questions and routes to the appropriate
data sources, then synthesizes a unified response integrating evidence
from all relevant pipeline stages.

Section 11.3.2 of the HCLS AI Factory architecture.
"""

from __future__ import annotations

import json
import logging
import os
import re
import time
import uuid
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum, unique
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

META_AGENT_SYSTEM_PROMPT = """\
You are a precision medicine clinical intelligence agent for the HCLS AI Factory.
You coordinate across genomic analysis, CAR-T cell therapy evaluation, and small
molecule drug discovery to provide unified patient treatment recommendations.

Your role:
1. Analyze clinical questions and determine which data sources are relevant.
2. Call the appropriate tools to gather evidence from each domain.
3. Synthesize findings into a coherent clinical assessment.
4. Provide confidence levels and cite your sources.

Pipeline context:
- Stage 1 (Genomics): Patient FASTQ -> VCF via NVIDIA Parabricks on DGX Spark.
  ~11.7M variants annotated with ClinVar (~2.7M entries) and AlphaMissense (71M predictions).
- Stage 2 (RAG/Chat): VCF -> Targets via Milvus vector DB (3.5M genomic evidence embeddings,
  BGE-small-en-v1.5) + LLM reasoning. Identifies drug targets from patient variants.
- Stage 3 (Drug Discovery): Targets -> Molecules via BioNeMo NIMs (MolMIM generation +
  DiffDock docking) with RDKit scoring. Produces ranked small molecule candidates.
- CAR-T Intelligence: Evaluates CAR-T cell therapy options for hematologic/oncologic targets
  using 11 specialized Milvus collections covering literature, clinical trials, constructs,
  safety profiles, manufacturing, and biomarkers.

Guidelines:
- Always ground your answers in retrieved evidence. Do not fabricate data.
- When multiple domains are relevant, query all of them before synthesizing.
- For cross-domain questions, explain how genomic findings connect to therapy options.
- Clearly distinguish between established evidence and inferred relationships.
- Use markdown formatting for readability.
- If a tool call fails, note the failure and proceed with available evidence.
- For drug target questions, always check both small molecule and CAR-T options.
"""

MAX_TOOL_ITERATIONS = 5
CLAUDE_TIMEOUT_SECONDS = 120
DEFAULT_MODEL = "claude-sonnet-4-20250514"

# CAR-T Milvus collection names
CART_COLLECTIONS = [
    "cart_literature",
    "cart_clinical_trials",
    "cart_constructs",
    "cart_safety_profiles",
    "cart_manufacturing",
    "cart_biomarkers",
    "cart_targets",
    "cart_toxicities",
    "cart_efficacy",
    "cart_resistance",
    "cart_combinations",
]


# ---------------------------------------------------------------------------
# Tool definitions for Claude tool-use
# ---------------------------------------------------------------------------

TOOL_DEFINITIONS: list[dict[str, Any]] = [
    {
        "name": "search_genomic_variants",
        "description": (
            "Search the Milvus genomic_evidence collection for patient variants. "
            "Can filter by gene symbol, genomic region (chromosome:start-end), "
            "clinical significance (pathogenic, likely_pathogenic, uncertain_significance, "
            "likely_benign, benign), or variant impact (HIGH, MODERATE, LOW, MODIFIER). "
            "Returns variant details including position, consequence, ClinVar annotation, "
            "and AlphaMissense pathogenicity scores."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "query": {
                    "type": "string",
                    "description": (
                        "Natural language query or gene name to search for in "
                        "the genomic evidence vector database."
                    ),
                },
                "gene": {
                    "type": "string",
                    "description": "Filter by gene symbol (e.g., VCP, BRCA1, EGFR).",
                },
                "chromosome": {
                    "type": "string",
                    "description": "Filter by chromosome (e.g., chr9, chrX).",
                },
                "clinical_significance": {
                    "type": "string",
                    "description": (
                        "Filter by ClinVar clinical significance: pathogenic, "
                        "likely_pathogenic, uncertain_significance, likely_benign, benign."
                    ),
                },
                "impact": {
                    "type": "string",
                    "enum": ["HIGH", "MODERATE", "LOW", "MODIFIER"],
                    "description": "Filter by variant impact level.",
                },
                "top_k": {
                    "type": "integer",
                    "description": "Maximum number of results to return (default 10).",
                    "default": 10,
                },
            },
            "required": ["query"],
        },
    },
    {
        "name": "get_variant_annotations",
        "description": (
            "Get detailed ClinVar and AlphaMissense annotations for a specific "
            "genomic variant identified by chromosome, position, reference allele, "
            "and alternate allele. Returns clinical significance, disease associations, "
            "AlphaMissense pathogenicity score and classification, review status, "
            "and HGVS notation."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "chromosome": {
                    "type": "string",
                    "description": "Chromosome (e.g., chr9).",
                },
                "position": {
                    "type": "integer",
                    "description": "Genomic position (1-based).",
                },
                "ref": {
                    "type": "string",
                    "description": "Reference allele.",
                },
                "alt": {
                    "type": "string",
                    "description": "Alternate allele.",
                },
                "variant_id": {
                    "type": "string",
                    "description": (
                        "Alternative: variant ID in chr_pos_ref_alt format. "
                        "If provided, chromosome/position/ref/alt are ignored."
                    ),
                },
            },
            "required": [],
        },
    },
    {
        "name": "search_cart_evidence",
        "description": (
            "Search the CAR-T Intelligence Agent's 11 specialized Milvus collections "
            "for CAR-T cell therapy evidence. Collections cover: literature, clinical "
            "trials, construct designs, safety profiles, manufacturing protocols, "
            "biomarkers, target antigens, toxicities, efficacy data, resistance "
            "mechanisms, and combination therapies. Best for hematologic malignancies "
            "and solid tumor immunotherapy questions."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "query": {
                    "type": "string",
                    "description": "Natural language query about CAR-T therapy.",
                },
                "collections": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": (
                        "Specific CAR-T collections to search. Options: "
                        "cart_literature, cart_clinical_trials, cart_constructs, "
                        "cart_safety_profiles, cart_manufacturing, cart_biomarkers, "
                        "cart_targets, cart_toxicities, cart_efficacy, "
                        "cart_resistance, cart_combinations. "
                        "If empty, searches all relevant collections."
                    ),
                },
                "target_antigen": {
                    "type": "string",
                    "description": (
                        "Filter by target antigen (e.g., CD19, BCMA, CD22, HER2)."
                    ),
                },
                "indication": {
                    "type": "string",
                    "description": (
                        "Filter by disease indication (e.g., ALL, DLBCL, multiple myeloma)."
                    ),
                },
                "top_k": {
                    "type": "integer",
                    "description": "Maximum number of results per collection (default 5).",
                    "default": 5,
                },
            },
            "required": ["query"],
        },
    },
    {
        "name": "get_cart_knowledge",
        "description": (
            "Retrieve structured knowledge graph context for CAR-T cell therapy. "
            "Returns information about target antigens, associated toxicity profiles, "
            "predictive biomarkers, approved products, and manufacturing considerations. "
            "Use this for factual lookups rather than semantic search."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "target_antigen": {
                    "type": "string",
                    "description": "Target antigen (e.g., CD19, BCMA, CD22).",
                },
                "category": {
                    "type": "string",
                    "enum": [
                        "targets",
                        "toxicities",
                        "biomarkers",
                        "products",
                        "manufacturing",
                    ],
                    "description": "Knowledge category to query.",
                },
                "disease": {
                    "type": "string",
                    "description": "Disease context for knowledge retrieval.",
                },
            },
            "required": [],
        },
    },
    {
        "name": "get_drug_candidates",
        "description": (
            "Retrieve drug discovery pipeline results for a specific target gene. "
            "Returns ranked small molecule candidates with SMILES, docking scores, "
            "QED (Quantitative Estimate of Drug-likeness), Lipinski properties, "
            "composite scores, and binding residue information. Reads from the "
            "pipeline output directory for completed runs."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "target_gene": {
                    "type": "string",
                    "description": "Gene symbol (e.g., VCP, EGFR, KRAS).",
                },
                "run_id": {
                    "type": "string",
                    "description": (
                        "Specific pipeline run ID. If not provided, uses the "
                        "most recent completed run for the target."
                    ),
                },
                "top_n": {
                    "type": "integer",
                    "description": "Number of top candidates to return (default 10).",
                    "default": 10,
                },
                "min_qed": {
                    "type": "number",
                    "description": "Minimum QED score filter (0-1).",
                },
                "max_lipinski_violations": {
                    "type": "integer",
                    "description": "Maximum Lipinski Rule of Five violations (default 1).",
                    "default": 1,
                },
            },
            "required": ["target_gene"],
        },
    },
    {
        "name": "get_protein_structures",
        "description": (
            "Fetch protein structure information for a drug target gene. "
            "Returns PDB IDs, experimental methods (X-ray, Cryo-EM, NMR), "
            "resolution, chain information, co-crystallized ligands, and "
            "binding site residues. Sources from the knowledge graph and "
            "pipeline structure manifest files."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "gene": {
                    "type": "string",
                    "description": "Gene symbol (e.g., VCP, EGFR).",
                },
                "pdb_id": {
                    "type": "string",
                    "description": (
                        "Specific PDB ID to look up (e.g., 5FTK). "
                        "If provided, returns details for that structure only."
                    ),
                },
                "method": {
                    "type": "string",
                    "enum": ["X-ray", "Cryo-EM", "NMR", "any"],
                    "description": "Filter by experimental method.",
                },
            },
            "required": [],
        },
    },
    {
        "name": "search_clinical_trials",
        "description": (
            "Search the cart_trials Milvus collection for relevant clinical trials. "
            "Returns trial IDs, phases, indications, interventions, enrollment status, "
            "and key outcome data. Useful for finding active trials for a target "
            "or therapeutic approach."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "query": {
                    "type": "string",
                    "description": "Natural language query about clinical trials.",
                },
                "target": {
                    "type": "string",
                    "description": "Target gene or antigen.",
                },
                "phase": {
                    "type": "string",
                    "enum": ["Phase I", "Phase II", "Phase III", "Phase IV", "any"],
                    "description": "Filter by trial phase.",
                },
                "status": {
                    "type": "string",
                    "enum": [
                        "recruiting",
                        "active",
                        "completed",
                        "terminated",
                        "any",
                    ],
                    "description": "Filter by enrollment status.",
                },
                "top_k": {
                    "type": "integer",
                    "description": "Maximum results to return (default 10).",
                    "default": 10,
                },
            },
            "required": ["query"],
        },
    },
    {
        "name": "trigger_drug_discovery",
        "description": (
            "Trigger a new drug discovery pipeline run for a target gene. "
            "Writes an event to the event bus that starts MolMIM molecule "
            "generation followed by DiffDock docking and RDKit scoring. "
            "This is an asynchronous operation -- it returns a run_id that "
            "can be used to check status later via get_drug_candidates. "
            "Use sparingly; only trigger when the user explicitly requests "
            "new molecule generation."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "target_gene": {
                    "type": "string",
                    "description": "Gene symbol to target (e.g., VCP).",
                },
                "reference_smiles": {
                    "type": "string",
                    "description": (
                        "Seed compound SMILES for MolMIM generation. "
                        "If not provided, uses default seed from knowledge base."
                    ),
                },
                "num_molecules": {
                    "type": "integer",
                    "description": "Number of molecules to generate (default 50).",
                    "default": 50,
                },
                "pdb_id": {
                    "type": "string",
                    "description": (
                        "PDB structure for docking. If not provided, uses the "
                        "best available structure from the knowledge base."
                    ),
                },
                "diversity": {
                    "type": "number",
                    "description": (
                        "MolMIM diversity parameter 0-1 (default 0.3). "
                        "Higher values produce more diverse molecules."
                    ),
                    "default": 0.3,
                },
            },
            "required": ["target_gene"],
        },
    },
]


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@unique
class ToolName(str, Enum):
    """Enumeration of available tools."""

    SEARCH_GENOMIC_VARIANTS = "search_genomic_variants"
    GET_VARIANT_ANNOTATIONS = "get_variant_annotations"
    SEARCH_CART_EVIDENCE = "search_cart_evidence"
    GET_CART_KNOWLEDGE = "get_cart_knowledge"
    GET_DRUG_CANDIDATES = "get_drug_candidates"
    GET_PROTEIN_STRUCTURES = "get_protein_structures"
    SEARCH_CLINICAL_TRIALS = "search_clinical_trials"
    TRIGGER_DRUG_DISCOVERY = "trigger_drug_discovery"


@dataclass
class ToolCallRecord:
    """Record of a single tool invocation within an agent turn."""

    tool_name: str
    tool_input: dict[str, Any]
    tool_output: Any = None
    success: bool = True
    error: str | None = None
    latency_ms: float = 0.0
    timestamp: str = field(default_factory=lambda: datetime.utcnow().isoformat())

    def to_dict(self) -> dict[str, Any]:
        return {
            "tool_name": self.tool_name,
            "tool_input": self.tool_input,
            "success": self.success,
            "error": self.error,
            "latency_ms": round(self.latency_ms, 1),
            "timestamp": self.timestamp,
        }


@dataclass
class AgentResponse:
    """Structured response from the meta-agent."""

    answer: str
    """Markdown-formatted clinical assessment."""

    sources: list[str] = field(default_factory=list)
    """List of data source identifiers used (e.g. 'genomic_evidence', 'cart_trials')."""

    confidence: float = 0.0
    """Estimated confidence 0-1, heuristic from source diversity and evidence quality."""

    follow_up_questions: list[str] = field(default_factory=list)
    """Suggested follow-up questions for the clinician."""

    tool_calls_made: list[ToolCallRecord] = field(default_factory=list)
    """Audit trail of all tool invocations."""

    latency_ms: float = 0.0
    """Total wall-clock time for the agent run in milliseconds."""

    model: str = ""
    """Claude model used for reasoning."""

    iterations: int = 0
    """Number of agentic loop iterations consumed."""

    request_id: str = field(default_factory=lambda: str(uuid.uuid4())[:12])
    """Unique identifier for this agent request."""

    def to_dict(self) -> dict[str, Any]:
        return {
            "answer": self.answer,
            "sources": self.sources,
            "confidence": round(self.confidence, 3),
            "follow_up_questions": self.follow_up_questions,
            "tool_calls_made": [tc.to_dict() for tc in self.tool_calls_made],
            "latency_ms": round(self.latency_ms, 1),
            "model": self.model,
            "iterations": self.iterations,
            "request_id": self.request_id,
        }


# ---------------------------------------------------------------------------
# Data source handlers
# ---------------------------------------------------------------------------

class GenomicDataSource:
    """
    Handler for genomic variant queries against Milvus genomic_evidence.

    Wraps the UnifiedMilvusClient and EvidenceEmbedder to provide
    semantic search and structured queries over patient variant data.
    """

    def __init__(
        self,
        milvus_host: str = "localhost",
        milvus_port: int = 19530,
        collection_name: str = "genomic_evidence",
        embedding_model: str = "BAAI/bge-small-en-v1.5",
    ):
        self.milvus_host = milvus_host
        self.milvus_port = milvus_port
        self.collection_name = collection_name
        self.embedding_model = embedding_model
        self._client = None
        self._embedder = None

    def _ensure_connected(self) -> None:
        """Lazy initialization of Milvus connection and embedder."""
        if self._client is not None:
            return

        try:
            from pymilvus import Collection, connections, utility

            connections.connect(
                alias="meta_agent_genomic",
                host=self.milvus_host,
                port=self.milvus_port,
            )
            if utility.has_collection(self.collection_name):
                self._client = Collection(self.collection_name)
                self._client.load()
                logger.info(
                    "GenomicDataSource connected to %s (%d entities)",
                    self.collection_name,
                    self._client.num_entities,
                )
            else:
                logger.warning(
                    "Collection %s does not exist in Milvus", self.collection_name
                )
        except Exception as exc:
            logger.error("Failed to connect to Milvus for genomic data: %s", exc)
            self._client = None

        try:
            from sentence_transformers import SentenceTransformer

            self._embedder = SentenceTransformer(self.embedding_model)
            logger.info("GenomicDataSource embedder loaded: %s", self.embedding_model)
        except Exception as exc:
            logger.error("Failed to load embedder: %s", exc)
            self._embedder = None

    def search_variants(
        self,
        query: str,
        gene: str | None = None,
        chromosome: str | None = None,
        clinical_significance: str | None = None,
        impact: str | None = None,
        top_k: int = 10,
    ) -> list[dict[str, Any]]:
        """
        Semantic search over genomic_evidence with optional structured filters.
        """
        self._ensure_connected()

        if self._client is None:
            return [{"error": "Milvus genomic_evidence collection not available"}]

        output_fields = [
            "chrom", "pos", "ref", "alt", "qual", "gene",
            "consequence", "impact", "genotype", "text_summary",
            "clinical_significance", "rsid", "am_pathogenicity", "am_class",
        ]

        # Build filter expression
        filters: list[str] = []
        if gene:
            safe_gene = re.sub(r'[^A-Za-z0-9_\-]', '', gene)
            filters.append(f'gene == "{safe_gene}"')
        if chromosome:
            safe_chrom = re.sub(r'[^a-zA-Z0-9]', '', chromosome)
            filters.append(f'chrom == "{safe_chrom}"')
        if clinical_significance:
            safe_cs = re.sub(r'[^a-zA-Z_]', '', clinical_significance)
            filters.append(f'clinical_significance == "{safe_cs}"')
        if impact:
            safe_impact = re.sub(r'[^A-Z]', '', impact.upper())
            filters.append(f'impact == "{safe_impact}"')

        filter_expr = " and ".join(filters) if filters else None

        # Try semantic search first
        if self._embedder is not None:
            try:
                embedding = self._embedder.encode(query, normalize_embeddings=True)
                search_params = {
                    "metric_type": "COSINE",
                    "params": {"nprobe": 16},
                }
                results = self._client.search(
                    data=[embedding.tolist()],
                    anns_field="embedding",
                    param=search_params,
                    limit=top_k,
                    expr=filter_expr,
                    output_fields=output_fields,
                )

                evidence = []
                for hits in results:
                    for hit in hits:
                        am_score = hit.entity.get("am_pathogenicity")
                        if am_score == -1.0:
                            am_score = None

                        evidence.append({
                            "id": hit.id,
                            "score": round(hit.score, 4),
                            "chrom": hit.entity.get("chrom"),
                            "pos": hit.entity.get("pos"),
                            "ref": hit.entity.get("ref"),
                            "alt": hit.entity.get("alt"),
                            "gene": hit.entity.get("gene"),
                            "consequence": hit.entity.get("consequence"),
                            "impact": hit.entity.get("impact"),
                            "clinical_significance": hit.entity.get("clinical_significance"),
                            "rsid": hit.entity.get("rsid"),
                            "text_summary": hit.entity.get("text_summary"),
                            "am_pathogenicity": am_score,
                            "am_class": hit.entity.get("am_class") or None,
                        })
                return evidence

            except Exception as exc:
                logger.warning("Semantic search failed, falling back to query: %s", exc)

        # Fallback: structured query only
        if filter_expr:
            try:
                results = self._client.query(
                    expr=filter_expr,
                    output_fields=output_fields,
                    limit=top_k,
                )
                for r in results:
                    r["score"] = 0.8  # Assign nominal score for direct query
                return results
            except Exception as exc:
                logger.error("Structured query failed: %s", exc)

        return [{"error": f"No results found for query: {query}"}]

    def get_variant_annotations(
        self,
        chromosome: str | None = None,
        position: int | None = None,
        ref: str | None = None,
        alt: str | None = None,
        variant_id: str | None = None,
    ) -> dict[str, Any]:
        """
        Look up detailed annotations for a specific variant.
        """
        self._ensure_connected()

        if self._client is None:
            return {"error": "Milvus genomic_evidence collection not available"}

        # Build lookup key
        if variant_id:
            parts = variant_id.split("_")
            if len(parts) >= 4:
                chromosome, position_str, ref, alt = parts[0], parts[1], parts[2], parts[3]
                position = int(position_str)

        if not all([chromosome, position, ref, alt]):
            return {"error": "Must provide chromosome+position+ref+alt or variant_id"}

        safe_chrom = re.sub(r'[^a-zA-Z0-9]', '', chromosome)
        safe_ref = re.sub(r'[^ACGTN]', '', ref.upper())
        safe_alt = re.sub(r'[^ACGTN]', '', alt.upper())

        output_fields = [
            "chrom", "pos", "ref", "alt", "qual", "gene",
            "consequence", "impact", "genotype", "text_summary",
            "clinical_significance", "rsid", "am_pathogenicity", "am_class",
        ]

        try:
            results = self._client.query(
                expr=f'chrom == "{safe_chrom}" and pos == {int(position)} and ref == "{safe_ref}" and alt == "{safe_alt}"',
                output_fields=output_fields,
                limit=1,
            )

            if results:
                variant = results[0]
                am_score = variant.get("am_pathogenicity")
                if am_score == -1.0:
                    am_score = None

                return {
                    "found": True,
                    "variant_id": f"{variant.get('chrom')}_{variant.get('pos')}_{variant.get('ref')}_{variant.get('alt')}",
                    "gene": variant.get("gene"),
                    "consequence": variant.get("consequence"),
                    "impact": variant.get("impact"),
                    "clinical_significance": variant.get("clinical_significance"),
                    "rsid": variant.get("rsid"),
                    "quality": variant.get("qual"),
                    "genotype": variant.get("genotype"),
                    "text_summary": variant.get("text_summary"),
                    "alphamissense": {
                        "pathogenicity_score": am_score,
                        "classification": variant.get("am_class") or None,
                    },
                }
            else:
                return {
                    "found": False,
                    "variant_id": f"{safe_chrom}_{position}_{safe_ref}_{safe_alt}",
                    "message": "Variant not found in genomic_evidence collection",
                }

        except Exception as exc:
            return {"error": f"Annotation lookup failed: {exc}"}


class CARTDataSource:
    """
    Handler for CAR-T cell therapy queries across the 11 specialized
    Milvus collections maintained by the CAR-T Intelligence Agent.
    """

    def __init__(
        self,
        milvus_host: str = "localhost",
        milvus_port: int = 19530,
        embedding_model: str = "BAAI/bge-small-en-v1.5",
    ):
        self.milvus_host = milvus_host
        self.milvus_port = milvus_port
        self.embedding_model = embedding_model
        self._connections: dict[str, Any] = {}
        self._embedder = None
        self._connected = False

    def _ensure_connected(self) -> None:
        """Lazy initialization of Milvus connections for CAR-T collections."""
        if self._connected:
            return

        try:
            from pymilvus import Collection, connections, utility

            connections.connect(
                alias="meta_agent_cart",
                host=self.milvus_host,
                port=self.milvus_port,
            )

            for coll_name in CART_COLLECTIONS:
                if utility.has_collection(coll_name):
                    coll = Collection(coll_name)
                    coll.load()
                    self._connections[coll_name] = coll
                    logger.debug("Loaded CAR-T collection: %s", coll_name)

            logger.info(
                "CARTDataSource connected: %d/%d collections available",
                len(self._connections),
                len(CART_COLLECTIONS),
            )
        except Exception as exc:
            logger.error("Failed to connect to Milvus for CAR-T data: %s", exc)

        try:
            from sentence_transformers import SentenceTransformer

            self._embedder = SentenceTransformer(self.embedding_model)
        except Exception as exc:
            logger.error("Failed to load embedder for CAR-T: %s", exc)

        self._connected = True

    def search_evidence(
        self,
        query: str,
        collections: list[str] | None = None,
        target_antigen: str | None = None,
        indication: str | None = None,
        top_k: int = 5,
    ) -> dict[str, list[dict[str, Any]]]:
        """
        Semantic search across CAR-T Milvus collections.

        Returns a dict mapping collection name -> list of results.
        """
        self._ensure_connected()

        target_collections = collections or list(self._connections.keys())

        # Validate requested collections exist
        target_collections = [
            c for c in target_collections if c in self._connections
        ]

        if not target_collections:
            return {"error": [{"message": "No CAR-T collections available"}]}

        all_results: dict[str, list[dict[str, Any]]] = {}

        # Build filter expression
        filters: list[str] = []
        if target_antigen:
            safe_antigen = re.sub(r'[^A-Za-z0-9_\-]', '', target_antigen)
            filters.append(f'target_antigen == "{safe_antigen}"')
        if indication:
            safe_indication = re.sub(r"[^A-Za-z0-9_ \-']", '', indication)
            filters.append(f'indication like "%{safe_indication}%"')

        filter_expr = " and ".join(filters) if filters else None

        if self._embedder is None:
            return {"error": [{"message": "Embedder not available for CAR-T search"}]}

        embedding = self._embedder.encode(query, normalize_embeddings=True)

        for coll_name in target_collections:
            coll = self._connections[coll_name]
            try:
                # Determine output fields based on collection schema
                schema_fields = [f.name for f in coll.schema.fields]
                output_fields = [
                    f for f in schema_fields
                    if f not in ("id", "embedding", "pk") and f != "vector"
                ]

                search_params = {
                    "metric_type": "COSINE",
                    "params": {"nprobe": 16},
                }

                # Find the vector field name
                vector_field = "embedding"
                for f_schema in coll.schema.fields:
                    if f_schema.dtype.name in ("FLOAT_VECTOR", "FLOAT16_VECTOR"):
                        vector_field = f_schema.name
                        break

                results = coll.search(
                    data=[embedding.tolist()],
                    anns_field=vector_field,
                    param=search_params,
                    limit=top_k,
                    expr=filter_expr,
                    output_fields=output_fields[:16],  # Milvus output field limit
                )

                coll_results = []
                for hits in results:
                    for hit in hits:
                        entry = {
                            "id": hit.id,
                            "score": round(hit.score, 4),
                            "collection": coll_name,
                        }
                        for field_name in output_fields[:16]:
                            val = hit.entity.get(field_name)
                            if val is not None:
                                entry[field_name] = val
                        coll_results.append(entry)

                if coll_results:
                    all_results[coll_name] = coll_results

            except Exception as exc:
                logger.warning("Search failed for collection %s: %s", coll_name, exc)
                all_results[coll_name] = [{"error": str(exc)}]

        return all_results

    def get_knowledge(
        self,
        target_antigen: str | None = None,
        category: str | None = None,
        disease: str | None = None,
    ) -> dict[str, Any]:
        """
        Structured knowledge graph lookup for CAR-T context.

        Returns factual information about CAR-T targets, toxicities,
        biomarkers, approved products, and manufacturing.
        """
        # In-memory knowledge base for common CAR-T targets
        cart_knowledge = {
            "CD19": {
                "targets": {
                    "antigen": "CD19",
                    "expression": "B-cell lineage (pan-B marker)",
                    "indications": [
                        "B-cell ALL",
                        "DLBCL",
                        "Follicular lymphoma",
                        "Mantle cell lymphoma",
                    ],
                    "approved_products": [
                        "Tisagenlecleucel (Kymriah)",
                        "Axicabtagene ciloleucel (Yescarta)",
                        "Lisocabtagene maraleucel (Breyanzi)",
                        "Brexucabtagene autoleucel (Tecartus)",
                    ],
                },
                "toxicities": {
                    "crs": "Cytokine Release Syndrome (Grade 3-4: 20-45%)",
                    "icans": "Immune effector cell-associated neurotoxicity (ICANS: 10-30%)",
                    "b_cell_aplasia": "Expected on-target off-tumor effect, managed with IgG replacement",
                    "management": "Tocilizumab (anti-IL6R) for CRS, corticosteroids for ICANS",
                },
                "biomarkers": {
                    "predictive": ["Tumor burden", "LDH levels", "Lymphodepletion intensity"],
                    "monitoring": ["CAR-T expansion peaks", "B-cell aplasia duration", "Serum IL-6"],
                },
                "manufacturing": {
                    "process": "Leukapheresis -> T-cell isolation -> Activation -> Transduction -> Expansion -> QC -> Infusion",
                    "vein_to_vein_time": "3-5 weeks",
                    "challenges": ["T-cell fitness in heavily pretreated patients", "Manufacturing failures ~5-10%"],
                },
            },
            "BCMA": {
                "targets": {
                    "antigen": "BCMA (B-cell maturation antigen / TNFRSF17)",
                    "expression": "Plasma cells, late-stage B cells",
                    "indications": ["Multiple myeloma", "Relapsed/refractory MM"],
                    "approved_products": [
                        "Idecabtagene vicleucel (Abecma)",
                        "Ciltacabtagene autoleucel (Carvykti)",
                    ],
                },
                "toxicities": {
                    "crs": "Cytokine Release Syndrome (Grade 3-4: 5-15%)",
                    "icans": "Neurotoxicity including movement/cognitive disorders (unique to BCMA CARs)",
                    "cytopenias": "Prolonged cytopenias common (40-60%)",
                    "management": "Tocilizumab for CRS, supportive care for cytopenias",
                },
                "biomarkers": {
                    "predictive": ["Soluble BCMA levels", "Tumor burden", "Extramedullary disease"],
                    "monitoring": ["CAR-T persistence", "sBCMA kinetics", "M-protein"],
                },
                "manufacturing": {
                    "process": "Similar to CD19 CARs, often with 4-1BB costimulatory domain",
                    "vein_to_vein_time": "4-6 weeks",
                    "challenges": ["Disease progression during manufacturing", "Prior BCMA-targeted therapy impact"],
                },
            },
            "CD22": {
                "targets": {
                    "antigen": "CD22",
                    "expression": "B-cell lineage (more restricted than CD19)",
                    "indications": ["B-cell ALL (CD19-relapsed)", "B-cell lymphomas"],
                    "approved_products": [],
                },
                "toxicities": {
                    "crs": "CRS generally milder than CD19 CARs",
                    "icans": "Lower ICANS rates compared to CD19 CARs",
                },
                "biomarkers": {
                    "predictive": ["CD22 surface density", "Prior CD19 CAR-T exposure"],
                    "monitoring": ["CD22 expression dynamics", "Antigen escape"],
                },
                "manufacturing": {
                    "process": "Standard CAR-T manufacturing",
                    "vein_to_vein_time": "3-5 weeks",
                    "challenges": ["CD22 downregulation/escape", "Dual-targeting strategies being explored"],
                },
            },
        }

        result: dict[str, Any] = {}

        if target_antigen:
            antigen_upper = target_antigen.upper()
            if antigen_upper in cart_knowledge:
                knowledge = cart_knowledge[antigen_upper]
                if category and category in knowledge:
                    result = {
                        "antigen": antigen_upper,
                        "category": category,
                        "data": knowledge[category],
                    }
                else:
                    result = {
                        "antigen": antigen_upper,
                        "data": knowledge,
                    }
            else:
                result = {
                    "antigen": antigen_upper,
                    "available": False,
                    "message": f"No structured knowledge available for {antigen_upper}. "
                    f"Available antigens: {list(cart_knowledge.keys())}",
                }
        elif disease:
            # Find relevant antigens for the disease
            matching = {}
            for antigen, data in cart_knowledge.items():
                indications = data.get("targets", {}).get("indications", [])
                if any(disease.lower() in ind.lower() for ind in indications):
                    matching[antigen] = data
            result = {
                "disease": disease,
                "matching_antigens": list(matching.keys()),
                "data": matching,
            }
        else:
            result = {
                "available_antigens": list(cart_knowledge.keys()),
                "total_knowledge_entries": sum(
                    len(v) for v in cart_knowledge.values()
                ),
            }

        return result


class DrugDiscoveryDataSource:
    """
    Handler for drug discovery pipeline results.

    Reads ranked candidate JSON from pipeline output directories and
    structure manifest files for protein structure information.
    """

    def __init__(
        self,
        output_base_dir: str | Path = "outputs",
        knowledge_module: Any = None,
    ):
        self.output_base_dir = Path(output_base_dir)
        self._knowledge = knowledge_module

    def _load_knowledge(self) -> dict[str, Any]:
        """Lazy-load the knowledge connections database."""
        if self._knowledge is not None:
            return self._knowledge

        try:
            # Try importing from the RAG pipeline knowledge module
            from importlib import import_module

            mod = import_module("rag_chat_pipeline.src.knowledge")
            self._knowledge = getattr(mod, "KNOWLEDGE_CONNECTIONS", {})
        except ImportError:
            try:
                # Fallback: try direct import from well-known path
                import sys

                rag_src = Path(__file__).resolve().parent.parent.parent / "core/engines/precision-intelligence" / "src"
                if rag_src.exists():
                    sys.path.insert(0, str(rag_src.parent))
                    from src.knowledge import KNOWLEDGE_CONNECTIONS  # type: ignore

                    self._knowledge = KNOWLEDGE_CONNECTIONS
                    sys.path.pop(0)
                else:
                    self._knowledge = {}
            except Exception:
                self._knowledge = {}

        return self._knowledge

    def get_drug_candidates(
        self,
        target_gene: str,
        run_id: str | None = None,
        top_n: int = 10,
        min_qed: float | None = None,
        max_lipinski_violations: int = 1,
    ) -> dict[str, Any]:
        """
        Retrieve ranked drug candidates from pipeline output.
        """
        safe_gene = re.sub(r'[^A-Za-z0-9_\-]', '', target_gene).upper()

        # Find pipeline output files
        candidate_files: list[Path] = []

        # Search patterns for pipeline output
        search_patterns = [
            self.output_base_dir / f"*{safe_gene}*" / "ranked_candidates.json",
            self.output_base_dir / f"*{safe_gene}*" / "results.json",
            self.output_base_dir / "ranked_candidates.json",
            self.output_base_dir / "results.json",
        ]

        for pattern in search_patterns:
            candidate_files.extend(
                sorted(pattern.parent.glob(pattern.name), key=lambda p: p.stat().st_mtime, reverse=True)
                if pattern.parent.exists()
                else []
            )

        # Also search recursively
        if self.output_base_dir.exists():
            for json_file in sorted(
                self.output_base_dir.rglob("ranked_candidates.json"),
                key=lambda p: p.stat().st_mtime,
                reverse=True,
            ):
                if json_file not in candidate_files:
                    candidate_files.append(json_file)

        if run_id:
            candidate_files = [
                f for f in candidate_files if run_id in str(f)
            ]

        if not candidate_files:
            # Return knowledge-base info if available
            knowledge = self._load_knowledge()
            gene_info = knowledge.get(safe_gene, {})

            if gene_info:
                return {
                    "target_gene": safe_gene,
                    "pipeline_results_available": False,
                    "message": (
                        f"No completed pipeline runs found for {safe_gene}. "
                        f"Knowledge base indicates: {gene_info.get('drugs', 'no drug info')}. "
                        f"Drug status: {gene_info.get('drug_status', 'unknown')}."
                    ),
                    "known_drugs": gene_info.get("drugs", []),
                    "drug_status": gene_info.get("drug_status"),
                    "pdb_ids": gene_info.get("pdb_ids", []),
                    "druggable": gene_info.get("druggable", False),
                }

            return {
                "target_gene": safe_gene,
                "pipeline_results_available": False,
                "message": f"No pipeline results or knowledge base entries found for {safe_gene}.",
            }

        # Load the most recent results file
        results_file = candidate_files[0]
        try:
            with open(results_file) as f:
                data = json.load(f)
        except Exception as exc:
            return {
                "target_gene": safe_gene,
                "error": f"Failed to load results from {results_file}: {exc}",
            }

        # Parse candidates -- handle both list and dict formats
        candidates = data if isinstance(data, list) else data.get("candidates", data.get("ranked_candidates", []))

        # Apply filters
        filtered = []
        for cand in candidates:
            if min_qed and cand.get("qed_score", 0) < min_qed:
                continue
            if cand.get("lipinski_violations", 0) > max_lipinski_violations:
                continue
            filtered.append(cand)

        # Take top N
        filtered = filtered[:top_n]

        return {
            "target_gene": safe_gene,
            "pipeline_results_available": True,
            "results_file": str(results_file),
            "total_candidates": len(candidates),
            "filtered_candidates": len(filtered),
            "top_candidates": filtered,
            "run_metadata": {
                "run_id": data.get("run_id", run_id or "unknown"),
                "completed_at": data.get("completed_at"),
                "pipeline_config": data.get("config", {}),
            },
        }

    def get_protein_structures(
        self,
        gene: str | None = None,
        pdb_id: str | None = None,
        method: str | None = None,
    ) -> dict[str, Any]:
        """
        Fetch protein structure information from knowledge base and manifests.
        """
        knowledge = self._load_knowledge()
        results: dict[str, Any] = {"structures": []}

        if pdb_id:
            safe_pdb = re.sub(r'[^A-Za-z0-9]', '', pdb_id).upper()
            # Search knowledge base for matching PDB
            for gene_name, info in knowledge.items():
                if safe_pdb in [p.upper() for p in info.get("pdb_ids", [])]:
                    results["gene"] = gene_name
                    results["protein"] = info.get("protein")
                    results["structures"].append({
                        "pdb_id": safe_pdb,
                        "gene": gene_name,
                        "protein": info.get("protein"),
                        "druggable": info.get("druggable", False),
                    })
                    break

        elif gene:
            safe_gene = re.sub(r'[^A-Za-z0-9_\-]', '', gene).upper()
            gene_info = knowledge.get(safe_gene, {})

            if gene_info:
                results["gene"] = safe_gene
                results["protein"] = gene_info.get("protein")
                results["function"] = gene_info.get("function")
                results["pathway"] = gene_info.get("pathway")
                results["druggable"] = gene_info.get("druggable", False)
                results["diseases"] = gene_info.get("diseases", [])

                for pid in gene_info.get("pdb_ids", []):
                    structure = {
                        "pdb_id": pid,
                        "gene": safe_gene,
                    }
                    results["structures"].append(structure)

            # Also check for structure manifests in pipeline output
            if self.output_base_dir.exists():
                for manifest_file in self.output_base_dir.rglob("structure_manifest.json"):
                    try:
                        with open(manifest_file) as f:
                            manifest = json.load(f)
                        if manifest.get("target_gene", "").upper() == safe_gene:
                            for struct in manifest.get("structures", []):
                                results["structures"].append(struct)
                            results["primary_structure"] = manifest.get("primary_structure")
                    except Exception:
                        pass

        # Filter by method if specified
        if method and method != "any":
            results["structures"] = [
                s for s in results["structures"]
                if s.get("method", "").lower() == method.lower()
            ]

        results["count"] = len(results["structures"])
        return results

    def search_clinical_trials(
        self,
        query: str,
        target: str | None = None,
        phase: str | None = None,
        status: str | None = None,
        top_k: int = 10,
    ) -> list[dict[str, Any]]:
        """
        Search for clinical trials. Currently queries the cart_clinical_trials
        collection if available, or returns knowledge-base trial references.
        """
        knowledge = self._load_knowledge()
        trial_results: list[dict[str, Any]] = []

        # Extract trial references from knowledge base
        query_lower = query.lower()
        for gene_name, info in knowledge.items():
            drugs = info.get("drugs", [])
            drug_status = info.get("drug_status", "")

            gene_relevant = (
                (target and gene_name.upper() == target.upper())
                or gene_name.lower() in query_lower
                or any(d.lower() in query_lower for d in (drugs if isinstance(drugs, list) else [drugs]))
            )

            if gene_relevant:
                trial_entry = {
                    "gene": gene_name,
                    "protein": info.get("protein"),
                    "drugs": drugs,
                    "drug_status": drug_status,
                    "diseases": info.get("diseases", []),
                    "druggable": info.get("druggable", False),
                    "source": "knowledge_base",
                }
                trial_results.append(trial_entry)

        # Filter by phase if specified
        if phase and phase != "any":
            trial_results = [
                t for t in trial_results
                if phase.lower() in str(t.get("drug_status", "")).lower()
            ]

        return trial_results[:top_k]

    def trigger_pipeline(
        self,
        target_gene: str,
        reference_smiles: str | None = None,
        num_molecules: int = 50,
        pdb_id: str | None = None,
        diversity: float = 0.3,
    ) -> dict[str, Any]:
        """
        Write an event to the event bus to trigger a drug discovery run.

        Returns the assigned run_id for tracking.
        """
        safe_gene = re.sub(r'[^A-Za-z0-9_\-]', '', target_gene).upper()
        run_id = f"{safe_gene}-{datetime.utcnow().strftime('%Y%m%d%H%M%S')}-{uuid.uuid4().hex[:6]}"

        # Look up defaults from knowledge base
        knowledge = self._load_knowledge()
        gene_info = knowledge.get(safe_gene, {})

        if not reference_smiles and gene_info.get("drugs"):
            # Note: actual SMILES would come from a compound database
            reference_smiles = None  # Leave for pipeline to resolve

        if not pdb_id and gene_info.get("pdb_ids"):
            pdb_id = gene_info["pdb_ids"][0]

        event = {
            "event_type": "drug_discovery_trigger",
            "run_id": run_id,
            "target_gene": safe_gene,
            "reference_smiles": reference_smiles,
            "num_molecules": num_molecules,
            "pdb_id": pdb_id,
            "diversity": diversity,
            "requested_at": datetime.utcnow().isoformat(),
            "status": "queued",
        }

        # Write event to event bus directory
        event_dir = self.output_base_dir / "events"
        event_dir.mkdir(parents=True, exist_ok=True)
        event_file = event_dir / f"trigger_{run_id}.json"

        try:
            with open(event_file, "w") as f:
                json.dump(event, f, indent=2)
            logger.info("Drug discovery trigger event written: %s", event_file)
        except Exception as exc:
            logger.error("Failed to write trigger event: %s", exc)
            event["error"] = str(exc)

        return {
            "run_id": run_id,
            "target_gene": safe_gene,
            "status": "queued",
            "pdb_id": pdb_id,
            "num_molecules": num_molecules,
            "message": (
                f"Drug discovery pipeline triggered for {safe_gene}. "
                f"Run ID: {run_id}. Use get_drug_candidates with this run_id "
                f"to check results once the pipeline completes."
            ),
        }


# ---------------------------------------------------------------------------
# MetaAgent
# ---------------------------------------------------------------------------

class MetaAgent:
    """
    Precision medicine meta-agent using Claude tool-use for multi-domain
    clinical intelligence.

    Receives clinical questions and orchestrates across genomic analysis,
    CAR-T evaluation, and drug discovery to provide unified treatment
    recommendations.

    Usage::

        agent = MetaAgent(anthropic_api_key="sk-ant-...")
        response = agent.ask("What therapeutic options exist for VCP mutations in FTD?")
        print(response.answer)
        print(f"Confidence: {response.confidence}")
        print(f"Sources: {response.sources}")
    """

    def __init__(
        self,
        anthropic_api_key: str | None = None,
        model: str = DEFAULT_MODEL,
        system_prompt: str = META_AGENT_SYSTEM_PROMPT,
        max_iterations: int = MAX_TOOL_ITERATIONS,
        timeout: int = CLAUDE_TIMEOUT_SECONDS,
        milvus_host: str = "localhost",
        milvus_port: int = 19530,
        embedding_model: str = "BAAI/bge-small-en-v1.5",
        output_base_dir: str | Path = "outputs",
        knowledge_module: Any = None,
    ):
        """
        Initialize the meta-agent.

        Args:
            anthropic_api_key: Anthropic API key. Falls back to ANTHROPIC_API_KEY env var.
            model: Claude model identifier.
            system_prompt: System prompt for Claude.
            max_iterations: Maximum agentic loop iterations before forced synthesis.
            timeout: Claude API call timeout in seconds.
            milvus_host: Milvus vector database host.
            milvus_port: Milvus vector database port.
            embedding_model: Sentence transformer model for embedding queries.
            output_base_dir: Base directory for drug discovery pipeline outputs.
            knowledge_module: Optional pre-loaded knowledge connections dict.
        """
        self.api_key = anthropic_api_key or os.getenv("ANTHROPIC_API_KEY")
        if not self.api_key:
            raise ValueError(
                "Anthropic API key required. Set ANTHROPIC_API_KEY or pass anthropic_api_key."
            )

        self.model = model
        self.system_prompt = system_prompt
        self.max_iterations = max_iterations
        self.timeout = timeout

        # Initialize Anthropic client
        try:
            import anthropic

            self._client = anthropic.Anthropic(
                api_key=self.api_key,
                timeout=self.timeout,
            )
        except ImportError as exc:
            raise ImportError(
                "anthropic package required. Install with: pip install anthropic"
            ) from exc

        # Initialize data source handlers
        self._genomic_source = GenomicDataSource(
            milvus_host=milvus_host,
            milvus_port=milvus_port,
            embedding_model=embedding_model,
        )
        self._cart_source = CARTDataSource(
            milvus_host=milvus_host,
            milvus_port=milvus_port,
            embedding_model=embedding_model,
        )
        self._drug_source = DrugDiscoveryDataSource(
            output_base_dir=output_base_dir,
            knowledge_module=knowledge_module,
        )

        logger.info(
            "MetaAgent initialized: model=%s, max_iterations=%d",
            self.model,
            self.max_iterations,
        )

    def ask(
        self,
        question: str,
        patient_id: str | None = None,
        context: str | None = None,
        conversation_history: list[dict[str, str]] | None = None,
    ) -> AgentResponse:
        """
        Main entry point: ask a clinical question and get a synthesized answer.

        Runs an agentic loop where Claude can call tools to gather evidence,
        up to max_iterations. Each iteration Claude decides whether to call
        another tool or provide a final answer.

        Args:
            question: Clinical question in natural language.
            patient_id: Optional patient identifier for report context.
            context: Optional additional context to prepend to the question.
            conversation_history: Optional prior conversation turns for continuity.

        Returns:
            AgentResponse with the synthesized answer and metadata.
        """
        start_time = time.monotonic()
        request_id = str(uuid.uuid4())[:12]

        logger.info(
            "MetaAgent.ask request_id=%s question=%s",
            request_id,
            question[:100],
        )

        # Build the initial message
        user_content = question
        if context:
            user_content = f"Context: {context}\n\nQuestion: {question}"
        if patient_id:
            user_content = f"Patient ID: {patient_id}\n\n{user_content}"

        # Build messages list
        messages: list[dict[str, Any]] = []
        if conversation_history:
            messages.extend(conversation_history)
        messages.append({"role": "user", "content": user_content})

        # Track tool calls
        tool_calls_made: list[ToolCallRecord] = []
        sources: set[str] = set()
        iterations = 0

        # Agentic loop
        for iteration in range(self.max_iterations):
            iterations = iteration + 1
            logger.debug(
                "MetaAgent iteration %d/%d (request_id=%s)",
                iterations,
                self.max_iterations,
                request_id,
            )

            # Call Claude with tool definitions
            response = self._call_claude(messages)

            if response is None:
                # API error -- synthesize from what we have
                logger.error("Claude API call failed on iteration %d", iterations)
                break

            # Check stop reason
            stop_reason = response.stop_reason

            # Process response content blocks
            assistant_content = []
            text_blocks: list[str] = []
            tool_use_blocks: list[dict[str, Any]] = []

            for block in response.content:
                if block.type == "text":
                    text_blocks.append(block.text)
                    assistant_content.append({"type": "text", "text": block.text})
                elif block.type == "tool_use":
                    tool_use_blocks.append({
                        "id": block.id,
                        "name": block.name,
                        "input": block.input,
                    })
                    assistant_content.append({
                        "type": "tool_use",
                        "id": block.id,
                        "name": block.name,
                        "input": block.input,
                    })

            # Add assistant response to messages
            messages.append({"role": "assistant", "content": assistant_content})

            # If Claude is done (no tool calls), we have our answer
            if stop_reason == "end_turn" or not tool_use_blocks:
                final_answer = "\n".join(text_blocks)
                elapsed_ms = (time.monotonic() - start_time) * 1000

                return AgentResponse(
                    answer=final_answer,
                    sources=sorted(sources),
                    confidence=self._estimate_confidence(tool_calls_made, sources),
                    follow_up_questions=self._extract_follow_ups(final_answer),
                    tool_calls_made=tool_calls_made,
                    latency_ms=elapsed_ms,
                    model=self.model,
                    iterations=iterations,
                    request_id=request_id,
                )

            # Execute tool calls and add results
            tool_results: list[dict[str, Any]] = []
            for tool_call in tool_use_blocks:
                tool_start = time.monotonic()
                tool_name = tool_call["name"]
                tool_input = tool_call["input"]
                tool_id = tool_call["id"]

                logger.info(
                    "Executing tool: %s (request_id=%s, iteration=%d)",
                    tool_name,
                    request_id,
                    iterations,
                )

                try:
                    result = self._execute_tool(tool_name, tool_input)
                    tool_elapsed = (time.monotonic() - tool_start) * 1000

                    # Track the source
                    source_name = self._tool_to_source(tool_name)
                    sources.add(source_name)

                    record = ToolCallRecord(
                        tool_name=tool_name,
                        tool_input=tool_input,
                        tool_output=_truncate_output(result),
                        success=True,
                        latency_ms=tool_elapsed,
                    )
                    tool_calls_made.append(record)

                    # Format result for Claude
                    result_str = json.dumps(result, indent=2, default=str)
                    if len(result_str) > 15000:
                        result_str = result_str[:15000] + "\n... (truncated)"

                    tool_results.append({
                        "type": "tool_result",
                        "tool_use_id": tool_id,
                        "content": result_str,
                    })

                except Exception as exc:
                    tool_elapsed = (time.monotonic() - tool_start) * 1000
                    logger.error("Tool execution failed: %s - %s", tool_name, exc)

                    record = ToolCallRecord(
                        tool_name=tool_name,
                        tool_input=tool_input,
                        success=False,
                        error=str(exc),
                        latency_ms=tool_elapsed,
                    )
                    tool_calls_made.append(record)

                    tool_results.append({
                        "type": "tool_result",
                        "tool_use_id": tool_id,
                        "content": json.dumps({
                            "error": str(exc),
                            "tool": tool_name,
                        }),
                        "is_error": True,
                    })

            # Add tool results to messages
            messages.append({"role": "user", "content": tool_results})

        # Max iterations exhausted -- forced synthesis
        logger.warning(
            "Max iterations (%d) exhausted for request_id=%s, forcing synthesis",
            self.max_iterations,
            request_id,
        )

        final_answer = self._synthesize(messages, tool_calls_made, question)
        elapsed_ms = (time.monotonic() - start_time) * 1000

        return AgentResponse(
            answer=final_answer,
            sources=sorted(sources),
            confidence=max(0.1, self._estimate_confidence(tool_calls_made, sources) - 0.1),
            follow_up_questions=self._extract_follow_ups(final_answer),
            tool_calls_made=tool_calls_made,
            latency_ms=elapsed_ms,
            model=self.model,
            iterations=iterations,
            request_id=request_id,
        )

    def _call_claude(
        self,
        messages: list[dict[str, Any]],
    ) -> Any:
        """
        Call Claude API with tool definitions.

        Args:
            messages: Conversation messages.

        Returns:
            Claude API response object, or None if the call fails.
        """
        try:
            response = self._client.messages.create(
                model=self.model,
                max_tokens=4096,
                temperature=0.3,
                system=self.system_prompt,
                tools=TOOL_DEFINITIONS,
                messages=messages,
            )
            return response

        except Exception as exc:
            logger.error("Claude API call failed: %s", exc)
            return None

    def _execute_tool(
        self,
        tool_name: str,
        tool_input: dict[str, Any],
    ) -> Any:
        """
        Route a tool call to the appropriate data source handler.

        Args:
            tool_name: Name of the tool to execute.
            tool_input: Input parameters for the tool.

        Returns:
            Tool execution result (JSON-serializable).

        Raises:
            ValueError: If tool_name is not recognized.
        """
        if tool_name == "search_genomic_variants":
            return self._genomic_source.search_variants(
                query=tool_input.get("query", ""),
                gene=tool_input.get("gene"),
                chromosome=tool_input.get("chromosome"),
                clinical_significance=tool_input.get("clinical_significance"),
                impact=tool_input.get("impact"),
                top_k=tool_input.get("top_k", 10),
            )

        elif tool_name == "get_variant_annotations":
            return self._genomic_source.get_variant_annotations(
                chromosome=tool_input.get("chromosome"),
                position=tool_input.get("position"),
                ref=tool_input.get("ref"),
                alt=tool_input.get("alt"),
                variant_id=tool_input.get("variant_id"),
            )

        elif tool_name == "search_cart_evidence":
            return self._cart_source.search_evidence(
                query=tool_input.get("query", ""),
                collections=tool_input.get("collections"),
                target_antigen=tool_input.get("target_antigen"),
                indication=tool_input.get("indication"),
                top_k=tool_input.get("top_k", 5),
            )

        elif tool_name == "get_cart_knowledge":
            return self._cart_source.get_knowledge(
                target_antigen=tool_input.get("target_antigen"),
                category=tool_input.get("category"),
                disease=tool_input.get("disease"),
            )

        elif tool_name == "get_drug_candidates":
            return self._drug_source.get_drug_candidates(
                target_gene=tool_input.get("target_gene", ""),
                run_id=tool_input.get("run_id"),
                top_n=tool_input.get("top_n", 10),
                min_qed=tool_input.get("min_qed"),
                max_lipinski_violations=tool_input.get("max_lipinski_violations", 1),
            )

        elif tool_name == "get_protein_structures":
            return self._drug_source.get_protein_structures(
                gene=tool_input.get("gene"),
                pdb_id=tool_input.get("pdb_id"),
                method=tool_input.get("method"),
            )

        elif tool_name == "search_clinical_trials":
            return self._drug_source.search_clinical_trials(
                query=tool_input.get("query", ""),
                target=tool_input.get("target"),
                phase=tool_input.get("phase"),
                status=tool_input.get("status"),
                top_k=tool_input.get("top_k", 10),
            )

        elif tool_name == "trigger_drug_discovery":
            return self._drug_source.trigger_pipeline(
                target_gene=tool_input.get("target_gene", ""),
                reference_smiles=tool_input.get("reference_smiles"),
                num_molecules=tool_input.get("num_molecules", 50),
                pdb_id=tool_input.get("pdb_id"),
                diversity=tool_input.get("diversity", 0.3),
            )

        else:
            raise ValueError(f"Unknown tool: {tool_name}")

    def _synthesize(
        self,
        messages: list[dict[str, Any]],
        tool_calls: list[ToolCallRecord],
        original_question: str,
    ) -> str:
        """
        Emergency fallback synthesis when max iterations are exhausted.

        Makes one final Claude call without tools to synthesize available evidence.
        """
        # Collect all successful tool outputs
        evidence_summary_parts: list[str] = []
        for tc in tool_calls:
            if tc.success and tc.tool_output:
                output_str = json.dumps(tc.tool_output, indent=1, default=str)
                if len(output_str) > 3000:
                    output_str = output_str[:3000] + "..."
                evidence_summary_parts.append(
                    f"## {tc.tool_name}\n```json\n{output_str}\n```"
                )

        evidence_text = "\n\n".join(evidence_summary_parts) or "No evidence gathered."

        synthesis_prompt = f"""\
You were asked: {original_question}

You gathered the following evidence from {len(tool_calls)} tool calls:

{evidence_text}

Please synthesize a comprehensive clinical assessment based on the available evidence.
Include confidence level and note any gaps in the evidence.
Format your response in markdown with clear sections.
"""

        try:
            response = self._client.messages.create(
                model=self.model,
                max_tokens=4096,
                temperature=0.3,
                system=self.system_prompt,
                messages=[{"role": "user", "content": synthesis_prompt}],
            )
            return response.content[0].text
        except Exception as exc:
            logger.error("Synthesis call failed: %s", exc)
            return (
                f"**Unable to fully synthesize response.**\n\n"
                f"The agent gathered evidence from {len(tool_calls)} tool calls but "
                f"could not complete synthesis. Evidence collected from: "
                f"{', '.join(tc.tool_name for tc in tool_calls if tc.success)}.\n\n"
                f"Please try rephrasing your question or breaking it into smaller parts."
            )

    def _estimate_confidence(
        self,
        tool_calls: list[ToolCallRecord],
        sources: set[str],
    ) -> float:
        """
        Heuristic confidence estimation based on source diversity and tool success.

        Factors:
        - Number of distinct data sources consulted (genomic, cart, drug_discovery)
        - Tool call success rate
        - Presence of specific evidence types (variants found, structures available)
        """
        if not tool_calls:
            return 0.1

        # Base confidence from tool success rate
        successful = sum(1 for tc in tool_calls if tc.success)
        total = len(tool_calls)
        success_rate = successful / total if total > 0 else 0.0

        # Source diversity bonus
        domain_sources = {
            "genomic_evidence",
            "variant_annotations",
            "cart_evidence",
            "cart_knowledge",
            "drug_candidates",
            "protein_structures",
            "clinical_trials",
        }
        source_diversity = len(sources.intersection(domain_sources)) / len(domain_sources)

        # Check for substantive results
        has_substantive = False
        for tc in tool_calls:
            if tc.success and tc.tool_output:
                output = tc.tool_output
                if isinstance(output, list) and len(output) > 0:
                    has_substantive = True
                elif isinstance(output, dict):
                    if output.get("found") or output.get("pipeline_results_available"):
                        has_substantive = True
                    if output.get("structures") or output.get("top_candidates"):
                        has_substantive = True

        # Compute weighted confidence
        confidence = (
            0.3 * success_rate
            + 0.3 * source_diversity
            + 0.2 * (1.0 if has_substantive else 0.2)
            + 0.2 * min(1.0, successful / 3)  # Bonus for multiple successful calls
        )

        return min(1.0, max(0.05, confidence))

    def _extract_follow_ups(self, answer: str) -> list[str]:
        """
        Extract suggested follow-up questions from the answer text.

        Looks for patterns like "You might also consider..." or bullet points
        under a "Follow-up" or "Next steps" heading.
        """
        follow_ups: list[str] = []

        # Look for explicit follow-up sections
        patterns = [
            r"(?:follow[- ]?up|next steps|additional questions|you (?:might|may|could) (?:also )?(?:ask|consider|explore|investigate)).*?:\s*\n((?:[-*]\s+.+\n?)+)",
            r"(?:suggested|recommended) (?:follow[- ]?up|questions).*?:\s*\n((?:[-*]\s+.+\n?)+)",
        ]

        for pattern in patterns:
            matches = re.findall(pattern, answer, re.IGNORECASE | re.MULTILINE)
            for match in matches:
                lines = match.strip().split("\n")
                for line in lines:
                    cleaned = re.sub(r'^[-*]\s+', '', line.strip())
                    if cleaned and len(cleaned) > 10:
                        follow_ups.append(cleaned)

        # If none found, generate domain-appropriate suggestions
        if not follow_ups:
            answer_lower = answer.lower()
            if "variant" in answer_lower or "genomic" in answer_lower:
                follow_ups.append(
                    "What are the functional consequences of the identified variants?"
                )
            if "drug" in answer_lower or "molecule" in answer_lower:
                follow_ups.append(
                    "What are the ADMET properties of the top drug candidates?"
                )
            if "car-t" in answer_lower or "cell therapy" in answer_lower:
                follow_ups.append(
                    "What is the expected toxicity profile for this CAR-T approach?"
                )
            if "target" in answer_lower:
                follow_ups.append(
                    "Are there alternative therapeutic modalities for this target?"
                )

        return follow_ups[:5]

    @staticmethod
    def _tool_to_source(tool_name: str) -> str:
        """Map tool name to data source identifier."""
        mapping = {
            "search_genomic_variants": "genomic_evidence",
            "get_variant_annotations": "variant_annotations",
            "search_cart_evidence": "cart_evidence",
            "get_cart_knowledge": "cart_knowledge",
            "get_drug_candidates": "drug_candidates",
            "get_protein_structures": "protein_structures",
            "search_clinical_trials": "clinical_trials",
            "trigger_drug_discovery": "drug_discovery_pipeline",
        }
        return mapping.get(tool_name, tool_name)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _truncate_output(result: Any, max_items: int = 20) -> Any:
    """
    Truncate large tool outputs for the audit trail.

    Keeps JSON structure but limits list lengths and string sizes.
    """
    if isinstance(result, list):
        truncated = result[:max_items]
        if len(result) > max_items:
            truncated.append({"_truncated": f"{len(result) - max_items} more items"})
        return truncated
    elif isinstance(result, dict):
        truncated = {}
        for key, value in result.items():
            if isinstance(value, list) and len(value) > max_items:
                truncated[key] = value[:max_items] + [
                    {"_truncated": f"{len(value) - max_items} more items"}
                ]
            elif isinstance(value, str) and len(value) > 2000:
                truncated[key] = value[:2000] + "... (truncated)"
            else:
                truncated[key] = value
        return truncated
    return result


# ---------------------------------------------------------------------------
# Factory function
# ---------------------------------------------------------------------------

def create_meta_agent(
    anthropic_api_key: str | None = None,
    model: str = DEFAULT_MODEL,
    milvus_host: str | None = None,
    milvus_port: int | None = None,
    output_base_dir: str | None = None,
    **kwargs: Any,
) -> MetaAgent:
    """
    Factory function to create a MetaAgent with settings from HCLSSettings.

    Reads defaults from the centralized configuration and allows overrides.

    Args:
        anthropic_api_key: Anthropic API key (falls back to env var).
        model: Claude model ID.
        milvus_host: Milvus host (falls back to HCLS_MILVUS_HOST).
        milvus_port: Milvus port (falls back to HCLS_MILVUS_PORT).
        output_base_dir: Drug discovery output directory.
        **kwargs: Additional arguments passed to MetaAgent.

    Returns:
        Configured MetaAgent instance.
    """
    # Try to load from centralized settings
    try:
        from hcls_common.config import get_settings

        settings = get_settings()
        milvus_host = milvus_host or settings.milvus_host
        milvus_port = milvus_port or settings.milvus_port
        embedding_model = kwargs.pop("embedding_model", settings.embedding_model)

        # Use Anthropic key from settings if configured
        if not anthropic_api_key and settings.llm_api_key:
            anthropic_api_key = settings.llm_api_key
    except Exception:
        milvus_host = milvus_host or "localhost"
        milvus_port = milvus_port or 19530
        embedding_model = kwargs.pop("embedding_model", "BAAI/bge-small-en-v1.5")

    return MetaAgent(
        anthropic_api_key=anthropic_api_key,
        model=model,
        milvus_host=milvus_host,
        milvus_port=milvus_port,
        embedding_model=embedding_model,
        output_base_dir=output_base_dir or "outputs",
        **kwargs,
    )
