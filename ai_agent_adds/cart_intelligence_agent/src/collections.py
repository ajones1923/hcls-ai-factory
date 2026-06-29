"""Milvus collection management for CAR-T Intelligence Agent.

Manages 11 CAR-T collections (10 domain-specific + 1 read-only genomic):
  - cart_literature    — Published research & patents
  - cart_trials        — ClinicalTrials.gov records
  - cart_constructs    — CAR construct designs
  - cart_assays        — In vitro / in vivo assay results
  - cart_manufacturing — Manufacturing / CMC records
  - cart_safety        — Pharmacovigilance & post-market safety
  - cart_biomarkers    — Predictive & pharmacodynamic biomarkers
  - cart_regulatory    — FDA regulatory milestones
  - cart_sequences     — Molecular & structural data
  - cart_realworld     — Real-world evidence & outcomes

Follows the same pymilvus pattern as:
  rag-chat-pipeline/src/milvus_client.py (MilvusClient)

Author: Adam Jones
Date: February 2026
"""

import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Dict, List, Optional

from loguru import logger
from pymilvus import (
    Collection,
    CollectionSchema,
    DataType,
    FieldSchema,
    connections,
    utility,
)

from src.models import (
    AssayResult,
    BiomarkerRecord,
    CARConstruct,
    CARTLiterature,
    ClinicalTrial,
    ManufacturingRecord,
    RealWorldRecord,
    RegulatoryRecord,
    SafetyRecord,
    SequenceRecord,
)


# ═══════════════════════════════════════════════════════════════════════
# COLLECTION SCHEMA DEFINITIONS
# ═══════════════════════════════════════════════════════════════════════

EMBEDDING_DIM = 384  # BGE-small-en-v1.5

# ── cart_literature ──────────────────────────────────────────────────

LITERATURE_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="PMID or patent number",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="title",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Paper or patent title",
    ),
    FieldSchema(
        name="text_chunk",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Text chunk used for embedding",
    ),
    FieldSchema(
        name="source_type",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="pubmed, pmc, patent, preprint, manual",
    ),
    FieldSchema(
        name="year",
        dtype=DataType.INT64,
        description="Publication year",
    ),
    FieldSchema(
        name="cart_stage",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="target_id, car_design, vector_eng, testing, clinical",
    ),
    FieldSchema(
        name="target_antigen",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Target antigen (e.g. CD19, BCMA)",
    ),
    FieldSchema(
        name="disease",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Disease or indication",
    ),
    FieldSchema(
        name="keywords",
        dtype=DataType.VARCHAR,
        max_length=1000,
        description="Comma-separated keywords / MeSH terms",
    ),
    FieldSchema(
        name="journal",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Journal name",
    ),
]

LITERATURE_SCHEMA = CollectionSchema(
    fields=LITERATURE_FIELDS,
    description="CAR-T published literature and patents",
)

# ── cart_trials ──────────────────────────────────────────────────────

TRIALS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=20,
        description="NCT number (e.g. NCT03958656)",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="title",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Official trial title",
    ),
    FieldSchema(
        name="text_summary",
        dtype=DataType.VARCHAR,
        max_length=3000,
        description="Brief summary for embedding",
    ),
    FieldSchema(
        name="phase",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="Trial phase",
    ),
    FieldSchema(
        name="status",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="Recruitment status",
    ),
    FieldSchema(
        name="sponsor",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Lead sponsor",
    ),
    FieldSchema(
        name="target_antigen",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Target antigen",
    ),
    FieldSchema(
        name="car_generation",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="1st, 2nd, 3rd, 4th, armored, universal",
    ),
    FieldSchema(
        name="costimulatory",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="CD28, 4-1BB, dual",
    ),
    FieldSchema(
        name="disease",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Disease or indication",
    ),
    FieldSchema(
        name="enrollment",
        dtype=DataType.INT64,
        description="Target enrollment count",
    ),
    FieldSchema(
        name="start_year",
        dtype=DataType.INT64,
        description="Study start year",
    ),
    FieldSchema(
        name="outcome_summary",
        dtype=DataType.VARCHAR,
        max_length=2000,
        description="Outcome summary if available",
    ),
]

TRIALS_SCHEMA = CollectionSchema(
    fields=TRIALS_FIELDS,
    description="CAR-T clinical trials from ClinicalTrials.gov",
)

# ── cart_constructs ──────────────────────────────────────────────────

CONSTRUCTS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Construct identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="name",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Product or construct name",
    ),
    FieldSchema(
        name="text_summary",
        dtype=DataType.VARCHAR,
        max_length=2000,
        description="Description for embedding",
    ),
    FieldSchema(
        name="target_antigen",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Target antigen",
    ),
    FieldSchema(
        name="scfv_origin",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="scFv antibody clone / origin",
    ),
    FieldSchema(
        name="costimulatory_domain",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Costimulatory domain(s)",
    ),
    FieldSchema(
        name="signaling_domain",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Intracellular signaling domain",
    ),
    FieldSchema(
        name="generation",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="CAR generation",
    ),
    FieldSchema(
        name="hinge_tm",
        dtype=DataType.VARCHAR,
        max_length=200,
        description="Hinge + transmembrane domain",
    ),
    FieldSchema(
        name="vector_type",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Lentiviral, retroviral, non-viral, etc.",
    ),
    FieldSchema(
        name="fda_status",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="FDA approval status",
    ),
    FieldSchema(
        name="known_toxicities",
        dtype=DataType.VARCHAR,
        max_length=500,
        description="Known toxicities (CRS, ICANS, etc.)",
    ),
]

CONSTRUCTS_SCHEMA = CollectionSchema(
    fields=CONSTRUCTS_FIELDS,
    description="CAR construct designs and FDA-approved products",
)

# ── cart_assays ──────────────────────────────────────────────────────

ASSAYS_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Assay record identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="text_summary",
        dtype=DataType.VARCHAR,
        max_length=2000,
        description="Assay description for embedding",
    ),
    FieldSchema(
        name="assay_type",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="cytotoxicity, cytokine, flow, proliferation, in_vivo, etc.",
    ),
    FieldSchema(
        name="construct_id",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Foreign key to cart_constructs",
    ),
    FieldSchema(
        name="target_antigen",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Target antigen tested",
    ),
    FieldSchema(
        name="cell_line",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Target cell line (e.g. Nalm-6, Raji, K562)",
    ),
    FieldSchema(
        name="effector_ratio",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="Effector:target ratio",
    ),
    FieldSchema(
        name="key_metric",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Primary metric name (e.g. % lysis, IFN-gamma pg/mL)",
    ),
    FieldSchema(
        name="metric_value",
        dtype=DataType.FLOAT,
        description="Numeric value for key_metric",
    ),
    FieldSchema(
        name="outcome",
        dtype=DataType.VARCHAR,
        max_length=20,
        description="success, partial, failure",
    ),
    FieldSchema(
        name="notes",
        dtype=DataType.VARCHAR,
        max_length=1000,
        description="Additional notes",
    ),
]

ASSAYS_SCHEMA = CollectionSchema(
    fields=ASSAYS_FIELDS,
    description="CAR-T in vitro and in vivo assay results",
)

# ── cart_manufacturing ───────────────────────────────────────────────

MANUFACTURING_FIELDS = [
    FieldSchema(
        name="id",
        dtype=DataType.VARCHAR,
        is_primary=True,
        max_length=100,
        description="Manufacturing record identifier",
    ),
    FieldSchema(
        name="embedding",
        dtype=DataType.FLOAT_VECTOR,
        dim=EMBEDDING_DIM,
        description="BGE-small-en-v1.5 text embedding",
    ),
    FieldSchema(
        name="text_summary",
        dtype=DataType.VARCHAR,
        max_length=2000,
        description="Process description for embedding",
    ),
    FieldSchema(
        name="process_step",
        dtype=DataType.VARCHAR,
        max_length=30,
        description="transduction, expansion, harvest, formulation, release, cryo",
    ),
    FieldSchema(
        name="vector_type",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Vector type used",
    ),
    FieldSchema(
        name="parameter",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Process parameter name",
    ),
    FieldSchema(
        name="parameter_value",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Parameter value",
    ),
    FieldSchema(
        name="target_spec",
        dtype=DataType.VARCHAR,
        max_length=100,
        description="Acceptance criteria / specification",
    ),
    FieldSchema(
        name="met_spec",
        dtype=DataType.VARCHAR,
        max_length=10,
        description="yes, no, borderline",
    ),
    FieldSchema(
        name="batch_id",
        dtype=DataType.VARCHAR,
        max_length=50,
        description="Manufacturing batch identifier",
    ),
    FieldSchema(
        name="notes",
        dtype=DataType.VARCHAR,
        max_length=1000,
        description="Additional notes",
    ),
]

MANUFACTURING_SCHEMA = CollectionSchema(
    fields=MANUFACTURING_FIELDS,
    description="CAR-T manufacturing / CMC process records",
)

# ── cart_safety ────────────────────────────────────────────────────

SAFETY_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, is_primary=True, max_length=100),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=EMBEDDING_DIM),
    FieldSchema(name="text_summary", dtype=DataType.VARCHAR, max_length=3000),
    FieldSchema(name="product", dtype=DataType.VARCHAR, max_length=200),
    FieldSchema(name="event_type", dtype=DataType.VARCHAR, max_length=30),
    FieldSchema(name="severity_grade", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="onset_timing", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="incidence_rate", dtype=DataType.VARCHAR, max_length=200),
    FieldSchema(name="management_protocol", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="outcome", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="reporting_source", dtype=DataType.VARCHAR, max_length=50),
    FieldSchema(name="year", dtype=DataType.INT64),
]

SAFETY_SCHEMA = CollectionSchema(
    fields=SAFETY_FIELDS,
    description="CAR-T pharmacovigilance and post-market safety data",
)

# ── cart_biomarkers ────────────────────────────────────────────────

BIOMARKER_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, is_primary=True, max_length=100),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=EMBEDDING_DIM),
    FieldSchema(name="text_summary", dtype=DataType.VARCHAR, max_length=3000),
    FieldSchema(name="biomarker_name", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="biomarker_type", dtype=DataType.VARCHAR, max_length=30),
    FieldSchema(name="assay_method", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="clinical_cutoff", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="predictive_value", dtype=DataType.VARCHAR, max_length=200),
    FieldSchema(name="associated_outcome", dtype=DataType.VARCHAR, max_length=200),
    FieldSchema(name="target_antigen", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="disease", dtype=DataType.VARCHAR, max_length=200),
    FieldSchema(name="evidence_level", dtype=DataType.VARCHAR, max_length=20),
]

BIOMARKER_SCHEMA = CollectionSchema(
    fields=BIOMARKER_FIELDS,
    description="CAR-T predictive and pharmacodynamic biomarkers",
)

# ── cart_regulatory ────────────────────────────────────────────────

REGULATORY_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, is_primary=True, max_length=100),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=EMBEDDING_DIM),
    FieldSchema(name="text_summary", dtype=DataType.VARCHAR, max_length=3000),
    FieldSchema(name="product", dtype=DataType.VARCHAR, max_length=200),
    FieldSchema(name="regulatory_event", dtype=DataType.VARCHAR, max_length=50),
    FieldSchema(name="date", dtype=DataType.VARCHAR, max_length=20),
    FieldSchema(name="agency", dtype=DataType.VARCHAR, max_length=20),
    FieldSchema(name="indication", dtype=DataType.VARCHAR, max_length=200),
    FieldSchema(name="decision", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="conditions", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="pivotal_trial", dtype=DataType.VARCHAR, max_length=100),
]

REGULATORY_SCHEMA = CollectionSchema(
    fields=REGULATORY_FIELDS,
    description="CAR-T FDA regulatory milestones and approvals",
)

# ── cart_sequences ─────────────────────────────────────────────────

SEQUENCE_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, is_primary=True, max_length=100),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=EMBEDDING_DIM),
    FieldSchema(name="text_summary", dtype=DataType.VARCHAR, max_length=3000),
    FieldSchema(name="construct_name", dtype=DataType.VARCHAR, max_length=200),
    FieldSchema(name="target_antigen", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="scfv_clone", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="binding_affinity_kd", dtype=DataType.VARCHAR, max_length=50),
    FieldSchema(name="heavy_chain_vregion", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="light_chain_vregion", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="framework", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="species_origin", dtype=DataType.VARCHAR, max_length=30),
    FieldSchema(name="immunogenicity_risk", dtype=DataType.VARCHAR, max_length=20),
    FieldSchema(name="structural_notes", dtype=DataType.VARCHAR, max_length=1000),
]

SEQUENCE_SCHEMA = CollectionSchema(
    fields=SEQUENCE_FIELDS,
    description="CAR-T molecular and structural data (scFv, binding affinity)",
)

# ── cart_realworld ─────────────────────────────────────────────────

REALWORLD_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, is_primary=True, max_length=100),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=EMBEDDING_DIM),
    FieldSchema(name="text_summary", dtype=DataType.VARCHAR, max_length=3000),
    FieldSchema(name="study_type", dtype=DataType.VARCHAR, max_length=30),
    FieldSchema(name="data_source", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="product", dtype=DataType.VARCHAR, max_length=200),
    FieldSchema(name="indication", dtype=DataType.VARCHAR, max_length=200),
    FieldSchema(name="population_size", dtype=DataType.INT64),
    FieldSchema(name="median_followup_months", dtype=DataType.FLOAT),
    FieldSchema(name="primary_endpoint", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="outcome_value", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="setting", dtype=DataType.VARCHAR, max_length=50),
    FieldSchema(name="special_population", dtype=DataType.VARCHAR, max_length=200),
]

REALWORLD_SCHEMA = CollectionSchema(
    fields=REALWORLD_FIELDS,
    description="CAR-T real-world evidence and outcomes",
)


# ── Genomic Evidence (read-only, created by rag-chat-pipeline) ──────

GENOMIC_EVIDENCE_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, is_primary=True, max_length=200),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=EMBEDDING_DIM),
    FieldSchema(name="chrom", dtype=DataType.VARCHAR, max_length=10),
    FieldSchema(name="pos", dtype=DataType.INT64),
    FieldSchema(name="ref", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="alt", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="qual", dtype=DataType.FLOAT),
    FieldSchema(name="gene", dtype=DataType.VARCHAR, max_length=50),
    FieldSchema(name="consequence", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="impact", dtype=DataType.VARCHAR, max_length=20),
    FieldSchema(name="genotype", dtype=DataType.VARCHAR, max_length=10),
    FieldSchema(name="text_summary", dtype=DataType.VARCHAR, max_length=2000),
    FieldSchema(name="clinical_significance", dtype=DataType.VARCHAR, max_length=200),
    FieldSchema(name="rsid", dtype=DataType.VARCHAR, max_length=20),
    FieldSchema(name="disease_associations", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="am_pathogenicity", dtype=DataType.FLOAT),
    FieldSchema(name="am_class", dtype=DataType.VARCHAR, max_length=30),
]

GENOMIC_EVIDENCE_SCHEMA = CollectionSchema(
    fields=GENOMIC_EVIDENCE_FIELDS,
    description="Genomic variant evidence (read-only, from rag-chat-pipeline)",
)


# ═══════════════════════════════════════════════════════════════════════
# COLLECTION REGISTRY
# ═══════════════════════════════════════════════════════════════════════

COLLECTION_SCHEMAS: Dict[str, CollectionSchema] = {
    "cart_literature": LITERATURE_SCHEMA,
    "cart_trials": TRIALS_SCHEMA,
    "cart_constructs": CONSTRUCTS_SCHEMA,
    "cart_assays": ASSAYS_SCHEMA,
    "cart_manufacturing": MANUFACTURING_SCHEMA,
    "cart_safety": SAFETY_SCHEMA,
    "cart_biomarkers": BIOMARKER_SCHEMA,
    "cart_regulatory": REGULATORY_SCHEMA,
    "cart_sequences": SEQUENCE_SCHEMA,
    "cart_realworld": REALWORLD_SCHEMA,
    "genomic_evidence": GENOMIC_EVIDENCE_SCHEMA,
}

# Maps collection names to their Pydantic model class for validation
# genomic_evidence is None because it's read-only (no inserts from this agent)
COLLECTION_MODELS: Dict[str, type] = {
    "cart_literature": CARTLiterature,
    "cart_trials": ClinicalTrial,
    "cart_constructs": CARConstruct,
    "cart_assays": AssayResult,
    "cart_manufacturing": ManufacturingRecord,
    "cart_safety": SafetyRecord,
    "cart_biomarkers": BiomarkerRecord,
    "cart_regulatory": RegulatoryRecord,
    "cart_sequences": SequenceRecord,
    "cart_realworld": RealWorldRecord,
    "genomic_evidence": None,
}


# ═══════════════════════════════════════════════════════════════════════
# COLLECTION MANAGER
# ═══════════════════════════════════════════════════════════════════════


class CARTCollectionManager:
    """Manages 11 CAR-T Milvus collections (10 owned + 1 read-only genomic).

    Provides create/drop/insert/search operations across the full set of
    CAR-T domain collections, following the same pymilvus patterns as
    rag-chat-pipeline/src/milvus_client.py.

    Usage:
        manager = CARTCollectionManager()
        manager.connect()
        manager.create_all_collections()
        stats = manager.get_collection_stats()
    """

    # IVF_FLAT index params shared across all collections
    INDEX_PARAMS = {
        "metric_type": "COSINE",
        "index_type": "IVF_FLAT",
        "params": {"nlist": 1024},
    }

    SEARCH_PARAMS = {
        "metric_type": "COSINE",
        "params": {"nprobe": 16},
    }

    def __init__(
        self,
        host: Optional[str] = None,
        port: Optional[int] = None,
        embedding_dim: int = EMBEDDING_DIM,
    ):
        """Initialize the collection manager.

        Args:
            host: Milvus server host. Defaults to MILVUS_HOST env var or localhost.
            port: Milvus server port. Defaults to MILVUS_PORT env var or 19530.
            embedding_dim: Embedding vector dimension (384 for BGE-small-en-v1.5).
        """
        self.host = host or os.environ.get("MILVUS_HOST", "localhost")
        self.port = port or int(os.environ.get("MILVUS_PORT", "19530"))
        self.embedding_dim = embedding_dim
        self._collections: Dict[str, Collection] = {}

    def connect(self) -> None:
        """Connect to the Milvus server."""
        logger.info(f"Connecting to Milvus at {self.host}:{self.port}")
        connections.connect(
            alias="default",
            host=self.host,
            port=self.port,
        )
        logger.info("Connected to Milvus")

    def disconnect(self) -> None:
        """Disconnect from the Milvus server."""
        connections.disconnect("default")
        self._collections.clear()
        logger.info("Disconnected from Milvus")

    # ── Collection lifecycle ─────────────────────────────────────────

    def create_collection(
        self,
        name: str,
        schema: CollectionSchema,
        drop_existing: bool = False,
    ) -> Collection:
        """Create a single collection with IVF_FLAT index on the embedding field.

        Args:
            name: Collection name (must be a recognized CAR-T or genomic collection).
            schema: The CollectionSchema defining the fields.
            drop_existing: If True, drop the collection first if it already exists.

        Returns:
            The pymilvus Collection object.
        """
        if drop_existing and utility.has_collection(name):
            logger.warning(f"Dropping existing collection: {name}")
            utility.drop_collection(name)

        if utility.has_collection(name):
            logger.info(f"Collection '{name}' already exists, loading reference")
            collection = Collection(name)
            self._collections[name] = collection
            return collection

        logger.info(f"Creating collection: {name}")
        collection = Collection(name=name, schema=schema)

        # Create IVF_FLAT index on the embedding field
        logger.info(f"Creating IVF_FLAT/COSINE index on '{name}.embedding'")
        collection.create_index(
            field_name="embedding",
            index_params=self.INDEX_PARAMS,
        )

        self._collections[name] = collection
        logger.info(f"Collection '{name}' created with index")
        return collection

    def create_all_collections(self, drop_existing: bool = False) -> Dict[str, Collection]:
        """Create all 11 CAR-T collections (10 domain + 1 read-only genomic).

        Args:
            drop_existing: If True, drop and recreate each collection.

        Returns:
            Dict mapping collection name to Collection object.
        """
        logger.info(f"Creating all {len(COLLECTION_SCHEMAS)} CAR-T collections")
        for name, schema in COLLECTION_SCHEMAS.items():
            self.create_collection(name, schema, drop_existing=drop_existing)
        logger.info(f"All {len(COLLECTION_SCHEMAS)} collections ready")
        return dict(self._collections)

    def drop_collection(self, name: str) -> None:
        """Drop a collection by name.

        Args:
            name: The collection name to drop.
        """
        if utility.has_collection(name):
            utility.drop_collection(name)
            self._collections.pop(name, None)
            logger.info(f"Collection '{name}' dropped")
        else:
            logger.warning(f"Collection '{name}' does not exist, nothing to drop")

    def get_collection(self, name: str) -> Collection:
        """Get a collection reference, creating it if needed.

        Args:
            name: The collection name.

        Returns:
            The pymilvus Collection object.

        Raises:
            ValueError: If the name is not a recognized CAR-T collection.
        """
        if name in self._collections:
            return self._collections[name]

        if utility.has_collection(name):
            collection = Collection(name)
            self._collections[name] = collection
            return collection

        if name in COLLECTION_SCHEMAS:
            return self.create_collection(name, COLLECTION_SCHEMAS[name])

        raise ValueError(
            f"Unknown collection '{name}'. "
            f"Valid collections: {list(COLLECTION_SCHEMAS.keys())}"
        )

    # ── Stats ────────────────────────────────────────────────────────

    def get_collection_stats(self) -> Dict[str, int]:
        """Get row counts for all 10 CAR-T collections.

        Returns:
            Dict mapping collection name to entity count.
            Collections that do not yet exist will show 0.
        """
        stats: Dict[str, int] = {}
        for name in COLLECTION_SCHEMAS:
            if utility.has_collection(name):
                collection = Collection(name)
                stats[name] = collection.num_entities
            else:
                stats[name] = 0
        return stats

    # ── Data operations ──────────────────────────────────────────────

    def _get_output_fields(self, collection_name: str) -> List[str]:
        """Return non-embedding field names for a given collection.

        Used to build the output_fields list for search results.
        Excludes the 'embedding' field since it is large and not
        needed in result payloads.

        Args:
            collection_name: The collection to get fields for.

        Returns:
            List of field name strings (e.g. ["id", "title", "text_chunk", ...]).

        Raises:
            ValueError: If the collection_name is not recognized.
        """
        if collection_name not in COLLECTION_SCHEMAS:
            raise ValueError(
                f"Unknown collection '{collection_name}'. "
                f"Valid collections: {list(COLLECTION_SCHEMAS.keys())}"
            )

        schema = COLLECTION_SCHEMAS[collection_name]
        return [
            field.name
            for field in schema.fields
            if field.name != "embedding"
        ]

    def insert_batch(
        self,
        collection_name: str,
        records: List[Dict[str, Any]],
    ) -> int:
        """Insert a batch of records into a collection.

        Each record dict must contain all required fields for the collection
        schema, including the pre-computed 'embedding' vector.

        Args:
            collection_name: Target collection name.
            records: List of dicts with field names matching the schema.

        Returns:
            Number of records successfully inserted.
        """
        try:
            collection = self.get_collection(collection_name)
            result = collection.insert(records)
            collection.flush()
            count = result.insert_count
            logger.info(f"Inserted {count} records into {collection_name}")
            return count
        except Exception as e:
            logger.error(f"Failed to insert batch into {collection_name}: {e}")
            raise

    def search(
        self,
        collection_name: str,
        query_embedding: List[float],
        top_k: int = 5,
        filter_expr: Optional[str] = None,
        score_threshold: float = 0.0,
    ) -> List[Dict[str, Any]]:
        """Search a single collection by vector similarity.

        Args:
            collection_name: The collection to search.
            query_embedding: 384-dim query vector (BGE-small-en-v1.5).
            top_k: Maximum number of results to return.
            filter_expr: Optional Milvus boolean filter expression
                (e.g. 'target_antigen == "CD19"').
            score_threshold: Minimum cosine similarity score (0.0-1.0).

        Returns:
            List of dicts with 'id', 'score', 'collection', and all output fields.
        """
        try:
            collection = self.get_collection(collection_name)
            collection.load()

            output_fields = self._get_output_fields(collection_name)

            results = collection.search(
                data=[query_embedding],
                anns_field="embedding",
                param=self.SEARCH_PARAMS,
                limit=top_k,
                output_fields=output_fields,
                expr=filter_expr,
            )

            # Convert results to list of dicts
            evidence_results: List[Dict[str, Any]] = []
            for hits in results:
                for hit in hits:
                    score = hit.score  # Cosine similarity (0-1)
                    if score < score_threshold:
                        continue

                    record: Dict[str, Any] = {
                        "id": hit.id,
                        "score": score,
                        "collection": collection_name,
                    }
                    for field_name in output_fields:
                        if field_name != "id":  # Already captured above
                            record[field_name] = hit.entity.get(field_name)

                    evidence_results.append(record)

            return evidence_results

        except Exception as e:
            logger.error(f"Search failed on {collection_name}: {e}")
            return []

    def search_all(
        self,
        query_embedding: List[float],
        top_k_per_collection: int = 5,
        filter_exprs: Optional[Dict[str, str]] = None,
        score_threshold: float = 0.0,
    ) -> Dict[str, List[Dict[str, Any]]]:
        """Search ALL CAR-T collections in parallel.

        Performs vector similarity search across every collection
        concurrently using a thread pool, then merges results.

        Args:
            query_embedding: 384-dim query vector (BGE-small-en-v1.5).
            top_k_per_collection: Max results per collection.
            filter_exprs: Optional dict of collection_name -> filter expression.
                Collections not in the dict get no filter.
            score_threshold: Minimum cosine similarity score (0.0-1.0).

        Returns:
            Dict mapping collection name -> list of result dicts.
        """
        collections = list(COLLECTION_SCHEMAS.keys())
        all_results: Dict[str, List[Dict[str, Any]]] = {}

        def _search_one(name: str) -> tuple:
            expr = (filter_exprs or {}).get(name)
            return name, self.search(
                collection_name=name,
                query_embedding=query_embedding,
                top_k=top_k_per_collection,
                filter_expr=expr,
                score_threshold=score_threshold,
            )

        with ThreadPoolExecutor(max_workers=len(collections)) as executor:
            futures = {
                executor.submit(_search_one, name): name
                for name in collections
            }
            for future in as_completed(futures):
                coll_name = futures[future]
                try:
                    name, hits = future.result()
                    all_results[name] = hits
                except Exception as e:
                    logger.warning(
                        f"Search failed for collection '{coll_name}': {e}"
                    )
                    all_results[coll_name] = []

        total = sum(len(v) for v in all_results.values())
        logger.info(
            f"Searched {len(collections)} collections, found {total} results"
        )
        return all_results
