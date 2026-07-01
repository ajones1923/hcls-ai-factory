"""
Milvus Collection Manager for the Oncology Intelligence Agent.

Defines schemas for 10 owned collections + 1 read-only genomic_evidence
collection (shared with the genomics pipeline). Each collection uses
BGE-small-en-v1.5 embeddings (dim=384) with IVF_FLAT / COSINE indexing.

Collections:
    onco_literature   - PubMed / PMC / preprint chunks tagged by cancer type
    onco_trials       - ClinicalTrials.gov summaries with biomarker criteria
    onco_variants     - Actionable somatic / germline variants (CIViC, OncoKB)
    onco_biomarkers   - Predictive / prognostic biomarkers and assays
    onco_therapies    - Approved & investigational therapies with MOA
    onco_pathways     - Signaling pathways, cross-talk, and druggable nodes
    onco_guidelines   - NCCN / ASCO / ESMO guideline recommendations
    onco_resistance   - Resistance mechanisms and bypass strategies
    onco_outcomes     - Real-world treatment outcome records
    onco_cases        - De-identified patient case snapshots
    genomic_evidence  - Read-only VCF-derived evidence from Stage 1

Author: Adam Jones
Date:   February 2026
"""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Dict, List, Optional, Type

from pymilvus import (
    Collection,
    CollectionSchema,
    DataType,
    FieldSchema,
    connections,
    utility,
)

from src.models import (
    CaseSnapshot,
    OncologyBiomarker,
    OncologyGuideline,
    OncologyLiterature,
    OncologyPathway,
    OncologyTherapy,
    OncologyTrial,
    OncologyVariant,
    OutcomeRecord,
    ResistanceMechanism,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

EMBEDDING_DIM = 384  # BGE-small-en-v1.5 output dimension

# ---------------------------------------------------------------------------
# Field definitions – one list per collection
# ---------------------------------------------------------------------------

# 1. onco_literature --------------------------------------------------------

ONCO_LITERATURE_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, max_length=100, is_primary=True),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=EMBEDDING_DIM),
    FieldSchema(name="title", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="text_chunk", dtype=DataType.VARCHAR, max_length=3000),
    FieldSchema(name="source_type", dtype=DataType.VARCHAR, max_length=20),
    FieldSchema(name="year", dtype=DataType.INT64),
    FieldSchema(name="cancer_type", dtype=DataType.VARCHAR, max_length=50),
    FieldSchema(name="gene", dtype=DataType.VARCHAR, max_length=50),
    FieldSchema(name="variant", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="keywords", dtype=DataType.VARCHAR, max_length=1000),
    FieldSchema(name="journal", dtype=DataType.VARCHAR, max_length=200),
]

ONCO_LITERATURE_SCHEMA = CollectionSchema(
    fields=ONCO_LITERATURE_FIELDS,
    description="PubMed / PMC / preprint literature chunks for precision oncology",
)

# 2. onco_trials ------------------------------------------------------------

ONCO_TRIALS_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, max_length=20, is_primary=True),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=EMBEDDING_DIM),
    FieldSchema(name="title", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="text_summary", dtype=DataType.VARCHAR, max_length=3000),
    FieldSchema(name="phase", dtype=DataType.VARCHAR, max_length=30),
    FieldSchema(name="status", dtype=DataType.VARCHAR, max_length=30),
    FieldSchema(name="sponsor", dtype=DataType.VARCHAR, max_length=200),
    FieldSchema(name="cancer_types", dtype=DataType.VARCHAR, max_length=200),
    FieldSchema(name="biomarker_criteria", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="enrollment", dtype=DataType.INT64),
    FieldSchema(name="start_year", dtype=DataType.INT64),
    FieldSchema(name="outcome_summary", dtype=DataType.VARCHAR, max_length=2000),
]

ONCO_TRIALS_SCHEMA = CollectionSchema(
    fields=ONCO_TRIALS_FIELDS,
    description="ClinicalTrials.gov summaries with biomarker eligibility criteria",
)

# 3. onco_variants ----------------------------------------------------------

ONCO_VARIANTS_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, max_length=100, is_primary=True),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=EMBEDDING_DIM),
    FieldSchema(name="gene", dtype=DataType.VARCHAR, max_length=50),
    FieldSchema(name="variant_name", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="variant_type", dtype=DataType.VARCHAR, max_length=30),
    FieldSchema(name="cancer_type", dtype=DataType.VARCHAR, max_length=50),
    FieldSchema(name="evidence_level", dtype=DataType.VARCHAR, max_length=20),
    FieldSchema(name="drugs", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="civic_id", dtype=DataType.VARCHAR, max_length=20),
    FieldSchema(name="vrs_id", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="text_summary", dtype=DataType.VARCHAR, max_length=3000),
    FieldSchema(name="clinical_significance", dtype=DataType.VARCHAR, max_length=200),
    FieldSchema(name="allele_frequency", dtype=DataType.FLOAT),
]

ONCO_VARIANTS_SCHEMA = CollectionSchema(
    fields=ONCO_VARIANTS_FIELDS,
    description="Actionable somatic / germline variants from CIViC and OncoKB",
)

# 4. onco_biomarkers --------------------------------------------------------

ONCO_BIOMARKERS_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, max_length=100, is_primary=True),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=EMBEDDING_DIM),
    FieldSchema(name="name", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="biomarker_type", dtype=DataType.VARCHAR, max_length=30),
    FieldSchema(name="cancer_types", dtype=DataType.VARCHAR, max_length=200),
    FieldSchema(name="predictive_value", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="testing_method", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="clinical_cutoff", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="text_summary", dtype=DataType.VARCHAR, max_length=3000),
    FieldSchema(name="evidence_level", dtype=DataType.VARCHAR, max_length=20),
]

ONCO_BIOMARKERS_SCHEMA = CollectionSchema(
    fields=ONCO_BIOMARKERS_FIELDS,
    description="Predictive and prognostic biomarkers with assay details",
)

# 5. onco_therapies ---------------------------------------------------------

ONCO_THERAPIES_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, max_length=100, is_primary=True),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=EMBEDDING_DIM),
    FieldSchema(name="drug_name", dtype=DataType.VARCHAR, max_length=200),
    FieldSchema(name="category", dtype=DataType.VARCHAR, max_length=30),
    FieldSchema(name="targets", dtype=DataType.VARCHAR, max_length=200),
    FieldSchema(name="approved_indications", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="resistance_mechanisms", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="evidence_level", dtype=DataType.VARCHAR, max_length=20),
    FieldSchema(name="text_summary", dtype=DataType.VARCHAR, max_length=3000),
    FieldSchema(name="mechanism_of_action", dtype=DataType.VARCHAR, max_length=500),
]

ONCO_THERAPIES_SCHEMA = CollectionSchema(
    fields=ONCO_THERAPIES_FIELDS,
    description="Approved and investigational therapies with mechanism of action",
)

# 6. onco_pathways ----------------------------------------------------------

ONCO_PATHWAYS_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, max_length=100, is_primary=True),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=EMBEDDING_DIM),
    FieldSchema(name="name", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="key_genes", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="therapeutic_targets", dtype=DataType.VARCHAR, max_length=300),
    FieldSchema(name="cross_talk", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="text_summary", dtype=DataType.VARCHAR, max_length=3000),
]

ONCO_PATHWAYS_SCHEMA = CollectionSchema(
    fields=ONCO_PATHWAYS_FIELDS,
    description="Signaling pathways, cross-talk, and druggable nodes",
)

# 7. onco_guidelines --------------------------------------------------------

ONCO_GUIDELINES_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, max_length=100, is_primary=True),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=EMBEDDING_DIM),
    FieldSchema(name="org", dtype=DataType.VARCHAR, max_length=20),
    FieldSchema(name="cancer_type", dtype=DataType.VARCHAR, max_length=50),
    FieldSchema(name="version", dtype=DataType.VARCHAR, max_length=50),
    FieldSchema(name="year", dtype=DataType.INT64),
    FieldSchema(name="key_recommendations", dtype=DataType.VARCHAR, max_length=3000),
    FieldSchema(name="text_summary", dtype=DataType.VARCHAR, max_length=3000),
    FieldSchema(name="evidence_level", dtype=DataType.VARCHAR, max_length=20),
]

ONCO_GUIDELINES_SCHEMA = CollectionSchema(
    fields=ONCO_GUIDELINES_FIELDS,
    description="NCCN / ASCO / ESMO guideline recommendations",
)

# 8. onco_resistance --------------------------------------------------------

ONCO_RESISTANCE_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, max_length=100, is_primary=True),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=EMBEDDING_DIM),
    FieldSchema(name="primary_therapy", dtype=DataType.VARCHAR, max_length=200),
    FieldSchema(name="gene", dtype=DataType.VARCHAR, max_length=50),
    FieldSchema(name="mechanism", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="bypass_pathway", dtype=DataType.VARCHAR, max_length=200),
    FieldSchema(name="alternative_therapies", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="text_summary", dtype=DataType.VARCHAR, max_length=3000),
]

ONCO_RESISTANCE_SCHEMA = CollectionSchema(
    fields=ONCO_RESISTANCE_FIELDS,
    description="Resistance mechanisms and bypass strategies for targeted therapies",
)

# 9. onco_outcomes ----------------------------------------------------------

ONCO_OUTCOMES_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, max_length=100, is_primary=True),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=EMBEDDING_DIM),
    FieldSchema(name="case_id", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="therapy", dtype=DataType.VARCHAR, max_length=200),
    FieldSchema(name="cancer_type", dtype=DataType.VARCHAR, max_length=50),
    FieldSchema(name="response", dtype=DataType.VARCHAR, max_length=20),
    FieldSchema(name="duration_months", dtype=DataType.FLOAT),
    FieldSchema(name="toxicities", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="biomarkers_at_baseline", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="text_summary", dtype=DataType.VARCHAR, max_length=3000),
]

ONCO_OUTCOMES_SCHEMA = CollectionSchema(
    fields=ONCO_OUTCOMES_FIELDS,
    description="Real-world treatment outcome records for precision oncology",
)

# 10. onco_cases ------------------------------------------------------------

ONCO_CASES_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, max_length=100, is_primary=True),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=EMBEDDING_DIM),
    FieldSchema(name="patient_id", dtype=DataType.VARCHAR, max_length=100),
    FieldSchema(name="cancer_type", dtype=DataType.VARCHAR, max_length=50),
    FieldSchema(name="stage", dtype=DataType.VARCHAR, max_length=20),
    FieldSchema(name="variants", dtype=DataType.VARCHAR, max_length=1000),
    FieldSchema(name="biomarkers", dtype=DataType.VARCHAR, max_length=1000),
    FieldSchema(name="prior_therapies", dtype=DataType.VARCHAR, max_length=500),
    FieldSchema(name="text_summary", dtype=DataType.VARCHAR, max_length=3000),
]

ONCO_CASES_SCHEMA = CollectionSchema(
    fields=ONCO_CASES_FIELDS,
    description="De-identified patient case snapshots for similarity search",
)

# 11. genomic_evidence (read-only) ------------------------------------------

GENOMIC_EVIDENCE_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, max_length=200, is_primary=True),
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
    description="Read-only VCF-derived genomic evidence from Stage 1 pipeline",
)

# ---------------------------------------------------------------------------
# Schema & model registries
# ---------------------------------------------------------------------------

COLLECTION_SCHEMAS: Dict[str, CollectionSchema] = {
    "onco_literature": ONCO_LITERATURE_SCHEMA,
    "onco_trials": ONCO_TRIALS_SCHEMA,
    "onco_variants": ONCO_VARIANTS_SCHEMA,
    "onco_biomarkers": ONCO_BIOMARKERS_SCHEMA,
    "onco_therapies": ONCO_THERAPIES_SCHEMA,
    "onco_pathways": ONCO_PATHWAYS_SCHEMA,
    "onco_guidelines": ONCO_GUIDELINES_SCHEMA,
    "onco_resistance": ONCO_RESISTANCE_SCHEMA,
    "onco_outcomes": ONCO_OUTCOMES_SCHEMA,
    "onco_cases": ONCO_CASES_SCHEMA,
    "genomic_evidence": GENOMIC_EVIDENCE_SCHEMA,
}

COLLECTION_MODELS: Dict[str, Optional[Type]] = {
    "onco_literature": OncologyLiterature,
    "onco_trials": OncologyTrial,
    "onco_variants": OncologyVariant,
    "onco_biomarkers": OncologyBiomarker,
    "onco_therapies": OncologyTherapy,
    "onco_pathways": OncologyPathway,
    "onco_guidelines": OncologyGuideline,
    "onco_resistance": ResistanceMechanism,
    "onco_outcomes": OutcomeRecord,
    "onco_cases": CaseSnapshot,
    "genomic_evidence": None,  # read-only — populated by Stage 1
}


# ---------------------------------------------------------------------------
# OncoCollectionManager
# ---------------------------------------------------------------------------


class OncoCollectionManager:
    """Manages Milvus collections for the Oncology Intelligence Agent.

    Provides helpers for connecting, creating / dropping collections,
    inserting data, and running vector similarity searches across one
    or all collections in parallel.
    """

    # IVF_FLAT index parameters
    INDEX_PARAMS = {
        "metric_type": "COSINE",
        "index_type": "IVF_FLAT",
        "params": {"nlist": 1024},
    }

    # Search parameters
    SEARCH_PARAMS = {
        "metric_type": "COSINE",
        "params": {"nprobe": 16},
    }

    def __init__(
        self,
        host: str = "localhost",
        port: int = 19530,
        alias: str = "default",
    ) -> None:
        self.host = host
        self.port = port
        self.alias = alias
        self._collections: Dict[str, Collection] = {}

    # ------------------------------------------------------------------
    # Connection lifecycle
    # ------------------------------------------------------------------

    def connect(self) -> None:
        """Establish a connection to the Milvus server."""
        logger.info("Connecting to Milvus at %s:%s (alias=%s)", self.host, self.port, self.alias)
        connections.connect(alias=self.alias, host=self.host, port=self.port)
        logger.info("Connected to Milvus successfully.")

    def disconnect(self) -> None:
        """Disconnect from the Milvus server and clear cached handles."""
        logger.info("Disconnecting from Milvus (alias=%s)", self.alias)
        connections.disconnect(alias=self.alias)
        self._collections.clear()
        logger.info("Disconnected from Milvus.")

    def is_connected(self) -> bool:
        """Return True if there is an active Milvus connection."""
        try:
            return utility.has_collection("__ping__") is not None
        except Exception:
            return False

    def list_collections(self) -> List[str]:
        """Return names of all collections in the current Milvus instance."""
        return utility.list_collections()

    def get_collection_count(self, name: str) -> int:
        """Return entity count for a collection."""
        col = self.get_collection(name)
        return col.num_entities

    # ------------------------------------------------------------------
    # Collection management
    # ------------------------------------------------------------------

    def create_collection(self, name: str) -> Collection:
        """Create a single collection by name and build its IVF_FLAT index.

        Args:
            name: Collection name (must exist in COLLECTION_SCHEMAS).

        Returns:
            The newly created (or already existing) pymilvus Collection.
        """
        if name not in COLLECTION_SCHEMAS:
            raise ValueError(f"Unknown collection: {name}")

        if utility.has_collection(name):
            logger.info("Collection '%s' already exists — skipping creation.", name)
            col = Collection(name=name)
        else:
            logger.info("Creating collection '%s' …", name)
            col = Collection(name=name, schema=COLLECTION_SCHEMAS[name])
            # Build the vector index on the embedding field
            col.create_index(
                field_name="embedding",
                index_params=self.INDEX_PARAMS,
            )
            logger.info("Collection '%s' created with IVF_FLAT index.", name)

        self._collections[name] = col
        return col

    def create_all_collections(self) -> Dict[str, Collection]:
        """Create all collections defined in COLLECTION_SCHEMAS.

        Returns:
            Dict mapping collection name to pymilvus Collection handle.
        """
        logger.info("Creating all %d collections …", len(COLLECTION_SCHEMAS))
        for name in COLLECTION_SCHEMAS:
            self.create_collection(name)
        logger.info("All collections created.")
        return dict(self._collections)

    def drop_collection(self, name: str) -> None:
        """Drop a collection by name.

        Args:
            name: The collection to drop.
        """
        if utility.has_collection(name):
            logger.warning("Dropping collection '%s'.", name)
            utility.drop_collection(name)
            self._collections.pop(name, None)
        else:
            logger.info("Collection '%s' does not exist — nothing to drop.", name)

    def get_collection(self, name: str) -> Collection:
        """Retrieve a Collection handle, loading it if necessary.

        Args:
            name: Collection name.

        Returns:
            pymilvus Collection (loaded into memory).
        """
        if name in self._collections:
            return self._collections[name]

        if not utility.has_collection(name):
            raise ValueError(f"Collection '{name}' does not exist in Milvus.")

        col = Collection(name=name)
        col.load()
        self._collections[name] = col
        return col

    def get_collection_stats(self, name: str) -> Dict[str, Any]:
        """Return basic statistics for a collection.

        Args:
            name: Collection name.

        Returns:
            Dict with entity count and field names.
        """
        col = self.get_collection(name)
        return {
            "name": name,
            "num_entities": col.num_entities,
            "fields": [f.name for f in col.schema.fields],
        }

    # ------------------------------------------------------------------
    # Data operations
    # ------------------------------------------------------------------

    def insert_batch(
        self,
        name: str,
        data: List[Dict[str, Any]],
    ) -> int:
        """Insert a batch of records into a collection.

        Args:
            name: Target collection name.
            data: List of dicts whose keys match the collection field names.

        Returns:
            Number of successfully inserted entities.
        """
        col = self.get_collection(name)
        if not data:
            logger.warning("insert_batch called with empty data for '%s'.", name)
            return 0

        # Transpose list-of-dicts into dict-of-lists for pymilvus
        field_names = [f.name for f in col.schema.fields]
        insert_data = {fn: [] for fn in field_names}
        for record in data:
            for fn in field_names:
                insert_data[fn].append(record.get(fn))

        res = col.insert([insert_data[fn] for fn in field_names])
        col.flush()
        inserted = res.insert_count
        logger.info("Inserted %d entities into '%s'.", inserted, name)
        return inserted

    def insert(
        self,
        collection_name: str = "",
        records: Optional[List[Dict[str, Any]]] = None,
        *,
        data: Any = None,
        name: str = "",
    ) -> int:
        """Flexible insert: accepts a single dict or list of dicts.

        Supports both calling conventions:
            insert(collection_name, [records])   — from ingest pipelines
            insert(collection_name=..., data={}) — from case manager
        """
        col_name = collection_name or name
        payload = records if records is not None else data
        if isinstance(payload, dict):
            payload = [payload]
        return self.insert_batch(name=col_name, data=payload)

    # ------------------------------------------------------------------
    # Search operations
    # ------------------------------------------------------------------

    def search(
        self,
        name: str = "",
        query_vector: Optional[List[float]] = None,
        top_k: int = 10,
        output_fields: Optional[List[str]] = None,
        expr: Optional[str] = None,
        *,
        collection: str = "",
        vector: Optional[List[float]] = None,
        filters: Optional[Dict[str, Any]] = None,
    ) -> List[Dict[str, Any]]:
        """Run a vector similarity search against a single collection.

        Supports two calling conventions:
            search(name, query_vector, ...)           — direct API
            search(collection=, vector=, filters=...) — from RAG engine
        """
        # Normalise parameter aliases
        name = name or collection
        query_vector = query_vector if query_vector is not None else vector
        if query_vector is None:
            raise ValueError("query_vector (or vector=) is required")

        # Convert numpy array to list if needed
        if hasattr(query_vector, "tolist"):
            query_vector = query_vector.tolist()

        # Build expr from filters dict if provided
        if filters and not expr:
            parts = []
            for k, v in filters.items():
                if "__gte" in k:
                    field = k.replace("__gte", "")
                    parts.append(f"{field} >= {int(v)}")
                elif "__lte" in k:
                    field = k.replace("__lte", "")
                    parts.append(f"{field} <= {int(v)}")
                else:
                    safe_v = str(v).replace('"', '\\"')
                    parts.append(f'{k} == "{safe_v}"')
            if parts:
                expr = " and ".join(parts)

        col = self.get_collection(name)
        col.load()

        # Default output fields: everything except the embedding vector
        if output_fields is None:
            output_fields = [
                f.name for f in col.schema.fields if f.dtype != DataType.FLOAT_VECTOR
            ]

        results = col.search(
            data=[query_vector],
            anns_field="embedding",
            param=self.SEARCH_PARAMS,
            limit=top_k,
            output_fields=output_fields,
            expr=expr,
        )

        hits: List[Dict[str, Any]] = []
        for hit in results[0]:
            record = {field: hit.entity.get(field) for field in output_fields}
            record["_distance"] = hit.distance
            record["_collection"] = name
            hits.append(record)

        logger.debug("Search in '%s' returned %d hits.", name, len(hits))
        return hits

    def search_all(
        self,
        query_vector: List[float],
        top_k: int = 5,
        collections: Optional[List[str]] = None,
        expr: Optional[str] = None,
        max_workers: int = 6,
    ) -> List[Dict[str, Any]]:
        """Search multiple collections in parallel and merge results.

        Args:
            query_vector:   Query embedding (dim=384).
            top_k:          Results per collection.
            collections:    Subset of collection names to search.
                            Defaults to all registered collections.
            expr:           Optional boolean filter expression (applied to
                            every collection — use with care).
            max_workers:    Thread pool size for parallel search.

        Returns:
            Combined results sorted by ascending distance (best first).
        """
        targets = collections or list(COLLECTION_SCHEMAS.keys())
        all_hits: List[Dict[str, Any]] = []

        def _search_one(col_name: str) -> List[Dict[str, Any]]:
            try:
                return self.search(col_name, query_vector, top_k=top_k, expr=expr)
            except Exception:
                logger.exception("search_all: error searching '%s'", col_name)
                return []

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(_search_one, c): c for c in targets}
            for future in as_completed(futures):
                all_hits.extend(future.result())

        # Sort by cosine distance (lower = more similar)
        all_hits.sort(key=lambda h: h.get("_distance", float("inf")))
        logger.info(
            "search_all across %d collections returned %d total hits.",
            len(targets),
            len(all_hits),
        )
        return all_hits
